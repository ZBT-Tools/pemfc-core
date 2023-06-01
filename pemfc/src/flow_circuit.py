# general imports
import numpy as np
from abc import ABC, abstractmethod

# local module imports
from . import interpolation as ip, global_functions as g_func, \
    channel as chl, output_object as oo
from .fluid import fluid as fluids


class ParallelFlowCircuit(ABC, oo.OutputObject):
    def __new__(cls, dict_flow_circuit, manifolds, channels,
                n_subchannels=1.0, **kwargs):
        circuit_type = dict_flow_circuit.get('type', 'Koh')
        if circuit_type == 'Koh':
            return super(ParallelFlowCircuit, cls).\
                __new__(KohFlowCircuit)
        elif circuit_type == 'ModifiedKoh':
            return super(ParallelFlowCircuit, cls).\
                __new__(ModifiedKohFlowCircuit)
        elif circuit_type == 'UpdatedKoh':
            return super(ParallelFlowCircuit, cls).\
                __new__(UpdatedKohFlowCircuit)
        elif circuit_type == 'Wang':
            return super(ParallelFlowCircuit, cls).\
                __new__(WangFlowCircuit)
        else:
            raise NotImplementedError

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0, **kwargs):
        name = dict_flow_circuit['name']
        super().__init__(name)
        assert isinstance(dict_flow_circuit, dict)
        assert isinstance(manifolds, (list, tuple))
        assert isinstance(channels,  (list, tuple))
        err_message = 'manifolds must be tuple or list with two objects of ' \
                      'class Channel'
        if len(manifolds) != 2:
            raise ValueError(err_message)
        elif not isinstance(manifolds[0], chl.Channel):
            raise TypeError(err_message)
        if not isinstance(channels[0], chl.Channel):
            raise TypeError(err_message)
        self.print_variables = \
            {
                'names': ['normalized_flow_distribution'],
                'units': ['-'],
                'sub_names': ['None']
            }
        self.combine_print_variables(self.print_variables,
                                     kwargs.get('print_variables', None))
        self.manifolds = manifolds
        self.manifolds[0].name = self.name + ': Inlet Manifold'
        self.manifolds[0].fluid.name = self.manifolds[0].name + ': ' \
            + self.manifolds[0].fluid.TYPE_NAME
        self.manifolds[1].name = self.name + ': Outlet Manifold'
        self.manifolds[1].fluid.name = self.manifolds[1].name + ': ' \
            + self.manifolds[1].fluid.TYPE_NAME
        self.channels = channels
        self.manifolds[0].flow_direction = 1
        self.shape = dict_flow_circuit.get('shape', 'U')
        if self.shape not in ('U', 'Z'):
            raise ValueError('shape of flow circuit must be either U or Z')
        if self.shape == 'U':
            self.manifolds[1].flow_direction = -1
        else:
            self.manifolds[1].flow_direction = 1

        if hasattr(self.manifolds[0].fluid, 'mass_fraction'):
            self.multi_component = True
        else:
            self.multi_component = False

        self.n_channels = len(self.channels)
        self.n_subchannels = n_subchannels
        self.tolerance = dict_flow_circuit.get('tolerance', 1e-6)
        self.max_iter = dict_flow_circuit.get('max_iter', 20)
        self.min_iter = dict_flow_circuit.get('min_iter', 3)
        self.calc_distribution = \
            dict_flow_circuit.get('calc_distribution', True)
        self.mass_flow_in = \
            self.manifolds[0].mass_flow_total[self.manifolds[0].id_in]
        self.vol_flow_in = 0.0
        self.channel_mass_flow = \
            np.ones(self.n_channels) * self.mass_flow_in / self.n_channels
        self.channel_vol_flow = np.zeros(self.channel_mass_flow.shape)
        self.channel_vol_flow_old = np.zeros(self.channel_vol_flow.shape)

        self.channel_length = \
            np.asarray([channel.length for channel in channels])
        self.channel_cross_area = \
            np.asarray([channel.cross_area for channel in channels])
        self.initialize = True
        self.update_channels(update_fluid=True)
        self.normalized_flow_distribution = \
            np.zeros(self.channel_vol_flow.shape)
        self.iteration = 0
        self.add_print_variables(self.print_variables)

    def update(self, inlet_mass_flow=None, calc_distribution=None,
               update_fluid=False, **kwargs):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.vol_flow_in = self.calc_volume_flow(mass_flow=inlet_mass_flow)
            if self.initialize:
                # homogeneous distribution
                self.channel_mass_flow[:] = inlet_mass_flow / self.n_channels
                self.channel_vol_flow[:] = self.vol_flow_in / self.n_channels
            else:
                # use previous distribution scaled to new mass flow
                self.channel_mass_flow[:] *= inlet_mass_flow / self.mass_flow_in
                self.channel_vol_flow[:] *= inlet_mass_flow / self.mass_flow_in
            # set new mass and volume flows
            self.mass_flow_in = inlet_mass_flow

        if self.initialize:
            self.update_channels(update_fluid=True)
            self.channel_vol_flow_old[:] = 1e8
        # channel_vol_flow_old = np.zeros(self.channel_vol_flow.shape)
        if calc_distribution is None:
            calc_distribution = self.calc_distribution
        if calc_distribution and self.n_channels > 1:
            for i in range(self.max_iter):
                # print(self.name + ' Iteration # ', str(i+1))
                self.iteration = i
                self.single_loop()
                if i == 0:
                    self.initialize = False

                error = \
                    np.sum(
                        np.divide(self.channel_vol_flow -
                                  self.channel_vol_flow_old[:],
                                  self.channel_vol_flow,
                                  where=self.channel_vol_flow != 0.0) ** 2.0)
                # print(channel_vol_flow_old)
                # print(self.channel_vol_flow)
                self.channel_vol_flow_old[:] = self.channel_vol_flow
                # print(error)
                if error < self.tolerance and i >= self.min_iter:
                    break
                if i == (self.max_iter - 1):
                    print('maximum number of iterations n = {} '
                          'with error = {} in update() of {} '
                          'reached'.format(self.max_iter, error, self))
        else:
            self.initialize = False
        # final channel updates within flow circuit iteration
        self.update_channels(update_fluid=True)
        try:
            self.normalized_flow_distribution[:] = \
                self.channel_mass_flow / np.average(self.channel_mass_flow)
        except FloatingPointError:
            self.normalized_flow_distribution[:] = 0.0

    @abstractmethod
    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        pass

    def update_channels(self, update_fluid=True):
        if self.initialize:
            channel_mass_flow_in = np.ones(self.n_channels) \
                * self.mass_flow_in / self.n_channels
            channel_mass_flow_out = channel_mass_flow_in
        else:
            channel_mass_flow_in = self.channel_mass_flow
            # channel_mass_flow_out = \
            #     np.array([channel.mass_flow_total[channel.id_out]
            #               for channel in self.channels])
            # channel_mass_flow_out *= self.n_subchannels
            channel_mass_flow_out = self.channel_mass_flow

        if self.multi_component:
            mass_fraction = \
                np.array([channel.fluid.mass_fraction[:, channel.id_out]
                          for channel in self.channels]).transpose()
        else:
            mass_fraction = 1.0

        mass_source = channel_mass_flow_out * mass_fraction
        # mass_source = self.channel_mass_flow * mass_fraction
        channel_enthalpy_out = \
            np.asarray([ch.g_fluid[ch.id_out] * ch.temperature[ch.id_out]
                        for ch in self.channels]) * self.n_subchannels

        p_channel_out = np.zeros(len(self.manifolds[1].pressure) - 1)
        if self.calc_distribution:
            self.manifolds[1].update(
                mass_flow_in=0.0, mass_source=mass_source,
                update_mass=True, update_flow=True,
                update_heat=False, update_fluid=update_fluid,
                enthalpy_source=channel_enthalpy_out)
            p_channel_out[:] = ip.interpolate_1d(self.manifolds[1].pressure)
        else:
            p_channel_out[:] = self.manifolds[1].p_out

        # Channel update
        for i, channel in enumerate(self.channels):
            channel.p_out = p_channel_out[i]
            channel.temperature[channel.id_in] = self.manifolds[0].temp_ele[i]
            channel.update(mass_flow_in=
                           channel_mass_flow_in[i] / self.n_subchannels,
                           update_mass=True, update_flow=True,
                           update_heat=False, update_fluid=update_fluid)

        # Inlet header update
        id_in = self.channels[-1].id_in
        self.manifolds[0].p_out = self.channels[-1].pressure[id_in]
        if self.multi_component:
            mass_fraction = self.manifolds[0].fluid.mass_fraction[:, :-1]
        else:
            mass_fraction = 1.0
        mass_source = -self.channel_mass_flow * mass_fraction
        if self.calc_distribution:
            self.manifolds[0].update(
                mass_flow_in=self.mass_flow_in,  # * 1.00000,
                mass_source=mass_source,
                update_mass=True, update_flow=True,
                update_heat=False, update_fluid=update_fluid)
        self.vol_flow_in = self.calc_volume_flow()


    def calc_volume_flow(self, mass_flow=None):
        if mass_flow is None:
            mass_flow = self.mass_flow_in
        if self.calc_distribution:
            id_in = self.manifolds[0].id_in
            density = self.manifolds[0].fluid.density[id_in]
        else:
            density = \
                np.average([channel.fluid.density[channel.id_in]
                            for channel in self.channels])
        return mass_flow / density


class KohFlowCircuit(ParallelFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)
        # Distribution factor
        self.alpha = np.ones(self.n_channels)
        id_in = self.channels[-1].id_in
        id_out = self.channels[-1].id_out
        self.dp_ref = \
            self.channels[-1].pressure[id_in] \
            - self.channels[-1].pressure[id_out]
        self.k_perm = np.zeros(self.n_channels)
        self.l_by_a = np.array([channel.length / channel.cross_area
                                for channel in self.channels])
        self.visc_channel = np.zeros(self.n_channels)
        self.dp_channel = np.zeros(self.n_channels)

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
        if update_channels:
            self.update_channels()
            self.dp_channel[:] = \
                np.array([channel.pressure[channel.id_in]
                          - channel.pressure[channel.id_out]
                          for channel in self.channels])
            self.channel_vol_flow[:] = \
                np.array([np.average(channel.vol_flow)
                          for channel in self.channels])
        # if np.min(np.abs(vol_flow_channel)) > g_par.SMALL:
            self.visc_channel[:] = \
                np.array([np.average(channel.fluid.viscosity)
                          for channel in self.channels])
        # velocity = np.array([np.average(channel.velocity)
        #                      for channel in self.channels])
        p_in = ip.interpolate_1d(self.manifolds[0].pressure)
        p_out = ip.interpolate_1d(self.manifolds[1].pressure)
        if np.any(self.channel_vol_flow == 0.0):
            raise ValueError('zero flow rates detected, '
                             'check boundary conditions')
        if self.initialize:
            self.k_perm[:] = self.channel_vol_flow / self.dp_channel \
                * self.visc_channel * self.l_by_a
        self.dp_ref = np.maximum(self.dp_channel[-1], 1e-3)
        self.alpha[:] = (p_in - p_out) / self.dp_ref
        try:
            self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a \
                * self.visc_channel[-1] / self.k_perm[-1] / self.n_subchannels
        except FloatingPointError:
            raise FloatingPointError('check if geomtries are adequate '
                                     'for flow conditions in {}'.
                                     format(self.name))

        p_in += self.dp_ref \
            + self.manifolds[1].pressure[self.manifolds[1].id_out] \
            - self.manifolds[0].p_out
        self.alpha[:] = (p_in - p_out) / self.dp_ref
        self.channel_vol_flow[:] = (p_in - p_out) * self.k_perm \
            / self.l_by_a * self.n_subchannels / self.visc_channel
        density = np.array([channel.fluid.density[channel.id_in]
                            for channel in self.channels])
        self.channel_mass_flow[:] = self.channel_vol_flow * density
        mass_flow_correction = \
            self.mass_flow_in / np.sum(self.channel_mass_flow)
        self.channel_mass_flow[:] *= mass_flow_correction


class ModifiedKohFlowCircuit(KohFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)
        self.urf = dict_flow_circuit.get('underrelaxation_factor', 0.5)

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
        if update_channels:
            self.update_channels()
            self.dp_channel[:] = \
                np.array([channel.pressure[channel.id_in]
                          - channel.pressure[channel.id_out]
                          for channel in self.channels])
            self.channel_vol_flow[:] = self.n_subchannels \
                * np.array([np.average(channel.vol_flow)
                            for channel in self.channels])
        # if np.min(np.abs(vol_flow_channel)) > g_par.SMALL:
            self.visc_channel[:] = \
                np.array([np.average(channel.fluid.viscosity)
                          for channel in self.channels])
        # velocity = np.array([np.average(channel.velocity)
        #                      for channel in self.channels])
        p_in = ip.interpolate_1d(self.manifolds[0].pressure)
        p_out = ip.interpolate_1d(self.manifolds[1].pressure)
        if np.any(self.channel_vol_flow == 0.0):
            raise ValueError('zero flow rates detected, '
                             'check boundary conditions')
        # if self.initialize:
        #     self.k_perm[:] = self.channel_vol_flow / self.dp_channel \
        #         * self.visc_channel * self.l_by_a
        self.alpha[:] = (p_in - p_out) / self.dp_channel
        self.channel_vol_flow[:] *= (self.urf + (1.0 - self.urf) * self.alpha)
        density = np.array([channel.fluid.density[channel.id_in]
                            for channel in self.channels])
        self.channel_mass_flow[:] = \
            self.channel_vol_flow * density
        mass_flow_correction = \
            self.mass_flow_in / np.sum(self.channel_mass_flow)
        self.channel_mass_flow[:] *= mass_flow_correction


class UpdatedKohFlowCircuit(KohFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)

        self.urf = dict_flow_circuit.get('underrelaxation_factor', 0.5)
        self.density_channel = np.zeros(self.visc_channel.shape)

    def update_manifolds(self, update_fluid=True):
        if self.initialize:
            channel_mass_flow_in = np.ones(self.n_channels) \
                * self.mass_flow_in / self.n_channels
            channel_mass_flow_out = channel_mass_flow_in
        else:
            channel_mass_flow_in = self.channel_mass_flow
            channel_mass_flow_out = \
                np.array([channel.mass_flow_total[channel.id_out]
                          for channel in self.channels])
            channel_mass_flow_out *= self.n_subchannels

        if self.multi_component:
            mass_fraction = \
                np.array([channel.fluid.mass_fraction[:, channel.id_out]
                          for channel in self.channels]).transpose()
        else:
            mass_fraction = 1.0

        mass_source = channel_mass_flow_out * mass_fraction
        # mass_source = self.channel_mass_flow * mass_fraction
        channel_enthalpy_out = \
            np.asarray([ch.g_fluid[ch.id_out] * ch.temperature[ch.id_out]
                        for ch in self.channels]) * self.n_subchannels
        self.manifolds[1].update(mass_flow_in=0.0, mass_source=mass_source,
                                 update_mass=True, update_flow=True,
                                 update_heat=False, update_fluid=update_fluid,
                                 enthalpy_source=channel_enthalpy_out)

        # Inlet header update
        self.manifolds[0].p_out = \
            self.manifolds[1].pressure[-1] + self.dp_ref
        if self.multi_component:
            mass_fraction = self.manifolds[0].fluid.mass_fraction[:, :-1]
        else:
            mass_fraction = 1.0
        mass_source = -self.channel_mass_flow * mass_fraction
        self.manifolds[0].update(mass_flow_in=self.mass_flow_in,  # * 1.00000,
                                 mass_source=mass_source,
                                 update_mass=True, update_flow=True,
                                 update_heat=False, update_fluid=update_fluid)
        id_in = self.manifolds[0].id_in
        self.vol_flow_in = \
            self.mass_flow_in / self.manifolds[0].fluid.density[id_in]

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow

        if self.iteration == 0:
            self.visc_channel[:] = \
                np.array([np.average(channel.fluid.viscosity)
                          for channel in self.channels])
            self.density_channel[:] = \
                np.array([channel.fluid.density[channel.id_in]
                          for channel in self.channels])
            self.dp_channel = \
                np.array([channel.pressure[channel.id_in]
                          - channel.pressure[channel.id_out]
                          for channel in self.channels])

            self.k_perm[:] = self.channel_vol_flow / self.dp_channel \
                * self.visc_channel * self.l_by_a / self.n_subchannels

        try:
            self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a[-1] \
                * self.visc_channel[-1] / self.k_perm[-1] / self.n_subchannels
        except FloatingPointError:
            raise FloatingPointError('check if geomtries are adequate '
                                     'for flow conditions in {}'.
                                     format(self.name))
        self.update_manifolds(update_fluid=True)
        p_in = ip.interpolate_1d(self.manifolds[0].pressure)
        p_out = ip.interpolate_1d(self.manifolds[1].pressure)

        alpha = (p_in - p_out) / self.dp_ref
        self.alpha[:] = self.alpha * self.urf + alpha * (1.0 - self.urf)
        # self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a \
        #     * self.visc_channel[-1] / self.k_perm[-1] / self.n_subchannels
        # p_in += self.dp_ref \
        #     + self.manifolds[1].pressure[self.manifolds[1].id_out] \
        #     - self.manifolds[0].p_out
        # self.alpha[:] = (p_in - p_out) / self.dp_ref
        channel_vol_flow = (p_in - p_out) * self.k_perm \
            / self.l_by_a * self.n_subchannels / self.visc_channel
        self.channel_vol_flow[:] = self.channel_vol_flow * self.urf \
            + channel_vol_flow * (1.0 - self.urf)
        # self.channel_vol_flow[:] = (p_in - p_out) * self.k_perm \
        #     / self.l_by_a * self.n_subchannels / self.visc_channel

        self.channel_mass_flow[:] = self.channel_vol_flow * self.density_channel
        mass_flow_correction = \
            self.mass_flow_in / np.sum(self.channel_mass_flow)
        self.channel_mass_flow[:] *= mass_flow_correction


class WangFlowCircuit(ParallelFlowCircuit):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)

        # self.zeta = np.zeros(self.n_channels)
        self.xsi = 1.0
        self.H = self.manifolds[0].cross_area / self.manifolds[1].cross_area
        F_c = np.array([np.average(channel.cross_area)
                        for channel in self.channels])
        # print('F_c: ', F_c)
        sum_Fc = g_func.add_source(np.copy(F_c), F_c[1:], direction=-1)
        # print('sum_Fc: ', sum_Fc)
        self.M = sum_Fc / np.average(self.manifolds[0].cross_area)
        # print('self.M: ', self.M)
        # self.M = np.sum(F_c) / np.average(self.manifolds[0].cross_area)
        self.E = self.manifolds[0].length / self.manifolds[0].d_h
        self.D_star = self.manifolds[0].d_h / self.manifolds[1].d_h
        self.sqr_M = self.M ** 2.0
        self.sqr_H = self.H ** 2.0

        # print('M = ', self.M)
        # print('E = ', self.E)

    def update_channels(self):
        super().update_channels()
        if self.initialize:
            self.f_in = np.copy(self.manifolds[0].friction_factor)
            self.f_out = np.copy(self.manifolds[1].friction_factor)
        # if self.initialize:
        self.zeta = np.array([channel.zeta_bends * channel.n_bends
                                 for channel in self.channels]) \
            + np.array([np.sum(channel.friction_factor * channel.dx /
                               channel.d_h)
                        for channel in self.channels])
        self.zeta[:] += 1.0 + self.manifolds[0].zeta_other \
            + self.manifolds[1].zeta_other
        # self.zeta[:] = 10.0
        self.initialize = False

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            id_in = self.manifolds[0].id_in
            self.vol_flow_in = self.mass_flow_in \
                / self.manifolds[0].fluid.density[id_in]
        self.update_channels()

        mfd_in = self.manifolds[0]
        mfd_out = self.manifolds[1]

        k_in_0 = 0.6
        k_out_0 = 1.0
        b_in = 0.01
        b_out = 0.01
        W_0 = self.vol_flow_in / mfd_in.cross_area
        # print('W_0: ', W_0)
        # print('Re_0:', W_0 * mfd_in.fluid.density[0] * mfd_in.d_h /
        #       mfd_in.fluid.viscosity[0])
        # print('mfd_in.velocity[:-1]: ', mfd_in.velocity[:-1])
        # mfd_in.velocity[0] = W_0

        # print('zeta = ', self.zeta)

        f_in = mfd_in.friction_factor
        f_out = mfd_out.friction_factor
        # f_in = self.f_in
        # f_out = self.f_out
        #print('f_in: ', f_in)
        #print('f_out: ', f_out)
        # f_in[:] = 0.038
        # f_out[:] = 0.038
        k_in = k_in_0 + b_in * np.log(mfd_in.velocity[:-1] / W_0)
        k_out = k_out_0 + b_out * np.log(mfd_out.velocity[:-1] / W_0)

        Q = 2.0 / (3.0 * self.zeta) * (k_in - k_out * self.sqr_H) \
            * self.sqr_M
        R = - 0.25 * self.E * self.xsi / self.zeta \
            * (f_in + f_out * self.D_star * self.sqr_H) * self.sqr_M
        avg_R = np.average(R)
        avg_Q = np.average(Q)

        cube_Q = np.power(Q, 3.0)
        condition = np.square(R) + cube_Q
        avg_condition = np.square(avg_R) + np.power(avg_Q, 3.0)
        condition_0 = np.square(R[0]) + np.power(Q[0], 3.0)
        x = mfd_in.x / mfd_in.length
        one_third = 1.0 / 3.0
        # print('avg_condition: ', avg_condition)
        # print('condition: ', condition)
        w = 1.0
        for i in range(self.n_channels):
            # print('w_i: ', w)
            # k_in_i = k_in_0 + b_in * np.log(w)
            # k_out_i = k_out_0 + b_out * np.log(w * self.H)
            # Q_i = 2.0 / (3.0 * self.zeta[i]) * (
            #             k_in_i - k_out_i * self.sqr_H) * self.sqr_M
            # R_i = - 0.25 * self.E * self.xsi / self.zeta[i] \
            #     * (f_in[i] + f_out[i] * self.D_star * self.sqr_H) * self.sqr_M
            # cube_Q_i = np.power(Q_i, 3.0)
            # square_R_i = np.square(R_i)
            # condition_i = square_R_i + cube_Q_i
            # print('cube_Q_i: ', cube_Q_i)
            # print('square_R_i: ', square_R_i)
            condition_i = condition[i]
            R_i = R[i]
            Q_i = Q[i]
            cube_Q_i = cube_Q[i]
            # print('condition: ', condition_i)

            if condition_i < 0.0:
                theta = np.arccos(R_i/np.sqrt(-cube_Q_i))
                sqrt_Q = np.sqrt(-Q_i)
                r_1 = 2.0 * sqrt_Q * np.cos(theta * one_third)
                r_2 = 2.0 * sqrt_Q * np.cos((theta + 2.0*np.pi) * one_third)
                w = (np.exp(r_1 + r_2 * x[i+1]) - np.exp(r_2 + r_1 * x[i+1])) \
                    / (np.exp(r_1) - np.exp(r_2))
                # print('i :', i, ', condition < 0,  w: ', w)
            elif condition_i == 0.0:
                r = - 0.5 * np.power(R_i, one_third)
                w = (1.0 - x[i+1]) * np.exp(r*x[i+1])
                # print('i :', i, ', condition == 0,  w: ', w)
            else:
                sqrt_condition = np.sqrt(condition_i)
                term_1 = np.cbrt(R_i + sqrt_condition)
                term_2 = np.cbrt(R_i - sqrt_condition)
                B = term_1 + term_2
                J = term_1 - term_2
                sqrt3_J_by_2 = np.sqrt(3.0) * J * 0.5
                w = np.exp(-B * x[i+1] * 0.5) \
                    * np.sin(sqrt3_J_by_2 * (1.0 - x[i+1])) \
                    / np.sin(sqrt3_J_by_2)
                # print('i :', i, ', condition > 0,  w: ', w)
            W = w * W_0
            mfd_in.velocity[i+1] = W
            mfd_out.velocity[i+1] = W * self.H \
                * mfd_in.fluid.density[i+1] / mfd_out.fluid.density[i+1]

        # print('condition: ', condition)
        mass_flow_in = \
            mfd_in.velocity * mfd_in.fluid.density * mfd_in.cross_area
        self.channel_mass_flow[:] = mass_flow_in[:-1] - mass_flow_in[1:]
        self.channel_vol_flow[:] = \
            self.channel_mass_flow / ip.interpolate_1d(mfd_in.fluid.density)
        # print('distribution: ', self.channel_vol_flow/(np.sum(
        #     self.channel_vol_flow)/self.n_channels))


def factory(dict_circuit, dict_in_manifold, dict_out_manifold,
            channels, channel_multiplier=1.0):
    if not isinstance(channels, (list, tuple)):
        raise TypeError('argument channels must be a list of type Channel')
    if not isinstance(channels[0], chl.Channel):
        raise TypeError('argument channels must be a list of type Channel')

    n_channels = len(channels)

    if hasattr(channels[0].fluid, 'dict') \
            and isinstance(channels[0].fluid, fluids.CanteraGasMixture):
        fluid_dict = channels[0].fluid.dict
        fluid_dict['nodes'] = n_channels + 1
        in_manifold_fluid = fluids.factory(fluid_dict)
        out_manifold_fluid = fluids.factory(fluid_dict)
    else:
        in_manifold_fluid = channels[0].fluid.copy()
        in_manifold_fluid.rescale(n_channels + 1)
        out_manifold_fluid = in_manifold_fluid.copy()

    manifolds = [chl.Channel(dict_in_manifold, in_manifold_fluid),
                 chl.Channel(dict_out_manifold, out_manifold_fluid)]

    return ParallelFlowCircuit(dict_circuit, manifolds, channels,
                               n_subchannels=channel_multiplier)
