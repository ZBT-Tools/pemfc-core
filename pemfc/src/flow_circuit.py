# general imports
import numpy as np
from abc import ABC, abstractmethod
import math
import copy

# local module imports
from . import interpolation as ip, global_functions as g_func, \
    channel as chl, output_object as oo, flow_resistance as fr
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
        elif circuit_type == 'VariableResistance':
            return super(ParallelFlowCircuit, cls).\
                __new__(VariableResistanceFlowCircuit)
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
        # self.update_channels(update_fluid=True)
        self.normalized_flow_distribution = \
            np.zeros(self.channel_vol_flow.shape)
        self.iteration = 0
        self.error = 1e8
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
                if self.mass_flow_in > 0.0:
                    self.channel_mass_flow[:] *= \
                        inlet_mass_flow / self.mass_flow_in
                    self.channel_vol_flow[:] *= \
                        inlet_mass_flow / self.mass_flow_in
            # set new mass and volume flows
            self.mass_flow_in = inlet_mass_flow

        if self.initialize:
            # self.update_channels(update_fluid=True)
            self.channel_vol_flow_old[:] = 1e8
        self.update_channels(update_fluid=True)
        # channel_vol_flow_old = np.zeros(self.channel_vol_flow.shape)
        if calc_distribution is None:
            calc_distribution = self.calc_distribution
        if calc_distribution and self.n_channels > 1:
            for i in range(self.max_iter):
                # print(self.name + ' Iteration # ', str(i+1))
                self.iteration = i
                self.single_loop(update_channels=True)
                if i == 0:
                    self.initialize = False

                self.error = g_func.calc_rrmse(
                    self.channel_vol_flow, self.channel_vol_flow_old)

                # print(channel_vol_flow_old)
                # print(self.channel_vol_flow)
                self.channel_vol_flow_old[:] = self.channel_vol_flow
                # print(error)
                if self.error < self.tolerance and i >= self.min_iter:
                    # print('number of iterations n = {} '
                    #       'with error = {} in update() of {} '
                    #       .format(self.iteration, self.error, self.name))
                    break

                if i == (self.max_iter - 1):
                    print('maximum number of iterations n = {} '
                          'with error = {} in update() of {} reached'
                          .format(self.max_iter, self.error, self.name))

        # else:
        #     self.initialize = False
            # Final channel updates within flow circuit iteration
            self.update_channels(update_fluid=True)

        try:
            self.normalized_flow_distribution[:] = \
                self.channel_mass_flow / np.average(self.channel_mass_flow)
        except FloatingPointError:
            self.normalized_flow_distribution[:] = 0.0

        # After first or any later update call, initialization is over
        self.initialize = False

    @abstractmethod
    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        pass

    def update_channels(self, update_fluid=True):
        if self.initialize:
            channel_mass_flow_in = np.ones(self.n_channels) \
                * self.mass_flow_in / self.n_channels
        else:
            channel_mass_flow_in = self.channel_mass_flow

        # Channel update
        for i, channel in enumerate(self.channels):
            # channel.p_out = p_channel_out[i]
            # channel.temperature[channel.id_in] = self.manifolds[0].temp_ele[i]
            channel.update(mass_flow_in=
                           channel_mass_flow_in[i] / self.n_subchannels,
                           update_mass=True, update_flow=True,
                           update_heat=False, update_fluid=update_fluid)

        self.vol_flow_in = self.calc_volume_flow()

    def update_manifolds(self, reference_pressure, update_fluid=True):
        if self.initialize:
            channel_mass_flow_in = np.ones(self.n_channels) \
                * self.mass_flow_in / self.n_channels
            channel_mass_flow_out = channel_mass_flow_in
        else:
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

        # Reference to inlet manifold
        mfd_in = self.manifolds[0]
        # Reference to outlet manifold
        mfd_out = self.manifolds[1]

        # Update outlet manifold
        mass_source = channel_mass_flow_out * mass_fraction
        # mass_source = self.channel_mass_flow * mass_fraction
        channel_enthalpy_out = \
            np.asarray([ch.g_fluid[ch.id_out] * ch.temperature[ch.id_out]
                        for ch in self.channels]) * self.n_subchannels
        mfd_out.update(mass_flow_in=0.0, mass_source=mass_source,
                       update_mass=True, update_flow=True,
                       update_heat=False, update_fluid=update_fluid,
                       enthalpy_source=channel_enthalpy_out)

        id_in = mfd_out.id_in - 1
        mfd_in.p_out = mfd_out.pressure[id_in] + reference_pressure

        # Update inlet manifold
        if self.multi_component:
            mass_fraction = self.manifolds[0].fluid.mass_fraction[:, :-1]
        else:
            mass_fraction = 1.0
        mass_source = -self.channel_mass_flow * mass_fraction
        mfd_in.update(mass_flow_in=self.mass_flow_in,  # * 1.00000,
                      mass_source=mass_source,
                      update_mass=True, update_flow=True,
                      update_heat=False, update_fluid=update_fluid)
        self.vol_flow_in = \
            self.mass_flow_in / mfd_in.fluid.density[mfd_in.id_in]

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

    def single_loop(self, inlet_mass_flow=None, update_channels=False):
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

        if np.any(self.channel_vol_flow == 0.0):
            raise ValueError('zero flow rates detected, '
                             'check boundary conditions')
        if self.initialize:
            self.k_perm[:] = self.channel_vol_flow / self.dp_channel \
                * self.visc_channel * self.l_by_a
        self.dp_ref = np.maximum(self.dp_channel[-1], 1e-3)

        try:
            self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a \
                * self.visc_channel[-1] / self.k_perm[-1] / self.n_subchannels
        except FloatingPointError:
            raise FloatingPointError('check if geomtries are adequate '
                                     'for flow conditions in {}'.
                                     format(self.name))
        self.update_manifolds(self.dp_ref)

        p_in = ip.interpolate_1d(self.manifolds[0].pressure)
        p_out = ip.interpolate_1d(self.manifolds[1].pressure)

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
        if update_channels:
            self.update_channels()
        try:
            self.dp_ref = self.vol_flow_in / np.sum(self.alpha) * self.l_by_a[-1] \
                * self.visc_channel[-1] / self.k_perm[-1] / self.n_subchannels
        except FloatingPointError:
            raise FloatingPointError('Check if geometries are adequate '
                                     'for flow conditions in {}'.
                                     format(self.name))
        self.update_manifolds(self.dp_ref, update_fluid=True)
        p_in = ip.interpolate_1d(self.manifolds[0].pressure)
        p_out = ip.interpolate_1d(self.manifolds[1].pressure)

        alpha = (p_in - p_out) / self.dp_ref
        self.alpha[:] = self.alpha * self.urf + alpha * (1.0 - self.urf)

        channel_vol_flow = (p_in - p_out) * self.k_perm \
            / self.l_by_a * self.n_subchannels / self.visc_channel
        self.channel_vol_flow[:] = self.channel_vol_flow * self.urf \
            + channel_vol_flow * (1.0 - self.urf)

        self.channel_mass_flow[:] = self.channel_vol_flow * self.density_channel
        mass_flow_correction = \
            self.mass_flow_in / np.sum(self.channel_mass_flow)
        self.channel_mass_flow[:] *= mass_flow_correction


class WangFlowCircuit(ParallelFlowCircuit):
    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):

        # manifold to channel resistance
        flow_res_dict = \
            {'type': 'Junction', 'coefficients': [1.02, 0.53, 33.67]}
        # flow_res_dict = {'type': 'Junction', 'coefficients': [0.0]}
        # flow_res_dict = {'type': 'Junction', 'coefficients': [0.02, -10.67]}
        self.mfd_to_chl_res = fr.FlowResistance(manifolds[0], flow_res_dict)
        # manifold to channel resistance
        # flow_res_dict['coefficients'] = [1.02, 0.53, 33.67]
        flow_res_dict['coefficients'] = [4.5]
        # flow_res_dict = {'type': 'Junction', 'coefficients': [-0.02, 10.67]}
        # flow_res_dict['coefficients'] = [1.0, 3.0]
        self.chl_to_mfd_res = fr.FlowResistance(manifolds[1], flow_res_dict)

        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)

        self.zeta = np.zeros(self.n_channels)
        self.xsi = 1.0
        self.H = self.manifolds[0].cross_area / self.manifolds[1].cross_area
        self.F_c = np.array([np.average(channel.cross_area)
                            for channel in self.channels])
        # print('F_c: ', F_c)
        # sum_Fc = g_func.add_source(np.copy(F_c), F_c[1:], direction=-1)
        # print('sum_Fc: ', sum_Fc)
        self.M = np.sum(self.F_c) / np.average(self.manifolds[0].cross_area)
        # print('self.M: ', self.M)
        # self.M = np.sum(F_c) / np.average(self.manifolds[0].cross_area)
        self.E = self.manifolds[0].length / self.manifolds[0].d_h
        self.D_star = self.manifolds[0].d_h / self.manifolds[1].d_h
        self.sqr_M = self.M ** 2.0
        self.sqr_H = self.H ** 2.0

        # Assign inlet and outlet manifold channel friction factors
        # to flow circuit model
        self.f_in = np.zeros(self.n_channels)
        self.f_out = np.zeros(self.n_channels)

    def update_manifold_friction_factors(self):
        for zeta in self.manifolds[0].zetas:
            if isinstance(zeta, fr.WallFrictionFlowResistance):
                self.f_in[:] = zeta.value[:-1]
                # self.f_in[:] = 0.0
                # self.f_in[:] = 0.038
        for zeta in self.manifolds[1].zetas:
            if isinstance(zeta, fr.WallFrictionFlowResistance):
                self.f_out[:] = zeta.value[:-1]
                # self.f_out[:] = 0.0
                # self.f_out[:] = 0.038

    def update_channels(self, **kwargs):
        super().update_channels(**kwargs)
        # if self.initialize:
        #     self.f_in = np.copy(self.manifolds[0].friction_factor)
        #     self.f_out = np.copy(self.manifolds[1].friction_factor)
        # if self.initialize:
        self.zeta = np.array([np.sum(channel.calculate_flow_resistance_sum())
                              for channel in self.channels])
        # Update t-junction branching flow resistances
        self.mfd_to_chl_res.update()
        self.chl_to_mfd_res.update()
        self.zeta[:] += \
            1.0 + self.mfd_to_chl_res.value + self.chl_to_mfd_res.value
        # # Manifold-to-Channel and Channel-to-Manifold resistance coefficients
        # self.zeta[:] += self.manifolds[0].zeta_other + self.manifolds[1].zeta_other
        # self.zeta[:] = 10.0
        # self.initialize = False

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
            id_in = self.manifolds[0].id_in
            self.vol_flow_in = self.mass_flow_in \
                / self.manifolds[0].fluid.density[id_in]
        # self.update_channels()

        self.update_manifold_friction_factors()
        mfd_in = self.manifolds[0]
        mfd_out = self.manifolds[1]

        k_in_0 = 0.6
        k_out_0 = 1.0
        b_in = 0.1
        b_out = 0.1
        W_0 = self.vol_flow_in / mfd_in.cross_area
        # print('W_0: ', W_0)
        # print('Re_0:', W_0 * mfd_in.fluid.density[0] * mfd_in.d_h /
        #       mfd_in.fluid.viscosity[0])
        Re_0 = W_0 * mfd_in.fluid.density[0] * mfd_in.d_h / mfd_in.fluid.viscosity[0]

        # print('mfd_in.velocity[:-1]: ', mfd_in.velocity[:-1])
        # mfd_in.velocity[0] = W_0

        # print('zeta = ', self.zeta)

        # f_in = mfd_in.friction_factor
        # f_out = mfd_out.friction_factor
        # f_in = self.f_in
        # f_out = self.f_out
        #print('f_in: ', f_in)
        #print('f_out: ', f_out)
        # f_in[:] = 0.038
        # f_out[:] = 0.038
        m = mfd_in.velocity[:-1]
        k_in = k_in_0 + b_in * g_func.np_log(mfd_in.velocity[:-1] / W_0)

        k_out = k_out_0 + b_out * g_func.np_log(mfd_out.velocity[:-1] / W_0)

        zeta = self.zeta
        # zeta = 10.0

        Q = g_func.np_div(2.0, (3.0 * zeta)) \
            * (k_in - k_out * self.sqr_H) * self.sqr_M

        R = - g_func.np_div(0.25 * self.E * self.xsi, zeta) \
            * (self.f_in + self.f_out * self.D_star * self.sqr_H) * self.sqr_M

        compare_paper = {'Re_0': Re_0, 'M': self.M, 'E': self.E,
                         'zeta': zeta}

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
        print('condition: ', condition)
        alpha = np.zeros(self.n_channels)
        if condition[0] < 0.0:
            switch_case = 1
        elif condition[0] == 0.0:
            switch_case = 2
        else:
            switch_case = 3
        for i in range(self.n_channels):
            condition_i = condition[i]
            R_i = R[i]
            Q_i = Q[i]
            cube_Q_i = cube_Q[i]
            # print('condition: ', condition_i)

            if condition_i < 0.0:
                theta = np.arccos(R_i/np.sqrt(-cube_Q_i))
                sqrt_Q = np.sqrt(-Q_i)
                r_1 = 2.0 * sqrt_Q * np.cos(theta * one_third)
                r_2 = 2.0 * sqrt_Q * np.cos((theta + 2.0 * np.pi) * one_third)
                w = (np.exp(r_1 + r_2 * x[i+1]) - np.exp(r_2 + r_1 * x[i+1])) \
                    / (np.exp(r_1) - np.exp(r_2))
                alpha[i] = - (r_2 * np.exp(r_1 + r_2 * x[i+1])
                              - r_1 * np.exp(r_2 + r_1 * x[i+1])) \
                    / (np.exp(r_1) - np.exp(r_2))
                # print('i :', i, ', condition < 0,  w: ', w)
            elif condition_i == 0.0:
                r = - 0.5 * np.power(R_i, one_third)
                w = (1.0 - x[i+1]) * np.exp(r*x[i+1])
                alpha[i] = - (r - 1.0 - r * x[i+1]) * np.exp(r * x[i+1])
                # print('i :', i, ', condition == 0,  w: ', w)
            else:
                sqrt_condition = np.sqrt(condition_i)
                term_1 = np.cbrt(R_i + sqrt_condition)
                term_2 = np.cbrt(R_i - sqrt_condition)
                B = term_1 + term_2
                J = term_1 - term_2
                sqrt3_J = math.sqrt(3.0) * J
                sqrt3_J_by_2 = math.sqrt(3.0) * J * 0.5
                w = np.exp(-B * x[i+1] * 0.5) \
                    * np.sin(sqrt3_J_by_2 * (1.0 - x[i+1])) \
                    / np.sin(sqrt3_J_by_2)

                alpha[i] = 0.5 * np.exp(-B * x[i+1] * 0.5) \
                    * ((B * np.sin(sqrt3_J_by_2 * (1.0 - x[i+1]))
                        + sqrt3_J * np.cos(sqrt3_J_by_2 * (1.0 - x[i+1])))
                       / np.sin(sqrt3_J_by_2))
                # print('i :', i, ', condition > 0,  w: ', w)
            # W = w * W_0
            # mfd_in.velocity[i+1] = W
            # mfd_out.velocity[i+1] = W * self.H \
            #     * mfd_in.fluid.density[i+1] / mfd_out.fluid.density[i+1]

        # print('condition: ', condition)
        mass_flow_in = \
            mfd_in.velocity * mfd_in.fluid.density * mfd_in.cross_area
        # self.channel_mass_flow[:] = mass_flow_in[:-1] - mass_flow_in[1:]

        self.channel_mass_flow[:] = alpha * self.mass_flow_in / self.n_channels
        self.channel_vol_flow[:] = \
            self.channel_mass_flow / ip.interpolate_1d(mfd_in.fluid.density)
        # print('distribution: ', self.channel_vol_flow/(np.sum(
        #     self.channel_vol_flow)/self.n_channels))


class VariableResistanceFlowCircuit(ParallelFlowCircuit):

    def __init__(self, dict_flow_circuit, manifolds, channels,
                 n_subchannels=1.0):
        super().__init__(dict_flow_circuit, manifolds, channels,
                         n_subchannels)
        self.urf = dict_flow_circuit.get('underrelaxation_factor', 0.3)

        self.turning_resistance = []
        self.manifold_wall_friction = []
        for mfd in manifolds:
            for zeta in mfd.zetas:
                if isinstance(zeta, fr.JunctionFlowResistance):
                    zeta_dict = copy.deepcopy(zeta.dict)
                    zeta_dict['type'] = zeta_dict['type'].replace('Main', 'Branch')
                    self.turning_resistance.append(
                        fr.FlowResistance(mfd, zeta_dict))
                if isinstance(zeta, fr.WallFrictionFlowResistance):
                    self.manifold_wall_friction.append(zeta)

    def calc_dp_ref(self, pressure_drop):
        """
        Calculates reference pressure drop between inlet and outlet manifold,
        which is defined at the last junction. Given the pressure drop in the
        channel in addition to turning losses, this function adds the
        dynamic pressure changes, which should only be non-zero for varying fluid
        properties and/or different cross-sectional areas between manifolds
        :param pressure_drop:
        :return:
        """
        # Reference to inlet manifold
        mfd_in = self.manifolds[0]
        # Reference to outlet manifold
        mfd_out = self.manifolds[1]
        ref_chl = self.channels[-1]
        id_in = mfd_out.id_in - 1
        id_out = mfd_in.id_out - 1
        id_out_chl = ref_chl.id_out
        id_in_chl = ref_chl.id_in
        dp_dyn = \
            (mfd_out.fluid.density[id_in] * mfd_out.velocity[id_in] ** 2.0
             - ref_chl.fluid.density[id_out_chl] * ref_chl.velocity[id_out_chl] ** 2.0
             + ref_chl.fluid.density[id_in_chl] * ref_chl.velocity[id_in_chl] ** 2.0
             - mfd_in.fluid.density[id_out] * mfd_in.velocity[id_out] ** 2.0) * 0.5
        return pressure_drop + dp_dyn

    def single_loop(self, inlet_mass_flow=None, update_channels=True):
        """
        Update the flow circuit
        """
        if inlet_mass_flow is not None:
            self.mass_flow_in = inlet_mass_flow
        if update_channels:
            self.update_channels(update_fluid=False)
        self.channel_vol_flow[:] = \
            np.array([np.average(channel.vol_flow)
                      for channel in self.channels])
        if np.any(self.channel_vol_flow == 0.0):
            raise ValueError('zero flow rates in channels detected, '
                             'check boundary conditions')
        # velocity = np.array([np.average(channel.velocity)
        #                      for channel in self.channels])

        # Update T-junction branching flow resistances
        for item in self.turning_resistance:
            assert isinstance(item, fr.JunctionFlowResistance)

        dp_turning = np.zeros((len(self.manifolds), self.n_channels))
        dp_wall = np.zeros((len(self.manifolds), self.n_channels))

        ref_chl = self.channels[-1]

        # Sub iteration loop for updating manifold pressure line including
        # last reference cell pressure drop
        error = 1e8
        while error > 1e-3:
            pressure_old = \
                np.asarray([mfd.pressure for mfd in self.manifolds]).flatten()

            for i, zeta_turn in enumerate(self.turning_resistance):
                # Update T-junction main-to-branch flow resistances
                zeta_turn.update()
                # and calculate the additional pressure drops based on current
                # flow distribution due to branching
                dp_turning[i, :] = zeta_turn.calc_pressure_drop()
                # and calculate pressure drop due to wall friction for the
                # distance from combined velocity reference point to the
                # branching point, which is half the nodal distance
                if self.manifold_wall_friction:
                    zeta_wall = \
                        self.manifold_wall_friction[i].value[zeta_turn.combined_ids]
                    dp_wall[i, :] = 0.5 * zeta_turn.density * zeta_wall \
                        * zeta_turn.velocity ** 2.0
                    # Only pressure drops for interior nodes need to be split
                    # in half, because outer nodes are already dx/2 long
                    dp_wall[i, 1:-2] *= 0.5

            dp_ref_chl = \
                ref_chl.pressure[ref_chl.id_in] - ref_chl.pressure[ref_chl.id_out] \
                + dp_turning[0, -1] + dp_turning[1, -1] \
                + dp_wall[0, -1] + dp_wall[1, -1]
            dp_ref = self.calc_dp_ref(dp_ref_chl)
            # Update manifolds with current reference pressure drop in last channel
            self.update_manifolds(dp_ref)

            pressure_new = \
                np.asarray([mfd.pressure for mfd in self.manifolds]).flatten()
            error = g_func.calc_rrmse(pressure_old, pressure_new)
            if not self.initialize:
                break

        if self.turning_resistance:
            # Calculate inlet pressures for channels according to inlet manifold
            # T-junctions, different equations for potentially combining or dividing
            # flow from manifold
            chl_in_density = \
                np.array([chl.fluid.density[chl.id_in] for chl in self.channels])
            chl_out_density = \
                np.array([chl.fluid.density[chl.id_out] for chl in self.channels])
            chl_in_velocity = \
                np.array([chl.velocity[chl.id_in] for chl in self.channels])
            chl_out_velocity = \
                np.array([chl.velocity[chl.id_out] for chl in self.channels])

            # Reference point is defined within junction flow resistance object as
            # the point with combined flow rates of main run
            p_branch_mfd = np.zeros(dp_turning.shape)
            for i, zeta_turn in enumerate(self.turning_resistance):

                id_div = np.nonzero(zeta_turn.mass_source_sign <= 0)
                id_com = np.nonzero(zeta_turn.mass_source_sign > 0)

                # Branch pressure for dividing tees from manifold
                p_branch_mfd[i, id_div] = zeta_turn.pressure[id_div] \
                    + 0.5 * zeta_turn.density[id_div] * zeta_turn.velocity[id_div] ** 2.0 \
                    - 0.5 * chl_in_density[id_div] * chl_in_velocity[id_div] ** 2.0 \
                    - dp_turning[i, id_div] - dp_wall[i, id_div]

                    # Branch pressure for combining tees with manifold
                p_branch_mfd[i, id_com] = zeta_turn.pressure[id_com] \
                    + 0.5 * zeta_turn.density[id_com] * zeta_turn.velocity[id_com] ** 2.0 \
                    - 0.5 * chl_out_density[id_com] * chl_out_velocity[id_com] ** 2.0 \
                    + dp_turning[i, id_com] + dp_wall[i, id_com]

                    # Assign outlet pressure to channels
            p_branch_out = np.min(p_branch_mfd, axis=0)
            for i, chl in enumerate(self.channels):
                chl.p_out = p_branch_out[i]

            dp_branches = p_branch_mfd[0] - p_branch_mfd[1]
        else:
            p_in = ip.interpolate_1d(self.manifolds[0].pressure)
            p_out = ip.interpolate_1d(self.manifolds[1].pressure)
            for i, chl in enumerate(self.channels):
                chl.p_out = p_out[i]
            dp_branches = p_in - p_out + np.sum(dp_turning, axis=0)

        dp_channel = \
            np.array([channel.pressure[channel.id_in]
                      - channel.pressure[channel.id_out]
                      for channel in self.channels])
        flow_correction = dp_branches / dp_channel

        channel_vol_flow_new = self.channel_vol_flow * flow_correction
        urf = self.urf
        # urf = 0.3
        channel_vol_flow = \
            self.channel_vol_flow * urf + channel_vol_flow_new * (1.0 - urf)
        density = np.array([channel.fluid.density[channel.id_in]
                            for channel in self.channels])
        channel_mass_flow = channel_vol_flow * density
        # Correct mass balance
        channel_mass_flow *= self.mass_flow_in / np.sum(channel_mass_flow)
        channel_vol_flow = channel_mass_flow / density
        # Assign to member variables
        self.channel_mass_flow[:] = channel_mass_flow
        self.channel_vol_flow[:] = channel_vol_flow


def add_flow_resistance(channel_dict, resistance_dict):
    """
    Add additional flow resistance input to channel dictionary, e.g. for T-Junctions.
    :param channel_dict: channel input dictionary for channel which holds resistance objects
    :param resistance_dict: input dictionary for flow resistance object
    :return: channel dictionary with additional flow resistance in 'flow_resistances' entry
    """
    if 'flow_resistances' in channel_dict:
        if isinstance(channel_dict['flow_resistances'], list):
            channel_dict['flow_resistances'].append(resistance_dict)
        else:
            raise TypeError('Entry "flow_resistances" in channel '
                            'dictionary must be a list')
    else:
        channel_dict['flow_resistances'] = [resistance_dict]
    return channel_dict


def create(dict_circuit, dict_in_manifold, dict_out_manifold,
           channels, channel_multiplier=1.0):
    """
    Factory method to create a concrete instantiation of the class
    ParallelFlowCircuit
    :param dict_circuit: flow circuit dictionary with inputs for the overall
    flow circuit properties
    :param dict_in_manifold: channel input dictionary for the inlet manifold
    :param dict_out_manifold: channel input dictionary for the outlet manifold
    :param channels: list of channel objects referencing the parallel channels
    :param channel_multiplier: mass flow multiplier for each representative
    channel such that header mass flows correspond to total channel mass flows
    :return: ParallelFlowCircuit object
    """
    """
    :param dict_circuit: 
    :param dict_in_manifold: 
    :param dict_out_manifold: 
    :param channels: 
    :param channel_multiplier: 
    :return: 
    """
    if not isinstance(channels, (list, tuple)):
        raise TypeError('argument channels must be a list of type Channel')
    if not isinstance(channels[0], chl.Channel):
        raise TypeError('argument channels must be a list of type Channel')

    n_channels = len(channels)

    if hasattr(channels[0].fluid, 'dict') \
            and isinstance(channels[0].fluid, fluids.CanteraGasMixture):
        fluid_dict = channels[0].fluid.dict
        fluid_dict['nodes'] = n_channels + 1
        in_manifold_fluid = fluids.create(fluid_dict)
        out_manifold_fluid = fluids.create(fluid_dict)
    else:
        in_manifold_fluid = channels[0].fluid.copy()
        in_manifold_fluid.rescale(n_channels + 1)
        out_manifold_fluid = in_manifold_fluid.copy()

    manifolds = [chl.Channel(dict_in_manifold, in_manifold_fluid),
                 chl.Channel(dict_out_manifold, out_manifold_fluid)]

    return ParallelFlowCircuit(dict_circuit, manifolds, channels,
                               n_subchannels=channel_multiplier)
