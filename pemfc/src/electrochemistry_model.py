# general imports
from abc import ABC, abstractmethod
import numpy as np
from scipy import optimize

# local module imports
from . import interpolation as ip, constants, channel as chl


class ElectrochemistryModel(ABC):

    def __init__(self, input_dict: dict, nodes: int):
        
        """voltage loss parameter, (Kulikovsky, 2013)"""

        self.n_ele = nodes - 1
        self.faraday = constants.FARADAY
        # index of fuel component in gas mixture
        self.id_fuel = input_dict['fuel_index']
        # charge number of reaction
        self.n_charge = input_dict['charge_number']
        # exchange current density
        vol_ex_cd = input_dict['vol_ex_cd']
        # proton conductivity of the catalyst layer
        self.prot_con_cl = input_dict['prot_con_cl']
        # diffusion coefficient of the reactant in the catalyst layer
        self.diff_coeff_cl = input_dict['diff_coeff_cl']
        # diffusion coefficient of the reactant in the gas diffusion layer
        self.diff_coeff_gdl = input_dict['diff_coeff_gdl']
        # tafel slope of the electrode
        self.tafel_slope = input_dict['tafel_slope']
        self.th_gdl = input_dict['thickness_gdl']
        self.th_cl = input_dict['thickness_cl']
        # could use a better name see (Kulikovsky, 2013) not sure if 2-D
        # exchange current density
        self.i_sigma = np.sqrt(2. * vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        # index of the first element with negative cell voltage
        self.index_cat = nodes - 1
        # characteristic current density, see (Kulikovsky, 2013)
        self.i_star = self.prot_con_cl * self.tafel_slope / self.th_cl
        # concentration at channel inlet
        self.conc_in = None
        # limiting current density due to diffusion through the gdl
        # at channel inlet (calculated when inlet concentration is known)
        self.i_lim_star = None
        # numerical parameter for tangent line extension at limiting current
        self.conc_eps = input_dict['c_eps']
        self.delta_i = input_dict['delta_i']
        # critical local current density where Kulikovsky model transitions
        # into linear tangent line near limiting current
        self.i_crit = np.zeros(self.n_ele)

        self.calc_act_loss = input_dict['calc_act_loss']
        self.calc_cl_diff_loss = input_dict['calc_cl_diff_loss']
        self.calc_gdl_diff_loss = input_dict['calc_gdl_diff_loss']

        # cell voltage loss
        self.v_loss = np.zeros(self.n_ele)
        self.updated_v_loss = False
        self.corrected_current_density = None

    def update(self, current_density: np.ndarray, channel: chl.Channel,
               current_control: bool = True):
        """
        This function coordinates the program sequence
        """
        # self.calc_temp_fluid_ele()
        # mole_flow_in, mole_source = self.calc_mass_balance(current_density)
        if np.any(current_density < 0.0):
            raise ValueError('current density became smaller 0')
        if not current_control and self.updated_v_loss:
            self.corrected_current_density = \
                self.calc_current_density(current_density, channel, self.v_loss)
        if current_control or self.corrected_current_density is None:
            corrected_current_density = current_density
        else:
            corrected_current_density = self.corrected_current_density
        self.update_voltage_loss(corrected_current_density, channel)

    def update_voltage_loss(self, current_density: np.ndarray,
                            channel: chl.Channel):
        eta = self.calc_electrode_loss(current_density, channel)
        self.v_loss[:] = eta # \
                         # + self.calc_plate_loss(current_density)
        self.updated_v_loss = True

    def calc_activation_loss(self, current_density: np.ndarray,
                             conc: np.ndarray):
        """
        Calculates the activation voltage loss,
        according to (Kulikovsky, 2013).
        """
        np.seterr(divide='ignore')
        try:
            v_loss_act = \
                np.where(np.logical_and(current_density > constants.SMALL,
                                        conc > constants.SMALL),
                         self.tafel_slope
                         * np.arcsinh((current_density / self.i_sigma) ** 2.
                                      / (2. * conc
                                         * (1. - np.exp(-current_density /
                                                        (2. * self.i_star))))),
                         0.0)
            np.seterr(divide='raise')
        except FloatingPointError:
            raise
        return v_loss_act

    def calc_transport_loss_catalyst_layer(self, current_density: np.ndarray,
                                           var: float, conc: np.ndarray):
        """
        Calculates the diffusion voltage loss in the catalyst layer
        according to (Kulikovsky, 2013).
        """
        try:
            i_hat = current_density / self.i_star
            short_save = np.sqrt(2. * i_hat)
            beta = \
                short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save)) \
                + np.pi * i_hat / (2. + i_hat)
        except FloatingPointError:
            test = np.any(current_density < 0.0)
            raise
        try:
            v_loss_cl_diff = \
                ((self.prot_con_cl * self.tafel_slope ** 2.)
                 / (4. * self.faraday * self.diff_coeff_cl * conc)
                 * (current_density / self.i_star
                    - np.log10(1. + np.square(current_density) /
                               (self.i_star ** 2. * beta ** 2.)))) / var
        except FloatingPointError:
            raise
        return v_loss_cl_diff

    def calc_transport_loss_diffusion_layer(self, var: float):
        """
        Calculates the diffusion voltage loss in the gas diffusion layer
        according to (Kulikovsky, 2013).
        """
        try:
            v_loss_gdl_diff = -self.tafel_slope * np.log10(var)
        except FloatingPointError:
            raise
        # nan_list = np.isnan(self.v_loss_gdl_diff)
        # if nan_list.any():
        #     v_loss_gdl_diff[np.argwhere(nan_list)[0, 0]:] = 1.e50
        return v_loss_gdl_diff

    def calc_electrode_loss(self, current_density: np.ndarray,
                            channel: chl.Channel):
        if hasattr(channel.fluid, 'gas'):
            conc = channel.fluid.gas.concentration[self.id_fuel]
        else:
            conc = channel.fluid.concentration[self.id_fuel]
        conc_ele = ip.interpolate_1d(conc)
        conc_ref = conc[channel.id_in]
        conc_star = conc_ele / conc_ref
        # if self.channel.flow_direction == 1:
        #     conc_in = conc[:-1]
        # else:
        #     conc_in = conc[1:]
        conc_in = conc[channel.id_in]
        if conc_in != self.conc_in:
            self.i_lim_star = self.n_charge * self.faraday * conc_in \
                              * self.diff_coeff_gdl / self.th_gdl
            self.conc_in = conc_in
        self.i_crit[:] = self.i_lim_star * (conc_ele - self.conc_eps) / conc_ref
        id_lin = np.argwhere(current_density >= self.i_crit)[:, 0]
        id_reg = np.argwhere(current_density < self.i_crit)[:, 0]
        if len(id_lin) > 0:
            i_crit = self.i_crit[id_lin]
            conc_crit = conc_ele[id_lin]
            conc_crit = \
                np.vstack((conc_crit, conc_crit, conc_crit))
            i_crit = np.vstack(
                (i_crit - self.delta_i, i_crit, i_crit + self.delta_i))
            conc_crit = conc_crit.transpose()
            i_crit = i_crit.transpose()
            # if np.any(i_crit < 0.0):
            #     raise ValueError
            eta_crit = \
                self.calc_electrode_loss_kulikovsky(i_crit, conc_crit, conc_ref,
                                                    update_members=False)

            grad_eta = np.gradient(eta_crit, self.delta_i, axis=-1)[:, 1]
            b = eta_crit[:, 1] - grad_eta * i_crit[:, 1]
            curr_den_lin = current_density[id_lin]
            eta_lin = grad_eta * curr_den_lin + b

            # curr_lin = current_density[id_lin[0]] \
            #     + current_density[id_lin] - self.i_crit[id_lin]
            # eta_lin = grad_eta * curr_lin + b

            eta_reg = \
                self.calc_electrode_loss_kulikovsky(current_density[id_reg],
                                                    conc_ele[id_reg], conc_ref,
                                                    update_members=False)
            eta = np.zeros(self.n_ele)
            eta[id_lin] = eta_lin
            eta[id_reg] = eta_reg
            return eta
        else:
            return self.calc_electrode_loss_kulikovsky(current_density,
                                                       conc_ele,
                                                       conc_ref)

    def calc_electrode_loss_kulikovsky(self, current_density: np.ndarray,
                                       conc: np.ndarray,
                                       conc_ref: float,
                                       update_members: bool = True):
        """
        Calculates the full voltage losses of the electrode
        """
        conc_star = conc / conc_ref
        var = 1. - current_density / (self.i_lim_star * conc_star)
        # var = np.where(var0 < 1e-4, 1e-4, var0)
        v_loss = np.zeros(current_density.shape)
        if self.calc_act_loss:
            v_loss_act = self.calc_activation_loss(current_density, conc_star)
            v_loss += v_loss_act
            # if update_members:
            #     self.v_loss_act[:] = v_loss_act
        if self.calc_gdl_diff_loss:
            v_loss_gdl_diff = self.calc_transport_loss_diffusion_layer(var)
            v_loss += v_loss_gdl_diff
            # if update_members:
            #     self.v_loss_gdl_diff[:] = v_loss_gdl_diff
        if self.calc_cl_diff_loss:
            v_loss_cl_diff = \
                self.calc_transport_loss_catalyst_layer(current_density,
                                                        var, conc)
            v_loss += v_loss_cl_diff
            # if update_members:
            #     self.v_loss_cl_diff[:] = v_loss_cl_diff
        return v_loss

    def calc_current_density(self, current_density: np.ndarray,
                             channel: chl.Channel, v_loss: np.ndarray):
        def func(curr_den: np.ndarray, over_pot: np.ndarray):
            return self.calc_electrode_loss(curr_den, channel) - over_pot
        return np.asarray(optimize.newton(func, current_density,
                                          args=(v_loss,)))
