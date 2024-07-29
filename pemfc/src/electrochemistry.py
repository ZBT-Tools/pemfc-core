# General imports
from abc import ABC
import copy
import numpy as np
from scipy import optimize

# local module imports
from . import interpolation as ip, constants, channel as chl
from . import discretization as dsct
from . import global_state as gs


class ElectrochemistryModel(ABC):

    def __init__(self, input_dict: dict, discretization: dsct.Discretization2D):
        
        """voltage loss parameter, (Kulikovsky, 2013)"""

        self.discretization = discretization
        shape = discretization.shape
        self.faraday = constants.FARADAY
        # Index of fuel component in gas mixture
        self.id_fuel = input_dict['fuel_index']
        # Charge number of reaction
        self.n_charge = input_dict['charge_number']

        self.th_gdl = input_dict['thickness_gdl']
        self.th_cl = input_dict['thickness_cl']

        # Proton conductivity of the catalyst layer
        self.prot_con_cl = input_dict['prot_con_cl']
        # Diffusion coefficient of the reactant in the catalyst layer
        self.diff_coeff_cl = input_dict['diff_coeff_cl']
        # Diffusion coefficient of the reactant in the gas diffusion layer
        self.diff_coeff_gdl = np.full(shape, input_dict['diff_coeff_gdl'])

        # Tafel slope of the electrode
        self.tafel_slope = input_dict['tafel_slope']

        # Could use a better name see (Kulikovsky, 2013) not sure if 2-D
        # exchange current density
        vol_ex_cd = input_dict['vol_ex_cd']
        self.i_sigma = np.sqrt(2. * vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        # Index of the first element with negative cell voltage
        self.index_cat = shape[0]
        # Characteristic current density, see (Kulikovsky, 2013)
        self.i_star = self.prot_con_cl * self.tafel_slope / self.th_cl

        # Numerical parameters for tangent line extension at limiting current
        self.conc_eps = input_dict['c_eps']
        self.delta_i = input_dict['delta_i']
        # Critical local current density where Kulikovsky model transitions
        # into linear tangent line near limiting current
        self.i_crit = np.zeros(shape)
        # Switch to include activation contribution to voltage losses
        self.calc_act_loss = input_dict['calc_act_loss']
        # Switch to include catalyst layer diffusion contribution to voltage
        # losses
        self.calc_cl_diff_loss = input_dict['calc_cl_diff_loss']
        # Switch to include gas diffusion layer diffusion contribution to
        # voltage losses
        self.calc_gdl_diff_loss = input_dict['calc_gdl_diff_loss']

        # Cell voltage loss
        self.v_loss = np.zeros(shape)
        self.updated_v_loss = False
        self.corrected_current_density = None

    def update(self, current_density: np.ndarray, concentration: np.ndarray,
               reference_concentration: float,
               current_control: bool = True, **kwargs):
        """
        This function coordinates the program sequence
        """
        # self.calc_temp_fluid_ele()
        # mole_flow_in, mole_source = self.calc_mass_balance(current_density)
        if np.any(current_density < 0.0):
            raise ValueError('current density became smaller 0')
        if not current_control and self.updated_v_loss:
            self.corrected_current_density = \
                self.calc_current_density(current_density,
                                          concentration,
                                          reference_concentration, self.v_loss)
        if current_control or self.corrected_current_density is None:
            corrected_current_density = current_density
        else:
            corrected_current_density = self.corrected_current_density
        self.update_voltage_loss(corrected_current_density, concentration,
                                 reference_concentration, **kwargs)

    def update_voltage_loss(self, current_density: np.ndarray,
                            concentration: np.ndarray,
                            reference_concentration: float, **kwargs):
        eta = self.calc_electrode_loss(current_density, concentration,
                                       reference_concentration, **kwargs)
        self.v_loss[:] = eta
        self.updated_v_loss = True

    def calc_electrode_loss(self, current_density: np.ndarray,
                            concentration: np.ndarray,
                            reference_concentration: (float, np.ndarray),
                            inlet_concentration: (float, np.ndarray) = None,
                            **kwargs):
        # conc_ele = ip.interpolate_1d(concentration)
        # conc_ele = np.asarray([
        #     conc_ele for i in range(current_density.shape[-1])]).transpose()
        # conc = concentration
        concentration[concentration <= 0.0] = 1e-3
            # i_lim_star = (self.n_charge * self.faraday * concentration
            #               * self.diff_coeff_cl / self.th_cl)
        if 'scaling_factors' in kwargs:
            scaling_factors = kwargs['scaling_factors']
            i_lim_star = np.divide(current_density, scaling_factors,
                                   out=np.zeros(current_density.shape),
                                   where=scaling_factors != 0)
            i_crit = i_lim_star * 1.1
            gdl_loss = False
            # diff_coeff_gdl_by_length = kwargs['scaling_factors']
            # i_lim_star = (self.n_charge * self.faraday * inlet_concentration
            #               * diff_coeff_gdl_by_length)
            # i_crit = (i_lim_star * (concentration - self.conc_eps)
            #           / reference_concentration)
            # gdl_loss = True
        else:
            # Additional in-plane diffusion resistance for
            # channel-land-discretization with channel located at index 1 and
            # land located at index 0
            diff_coeff_gdl_by_length = self.diff_coeff_gdl / self.th_gdl
            if current_density.shape[-1] == 2:
                diff_coeff_gdl_by_length[:, 0] = (
                        1.0 / (self.th_gdl + self.discretization.dx[-1, :, 0] /
                               self.diff_coeff_gdl[:, 0])
                )
                # diff_coeff_gdl_by_length[:, 0] = (
                #         1.0 / (1.0 * self.th_gdl /
                #                self.diff_coeff_gdl[:, 0])
                # )
            # Limiting current density due to diffusion through the gdl
            # at inlet (calculated when inlet concentration is known)
            i_lim_star = (self.n_charge * self.faraday * inlet_concentration
                          * diff_coeff_gdl_by_length)
            i_crit = (i_lim_star * (concentration - self.conc_eps)
                      / reference_concentration)
            # i_crit = 0.95 * i_lim_star
            gdl_loss = True

        # TODO: Algorithm needs to be stabilized or dampened for convergence,
        #  otherwise just fluctuates between minimal and maximal currents each
        #  iteration when approaching limiting currents
        id_lin = np.nonzero(current_density > i_crit)
        id_reg = np.nonzero(current_density <= i_crit)
        if id_lin[0].size:
            i_crit_lin = np.copy(i_crit[id_lin])
            conc_crit = concentration[id_lin]
            conc_crit_stack = np.stack((conc_crit, conc_crit),
                                       axis=conc_crit.ndim)
            delta_i = self.delta_i * np.average(i_crit_lin)
            i_crit_stack = np.stack(
                (i_crit_lin - delta_i, i_crit_lin + delta_i),
                axis=i_crit_lin.ndim)
            i_lim_star_lin = np.copy(i_lim_star[id_lin])
            i_lim_star_stack = np.stack((i_lim_star_lin, i_lim_star_lin),
                                        axis=i_lim_star_lin.ndim)
            eta_crit_stack = self.calc_electrode_loss_kulikovsky(
                i_crit_stack, conc_crit_stack, reference_concentration,
                i_lim_star_stack, gdl_loss=gdl_loss)

            grad_eta = np.moveaxis(
                np.gradient(eta_crit_stack, delta_i, axis=-1), -1, 0)[0]
            eta_crit = np.average(eta_crit_stack, axis=-1)
            b = eta_crit - grad_eta * i_crit_lin
            eta_lin = grad_eta * current_density[id_lin] + b

            # curr_lin = current_density[id_lin[0]] \
            #     + current_density[id_lin] - self.i_crit[id_lin]
            # eta_lin = grad_eta * curr_lin + b
            i_lim_star_reg = i_lim_star[id_reg]
            eta_reg = self.calc_electrode_loss_kulikovsky(
                current_density[id_reg], concentration[id_reg],
                reference_concentration, i_lim_star_reg, gdl_loss=gdl_loss)
            eta = np.zeros(current_density.shape)
            eta[id_lin] = eta_lin
            eta[id_reg] = eta_reg
            return eta
        else:
            eta = self.calc_electrode_loss_kulikovsky(
                current_density, concentration, reference_concentration,
                i_lim_star, gdl_loss=gdl_loss)
            # if gs.global_state.iteration == 50:
            #     print('test')
            return eta

    def calc_electrode_loss_kulikovsky(
            self, current_density: np.ndarray, conc: np.ndarray,
            conc_ref: float, i_lim_star: np.ndarray, gdl_loss=True):
        """
        Calculates the full voltage losses of the electrode
        """
        conc_star = conc / conc_ref
        if gdl_loss:
            var = np.maximum(1. - current_density / (i_lim_star * conc_star),
                             1e-6)
        else:
            self.calc_gdl_diff_loss = False
            var = 1.0

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
                                           var: np.ndarray, conc: np.ndarray):
        """
        Calculates the diffusion voltage loss in the catalyst layer
        according to (Kulikovsky, 2013).
        """
        # try:
        i_hat = current_density / self.i_star
        i_hat[i_hat <= 0.0] = 1e-3
        short_save = np.sqrt(2. * i_hat)
        beta = \
            short_save / (1. + np.sqrt(1.12 * i_hat) * np.exp(short_save)) \
            + np.pi * i_hat / (2. + i_hat)
        # except FloatingPointError:
        #     test = np.any(current_density < 0.0)
        #     raise
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

    def calc_transport_loss_diffusion_layer(self, var: np.ndarray):
        """
        Calculates the diffusion voltage loss in the gas diffusion layer
        according to (Kulikovsky, 2013).
        """
        if np.all(var > 0.0):
            v_loss_gdl_diff = - self.tafel_slope * np.log10(var)
        else:
            raise ValueError('value must not become equal or smaller zero')
        # nan_list = np.isnan(self.v_loss_gdl_diff)
        # if nan_list.any():
        #     v_loss_gdl_diff[np.argwhere(nan_list)[0, 0]:] = 1.e50
        return v_loss_gdl_diff

    def calc_current_density(self, current_density: np.ndarray,
                             concentration: np.ndarray,
                             reference_concentration: float,
                             v_loss: np.ndarray):
        def func(curr_den: np.ndarray, over_pot: np.ndarray):
            return self.calc_electrode_loss(curr_den, concentration,
                                            reference_concentration) - over_pot
        return np.asarray(optimize.newton(func, current_density,
                                          args=(v_loss,)))


