# general imports
import warnings
import numpy as np
from scipy import optimize

# local module imports
from . import interpolation as ip, layers as layers, constants, \
    global_functions as g_func, flow_field as ff, \
    channel as chl
from .fluid import fluid as fluids

warnings.filterwarnings("ignore")


class HalfCell:

    # Class variables constant across all instances of the class
    # (under construction)

    def __init__(self, halfcell_dict, cell_dict, channel, number=None):
        self.number = number
        self.name = halfcell_dict['name']
        self.n_nodes = channel.n_nodes
        n_ele = self.n_nodes - 1
        self.n_ele = n_ele
        # Discretization in elements and nodes along the x-axis (flow axis)

        """half cell geometry parameter"""
        self.width = cell_dict["width"]
        self.length = cell_dict["length"]

        # Reference to channel object
        self.channel = channel
        self.channel.name = self.name + ' Channel'
        self.channel.fluid.name = \
            self.name + ' Fluid'  # + self.channel.fluid.TYPE_NAME

        # number of channels of each half cell
        self.n_channel = halfcell_dict['channel_number']

        flowfield_dict = {**cell_dict, **halfcell_dict}

        self.flow_field = \
            ff.FlowField(self.name + 'Flow Field', flowfield_dict, self.channel)

        # fuel must be at first position
        self.id_fuel = 0
        if isinstance(self.channel, chl.TwoPhaseMixtureChannel):
            self.id_h2o = self.channel.fluid.id_pc

        self.n_charge = halfcell_dict['charge_number']
        self.n_stoi = np.asarray(halfcell_dict['reaction_stoichiometry'])

        self.faraday = constants.FARADAY

        self.is_cathode = halfcell_dict['is_cathode']

        self.calc_act_loss = halfcell_dict['calc_act_loss']
        self.calc_cl_diff_loss = halfcell_dict['calc_cl_diff_loss']
        self.calc_gdl_diff_loss = halfcell_dict['calc_gdl_diff_loss']

        # thickness of the gas diffusion layer
        self.th_gdl = halfcell_dict['thickness_gdl']
        # thickness of the catalyst layer
        self.th_cl = halfcell_dict['thickness_cl']

        bpp_layer_dict = \
            {'thickness': halfcell_dict['thickness_bpp'],
             'width': self.flow_field.width_straight_channels,
             'length': self.flow_field.length_straight_channels,
             'electrical_conductivity':
                 halfcell_dict['electrical_conductivity_bpp'],
             'thermal_conductivity':
                 halfcell_dict['thermal_conductivity_bpp']}
        # 'porosity': self.channel.cross_area * self.n_channel / (
        #             self.th_bpp * self.width)}
        self.bpp = layers.SolidLayer(bpp_layer_dict, self.channel.dx)
        gde_layer_dict = \
            {'thickness': halfcell_dict['thickness_gdl']
                + halfcell_dict['thickness_cl'],
             'width': self.flow_field.width_straight_channels,
             'length': self.flow_field.length_straight_channels,
             'electrical_conductivity':
                 halfcell_dict['electrical_conductivity_gde'],
             'thermal_conductivity':
                 halfcell_dict['thermal_conductivity_gde']}
        # 'porosity':
        #    (self.th_gdl * halfcell_dict['porosity gdl']
        #     + self.th_cl * halfcell_dict['porosity cl'])
        #    / (self.th_gde + self.th_cl)}
        self.gde = layers.SolidLayer(gde_layer_dict, self.channel.dx)
        self.thickness = self.bpp.thickness + self.gde.thickness

        """voltage loss parameter, (Kulikovsky, 2013)"""
        # exchange current density
        vol_ex_cd = halfcell_dict['vol_ex_cd']
        # proton conductivity of the catalyst layer
        self.prot_con_cl = halfcell_dict['prot_con_cl']
        # diffusion coefficient of the reactant in the catalyst layer
        self.diff_coeff_cl = halfcell_dict['diff_coeff_cl']
        # diffusion coefficient of the reactant in the gas diffusion layer
        self.diff_coeff_gdl = halfcell_dict['diff_coeff_gdl']
        # tafel slope of the electrode
        self.tafel_slope = halfcell_dict['tafel_slope']
        # could use a better name see (Kulikovsky, 2013) not sure if 2-D
        # exchange current density
        self.i_sigma = np.sqrt(2. * vol_ex_cd * self.prot_con_cl
                               * self.tafel_slope)
        # index of the first element with negative cell voltage
        self.index_cat = self.n_nodes - 1
        # characteristic current density, see (Kulikovsky, 2013)
        self.i_star = self.prot_con_cl * self.tafel_slope / self.th_cl
        # concentration at channel inlet
        self.conc_in = None
        # limiting current density due to diffusion through the gdl
        # at channel inlet (calculated when inlet concentration is known)
        self.i_lim_star = None
        # numerical parameter for tangent line extension at limiting current
        self.conc_eps = halfcell_dict['c_eps']
        self.delta_i = halfcell_dict['delta_i']
        # critical local current density where Kulikovsky model transitions
        # into linear tangent line near limiting current
        self.i_crit = np.zeros(n_ele)

        # cell voltage loss
        self.v_loss = np.zeros(n_ele)
        self.updated_v_loss = False
        # boolean to hint if the cell voltage runs below zero
        # if HT-PEMFC True; if NT-PEMFC False
        self.break_program = False

        self.target_stoi = halfcell_dict['stoichiometry']
        # stoichiometry of the reactant at the channel inlet
        self.inlet_stoi = 0.0
        # cross water flux through the membrane
        self.w_cross_flow = np.zeros(n_ele)

        self.corrected_current_density = None

    def update(self, current_density, update_channel=False,
               current_control=True):
        """
        This function coordinates the program sequence
        """
        # self.calc_temp_fluid_ele()
        # mole_flow_in, mole_source = self.calc_mass_balance(current_density)
        if np.any(current_density < 0.0):
            raise ValueError('current density became smaller 0')
        if not current_control and self.updated_v_loss:
            self.corrected_current_density = \
                self.calc_current_density(current_density, self.v_loss)
        if current_control or self.corrected_current_density is None:
            corrected_current_density = current_density
        else:
            corrected_current_density = self.corrected_current_density
        if not self.break_program:
            # self.channel.update(mole_flow_in, mole_source)
            # self.channel.mole_flow[:] = mole_flow_in
            self.channel.mass_source[:], self.channel.mole_source[:] = \
                self.calc_mass_source(current_density)
            if update_channel:
                self.channel.update(update_mass=True, update_flow=False,
                                    update_heat=False, update_fluid=True)
            self.update_voltage_loss(corrected_current_density)

            # calculate stoichiometry
            current = np.sum(current_density * self.flow_field.active_area_dx)
            self.inlet_stoi = \
                self.channel.mole_flow[self.id_fuel, self.channel.id_in] \
                * self.faraday * self.n_charge \
                / (current * abs(self.n_stoi[self.id_fuel]))
            if current_control and self.inlet_stoi < 1.0:
                raise ValueError('stoichiometry of cell {0} '
                                 'becomes smaller than one: {1:0.3f}'
                                 .format(self.number, self.inlet_stoi))

    # def calc_mass_balance(self, current_density, stoi=None):
    #     n_species = self.channel.fluid.n_species
    #     mole_flow_in = np.zeros((n_species, self.n_nodes))
    #     mole_source = np.zeros((n_species, self.n_ele))
    #     mole_flow_in[self.id_fuel, :], mole_source[self.id_fuel, :] = \
    #         self.calc_fuel_flow(current_density, stoi)
    #     mole_flow_in[self.id_inert, :] = \
    #         mole_flow_in[self.id_fuel, self.channel.id_in] \
    #         * self.inert_reac_ratio
    #     air_flow_in = np.sum(mole_flow_in[:, self.channel.id_in])
    #     mole_flow_in[self.id_h2o, :], mole_source[self.id_h2o, :] = \
    #         self.calc_water_flow(current_density, air_flow_in)
    #     return mole_flow_in, mole_source

    def calc_mass_balance(self, current_density, stoi=None):
        avg_current_density = \
            np.average(current_density, weights=self.flow_field.active_area_dx)
        mass_flow_in, mole_flow_in = \
            self.calc_inlet_flow(avg_current_density, stoi)
        mass_flow_in = g_func.fill_transposed(mass_flow_in,
                                              self.channel.mass_flow.shape)
        mole_flow_in = g_func.fill_transposed(mole_flow_in,
                                              self.channel.mole_flow.shape)
        mass_source, mole_source = self.calc_mass_source(current_density)
        return mass_flow_in, mole_flow_in, mass_source, mole_source

    def calc_inlet_flow(self, current_density, stoi=None):
        if stoi is None:
            stoi = self.target_stoi
        if np.ndim(current_density) > 0:
            raise ValueError('current_density must be scalar')
        mole_flow_in = np.zeros(self.channel.fluid.n_species)
        mole_flow_in[self.id_fuel] = \
            current_density * self.flow_field.active_area \
            * stoi * abs(self.n_stoi[self.id_fuel]) \
            / (self.n_charge * self.faraday)
        inlet_composition = \
            self.channel.fluid.mole_fraction[:, self.channel.id_in]
        for i in range(len(mole_flow_in)):
            if i != self.id_fuel:
                mole_flow_in[i] = mole_flow_in[self.id_fuel] \
                    * inlet_composition[i] / inlet_composition[self.id_fuel]
        mass_flow_in = mole_flow_in * self.channel.fluid.species_mw
        return mass_flow_in, mole_flow_in

    def calc_mass_source(self, current_density):
        mole_source = np.zeros((self.channel.fluid.n_species, self.n_ele))

        for i in range(len(mole_source)):
            mole_source[i] = \
                current_density * self.flow_field.active_area_dx \
                * self.n_stoi[i] / (self.n_charge * self.faraday)

        # water cross flow
        water_cross_flow = self.flow_field.active_area_dx * self.w_cross_flow
        mole_source[self.id_h2o] += \
            self.flow_field.active_area_dx * self.w_cross_flow
        # self.channel.flow_direction
        mass_source = (mole_source.transpose()
                       * self.channel.fluid.species_mw).transpose()
        return mass_source, mole_source

    def calc_fuel_flow(self, current_density, stoi=None):
        """
        Calculates the reactant molar flow [mol/s]
        """
        if stoi is None:
            stoi = self.target_stoi
        curr_den = \
            np.average(current_density, weights=self.flow_field.active_area_dx)
        # curr_den = self.target_cd
        mol_flow_in = curr_den * self.flow_field.active_area * stoi \
            * abs(self.n_stoi[self.id_fuel]) / (self.n_charge * self.faraday)
        dmol = current_density * self.flow_field.active_area_dx \
            * self.n_stoi[self.id_fuel] / (self.n_charge * self.faraday)
        # g_func.add_source(self.mol_flow[self.id_fuel], dmol,
        #                   self.flow_direction)
        return mol_flow_in, dmol

    def calc_water_flow(self, current_density, air_flow_in):
        """"
        Calculates the water molar flow [mol/s]
        """
        if not isinstance(self.channel.fluid, fluids.TwoPhaseMixture):
            raise TypeError('Fluid in channel must be of type TwoPhaseMixture')
        id_in = self.channel.id_in
        humidity_in = self.channel.fluid.humidity[id_in]
        sat_p = self.channel.fluid.saturation_pressure[id_in]
        mol_flow_in = air_flow_in * sat_p * humidity_in / \
            (self.channel.pressure[id_in] - humidity_in * sat_p)
        dmol = np.zeros_like(current_density)
        h2o_prod = self.flow_field.active_area_dx * self.n_stoi[self.id_h2o] \
            * current_density / (self.n_charge * self.faraday)
        dmol += h2o_prod
        h2o_cross = self.flow_field.active_area_dx * self.w_cross_flow
        # * self.channel.flow_direction
        dmol += h2o_cross
        return mol_flow_in, dmol

    def update_voltage_loss(self, current_density):
        eta = self.calc_electrode_loss(current_density)
        self.v_loss[:] = eta \
            + self.calc_plate_loss(current_density)
        self.updated_v_loss = True

    def calc_plate_loss(self, current_density):
        current = current_density * self.flow_field.active_area_dx
        v_loss_bpp = current / self.bpp.electrical_conductance[0]
        # self.v_loss_bpp[:] = current / self.bpp.electrical_conductance[0]
        return v_loss_bpp

    def calc_activation_loss(self, current_density, conc):
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

    def calc_transport_loss_catalyst_layer(self, current_density, var, conc):
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

    def calc_transport_loss_diffusion_layer(self, var):
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

    def calc_electrode_loss(self, current_density):
        if hasattr(self.channel.fluid, 'gas'):
            conc = self.channel.fluid.gas.concentration[self.id_fuel]
        else:
            conc = self.channel.fluid.concentration[self.id_fuel]
        conc_ele = ip.interpolate_1d(conc)
        conc_ref = conc[self.channel.id_in]
        conc_star = conc_ele / conc_ref
        # if self.channel.flow_direction == 1:
        #     conc_in = conc[:-1]
        # else:
        #     conc_in = conc[1:]
        conc_in = conc[self.channel.id_in]
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

    def calc_electrode_loss_kulikovsky(self, current_density, conc, conc_ref,
                                       update_members=True):
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

    def calc_current_density(self, current_density, v_loss):
        def func(curr_den, over_pot):
            return self.calc_electrode_loss(curr_den) \
                   + self.calc_plate_loss(curr_den) - over_pot
        return optimize.newton(func, current_density, args=(v_loss, ))
