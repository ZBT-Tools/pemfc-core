# general imports
import warnings
import numpy as np

# local module imports
from . import solid_layer as sl, constants, \
    global_functions as g_func, flow_field as ff, \
    channel as chl
from .fluid import fluid as fluids
from . import electrochemistry as electrochem

warnings.filterwarnings("ignore")


class HalfCell:

    def __init__(self, halfcell_dict, cell_dict, channel, number=None):
        self.number = number
        self.name = halfcell_dict['name']
        # Discretization in elements and nodes along the x-axis (flow axis)
        self.n_nodes = channel.n_nodes
        n_ele = self.n_nodes - 1
        self.n_ele = n_ele

        # Additional resolution in width direction in regions
        # "under land" and "under channel" if True
        self.channel_land_discretization = cell_dict['channel_land_discretization']

        # Half cell geometry parameter
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

        self.faraday = constants.FARADAY

        # Initialize electrochemistry model
        electrochemistry_dict = halfcell_dict['electrochemistry']
        electrochemistry_dict['fuel_index'] = self.id_fuel
        self.electrochemistry = electrochem.ElectrochemistryModel(
            electrochemistry_dict, self.n_nodes)

        # Initialize bipolar plate (bpp)
        bpp_dict = halfcell_dict['bpp']
        bpp_dict.update(
            {'name': self.name + ' BPP',
             'width': self.flow_field.width_straight_channels,
             'length': self.flow_field.length_straight_channels})
        # 'porosity': self.channel.cross_area * self.n_channel / (
        #             self.th_bpp * self.width)}

        layer_discretization = self.channel.dx.shape
        if self.channel_land_discretization is True:
            layer_discretization += (2, )
        self.bpp = sl.SolidLayer(bpp_dict, layer_discretization)

        # initialize gas diffusion electrode (gde: gdl + cl
        gde_dict = halfcell_dict['gde']
        gde_dict.update(
            {'name': self.name + ' GDE',
             'thickness': electrochemistry_dict['thickness_gdl']
                + electrochemistry_dict['thickness_cl'],
             'width': self.flow_field.width_straight_channels,
             'length': self.flow_field.length_straight_channels})
        # 'porosity':
        #    (self.th_gdl * halfcell_dict['porosity gdl']
        #     + self.th_cl * halfcell_dict['porosity cl'])
        #    / (self.th_gde + self.th_cl)}
        self.gde = sl.SolidLayer(gde_dict, layer_discretization)
        self.thickness = self.bpp.thickness + self.gde.thickness

        self.n_charge = self.electrochemistry.n_charge
        self.n_stoi = \
            np.asarray(halfcell_dict['reaction_stoichiometry'])

        # boolean to hint if the cell voltage runs below zero
        # if HT-PEMFC True; if NT-PEMFC False
        self.break_program = False

        self.target_stoi = halfcell_dict['stoichiometry']
        # stoichiometry of the reactant at the channel inlet
        self.inlet_stoi = 0.0
        # cross water flux through the membrane
        self.w_cross_flow = np.zeros(layer_discretization)
        # voltage loss
        self.v_loss = np.zeros(layer_discretization)

    def update(self, current_density, update_channel=False,
               current_control=True):
        """
        This function coordinates the program sequence
        """
        if not self.break_program:
            self.electrochemistry.update(current_density, self.channel)
            # self.channel.update(mole_flow_in, mole_source)
            # self.channel.mole_flow[:] = mole_flow_in
            self.channel.mass_source[:], self.channel.mole_source[:] = \
                self.calc_mass_source(current_density)
            if update_channel:
                self.channel.update(update_mass=True, update_flow=False,
                                    update_heat=False, update_fluid=True)
            self.update_voltage_loss(current_density)

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

    def update_voltage_loss(self, current_density: np.ndarray):
        area = self.flow_field.active_area_dx
        bpp_loss = self.bpp.calc_voltage_loss(current_density, area)
        gde_loss = self.gde.calc_voltage_loss(current_density, area)
        self.v_loss[:] = self.electrochemistry.v_loss \
            + bpp_loss \
            + gde_loss

    @staticmethod
    def calc_faraday_flow(fluid, current, stoichiometry, reaction_stoichiometry,
                          charge_number, reactant_index=0):
        """
        Calculates the corresponding species mass and molar flows
        :param fluid: object of type GasMixture, CanteraGasMixture, TwoPhaseMixture,
        or CanteraTwoPhaseMixture from module pemfc.fluid.fluid
        :param current: scalar value providing electrical current (A)
        :param stoichiometry: scalar for flow stoichiometry
        :param reaction_stoichiometry: scalar for reaction stoichiometry,
        moles of reactants used in reaction balance
        :param charge_number: number of electron charges transferred in reaction balance
        :param reactant_index: index in species array for the reactant species
        :return: mass_flow, mole_flow (species array according to fluid object)
        """
        if not isinstance(fluid, (fluids.GasMixture, fluids.TwoPhaseMixture,
                                  fluids.CanteraGasMixture, fluids.CanteraTwoPhaseMixture)):
            raise TypeError('Parameter "fluid" must be object of type GasMixture, '
                            ' TwoPhaseMixture, CanteraGasMixture, or CanteraTwoPhaseMixture '
                            'from module pemfc.fluid.fluid')
        mole_flow = np.zeros(fluid.n_species)

        mole_flow[reactant_index] = current * stoichiometry \
                                    * abs(reaction_stoichiometry) / (
                                                charge_number * constants.FARADAY)
        composition = fluid.mole_fraction[:, 0]
        for i in range(len(mole_flow)):
            if i != reactant_index:
                mole_flow[i] = mole_flow[reactant_index] \
                               * composition[i] / composition[reactant_index]
        mass_flow = mole_flow * fluid.species_mw
        return mass_flow, mole_flow

