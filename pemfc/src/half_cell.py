# General imports
import warnings
import numpy as np

# Local module imports
from . import transport_layer as tl, constants, flow_field as ff, channel as chl
from .fluid import fluid as fluids
from . import electrochemistry as electrochem
from . import interpolation as ip
from . import discretization as dsct
from . import diffusion_transport as diff

warnings.filterwarnings("ignore")


class HalfCell:

    def __init__(self, halfcell_dict, cell_dict, channel, number=None):
        self.number = number
        self.name = halfcell_dict['name']
        # Discretization in elements and nodes along the x-axis (flow axis)
        self.n_nodes = channel.n_nodes
        n_ele = self.n_nodes - 1
        self.n_ele = n_ele

        # Half cell geometry parameter
        self.width = cell_dict["width"]
        self.length = cell_dict["length"]

        # Reference to channel object
        self.channel = channel
        self.channel.name = self.name + ' Channel'
        self.channel.fluid.name = \
            self.name + ' Fluid'  # + self.channel.fluid.TYPE_NAME

        # Number of channels of each half cell
        self.n_channels = halfcell_dict['channel_number']

        # Fuel must be at first position
        self.id_fuel = 0
        if isinstance(self.channel, chl.TwoPhaseMixtureChannel):
            self.id_h2o = self.channel.fluid.id_pc
        self.faraday = constants.FARADAY

        # Initialize flow field geometry
        flowfield_dict = {
            'channel_number': self.n_channels,
            'rib_width': halfcell_dict['rib_width'],
            'width': self.width,
            'length': self.length
        }
        self.flow_field = \
            ff.FlowField(self.name + 'Flow Field', flowfield_dict, self.channel)

        # Discretization settings
        # Additional resolution in width direction in regions
        # "under land" and "under channel" if True
        self.channel_land_discretization = \
            cell_dict['channel_land_discretization']
        discretization_shape = self.channel.dx.shape
        if self.channel_land_discretization is True:
            discretization_shape += (2, )
            land_channel_ratio = self.flow_field.rib_width / self.channel.width
        else:
            land_channel_ratio = 1.0
        discretization_dict = {
            'shape': discretization_shape,
            'ratio': (1.0, land_channel_ratio),
            'width': self.flow_field.width_straight_channels / self.n_channels,
            'length': self.flow_field.length_straight_channels,
        }
        self.discretization = dsct.Discretization2D(
            discretization_dict=discretization_dict)

        # Initialize electrochemistry model
        electrochemistry_dict = halfcell_dict['electrochemistry']
        electrochemistry_dict['fuel_index'] = self.id_fuel
        self.electrochemistry = electrochem.ElectrochemistryModel(
            electrochemistry_dict, self.discretization)

        # Initialize bipolar plate (bpp)
        bpp_dict = halfcell_dict['bpp']
        bpp_dict['name'] = self.name + ' BPP'
        # 'porosity': self.channel.cross_area * self.n_channel / (
        #             self.th_bpp * self.width)}
        bpp_transport_properties = {
            'thermal': bpp_dict['thermal_conductivity'],
            'electrical': bpp_dict['electrical_conductivity']}
        self.bpp = tl.TransportLayer2D(bpp_dict, bpp_transport_properties,
                                       self.discretization)

        # Initialize gas diffusion electrode (gde: gdl + cl)
        gde_dict = halfcell_dict['gde']
        gde_dict.update(
            {'name': self.name + ' GDE',
             'thickness': electrochemistry_dict['thickness_gdl']
                + electrochemistry_dict['thickness_cl']})
        gde_transport_properties = {
            'thermal': gde_dict['thermal_conductivity'],
            'electrical': gde_dict['electrical_conductivity']}
        # 'porosity':
        #    (self.th_gdl * halfcell_dict['porosity gdl']
        #     + self.th_cl * halfcell_dict['porosity cl'])
        #    / (self.th_gde + self.th_cl)}
        self.gde = tl.TransportLayer2D(gde_dict, gde_transport_properties,
                                       self.discretization)

        # Initialize diffusion transport model
        gdl_dict = gde_dict.copy()
        gdl_dict.update(
            {'name': self.name + ' GDL',
             'thickness': electrochemistry_dict['thickness_gdl']})

        self.gdl_diffusion = diff.DiffusionTransport.create(
            gdl_dict, self.channel.fluid,
            dsct.Discretization3D.create_from(self.discretization,
                                              gdl_dict['thickness'], 2))

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
        self.w_cross_flow = np.zeros(self.gde.discretization.shape)
        # voltage loss
        self.voltage_loss = np.zeros(self.gde.discretization.shape)

    def update(self, current_density: np.ndarray,
               temperature: np.ndarray, update_channel=False,
               current_control=True):
        """
        This function coordinates the program sequence
        """
        if not self.break_program:
            # self.channel.update(mole_flow_in, mole_source)
            # self.channel.mole_flow[:] = mole_flow_in
            current = self.surface_flux_to_channel_source(current_density)
            mass_source, mole_source = self.calc_mass_source(current_density)
            concentration = self.channel.fluid.concentration
            self.gdl_diffusion.update(temperature, self.channel.pressure,
                                      concentration, mole_source)

            self.electrochemistry.update(current_density, self.channel)
            # self.channel.mass_source[:], self.channel.mole_source[:] = (
            #     mass_source, mole_source)
            if update_channel:
                self.channel.update(update_mass=True, update_flow=False,
                                    update_heat=False, update_fluid=True)
            self.update_voltage_loss(current_density)

            # Calculate stoichiometry
            self.inlet_stoi = (
                self.channel.mole_flow[self.id_fuel, self.channel.id_in]
                * self.faraday * self.n_charge
                / (np.sum(current) * abs(self.n_stoi[self.id_fuel])))
            if current_control and self.inlet_stoi < 1.0:
                raise ValueError('stoichiometry of cell {0} '
                                 'becomes smaller than one: {1:0.3f}'
                                 .format(self.number, self.inlet_stoi))

    def calc_inlet_flow(self, current_density, stoi=None):
        if stoi is None:
            stoi = self.target_stoi
        if np.ndim(current_density) > 0:
            raise ValueError('current_density must be scalar')
        mole_flow_in = np.zeros(self.channel.fluid.n_species)
        current = np.sum(self.surface_flux_to_channel_source(current_density))
        mole_flow_in[self.id_fuel] = \
            current * stoi * abs(self.n_stoi[self.id_fuel]) \
            / (self.n_charge * self.faraday)
        inlet_composition = \
            self.channel.fluid.mole_fraction[:, self.channel.id_in]
        for i in range(len(mole_flow_in)):
            if i != self.id_fuel:
                mole_flow_in[i] = mole_flow_in[self.id_fuel] \
                    * inlet_composition[i] / inlet_composition[self.id_fuel]
        mass_flow_in = mole_flow_in * self.channel.fluid.species_mw
        return mass_flow_in, mole_flow_in

    def calc_mass_source(self, current_density: np.ndarray):
        if not isinstance(current_density, np.ndarray):
            current_density = np.asarray(current_density)
        mole_source = np.zeros((self.channel.fluid.n_species,
                                *current_density.shape))

        for i in range(len(mole_source)):
            mole_source[i] = (current_density * self.n_stoi[i]
                              / (self.n_charge * self.faraday))

        # Water cross flow
        # water_cross_source = self.surface_flux_to_channel_source(
        #     self.w_cross_flow)
        mole_source[self.id_h2o] += self.w_cross_flow
        # self.channel.flow_direction
        mass_source = (mole_source.transpose()
                       * self.channel.fluid.species_mw).transpose()
        return mass_source, mole_source

    def surface_flux_to_channel_source(self, flux: np.ndarray):
        if np.isscalar(flux) or flux.ndim in (0, 1):
            return np.sum(self.discretization.d_area, axis=0) * flux
        elif flux.ndim == 2:
            return np.sum(self.discretization.d_area * flux, axis=1)
        else:
            raise ValueError('flux variable must be one- or'
                             'two-dimensional numpy array')

    def update_voltage_loss(self, current_density: np.ndarray):
        """
        Calculates voltage loss at electrode. Solid layer losses are set to zero
        at the current implementation, due to its inclusion in the overall
        electrical system calculation at stack level.
        Args:
            current_density: numpy array of local area-specific current density
                             distribution
        Returns: None
        """
        # area = self.discretization.d_area
        bpp_loss = 0.0  # self.bpp.calc_voltage_loss(current_density, area)
        gde_loss = 0.0  # self.gde.calc_voltage_loss(current_density, area)
        self.voltage_loss[:] = self.electrochemistry.v_loss + bpp_loss + gde_loss

    @staticmethod
    def calc_faraday_flow(fluid, current, stoichiometry, reaction_stoichiometry,
                          charge_number, reactant_index=0):
        """
        Calculates the corresponding species mass and molar flows

        Args:
            fluid: object of type GasMixture, CanteraGasMixture, TwoPhaseMixture,
                   or CanteraTwoPhaseMixture from module pemfc.fluid.fluid
            current: scalar value providing electrical current (A)
            stoichiometry: scalar for flow stoichiometry
            reaction_stoichiometry: scalar for reaction stoichiometry,
                                    moles of reactants used in reaction balance
            charge_number: number of electron charges transferred in reaction
                           balance
            reactant_index: index in species array for the reactant species
        Returns: mass_flow, mole_flow (species array according to fluid object)
        """
        if not isinstance(fluid, (fluids.GasMixture, fluids.TwoPhaseMixture,
                                  fluids.CanteraGasMixture,
                                  fluids.CanteraTwoPhaseMixture)):
            raise TypeError(
                'Parameter "fluid" must be object of type GasMixture, '
                'TwoPhaseMixture, CanteraGasMixture, or CanteraTwoPhaseMixture '
                'from module pemfc.fluid.fluid')
        mole_flow = np.zeros(fluid.n_species)

        mole_flow[reactant_index] = (
                current * stoichiometry * abs(reaction_stoichiometry)
                / (charge_number * constants.FARADAY))
        composition = fluid.mole_fraction[:, 0]
        for i in range(len(mole_flow)):
            if i != reactant_index:
                mole_flow[i] = mole_flow[reactant_index] \
                               * composition[i] / composition[reactant_index]
        mass_flow = mole_flow * fluid.species_mw
        return mass_flow, mole_flow

    def calc_humidity(self):
        """
        For now this is a helper function to map the channel humidity to the
        gde-membrane interface in the correct discretization. In the future, a
        approximation according to composition gradients and thermodynamics
        should be implemented
        """
        humidity = ip.interpolate_1d(self.channel.fluid.humidity)
        return np.asarray([humidity for i in
                           range(self.discretization.shape[1])]).transpose()
