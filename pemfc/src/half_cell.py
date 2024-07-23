# General imports
import warnings
import numpy as np
from scipy import ndimage

# Local module imports
from . import transport_layer as tl, constants, flow_field as ff, channel as chl
from .fluid import fluid as fluids
from . import electrochemistry as electrochem
from . import interpolation as ip
from . import discretization as dsct
from . import diffusion_transport as diff
from . import global_functions as gf
from . import global_state as gs
from . import porous_two_phase_flow as p2pf

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
        else:
            raise NotImplementedError('fluid for HalfCell class must be a '
                                      'TwoPhaseMixture at the moment')

        inert_ids = [i for i in range(self.channel.fluid.n_species)
                     if i not in (self.id_fuel, self.id_h2o)]
        # self.id_inert = None
        if inert_ids:
            self.id_inert = inert_ids[-1]
        else:
            raise ValueError('species in fluid the HalfCell class must '
                             'include an inert species')

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
        electrochemistry_dict['thickness_gdl'] = (
            halfcell_dict)['gde']['thickness']
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
        gde_dict = halfcell_dict['gde'].copy()
        gde_dict.update(
            {'name': self.name + ' GDE',
             'thickness': gde_dict['thickness']
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

        self.calc_gdl_diffusion = True
        self.calc_two_phase_flow = True

        if self.calc_gdl_diffusion:
            gdl_diffusion_dict = gde_dict.copy()
            gdl_diffusion_dict.update(
                {'name': self.name + ' GDL',
                 'type': 'GasMixture',
                 'thickness': halfcell_dict['gde']['thickness'],
                 'boundary_patches': {
                     'Neumann': {'axes': (0,),
                                 'indices': (-1,)}}})
            # gdl_diffusion_dict['diffusion_coefficient'] = (
            #     electrochemistry_dict['diff_coeff_gdl'])
            if self.channel_land_discretization:
                # Initialize diffusion transport model
                nx = 2
                ny = self.discretization.shape[0]
                nz = 2
                gdl_diffusion_dict['boundary_patches']['Dirichlet'] = {
                    'axes': (0, 2), 'indices': (0, list(range(int(nz/2), nz)))}
            else:
                nx = 2
                ny = self.discretization.shape[0]
                nz = 1
                gdl_diffusion_dict['boundary_patches']['Dirichlet'] = {
                    'axes': (0,), 'indices': (0,)}
            discretization_dict = {
                'shape': (nx, ny, nz),
                'depth': gdl_diffusion_dict['thickness'],
                'length': self.discretization.length[0],
                'width': self.discretization.length[1] * 0.5}
            gdl_discretization_2d = dsct.Discretization(discretization_dict)
            gdl_discretization = dsct.Discretization3D.create_from(
                    gdl_discretization_2d, gdl_diffusion_dict['thickness'])
            self.gdl_diffusion = diff.DiffusionTransport.create(
                gdl_diffusion_dict, gdl_discretization,
                self.channel.fluid, self.id_inert)
            self.explicit_source_terms = np.zeros(
                self.gdl_diffusion.solution_array.shape)
            self.implicit_source_terms = np.zeros(
                self.gdl_diffusion.solution_array.shape)
            if self.calc_two_phase_flow:
                gdl_two_phase_flow_dict = gdl_diffusion_dict.copy()
                gdl_two_phase_flow_dict['type'] = 'DarcyFlow'
                self.two_phase_flow = p2pf.TwoPhaseMixtureDiffusionTransport(
                    gdl_two_phase_flow_dict, gdl_discretization,
                    self.gdl_diffusion.fluid)

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
               temperature: np.ndarray, update_channel=True,
               current_control=True):
        """
        This function coordinates the program sequence
        """

        if not self.break_program:
            # self.channel.update(mole_flow_in, mole_flux)
            # self.channel.mole_flow[:] = mole_flow_in
            mass_flux, mole_flux = self.calc_species_flux(current_density)
            channel_concentration = self.channel.fluid.concentration
            reference_fuel_concentration = np.max(
                channel_concentration[self.id_fuel, self.channel.id_in])
            if self.calc_gdl_diffusion:
                # Explicitly calculate diffusion in GDL and use calculated
                # concentrations at CL-GDL-Interface for the electrochemistry
                # model, therefore the electrochemistry model does not need
                # to account for GDL losses anymore
                if not self.calc_two_phase_flow:
                    self.gdl_diffusion.update(
                        temperature, self.channel.pressure,
                        channel_concentration, mole_flux,
                        self.explicit_source_terms)
                if self.calc_two_phase_flow:
                    # Adjust two-phase flow boundary conditions
                    liquid_flux_ratio = 0.5
                    ch_gdl_sat = (np.ones(self.two_phase_flow.saturation.shape)
                                  * 0.01)
                    liquid_water_flux = (mass_flux[self.id_h2o]
                                         * liquid_flux_ratio)
                    mole_flux[self.id_h2o] *= (1.0 - liquid_flux_ratio)
                    iteration = 0
                    molar_evap_rate = np.zeros(
                        self.two_phase_flow.evaporation_rate.shape)
                    molar_cond_coeff = np.zeros(
                        self.two_phase_flow.implicit_condensation_coeff.shape)

                    # Adjust numerical settings
                    error = np.inf
                    max_iterations = 10
                    min_iterations = 1
                    urf = 0.5
                    error_tolerance = 1e-5
                    while (all((error > error_tolerance,
                                iteration < max_iterations))
                           or iteration < min_iterations):
                        # Update gas diffusion in GDL
                        self.gdl_diffusion.update(
                            temperature, self.channel.pressure,
                            channel_concentration, mole_flux,
                            self.explicit_source_terms,
                            self.implicit_source_terms)

                        # Update two-phase flow in GDL
                        temp = self.gdl_diffusion.fluid.temperature
                        press = self.gdl_diffusion.fluid.pressure
                        mol_comp = self.gdl_diffusion.solution_array
                        self.two_phase_flow.update(
                            temp, press, mol_comp, ch_gdl_sat,
                            liquid_water_flux, update_fluid=True)
                        fluid = self.two_phase_flow.fluid
                        molar_mass = fluid.species_mw[fluid.id_pc]
                        molar_evap_rate_old = np.copy(molar_evap_rate)
                        molar_cond_coeff_old = np.copy(molar_cond_coeff)
                        molar_evap_rate[:] = (
                                self.two_phase_flow.evaporation_rate /
                                molar_mass)
                        molar_evap_rate = (urf * molar_evap_rate
                                           + (1.0 - urf) * molar_evap_rate_old)
                        molar_cond_coeff[:] = (
                                self.two_phase_flow.implicit_condensation_coeff
                                / molar_mass)
                        molar_cond_coeff = (
                            urf * molar_cond_coeff
                            + (1.0 - urf) * molar_cond_coeff_old)
                        self.explicit_source_terms[fluid.id_pc] = (
                            1.0 * molar_evap_rate
                                )
                        self.implicit_source_terms[fluid.id_pc] = (
                            -1.0 * molar_cond_coeff)
                        # Error calculation
                        diff_e = molar_evap_rate - molar_evap_rate_old
                        diff_e[:] = np.divide(diff_e, molar_evap_rate,
                                              where=molar_evap_rate != 0.0)
                        diff_e = diff_e.ravel()
                        error = (np.dot(diff_e.transpose(), diff_e)
                                 / (2.0 * len(diff_e)))
                        iteration += 1
                # Reshape 3D-concentration fields from the GDL-Diffusion
                # sub-model to the reduced discretization in this model
                # Take only the concentration at the CL-Interface
                # (axis: 1, last index)
                fuel_cl_concentration = (
                    self.gdl_diffusion.solution_array[self.id_fuel, -1, :, :])
                reduced_shape = mole_flux[self.id_fuel].shape
                # Average along z-axis according to the reduced discretization
                filter_size = (
                    int(fuel_cl_concentration.shape[-1] / reduced_shape[-1]))
                fuel_cl_concentration = ndimage.uniform_filter(
                    fuel_cl_concentration, size=filter_size, axes=(-1,))
                flux_scaling_factors = ndimage.uniform_filter(
                    self.gdl_diffusion.flux_scaling_factors[self.id_fuel],
                    size=filter_size, axes=(-1,))
                diff_coeff_by_length = ndimage.uniform_filter(
                    self.gdl_diffusion.diff_coeff_by_length[self.id_fuel],
                    size=filter_size, axes=(-1,))
                # Rescale to reduced discretization
                fuel_cl_concentration = gf.rescale(fuel_cl_concentration,
                                                   reduced_shape)
                flux_scaling_factors = gf.rescale(flux_scaling_factors,
                                                  reduced_shape)
                diff_coeff_by_length = gf.rescale(diff_coeff_by_length,
                                                  reduced_shape)
                # if gs.global_state.iteration == 20:
                #     print('test')
                fuel_gdl_concentration = np.asarray([
                    ip.interpolate_1d(channel_concentration[self.id_fuel])
                    for i in range(current_density.shape[-1])]).transpose()
                self.electrochemistry.update(
                    current_density, fuel_cl_concentration,
                    reference_fuel_concentration,
                    scaling_factors=flux_scaling_factors,
                    inlet_concentration=reference_fuel_concentration)

                mass_flux_old = np.copy(mass_flux)
                mole_flux_old = np.copy(mole_flux)
                # mole_flux_avg_old = np.average(np.average(mole_flux, axis=-1),
                #                                axis=-1)
                # Recalculate mass flux at GDL/Channel interface
                bc_type = 'Dirichlet'
                mole_flux = - self.gdl_diffusion.calc_boundary_flux(bc_type)
                mole_flux_avg = np.average(np.average(mole_flux, axis=-1),
                                           axis=-1)
                mass_flux = (mole_flux.transpose() *
                             self.channel.fluid.species_mw).transpose()
                axes, indices, _ = self.gdl_diffusion.get_boundary_indices(
                    bc_type)
                flux_interface_area = (
                    self.gdl_diffusion.get_values(
                        self.gdl_diffusion.discretization.d_area[axes[0]],
                        axes,
                        indices))

                # # Recalculate mass flux at GDL/Channel interface
                # bc_type = 'Neumann'
                # mole_flux_cl = - self.gdl_diffusion.calc_boundary_flux(
                #     bc_type)
                # mole_flux_avg_cl = np.average(np.average(mole_flux_cl, axis=-1),
                #                               axis=-1)
                #
                # mass_flux_cl = (mole_flux_cl.transpose() *
                #                 self.channel.fluid.species_mw).transpose()
                # axes, indices, _ = self.gdl_diffusion.get_boundary_indices(
                #     bc_type)
                # flux_interface_area_cl = (
                #     self.gdl_diffusion.get_values(
                #         self.gdl_diffusion.discretization.d_area[axes[0]],
                #         axes,
                #         indices))
                # Factor 2 due to symmetry in gdl_diffusion model
                flux_interface_area *= 2.0

            else:
                fuel_gdl_concentration = np.asarray([
                    ip.interpolate_1d(channel_concentration[self.id_fuel])
                    for i in range(current_density.shape[-1])]).transpose()
                self.electrochemistry.update(
                    current_density, fuel_gdl_concentration,
                    reference_fuel_concentration,
                    inlet_concentration=reference_fuel_concentration)
                flux_interface_area = self.discretization.d_area

            if update_channel:
                # Calculate mole and mass source from fluxes
                mass_source = self.surface_flux_to_channel_source(
                    -mass_flux, area=flux_interface_area)
                mass_source_old = self.surface_flux_to_channel_source(
                    mass_flux_old)
                mole_source = self.surface_flux_to_channel_source(
                    -mole_flux, flux_interface_area)
                self.channel.mass_source[:], self.channel.mole_source[:] = (
                    mass_source, mole_source)
                self.channel.update(update_mass=True, update_flow=False,
                                    update_heat=False, update_fluid=True)
            # if gs.global_state.iteration == 50:
            #     print('test')
            self.update_voltage_loss(current_density)

            # Calculate stoichiometry
            current = self.surface_flux_to_channel_source(current_density)
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

    def calc_species_flux(self, current_density: np.ndarray):
        if not isinstance(current_density, np.ndarray):
            current_density = np.asarray(current_density)
        mole_flux = np.zeros((self.channel.fluid.n_species,
                             *current_density.shape))

        for i in range(len(mole_flux)):
            mole_flux[i] = (current_density * self.n_stoi[i]
                            / (self.n_charge * self.faraday))

        # Water cross flow
        # water_cross_source = self.surface_flux_to_channel_source(
        #     self.w_cross_flow)
        mole_flux[self.id_h2o] += self.w_cross_flow
        # self.channel.flow_direction
        mass_flux = (mole_flux.transpose()
                     * self.channel.fluid.species_mw).transpose()
        return mass_flux, mole_flux

    def surface_flux_to_channel_source(self, flux: np.ndarray, area=None):
        if area is None:
            area = self.discretization.d_area
        if np.isscalar(flux) or flux.ndim in (0, 1):
            return np.sum(area, axis=0) * flux
        elif flux.ndim == 2:
            return np.sum(area * flux, axis=1)
        elif flux.ndim == 3 and flux.shape[0] == self.channel.fluid.n_species:
            return np.asarray(
                [self.surface_flux_to_channel_source(flux[i], area)
                 for i in range(self.channel.fluid.n_species)])
        else:
            raise ValueError('flux variable must be one- or '
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
