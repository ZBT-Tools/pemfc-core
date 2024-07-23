import numpy as np
import json
import os
from porous_two_phase_flow import porous_layer as pl
from porous_two_phase_flow import helper_functions as hf

from .fluid import (fluid as fl, diffusion_coefficient as dc,
                    evaporation_model as evap)
from . import transport_layer as tl
from . import discretization as dsct
from . import linear_system as ls
from . import global_functions as gf
from . import matrix_functions as mf
from . import global_state as gs
from . import diffusion_transport as dt


class TwoPhaseMixtureDiffusionTransport:
    """
    Class to describe the diffusion in a porous layer (here defined through the
    SolidLayer class
    """
    def __init__(self, input_dict: dict,
                 discretization: dsct.Discretization3D,
                 fluid: (fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 *args, **kwargs):
        if not isinstance(fluid, fl.TwoPhaseMixture):
            raise TypeError('fluid must be of type TwoPhaseMixture')

        with open(os.path.join(gs.global_state.base_directory, 'settings',
                               'two_phase_settings.json')) as file:
            input_dict.update(json.load(file))
        input_dict['volume_fraction'] = (
            input_dict['porosity_model']['porosity'])
        input_dict['effective'] = True
        input_dict['permeability'] = (
            input_dict['porosity_model']['permeability'])
        self.dict = input_dict
        self.transport = dt.DiffusionTransport.create(
            input_dict, discretization)

        # self.fluid = fluid.copy(self.liquid_transport.base_shape, plot_axis=-2)
        if fluid.array_shape != self.transport.base_shape:
            self.fluid = fluid.copy(self.transport.base_shape, plot_axis=-2)
            self.fluid_reference = False
        else:
            self.fluid = fluid
            self.fluid_reference = True

        self.dict = input_dict
        # Create porosity model
        self.porosity_model = pl.PorousLayer(self.dict['porosity_model'],
                                             fluid)
        if not isinstance(self.porosity_model, pl.PorousTwoPhaseLayer):
            raise TypeError('attribute "porosity_model" must be of type '
                            'PorousTwoPhaseLayer')
        self.saturation_model = self.porosity_model.saturation_model
        # Create evaporation model
        self.evaporation_model = evap.EvaporationModel(
            input_dict['evaporation_model'], fluid)

        # Setup variables
        self.saturation_min = 1e-6
        self.saturation = (np.ones(self.transport.base_shape) *
                           self.saturation_min)
        self.capillary_pressure = np.zeros(self.transport.base_shape)
        self.net_evaporation_rate = np.zeros(self.transport.base_shape)
        self.condensation_rate = np.zeros(self.transport.base_shape)
        self.evaporation_rate = np.zeros(self.transport.base_shape)
        self.implicit_condensation_coeff = np.zeros(
            self.transport.base_shape)

        self.specific_liquid_surface_area = np.zeros(self.transport.base_shape)

        # Inner looping numerical settings
        self.error = np.inf
        self.max_iterations = 1
        self.min_iterations = 1
        self.urf = 0.5
        self.error_tolerance = 1e-5

    def update(self, temperature: np.ndarray,
               gas_pressure: np.ndarray,
               molar_composition: np.ndarray,
               liquid_saturation_value: np.ndarray,
               liquid_mass_flux: np.ndarray,
               update_fluid: bool = True,
               *args, **kwargs):
        # Rescale input data
        if temperature.shape != self.transport.base_shape:
            temperature = self.transport.rescale_input(temperature)
        if gas_pressure.shape != self.transport.base_shape:
            gas_pressure = self.transport.rescale_input(gas_pressure)
        if molar_composition[0].shape != self.transport.base_shape:
            molar_composition = self.transport.rescale_input(
                molar_composition, except_first_axis=True)
        if liquid_saturation_value.shape != self.transport.base_shape:
            liquid_saturation_value = self.transport.rescale_input(
                liquid_saturation_value)
        if liquid_mass_flux.shape != self.transport.base_shape:
            liquid_mass_flux = self.transport.rescale_input(liquid_mass_flux)
        # if update_fluid:
        #     self.fluid.update(temperature, gas_pressure,
        #                       gas_mole_composition=molar_composition)

        # Boundary conditions for liquid pressure
        sigma_liquid = \
            self.fluid.phase_change_species.calc_surface_tension(temperature)[0]
        humidity = self.fluid.humidity
        # humidity[humidity < 0.5] = 0.5
        # humidity[humidity > 1.2] = 1.2
        show_humidity = np.moveaxis(humidity, (0, 1, 2), (1, 0, 2))
        capillary_pressure_boundary = \
            self.saturation_model.calc_capillary_pressure(
                liquid_saturation_value, surface_tension=sigma_liquid,
                humidity=humidity)
        liquid_pressure_boundary = capillary_pressure_boundary + gas_pressure

        # Calculate constant factor of transport property from absolute
        # permeability
        absolute_permeability = np.asarray(
            self.porosity_model.dict['permeability'])
        dens_by_visc = self.fluid.liquid.density / self.fluid.liquid.viscosity
        constant_transport_property = [value * dens_by_visc for value in
                                       absolute_permeability]

        # capillary_pressure = np.copy(gas_pressure)
        iteration = 0

        # Setup variable containers for previous iteration
        saturation_old = np.copy(self.saturation)
        capillary_pressure_old = np.copy(self.capillary_pressure)
        capillary_pressure = np.copy(self.capillary_pressure)
        vol_net_evap_rate = np.zeros(self.net_evaporation_rate.shape)
        saturation = np.copy(self.saturation)
        vol_cond_rate = np.zeros(self.condensation_rate.shape)
        vol_evap_rate = np.zeros(self.evaporation_rate.shape)
        vol_cond_coeff = np.zeros(self.implicit_condensation_coeff.shape)

        specific_area = np.zeros(self.specific_liquid_surface_area.shape)
        while (all((self.error > self.error_tolerance,
                    iteration < self.max_iterations)) or
               iteration < self.min_iterations):
            # Calculate relative permeability from saturation
            relative_permeability = (
                self.porosity_model.calc_relative_permeability(saturation))
            transport_property = (constant_transport_property
                                  * relative_permeability)

            # Calculate source terms from evaporation/condensation
            specific_area = (
                self.porosity_model.calc_specific_interfacial_area(saturation)
                * 1.0)
            # specific_area = 1000.0
            # evaporation_rate = evap_model.calc_evaporation_rate(
            #     temperature=t, pressure=p, capillary_pressure=p_cap)
            net_evap_rate, evap_rate, cond_rate, cond_coeffs = (
                self.evaporation_model.calc_evaporation_rate(
                    saturation=saturation, temperature=temperature,
                    pressure=gas_pressure,
                    capillary_pressure=None,
                    porosity=self.porosity_model.porosity))
            # specific_area = 1.0
            vol_net_evap_rate = specific_area * net_evap_rate
            vol_evap_rate = specific_area * evap_rate
            vol_cond_rate = specific_area * cond_rate
            vol_cond_coeff = specific_area * cond_coeffs
            show_volumetric_evap_rate = np.moveaxis(vol_net_evap_rate,
                                                    (0, 1, 2), (1, 0, 2))

            self.transport.update(
                liquid_pressure_boundary,
                liquid_mass_flux,
                -vol_net_evap_rate,
                transport_property)

            # Save old capillary pressure
            capillary_pressure_old[:] = np.copy(capillary_pressure)

            # Calculate and constrain capillary pressure
            liquid_pressure = self.transport.solution_array
            capillary_pressure = liquid_pressure - gas_pressure
            show_capillary_pressure = np.round(np.moveaxis(capillary_pressure,
                                               (0, 1, 2), (1, 0, 2)), 2)
            min_value = self.saturation_model.calc_capillary_pressure(
                self.saturation_min,
                surface_tension=np.min(self.fluid.surface_tension),
                humidity=np.max(self.fluid.humidity)
            )
            max_value = self.saturation_model.calc_capillary_pressure(
                0.99,
                surface_tension=np.max(self.fluid.surface_tension),
                humidity=np.min(self.fluid.humidity)
            )
            capillary_pressure[capillary_pressure < min_value] = min_value
            capillary_pressure[capillary_pressure > max_value] = max_value
            show_capillary_pressure = np.round(np.moveaxis(capillary_pressure,
                                               (0, 1, 2), (1, 0, 2)), 2)
            # Save old saturation
            saturation_old[:] = np.copy(saturation)

            # Calculate new saturation
            saturation = self.saturation_model.calc_saturation(
                capillary_pressure)
            show_saturation = np.moveaxis(saturation,
                                          (0, 1, 2), (1, 0, 2))
            # Apply underrelaxation
            saturation[:] = (self.urf * saturation
                             + (1.0 - self.urf) * saturation_old)
            # Constrain solution
            # saturation[saturation < self.saturation_min] = self.saturation_min
            # saturation[saturation > 1.0] = 1.0
            show_saturation = np.moveaxis(saturation,
                                          (0, 1, 2), (1, 0, 2))

            # Error calculation
            s_diff = saturation - saturation_old
            s_diff[:] = np.divide(s_diff, saturation, where=saturation != 0.0)
            p_diff = capillary_pressure - capillary_pressure_old
            p_diff[:] = np.divide(p_diff, capillary_pressure,
                                  where=capillary_pressure != 0.0)
            s_diff = s_diff.ravel()
            p_diff = p_diff.ravel()
            error_s = np.dot(s_diff.transpose(), s_diff) / (2.0 * len(s_diff))
            error_p = np.dot(p_diff.transpose(), p_diff) / (2.0 * len(p_diff))
            self.error = error_s + error_p
            iteration += 1

        show_concentration = np.moveaxis(self.fluid.gas.concentration,
                                         (0, 1, 2, 3), (0, 2, 1, 3))
        self.net_evaporation_rate[:] = vol_net_evap_rate
        self.condensation_rate[:] = vol_cond_rate
        self.evaporation_rate[:] = vol_evap_rate
        self.implicit_condensation_coeff[:] = vol_cond_coeff
        self.specific_liquid_surface_area[:] = specific_area
        self.capillary_pressure[:] = capillary_pressure
        self.saturation[:] = saturation

