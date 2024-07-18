import numpy as np
import json
import os
from porous_two_phase_flow import porous_layer as pl

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
                 discretization: dsct.Discretization,
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
        # Create evaporation model
        self.evaporation_model = evap.EvaporationModel(
            input_dict['evaporation_model'], fluid)

        # Setup variables
        self.saturation = np.ones(self.transport.base_shape) * 1e-2
        self.capillary_pressure = np.zeros(self.transport.base_shape)

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
        if update_fluid:
            self.fluid.update(temperature, gas_pressure, molar_composition)

        # Boundary conditions for liquid pressure
        sigma_liquid = \
            self.fluid.phase_change_species.calc_surface_tension(temperature)[0]
        humidity = self.fluid.humidity
        capillary_pressure_boundary = \
            self.porosity_model.saturation_model.calc_capillary_pressure(
                liquid_saturation_value, surface_tension=sigma_liquid,
                humidity=humidity)
        liquid_pressure_boundary = capillary_pressure_boundary + gas_pressure

        # Calculate constant factor of transport property from absolute
        # permeability
        absolute_permeability = np.asarray(
            self.porosity_model.dict['permeability'])
        dens_by_visc = self.fluid.liquid.density / self.fluid.liquid.viscosity
        transport_property = [value * dens_by_visc for value in
                              absolute_permeability]

        capillary_pressure = np.copy(gas_pressure)

        # Calculate relative permeability from saturation
        relative_permeability = self.porosity_model.calc_relative_permeability(
            self.saturation)


        # Calculate source terms from evaporation/condensation
        interfacial_area = (
            self.porosity_model.calc_two_phase_interfacial_area(self.saturation))
        # evaporation_rate = evap_model.calc_evaporation_rate(
        #     temperature=t, pressure=p, capillary_pressure=p_cap)
        evaporation_rate = self.evaporation_model.calc_evaporation_rate(
            saturation=self.saturation, temperature=temperature,
            pressure=gas_pressure,
            capillary_pressure=ca, porosity=self.porosity_model.porosity)
        specific_area = (interfacial_area / self.transport.transport_layers[0]
                         .discretization.d_volume)
        specific_area = 1.0
        volumetric_evap_rate = specific_area * evaporation_rate
        self.transport.update(
            liquid_pressure_boundary, liquid_mass_flux, transport_property)
        # Create liquid pressure diffusion transport system

