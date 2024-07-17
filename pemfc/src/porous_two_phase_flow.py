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
        self.initialize = True
        self.fluid = fluid.copy(discretization.shape, plot_axis=-2)
        with open(os.path.join(gs.global_state.base_directory, 'settings',
                               'two_phase_settings.json')) as file:
            input_dict.update(json.load(file))
        input_dict['volume_fraction'] = (
            input_dict)['porosity_model']['porosity']
        input_dict['effective'] = True
        self.dict = input_dict
        # Create porosity model
        self.porosity_model = pl.PorousLayer(self.dict['porosity_model'],
                                             fluid)
        # Create liquid pressure diffusion transport system
        transport_property = (
                fluid.liquid.density / fluid.liquid.viscosity *
                np.asarray(self.dict['porosity_model']['permeability']))
        transport_layers = [tl.TransportLayer.create(
                input_dict, {'diffusion': transport_property}, discretization)]
        self.liquid_transport = dt.DiffusionTransport.create(input_dict, )


    def update(self, temperature: np.ndarray, pressure: np.ndarray,
               molar_composition: np.ndarray, flux: np.ndarray,
               *args, **kwargs):
        temperature = gf.rescale(temperature, self.fluid.temperature.shape)
        pressure = self.reshape_1d_input(pressure)
        pressure = gf.rescale(pressure, self.fluid.pressure.shape)

        molar_composition = np.asarray(
            [self.reshape_input(mol_comp) for mol_comp in molar_composition])
        molar_composition = gf.rescale(
            molar_composition, self.fluid.gas.mole_fraction.shape)

        flux = np.asarray([self.reshape_input(item) for item in flux])
        flux = gf.rescale(flux, self.fluid.gas.mole_fraction.shape)
        if self.initialize:
            self.fluid.update(temperature, pressure, molar_composition)
        else:
            self.fluid.update(temperature, pressure)

        # Update evaporation model
        self.evaporation_model.update(temperature, pressure)
        # Update species transport for active species
        boundary_composition = [
            mol_comp for i, mol_comp in enumerate(molar_composition)
            if i != self.id_inert]
        boundary_flux = [flux[i] for i in range(len(flux)) if i !=
                         self.id_inert]
        self.gas_diff_coeff.update()
        for i, trans_layer in enumerate(self.transport_layers):
            trans_layer.update(
                {self.transport_type: [self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i],
                                       self.gas_diff_coeff.d_eff[i]]})
            # trans_layer.conductance['diffusion'][:, 0:5, :, 0:10] = 1e-16

        for i, lin_sys in enumerate(self.linear_systems):
            boundary_concentration = mf.get_axis_values(
                boundary_composition[i],
                axes=self.dirichlet_bc.axes,
                indices=self.dirichlet_bc.indices)
            self.neumann_bc.values = boundary_flux[i]
            self.dirichlet_bc.values = boundary_concentration
            rhs_input = ls.RHSInput(neumann_bc=self.neumann_bc,
                                    dirichlet_bc=self.dirichlet_bc)
            lin_sys.update(rhs_values=rhs_input)

        concentrations = [lin_sys.solution_array for
                          lin_sys in self.linear_systems]
        conc_array = np.array(concentrations)
        show_concentrations = np.round(np.moveaxis(conc_array,
                                                   (0, 1, 2, 3), (0, 2, 1, 3)),
                                       2)
        self.flux_scaling_factors[:], self.diff_coeff_by_length[:] = (
            self.calc_flux_scaling_factors(
                concentrations, boundary_composition, boundary_flux))
        concentrations = [np.where(conc < 0.0, 0.0, conc) for conc in
                          concentrations]
        concentrations = self.calculate_inert_concentration(concentrations)
        self.solution_array[:] = concentrations
        show_concentrations = np.round(np.moveaxis(concentrations,
                                                   (0, 1, 2, 3), (0, 2, 1, 3)),
                                       2)
        conductances = [lin_sys.conductance for lin_sys in self.linear_systems]
        # show_conductances = np.round(np.)
        self.initialize = False

    def calculate_inert_concentration(self, concentrations: list) -> np.ndarray:
        # Update species transport for inert species
        if isinstance(self.fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            gas_constant = self.fluid.gas_constant
        else:
            gas_constant = self.fluid.gas.gas_constant
        total_gas_concentration = (
                self.fluid.pressure / (gas_constant * self.fluid.temperature))
        active_concentrations = np.sum(concentrations, axis=0)
        total_gas_concentration = gf.rescale(total_gas_concentration,
                                             active_concentrations.shape)
        inert_concentration = total_gas_concentration - active_concentrations
        concentrations.insert(self.id_inert, inert_concentration)
        return np.asarray(concentrations)