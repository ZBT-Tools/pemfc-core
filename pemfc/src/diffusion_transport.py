from abc import ABC, abstractmethod
import numpy as np
import json
import os

from .fluid import (fluid as fl, diffusion_coefficient as dc,
                    evaporation_model as evap)
from . import transport_layer as tl
from . import discretization as dsct
from . import linear_system as ls
from . import global_functions as gf
from . import matrix_functions as mf
from . import global_state as gs
from . import porous_two_phase_flow as p2pf
from porous_two_phase_flow import porous_layer as pl


class DiffusionTransport(ABC):
    """
    Abstract base class to describe the calculate the diffusion transport in
    a solid or porous layer (here defined through the TransportLayer class).
    """
    def __init__(self, input_dict: dict,
                 transport_layers: list[tl.TransportLayer],
                 *args, **kwargs):
        if isinstance(transport_layers, tl.TransportLayer):
            transport_layers = [transport_layers]
        self.dict = input_dict
        self.boundary_patches = input_dict['boundary_patches']
        self.transport_layers = transport_layers
        self.transport_type = 'diffusion'

        self.linear_systems = [
            ls.BasicLinearSystem.create(item, self.transport_type)
            for item in self.transport_layers]
        self.base_shape = self.linear_systems[0].solution_array.shape
        if len(self.linear_systems) == 1:
            solution_shape = self.base_shape
        else:
            solution_shape = (len(self.linear_systems), *self.base_shape)
        self.solution_array = np.zeros(solution_shape)
        self.neumann_bc = self.create_boundary_condition('Neumann')
        self.dirichlet_bc = self.create_boundary_condition('Dirichlet')
        self.flux_scaling_factors = np.zeros(self.neumann_bc.values.shape)
        self.diff_coeff_by_length = np.zeros(self.neumann_bc.values.shape)

    @classmethod
    def create(cls, input_dict: dict,
               discretization: dsct.Discretization,
               fluid: (fl.GasMixture, fl.CanteraTwoPhaseMixture,
                       fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture) = None,
               id_inert: int = None, *args,
               **kwargs):
        model_type = input_dict['type']
        input_dict['shift_axis'] = 0
        if model_type == 'Constant':
            input_dict['effective'] = False
            transport_layers = [
                tl.TransportLayer.create(
                    input_dict,
                    {'diffusion': input_dict['diffusion_coefficient']},
                    discretization)]

            return ConstantDiffusionTransport(
                input_dict, transport_layers, **kwargs)
        elif model_type == 'GasMixture' and fluid is not None:
            if id_inert is None:
                raise ValueError('for type GasMixtureDiffusionTransport the '
                                 'integer value "id_inert" must be specified')
            input_dict['volume_fraction'] = input_dict['porosity']
            input_dict['effective'] = True
            species_names = fluid.species_names
            transport_layers = [
                tl.TransportLayer.create(
                    input_dict, {'diffusion': 0.0}, discretization)
                for i in range(len(species_names)) if i != id_inert]

            return GasMixtureDiffusionTransport(
                input_dict, transport_layers, fluid, id_inert=id_inert,
                **kwargs)
        elif model_type == 'DarcyFlow':
            input_dict['volume_fraction'] = input_dict['porosity']
            input_dict['effective'] = True
            permeability = input_dict.get('permeability', 0.0)
            transport_layers = [
                tl.TransportLayer.create(
                    input_dict, {'diffusion': permeability}, discretization)]
            return DarcyFlowDiffusionTransport(
                input_dict, transport_layers, **kwargs)
        # elif (model_type == 'TwoPhaseMixture'
        #       and isinstance(fluid, fl.TwoPhaseMixture)):
        #     input_dict['volume_fraction'] = input_dict['porosity']
        #     input_dict['effective'] = True
        #     permeability = input_dict.get('permeability', 0.0)
        #     transport_layers = [
        #         tl.TransportLayer.create(
        #             input_dict, {'diffusion': permeability}, discretization)]
        #     return p2pf.TwoPhaseMixtureDiffusionTransport(
        #         input_dict, transport_layers, fluid, **kwargs)
        else:
            raise NotImplementedError(
                'given type {} not implemented as concrete '
                'implementation of abstract class DiffusionTransport'.format(
                    model_type))

    def create_boundary_condition(self, name: str):
        axes = self.boundary_patches[name]['axes']
        indices = self.boundary_patches[name]['indices']
        shape = self.calc_bc_array_shape(axes, indices)
        return ls.BoundaryData(values=np.zeros(shape), axes=axes,
                               indices=indices)

    @abstractmethod
    def update(self, *args, **kwargs):
        pass

    @staticmethod
    def reshape_1d_input(array: np.ndarray):
        return np.moveaxis(np.asarray([[array,],]), (0, 1, 2), (0, 2, 1))

    @staticmethod
    def reshape_2d_input(array: np.ndarray):
        return np.moveaxis(np.asarray([array,]), (0, 1, 2), (0, 1, 2))

    @classmethod
    def reshape_input(cls, array: np.ndarray, except_first_axis=False):
        array = np.asarray(array)
        if except_first_axis:
            return np.asarray([cls.reshape_input(sub_array)
                               for sub_array in array])
        else:
            if array.ndim == 1:
                return cls.reshape_1d_input(array)
            elif array.ndim == 2:
                return cls.reshape_2d_input(array)
            elif array.ndim == 3:
                return array
            else:
                raise ValueError('only up to three-dimensional arrays allowed')

    def rescale_input(self, array: np.ndarray, except_first_axis=False,
                      shape=None):
        if shape is None:
            shape = self.base_shape
        reshaped_array = self.reshape_input(
            array, except_first_axis=except_first_axis)
        if except_first_axis:
            return np.asarray([gf.rescale(arr, shape)
                               for arr in reshaped_array])
        else:
            return gf.rescale(reshaped_array, shape)

    def calc_bc_array_shape(self, axes: tuple, indices: tuple):
        if len(self.linear_systems) > 1:
            return np.asarray([mf.get_axis_values(item, axes, indices)
                               for item in self.solution_array]).shape
        else:
            return mf.get_axis_values(self.solution_array, axes, indices).shape

    def get_boundary_values(self, values: (list, np.ndarray), axes: tuple,
                            indices: tuple) \
            -> np.ndarray:
        if isinstance(values, list):
            return np.asarray([mf.get_axis_values(item, axes, indices)
                               for item in values])
        elif isinstance(values, np.ndarray):
            if values.ndim == len(self.base_shape) + 1:
                return np.asarray([mf.get_axis_values(item, axes, indices)
                                   for item in values])
            elif values.ndim == len(self.base_shape):
                return mf.get_axis_values(values, axes, indices)
            else:
                raise ValueError('values must be of same shape as '
                                 'attribute solution_array')
        else:
            raise TypeError('values must be of types (list, np.ndarray)')

    def calc_boundary_flux(self, boundary_type: str):
        if boundary_type == 'Neumann':
            axes = self.neumann_bc.axes
            indices = self.neumann_bc.indices
        elif boundary_type == 'Dirichlet':
            axes = self.dirichlet_bc.axes
            indices = self.dirichlet_bc.indices
        else:
            raise ValueError('"boundary_type" must be "Neumann" or "Dirichlet"')
        index_next = 1 if axes[0] == 0 else -2
        indices_next = list(indices)
        indices_next[0] = index_next
        indices_next = tuple(indices_next)
        boundary_values = (self.get_boundary_values(self.solution_array, axes,
                           indices))
        adjacent_values = (self.get_boundary_values(self.solution_array, axes,
                           indices_next))
        # Calculate boundary flux: defined positive if going into domain,
        # negative if going out of domain
        d_volume = self.get_boundary_values(
            self.transport_layers[0].d_volume, axes, indices)
        d_area = self.get_boundary_values(
            self.transport_layers[0].discretization.d_area, axes, indices)
        diff_coeffs = [item.transport_properties[self.transport_type][axes[0]]
                       for item in self.transport_layers]
        diff_coeff = self.get_boundary_values(diff_coeffs, axes, indices)
        # TODO: Include inert species diffusion coefficient

        diff_coeff_by_length = diff_coeff * d_area[axes[0]] / d_volume
        boundary_flux = - diff_coeff_by_length * (adjacent_values -
                                                       boundary_values)
        return boundary_flux

    def calc_flux_scaling_factors(
            self, solution_values: (list, np.ndarray),
            dirichlet_bc_values: (list, np.ndarray),
            neumann_bc_flux: (list, np.ndarray)):
        axes, indices = self.neumann_bc.axes, self.neumann_bc.indices

        flux_boundary_concentrations = self.get_boundary_values(
            solution_values, axes, indices)
        boundary_concentration = self.get_boundary_values(
            dirichlet_bc_values, axes, indices)
        flux_boundary_flux = self.get_boundary_values(
            neumann_bc_flux, axes, indices)

        delta_conc = boundary_concentration - flux_boundary_concentrations
        flux_scaling_factors = (delta_conc / boundary_concentration)

        diff_coeff_by_length = np.abs(
            np.divide(flux_boundary_flux, delta_conc,
                      out=np.zeros(delta_conc.shape),
                      where=delta_conc != 0.0))

        return flux_scaling_factors, diff_coeff_by_length


class ConstantDiffusionTransport(DiffusionTransport):
    """
    Abstract base class to describe the calculate the diffusion transport in
    a solid or porous layer (here defined through the TransportLayer class).
    """
    def __init__(self, input_dict: dict,
                 transport_layers: list[tl.TransportLayer], *args, **kwargs):
        super().__init__(input_dict, transport_layers, *args, **kwargs)

    def update(self, *args, **kwargs):
        raise NotImplementedError


class DarcyFlowDiffusionTransport(DiffusionTransport):
    """
    Abstract base class to describe the calculate the diffusion transport in
    a solid or porous layer (here defined through the TransportLayer class).
    """
    def __init__(self, input_dict: dict,
                 transport_layers: list[tl.TransportLayer], *args, **kwargs):
        super().__init__(input_dict, transport_layers, *args, **kwargs)

    def update(self, dirichlet_values: np.ndarray,
               neumann_values: np.ndarray,
               source_values: np.ndarray,
               transport_property: (list, np.ndarray), *args, **kwargs):
        dirichlet_values = self.rescale_input(dirichlet_values)
        neumann_values = self.rescale_input(neumann_values)
        for i, trans_layer in enumerate(self.transport_layers):
            trans_layer.update({self.transport_type: transport_property})
            # trans_layer.conductance['diffusion'][:, 0:5, :, 0:10] = 1e-16

        for i, lin_sys in enumerate(self.linear_systems):
            boundary_values = mf.get_axis_values(
                    dirichlet_values,
                    axes=self.dirichlet_bc.axes,
                    indices=self.dirichlet_bc.indices)
            self.neumann_bc.values = neumann_values
            self.dirichlet_bc.values = boundary_values
            volumetric_source = ls.SourceData(values=source_values)
            rhs_input = ls.RHSInput(neumann_bc=self.neumann_bc,
                                    dirichlet_bc=self.dirichlet_bc,
                                    volumetric_sources=volumetric_source)
            lin_sys.update(rhs_input=rhs_input)

        solution_list = [lin_sys.solution_array for lin_sys
                         in self.linear_systems]
        solution_array = np.asarray(solution_list[0])
        show_solution = np.round(np.moveaxis(solution_array,
                                 (0, 1, 2), (1, 0, 2)), 2)
        self.flux_scaling_factors, self.diff_coeff_by_length = (
            self.calc_flux_scaling_factors(
                solution_array, dirichlet_values, neumann_values))
        self.solution_array[:] = solution_array
        # # self.solution_array[self.solution_array < 0.0] = 0.0
        # show_solution = np.round(np.moveaxis(solution_array,
        #                          (0, 1, 2), (1, 0, 2)), 2)


class GasMixtureDiffusionTransport(DiffusionTransport):
    """
    Abstract base class to describe the calculate the diffusion transport in
    a solid or porous layer (here defined through the TransportLayer class).
    """
    def __init__(self, input_dict: dict,
                 transport_layers: list[tl.TransportLayer],
                 fluid: (fl.GasMixture, fl.CanteraGasMixture,
                         fl.TwoPhaseMixture, fl.CanteraTwoPhaseMixture),
                 id_inert: int, **kwargs):
        super().__init__(input_dict, transport_layers)
        self.id_inert = id_inert
        self.initialize = True
        self.fluid = fluid.copy(self.base_shape, plot_axis=-2)

        if isinstance(fluid, (fl.GasMixture, fl.CanteraGasMixture)):
            self.gas_diff_coeff = (
                dc.MixtureAveragedDiffusionCoefficient(self.fluid))
        elif isinstance(fluid, (fl.TwoPhaseMixture,
                                fl.CanteraTwoPhaseMixture)):
            self.gas_diff_coeff = (
                dc.MixtureAveragedDiffusionCoefficient(self.fluid.gas))
        else:
            raise TypeError(
                'fluid must be of types (GasMixture, CanteraGasMixture, '
                'TwoPhaseMixture, or CanteraTwoPhaseMixture)')

    def update(self, temperature: np.ndarray, pressure: np.ndarray,
               molar_composition: np.ndarray, flux: np.ndarray,
               *args, **kwargs):
        # temperature = gf.rescale(temperature, self.fluid.temperature.shape)
        # pressure = self.reshape_input(pressure)
        # pressure = gf.rescale(pressure, self.fluid.pressure.shape)
        # molar_composition = np.asarray(
        #     [self.reshape_input(mol_comp) for mol_comp in molar_composition])
        # molar_composition = gf.rescale(
        #     molar_composition, self.fluid.gas.mole_fraction.shape)
        # flux = np.asarray([self.reshape_input(item) for item in flux])
        # flux = gf.rescale(flux, self.fluid.gas.mole_fraction.shape)
        temperature = self.rescale_input(temperature)
        pressure = self.rescale_input(pressure)
        molar_composition = self.rescale_input(molar_composition,
                                               except_first_axis=True)
        flux = self.rescale_input(flux, except_first_axis=True)
        if self.initialize:
            self.fluid.update(temperature, pressure, molar_composition)
        else:
            self.fluid.update(temperature, pressure, self.solution_array)

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
            lin_sys.update(rhs_input=rhs_input)

        concentrations = [lin_sys.solution_array for
                          lin_sys in self.linear_systems]
        conc_array = np.array(concentrations)
        show_concentrations = np.round(np.moveaxis(conc_array,
                                       (0, 1, 2, 3), (0, 2, 1, 3)), 2)
        self.flux_scaling_factors, self.diff_coeff_by_length = (
            self.calc_flux_scaling_factors(
                concentrations, boundary_composition, boundary_flux))
        concentrations = [np.where(conc < 0.0, 0.0, conc) for conc in
                          concentrations]
        concentrations = self.calculate_inert_concentration(concentrations)
        self.solution_array = concentrations
        show_concentrations = np.round(np.moveaxis(concentrations,
                                       (0, 1, 2, 3), (0, 2, 1, 3)), 2)
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
