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
            ls.BasicLinearSystem.create(item, self.transport_type,
                                        previously_shifted_axis=item.shift_axis)
            for item in self.transport_layers]
        # for lin_sys in self.linear_systems:
        #     lin_sys.sparse_solve = False
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
        self.discretization = self.transport_layers[0].discretization

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
            input_dict['effective'] = False
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

    def get_values(self, values: (list, np.ndarray), axes: tuple,
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

    def get_boundary_indices(self, boundary_type: str):
        if boundary_type == 'Neumann':
            axes = self.neumann_bc.axes
            indices = self.neumann_bc.indices
        elif boundary_type == 'Dirichlet':
            axes = self.dirichlet_bc.axes
            indices = self.dirichlet_bc.indices
        else:
            raise ValueError('"boundary_type" must be "Neumann" or "Dirichlet"')
        index_adjacent = 1 if indices[0] == 0 else -2
        indices_adjacent = list(indices)
        indices_adjacent[0] = index_adjacent
        indices_adjacent = tuple(indices_adjacent)
        return axes, indices, indices_adjacent

    @abstractmethod
    def calc_boundary_flux(self, boundary_type: str):
        pass

    def calc_flux_scaling_factors(
            self, solution_values: (list, np.ndarray),
            dirichlet_bc_values: (list, np.ndarray),
            neumann_bc_flux: (list, np.ndarray)):
        axes, indices = self.neumann_bc.axes, self.neumann_bc.indices

        flux_boundary_concentrations = self.get_values(
            solution_values, axes, indices)
        boundary_concentration = self.get_values(
            dirichlet_bc_values, axes, indices)
        flux_boundary_flux = self.get_values(
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

    def calc_boundary_flux(self, boundary_type: str):
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
               transport_property: (list, np.ndarray),
               *args, **kwargs):
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
                                    explicit_sources=volumetric_source)
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

    def calc_boundary_flux(self, boundary_type: str):
        raise NotImplementedError


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

        self.ids_active = tuple(i for i in range(fluid.n_species)
                                if i != self.id_inert)
        self.initialize = True
        self.fluid = fluid.copy(self.base_shape, plot_axis=-2)
        if self.id_inert is not None:
            self.solution_array = np.zeros((len(self.linear_systems) + 1,
                                            *self.base_shape))
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
               dirichlet_molar_composition: np.ndarray,
               neumann_flux: np.ndarray,
               explicit_source: np.ndarray = None,
               implicit_source: np.ndarray = None,
               volume_fraction: np.ndarray = None,
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
        dirichlet_molar_composition = self.rescale_input(
            dirichlet_molar_composition, except_first_axis=True)
        neumann_flux = self.rescale_input(neumann_flux, except_first_axis=True)
        if explicit_source is not None:
            explicit_source = self.rescale_input(explicit_source,
                                                 except_first_axis=True)
            explicit_source = [explicit_source[i] for i in range(len(
                explicit_source)) if i != self.id_inert]
        if implicit_source is not None:
            implicit_source = self.rescale_input(implicit_source,
                                                 except_first_axis=True)
            implicit_source = [implicit_source[i] for i in range(len(
                implicit_source)) if i != self.id_inert]
        if self.initialize:
            self.fluid.update(temperature, pressure,
                              dirichlet_molar_composition)

        # Update species transport for active species
        boundary_composition = [
            dirichlet_molar_composition[i] for i in
            range(len(dirichlet_molar_composition)) if i != self.id_inert]
        boundary_flux = [neumann_flux[i] for i in range(len(neumann_flux))
                         if i != self.id_inert]
        self.gas_diff_coeff.update()

        for i, trans_layer in enumerate(self.transport_layers):
            trans_layer.update(
                {self.transport_type:
                    [self.gas_diff_coeff.d_eff[self.ids_active[i]],
                     self.gas_diff_coeff.d_eff[self.ids_active[i]] * 0.0,
                     self.gas_diff_coeff.d_eff[self.ids_active[i]] * 0.0]},
                volume_fraction=volume_fraction)
            # trans_layer.conductance['diffusion'][:, 0:5, :, 0:10] = 1e-16

        for i, lin_sys in enumerate(self.linear_systems):
            boundary_concentration = mf.get_axis_values(
                    boundary_composition[i],
                    axes=self.dirichlet_bc.axes,
                    indices=self.dirichlet_bc.indices)
            self.neumann_bc.values = boundary_flux[i]
            self.dirichlet_bc.values = boundary_concentration
            if explicit_source is not None:
                explicit_source_data = ls.SourceData(values=explicit_source[i])
            else:
                explicit_source_data = None
            if implicit_source is not None:
                implicit_source_data = ls.SourceData(values=implicit_source[i])
            else:
                implicit_source_data = None
            rhs_input = ls.RHSInput(neumann_bc=self.neumann_bc,
                                    dirichlet_bc=self.dirichlet_bc,
                                    explicit_sources=explicit_source_data,
                                    implicit_sources=implicit_source_data)
            lin_sys.update(rhs_input=rhs_input)

        concentrations = [lin_sys.solution_array for
                          lin_sys in self.linear_systems]
        conc_array = np.array(concentrations)
        show_concentrations = np.moveaxis(conc_array,
                                          (0, 1, 2, 3), (0, 2, 1, 3))
        self.flux_scaling_factors, self.diff_coeff_by_length = (
            self.calc_flux_scaling_factors(
                concentrations, boundary_composition, boundary_flux))
        # concentrations = [np.where(conc < 0.0, 0.0, conc) for conc in
        #                   concentrations]
        concentrations = self.calculate_inert_concentration(concentrations)
        # concentrations = np.asarray([np.where(conc < 0.0, 0.0, conc)
        #                              for conc in concentrations])
        self.solution_array = concentrations
        show_concentrations = np.moveaxis(concentrations,
                                          (0, 1, 2, 3), (0, 2, 1, 3))
        conductances = [lin_sys.conductance for lin_sys in self.linear_systems]
        # if isinstance(self.fluid, (fl.GasMixture, fl.CanteraGasMixture)):
        #     self.fluid.update(temperature, pressure,
        #                       mole_composition=self.solution_array)
        # elif isinstance(self.fluid, (fl.TwoPhaseMixture,
        #                 fl.CanteraTwoPhaseMixture)):
        #     self.fluid.update(temperature, pressure,
        #                       gas_mole_composition=self.solution_array)
        diff_conc = np.diff(conc_array, axis=1)
        self.initialize = False

    def calc_boundary_flux(self, boundary_type: str):

        # Calculate boundary flux: defined positive if going into domain,
        # negative if going out of domain
        axes, indices_boundary, indices_adjacent = self.get_boundary_indices(
            boundary_type)
        boundary_values = (self.get_values(self.solution_array, axes,
                                           indices_boundary))
        adjacent_values = (self.get_values(self.solution_array, axes,
                                           indices_adjacent))
        # Be careful to use the shifted node d_volume from one of the
        # TransportLayer attributes instead of the original discrete volumes
        # from the Discretization attributes
        trans_layer = self.transport_layers[0]
        dx = self.get_values(trans_layer.dx[axes[0]], axes, indices_adjacent)
        # d_volume = self.get_values(
        #     trans_layer.d_volume, axes, indices_boundary)
        # d_area = self.get_values(
        #     self.discretization.d_area, axes, indices_boundary)
        # d_area = d_area[axes[0]]
        # area_sum = np.sum(d_area)
        diff_coeffs_boundary = self.get_diff_coeffs(axes, indices_boundary)
        # diff_coeffs_adjacent = self.get_diff_coeffs(axes, indices_adjacent)
        diff_coeffs = np.stack([self.get_diff_coeffs(axes, indices_boundary),
                                self.get_diff_coeffs(axes, indices_adjacent)],
                               axis=diff_coeffs_boundary.ndim)
        diff_coeffs = np.average(diff_coeffs, axis=-1)

        diff_coeff_by_length = diff_coeffs / dx
        diff_conc = boundary_values - adjacent_values
        # diff_conc_avg = np.average(np.average(diff_conc, axis=-1),
        #                            axis=-1)
        # diff_coeff_avg = np.average(np.average(diff_coeffs, axis=-1),
        #                             axis=-1)
        boundary_flux = - diff_coeff_by_length * (boundary_values -
                                                  adjacent_values)
        return boundary_flux

    def get_diff_coeffs(self, axes, indices):
        diff_coeffs = self.get_values(self.gas_diff_coeff.d_eff, axes, indices)
        trans_layer = self.transport_layers[0]
        if trans_layer.effective:
            diff_coeffs *= (trans_layer.volume_fraction
                            ** trans_layer.bruggeman_exponent)
        return diff_coeffs

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
