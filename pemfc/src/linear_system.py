# General imports
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING
import math
import numpy as np
# from numpy import linalg
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve
# import pemfc.src.matrix_functions as mtx_func
# Local module imports
from . import (
    stack as stack_module, channel as chl, matrix_functions as mtx_func,
    global_functions as g_func, transport_layer as tl, discretization as dsct)

if TYPE_CHECKING:
    from .stack import Stack
    from .cell import Cell

# import pandas as pd
# from numba import jit

np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class LinearSystem(ABC):

    def __init__(self, shape: tuple[int, ...]):

        # self.cells = stack.cells
        # if not isinstance(stack.cells, (list, tuple)):
        #     raise TypeError
        # if not isinstance(self.cells[0], cell_module.Cell):
        #     raise TypeError
        # self.n_cells = stack.n_cells
        self.shape = shape
        self.n = np.prod(shape)

        # Use SciPy sparse solver, efficient for larger sparse matrices
        self.sparse_solve = True

        # Instead of solving the completely coupled temperature src at once
        # solve the decoupled cell-wise temperature systems and iterate
        # however not working yet!!!
        self.solve_individual_cells = False

        # Conductance Matrix as ndarray
        self.mtx = np.zeros((self.n, self.n))
        # Solution vector
        self.solution_vector = np.zeros(self.n)
        # Right side of the matrix src: mtx * solution_vector = rhs,
        # contains the power sources and explicit coupled terms
        self.rhs = np.zeros(self.solution_vector.shape)

        self.solution_array = np.reshape(
            self.solution_vector, self.shape, order='F')

    def update(self, *args, **kwargs):
        """
        Wrapper to update the overall linear system.
        """
        self.update_and_solve_linear_system()

    def update_and_solve_linear_system(self, *args, **kwargs):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix(*args, **kwargs)
        self.update_rhs(*args, **kwargs)
        self.solve_system()

    @abstractmethod
    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        pass

    @abstractmethod
    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        pass

    def solve_system(self):
        """
        Solves for the solution vector.
        """
        if self.sparse_solve:
            cond_number = np.linalg.cond(self.mtx)
            mtx = sparse.csr_matrix(self.mtx)
            self.solution_vector[:] = spsolve(mtx, self.rhs)
        else:
            self.solution_vector[:] = np.linalg.tensorsolve(self.mtx, self.rhs)
        self.solution_array[:] = np.reshape(
            self.solution_vector, self.shape, order='F')

    @abstractmethod
    def add_explicit_source(self, *args, **kwargs):
        pass

    @abstractmethod
    def add_implicit_source(self, *args, **kwargs):
        pass

    # @staticmethod
    # def shift_nodes(conductance_list: list, axis=0, **kwargs):
    #     conductance_list = [np.asarray(item) for item in conductance_list]
    #     shift_axes = tuple(i for i in range(len(conductance_list)) if i != axis)
    #     for i in shift_axes:
    #         conductance_list[i] = (
    #                 (conductance_list[i]
    #                  + np.roll(conductance_list[i], 1, axis=axis)) * 0.5)
    #         conductance_list[i] = np.concatenate(
    #             (conductance_list[i], [conductance_list[i][0]]), axis=axis)
    #         conductance_list[i][0] *= 0.5
    #         conductance_list[i][-1] *= 0.5
    #     return conductance_list

    @staticmethod
    def shift_nodes(conductance: list, axis=0, **kwargs):
        conductance = [np.asarray(item) for item in conductance]
        if len(conductance) != conductance[0].ndim:
            raise ValueError('conductance values must be given for '
                             'all dimensions')
        shift_axes = tuple(i for i in range(len(conductance)) if i != axis)
        if axis == -1:
            axis = conductance[0].ndim - 1
        for i in shift_axes:
            old_shape = conductance[i].shape
            new_shape = tuple(size + 1 if j == axis else size
                              for j, size in enumerate(old_shape))
            new_conductance = np.zeros(new_shape)
            half_values = 0.5 * np.moveaxis(conductance[i], axis, 0)
            np.moveaxis(new_conductance, axis, 0)[:-1] += half_values
            np.moveaxis(new_conductance, axis, 0)[1:] += half_values
            conductance[i] = new_conductance
        return conductance


class BasicLinearSystem(LinearSystem):
    def __init__(self, transport_layer: tl.TransportLayer, transport_type: str,
                 init_value=0.0):
        self.transport_layer = transport_layer
        self.type = transport_type
        conductance = self.transport_layer.conductance[self.type]
        # Shift conductance nodes along first axis to the outer edges of layer
        self.conductance = self.shift_nodes(conductance, axis=0)
        shape = self.conductance[1].shape
        super().__init__(shape)

        # inter_node_conductance = (
        #     [self.transport_layer.calc_inter_node_conductance(
        #         self.conductance[i], axis=i)
        #      for i in range(len(self.conductance))])
        # self.mtx_const = mtx_func.build_cell_conductance_matrix(
        #     inter_node_conductance)

        # With shifted nodes for axis 0, only the inter-nodal conductance for
        # the remaining axes must be recalculated
        self.mtx_const = mtx_func.build_cell_conductance_matrix(
                [self.conductance[0],
                 tl.TransportLayer.calc_inter_node_conductance(
                    self.conductance[1], axis=1),
                 tl.TransportLayer.calc_inter_node_conductance(
                    self.conductance[2], axis=2)])
        self.mtx_dyn = np.zeros(self.mtx_const.shape)
        self.rhs_const = np.zeros(self.rhs.shape)
        self.rhs_dyn = np.zeros(self.rhs.shape)

        self.index_array = self.create_index_array()

        self.boundary_planes = [self.get_boundary_planes(axis=i) for i in
                                range(len(self.shape))]

        # self.set_neumann_boundary_conditions(200.0, axis=(0, 2),
        #                                      indices=(0, 1))
        # self.set_dirichlet_boundary_conditions(15.0, axis=(0, 2),
        #                                        indices=(0, 0))
        # Initialize solution values
        self.solution_vector[:] = init_value
        self.solution_array[:] = init_value

    @classmethod
    def create(cls, transport_layer: tl.TransportLayer, transport_type: str,
               init_value=0.0):
        if isinstance(transport_layer, tl.TransportLayer2D):
            return BasicLinearSystem2D(transport_layer, transport_type,
                                       init_value)
        elif isinstance(transport_layer, tl.TransportLayer3D):
            return BasicLinearSystem3D(transport_layer, transport_type,
                                       init_value)
        else:
            raise NotImplementedError(
                'argument "transport_layer" must be of '
                'types (TransportLayer2D, TransportLayer3D)')

    @staticmethod
    def add_explicit_source(
            rhs_vector, source_term, index_array=None, replace=False):
        if isinstance(source_term, np.ndarray):
            if rhs_vector.ndim != source_term.ndim:
                raise ValueError('source_term must have same dimensions as '
                                 'rhs_vector')
        if index_array is None:
            if np.isscalar(source_term):
                source_vector = np.full_like(rhs_vector, -source_term)
            else:
                source_vector = np.asarray(-source_term)
        else:
            if replace is True:
                source_vector = np.copy(rhs_vector)
                np.put(source_vector, index_array, -source_term)
            else:
                source_vector = np.zeros(rhs_vector.shape)
                np.put(source_vector, index_array, -source_term)
        if replace is True:
            rhs_vector[:] = source_vector
        else:
            rhs_vector += source_vector
        return rhs_vector, source_vector

    def add_implicit_source(self, matrix, coefficients,
                            index_array=None, replace=False):
        # matrix_size = matrix.shape[0]
        diag_vector = np.diagonal(matrix).copy()
        new_diag_vector, source_vector = self.add_explicit_source(
            diag_vector, -coefficients.flatten(order='F'), index_array,
            replace=replace)
        np.fill_diagonal(matrix, new_diag_vector)
        return matrix, source_vector

    @staticmethod
    def set_implicit_fixed(matrix, index_array):
        index_array = index_array.flatten(order='F')
        row_length = matrix.shape[0]
        for row_id in index_array:
            row = np.zeros(row_length)
            row[row_id] = 1.0
            matrix[row_id] = row
        return matrix

    def set_neumann_boundary_conditions(
            self, flux_value: float, axes: tuple, indices: tuple):

        index_array = mtx_func.get_axis_values(self.index_array, axes, indices)
        d_area = mtx_func.get_axis_values(
            self.transport_layer.discretization.d_area[axes[0]], axes, indices)
        source = flux_value * d_area
        rhs_vector, _ = self.add_explicit_source(
            self.rhs_const, source.flatten(order='F'),
            index_array.flatten(order='F'), replace=True)

    def set_dirichlet_boundary_conditions(
            self, fixed_value: float, axes: tuple, indices: tuple):
        index_array = mtx_func.get_axis_values(self.index_array, axes, indices)
        index_vector = index_array.flatten(order='F')
        self.set_implicit_fixed(self.mtx_const, index_vector)
        self.add_explicit_source(self.rhs_const, -fixed_value, index_vector,
                                 replace=True)

    def get_boundary_planes(self, axis):
        first_plane = np.moveaxis(self.index_array, axis, 0)[0]
        last_plane = np.moveaxis(self.index_array, axis, 0)[-1]
        return first_plane, last_plane

    def create_index_array(self):
        """
        Create array with shape of linear systems shape containing the ids of
        the corresponding entries in the flattened solution vector or the
        symmetric coefficient matrix, respectively.
        """
        ids_flat = np.arange(self.solution_vector.shape[0])
        return np.reshape(ids_flat, self.solution_array.shape, order='F')

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        self.rhs[:] = self.rhs_const + self.rhs_dyn

    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        self.mtx[:] = self.mtx_const + self.mtx_dyn


class BasicLinearSystem2D(BasicLinearSystem):
    def __init__(self, transport_layer: tl.TransportLayer, transport_type: str,
                 init_value=0.0):

        super().__init__(transport_layer, transport_type, init_value)


class BasicLinearSystem3D(BasicLinearSystem):
    def __init__(self, transport_layer: tl.TransportLayer, transport_type: str,
                 init_value=0.0):
        super().__init__(transport_layer, transport_type, init_value)


class StackedLayerLinearSystem(LinearSystem):

    def __init__(self, layers: list[tl.TransportLayer2D], transport_type: str,
                 init_value=0.0):
        self.conductance = self.calculate_conductance(layers, transport_type)
        shape = self.conductance[1].shape
        super().__init__(shape)
        self.type = transport_type
        self.layers = layers

        self.mtx_const = mtx_func.build_cell_conductance_matrix(
                [self.conductance[0],
                 tl.TransportLayer2D.calc_inter_node_conductance(
                    self.conductance[1], axis=1),
                 tl.TransportLayer2D.calc_inter_node_conductance(
                    self.conductance[2], axis=2)])

        self.mtx_dyn = np.zeros(self.mtx_const.shape)
        self.rhs_const = np.zeros(self.rhs.shape)
        self.rhs_dyn = np.zeros(self.rhs.shape)
        self.index_array = self.create_cell_index_list(self.shape)

        # Initialize solution values
        self.solution_vector[:] = init_value
        self.solution_array[:] = init_value

    def calculate_conductance(self, layers: list[tl.TransportLayer2D],
                              transport_type: str):
        if transport_type not in ('electrical', 'thermal', 'diffusion'):
            raise ValueError("transport_type argument must be either "
                             "'electrical', 'thermal' or 'diffusion'")

        # Stack thermal conductances along through-plane direction,
        # i.e. x-coordinate
        dims = len(layers[0].conductance[transport_type])
        conductance = [
            [layer.conductance[transport_type][i] for layer in layers]
            for i in range(dims)]
        return self.stack_conductance_layers(conductance)

    def stack_conductance_layers(self, conductance_list: list, axis: int = 0):
        return super().shift_nodes(conductance_list, axis=axis)

    @staticmethod
    def add_explicit_source(rhs_vector, source_term, index_array,
                            layer_id=None, replace=False):
        if isinstance(source_term, np.ndarray):
            if rhs_vector.ndim != source_term.ndim:
                raise ValueError('source_term must have same dimensions as '
                                 'rhs_vector')
        if layer_id is None:
            if np.isscalar(source_term):
                source_vector = np.full_like(rhs_vector, -source_term)
            else:
                source_vector = np.asarray(-source_term)
        else:
            if replace is True:
                source_vector = np.copy(rhs_vector)
                np.put(source_vector, index_array[layer_id], -source_term)
            else:
                source_vector = np.zeros(rhs_vector.shape)
                np.put(source_vector, index_array[layer_id], -source_term)
        if replace is True:
            rhs_vector[:] = source_vector
        else:
            rhs_vector += source_vector
        return rhs_vector, source_vector

    def add_implicit_source(self, matrix, coefficients, index_array,
                            layer_id=None, replace=False):
        # matrix_size = matrix.shape[0]
        diag_vector = np.diagonal(matrix).copy()
        new_diag_vector, source_vector = self.add_explicit_source(
            diag_vector, -coefficients.flatten(order='F'), index_array,
            layer_id=layer_id, replace=replace)
        np.fill_diagonal(matrix, new_diag_vector)
        return matrix, source_vector

    @staticmethod
    def set_implicit_layer_fixed(matrix, index_array, layer_id):
        row_length = matrix.shape[0]
        for row_id in index_array[layer_id]:
            row = np.zeros(row_length)
            row[row_id] = 1.0
            matrix[row_id] = row
        return matrix

    def set_neumann_boundary_conditions(self, flux_value, layer_id):
        source = flux_value * self.layers[layer_id].discretization.d_area
        self.add_explicit_source(
            self.rhs_const, source.flatten(order='F'),
            self.index_array, layer_id)

    def set_dirichlet_boundary_conditions(self, fixed_value, layer_id):
        self.set_implicit_layer_fixed(
            self.mtx_const, self.index_array, layer_id)
        self.rhs_const[:], _ = self.add_explicit_source(
            self.rhs_const, -fixed_value, self.index_array, layer_id,
            replace=True)

    @staticmethod
    def create_cell_index_list(shape: tuple[int, ...]):
        """
        Create list of lists with each list containing flattened order of indices
        for a continuous functional layer (for the conductance matrix).
        Functional layers a discretized in x-direction (z-direction in old version),
        the remaining flattened order is equal to the order in the overall matrix
        and the corresponding right-hand side vector
        (typically first y-, then z-direction within a layer a.k.a x-plane)

        Args:
            shape: tuple of size 3 containing the number of layers (index 0),
            the number of y-elements (index 1) and the number z-elements (index 2)
        """
        index_list = []
        for i in range(shape[0]):
            index_list.append(
                [(j * shape[0]) + i for j in range(shape[1] * shape[2])])
        return index_list

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        pass

    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        pass


class CellLinearSystem(StackedLayerLinearSystem):

    def __init__(self, cell: Cell, transport_type: str, init_value=0.0):
        self.cell = cell

        super().__init__(cell.layers, transport_type, init_value)

        if transport_type == 'thermal':
            self.add_ambient_convection()

    def set_layer_boundary_conditions(self, layer_id):
        if 'flux_endplate' in self.cell.cell_dict:
            flux_value = self.cell.cell_dict['flux_endplate']
            self.set_neumann_boundary_conditions(flux_value, layer_id)

        elif 'value_endplate' in self.cell.cell_dict:
            fixed_value = self.cell.cell_dict['value_endplate']
            self.set_dirichlet_boundary_conditions(fixed_value, layer_id)
        else:
            raise KeyError('either values for "flux_endplate" or '
                           '"value_endplate" must be provided')

    def stack_conductance_layers(self, conductance_list: list,
                                 modify_values=True, axis=0, **kwargs):

        # Split bipolar plate in two elements among x-direction if
        # channel-land-discretization is applied
        if self.cell.additional_layer:
            exp = kwargs.get('exponents', (-1.0, 1.0, 1.0))
            # exp = (-1.0, 1.0, 1.0)
            cat_bpp_split_ratio = (self.cell.cathode.channel.height
                                   / self.cell.cathode.bpp.thickness)
            ano_bpp_split_ratio = (self.cell.anode.channel.height
                                   / self.cell.anode.bpp.thickness)
            for i in range(len(conductance_list)):
                value = np.copy(conductance_list[i][0])
                conductance_list[i].insert(0, value * math.pow(
                    1.0 - cat_bpp_split_ratio, exp[i]))
                conductance_list[i][1] = value * math.pow(
                    cat_bpp_split_ratio, exp[i])
                value = np.copy(conductance_list[i][-1])
                conductance_list[i].append(value * math.pow(
                    1.0 - ano_bpp_split_ratio, exp[i]))
                conductance_list[i][-2] = (
                        value * math.pow(ano_bpp_split_ratio, exp[i]))
        conductance_list = [np.asarray(item) for item in conductance_list]

        # Set solid transport property to zero at channel domain
        # Land: index 0, Channel: index 1
        if (self.cell.channel_land_discretization
                and modify_values):
            for i in range(len(conductance_list)):
                conductance_list[i][[1, -2], :, 1] = 0.0
        conductance_list = self.shift_nodes(conductance_list, axis=axis)
        return conductance_list

    def add_ambient_convection(self):
        """
        Add coefficients to matrix and rhs for each cell
        due to convection to ambient
        """
        if 'alpha_amb' in self.cell.cell_dict:
            alpha_amb = self.cell.cell_dict['alpha_amb']
            temp_amb = self.cell.cell_dict['temp_amb']

            th_layer_amb = self.stack_conductance_layers(
                [self.cell.th_layer], axis=-1,
                modify_values=False,
                exponents=(1.0, ))
            k_amb = self.cell.calc_ambient_conductance(alpha_amb, th_layer_amb)
            k_amb_vector = k_amb.flatten(order='F')

            self.add_implicit_source(
                self.mtx_const, -k_amb_vector, self.index_array)
            self.add_explicit_source(
                self.rhs_const, k_amb_vector * temp_amb, self.index_array)

    @staticmethod
    def create_cell_index_list(shape: tuple[int, ...]):
        """
        Create list of lists with each list containing flattened order of indices
        for a continuous functional layer (for the conductance matrix).
        Functional layers a discretized in x-direction (z-direction in old version),
        the remaining flattened order is equal to the order in the overall matrix
        and the corresponding right-hand side vector
        (typically first y-, then z-direction within a layer a.k.a x-plane)

        Args:
            shape: tuple of size 3 containing the number of layers (index 0),
            the number of y-elements (index 1) and the number z-elements (index 2)
        """
        index_list = []
        for i in range(shape[0]):
            index_list.append(
                [(j * shape[0]) + i for j in range(shape[1] * shape[2])])
        return index_list

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        pass

    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        pass


class StackLinearSystem(LinearSystem, ABC):

    def __init__(self, shape: tuple[int, ...], stack: stack_module.Stack):
        super().__init__(shape)
        self.cells = stack.cells
        # self.transport_type = transport_type
        # Setup constant part of conductance matrix
        if isinstance(self, TemperatureSystem):
            self.transport_type = 'thermal'
        elif isinstance(self, ElectricalSystem):
            self.transport_type = 'electrical'
        else:
            raise NotImplementedError('subclass {} of class StackLinearSystem '
                                      'is not implemented'.format(type(self)))

        self.cell_systems = [cell.linear_systems[self.transport_type]
                             for cell in self.cells]
        self.index_list, self.layer_index_list = (
            self.create_stack_index_list(self.cells, self.transport_type))

        self.mtx_const = self.connect_cells(self.transport_type)
        self.rhs_const = np.hstack(
            [cell_sys.rhs_const for cell_sys in self.cell_systems])
        # if self.sparse_solve:
        #     self.mtx_const = sparse.csr_matrix(self.mtx_const)

    def add_explicit_source(self, rhs_vector, source_term, index_array,
                            layer_id=None, replace=False):
        return self.cell_systems[0].add_explicit_source(
            rhs_vector, source_term, index_array, layer_id=layer_id,
            replace=replace)

    def add_implicit_source(self, matrix, coefficients, index_array,
                            layer_id=None, replace=False):
        # matrix_size = matrix.shape[0]
        return self.cell_systems[0].add_implicit_source(
            matrix, coefficients, index_array,
            layer_id=layer_id, replace=replace)

    def connect_cells(self, transport_type: str):
        matrix = sp_la.block_diag(
            *[cell_sys.mtx_const for cell_sys in self.cell_systems])
        n_cells = len(self.cell_systems)
        cell_ids = np.asarray([list(range(n_cells-1)),
                               list(range(1, n_cells))]).transpose()
        layer_ids = np.asarray([(-1, 0) for i in range(n_cells-1)])
        conductance = np.asarray(
            [self.cell_systems[i].conductance[0][layer_ids[i][0]]
             for i in range(n_cells-1)])
        # old_matrix = np.copy(matrix)
        mtx_func.connect_cells(matrix, cell_ids, layer_ids,
                               conductance, self.index_list)
        # diff = matrix - old_matrix
        return matrix

    def set_matrix_dirichlet_bc(self, cell_ids, layer_ids):
        if not len(cell_ids) == len(layer_ids):
            raise ValueError('cell and layer index lists '
                             'must have equal length')
        for i in range(len(cell_ids)):
            row_ids = self.index_list[cell_ids[i]][:][layer_ids[i]]
            self.mtx[row_ids, :] = 0.0
            self.mtx[row_ids, row_ids] = 1.0

    @abstractmethod
    def update(self, *args, **kwargs):
        """
        This function coordinates the program sequence
        """
        pass

    @staticmethod
    def create_stack_index_list(cells: list[Cell], transport_type: str):
        n_cells = len(cells)
        index_list = []
        layer_ids = [[] for _ in range(cells[-1].nx)]
        for i in range(n_cells):
            index_array = \
                (np.prod(cells[i - 1].membrane.discretization.shape, dtype=np.int32)
                 * cells[i - 1].nx) * i \
                + cells[i].linear_systems[transport_type].index_array
            index_list.append(index_array.tolist())

        for i in range(n_cells):
            for j in range(cells[i].nx):
                layer_ids[j].append(index_list[i][j])
        layer_index_list = []
        for sub_list in layer_ids:
            layer_index_list.append(np.hstack(sub_list))
        return index_list, layer_index_list

    def update_cell_solution(self):
        """
        Transfer the general solution vector (1D) into the
        single cell solution arrays (3D)
        """

        cell_solution_vectors = np.array_split(self.solution_vector,
                                               len(self.cells))
        cell_solution_list = []
        for i, cell_sys in enumerate(self.cell_systems):
            new_cell_solution_array = np.reshape(
                cell_solution_vectors[i], cell_sys.shape, order='F')
            cell_solution_array = cell_sys.solution_array
            cell_solution_array[:] = new_cell_solution_array
            cell_solution_list.append(cell_solution_array)
        return np.asarray(cell_solution_list)


class TemperatureSystem(StackLinearSystem):

    def __init__(self, stack: stack_module.Stack, input_dict: dict):

        shape = (len(stack.cells), *stack.cells[0].thermal_system.shape)

        super().__init__(shape, stack)

        self.input_dict = input_dict

        # Sub channel ratios
        self.n_cat_channels = stack.fuel_circuits[0].n_subchannels
        self.n_ano_channels = stack.fuel_circuits[1].n_subchannels

        # Coolant flow settings
        self.cool_flow = False
        if stack.coolant_circuit is not None:
            self.cool_flow = True
            self.cool_channels = stack.coolant_circuit.channels
            self.n_cool = len(self.cool_channels)
            if not isinstance(self.cool_channels[0], chl.Channel):
                raise TypeError
            if self.n_cool == (len(self.cells) + 1):
                self.cool_ch_bc = True
            else:
                self.cool_ch_bc = False
            self.n_cell_cool_channels = stack.coolant_circuit.n_subchannels

        # Thermo-neutral cell potential
        self.e_tn = self.cells[0].e_tn
        # Open circuit potential
        self.e_0 = self.cells[0].e_0

    def update(self, *args, **kwargs):
        """
        This function coordinates the program sequence
        """
        if self.cool_flow:
            self.update_coolant_channel()
        self.update_gas_channel()

        self.update_and_solve_linear_system(gas_transfer=True,
                                            electrochemical_heat=True)
        A = self.mtx
        b = self.rhs
        x = self.solution_vector

        if np.any(self.solution_vector < 200.0):
            raise ValueError('temperature too low, check boundary conditions')
        if np.any(self.solution_vector > 1000.0):
            raise ValueError('temperature too high, check boundary conditions')
        self.update_cell_solution()

    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients, add implicit sources here.
        """

        gas_transfer = kwargs.get('gas_transfer', True)
        source_vectors = []
        for i, cell_sys in enumerate(self.cell_systems):
            cell_sys.mtx_dyn[:, :] = 0.0
            source_vectors.append(np.zeros(cell_sys.rhs_dyn.shape))

            if gas_transfer:
                # Add thermal conductance for heat transfer to cathode gas
                idx = self.cells[i].interface_id['cathode_bpp_gde']
                channel = self.cells[i].cathode.channel
                source = - channel.k_coeff  # * self.n_cat_channels

                # k_test = np.copy(channel.k_coeff)
                # k_test[:] = 0.5  # channel.k_coeff[0]
                # source = - k_test

                if self.cells[i].channel_land_discretization:
                    source = np.asarray(
                        [np.zeros(source.shape), source]).transpose()
                # source /= cell.thermal_conductance[0].shape[-1]

                matrix, source_vec_1 = self.add_implicit_source(
                    cell_sys.mtx_dyn, source, cell_sys.index_array,
                    layer_id=idx)
                source_vectors[i][:] += source_vec_1

                # Add thermal conductance for heat transfer to anode gas
                idx = self.cells[i].interface_id['anode_gde_bpp']
                channel = self.cells[i].anode.channel
                source = - channel.k_coeff  # * self.n_ano_channels

                # k_test = np.copy(channel.k_coeff)
                # k_test[:] = 0.5  # channel.k_coeff[0]
                # # k_test[:] += np.linspace(0, 0.1, k_test.shape[0])
                # source = - k_test

                if self.cells[i].channel_land_discretization:
                    source = np.asarray(
                        [np.zeros(source.shape), source]).transpose()
                # source /= cell.thermal_conductance[0].shape[-1]

                matrix, source_vec_2 = self.add_implicit_source(
                    cell_sys.mtx_dyn, source, cell_sys.index_array,
                    layer_id=idx)
                source_vectors[i][:] += source_vec_2

            # Add thermal conductance for heat transfer to coolant
            if self.cool_flow:
                source_vec_3 = np.zeros(cell_sys.rhs_dyn.shape)
                if self.cool_ch_bc:
                    cell_cool_channels = [self.cool_channels[i],
                                          self.cool_channels[i + 1]]
                    n_gas_channels = [self.n_cat_channels, self.n_ano_channels]
                    factors = [0.5, 0.5]
                    if self.cells[i].first_cell:
                        factors[0] += 0.5
                    if self.cells[i].last_cell:
                        factors[1] += 0.5
                    layer_ids = [0, -1]
                else:
                    cell_cool_channels = [self.cool_channels[i - 1],
                                          self.cool_channels[i]]
                    n_gas_channels = [self.n_ano_channels, self.n_cat_channels]
                    factors = [0.5, 0.5]
                    layer_ids = [0, -1]
                    if self.cells[i].first_cell:
                        cell_cool_channels.pop(0)
                        factors.pop(0)
                        layer_ids.pop(0)
                        n_gas_channels.pop(0)
                    if self.cells[i].last_cell:
                        cell_cool_channels.pop(1)
                        factors.pop(1)
                        layer_ids.pop(1)
                        n_gas_channels.pop(1)
                for j in range(len(cell_cool_channels)):
                    source = - cell_cool_channels[j].k_coeff
                    source *= self.n_cell_cool_channels / n_gas_channels[j]
                    source *= factors[j]

                    # if cell.channel_land_discretization:
                    #     source = np.asarray(
                    #         [np.zeros(source.shape), source]).transpose()
                    source /= cell_sys.conductance[0].shape[-1]

                    matrix, source_vec = self.add_implicit_source(
                        cell_sys.mtx_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=layer_ids[j])
                    source_vec_3[:] += source_vec
                source_vectors[i][:] += source_vec_3

        dyn_vec = np.hstack(source_vectors)
        mtx_dyn = np.diag(dyn_vec)
        self.mtx[:] = self.mtx_const + mtx_dyn

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the src to the src must be defined negative.
        """
        gas_transfer = kwargs.get('gas_transfer', True)
        electrochemical_heat = kwargs.get('electrochemical_heat', True)

        for i, cell_sys in enumerate(self.cell_systems):
            cell_sys.rhs_dyn[:] = 0.0
            if gas_transfer:
                # Cathode bpp-gde source
                idx = self.cells[i].interface_id['cathode_bpp_gde']
                # h_vap = w_prop.water.calc_h_vap(cell.cathode.channel.temp[:-1])
                channel = self.cells[i].cathode.channel
                source = channel.k_coeff * channel.temp_ele  # * self.n_cat_channels

                # temp_test = np.copy(channel.temp_ele)
                # temp_test[:] = 343.15  # channel.temp_ele[0] + 100.0
                # k_test = np.copy(channel.k_coeff)
                # # k_test[:] = 0.5  # channel.k_coeff[0]
                # source = k_test * temp_test

                source += getattr(channel, 'condensation_heat', 0.0)

                # Heat transport to reactant gas only in the part of layer at
                # channel (index 1 of z-direction)
                if self.cells[i].channel_land_discretization:
                    source = np.asarray(
                        [np.zeros(source.shape), source]).transpose()
                # source /= cell.thermal_conductance[0].shape[-1]

                cell_sys.rhs_dyn[:], _ = (
                    self.add_explicit_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=idx))

                # Anode bpp-gde
                idx = self.cells[i].interface_id['anode_gde_bpp']
                # h_vap = w_prop.water.calc_h_vap(cell.anode.temp_fluid[:-1])
                channel = self.cells[i].anode.channel
                source = channel.k_coeff * channel.temp_ele  # * self.n_ano_channels

                # temp_test = np.copy(channel.temp_ele)
                # temp_test[:] = 343.15  # channel.temp_ele[0] + 100.0
                # k_test = np.copy(channel.k_coeff)
                # # k_test[:] = 0.5  # channel.k_coeff[0]
                # # k_test[:] += np.linspace(0, 0.1, k_test.shape[0])
                # source = k_test * temp_test

                source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0

                if self.cells[i].channel_land_discretization:
                    source = np.asarray(
                        [np.zeros(source.shape), source]).transpose()
                # source /= cell.thermal_conductance[0].shape[-1]

                cell_sys.rhs_dyn[:], _ = (
                    self.add_explicit_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=idx))

            if electrochemical_heat:
                # Cathode gde-mem source
                idx = self.cells[i].interface_id['cathode_gde_mem']
                current = (self.cells[i].current_density[idx]
                           * self.cells[i].d_area)

                cathode_ohmic_heat_membrane = (
                        0.5 * self.cells[i].membrane.omega * np.square(current))
                source = np.copy(cathode_ohmic_heat_membrane)

                # test_current_density = np.zeros(current.shape)
                # test_current_density[:] = 30000.0
                # test_current = test_current_density * cell.d_area
                # test_omega = np.copy(cell.membrane.omega)
                # test_omega[:] = test_omega[0]
                # test_ohmic_heat = 0.5 * test_omega * np.square(test_current)
                # source = np.copy(test_ohmic_heat)

                v_loss = np.minimum(
                    self.e_0, self.cells[i].cathode.voltage_loss)
                v_loss[v_loss < 0.0] = 0.0
                reaction_heat = (self.e_tn - self.e_0 + v_loss) * current
                # reaction_heat_sum = np.sum(reaction_heat)

                # test_v_loss = np.copy(v_loss)
                # test_v_loss[:] = 0.3
                # test_reaction_heat = (
                #         (self.e_tn - self.e_0 + test_v_loss) * test_current)
                # test_reaction_heat_sum = np.sum(test_reaction_heat)

                source += reaction_heat
                # source *= 1.0
                cell_sys.rhs_dyn[:], _ = (
                    self.add_explicit_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=idx))

                # Anode gde-mem source
                idx = self.cells[i].interface_id['anode_mem_gde']
                current = (self.cells[i].current_density[idx]
                           * self.cells[i].d_area)
                anode_ohmic_heat_membrane = (
                        0.5 * self.cells[i].membrane.omega * np.square(current))
                source = anode_ohmic_heat_membrane
                v_loss = np.minimum(self.e_0, self.cells[i].anode.voltage_loss)
                v_loss[v_loss < 0.0] = 0.0
                reaction_heat = v_loss * current
                source += reaction_heat
                source *= 1.0
                cell_sys.rhs_dyn[:], _ = self.add_explicit_source(
                    cell_sys.rhs_dyn, source.flatten(order='F'),
                    cell_sys.index_array, layer_id=idx)

            # Coolant source coefficients
            if self.cool_flow:
                if self.cool_ch_bc:
                    cell_cool_channels = [self.cool_channels[i],
                                          self.cool_channels[i + 1]]
                    n_gas_channels = [self.n_cat_channels, self.n_ano_channels]
                    factors = [0.5, 0.5]
                    if self.cells[i].first_cell:
                        factors[0] += 0.5
                    if self.cells[i].last_cell:
                        factors[1] += 0.5
                    layer_ids = [0, -1]
                else:
                    cell_cool_channels = [self.cool_channels[i - 1],
                                          self.cool_channels[i]]
                    n_gas_channels = [self.n_ano_channels, self.n_cat_channels]
                    factors = [0.5, 0.5]
                    layer_ids = [0, -1]
                    if self.cells[i].first_cell:
                        cell_cool_channels.pop(0)
                        factors.pop(0)
                        layer_ids.pop(0)
                        n_gas_channels.pop(0)
                    if self.cells[i].last_cell:
                        cell_cool_channels.pop(1)
                        factors.pop(1)
                        layer_ids.pop(1)
                        n_gas_channels.pop(1)
                for j, cool_chl in enumerate(cell_cool_channels):
                    source = cool_chl.k_coeff * cool_chl.temp_ele
                    source *= self.n_cell_cool_channels / n_gas_channels[j]
                    # if cell.channel_land_discretization:
                    #     source = np.asarray(
                    #         [np.zeros(source.shape), source]).transpose()
                    source /= cell_sys.conductance[0].shape[-1]
                    source *= factors[j]
                    cell_sys.rhs_dyn[:], _ = (
                        self.add_explicit_source(
                            cell_sys.rhs_dyn, source.flatten(order='F'),
                            cell_sys.index_array, layer_id=layer_ids[j]))

        rhs_dyn = np.hstack([cell_sys.rhs_dyn for cell_sys in
                             self.cell_systems])
        self.rhs[:] = self.rhs_const + rhs_dyn

    def update_gas_channel(self):
        """
        Calculates the fluid temperatures in the anode and cathode channels
        """
        for i, cell in enumerate(self.cells):
            idx = cell.interface_id['cathode_bpp_gde']
            cell.cathode.channel.update_heat(
                wall_temp=g_func.reduce_dimension(cell.temp_layer[idx]),
                update_fluid=False)
            idx = cell.interface_id['anode_gde_bpp']
            cell.anode.channel.update_heat(
                wall_temp=g_func.reduce_dimension(cell.temp_layer[idx]),
                update_fluid=False)

    def update_coolant_channel(self):
        """
        Calculates the coolant channel temperatures.
        """
        for i, cool_chl in enumerate(self.cool_channels):
            if self.cool_ch_bc:
                if i == self.n_cool - 1:
                    idx = self.cells[i - 1].interface_id['anode_bpp_bc']
                    wall_temp = self.cells[i - 1].temp_layer[idx]
                else:
                    idx = self.cells[i - 1].interface_id['cathode_bc_bpp']
                    wall_temp = self.cells[i].temp_layer[idx]
            else:
                idx = self.cells[i + 1].interface_id['cathode_bc_bpp']
                wall_temp = self.cells[i + 1].temp_layer[idx]
            wall_temp = g_func.reduce_dimension(wall_temp)
            cool_chl.update_heat(wall_temp=wall_temp, update_fluid=False)


class ElectricalSystem(StackLinearSystem):

    def __init__(self, stack: stack_module.Stack, input_dict: dict):

        shape = (len(stack.cells), *stack.cells[0].electrical_system.shape)
        super().__init__(shape, stack)

        """Building up the base conductance matrix"""
        # Vector for dynamically changing heat conductance values
        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        self.input_dict = input_dict

        # Variables
        self.current_control = stack.current_control
        if self.current_control:
            self.current_density_target = stack.current_density_target
        else:
            self.v_tar = stack.voltage_target
            self.e_0_stack = np.sum([cell.e_0 for cell in self.cells])
            self.v_loss_tar = self.e_0_stack - self.v_tar

        # Current density of the elements in trough-plane- / x-direction
        diff_shape = (len(self.cells),
                      *self.cell_systems[0].conductance[0].shape)
        self.current_density = np.zeros(diff_shape)

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        # raise NotImplementedError
        cell_rhs_list = []

        # Current boundary conditions are applied at outer plate
        # (layer id 0) of cell 0
        index = 0
        cell = self.cells[index]
        cell_sys = self.cell_systems[index]
        cell_rhs = np.zeros(cell.voltage_layer.flatten().shape)
        if self.current_control:
            bc_current = self.current_density_target * cell.d_area
            if np.sum(cell_sys.solution_array[0]) == 0.0:
                rhs_bc_values = self.current_density_target * cell.d_area
            else:
                inlet_current = (
                    np.abs(cell_sys.solution_array[0]
                           - cell_sys.solution_array[1])
                    * self.cell_systems[index].conductance[0][0])
                correction_factors = bc_current / inlet_current
                rhs_bc_values = inlet_current * correction_factors

        else:
            rhs_bc_values = self.v_loss_tar
        self.add_explicit_source(
            cell_rhs, rhs_bc_values.flatten(order='F'),
            cell_sys.index_array, layer_id=0)
        cell_rhs_list.append(cell_rhs)
        for i in range(1, len(self.cells)):
            cell_rhs_list.append(
                np.zeros(cell_sys.solution_array.flatten().shape))
        self.rhs[:] = np.hstack(cell_rhs_list)

    def update_matrix(self, electrochemical_conductance, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        mtx_list = \
            [mtx_func.build_cell_conductance_matrix([cond])
             for cond in electrochemical_conductance]
        mtx_dyn = sp_la.block_diag(*mtx_list)
        self.mtx[:] = self.mtx_const + mtx_dyn
        # Modify matrix for dirichlet boundary conditions in last cell
        # (and first cell if voltage is given as BC)
        self.set_matrix_dirichlet_bc([-1], [-1])
        if not self.current_control:
            self.set_matrix_dirichlet_bc([0], [0])

    def update(self, current_density=None, voltage=None):
        """
        Wrapper function to solve underlying linear system with updated input
        and update every member results accordingly
        """
        if current_density is not None:
            self.current_density_target = current_density
        if voltage is not None:
            self.v_tar = voltage
            self.v_loss_tar = self.e_0_stack - self.v_tar

        dynamic_conductances = self.create_dynamic_conductance_list()
        self.update_and_solve_linear_system(dynamic_conductances)
        voltage_array = self.update_cell_voltage()
        constant_conductances = np.asarray([cell_sys.conductance[0]
                                            for cell_sys in self.cell_systems])
        conductance_array = (constant_conductances
                             + np.asarray(dynamic_conductances))
        self.update_current_density(voltage_array, conductance_array)

    def create_dynamic_conductance_list(self):
        dynamic_cell_conductance_list = []
        for i, cell_sys in enumerate(self.cell_systems):
            shape = cell_sys.conductance[0].shape
            dynamic_cell_conductance = np.zeros(shape)
            layer_id = (cell_sys.conductance[0].shape[0] // 2)
            dynamic_cell_conductance[layer_id] = (
                self.cells[i].electrochemical_conductance)
            dynamic_cell_conductance_list.append(dynamic_cell_conductance)
        return dynamic_cell_conductance_list

    def update_cell_voltage(self):
        """
        Wrapper for the update_cell_solution-function which the cell-wise parts
        of the flattened solution_vector into the corresponding cell member (
        voltage_layer). Furthermore, the area-distribution of the cell voltage
        loss is calculated from cathode to anode bipolar plate for each cell.
        Returns: the reshaped voltage_array containing the over-potential in
        each cell layer discretized according to the y-z-discretization of
        each cell.
        """
        voltage_array = self.update_cell_solution()
        for i, cell in enumerate(self.cells):
            # cell.voltage_layer[:] = cell.e_0 - cell.voltage_layer
            cell.voltage_loss[:] = np.abs(
                voltage_array[i][0] - voltage_array[i][-1])
            # cell.update_voltage_loss(v_diff[i])
        return voltage_array

    def update_current_density(self, voltage_array, conductance):
        """
        Calculates the current density in each layer of each in through-plane
        direction, i.e. the x-direction
        """
        v_diff = np.abs(np.diff(voltage_array, axis=1))
        active_area = np.array([cell.d_area for cell in self.cells])
        active_area_array = np.asarray(
            [np.stack([cell.d_area for i in range(v_diff.shape[1])], axis=0)
             for cell in self.cells])
        current_density = v_diff * conductance / active_area_array
        self.current_density[:] = current_density
        for i, cell in enumerate(self.cells):
            cell.current_density[:] = current_density[i]
        return current_density


# disc_dict = {
#     'length': 1.0,
#     'width': 1.0,
#     'depth': 1.0,
#     'shape': (10, 2),
#     'ratio': (1, 1),
#     'direction': (1, 1),
# }
#
# discretization_2d = dsct.Discretization(disc_dict)
# discretization_3d = dsct.Discretization3D.create_from(
#     discretization_2d, 0.4, 3)
#
# trans_layer = tl.TransportLayer.create({}, {'diffusion': [1.0, 2.0, 3.0]},
#                                        discretization_3d)
# lin_sys = BasicLinearSystem.create(trans_layer, trans_layer.types[0])
