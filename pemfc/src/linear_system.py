# General imports
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING
import math
import numpy as np
from numpy import linalg
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve
# Local module imports
from . import (
    matrix_functions as mtx_func, stack as stack_module, channel as chl,
    global_functions as g_func, transport_layer as tl)

if TYPE_CHECKING:
    from pemfc.src.stack import Stack
    from pemfc.src.cell import Cell

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


class CellLinearSystem(LinearSystem):

    def __init__(self, cell: Cell, transport_type: str, init_value=0.0):
        self.cell = cell
        self.conductance = self.calculate_conductance(
            self.cell.layers, transport_type)
        shape = self.conductance[1].shape
        super().__init__(shape)
        self.type = transport_type

        self.mtx_const = mtx_func.build_cell_conductance_matrix(
                [self.conductance[0],
                 tl.TransportLayer.calc_inter_node_conductance(
                    self.conductance[1], axis=1) * 1.0,
                 tl.TransportLayer.calc_inter_node_conductance(
                    self.conductance[2], axis=2) * 1.0])

        self.mtx_dyn = np.zeros(self.mtx_const.shape)
        self.rhs_const = np.zeros(self.rhs.shape)
        self.rhs_dyn = np.zeros(self.rhs.shape)
        self.index_array = self.create_cell_index_list(self.shape)

        # Initialize solution values
        self.solution_vector[:] = init_value
        self.solution_array[:] = init_value

    def calculate_conductance(self, layers: list[tl.TransportLayer],
                              transport_type: str):
        if transport_type not in ('electrical', 'thermal'):
            raise ValueError("transport_type argument must be either "
                             "'electrical' or 'thermal'")

        # Stack thermal conductances along through-plane direction,
        # i.e. x-coordinate
        dims = len(layers[0].conductance[transport_type])
        conductance = [
            [layer.conductance[transport_type][i] for layer in layers]
            for i in range(dims)]
        return self.stack_cell_property(conductance, exp=(-1.0, 1.0, 1.0),
                                        stacking_axis=0, modify_values=True)

    def set_layer_boundary_conditions(self, layer_id):
        if 'flux_endplate' in self.cell.cell_dict:
            flux_endplate = self.cell.cell_dict['flux_endplate']
            source = flux_endplate * self.cell.cathode.discretization.d_area
            mtx_func.add_explicit_layer_source(
                self.rhs_const, source.flatten(order='F'),
                self.index_array, layer_id)

        elif 'temp_endplate' in self.cell.cell_dict:
            mtx_func.set_implicit_layer_fixed(self.mtx_const,
                                         self.index_array, layer_id)
            source = self.cell.cell_dict['temp_endplate']
            self.rhs_const[:], _ = mtx_func.add_explicit_layer_source(
                self.rhs_const, -source, self.index_array, layer_id,
                replace=True)
        else:
            raise KeyError('either values for "flux_endplate" or '
                           '"temp_endplate" must be provided')

    def stack_cell_property(self, cell_property: list, stacking_axis,
                            exp: tuple,
                            modify_values=False,
                            shift_along_axis=(False, True, True)):

        # Split bipolar plate in two elements among x-direction if
        # channel-land-discretization is applied
        if self.cell.additional_layer:
            cat_bpp_split_ratio = (self.cell.cathode.channel.height
                                   / self.cell.cathode.bpp.thickness)
            ano_bpp_split_ratio = (self.cell.anode.channel.height
                                   / self.cell.anode.bpp.thickness)
            for i in range(len(cell_property)):
                value = np.copy(cell_property[i][0])
                cell_property[i].insert(0, value * math.pow(
                    1.0 - cat_bpp_split_ratio, exp[i]))
                cell_property[i][1] = value * math.pow(
                    cat_bpp_split_ratio, exp[i])
                value = np.copy(cell_property[i][-1])
                cell_property[i].append(value * math.pow(
                    1.0 - ano_bpp_split_ratio, exp[i]))
                cell_property[i][-2] = (
                        value * math.pow(ano_bpp_split_ratio, exp[i]))
        cell_property = [np.asarray(item) for item in cell_property]

        # Set solid transport property to zero at channel domain
        # Channel: index 1, Land: index 0
        if self.cell.channel_land_discretization and modify_values:
            for i in range(len(cell_property)):
                cell_property[i][[1, -2], :, 1] = 0.0

        for i in range(len(cell_property)):
            if shift_along_axis[i]:
                cell_property[i] = (
                        (cell_property[i]
                         + np.roll(cell_property[i], 1, axis=0)) * 0.5)
                cell_property[i] = np.concatenate(
                    (cell_property[i], [cell_property[i][0]]),
                    axis=stacking_axis)
                cell_property[i][0] *= 0.5
                cell_property[i][-1] *= 0.5
        return cell_property

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
        # if self.sparse_solve:
        #     self.mtx_const = sparse.csr_matrix(self.mtx_const)

    def connect_cells(self, transport_type: str):
        matrix = sp_la.block_diag(
            *[lin_sys.mtx_const for lin_sys in self.cell_systems])
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
                (np.prod(cells[i - 1].membrane.dsct.shape, dtype=np.int32)
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

        """Building up the base conductance matrix"""
        # Vector for dynamically changing heat conductance values
        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        self.input_dict = input_dict

        self.rhs_const = np.hstack(
            [lin_sys.rhs_const for lin_sys in self.cell_systems])

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

    def add_ambient_convection(self):
        """
        Add coefficients to matrix and rhs for each cell
        due to convection to ambient
        """
        alpha_amb = self.input_dict['alpha_amb']
        temp_amb = self.input_dict['temp_amb']
        for cell in self.cells:
            cell.k_amb = cell.calc_ambient_conductance(alpha_amb)
            # if cell.last_cell:
            k_amb_vector = cell.k_amb.flatten(order='F')
            # else:
            #     k_amb_vector = cell.k_amb[:-1].flatten(order='F')

            mtx_func.add_implicit_layer_source(
                cell.thermal_mtx_const, -k_amb_vector, cell.index_array)
            mtx_func.add_explicit_layer_source(
                cell.thermal_rhs_const, k_amb_vector * temp_amb,
                cell.index_array)

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
                layer_id = self.cells[i].layer_id['cathode_gde']
                channel = self.cells[i].cathode.channel
                source = - channel.k_coeff  # * self.n_cat_channels

                # k_test = np.copy(channel.k_coeff)
                # k_test[:] = 0.5  # channel.k_coeff[0]
                # source = - k_test

                if self.cells[i].channel_land_discretization:
                    source = np.asarray(
                        [np.zeros(source.shape), source]).transpose()
                # source /= cell.thermal_conductance[0].shape[-1]

                matrix, source_vec_1 = mtx_func.add_implicit_layer_source(
                    cell_sys.mtx_dyn, source, cell_sys.index_array,
                    layer_id=layer_id)
                source_vectors[i][:] += source_vec_1

                # Add thermal conductance for heat transfer to anode gas
                layer_id = self.cells[i].layer_id['anode_gde'] + 1
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

                matrix, source_vec_2 = mtx_func.add_implicit_layer_source(
                    cell_sys.mtx_dyn, source, cell_sys.index_array,
                    layer_id=layer_id)
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

                    matrix, source_vec = mtx_func.add_implicit_layer_source(
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
                layer_id = self.cells[i].layer_id['cathode_gde']
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
                    mtx_func.add_explicit_layer_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=layer_id))

                # Anode bpp-gde
                layer_id = self.cells[i].layer_id['anode_gde'] + 1
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
                    mtx_func.add_explicit_layer_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=layer_id))

            if electrochemical_heat:
                # Cathode gde-mem source
                layer_id = self.cells[i].layer_id['cathode_gde'] + 1
                current = (self.cells[i].current_density[layer_id]
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
                    mtx_func.add_explicit_layer_source(
                        cell_sys.rhs_dyn, source.flatten(order='F'),
                        cell_sys.index_array, layer_id=layer_id))

                # Anode gde-mem source
                layer_id = self.cells[i].layer_id['anode_gde']
                current = (self.cells[i].current_density[layer_id]
                           * self.cells[i].d_area)
                anode_ohmic_heat_membrane = (
                        0.5 * self.cells[i].membrane.omega * np.square(current))
                source = anode_ohmic_heat_membrane
                v_loss = np.minimum(self.e_0, self.cells[i].anode.voltage_loss)
                v_loss[v_loss < 0.0] = 0.0
                reaction_heat = v_loss * current
                source += reaction_heat
                source *= 1.0
                cell_sys.rhs_dyn[:], _ = mtx_func.add_explicit_layer_source(
                    cell_sys.rhs_dyn, source.flatten(order='F'),
                    cell_sys.index_array, layer_id=layer_id)

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
                        mtx_func.add_explicit_layer_source(
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
            layer_id = cell.layer_id['cathode_gde']
            cell.cathode.channel.update_heat(
                wall_temp=g_func.reduce_dimension(cell.temp_layer[layer_id]),
                update_fluid=False)
            layer_id = cell.layer_id['anode_gde']
            cell.anode.channel.update_heat(
                wall_temp=g_func.reduce_dimension(cell.temp_layer[layer_id]),
                update_fluid=False)

    def update_coolant_channel(self):
        """
        Calculates the coolant channel temperatures.
        """
        for i, cool_chl in enumerate(self.cool_channels):
            if self.cool_ch_bc:
                if i == self.n_cool - 1:
                    wall_temp = self.cells[i - 1].temp_layer[-1]
                else:
                    wall_temp = self.cells[i].temp_layer[0]
            else:
                wall_temp = self.cells[i + 1].temp_layer[0]
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
                    np.abs(cell_sys.solution_array[0] - cell_sys.solution_array[1])
                    * self.cell_systems[index].conductance[0][0])
                correction_factors = bc_current / inlet_current
                rhs_bc_values = inlet_current * correction_factors

        else:
            rhs_bc_values = self.v_loss_tar
        mtx_func.add_explicit_layer_source(
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
