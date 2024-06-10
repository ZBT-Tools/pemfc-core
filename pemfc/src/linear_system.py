# General imports
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING
import numpy as np
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve
# Local module imports
from . import (
    matrix_functions as mtx_func, stack as stack_module, channel as chl,
    global_functions as g_func)

if TYPE_CHECKING:
    from pemfc.src.stack import Stack

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
            mtx = sparse.csr_matrix(self.mtx)
            self.solution_vector[:] = spsolve(mtx, self.rhs)
        else:
            self.solution_vector[:] = np.linalg.tensorsolve(self.mtx, self.rhs)


class StackLinearSystem(LinearSystem, ABC):

    def __init__(self, shape: tuple[int, ...], stack: stack_module.Stack):
        super().__init__(shape)
        self.cells = stack.cells
        # self.transport_type = transport_type
        self.index_list, self.layer_index_list = (
            mtx_func.create_stack_index_list(self.cells))

        # Setup constant part of conductance matrix
        if isinstance(self, TemperatureSystem):
            transport_type = 'thermal'
        elif isinstance(self, ElectricalSystem):
            transport_type = 'electrical'
        else:
            raise NotImplementedError('subclass {} of class StackLinearSystem '
                                      'is not implemented'.format(type(self)))
        self.mtx_const = self.connect_cells(transport_type)
        # if self.sparse_solve:
        #     self.mtx_const = sparse.csr_matrix(self.mtx_const)

    def connect_cells(self, transport_type: str):
        matrix = sp_la.block_diag(*[cell.mtx_const[transport_type]
                                    for cell in self.cells])
        n_cells = len(self.cells)
        cell_ids = np.asarray([list(range(n_cells-1)),
                               list(range(1, n_cells))]).transpose()
        layer_ids = np.asarray([(-1, 0) for i in range(n_cells-1)])
        conductance = np.asarray(
            [self.cells[i].conductance[transport_type][0][layer_ids[i][0]]
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

    @abstractmethod
    def update_cell_solution(self):
        """
        This function coordinates the program sequence

        Returns:

        """
        pass

    def update_cell_solution_general(self, cell_solution_name: str,
                                     cell_shape_name: str):
        """
        Transfer the general solution vector (1D) into the single cell solution arrays (3D)
        """

        cell_solution_vectors = np.array_split(self.solution_vector,
                                               len(self.cells))
        cell_solution_list = []
        for i, cell in enumerate(self.cells):
            # cell_shape = getattr(cell, cell_shape_name)
            cell_solution_array = np.reshape(cell_solution_vectors[i],
                                             self.cells[i].voltage_shape,
                                             order='F')
            cell_solution = getattr(cell, cell_solution_name)
            cell_solution[:] = cell_solution_array
            cell_solution_list.append(cell_solution_array)
        return np.asarray(cell_solution_list)


class TemperatureSystem(StackLinearSystem):

    def __init__(self, stack: stack_module.Stack, input_dict: dict):

        shape = (len(stack.cells), *stack.cells[0].temp_shape)
        super().__init__(shape, stack)

        """Building up the base conductance matrix"""
        # Vector for dynamically changing heat conductance values
        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        self.input_dict = input_dict

        self.rhs_const = np.hstack(
            [cell.thermal_rhs_const for cell in self.cells])

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
            self.n_cool_sub_channels = stack.coolant_circuit.n_subchannels

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
        self.update_and_solve_linear_system()
        if np.any(self.solution_vector < 200.0):
            raise ValueError('temperature too low, check boundary conditions')
        if np.any(self.solution_vector > 1000.0):
            raise ValueError('temperature too high, check boundary conditions')
        self.update_cell_solution()

    def update_gas_channel(self):
        """
        Calculates the fluid temperatures in the anode and cathode channels
        """
        for i, cell in enumerate(self.cells):
            cell.cathode.channel.update_heat(
                wall_temp=g_func.retrieve_dimension(cell.temp_layer[1]),
                update_fluid=False)
            cell.anode.channel.update_heat(
                wall_temp=g_func.retrieve_dimension(cell.temp_layer[4]),
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
            wall_temp = g_func.retrieve_dimension(wall_temp)
            cool_chl.update_heat(wall_temp=wall_temp, update_fluid=False)

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the src to the src must be defined negative.
        """
        for i, cell in enumerate(self.cells):
            cell.thermal_rhs_dyn[:] = 0.0

            # Cathode bpp-gde source
            layer_id = cell.layer_id['cathode_gde']
            # h_vap = w_prop.water.calc_h_vap(cell.cathode.channel.temp[:-1])
            channel = cell.cathode.channel
            source = channel.k_coeff * channel.temp_ele  # * self.n_cat_channels
            source += getattr(channel, 'condensation_heat', 0.0)
            # Heat transport to reactant gas only in the part of layer at
            # channel (index 1 of z-direction)
            if cell.channel_land_discretization:
                source = np.asarray([np.zeros(source.shape), source])
            mtx_func.add_explicit_layer_source(
                cell.thermal_rhs_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Cathode gde-mem source
            layer_id = cell.layer_id['cathode_gde'] + 1
            current = (cell.current_density[layer_id] * cell.d_area)
            cathode_ohmic_heat_membrane = \
                0.5 * cell.membrane.omega * np.square(current)
            source = cathode_ohmic_heat_membrane
            v_loss = np.minimum(self.e_0, cell.cathode.voltage_loss)
            v_loss[v_loss < 0.0] = 0.0
            reaction_heat = (self.e_tn - self.e_0 + v_loss) * current
            source += reaction_heat
            mtx_func.add_explicit_layer_source(
                cell.thermal_rhs_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Anode gde-mem source
            layer_id = cell.layer_id['anode_gde']
            current = (cell.current_density[layer_id] * cell.d_area)
            anode_ohmic_heat_membrane = (
                    0.5 * cell.membrane.omega * np.square(current))
            source = anode_ohmic_heat_membrane
            v_loss = np.minimum(self.e_0, cell.anode.voltage_loss)
            v_loss[v_loss < 0.0] = 0.0
            reaction_heat = v_loss * current
            source += reaction_heat
            mtx_func.add_explicit_layer_source(
                cell.thermal_rhs_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Anode bpp-gde
            layer_id = cell.layer_id['anode_gde'] + 1
            # h_vap = w_prop.water.calc_h_vap(cell.anode.temp_fluid[:-1])
            channel = cell.anode.channel
            source = channel.k_coeff * channel.temp_ele  # * self.n_ano_channels
            source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0
            if cell.channel_land_discretization:
                source = np.asarray([np.zeros(source.shape), source])
            mtx_func.add_explicit_layer_source(
                cell.thermal_rhs_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Cooling channels
            if self.cool_flow:
                if self.cool_ch_bc:
                    cool_chl = self.cool_channels[i]
                    source = cool_chl.k_coeff * cool_chl.temp_ele
                    source *= self.n_cool_sub_channels
                    mtx_func.add_explicit_layer_source(
                        cell.thermal_rhs_dyn, source,
                        cell.index_array, layer_id=0)
                    if cell.last_cell:
                        cool_chl = self.cool_channels[i + 1]
                        source = cool_chl.k_coeff * cool_chl.temp_ele
                        source *= self.n_cool_sub_channels
                        mtx_func.add_explicit_layer_source(
                            cell.thermal_rhs_dyn, source,
                            cell.index_array, layer_id=-1)
                else:
                    if not cell.first_cell:
                        cool_chl = self.cool_channels[i - 1]
                        source = cool_chl.k_coeff * cool_chl.temp_ele
                        source *= self.n_cool_sub_channels
                        mtx_func.add_explicit_layer_source(
                            cell.thermal_rhs_dyn, source,
                            cell.index_array, layer_id=0)

        rhs_dyn = np.hstack([cell.thermal_rhs_dyn for cell in self.cells])
        self.rhs[:] = self.rhs_const + rhs_dyn

    def update_matrix(self, *args, **kwargs):
        """
        Updates matrix coefficients, add implicit sources here.
        """
        source_vectors = []
        for i, cell in enumerate(self.cells):
            cell.thermal_mtx_dyn[:, :] = 0.0
            source_vectors.append(np.zeros(cell.thermal_rhs_dyn.shape))

            # Add thermal conductance for heat transfer to cathode gas
            layer_id = cell.layer_id['cathode_gde']
            source = -cell.cathode.channel.k_coeff  # * self.n_cat_channels
            if cell.channel_land_discretization:
                source = np.asarray([np.zeros(source.shape), source])
            matrix, source_vec_1 = mtx_func.add_implicit_layer_source(
                cell.thermal_mtx_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Add thermal conductance for heat transfer to anode gas
            layer_id = cell.layer_id['anode_gde'] + 1
            source = -cell.anode.channel.k_coeff  # * self.n_ano_channels
            if cell.channel_land_discretization:
                source = np.asarray([np.zeros(source.shape), source])
            matrix, source_vec_2 = mtx_func.add_implicit_layer_source(
                cell.thermal_mtx_dyn, source, cell.index_array,
                layer_id=layer_id)

            # Add thermal conductance for heat transfer to coolant
            source_vec_3 = np.zeros(source_vec_1.shape)
            # TODO: implement correct cooling with new discretization
            raise NotImplementedError('implement correct cooling')
            if self.cool_flow:
                if self.cool_ch_bc:
                    source = - self.cool_channels[i].k_coeff
                    source *= self.n_cool_sub_channels / self.n_cat_channels
                    matrix, source_vec = mtx_func.add_implicit_layer_source(
                        cell.thermal_mtx_dyn, source,
                        cell.index_array, layer_id=0)
                    source_vec_3[:] = source_vec
                    if cell.last_cell:
                        source = - self.cool_channels[i + 1].k_coeff
                        source *= self.n_cool_sub_channels / self.n_ano_channels
                        matrix, source_vec = mtx_func.add_implicit_layer_source(
                            cell.thermal_mtx_dyn, source,
                            cell.index_array, layer_id=-1)
                        source_vec_3[:] += source_vec
                else:
                    if not cell.first_cell:
                        source = - self.cool_channels[i - 1].k_coeff
                        source *= self.n_cool_sub_channels / self.n_cat_channels
                        matrix, source_vec = mtx_func.add_implicit_layer_source(
                            cell.thermal_mtx_dyn, source,
                            cell.index_array, layer_id=0)
                        source_vec_3[:] = source_vec

            source_vectors[i][:] = source_vec_1 + source_vec_2 + source_vec_3

        dyn_vec = np.hstack(source_vectors)
        # if self.sparse_solve:
        #     self.mtx = \
        #         self.mtx_const + sparse.diags([dyn_vec], [0], format='csr')
        # else:
        self.mtx[:] = self.mtx_const + np.diag(dyn_vec)

    def update_cell_solution(self):
        temp_solution_name = 'temp_layer'
        if not hasattr(self.cells[0], temp_solution_name):
            raise ValueError('attribute {} is not found in {}'.format(
                temp_solution_name, self.cells[0].__name__))
        temp_shape_name = 'temp_shape'
        if not hasattr(self.cells[0], temp_shape_name):
            raise ValueError('attribute {} is not found in {}'.format(
                temp_shape_name, self.cells[0].__name__))
        return self.update_cell_solution_general(temp_solution_name,
                                                 temp_shape_name)


class ElectricalSystem(StackLinearSystem):

    def __init__(self, stack: stack_module.Stack, input_dict: dict):

        shape = (len(stack.cells), *stack.cells[0].temp_shape)
        super().__init__(shape, stack)

        """Building up the base conductance matrix"""
        # Vector for dynamically changing heat conductance values
        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        self.input_dict = input_dict

        # self.rhs_const = np.hstack([cell.heat_rhs_const for cell in self.cells])

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
                      *self.cells[0].electrical_conductance[0].shape)
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
        cell_0 = self.cells[0]
        cell_rhs = np.zeros(cell_0.voltage_layer.flatten().shape)
        if self.current_control:
            bc_current = self.current_density_target * cell_0.d_area
            if np.sum(cell_0.voltage_layer[0]) == 0.0:
                rhs_bc_values = self.current_density_target * cell_0.d_area
            else:
                inlet_current = (
                    np.abs(cell_0.voltage_layer[0] - cell_0.voltage_layer[1])
                    * cell_0.electrical_conductance[0][0])
                correction_factors = bc_current / inlet_current
                rhs_bc_values = inlet_current * correction_factors

        else:
            rhs_bc_values = self.v_loss_tar
        mtx_func.add_explicit_layer_source(cell_rhs, rhs_bc_values,
                                           cell_0.index_array, layer_id=0)
        cell_rhs_list.append(cell_rhs)
        for i in range(1, len(self.cells)):
            cell_rhs_list.append(
                np.zeros(self.cells[i].voltage_layer.flatten().shape))
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

    def update_cell_solution(self):
        solution_name = 'voltage_layer'
        if not hasattr(self.cells[0], solution_name):
            raise ValueError('attribute {} is not found in {}'.format(
                solution_name, self.cells[0].__name__))
        shape_name = 'voltage_shape'
        if not hasattr(self.cells[0], shape_name):
            raise ValueError('attribute {} is not found in {}'.format(
                shape_name, self.cells[0].__name__))
        return self.update_cell_solution_general(solution_name, shape_name)

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
        constant_conductances = np.asarray([cell.electrical_conductance[0]
                                            for cell in self.cells])
        conductance_array = (constant_conductances
                             + np.asarray(dynamic_conductances))
        self.update_current_density(voltage_array, conductance_array)

    def create_dynamic_conductance_list(self):
        dynamic_cell_conductance_list = []
        for cell in self.cells:
            shape = self.cells[0].electrical_conductance[0].shape

            dynamic_cell_conductance = np.zeros(shape)
            layer_id = (
                self.cells[0].electrical_conductance[0].shape[0] // 2)
            dynamic_cell_conductance[layer_id] = (
                cell.electrochemical_conductance)
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
            cell.voltage_loss[:] = np.abs(voltage_array[i][0] - voltage_array[i][-1])
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
