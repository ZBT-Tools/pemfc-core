# general imports
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve

if TYPE_CHECKING:
    from stack import Stack

# local module imports
from . import matrix_functions as mtx_func, stack as stack_module, channel as chl

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
        # TODO: Check shape of mtx, not conform to mtx_const and mtx_dyn
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

    def update_and_solve_linear_system(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
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
        conductance = \
            np.asarray([self.cells[i].conductance[transport_type][0][layer_ids[i][0]]
                        for i in range(n_cells-1)])
        # old_matrix = np.copy(matrix)
        mtx_func.connect_cells(matrix, cell_ids, layer_ids,
                               conductance, self.index_list)
        # diff = matrix - old_matrix
        return matrix

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

    def update_cell_solution_general(self, cell_solution_name: str, cell_shape_name: str):
        """
        Transfer the general solution vector (1D) into the single cell solution arrays (3D)
        """
        # TODO: Rework for 3D required
        raise NotImplementedError
        for i, cell in enumerate(self.cells):
            n_layer = getattr(cell, cell_shape_name)[0]
            for j in range(n_layer):
                index_vector = self.index_list[i][j]
                cell_solution_vector = getattr(cell, cell_solution_name)
                cell_solution_vector[j] = self.solution_vector[index_vector]


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

        # coolant flow settings
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

        # thermodynamic neutral cell potential
        self.e_tn = self.cells[0].e_tn
        # open circuit potential
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
            cell.cathode.channel.update_heat(wall_temp=cell.temp_layer[1],
                                             update_fluid=False)
            cell.anode.channel.update_heat(wall_temp=cell.temp_layer[4],
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
            cool_chl.update_heat(wall_temp=wall_temp, update_fluid=False)

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the src to the src must be defined negative.
        """
        for i, cell in enumerate(self.cells):
            cell.thermal_rhs_dyn[:] = 0.0

            source = np.zeros(cell.temp_layer[0].shape)
            # Cathode bpp-gde source
            # h_vap = w_prop.water.calc_h_vap(cell.cathode.channel.temp[:-1])
            channel = cell.cathode.channel
            source += channel.k_coeff * channel.temp_ele * self.n_cat_channels
            source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0
            mtx_func.add_explicit_layer_source(cell.thermal_rhs_dyn, source,
                                               cell.index_array, layer_id=1)

            current = cell.i_cd * cell.d_area
            half_ohmic_heat_membrane = \
                0.5 * cell.membrane.omega * np.square(current)

            # Cathode gde-mem source
            source[:] = 0.0
            source += half_ohmic_heat_membrane
            v_loss = np.minimum(self.e_0, cell.cathode.v_loss)
            v_loss[v_loss < 0.0] = 0.0
            reaction_heat = \
                (self.e_tn - self.e_0 + v_loss) * current
            source += reaction_heat
            mtx_func.add_explicit_layer_source(cell.thermal_rhs_dyn, source,
                                               cell.index_array, layer_id=2)

            # Anode gde-mem source
            source[:] = 0.0
            source += half_ohmic_heat_membrane
            v_loss = np.minimum(self.e_0, cell.anode.v_loss)
            v_loss[v_loss < 0.0] = 0.0
            reaction_heat = v_loss * current
            source += reaction_heat
            mtx_func.add_explicit_layer_source(cell.thermal_rhs_dyn, source,
                                               cell.index_array, layer_id=3)

            # Anode bpp-gde source
            source[:] = 0.0
            # h_vap = w_prop.water.calc_h_vap(cell.anode.temp_fluid[:-1])
            channel = cell.anode.channel
            source = channel.k_coeff * channel.temp_ele * self.n_ano_channels
            source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0
            mtx_func.add_explicit_layer_source(cell.thermal_rhs_dyn, source,
                                               cell.index_array, layer_id=4)

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

            # add thermal conductance for heat transfer to cathode gas
            source = -cell.cathode.channel.k_coeff * self.n_cat_channels
            matrix, source_vec_1 = mtx_func.add_implicit_layer_source(
                    cell.thermal_mtx_dyn, source, cell.index_array, layer_id=1)

            # add thermal conductance for heat transfer to anode gas
            source = -cell.anode.channel.k_coeff * self.n_ano_channels
            matrix, source_vec_2 = mtx_func.add_implicit_layer_source(
                cell.thermal_mtx_dyn, source, cell.index_array, 4)

            # add thermal conductance for heat transfer to coolant
            source_vec_3 = np.zeros(source_vec_1.shape)
            if self.cool_flow:
                if self.cool_ch_bc:
                    source = - self.cool_channels[i].k_coeff
                    source *= self.n_cool_sub_channels
                    matrix, source_vec = mtx_func.add_implicit_layer_source(
                        cell.thermal_mtx_dyn, source,
                        cell.index_array, layer_id=0)
                    source_vec_3[:] = source_vec
                    if cell.last_cell:
                        source = - self.cool_channels[i + 1].k_coeff
                        source *= self.n_cool_sub_channels
                        matrix, source_vec = mtx_func.add_implicit_layer_source(
                            cell.thermal_mtx_dyn, source,
                            cell.index_array, layer_id=-1)
                        source_vec_3[:] += source_vec
                else:
                    if not cell.first_cell:
                        source = - self.cool_channels[i - 1].k_coeff
                        source *= self.n_cool_sub_channels
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
        self.update_cell_solution_general(temp_solution_name, temp_shape_name)


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
            self.i_cd_tar = stack.i_cd_target
        else:
            self.v_tar = stack.v_target
            self.e_0_stack = np.sum([cell.e_0 for cell in self.cells])
            self.v_loss_tar = self.e_0_stack - self.v_tar

        # Current density of the elements in z-direction
        self.i_cd = np.zeros(self.shape)

        # Accumulated voltage loss over the stack at the lower end plate
        self.v_end_plate = np.zeros(self.shape)

    def update_rhs(self, *args, **kwargs):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        # raise NotImplementedError

        cell_0 = self.cells[0]
        if self.current_control:
            v_loss, v_loss_total = self.calc_voltage_loss()
            i_bc = v_loss[0] * cell_0.electrochemical_conductance
            i_target = self.i_cd_tar * cell_0.d_area
            i_correction_factor = i_target \
                                  / np.average(i_bc, weights=cell_0.d_area)
            v_loss_total *= - 1.0 * i_correction_factor
            self.rhs[:] = v_loss_total * cell_0.electrochemical_conductance
        else:
            self.rhs[:] = - self.v_loss_tar * cell_0.electrochemical_conductance

    def update_matrix(self, electrochemical_conductance, *args, **kwargs):
        """
        Updates matrix coefficients
        """
        # raise NotImplementedError
        mtx_list = \
            [mtx_func.build_cell_conductance_matrix([cond])
             for cond in electrochemical_conductance]
        mtx_dyn = sp_la.block_diag(*mtx_list)
        self.mtx[:] = self.mtx_const + mtx_dyn

    def update_cell_solution(self):
        raise NotImplementedError
        temp_solution_name = 'temp_layer'
        if not hasattr(self.cells[0], temp_solution_name):
            raise ValueError('attribute {} is not found in {}'.format(
                temp_solution_name, self.cells[0].__name__))
        temp_shape_name = 'temp_shape'
        if not hasattr(self.cells[0], temp_shape_name):
            raise ValueError('attribute {} is not found in {}'.format(
                temp_shape_name, self.cells[0].__name__))
        self.update_cell_solution_general(temp_solution_name, temp_shape_name)

    def update(self, current_density=None, voltage=None):
        """
        Coordinates the program sequence
        """
        # raise NotImplementedError
        # resistance = \
        #     np.asarray([cell.resistance for cell in self.cells])
        # self.resistance[:] = resistance.flatten()
        if current_density is not None:
            self.i_cd_tar = current_density
        if voltage is not None:
            self.v_tar = voltage
            self.v_loss_tar = self.e_0_stack - self.v_tar
        electrochemical_conductance = [
            cell.electrochemical_conductance for cell in self.cells]
        active_area = [cell.d_area for cell in self.cells]
        # if self.n_cells > 1:
        self.update_matrix(electrochemical_conductance)
        self.rhs[:self.n_ele] = self.calc_boundary_condition()
        self.i_cd[:] = self.calc_i(conductance_z, active_area)

        # else:
        #     i_bc = self.calc_boundary_condition()
        #     self.i_cd[:] = - i_bc / active_area
        #     v_diff = - i_bc / np.array([cell.electrochemical_conductance
        #                                 for cell in self.cells]).flatten()
        #     v_diff = v_diff.reshape((self.n_cells, self.n_ele))
        #     self.update_cell_voltage(v_diff)

        self.update_and_solve_linear_system()
        self.update_cell_solution()

    def calc_voltage_loss(self):
        v_loss = \
            np.asarray([np.average(cell.v_loss, weights=cell.d_area)
                        for cell in self.cells])
        v_loss_total = np.sum(v_loss)
        return v_loss, v_loss_total

    def calc_i(self, conductance, active_area):
        """
        Calculates the current density
        of the elements in z-direction.
        """
        raise NotImplementedError
        self.v_end_plate[:] = - self.rhs[:self.n_ele] / conductance[:self.n_ele]
        if self.solve_sparse:
            v_new = spsolve(self.mat, self.rhs)
            # mat_const = self.mat_const.toarray()
            # mat = self.mat.toarray()
        else:
            v_new = np.linalg.tensorsolve(self.mat, self.rhs)
        v_new = np.hstack((self.v_end_plate, v_new, np.zeros(self.n_ele)))
        v_diff = v_new[:-self.n_ele] - v_new[self.n_ele:]
        i_ca_vec = v_diff * conductance / active_area
        v_diff = v_diff.reshape((self.n_cells, self.n_ele))
        self.update_cell_voltage(v_diff)
        return np.reshape(i_ca_vec.flatten(), (self.n_cells, self.n_ele))

    def update_cell_voltage(self, v_diff):
        raise NotImplementedError
        for i, cell in enumerate(self.cells):
            cell.v[:] = cell.e_0 - v_diff[i]
            # cell.v_loss[:] = v_diff[i]
            cell.update_voltage_loss(v_diff[i])