# general imports
import numpy as np
import copy
from abc import ABC, abstractmethod
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve

# local module imports
from . import matrix_functions as mtx_func, stack as stack_module, cell as cell_module, \
    global_functions as g_func, channel as chl

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
        self.n = np.prod(*shape)

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

    def update_and_solve_linear_system(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
        self.solve_system()

    @abstractmethod
    def update_rhs(self):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        pass

    @abstractmethod
    def update_matrix(self):
        """
        Updates matrix coefficients
        """
        pass

    def solve_system(self):
        """
        Solves the layer temperatures.
        """
        if self.sparse_solve:
            self.solution_vector[:] = spsolve(self.mtx, self.rhs)
        else:
            self.solution_vector[:] = np.linalg.tensorsolve(self.mtx, self.rhs)


class StackLinearSystem(LinearSystem, ABC):

    def __init__(self, shape: tuple[int, ...], stack: stack_module.Stack):
        super().__init__(shape)
        self.cells = stack.cells

        self.index_list, self.layer_index_list = mtx_func.create_stack_index_list(self.cells)

        # constant part of conductance matrix
        self.mtx_const = self.connect_cells()
        if self.sparse_solve:
            self.mtx_const = sparse.csr_matrix(self.mtx_const)

    @abstractmethod
    def connect_cells(self):
        pass

    @abstractmethod
    def update(self):
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

        self.rhs_const = np.hstack([cell.heat_rhs_const for cell in self.cells])

        # TODO: Update TemperatureSystem for 3D

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
        # Add coefficients to matrix and rhs for each cell due to convection to ambient
        alpha_amb = self.input_dict['alpha_amb']
        temp_amb = self.input_dict['temp_amb']
        for cell in self.cells:
            cell.k_amb = cell.calc_ambient_conductance(alpha_amb)
            # if cell.last_cell:
            k_amb_vector = cell.k_amb.flatten(order='F')
            # else:
            #     k_amb_vector = cell.k_amb[:-1].flatten(order='F')

            mtx_func.add_implicit_layer_source(
                cell.heat_mtx_const, -k_amb_vector, cell.index_array)
            mtx_func.add_explicit_layer_source(
                cell.heat_rhs_const, k_amb_vector * temp_amb, cell.index_array)

    def connect_cells(self):
        matrix = sp_la.block_diag(*[cell.heat_mtx_const for cell in self.cells])
        n_cells = len(self.cells)
        cell_ids = np.asarray([list(range(n_cells-1)),
                               list(range(1, n_cells))]).transpose()
        layer_ids = np.asarray([(-1, 0) for i in range(n_cells-1)])
        conductance = \
            np.asarray([self.cells[i].thermal_conductance[0][layer_ids[i][0]]
                        for i in range(n_cells-1)])
        # old_matrix = np.copy(matrix)
        mtx_func.connect_cells(matrix, cell_ids, layer_ids,
                               conductance, self.index_list)
        # diff = matrix - old_matrix
        return matrix

    def update(self):
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

    def update_and_solve_linear_system(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
        self.solve_system()

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

    def update_rhs(self):
        """
        Create vector with the right hand side entries,
        Sources from outside the src to the src must be defined negative.
        """
        pass

    def update_matrix(self):
        """
        Updates matrix coefficients
        """
        pass

    def update_cell_solution(self):
        temp_solution_name = 'temp_layer'
        if not hasattr(self.cells[0], temp_solution_name):
            raise ValueError('attribute {} is not found in {}'.format(temp_solution_name,
                                                                      self.cells[0].__name__))
        temp_shape_name = 'temp_shape'
        if not hasattr(self.cells[0], temp_shape_name):
            raise ValueError('attribute {} is not found in {}'.format(temp_shape_name,
                                                                      self.cells[0].__name__))
        self.update_cell_solution_general(temp_solution_name, temp_shape_name)
