# general imports
import numpy as np
import copy
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve

# local module imports
from . import matrix_functions as mtx, cell as fcell, \
    global_functions as g_func, channel as chl

# import pandas as pd
# from numba import jit

np.set_printoptions(linewidth=10000, threshold=None, precision=2)


class TemperatureSystem:

    def __init__(self, stack, temp_dict):
        self.cells = stack.cells
        if not isinstance(stack.cells, (list, tuple)):
            raise TypeError
        if not isinstance(self.cells[0], fcell.Cell):
            raise TypeError
        self.n_cells = stack.n_cells

        # use SciPy sparse solver, efficient for larger sparse matrices
        self.sparse_solve = True

        # instead of solving the completely coupled temperature src at once
        # solve the decoupled cell-wise temperature systems and iterate
        # however not working yet!!!
        self.solve_individual_cells = False

        # sub channel ratios
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
            if self.n_cool == (self.n_cells + 1):
                self.cool_ch_bc = True
            else:
                self.cool_ch_bc = False
            self.n_cool_sub_channels = stack.coolant_circuit.n_subchannels

        # thermodynamic neutral cell potential
        self.e_tn = self.cells[0].e_tn
        # open circuit potential
        self.e_0 = self.cells[0].e_0
        # conductance matrix
        self.mtx = None
        # solid layer temperature vector
        self.temp_layer_vec = \
            np.hstack([cell.heat_rhs for cell in self.cells])
        # right side of the matrix src: mat T = rhs,
        # contains the power sources and explicit coupled terms
        self.rhs = np.zeros(np.shape(self.temp_layer_vec))

        """Building up the base conductance matrix"""
        # vector for dynamically changing heat conductance values
        # self.dyn_vec = np.zeros(np.shape(self.temp_layer_vec))
        # Add constant source and sink coefficients to heat conductance matrix
        # Heat transfer to ambient
        alpha_amb = temp_dict['alpha_amb']
        temp_amb = temp_dict['temp_amb']
        # ambient temperature
        for cell in self.cells:
            cell.k_amb = cell.calc_ambient_conductance(alpha_amb)
            if cell.last_cell:
                k_amb_vector = cell.k_amb.transpose().flatten()
            else:
                k_amb_vector = cell.k_amb[:-1].transpose().flatten()

            cell.add_implicit_layer_source(cell.heat_mtx_const, -k_amb_vector)
            cell.add_explicit_layer_source(cell.heat_rhs_const,
                                           k_amb_vector * temp_amb)

        self.rhs_const = \
            np.hstack([cell.heat_rhs_const for cell in self.cells])

        self.index_list, self.layer_index_list = \
            mtx.create_index_lists(self.cells)

        # constant part of conductance matrix
        self.mtx_const = self.connect_cells()
        if self.sparse_solve:
            self.mtx_const = sparse.csr_matrix(self.mtx_const)

    def connect_cells(self):
        matrix = sp_la.block_diag(*[cell.heat_mtx_const for cell in self.cells])
        cell_ids = np.asarray([list(range(self.n_cells-1)),
                               list(range(1, self.n_cells))]).transpose()
        layer_ids = np.asarray([(-1, 0) for i in range(self.n_cells-1)])
        conductance = \
            np.asarray([self.cells[i].thermal_conductance_x[layer_ids[i][0]]
                        for i in range(self.n_cells-1)])
        mtx.connect_cells(matrix, cell_ids, layer_ids,
                          conductance, self.index_list)
        return matrix

    def update(self):
        """
        This function coordinates the program sequence
        """
        self.update_gas_channel()
        if self.cool_flow:
            self.update_coolant_channel()
        self.update_temp_layer()
        if np.any(self.temp_layer_vec < 200.0):
            raise ValueError('temperature too low, check boundary conditions')
        if np.any(self.temp_layer_vec > 1000.0):
            raise ValueError('temperature too high, check boundary conditions')

    def update_temp_layer(self):
        """
        This function coordinates the temp_layer program sequence
        """
        self.update_matrix()
        self.update_rhs()
        self.solve_system()
        self.update_cell_layer_temperatures()

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
        Creates a vector with the right hand side entries,
        add explicit heat sources here.
        Sources from outside the src
        to the src must be defined negative.
        """
        for i, cell in enumerate(self.cells):
            cell.heat_rhs_dyn[:] = 0.0

            source = np.zeros(cell.temp_layer[0].shape)
            # Cathode bpp-gde source
            # h_vap = w_prop.water.calc_h_vap(cell.cathode.channel.temp[:-1])
            channel = cell.cathode.channel
            source += channel.k_coeff * channel.temp_ele * self.n_cat_channels
            source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 1)

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
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 2)

            # Anode gde-mem source
            source[:] = 0.0
            source += half_ohmic_heat_membrane
            v_loss = np.minimum(self.e_0, cell.anode.v_loss)
            v_loss[v_loss < 0.0] = 0.0
            reaction_heat = v_loss * current
            source += reaction_heat
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 3)

            # Anode bpp-gde source
            source[:] = 0.0
            # h_vap = w_prop.water.calc_h_vap(cell.anode.temp_fluid[:-1])
            channel = cell.anode.channel
            source = channel.k_coeff * channel.temp_ele * self.n_ano_channels
            source += getattr(channel, 'condensation_heat', 0.0)  # * 0.0
            cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 4)

            # Cooling channels
            if self.cool_flow:
                if self.cool_ch_bc:
                    cool_chl = self.cool_channels[i]
                    source = cool_chl.k_coeff * cool_chl.temp_ele
                    source *= self.n_cool_sub_channels
                    cell.add_explicit_layer_source(cell.heat_rhs_dyn,
                                                   source, layer_id=0)
                    if cell.last_cell:
                        cool_chl = self.cool_channels[i + 1]
                        source = cool_chl.k_coeff * cool_chl.temp_ele
                        source *= self.n_cool_sub_channels
                        cell.add_explicit_layer_source(cell.heat_rhs_dyn,
                                                       source, layer_id=-1)
                else:
                    if not cell.first_cell:
                        cool_chl = self.cool_channels[i - 1]
                        source = cool_chl.k_coeff * cool_chl.temp_ele
                        source *= self.n_cool_sub_channels
                        cell.add_explicit_layer_source(cell.heat_rhs_dyn,
                                                       source, layer_id=0)

        if not self.solve_individual_cells:
            rhs_dyn = np.hstack([cell.heat_rhs_dyn for cell in self.cells])
            self.rhs = self.rhs_const + rhs_dyn

    def update_matrix(self):
        """
        Updates the thermal conductance matrix
        """
        source_vectors = []
        for i, cell in enumerate(self.cells):
            cell.heat_mtx_dyn[:, :] = 0.0
            source_vectors.append(np.zeros(cell.heat_rhs_dyn.shape))

            # add thermal conductance for heat transfer to cathode gas
            source = -cell.cathode.channel.k_coeff * self.n_cat_channels
            matrix, source_vec_1 = \
                cell.add_implicit_layer_source(cell.heat_mtx_dyn, source, 1)

            # add thermal conductance for heat transfer to anode gas
            source = -cell.anode.channel.k_coeff * self.n_ano_channels
            matrix, source_vec_2 = \
                cell.add_implicit_layer_source(cell.heat_mtx_dyn, source, 4)

            # add thermal conductance for heat transfer to coolant
            source_vec_3 = np.zeros(source_vec_1.shape)
            if self.cool_flow:
                if self.cool_ch_bc:
                    source = - self.cool_channels[i].k_coeff
                    source *= self.n_cool_sub_channels
                    matrix, source_vec = \
                        cell.add_implicit_layer_source(cell.heat_mtx_dyn,
                                                       source, layer_id=0)
                    source_vec_3[:] = source_vec
                    if cell.last_cell:
                        source = - self.cool_channels[i + 1].k_coeff
                        source *= self.n_cool_sub_channels
                        matrix, source_vec = \
                            cell.add_implicit_layer_source(cell.heat_mtx_dyn,
                                                           source, layer_id=-1)
                        source_vec_3[:] += source_vec
                else:
                    if not cell.first_cell:
                        source = - self.cool_channels[i - 1].k_coeff
                        source *= self.n_cool_sub_channels
                        matrix, source_vec = \
                            cell.add_implicit_layer_source(cell.heat_mtx_dyn,
                                                           source, layer_id=0)
                        source_vec_3[:] = source_vec

            source_vectors[i][:] = source_vec_1 + source_vec_2 + source_vec_3

        dyn_vec = np.hstack(source_vectors)
        if self.sparse_solve:
            self.mtx = \
                self.mtx_const + sparse.diags([dyn_vec], [0], format='csr')
        else:
            self.mtx = self.mtx_const + np.diag(dyn_vec)

    def solve_system(self):
        if self.solve_individual_cells:
            self.solve_cells()
        else:
            self.solve_implicit_system()

    def solve_implicit_system(self):
        """
        Solves the layer temperatures.
        """
        # tx_df = pd.DataFrame(self.mtx)
        # mtx_df.to_clipboard(index=False, header=False, sep=' ')
        # rhs_df = pd.DataFrame(self.rhs)
        # rhs_df.to_clipboard(index=False, header=False, sep=' ')
        if self.sparse_solve:
            self.temp_layer_vec[:] = spsolve(self.mtx, self.rhs)
        else:
            self.temp_layer_vec[:] = np.linalg.tensorsolve(self.mtx, self.rhs)

    def connect_to_next_cell(self, i):
        cell = self.cells[i]
        conductance = cell.thermal_conductance_x[-1]
        source = conductance * self.cells[i + 1].temp_layer[0]
        cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, -1)
        coeff = - conductance
        cell.add_implicit_layer_source(cell.heat_mtx_dyn, coeff, -1)

    def connect_to_previous_cell(self, i):
        cell = self.cells[i]
        conductance = self.cells[i - 1].thermal_conductance_x[-1]
        source = - conductance * self.cells[i - 1].temp_layer[-1]
        # source = conductance * (self.cells[i - 1].temp_layer[-1] -
        #                        self.cells[i].temp_layer[0])
        cell.add_explicit_layer_source(cell.heat_rhs_dyn, source, 0)
        source = conductance
        cell.add_implicit_layer_source(cell.heat_mtx_dyn, source, 0)

    def solve_cells(self):
        """
        Solves the layer temperatures solving each cell individually and
        explicit coupling between cells (self.solve_individual_cells = True),
        however not working properly at the moment!!!
        """
        tolerance = 1e-10
        error = 1e8
        under_relaxation = 0.0
        temp_old = g_func.full(self.temp_layer_vec.shape, 1e8)
        temp_new = np.zeros(self.temp_layer_vec.shape)
        counter = 0
        while (error > tolerance) or (counter < 3):
            temp_new_list = []
            for i, cell in enumerate(self.cells):
                if counter > 0:
                    if self.n_cells > 1:
                        if cell.first_cell:
                            self.connect_to_next_cell(i)
                        elif cell.last_cell:
                            self.connect_to_previous_cell(i)
                        else:
                            self.connect_to_previous_cell(i)
                            self.connect_to_next_cell(i)
                cell.heat_mtx[:] = cell.heat_mtx_const + cell.heat_mtx_dyn
                cell.heat_rhs[:] = cell.heat_rhs_const + cell.heat_rhs_dyn
                temp_layer_vec = \
                    np.linalg.tensorsolve(cell.heat_mtx, cell.heat_rhs)
                if counter > 0:
                    temp_layer_vec = (1.0 - under_relaxation) * temp_layer_vec \
                        + under_relaxation * temp_old_array[i]
                temp_new_list.append(temp_layer_vec)
                for j in range(cell.n_layer):
                    index_vector = cell.index_array[j]
                    cell.temp_layer[j] = temp_layer_vec[index_vector]
            temp_old_array = copy.deepcopy(temp_new_list)
            temp_new[:] = np.concatenate(temp_new_list, axis=0)
            error = np.abs(np.sum(((temp_old - temp_new) / temp_new) ** 2.0))
            temp_old[:] = temp_new
            counter += 1
        self.temp_layer_vec[:] = temp_new

    def update_cell_layer_temperatures(self):
        """
        From 1D temperature vector to 2D cell temperature arrays
        """
        for i, cell in enumerate(self.cells):
            for j in range(cell.n_layer):
                # index_vector = cell.index_array[j]
                index_vector = self.index_list[i][j]
                cell.temp_layer[j] = self.temp_layer_vec[index_vector]

