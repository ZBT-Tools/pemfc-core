# general imports
import numpy as np
from scipy import linalg as sp_la
from scipy import sparse
from scipy.sparse.linalg import spsolve

# local module imports
from . import matrix_functions as mtx


class ElectricalCoupling:

    def __init__(self, stack):
        # Reference to stack object
        self.stack = stack
        # Reference to list of cells
        self.cells = stack.cells
        # Number of the stack cells
        self.n_cells = self.stack.n_cells

        # Variables
        self.current_control = stack.current_control
        if self.current_control:
            self.i_cd_tar = stack.i_cd_target
        else:
            self.v_tar = stack.v_target
            self.e_0_stack = np.sum([cell.e_0 for cell in self.cells])
            self.v_loss_tar = self.e_0_stack - self.v_tar

        # TODO: Update ElectricalCoupling for 3D
        # number of the nodes along the channel
        # self.n_ele = self.cells[0].n_ele
        self.shape = self.cells[0].membrane.dsct.shape
        # number of the elements along the channel
        self.i_cd = np.zeros((self.n_cells, *self.shape))
        # current density of the elements in z-direction
        # self.resistance = np.zeros((self.n_cells, self.n_ele)).flatten()
        # combined cell & bipolar plate resistance vector in z-direction
        self.v_end_plate = np.zeros(self.shape)
        # accumulated voltage loss over the stack at the lower end plate
        if self.n_cells > 1:
            self.mat = None
            # electrical conductance matrix
            self.rhs = np.zeros((self.n_cells - 1) * np.prod(self.shape))
            # right hand side terms, here the current
            # self.c_width = \
            #     self.cells[0].cathode.rib_width \
            #     * (self.cells[0].cathode.n_channel + 1)
            # self.cells[0].width_straight_channels
            # width of the channel

            self.solve_sparse = True

            cell_mat_y_list = [cell.elec_mtx_const for cell in self.cells]

            # TODO: Update stack conductance matrix: in 3D the overlap does not seem suitable,
            #  just additional x-conductance between bipolar half plates seems correct

            stack_y_elec_cond_mat = sp_la.block_diag(*cell_mat_y_list)
            self.mat_const = mtx.block_diag_overlap(cell_mat_y_list,
                                                    (self.n_ele, self.n_ele))
            self.mat_const = \
                self.mat_const[self.n_ele:-self.n_ele, self.n_ele:-self.n_ele]
            if self.solve_sparse:
                self.mat_const = sparse.csr_matrix(self.mat_const)

    def update(self, current_density=None, voltage=None):
        """
        Coordinates the program sequence
        """
        # resistance = \
        #     np.asarray([cell.resistance for cell in self.cells])
        # self.resistance[:] = resistance.flatten()
        if current_density is not None:
            self.i_cd_tar = current_density
        if voltage is not None:
            self.v_tar = voltage
            self.v_loss_tar = self.e_0_stack - self.v_tar
        conductance_z = \
            np.asarray([cell.electrochemical_conductance for cell in self.cells]).flatten()
        # conductance = (self.c_width * self.dx / resistance).flatten()
        # conductance = 1.0 / self.resistance
        active_area = \
            np.array([cell.d_area for cell in self.cells]).flatten()
        if self.n_cells > 1:
            self.update_mat(conductance_z)
            self.rhs[:self.n_ele] = self.calc_boundary_condition()
            self.i_cd[:] = self.calc_i(conductance_z, active_area)

        else:
            i_bc = self.calc_boundary_condition()
            self.i_cd[:] = - i_bc / active_area
            v_diff = - i_bc / np.array([cell.electrochemical_conductance
                                        for cell in self.cells]).flatten()
            v_diff = v_diff.reshape((self.n_cells, self.n_ele))
            self.update_cell_voltage(v_diff)

    def update_mat(self, conductance):
        """
        Updates the conductance matrix
        """
        cell_c_mid = \
            np.hstack((conductance[:-self.n_ele] + conductance[self.n_ele:]))
        mat_dyn = \
            - np.diag(cell_c_mid, 0) \
            + np.diag(conductance[:-self.n_ele][self.n_ele:], self.n_ele) \
            + np.diag(conductance[:-self.n_ele][self.n_ele:], -self.n_ele)
        if self.solve_sparse:
            mat_dyn = sparse.csr_matrix(mat_dyn)
        self.mat = self.mat_const + mat_dyn

    def calc_voltage_loss(self):
        v_loss = \
            np.asarray([np.average(cell.v_loss, weights=cell.d_area)
                        for cell in self.cells])
        v_loss_total = np.sum(v_loss)
        return v_loss, v_loss_total

    def calc_boundary_condition(self):
        """
        Updates the right hand side of the linear src.
        """
        cell_0 = self.cells[0]
        if self.current_control:
            v_loss, v_loss_total = self.calc_voltage_loss()
            i_bc = v_loss[0] * cell_0.electrochemical_conductance
            i_target = self.i_cd_tar * cell_0.d_area
            i_correction_factor = i_target \
                / np.average(i_bc, weights=cell_0.d_area)
            v_loss_total *= - 1.0 * i_correction_factor
            return v_loss_total * cell_0.electrochemical_conductance
        else:
            return - self.v_loss_tar * cell_0.electrochemical_conductance

    def calc_i(self, conductance, active_area):
        """
        Calculates the current density
        of the elements in z-direction.
        """
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
        for i, cell in enumerate(self.cells):
            cell.v[:] = cell.e_0 - v_diff[i]
            # cell.v_loss[:] = v_diff[i]
            cell.update_voltage_loss(v_diff[i])
