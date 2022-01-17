# general imports
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

# local module imports
from . import matrix_functions as mtx


class ElectricalCoupling:

    def __init__(self, stack):
        self.stack = stack
        self.cells = stack.cells
        # Handover
        self.n_cells = self.stack.n_cells
        # number of the stack cells
        self.dx = self.cells[0].dx
        # length of an element
        self.th_plate = self.cells[0].cathode.bpp.thickness
        # thickness of the bipolar plate

        # Variables
        self.current_control = stack.current_control
        if self.current_control:
            self.i_cd_tar = stack.i_cd_target
        else:
            self.v_tar = stack.v_target
            self.e_0_stack = np.sum([cell.e_0 for cell in self.cells])
            self.v_loss_tar = self.e_0_stack - self.v_tar
        # number of the nodes along the channel
        self.n_ele = self.cells[0].n_ele
        # number of the elements along the channel
        self.i_cd = np.zeros((self.n_cells, self.n_ele))
        # current density of the elements in z-direction
        # self.resistance = np.zeros((self.n_cells, self.n_ele)).flatten()
        # combined cell & bipolar plate resistance vector in z-direction
        self.v_end_plate = np.zeros(self.n_ele)
        # accumulated voltage loss over the stack at the lower end plate
        if self.n_cells > 1:
            self.mat = None
            # electrical conductance matrix
            self.rhs = np.zeros((self.n_cells - 1) * self.n_ele)
            # right hand side terms, here the current
            # self.c_width = \
            #     self.cells[0].cathode.rib_width \
            #     * (self.cells[0].cathode.n_channel + 1)
            # self.cells[0].width_straight_channels
            # width of the channel

            self.solve_sparse = True
            cell_mat_x_list = [cell.elec_x_mat_const for cell in self.cells]

            self.mat_const = mtx.block_diag_overlap(cell_mat_x_list,
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
            np.asarray([cell.conductance_z for cell in self.cells]).flatten()
        # conductance = (self.c_width * self.dx / resistance).flatten()
        # conductance = 1.0 / self.resistance
        active_area = \
            np.array([cell.active_area_dx for cell in self.cells]).flatten()
        if self.n_cells > 1:
            self.update_mat(conductance_z)
            self.rhs[:self.n_ele] = self.calc_boundary_condition()
            self.i_cd[:] = self.calc_i(conductance_z, active_area)

        else:
            i_bc = self.calc_boundary_condition()
            self.i_cd[:] = - i_bc / active_area
            v_diff = - i_bc / np.array([cell.conductance_z
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
            np.asarray([np.average(cell.v_loss, weights=cell.active_area_dx)
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
            i_bc = v_loss[0] * cell_0.conductance_z
            i_target = self.i_cd_tar * cell_0.active_area_dx
            i_correction_factor = i_target \
                / np.average(i_bc, weights=cell_0.active_area_dx)
            v_loss_total *= - 1.0 * i_correction_factor
            return v_loss_total * cell_0.conductance_z
        else:
            return - self.v_loss_tar * cell_0.conductance_z

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
