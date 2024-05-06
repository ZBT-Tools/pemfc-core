# general imports
import numpy as np
import copy

# local modul imports
from . import interpolation as ip, matrix_functions as mtx, half_cell as h_c, \
    global_functions as g_func, membrane as membrane
from .output_object import OutputObject


class Cell(OutputObject):

    def __init__(self, cell_dict, membrane_dict, half_cell_dicts,
                 channels, number=None):
        name = 'Cell'  # + str(number)
        self.number = number
        super().__init__(name)
        self.cell_dict = cell_dict

        # print('Initializing: ', self.name)
        self.n_layer = 5
        self.n_electrodes = 2

        self.first_cell = cell_dict['first_cell']
        self.last_cell = cell_dict['last_cell']

        if self.last_cell:
            self.n_layer += 1
        n_nodes = channels[0].n_nodes
        # Number of nodes along the channel
        self.n_ele = n_nodes - 1

        # Additional resolution in width direction in regions
        # "under land" and "under channel" if True
        self.channel_land_discretization = cell_dict['channel_land_discretization']

        # Underrelaxation factor
        self.urf = cell_dict['underrelaxation_factor']

        self.e_0 = cell_dict['open_circuit_voltage']
        self.e_tn = cell_dict['thermoneutral_voltage']

        self.width = self.cell_dict['width']
        self.length = self.cell_dict['length']
        self.active_area = self.width * self.length

        half_cell_dicts = copy.deepcopy(half_cell_dicts)
        # Create half cell objects
        for i in range(len(half_cell_dicts)):
            # name = self.name + ': ' + half_cell_dicts[i]['name']
            name = half_cell_dicts[i]['name']
            half_cell_dicts[i]['name'] = name
        self.half_cells = [h_c.HalfCell(half_cell_dicts[i], cell_dict,
                                        channels[i], number=self.number)
                           for i in range(len(half_cell_dicts))]
        self.cathode = self.half_cells[0]
        self.anode = self.half_cells[1]

        # Append half cell names to output data
        for half_cell in self.half_cells:
            half_cell.channel.extend_data_names(half_cell.channel.name)
            half_cell.channel.fluid.extend_data_names(
                half_cell.channel.fluid.name)

        # self.dx = self.cathode.channel.dx

        """heat conductivity along and through the cell layers"""
        # self.width_straight_channels = \
        #     self.cathode.flow_field.width_straight_channels
        self.d_area = self.cathode.gde.dsct.d_area

        # Setup membrane
        self.membrane = membrane.Membrane(membrane_dict, self.cathode.gde.dsct)

        # Overall cell thickness (cell pitch)
        self.thickness = self.cathode.thickness + self.membrane.thickness \
            + self.anode.thickness

        # Cell coordinates in z-direction (stacking/current direction)
        # will be initialized correctly through stack class
        self.coords = [0.0, 0.0]

        # self.is_ht_pem = self.cell_dict['is_ht_pem']

        # heat conductivity along the gas diffusion electrode and membrane
        self.th_layer = \
            np.asarray([self.cathode.bpp.thickness,
                        self.cathode.gde.thickness,
                        self.membrane.thickness,
                        self.anode.gde.thickness,
                        self.anode.bpp.thickness])

        self.thermal_conductance_z = \
            np.asarray([self.cathode.bpp.thermal_conductance[0],
                        self.cathode.gde.thermal_conductance[0],
                        self.membrane.thermal_conductance[0],
                        self.anode.gde.thermal_conductance[0],
                        self.anode.bpp.thermal_conductance[0]])

        self.thermal_conductance_x = \
            np.asarray([self.cathode.bpp.thermal_conductance[1],
                        self.cathode.gde.thermal_conductance[1],
                        self.membrane.thermal_conductance[1],
                        self.anode.gde.thermal_conductance[1],
                        self.anode.bpp.thermal_conductance[1]])

        # self.thermal_conductance_y = \
        #     np.asarray([self.cathode.bpp.thermal_conductance[1],
        #                 self.cathode.gde.thermal_conductance[1],
        #                 self.membrane.thermal_conductance[1],
        #                 self.anode.gde.thermal_conductance[1],
        #                 self.anode.bpp.thermal_conductance[1]])

        self.thermal_conductance_x = \
            (self.thermal_conductance_x
             + np.roll(self.thermal_conductance_x, 1, axis=0)) * 0.5
        self.thermal_conductance_x = \
            np.vstack((self.thermal_conductance_x,
                       [self.thermal_conductance_x[0]]))
        if self.first_cell:
            self.thermal_conductance_x[0] *= 0.5
        if self.last_cell:
            self.thermal_conductance_x[-1] *= 0.5

        if self.last_cell:
            heat_cond_mtx = \
                mtx.build_cell_conductance_matrix(self.thermal_conductance_x,
                                                  self.thermal_conductance_z,
                                                  self.n_ele)
        else:
            heat_cond_mtx = \
                mtx.build_cell_conductance_matrix(
                    self.thermal_conductance_x[:-1],
                    self.thermal_conductance_z[:-1],
                    self.n_ele)
        self.heat_mtx_const = heat_cond_mtx
        # self.heat_mtx_const = np.zeros(self.heat_cond_mtx.shape)
        self.heat_mtx_dyn = np.zeros(self.heat_mtx_const.shape)
        self.heat_mtx = np.zeros(self.heat_mtx_dyn.shape)

        self.heat_rhs_const = np.zeros(self.n_layer * self.n_ele)
        self.heat_rhs_dyn = np.zeros(self.heat_rhs_const.shape)
        self.heat_rhs = np.zeros(self.heat_rhs_dyn.shape)

        # Create array for each thermal layer with indices according to
        # corresponding position in center diagonal of conductance matrix and
        # right hand side vector
        index_list = []
        for i in range(self.n_layer):
            index_list.append([(j * self.n_layer) + i
                               for j in range(self.n_ele)])
        self.index_array = np.asarray(index_list)

        # Set constant thermal boundary conditions
        if self.first_cell:
            end_plate_heat = cell_dict['heat_flux']
            heat_dx = end_plate_heat * self.cathode.discretization.d_area
            self.add_explicit_layer_source(self.heat_rhs_const, heat_dx, 0)
        if self.last_cell:
            end_plate_heat = cell_dict['heat_flux']
            heat_dx = end_plate_heat * self.anode.discretization.d_area
            self.add_explicit_layer_source(self.heat_rhs_const, heat_dx, -1)

        # Create electric conductance matrix
        self.elec_cond = \
            np.asarray([self.cathode.bpp.electrical_conductance[1],
                        self.anode.bpp.electrical_conductance[1]])
        self.elec_cond = \
            (self.elec_cond + np.roll(self.elec_cond, 1, axis=1)) * 0.5
        self.elec_cond = self.elec_cond[:, :-1]
        self.elec_x_mat_const = \
            mtx.build_z_cell_conductance_matrix(self.elec_cond.transpose())
        # print(self.elec_x_mat_const)

        # boolearn alarm values
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        # cell height
        self.height = self.membrane.thickness \
            + self.cathode.bpp.thickness \
            + self.cathode.gde.thickness \
            + self.anode.bpp.thickness \
            + self.anode.gde.thickness

        # Initializing temperatures with average channel fluid temperature
        temp_init = np.average([hc.channel.fluid.temperature
                                for hc in self.half_cells])
        # membrane temperature
        self.temp_mem = np.zeros(self.membrane.temp.shape)
        self.temp_layer = \
            g_func.full((self.n_layer,) + self.temp_mem.shape, temp_init)
        # interface names according to temperature array
        self.temp_names = ['Cathode BPP-BPP',
                           'Cathode BPP-GDE',
                           'Cathode GDE-MEM',
                           'Anode MEM-GDE',
                           'Anode GDE-BPP',
                           'Anode BPP-BPP']

        # current density
        self.i_cd = np.zeros(self.temp_mem.shape)
        # cell voltage
        self.v = np.zeros(self.i_cd.shape)
        # voltage loss
        self.v_loss = np.zeros(self.v.shape)
        # self.resistance_z = np.zeros(n_ele)
        # through-plane cell resistance
        self.conductance_z = np.zeros(self.i_cd.shape)

        self.add_print_data(self.i_cd, 'Current Density', 'A/m²')
        self.add_print_data(self.temp_layer, 'Temperature', 'K',
                            self.temp_names[:self.n_layer])
        self.add_print_data(self.v, 'Cell Voltage', 'V')

    def calc_ambient_conductance(self, alpha_amb):
        """
        :param alpha_amb: heat transfer coefficient for free or forced
        convection to the ambient
        :return: discretized conductance for each layer of cell based on
        external surface
        """
        # geometry model
        # ext_surface_factor = (self.width + self.length) \
        #     / (self.cathode.channel.length * self.width_straight_channels)
        # k_amb = np.full((self.n_layer, self.n_ele), 0.)
        # convection conductance to the environment
        # ToDo: Check thickness calculation
        th_layer_amb = (self.th_layer + np.roll(self.th_layer, 1)) * 0.5
        # if self.last_cell:
        th_layer_amb = np.hstack((th_layer_amb, th_layer_amb[0]))
        # ToDo: Check conductance calculation
        k_amb = np.outer(th_layer_amb, self.cathode.discretization.dx[0]) \
            * alpha_amb * self.cathode.flow_field.external_surface_factor
        if self.first_cell:
            k_amb[0] *= 0.5
        if self.last_cell:
            k_amb[-1] *= 0.5
        return k_amb

    def add_explicit_layer_source(self, rhs_vector, source_term,
                                  layer_id=None):
        if layer_id is None:
            if np.isscalar(source_term):
                source_vector = np.full_like(rhs_vector, -source_term)
            else:
                source_vector = np.asarray(-source_term)
        else:
            source_vector = np.zeros(rhs_vector.shape)
            np.put(source_vector, self.index_array[layer_id], -source_term)
        rhs_vector += source_vector
        return rhs_vector, source_vector

    def add_implicit_layer_source(self, matrix, coefficient, layer_id=None):
        matrix_size = matrix.shape[0]
        if layer_id is None:
            if np.isscalar(coefficient):
                source_vector = g_func.full(matrix_size, coefficient)
            else:
                source_vector = np.asarray(coefficient)
        else:
            source_vector = np.zeros(matrix_size)
            np.put(source_vector, self.index_array[layer_id], coefficient)
        matrix += np.diag(source_vector)
        return matrix, source_vector

    def update(self, current_density, update_channel=False,
               current_control=True, urf=None):
        """
        This function coordinates the program sequence
        """
        if urf is None:
            urf = self.urf
        current_density = (1.0 - urf) * current_density + urf * self.i_cd
        # if g_par.iteration > 50:
        #     self.urf *= 0.99
        # self.urf = max(self.urf, 0.8)
        # self.temp_mem[:] = .5 * (self.temp_layer[2] + self.temp_layer[3])
        self.membrane.temp[:] = 0.5 * (self.temp_layer[2] + self.temp_layer[3])
        if isinstance(self.membrane, membrane.WaterTransportMembrane):
            self.cathode.w_cross_flow[:] = self.membrane.water_flux * -1.0
            self.anode.w_cross_flow[:] = self.membrane.water_flux
        # self.cathode.set_layer_temperature([self.temp[2], self.temp[3],
        #                                     self.temp[4]])
        # self.anode.set_layer_temperature([self.temp[0], self.temp[1]])
        self.cathode.update(current_density, update_channel=update_channel,
                            current_control=current_control)
        self.anode.update(current_density, update_channel=update_channel,
                          current_control=True)
        if self.cathode.electrochemistry.corrected_current_density is not None:
            corrected_current_density = \
                self.cathode.electrochemistry.corrected_current_density
        else:
            corrected_current_density = current_density
        if self.anode.break_program or self.cathode.break_program:
            self.break_program = True
        else:
            humidity = np.asarray([self.cathode.channel.fluid.humidity,
                                   self.anode.channel.fluid.humidity])
            humidity_ele = \
                np.array([ip.interpolate_1d(humidity[0]),
                          ip.interpolate_1d(humidity[1])])
            self.membrane.update(corrected_current_density, humidity_ele)
            self.calc_voltage_loss()
            self.calc_conductance(corrected_current_density)
            # if np.any(self.v_alarm) and current_control:
            #     self.correct_voltage_loss()
                # raise ValueError('voltage losses greater than '
                #                  'open circuit voltage')
            self.i_cd[:] = current_density

    def calc_voltage_loss(self):
        """
        Calculates the cell voltage loss. If the cell voltage loss is larger
        than the open circuit cell voltage, the cell voltage is set to zero.
        """
        self.v_loss[:] = \
            self.membrane.v_loss + self.cathode.v_loss + self.anode.v_loss
        self.v_alarm = self.v_loss >= self.e_0
        # self.v_loss[:] = np.minimum(self.v_loss, self.e_0)

    def update_voltage_loss(self, v_loss):
        v_loss_factor = v_loss / self.v_loss
        self.v_loss[:] *= v_loss_factor
        for electrode in self.half_cells:
            electrode.v_loss[:] *= v_loss_factor
        self.membrane.v_loss[:] *= v_loss_factor

    # def correct_voltage_loss(self):
    #     try:
    #         id = np.argwhere(np.invert(self.v_alarm))[-1, -1]
    #     except IndexError:
    #         id = 0
    #     np.seterr(divide='ignore')
    #     # gdl_loss_ratio = \
    #     #     self.cathode.gdl_diff_loss[id] / self.anode.gdl_diff_loss[id]
    #     cl_loss_ratio = \
    #         np.where(self.anode.v_loss_cl_diff[id] == 0.0, 0.0,
    #                  self.cathode.v_loss_cl_diff[id]
    #                  / self.anode.v_loss_cl_diff[id])
    #     cat_loss_ratio = \
    #         np.where(self.cathode.v_loss_cl_diff[id] == 0.0, 0.0,
    #                  self.cathode.v_loss_gdl_diff[id]
    #                  / self.cathode.v_loss_cl_diff[id])
    #     ano_loss_ratio = \
    #         np.where(self.anode.v_loss_cl_diff[id] == 0.0, 0.0,
    #                  self.anode.v_loss_gdl_diff[id]
    #                  / self.anode.v_loss_cl_diff[id])
    #     v_loss_max = self.e_0 - 0.05
    #     v_loss_diff_max = v_loss_max - self.membrane.v_loss[id] \
    #                       - self.cathode.v_loss_bpp[id] - self.cathode.v_loss_act[id] \
    #                       - self.anode.v_loss_bpp[id] - self.anode.v_loss_act[id]
    #     v_loss_cat_cl = v_loss_diff_max * cl_loss_ratio \
    #         / (1.0 + ano_loss_ratio
    #            + cat_loss_ratio * cl_loss_ratio + cl_loss_ratio)
    #     v_loss_ano_cl = \
    #         np.where(cl_loss_ratio == 0.0, 0.0,
    #                  v_loss_cat_cl / cl_loss_ratio)
    #     v_loss_cat_gdl = v_loss_cat_cl * cat_loss_ratio
    #     v_loss_ano_gdl = v_loss_ano_cl * ano_loss_ratio
    #     np.seterr(divide='raise')
    #     self.cathode.v_loss_gdl_diff[:] = \
    #         np.where(self.v_alarm, v_loss_cat_gdl, self.cathode.v_loss_gdl_diff)
    #     self.cathode.v_loss_cl_diff[:] = \
    #         np.where(self.v_alarm, v_loss_cat_cl, self.cathode.v_loss_cl_diff)
    #     self.anode.v_loss_gdl_diff[:] = \
    #         np.where(self.v_alarm, v_loss_ano_gdl, self.anode.v_loss_gdl_diff)
    #     self.anode.v_loss_cl_diff[:] = \
    #         np.where(self.v_alarm, v_loss_ano_cl, self.anode.v_loss_cl_diff)
    #
    #     self.cathode.v_loss[:] = \
    #         self.cathode.v_loss_act + self.cathode.v_loss_cl_diff \
    #         + self.cathode.v_loss_gdl_diff + self.cathode.v_loss_bpp
    #     self.anode.v_loss[:] = \
    #         self.anode.v_loss_act + self.anode.v_loss_cl_diff \
    #         + self.anode.v_loss_gdl_diff + self.anode.v_loss_bpp
    #     self.calc_voltage_loss()
    #     # raise ValueError

    def calc_conductance(self, current_density):
        """
        Calculates the area-specific electrical resistance of the element in
        z-direction
        """
        current = current_density * self.membrane.dsct.d_area
        resistance_z = self.v_loss / current
        self.conductance_z[:] = 1.0 / resistance_z
