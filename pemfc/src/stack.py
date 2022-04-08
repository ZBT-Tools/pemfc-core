# general imports
import numpy as np

# local module imports
from . import electrical_coupling as el_cpl, flow_circuit as flow_circuit, \
    cell as cl, temperature_system as therm_cpl, fluid as fluid, channel as chl
# from data import input_dicts
# from ..gui import data_transfer

# gui_data = True


class Stack:

    def __init__(self, settings, n_nodes, current_control=False):

        # Read settings dictionaries
        stack_dict = settings['stack']

        self.n_cells = stack_dict['cell_number']
        # number of cells of the stack
        n_ele = n_nodes - 1
        # node points/elements along the x-axis
        self.calc_temp = stack_dict['calc_temperature']
        # switch to calculate the temperature distribution
        self.calc_electric = stack_dict['calc_current_density']
        # switch to calculate the current density distribution
        # self.calc_flow_dis = stack_dict['calc_flow_distribution']
        # switch to calculate the flow distribution

        # decompose input dict the individual objects
        cell_dict = settings['cell']
        membrane_dict = settings['membrane']
        anode_dict = settings['anode']
        cathode_dict = settings['cathode']
        ano_channel_dict = anode_dict['channel']
        cat_channel_dict = cathode_dict['channel']
        ano_fluid_dict = ano_channel_dict['fluid']
        cat_fluid_dict = cat_channel_dict['fluid']
        ano_flow_circuit_dict = anode_dict['flow_circuit']
        ano_in_manifold_dict = ano_flow_circuit_dict['inlet_manifold']
        ano_out_manifold_dict = ano_flow_circuit_dict['outlet_manifold']
        cat_flow_circuit_dict = cathode_dict['flow_circuit']
        cat_in_manifold_dict = cat_flow_circuit_dict['inlet_manifold']
        cat_out_manifold_dict = cat_flow_circuit_dict['outlet_manifold']
        temperature_dict = settings['temperature_system']

        half_cell_dicts = [cathode_dict, anode_dict]
        channel_dicts = [cat_channel_dict, ano_channel_dict]
        fluid_dicts = [cat_fluid_dict, ano_fluid_dict]
        manifold_in_dicts = [cat_in_manifold_dict, ano_in_manifold_dict]
        manifold_out_dicts = [cat_out_manifold_dict, ano_out_manifold_dict]
        flow_circuit_dicts = [cat_flow_circuit_dict, ano_flow_circuit_dict]

        # switch for cell discretizsation
        cell_dict['channel_land_discretization'] = \
            settings['simulation']['channel_land_discretization']

        # Initialize fluid channels
        fluids, channels = [], []
        for i in range(len(half_cell_dicts)):
            temp_in = manifold_in_dicts[i]['temp_in']
            p_out = manifold_out_dicts[i]['p_out']
            fluid_dicts[i]['temperature'] = temp_in
            fluid_dicts[i]['pressure'] = p_out
            channel_dicts[i]['temp_in'] = temp_in
            channel_dicts[i]['p_out'] = p_out
            fluid_dicts[i]['nodes'] = n_nodes
            # fluids.append(fluid.factory2(fluid_dicts[i]))
            channels.append([chl.Channel(channel_dicts[i],
                                         fluid.factory(fluid_dicts[i]))
                             for j in range(self.n_cells)])

        # Initialize fuel cells
        self.cells = []
        endplate_heat_flux = temperature_dict['heat_flux']
        for i in range(self.n_cells):
            if self.n_cells == 1:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = True
                cell_dict['heat_flux'] = endplate_heat_flux
            elif i == 0:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = False
                cell_dict['heat_flux'] = endplate_heat_flux
            elif i == self.n_cells-1:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = True
            else:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = False
                cell_dict['heat_flux'] = endplate_heat_flux

            cell_channels = [channels[0][i], channels[1][i]]
            # Cell constructor
            cell = cl.Cell(cell_dict, membrane_dict, half_cell_dicts,
                           cell_channels, number=i)
            if i == 0:
                cell.coords[0] = 0.0
                cell.coords[1] = cell.thickness
            else:
                cell.coords[0] = self.cells[i - 1].coords[1]
                cell.coords[1] = cell.coords[0] + cell.thickness
            self.cells.append(cell)

        # Initialize flow circuits
        manifold_length = \
            self.cells[-1].coords[-1] - self.cells[0].coords[0]
        self.fuel_circuits = []
        for i in range(len(half_cell_dicts)):
            manifold_in_dicts[i]['length'] = manifold_length
            manifold_out_dicts[i]['length'] = manifold_length
            sub_channel_number = self.cells[0].half_cells[i].n_channel
            self.fuel_circuits.append(
                flow_circuit.factory(flow_circuit_dicts[i],
                                     manifold_in_dicts[i],
                                     manifold_out_dicts[i],
                                     channels[i], sub_channel_number))

        cool_flow = stack_dict['cool_flow']
        if cool_flow:
            coolant_channel_dict = settings['coolant_channel']
            coolant_dict = coolant_channel_dict['fluid']
            coolant_dict['nodes'] = n_nodes
            dict_coolant_flow_circuit = \
                settings['coolant_flow_circuit']
            dict_coolant_in_manifold = \
                dict_coolant_flow_circuit['inlet_manifold']
            dict_coolant_out_manifold = \
                dict_coolant_flow_circuit['outlet_manifold']
            dict_coolant_in_manifold['length'] = manifold_length
            dict_coolant_out_manifold['length'] = manifold_length
            coolant_dict['temperature'] = dict_coolant_in_manifold['temp_in']
            coolant_dict['pressure'] = dict_coolant_out_manifold['p_out']

            cool_bc = temperature_dict['cool_ch_bc']
            if cool_bc:
                n_cool = self.n_cells + 1
            else:
                n_cool = self.n_cells - 1

            n_cool_cell = temperature_dict['cool_ch_numb']
            cool_channels = []
            for i in range(n_cool):
                cool_channels.append(
                    chl.Channel(coolant_channel_dict,
                                fluid.factory(coolant_dict),
                                number=str(i)))
                cool_channels[i].extend_data_names(cool_channels[i].name)
                cool_channels[i].fluid.name = \
                    cool_channels[i].name + ' Fluid'
                # + cool_channels[i].fluid.TYPE_NAME
                cool_channels[i].fluid.extend_data_names(
                    cool_channels[i].fluid.name)
            if cool_bc:
                cool_channels[0].height *= 1.0
                cool_channels[-1].height *= 1.0
            if n_cool > 0:
                self.coolant_circuit = \
                    flow_circuit.factory(dict_coolant_flow_circuit,
                                         dict_coolant_in_manifold,
                                         dict_coolant_out_manifold,
                                         cool_channels, n_cool_cell)
            else:
                self.coolant_circuit = None
        else:
            self.coolant_circuit = None
        self.coolant_temp_diff = \
            temperature_dict.get('cool_temp_diff', None)
        self.coolant_mass_flow = \
            temperature_dict.get('cool_mass_flow', None)

        if 'coolant_temp_control' in temperature_dict:
            self.coolant_temp_control = temperature_dict['cool_temp_control']
        elif 'coolant_control' in temperature_dict:
            if temperature_dict['coolant_control'].lower() \
                    == 'temperature difference':
                self.coolant_temp_control = True
            elif temperature_dict['coolant_control'].lower() == 'mass flow':
                self.coolant_temp_control = False
        elif self.coolant_temp_diff is not None \
                and self.coolant_mass_flow is None:
            self.coolant_temp_control = True
        elif self.coolant_temp_diff is None \
                and self.coolant_mass_flow is not None:
            self.coolant_temp_control = False
        else:
            raise ValueError('either mass flow or temperature difference '
                             'must be provided for coolant flow')

        self.flow_circuits = \
            [self.fuel_circuits[0], self.fuel_circuits[1], self.coolant_circuit]

        self.current_control = current_control

        i_cd_target = np.asarray(stack_dict['init_current_density'])
        if i_cd_target.ndim > 0:
            self.i_cd_target = i_cd_target[0]
        else:
            self.i_cd_target = i_cd_target
        v_target = np.asarray(stack_dict['init_current_density'])
        if v_target.ndim > 0:
            self.v_target = v_target[0]
        else:
            self.v_target = v_target

        # Initialize the electrical coupling
        self.elec_sys = el_cpl.ElectricalCoupling(self)
        # Initialize temperature src
        self.temp_sys = therm_cpl.TemperatureSystem(self, temperature_dict)

        """boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        # target current density

        # current density array
        self.i_cd = np.zeros((self.n_cells, n_ele))
        self.i_cd[:] = self.i_cd_target
        # for i in range(self.n_cells):
        #     self.i_cd[i, :] = \
        #         g_func.exponential_distribution(self.i_target, n_ele,
        #                                         a=0.5, b=0.0)
        # current density array of previous iteration step
        self.i_cd_old = np.copy(self.i_cd)
        self.i_cd_avg = self.i_cd_target
        # voltage array
        self.v = np.zeros(self.n_cells)
        self.v_stack = None
        self.v_loss = None
        self.e_0 = self.n_cells * self.cells[0].e_0

        # old temperature for convergence calculation
        self.temp_old = np.zeros(self.temp_sys.temp_layer_vec.shape)
        self.temp_old[:] = self.temp_sys.temp_layer_vec

    def update(self, current_density=None, voltage=None):
        """
        This function coordinates the program sequence
        """
        update_inflows = False
        if current_density is not None:
            self.i_cd[:] = current_density
            self.i_cd_avg = current_density
            update_inflows = True
        elif voltage is not None:
            self.v_stack = voltage
            update_inflows = True
        if self.current_control is False:
            update_inflows = True
        self.update_flows(update_inflows,
                          coolant_temp_diff=self.coolant_temp_diff,
                          coolant_mass_flow=self.coolant_mass_flow)
        for i, cell in enumerate(self.cells):
            cell.update(self.i_cd[i, :], current_control=self.current_control,
                        update_channel=False)
            if cell.break_program:
                self.break_program = True
                break
        self.i_cd_old[:] = self.elec_sys.i_cd
        self.temp_old[:] = self.temp_sys.temp_layer_vec
        if not self.break_program:
            if self.calc_temp:
                self.temp_sys.update()
            if self.calc_electric:
                self.elec_sys.update(current_density=current_density,
                                     voltage=voltage)
                self.i_cd[:] = self.elec_sys.i_cd
            self.v[:] = \
                np.asarray([np.average(cell.v, weights=cell.active_area_dx)
                            for cell in self.cells])
            if self.current_control:
                self.v_stack = np.sum(self.v)
                self.v_loss = self.e_0 - self.v_stack
            self.i_cd_avg = np.average(self.i_cd[0],
                                       weights=self.cells[0].active_area_dx)

    def update_flows(self, update_inflows=False,
                     coolant_temp_diff=None, coolant_mass_flow=None):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        mass_flows_in = [None, None]
        if update_inflows:
            mass_flows_in[:] = self.calc_mass_flows()
        for i in range(len(self.fuel_circuits)):
            self.fuel_circuits[i].update(mass_flows_in[i])
        if self.coolant_circuit is not None:
            cool_mass_flow = None
            if self.current_control or update_inflows:
                if not self.coolant_temp_control:
                    cool_mass_flow = coolant_mass_flow
                else:
                    cool_mass_flow = self.calc_cool_mass_flow(coolant_temp_diff)
            self.coolant_circuit.update(cool_mass_flow)

    def calc_cool_mass_flow(self, coolant_temp_diff):
        n_cool_cell = self.coolant_circuit.n_subchannels
        if self.v_loss is None:
            v_loss = 0.5 * self.n_cells
        else:
            v_loss = self.v_loss
        heat = self.i_cd_avg * self.cells[0].active_area * v_loss
        cp_cool = \
            np.average([np.average(channel.fluid.specific_heat)
                        for channel in self.coolant_circuit.channels])
        return heat / (cp_cool * coolant_temp_diff)  # * n_cool_cell

    def calc_mass_flows(self):
        mass_flows_in = []
        for i in range(len(self.cells[0].half_cells)):
            cell_mass_flow, cell_mole_flow = \
                self.cells[0].half_cells[i].calc_inlet_flow(self.i_cd_avg)
            cell_mass_flow = np.sum(cell_mass_flow, axis=0)
            mass_flow = cell_mass_flow \
                * self.cells[0].half_cells[i].n_channel * self.n_cells
            mass_flows_in.append(mass_flow)
        return mass_flows_in

