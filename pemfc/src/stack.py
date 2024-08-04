# General imports
import numpy as np
# Local module imports
from . import (flow_circuit as flow_circuit, cell as cl, channel as chl,
               linear_system as lin_sys)
from .fluid import fluid as fluid
from .output_object import OutputObject1D


class Stack(OutputObject1D):

    def __init__(self, settings, n_nodes, current_control=False):

        super().__init__('Stack')
        # Read settings dictionaries
        stack_dict = settings['stack']
        # Number of cells of the stack
        self.n_cells = stack_dict['cell_number']
        # Switch to calculate the temperature distribution
        self.calc_temp = stack_dict['calc_temperature']
        # Switch to calculate the current density distribution
        self.calc_electric = stack_dict['calc_current_density']
        # Switch to calculate the flow distribution
        # self.calc_flow_dis = stack_dict['calc_flow_distribution']

        # Decompose input dict the individual objects
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

        # Add evaporation model settings to channel dictionary
        keys = ('two_phase_flow', 'evaporation_model')
        for i in range(len(half_cell_dicts)):
            two_phase_settings = half_cell_dicts[i].get(keys[0], None)
            if two_phase_settings is not None:
                evaporation_settings = two_phase_settings.get(keys[1], None)
                if evaporation_settings is not None :
                    channel_dicts[i][keys[1]] = evaporation_settings

        # Switch for cell discretizsation
        cell_dict['channel_land_discretization'] = \
            settings['simulation']['channel_land_discretization']
        # Add underrelaxation factor to cell settings
        cell_dict['underrelaxation_factor'] = \
            settings['simulation']['underrelaxation_factor']
        # Add GDL diffusion and two-phase switches to half-cell dicts
        for half_cell_dict in half_cell_dicts:
            half_cell_dict['calc_gdl_diffusion'] = (
                settings['simulation'].get('calc_gdl_diffusion', False))
            half_cell_dict['calc_two_phase_flow'] = (
                settings['simulation'].get('calc_two_phase_flow', False))

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
                                         fluid.create(fluid_dicts[i]))
                             for j in range(self.n_cells)])

        # Initialize fuel cells
        self.cells = []
        if temperature_dict.get('bc_endplate') == 'fixed':
            try:
                cell_dict['value_endplate'] = temperature_dict['value_endplate']
            except KeyError:
                raise KeyError('value for "value_endplate" must be provided, '
                               'when using "fixed" end-plate boundary '
                               'conditions')
        else:
            cell_dict['flux_endplate'] = temperature_dict['flux_endplate']

        if 'alpha_amb' in temperature_dict:
            cell_dict['alpha_amb'] = temperature_dict['alpha_amb']
            try:
                cell_dict['temp_amb'] = temperature_dict['temp_amb']
            except KeyError:
                raise KeyError("parameter 'temp_amb' must be provided if "
                               "'allpha_amb' is given")
        for i in range(self.n_cells):
            if self.n_cells == 1:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = True
            elif i == 0:
                cell_dict['first_cell'] = True
                cell_dict['last_cell'] = False
            elif i == self.n_cells-1:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = True
            else:
                cell_dict['first_cell'] = False
                cell_dict['last_cell'] = False

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
            sub_channel_number = self.cells[0].half_cells[i].n_channels
            self.fuel_circuits.append(
                flow_circuit.create(flow_circuit_dicts[i],
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
                                fluid.create(coolant_dict),
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
                    flow_circuit.create(dict_coolant_flow_circuit,
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

        # # prepend names of output values with flow circuit name
        # for item in self.flow_circuits:
        #     item.extend_data_names(item.name)

        self.current_control = current_control

        current_density_target = np.asarray(stack_dict['init_current_density'])
        if current_density_target.ndim > 0:
            self.current_density_target = current_density_target[0]
        else:
            self.current_density_target = current_density_target
        voltage_target = np.asarray(stack_dict['init_current_density'])
        if voltage_target.ndim > 0:
            self.voltage_target = voltage_target[0]
        else:
            self.voltage_target = voltage_target

        # Initialize temperature system
        # self.temp_sys = therm_cpl.TemperatureSystem(self, temperature_dict)
        self.temp_sys = lin_sys.TemperatureSystem(self, temperature_dict)
        # Initialize the electrical coupling
        self.elec_sys = lin_sys.ElectricalSystem(self, {})

        """Boolean alarms"""
        self.v_alarm = False
        # True if :voltage loss > cell voltage
        self.break_program = False
        # True if the program aborts because of some critical impact

        # Stack current density array
        self.current_density = np.zeros(self.elec_sys.current_density.shape)
        self.current_density[:] = self.current_density_target
        # for i in range(self.n_cells):
        #     self.i_cd[i, :] = \
        #         g_func.exponential_distribution(self.i_target, n_ele,
        #                                         a=0.5, b=0.0)
        # Current density array of previous iteration step
        self.current_density_old = np.copy(self.current_density)
        self.current_density_avg = self.current_density_target
        # Voltage array
        self.voltage_cells = np.zeros(self.n_cells)
        self.voltage_stack = None
        self.voltage_loss = None
        self.e_0 = self.n_cells * self.cells[0].e_0

        # old temperature for convergence calculation
        self.temp_old = np.zeros(self.temp_sys.solution_vector.shape)
        self.temp_old[:] = self.temp_sys.solution_vector

        # Add data container for output
        self.add_print_data(self.voltage_cells, 'Cell Voltage', 'V')

    def update(self, current_density=None, voltage=None):
        """
        This function coordinates the program sequence
        """
        update_inflows = False
        if current_density is not None:
            self.current_density_avg = current_density
        if any((current_density, voltage)):
            update_inflows = True
        if self.current_control is False:
            update_inflows = True
        self.update_flows(update_inflows,
                          coolant_temp_diff=self.coolant_temp_diff,
                          coolant_mass_flow=self.coolant_mass_flow)
        for i, cell in enumerate(self.cells):
            cell.update(self.current_density[i, :],
                        current_control=self.current_control,
                        update_channel=True)
            if cell.break_program:
                self.break_program = True
                break
        self.temp_old[:] = self.temp_sys.solution_vector
        if not self.break_program:
            if self.calc_temp:
                self.temp_sys.update()
            if self.calc_electric:
                self.update_electric_system(current_density, voltage)

    def update_electric_system(self, current_density, voltage):
        self.current_density_old[:] = self.elec_sys.current_density
        self.elec_sys.update(current_density=current_density,
                             voltage=voltage)
        self.current_density[:] = self.elec_sys.current_density

        self.voltage_cells[:] = np.asarray(
            [np.average(cell.e_0 - cell.voltage_loss, weights=cell.d_area)
             for cell in self.cells])

        if self.current_control:
            self.voltage_stack = np.sum(self.voltage_cells)
            self.voltage_loss = self.e_0 - self.voltage_stack

        # Calculate average current density from first cell
        ref_cell = self.cells[-1]
        ref_current_density = (
            ref_cell.current_density)[ref_cell.layer_id['membrane']]
        self.current_density_avg = np.average(ref_current_density,
                                              weights=ref_cell.d_area)

    def update_flows(self, update_inflows=False,
                     coolant_temp_diff=None, coolant_mass_flow=None):
        """
        This function updates the flow distribution of gas over the stack cells
        """
        mass_flows_in = [None, None]
        if update_inflows:
            mass_flows_in[:] = self.calc_fuel_mass_flows()
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
        if self.voltage_loss is None:
            v_loss = 0.5 * self.n_cells
        else:
            v_loss = self.voltage_loss
        v_heat = self.temp_sys.e_tn * self.n_cells - self.e_0 + v_loss
        heat = self.current_density_avg * self.cells[0].active_area * v_heat
        cp_cool = \
            np.average([np.average(channel.fluid.specific_heat)
                        for channel in self.coolant_circuit.channels])
        return heat / (cp_cool * coolant_temp_diff)  # * n_cool_cell

    def calc_fuel_mass_flows(self):
        mass_flows_in = []
        for i in range(len(self.cells[0].half_cells)):
            cell_mass_flow, cell_mole_flow = \
                self.cells[0].half_cells[i].calc_inlet_flow(
                    self.current_density_avg)
            cell_mass_flow = np.sum(cell_mass_flow, axis=0)
            mass_flow = cell_mass_flow \
                * self.cells[0].half_cells[i].n_channels * self.n_cells
            mass_flows_in.append(mass_flow)
        return mass_flows_in

