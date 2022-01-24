# general imports
import os
import numpy as np
import cProfile
import timeit

# local module imports
from . import stack
from . import output
from data import input_dicts
# from ..gui import data_transfer


def do_c_profile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats('cumtime')
    return profiled_func


class Simulation:

    def __init__(self, ):
        # dict_simulation = input_dicts.dict_simulation
        dict_simulation = input_dicts.sim_dict['simulation']
        # set convergence criteria for current and temperature distribution
        self.it_crit = dict_simulation['convergence_criteria']

        # maximum and minimum number of iterations for global algorithm
        self.max_it = dict_simulation['maximum_iteration']
        self.min_it = dict_simulation['minimum_iteration']

        self.timing = {'start': 0.0,
                       'initialization': 0.0,
                       'simulation': 0.0,
                       'output': 0.0}

        """General variables"""
        # initialize stack object
        n_nodes = dict_simulation['elements'] + 1
        self.current_control = dict_simulation.get('current_control', True)
        if 'operation_control' in dict_simulation:
            if dict_simulation['operation_control'].lower() == 'current':
                self.current_control = True
            elif dict_simulation['operation_control'].lower() == 'voltage':
                self.current_control = False

        self.current_density = dict_simulation['current_density']
        self.average_cell_voltage = dict_simulation['average_cell_voltage']
        if self.current_control and self.current_density is None:
            raise ValueError('parameter current_density must be provided')
        elif not self.current_density and self.average_cell_voltage is None:
            raise ValueError('parameter average_cell_voltage must be provided')
        stack_dict = input_dicts.sim_dict['stack']
        cell_number = stack_dict['cell_number']
        if self.current_control:
            stack_dict['init_current_density'] = self.current_density
        else:
            stack_dict['voltage'] = self.average_cell_voltage * cell_number

        self.stack = stack.Stack(n_nodes, current_control=self.current_control)

        # initialize output object
        output_dict = input_dicts.sim_dict['output']
        self.output = output.Output(output_dict)

    # @do_c_profile
    def run(self):
        """
        This function coordinates the program sequence
        """
        if self.current_control:
            target_value = self.current_density
        else:
            target_value = self.average_cell_voltage * self.stack.n_cells
        if not isinstance(target_value, (list, tuple, np.ndarray)):
            target_value = [target_value]
        cell_voltages = []
        current_densities = []
        local_data_list = []
        global_data_list = []
        for i, tar_value in enumerate(target_value):
            current_errors = []
            temp_errors = []
            self.timing['simulation'] = 0.0
            self.timing['output'] = 0.0
            simulation_start_time = timeit.default_timer()
            counter = 0
            while True:
                if counter == 0:
                    if self.current_control:
                        self.stack.update(current_density=tar_value)
                    else:
                        self.stack.update(voltage=tar_value)
                else:
                    self.stack.update()
                if self.stack.break_program:
                    break
                current_error, temp_error = self.calc_convergence_criteria()
                current_errors.append(current_error)
                temp_errors.append(temp_error)
                if len(target_value) < 1:
                    print(counter)
                counter += 1
                if ((current_error < self.it_crit and temp_error < self.it_crit)
                        and counter > self.min_it) or counter > self.max_it:
                    break
            simulation_stop_time = timeit.default_timer()
            simulation_time = simulation_stop_time - simulation_start_time
            self.timing['simulation'] += simulation_time
            output_start_time = timeit.default_timer()

            if not self.stack.break_program:
                # voltage_loss = self.get_voltage_losses(self.stack)
                cell_voltages.append(np.average([cell.v for cell in
                                                 self.stack.cells]))
                current_densities.append(self.stack.i_cd_avg)

                case_name = 'Case'+str(i)
                self.output.save(case_name, self.stack)
                local_data_list.append(self.output.get_data(self.stack))
                if self.output.save_plot:
                    path = os.path.join(self.output.output_dir, case_name,
                                        'plots', 'Convergence.png')
                    self.output.create_figure(path, list(range(counter)),
                                              [current_errors, temp_errors],
                                              xlabels='Iteration',
                                              ylabels='Error',
                                              yscale='log',
                                              legend=['Current Density',
                                                      'Temperature'])
            else:
                target_value = target_value[0:-i]
                break

            average_current_density = \
                np.average([np.average(cell.i_cd, weights=cell.active_area_dx)
                            for cell in self.stack.cells])
            global_data = \
                {'Stack Voltage': {'value': self.stack.v_stack, 'units': 'V'},
                 'Average Cell Voltage':
                     {'value': self.stack.v_stack / self.stack.n_cells,
                      'units': 'V'},
                 'Minimum Cell Voltage':
                     {'value': np.min(self.stack.v), 'units': 'V'},
                 'Maximum Cell Voltage':
                     {'value': np.max(self.stack.v), 'units': 'V'},
                 'Average Current Density':
                     {'value': average_current_density, 'units': 'A/m²'},
                 'Stack Power Density':
                     {'value': self.stack.v_stack * average_current_density,
                      'units': 'W/m²'},
                 'Stack Power':
                     {'value': self.stack.v_stack * average_current_density
                      * self.stack.cells[0].active_area,
                      'units': 'W'},
                 'Cooling Mass Flow Rate:':
                     {'value': self.stack.coolant_circuit.mass_flow_in,
                      'units': 'kg/s', 'format': '.4E'},
                 'Cathode Mass Flow Rate:':
                     {'value': self.stack.fuel_circuits[0].mass_flow_in,
                      'units': 'kg/s', 'format': '.4E'},
                 'Anode Mass Flow Rate:':
                     {'value': self.stack.fuel_circuits[1].mass_flow_in,
                      'units': 'kg/s', 'format': '.4E'}
                 }
            global_data_list.append(global_data)
            output_stop_time = timeit.default_timer()
            self.timing['output'] += output_stop_time - output_start_time
            self.output.print_global_data(self, global_data)

        output_start_time = timeit.default_timer()

        # print polarization curve, if more than one current_density value
        # was provided
        if len(target_value) > 1:
            self.output.write_data(current_densities, cell_voltages,
                                   'Current Density [A/m²]', 'Cell Voltage',
                                   units='V', directory=self.output.output_dir,
                                   save_csv=True, save_plot=True,
                                   write_mode='w')

        output_stop_time = timeit.default_timer()
        self.timing['output'] += output_stop_time - output_start_time
        return global_data_list, local_data_list

    @staticmethod
    def get_voltage_losses(fc_stack):
        """
        Saves the average voltage losses of the stack
        """
        voltage_loss = \
            {'activation': {'anode': {}, 'cathode': {}},
             'diffusion': {'CL': {'anode': {}, 'cathode': {}},
                           'GDL': {'anode': {}, 'cathode': {}}},
             'membrane': {}}

        voltage_loss['activation']['anode']['cells'] = \
            np.asarray([cell.cathode.v_loss_act for cell in fc_stack.cells])
        voltage_loss['activation']['cathode']['cells'] = \
            np.asarray([cell.anode.v_loss_act for cell in fc_stack.cells])
        voltage_loss['diffusion']['CL']['anode']['cells'] = \
            np.asarray([cell.anode.v_loss_cl_diff for cell in fc_stack.cells])
        voltage_loss['diffusion']['CL']['cathode']['cells'] = \
            np.asarray([cell.cathode.v_loss_cl_diff for cell in fc_stack.cells])
        voltage_loss['diffusion']['GDL']['anode']['cells'] = \
            np.asarray([cell.anode.v_loss_gdl_diff for cell in fc_stack.cells])
        voltage_loss['diffusion']['GDL']['cathode']['cells'] = \
            np.asarray([cell.cathode.v_loss_gdl_diff
                        for cell in fc_stack.cells])
        voltage_loss['membrane']['cells'] = \
            np.asarray([cell.membrane.v_loss for cell in fc_stack.cells])

        voltage_loss['activation']['anode']['average'] = \
            np.average(voltage_loss['activation']['anode']['cells'])
        voltage_loss['activation']['cathode']['average'] = \
            np.average(voltage_loss['activation']['cathode']['cells'])
        voltage_loss['diffusion']['CL']['anode']['average'] = \
            np.average(voltage_loss['diffusion']['CL']['anode']['cells'])
        voltage_loss['diffusion']['CL']['cathode']['average'] = \
            np.average(voltage_loss['diffusion']['CL']['cathode']['cells'])
        voltage_loss['diffusion']['GDL']['anode']['average'] = \
            np.average(voltage_loss['diffusion']['GDL']['anode']['cells'])
        voltage_loss['diffusion']['GDL']['cathode']['average'] = \
            np.average(voltage_loss['diffusion']['GDL']['cathode']['cells'])
        voltage_loss['membrane']['average'] = \
            np.average(voltage_loss['membrane']['cells'])

        return voltage_loss

    def calc_convergence_criteria(self):
        """
        Calculates the convergence criteria according to (Koh, 2003)
        """
        i_cd_vec = self.stack.i_cd.flatten()
        current_error = \
            np.abs(np.sum(((i_cd_vec - self.stack.i_cd_old.flatten())
                           / i_cd_vec) ** 2.0))
        # self.temp_criteria =\
        #     np.abs(np.sum(((self.temp_old
        #                     - self.stack.temp_sys.temp_layer[0][0, 0]))
        #                   / self.stack.temp_sys.temp_layer[0][0, 0]))
        if self.stack.calc_temp:
            temp_error =\
                np.abs(np.sum(((self.stack.temp_old
                                - self.stack.temp_sys.temp_layer_vec)
                              / self.stack.temp_sys.temp_layer_vec) ** 2.0))
        else:
            temp_error = 0.0
        return current_error, temp_error
