channel_discretization = \
    {'label': 'Channel Discretization:', 'value': 10,
     'sim_name': ['simulation', 'elements'], 'dtype': 'int', 'type': 'EntrySet'}

convergence_criteria = \
    {'label': 'Convergence Tolerance:', 'value': 1e-6,
     'sim_name': ['simulation', 'convergence_criteria'], 'dtype': 'float',
     'dimensions': '-', 'type': 'EntrySet'}

minimum_iteration = \
    {'label': 'Minimum Number of Iterations:', 'value': 3,
     'sim_name': ['simulation', 'minimum_iteration'], 'dtype': 'int',
     'type': 'EntrySet'}

maximum_iteration = \
    {'label': 'Maximum Number of Iterations:', 'value': 50,
     'sim_name': ['simulation', 'maximum_iteration'], 'dtype': 'int',
     'type': 'EntrySet'}

underrelaxation_factor = \
    {'label': 'Under-relaxation Factor:', 'value': 0.5,
     'sim_name': ['cell', 'underrelaxation_factor'], 'dtype': 'float',
     'dimensions': '-', 'type': 'EntrySet'}

main_numerical_settings_label = \
    {'label': 'Main Numerical Settings',
     'font': 'Arial 10 bold', 'type': 'Label', 'sticky': 'WNS'}

main_numerical_settings_frame_dict = \
    {'title': 'Main Numerical Settings', 'show_title': False,
     'widget_dicts': [main_numerical_settings_label,
                      channel_discretization,
                      convergence_criteria,
                      minimum_iteration,
                      maximum_iteration,
                      underrelaxation_factor],
     'size_label': 'xl', 'size_unit': 's', 'font': 'Arial 10', 'sticky': 'WEN'}
# , 'highlightbackground': 'grey', 'highlightthickness': 1}

anode_label = {'label': 'Anode', 'row': 1, 'column': 1,
                        'type': 'Label', 'sticky': 'WENS'}

cathode_label = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}

cooling_label = \
    {'label': 'Cooling', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}

convergence_criteria_flow = \
    {'label': 'Convergence Tolerance:',
     'value': [1e-6, 1e-6, 1e-6],
     'sim_name': [['anode', 'flow_circuit', 'tolerance'],
                  ['cathode', 'flow_circuit', 'tolerance'],
                  ['coolant_flow_circuit', 'tolerance']],
     'sticky': ['NW', 'NWE'], 'dimensions': '-',
     'dtype': 'float', 'type': 'EntrySet'}

min_iter_flow = \
    {'label': 'Minimum Iteration Number:',
     'value': [3, 3, 3],
     'sim_name': [['anode', 'flow_circuit', 'min_iter'],
                  ['cathode', 'flow_circuit', 'min_iter'],
                  ['coolant_flow_circuit', 'min_iter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'int', 'type': 'EntrySet'}

max_iter_flow = \
    {'label': 'Maximum Iteration Number:',
     'value': [20, 20, 20],
     'sim_name': [['anode', 'flow_circuit', 'max_iter'],
                  ['cathode', 'flow_circuit', 'max_iter'],
                  ['coolant_flow_circuit', 'max_iter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'int', 'type': 'EntrySet'}

underrelaxation_factor_flow = \
    {'label': 'Under-relaxation Factor:',
     'value': [0.5, 0.5, 0.5],
     'sim_name': [['anode', 'flow_circuit', 'underrelaxation_factor'],
                  ['cathode', 'flow_circuit', 'underrelaxation_factor'],
                  ['coolant_flow_circuit', 'underrelaxation_factor']],
     'sticky': ['NW', 'NWE'], 'dimensions': '-',
     'dtype': 'float', 'type': 'EntrySet'}

flow_numerical_settings_frame_dict = \
    {'title': 'Flow Circuit Algorithm', 'show_title': True,
     'font': 'Arial 10', 'size_label': 'm', 'size_unit': 's',
     'widget_dicts': [anode_label, cathode_label, cooling_label,
                      convergence_criteria_flow,
                      min_iter_flow,
                      max_iter_flow,
                      underrelaxation_factor_flow],
     'sticky': 'WEN'}
# , 'highlightbackground': 'grey', 'highlightthickness': 1}


# anode_label_2 = {'label': 'Anode', 'row': 1, 'column': 1,
#                         'type': 'Label', 'sticky': 'WENS'}
#
# cathode_label_2 = \
#     {'label': 'Cathode', 'row': 1, 'column': 2,
#      'type': 'Label', 'sticky': 'WENS'}

c_eps = \
    {'label': 'Critical Concentration Value:', 'value': [0.02, 0.02],
     'sim_name': [['anode', 'c_eps'], ['cathode', 'c_eps']],
     'dtype': 'float', 'dimensions': 'mol/m³', 'type': 'EntrySet'}

delta_i = \
    {'label': 'Numerical Current Density Difference:', 'value': [5.0, 5.0],
     'sim_name': [['anode', 'delta_i'], ['cathode', 'delta_i']],
     'dtype': 'float', 'dimensions': 'A/m²', 'type': 'EntrySet'}

current_linearization_frame_dict = \
    {'title': 'Limiting Current Linearization', 'show_title': True,
     'font': 'Arial 10', 'size_label': 'l', 'size_unit': 'm',
     'widget_dicts': [anode_label, cathode_label,
                      c_eps,
                      delta_i],
     'sticky': 'WEN'}
# , 'highlightbackground': 'grey', 'highlightthickness': 1}

numerical_settings_frame_dict = \
    {'title': 'Numerical Settings', 'show_title': True,
     'font': 'Arial 10 bold',
     'sub_frame_dicts': [main_numerical_settings_frame_dict,
                         flow_numerical_settings_frame_dict,
                         current_linearization_frame_dict],
     'sticky': 'WEN', 'highlightbackground': 'grey', 'highlightthickness': 1}

calc_temperature = \
    {'label': 'Calculate Temperature Distribution:', 'value': True,
     'sim_name': ['stack', 'calc_temperature'],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

calc_current_density = \
    {'label': 'Calculate Current Density Distribution:', 'value': True,
     'sim_name': ['stack', 'calc_current_density'],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

calc_mem_loss = \
    {'label': 'Calculate Membrane Loss:', 'value': True,
     'sim_name': ['membrane', 'calc_loss'],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

main_physics_switch_frame_dict = \
    {'title': 'Main Physics Model Switches', 'show_title': False,
     'font': 'Arial 10 bold', 'size_label': 'xl', 'size_unit': 'null',
     'widget_dicts': [calc_temperature,
                      calc_current_density,
                      calc_mem_loss],
     'sticky': 'WEN'}
# , 'highlightbackground': 'grey', 'highlightthickness': 1}


calc_act_loss = \
    {'label': 'Calculate Activation Loss:', 'value': [True, True],
     'sim_name': [['anode', 'calc_act_loss'],
                  ['cathode', 'calc_act_loss']],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

calc_cl_diff_loss = \
    {'label': 'Calculate Catalyst Diffusion Loss:', 'value': [True, True],
     'sim_name': [['anode', 'calc_cl_diff_loss'],
                  ['cathode', 'calc_cl_diff_loss']],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

calc_gdl_diff_loss = \
    {'label': 'Calculate Catalyst Diffusion Loss:', 'value': [True, True],
     'sim_name': [['anode', 'calc_gdl_diff_loss'],
                  ['cathode', 'calc_gdl_diff_loss']],
     'dtype': 'bool', 'type': 'CheckButtonSet'}

electrode_loss_switch_frame_dict = \
    {'title': 'Physics Model Switches', 'show_title': False,
     'font': 'Arial 10 bold', 'size_label': 'l', 'size_unit': 'null',
     'widget_dicts': [anode_label, cathode_label,
                      calc_act_loss,
                      calc_cl_diff_loss,
                      calc_gdl_diff_loss],
     'sticky': 'WEN'}
     # , 'highlightbackground': 'grey', 'highlightthickness': 1}

physics_switch_frame_dict = \
    {'title': 'Physics Model Switches', 'show_title': True,
     'font': 'Arial 10 bold',
     'sub_frame_dicts': [main_physics_switch_frame_dict,
                         electrode_loss_switch_frame_dict],
     'sticky': 'WEN', 'highlightbackground': 'grey', 'highlightthickness': 1}

load_settings_button_dict = \
    {'label': 'Load Settings', 'takefocus': 0, 'row': 0, 'column': 0,
     'sticky': '', 'width': 20, 'weights': [1, 1],
     'type': 'OpenFileButton'}
save_settings_button_dict = \
    {'label': 'Save Settings', 'takefocus': 0, 'row': 0, 'column': 1,
     'sticky': '', 'width': 20, 'weights': [1, 1],
     'type': 'SaveFileButton'}

load_settings_frame_dict = \
    {'title': 'Settings IO', 'show_title': False, 'font': 'Arial 10 bold',
     'widget_dicts': [load_settings_button_dict,
                      save_settings_button_dict],
     'sticky': 'WEN'}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}

output_dir_button_dict = \
    {'label': 'Open', 'type': 'OpenDirectoryButton', 'width': 10}
output_dir = \
    {'label': 'Output Directory:', 'button_dict': output_dir_button_dict,
     'sim_name': ['output', 'directory'], 'width': 40,
     'dtype': 'string', 'type': 'EntryButtonSet'}


output_frame_dict = \
    {'title': 'Output Results', 'show_title': False, 'font': 'Arial 10 bold',
     'widget_dicts': [output_dir],
     'sticky': 'WENS'}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}

run_button_dict = {'label': 'Run Simulation', 'type': 'RunButton',
                   'columnspan': 3, 'width': 20, 'sticky': 'WNE'}

empty_row = {'label': ' ',  'font': 'Arial 1',  'height': 100,
             'type': 'Label', 'sticky': 'WENS'}

run_button_frame_dict = \
    {'title': 'Run Simulation', 'show_title': False,
     'widget_dicts': [run_button_dict], 'sticky': ''}

simulation_frame_dict = \
    {'title': 'Simulation Settings', 'show_title': False,
     'sub_frame_dicts': [numerical_settings_frame_dict,
                         physics_switch_frame_dict,
                         output_frame_dict,
                         load_settings_frame_dict,
                         run_button_frame_dict],
     # 'widget_dicts': [cell_number, cell_length, cell_width],
     'sticky': 'WEN'}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}


tab_dict = \
    {'title': 'Simulation', 'show_title': False,
     'sub_frame_dicts': [simulation_frame_dict]}
