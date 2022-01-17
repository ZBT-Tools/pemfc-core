
# Cooling Settings
cool_circuit = \
    {'label': 'Activate Cooling:', 'value': True, 'sticky': ['NW', 'NWE'],
     'sim_name': ['stack', 'cool_flow'], 'type': 'CheckButtonSet',
     'specifier': 'checklist_activate_cooling',
     'command': {'function': 'set_status',
                 'args': [[[1, 0], [1, 1], [1, 2],
                           [2, 0], [2, 1], [2, 2],
                           [3, 0], [3, 1], [3, 2],
                           [4, 0], [4, 1], [4, 2],
                           [5, 0], [5, 1], [5, 2],
                           [6, 0], [6, 1], [6, 2],
                           [7, 0], [7, 1], [7, 2]]],
                 'args2': [1, 2, 3, 4, 5, 6, 7]}}

cool_channel_length = \
    {'label': 'Coolant Channel Length:', 'value': 0.4,
     'sim_name': ['coolant_channel', 'length'], 'specifier': 'disabled_cooling',
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_channel_height = \
    {'label': 'Coolant Channel Height:', 'value': 1e-3,
     'sim_name': ['coolant_channel', 'height'], 'specifier': 'disabled_cooling',
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}
cool_channel_width = \
    {'label': 'Coolant Channel Width:', 'value': 1e-3,
     'sim_name': ['coolant_channel', 'width'], 'specifier': 'disabled_cooling',
     'dtype': 'float', 'dimensions': 'm', 'type': 'EntrySet'}

cool_channel_number = \
    {'label': 'Coolant Channel Number:', 'value': 2,
     'sim_name': ['temperature_system', 'cool_ch_numb'],
     'specifier': 'disabled_cooling', 'type': 'EntrySet'}

cool_channel_bends = \
    {'label': 'Number of Coolant Channel Bends:', 'value': 0,
     'sim_name': ['coolant_channel', 'bend_number'],
     'specifier': 'disabled_cooling', 'type': 'EntrySet'}

cool_bend_pressure_loss_coefficient = \
    {'label': 'Pressure Loss Coefficient for Coolant Channel Bend:',
     'sim_name': ['coolant_channel', 'bend_friction_factor'],
     'specifier': 'disabled_cooling', 'value': 0.5, 'dimensions': '-',
     'type': 'EntrySet'}

cool_flow_end_cells = \
    {'label': 'Activate Cooling Flow at End Plates:', 'value': False,
     'sim_name': ['temperature_system', 'cool_ch_bc'],
     'sticky': ['NW', 'NWE'], 'type': 'CheckButtonSet'}

channel_flow_direction = \
    {'label': 'Channel Flow Direction (1 or -1):', 'value': 1,
     'sim_name': ['coolant_channel', 'flow_direction'],
     'specifier': 'disabled_cooling', 'dtype': 'int', 'type': 'EntrySet',
     'sticky': ['NW', 'NWE']}

cool_frame_dict = \
    {'title': 'Cooling Settings', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN', 'size_label': 'xl', 'size_unit': 's',
     'widget_dicts': [cool_circuit, cool_channel_number,
                      cool_channel_length, cool_channel_height,
                      cool_channel_width, cool_channel_bends,
                      cool_bend_pressure_loss_coefficient,
                      channel_flow_direction,
                      cool_flow_end_cells],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Cooling', 'show_title': False,
            'sub_frame_dicts': [cool_frame_dict]}

# geometry_frame_dict = \
#     {'title': 'Geometry', 'show_title': False, 'font': 'Arial 10 bold',
#      'sub_frame_dicts': [cool_frame_dict, manifold_frame_dict,
#                          cell_frame_dict],
#      'highlightbackground': 'grey', 'highlightthickness': 1}

# main_frame_dicts = [geometry_frame_dict, simulation_frame_dict]







