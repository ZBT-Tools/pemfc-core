# Manifold Settings
calc_distribution = \
    {'label': 'Activate Calculation:', 'number': 3,
     'value': [False, False, False],
     'specifier': 'checklist_activate_calculation',
     'sim_name': [['anode', 'flow_circuit', 'calc_distribution'],
                  ['cathode', 'flow_circuit', 'calc_distribution'],
                  ['coolant_flow_circuit', 'calc_distribution']],
     'sticky': ['NW', 'NWE'], 'type': 'CheckButtonSet',
     'command': {'function': 'set_status',
                 'args': [[[3, 1], [6, 1], [7, 1], [8, 1], [9, 1], [10, 1],
                           [13, 1], [14, 1], [15, 1], [16, 1], [17, 1]],
                          [[3, 2], [6, 2], [7, 2], [8, 2], [9, 2], [10, 2],
                           [13, 2], [14, 2], [15, 2], [16, 2], [17, 2]],
                          [[3, 3], [6, 3], [7, 3], [8, 3], [9, 3], [10, 3],
                           [13, 3], [14, 3], [15, 3], [16, 3], [17, 3]]],
                 'args2': [4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18]}}
                 # 'args': [[[3, 1], [4, 1], [5, 1], [6, 1], [7, 1], [8, 1],
                 #           [9, 1], [10, 1], [11, 1], [12, 1], [13, 1]],
                 #          [[3, 2], [4, 2], [5, 2], [6, 2], [7, 2], [8, 2],
                 #           [9, 2], [10, 2], [11, 2], [12, 2], [13, 2]],
                 #          [[3, 3], [4, 3], [5, 3], [6, 3], [7, 3], [8, 3],
                 #           [9, 3], [10, 3], [11, 3], [12, 3], [13, 3]]]}}


anode_label_manifold = {'label': 'Anode', 'row': 1, 'column': 1,
                        'type': 'Label', 'sticky': 'WENS'}
cathode_label_manifold = \
    {'label': 'Cathode', 'row': 1, 'column': 2,
     'type': 'Label', 'sticky': 'WENS'}
cooling_label_manifold = \
    {'label': 'Cooling', 'row': 1, 'column': 3,
     'type': 'Label', 'sticky': 'WENS'}
empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
             'type': 'Label', 'sticky': 'WENS'}
inlet_manifold_label = {'label': 'Inlet Manifold', 'font': 'Arial 10 bold',
                        'type': 'Label', 'sticky': 'WNS'}
outlet_manifold_label = {'label': 'Outlet Manifold', 'font': 'Arial 10 bold',
                         'type': 'Label', 'sticky': 'WNS'}
manifold_configuration = \
    {'label': 'Flow Configuration:', 'number': 3,
     'specifier': 'disabled_manifolds',
     'sim_name': [['anode', 'flow_circuit', 'shape'],
                  ['cathode', 'flow_circuit', 'shape'],
                  ['coolant_flow_circuit', 'shape']],
     'value': ['U', 'Z'], 'type': 'ComboboxSet'}


# inlet_manifold_cross_section = \
#     {'label': 'Inlet Manifold Cross-Section:', 'number': 3,
#      'sim_name': [['anode', 'flow_circuit', 'inlet_manifold',
#                    'cross_sectional_shape'],
#                   ['cathode', 'flow_circuit', 'inlet_manifold',
#                    'cross_sectional_shape'],
#                   ['coolant_flow_circuit', 'inlet_manifold',
#                    'cross_sectional_shape']],
#      'value': ['circular', 'rectangular'], 'type': 'ComboboxSet',
#      'command': {'function': 'show_connected_widgets',
#                  'args': [[[[[6, 1]], [[8, 1], [10, 1]]],
#                            [[[8, 1], [10, 1]], [[6, 1]]]],
#                           [[[[6, 2], [7, 2], [8, 2]], [[9, 2]]],
#                            [[[6, 2], [7, 2], [8, 2], [9, 2]], []]],
#                           [[[[6, 3], [7, 3], [8, 3]], [[9, 3]]],
#                            [[[6, 3], [7, 3], [8, 3], [9, 3]], []]]]
#                 }
#      }

inlet_manifold_cross_section = \
    {'label': 'Cross-Section:', 'number': 3, 'specifier': 'disabled_manifolds',
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape'],
                  ['coolant_flow_circuit', 'inlet_manifold',
                   'cross_sectional_shape']],
     'value': ['circular', 'rectangular'], 'type': 'ComboboxSet',
     'command': {'function': 'set_status',
                 # 'args': [[[[[5, 1]], [[6, 1], [7, 1]]],
                 #           [[[6, 1], [7, 1]], [[5, 1]]]],
                 #          [[[[5, 2]], [[6, 2], [7, 2]]],
                 #           [[[6, 2], [7, 2]], [[5, 2]]]],
                 #          [[[[5, 3]], [[6, 3], [7, 3]]],
                 #           [[[6, 3], [7, 3]], [[5, 3]]]]]
                 'args': [[[[[7, 1]], [[8, 1], [9, 1]]],
                           [[[8, 1], [9, 1]], [[7, 1]]]],
                          [[[[7, 2]], [[8, 2], [9, 2]]],
                           [[[8, 2], [9, 2]], [[7, 2]]]],
                          [[[[7, 3]], [[8, 3], [9, 3]]],
                           [[[8, 3], [9, 3]], [[7, 3]]]]]
                }
     }

outlet_manifold_cross_section = \
    {'label': 'Cross-Section:', 'number': 3, 'specifier': 'disabled_manifolds',
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape'],
                  ['cathode', 'flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape'],
                  ['coolant_flow_circuit', 'outlet_manifold',
                   'cross_sectional_shape']],
     'value': ['circular', 'rectangular'], 'type': 'ComboboxSet',
     'command': {'function': 'set_status',
                 # 'args': [[[[[10, 1]], [[11, 1], [12, 1]]],
                 #           [[[11, 1], [12, 1]], [[10, 1]]]],
                 #          [[[[10, 2]], [[11, 2], [12, 2]]],
                 #           [[[11, 2], [12, 2]], [[10, 2]]]],
                 #          [[[[10, 3]], [[11, 3], [12, 3]]],
                 #           [[[11, 3], [12, 3]], [[10, 3]]]]]
                 'args': [[[[[14, 1]], [[15, 1], [16, 1]]],
                           [[[15, 1], [16, 1]], [[14, 1]]]],
                          [[[[14, 2]], [[15, 2], [16, 2]]],
                           [[[15, 2], [16, 2]], [[14, 2]]]],
                          [[[[14, 3]], [[15, 3], [16, 3]]],
                           [[[15, 3], [16, 3]], [[14, 3]]]]]
                 }
     }

inlet_manifold_diameter = \
    {'label': 'Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'diameter'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'diameter'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'diameter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

outlet_manifold_diameter = \
    {'label': 'Diameter:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'diameter'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'diameter'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'diameter']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

inlet_manifold_height = \
    {'label': 'Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'height'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'height'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'height']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

outlet_manifold_height = \
    {'label': 'Height:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'height'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'height'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'height']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

inlet_manifold_width = \
    {'label': 'Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold', 'width'],
                  ['cathode', 'flow_circuit', 'inlet_manifold', 'width'],
                  ['coolant_flow_circuit', 'inlet_manifold', 'width']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

outlet_manifold_width = \
    {'label': 'Width:', 'value': [1e-2, 1e-2, 1e-2],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold', 'width'],
                  ['cathode', 'flow_circuit', 'outlet_manifold', 'width'],
                  ['coolant_flow_circuit', 'outlet_manifold', 'width']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': 'm',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

inlet_pressure_loss_coeff = \
    {'label': 'T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow_circuit', 'inlet_manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow_circuit', 'inlet_manifold',
                   'constant_friction_factor'],
                  ['coolant_flow_circuit', 'inlet_manifold',
                   'constant_friction_factor']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': '-',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

outlet_pressure_loss_coeff = \
    {'label': 'T-Junction Loss Coefficient:',
     'value': [0.4, 0.4, 0.4],
     'sim_name': [['anode', 'flow_circuit', 'outlet_manifold',
                   'constant_friction_factor'],
                  ['cathode', 'flow_circuit', 'outlet_manifold',
                   'constant_friction_factor'],
                  ['coolant_flow_circuit', 'outlet_manifold',
                   'constant_friction_factor']],
     'sticky': ['NW', 'NWE'], 'dtype': 'float', 'dimensions': '-',
     'specifier': 'disabled_manifolds', 'type': 'EntrySet'}

manifold_frame_dict = \
    {'title': 'Manifold Settings', 'show_title': False, 'font': 'Arial 10 bold',
     'sticky': 'WEN',
     'widget_dicts': [anode_label_manifold, cathode_label_manifold,
                      cooling_label_manifold, calc_distribution,
                      manifold_configuration,
                      empty_row,
                      inlet_manifold_label,
                      inlet_manifold_cross_section, inlet_manifold_diameter,
                      inlet_manifold_width, inlet_manifold_height,
                      inlet_pressure_loss_coeff,
                      empty_row,
                      outlet_manifold_label,
                      outlet_manifold_cross_section, outlet_manifold_diameter,
                      outlet_manifold_width, outlet_manifold_height,
                      outlet_pressure_loss_coeff]
     }
     # 'highlightbackground': 'grey', 'highlightthickness': 1}
command_order = list(range(len(manifold_frame_dict['widget_dicts'])))
command_order.insert(0, command_order.pop(7))
command_order.insert(1, command_order.pop(14))
manifold_frame_dict['command_order'] = command_order

tab_dict = {'title': 'Manifolds', 'show_title': False,
            'sub_frame_dicts': [manifold_frame_dict]}
