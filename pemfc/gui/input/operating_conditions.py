# Electrochemistry Settings
# current_control = \
#     {'label': 'Current Control:', 'value': True, 'sticky': ['NW', 'NWE'],
#      'sim_name': ['simulation', 'current_control'], 'type': 'CheckButtonSet',
#      # 'command': {'function': 'set_visibility',
#      #             'args': [[[1, 0], [1, 1], [1, 2]],
#      #                      [[2, 0], [2, 1], [2, 2]]]}
#      }

current_control = \
    {'label': 'Operation Control:', 'number': 1,
     'sim_name': ['simulation', 'operation_control'],
     'value': ['Current', 'Voltage'],
     'type': 'ComboboxSet', 'specifier': 'dropdown_activate',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[1, 0], [1, 1], [1, 2]],
                            [[2, 0], [2, 1], [2, 2]]],
                           [[[2, 0], [2, 1], [2, 2]],
                            [[1, 0], [1, 1], [1, 2]]]]],
                 }
     }

current_density = \
    {'label': 'Current Density:', 'value': 1000.0,
     'sim_name': ['simulation', 'current_density'], 'specifier': 'visibility',
     'dtype': 'float', 'dimensions': 'A/m²', 'type': 'EntrySet'}

average_cell_voltage = \
    {'label': 'Cell Voltage:', 'value': 0.5,
     'sim_name': ['simulation', 'average_cell_voltage'],
     'specifier': 'visibility', 'dtype': 'float', 'dimensions': 'V',
     'type': 'EntrySet'}

electrochem_frame_dict = \
    {'title': 'Electrochemistry Settings', 'show_title': False, 'sticky': 'NWE',
     'widget_dicts': [current_control,
                      current_density,
                      average_cell_voltage],
     'size_label': 'l', 'size_unit': 'm',
     'highlightbackground': 'grey', 'highlightthickness': 1}

anode_label = {'label': 'Anode', 'row': 1, 'column': 1, 'columnspan': 2,
               'pady': 0, 'type': 'Label', 'sticky': 'WENS'}
cathode_label = {'label': 'Cathode', 'row': 1, 'column': 3, 'columnspan': 2,
                 'pady': 0, 'type': 'Label', 'sticky': 'WENS'}
hydrogen_label = {'label': 'H2', 'row': 2, 'column': 1, 'pady': 0,
                  'type': 'Label', 'sticky': 'WENS'}
nitrogen_label_1 = {'label': 'N2', 'row': 2, 'column': 2, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
oxygen_label = {'label': 'O2', 'row': 2, 'column': 3, 'pady': 0,
                'type': 'Label', 'sticky': 'WENS'}
nitrogen_label_2 = {'label': 'N2', 'row': 2, 'column': 4, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
dry_gas_composition_dict = \
    {'label': 'Dry Molar Gas Fractions:',
     'number': 4, 'value': [1.0, 0.0, 0.21, 0.79], 'width': 5,
     'sim_name':
         [['anode', 'channel', 'fluid', 'components', 'H2', 'molar_fraction'],
          ['anode', 'channel', 'fluid', 'components', 'N2', 'molar_fraction'],
          ['cathode', 'channel', 'fluid', 'components', 'O2',
           'molar_fraction'],
          ['cathode', 'channel', 'fluid', 'components', 'N2',
           'molar_fraction']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

relative_humidity = \
    {'label': 'Relative Humidity:', 'number': 2, 'value': 0.5,
     'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name':
         [['anode', 'channel', 'fluid', 'humidity'],
          ['cathode', 'channel', 'fluid', 'humidity']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

gas_inlet_temperature = \
    {'label': 'Gas Inlet Temperature:', 'number': 2, 'value': 343.15,
     'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name':
         [['anode', 'channel', 'temp_in'],
          ['cathode', 'channel', 'temp_in']],
     'dtype': 'float', 'dimensions': 'K', 'type': 'EntrySet'}

gas_outlet_pressure = \
    {'label': 'Gas Outlet Pressure:', 'number': 2, 'value': 101325.0,
     'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name':
         [['anode', 'channel', 'p_out'],
          ['cathode', 'channel', 'p_out']],
     'dtype': 'float', 'dimensions': 'Pa', 'type': 'EntrySet'}

gas_frame_dict = \
    {'title': 'Gas Supply Settings', 'show_title': False, 'sticky': 'NWE',
     'widget_dicts': [
         anode_label,
         cathode_label,
         hydrogen_label,
         nitrogen_label_1,
         oxygen_label,
         nitrogen_label_2,
         dry_gas_composition_dict,
         relative_humidity,
         gas_inlet_temperature,
         gas_outlet_pressure], 'size_label': 'l', 'size_unit': 's',
     'highlightbackground': 'grey', 'highlightthickness': 1}

coolant_control = \
    {'label': 'Coolant Control:', 'number': 1, 'width': 25,
     'sim_name': ['temperature_system', 'coolant_control'],
     'value': ['Temperature Difference', 'Mass Flow'],
     # 'columnspan': [1, 2],
     'type': 'ComboboxSet', 'specifier': 'dropdown_activate',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[1, 0]], [[2, 0]]],
                           [[[2, 0]], [[1, 0]]]]],
                 'args2': [1, 2]}
     }

coolant_mass_flow = \
    {'label': 'Coolant Mass Flow:', 'value': 1e-3,
     'sim_name': ['temperature_system', 'cool_mass_flow'],
     'dtype': 'float', 'dimensions': 'kg/s', 'type': 'EntrySet'}

coolant_mass_flow_frame = \
    {'title': 'Coolant Mass Flow:', 'show_title': False, 'sticky': 'NWE',
     'widget_dicts': [coolant_mass_flow],  'columnspan': 3,
     'specifier': 'visibility', 'padx': 0, 'pady': 0}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}

coolant_temp_diff = \
    {'label': 'Coolant Temperature Difference:', 'value': 10.0,
     'sim_name': ['temperature_system', 'cool_temp_diff'],
     'dtype': 'float', 'dimensions': 'K', 'type': 'EntrySet'}

coolant_temp_diff_frame = \
    {'title': 'Coolant Temperature Difference:', 'show_title': False,
     'sticky': 'NWE', 'columnspan': 3, 'specifier': 'visibility',
     'widget_dicts': [coolant_temp_diff],
     'padx': 0, 'pady': 0}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}

coolant_control_frame = \
    {'title': 'Coolant Control Settings', 'show_title': False, 'sticky': 'NWE',
     'columnspan': 3, 'size_label': 'l', 'size_unit': 's',
     'widget_dicts': [coolant_control,
                      coolant_temp_diff_frame,
                      coolant_mass_flow_frame],
     'highlightbackground': 'grey', 'highlightthickness': 1}

coolant_inlet_temperature = \
    {'label': 'Coolant Inlet Temperature:', 'number': 1, 'value': 343.15,
     'width': 11,  # 'columnspan': [1, 2, 2, 1],
     'sim_name': ['coolant_flow_circuit', 'inlet_manifold', 'temp_in'],
     'dtype': 'float', 'dimensions': 'K', 'type': 'EntrySet'}

coolant_outlet_pressure = \
    {'label': 'Coolant Outlet Pressure:', 'number': 1, 'value': 101325.0,
     'width': 11,  # 'columnspan': [1, 2, 2, 1],
     'sim_name':
         ['coolant_channel', 'p_out'],
     'dtype': 'float', 'dimensions': 'Pa', 'type': 'EntrySet'}

endplate_heat_flux = \
    {'label': 'Heat Flux at End Plates:', 'number': 1, 'value': 0.0,
     'sim_name': ['temperature_system', 'heat_flux'],
     'width': 11, 'dtype': 'float', 'dimensions': 'W/m²', 'type': 'EntrySet'}

convection_coefficient_environment = \
    {'label': 'Convection Coefficient to Environment:',
     'number': 1, 'value': 0.0, 'width': 11,
     'sim_name': ['temperature_system', 'alpha_amb'],
     'dtype': 'float', 'dimensions': 'W/(m²-K)', 'type': 'EntrySet'}

ambient_temperature = \
    {'label': 'Ambient Temperature:',
     'number': 1, 'value': 293.15, 'width': 11,
     'sim_name': ['temperature_system', 'temp_amb'],
     'dtype': 'float', 'dimensions': 'K', 'type': 'EntrySet'}

cooling_frame_dict = \
    {'title': 'Thermal Settings', 'show_title': False, 'sticky': 'NWE',
     'size_label': 'xl', 'size_unit': 'l',
     'widget_dicts': [
         coolant_control_frame,
         coolant_inlet_temperature,
         coolant_outlet_pressure,
         endplate_heat_flux,
         convection_coefficient_environment,
         ambient_temperature],
     'highlightbackground': 'grey', 'highlightthickness': 1}

operating_frame_dict = \
    {'title': 'Operating Conditions', 'show_title': False,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'widget_dicts': [electrochem_frame_dict, gas_frame_dict,
                      cooling_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = \
    {'title': 'Operating Conditions', 'show_title': False,
     'sub_frame_dicts': [operating_frame_dict]}
