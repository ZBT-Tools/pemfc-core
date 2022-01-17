from pemfc.src.global_functions import sub, sup

# Thermal conductivities
anode_label = {'label': 'Anode', 'row': 1, 'column': 1, 'columnspan': 2,
               'pady': 0, 'type': 'Label', 'sticky': 'WENS'}
cathode_label = {'label': 'Cathode', 'row': 1, 'column': 3, 'columnspan': 2,
                 'pady': 0, 'type': 'Label', 'sticky': 'WENS'}
in_plane_label_1 = {'label': 'ip', 'row': 2, 'column': 1, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
through_plane_label_1 = {'label': 'tp', 'row': 2, 'column': 2, 'pady': 0,
                         'type': 'Label', 'sticky': 'WENS'}
in_plane_label_2 = {'label': 'ip', 'row': 2, 'column': 3, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
through_plane_label_2 = {'label': 'tp', 'row': 2, 'column': 4, 'pady': 0,
                         'type': 'Label', 'sticky': 'WENS'}
electrical_conductivity_bpp = \
    {'label': 'BPP Electrical Conductivity:',
     'number': 4, 'value': 6e4, 'width': 5, 'types': 'multiinput',
     'sim_name': [['anode', 'electrical_conductivity_bpp', [0, 1]],
                  ['cathode', 'electrical_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

electrical_conductivity_gde = \
    {'label': 'GDE Electrical Conductivity:',
     'number': 4, 'value': 500.0, 'width': 5, 'types': 'multiinput',
     'sim_name': [['anode', 'electrical_conductivity_gde', [0, 1]],
                  ['cathode', 'electrical_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

thermal_conductivity_bpp = \
    {'label': 'BPP Thermal Conductivity:', 'types': 'multiinput',
     'value': [[100.0, 100.0], [100.0, 100.0]], 'width': 5,
     'sim_name': [['anode', 'thermal_conductivity_bpp', [0, 1]],
                  ['cathode', 'thermal_conductivity_bpp', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

thermal_conductivity_gde = \
    {'label': 'GDE Thermal Conductivity:', 'types': 'multiinput',
     'value': [[1.0, 2.0], [3.0, 4.0]], 'width': 5,
     'sim_name': [['anode', 'thermal_conductivity_gde', [0, 1]],
                  ['cathode', 'thermal_conductivity_gde', [2, 3]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

# empty_row = {'label': ' ',  'font': 'Arial 1',  # 'row': 1, 'column': 1,
#              'type': 'Label', 'sticky': 'WENS'}

porosity_gdl = \
    {'label': 'GDL Porosity:', 'number': 2, 'value': 0.8,
     'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'porosity_gdl'], ['cathode', 'porosity_gdl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}
porosity_cl = \
    {'label': 'Catalyst Porosity:', 'number': 2, 'value': 0.5,
     'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'porosity_cl'], ['cathode', 'porosity_cl']],
     'dtype': 'float', 'dimensions': '-', 'type': 'EntrySet'}

porous_frame_dict = \
    {'title': 'Porous Layers', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'size_label': 's', 'size_unit': 'l',
     'widget_dicts': [anode_label,
                      cathode_label,
                      in_plane_label_1, through_plane_label_1,
                      in_plane_label_2, through_plane_label_2,
                      electrical_conductivity_bpp,
                      electrical_conductivity_gde,
                      thermal_conductivity_bpp,
                      thermal_conductivity_gde,
                      # empty_row,
                      porosity_gdl, porosity_cl]}
     #'highlightbackground': 'grey', 'highlightthickness': 1}


membrane_model = \
    {'label': 'Membrane Model:', 'number': 1,
     'sim_name': ['membrane', 'type'],
     'value': ['Constant', 'Springer', 'Linear'],
     'type': 'ComboboxSet', 'specifier': 'dropdown_activate',
     'command': {'function': 'show_connected_widgets',
                 'args': [[[[[1, 0]], [[2, 0], [3, 0]]],
                           [[[2, 0]], [[1, 0], [3, 0]]],
                           [[[3, 0]], [[1, 0], [2, 0]]]]],
                 }
     }

mem_ionic_conductivity = \
    {'label': 'Ionic Conductivity:', 'value': 5.0,
     'sim_name': ['membrane', 'ionic_conductivity'],
     'dtype': 'float', 'dimensions': 'S/m', 'type': 'EntrySet'}

empty_row = {'label': ' ', 'type': 'Label', 'sticky': 'WENS'}

constant_frame = \
    {'title': 'Constant Ionic Conductivity', 'specifier': 'visibility',
     'widget_dicts': [mem_ionic_conductivity,
                      empty_row],
     'sticky': 'WEN', 'columnspan': 2, 'padx': 0, 'pady': 0}

mem_constant_resistance = \
    {'label': 'Constant Resistance Coefficient:', 'value': 4.3e-5,
     'sim_name': ['membrane', 'basic_resistance'],
     'dtype': 'float', 'dimensions': 'Ohm-m²', 'type': 'EntrySet'}

mem_temp_coefficient = \
    {'label': 'Linear Temperature Coefficient:', 'value': 7e-8,
     'sim_name': ['membrane', 'temperature_coefficient'],
     'dtype': 'float', 'dimensions': 'Ohm-m²/K', 'type': 'EntrySet'}

linear_frame = \
    {'title': 'Linear Ionic Conductivity', 'specifier': 'visibility',
     'widget_dicts': [mem_constant_resistance, mem_temp_coefficient],
     'sticky': 'WEN', 'columnspan': 2, 'padx': 0, 'pady': 0}

mem_acid_group_conc = \
    {'label': 'Acid Group Concentration:', 'value': 1.2e3,
     'sim_name': ['membrane', 'acid_group_concentration'],
     'dtype': 'float', 'dimensions': 'mol/m³', 'type': 'EntrySet'}

mem_vapour_transport_coefficient = \
    {'label': 'Vapour Transport Coefficient:', 'value': 6.2e-6,
     'sim_name': ['membrane', 'vapour_transport_coefficient'],
     'dtype': 'float', 'dimensions': 'm/s', 'type': 'EntrySet'}

springer_frame = \
    {'title': 'Springer Ionic Conductivity', 'specifier': 'visibility',
     'widget_dicts': [empty_row, empty_row], 'padx': 0, 'pady': 0,
     'sticky': 'WEN', 'columnspan': 2}

membrane_model_frame_dict = \
    {'title': 'Membrane Model Settings', 'show_title': False,
     'font': 'Arial 10 bold', 'sticky': 'WEN', 'padx': 0, 'pady': 0,  #'row': 1,
     'size_label': 'l', 'size_unit': 'l',
     'widget_dicts': [
                      membrane_model,
                      constant_frame,
                      springer_frame,
                      linear_frame,
                      ]}
     # 'highlightbackground': 'grey', 'highlightthickness': 1}

in_plane_label_3 = {'label': 'ip', 'row': 2, 'column': 1, 'pady': 0,
                    'type': 'Label', 'sticky': 'WENS'}
through_plane_label_3 = {'label': 'tp', 'row': 2, 'column': 2, 'pady': 0,
                         'type': 'Label', 'sticky': 'WENS'}
thermal_conductivity_mem = \
    {'label': 'Membrane Thermal Conductivity:',
     'value': [[0.26, 0.26]], 'types': 'multiinput',
     #'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['membrane', 'thermal_conductivity', [0, 1]]],
     'dtype': 'float', 'dimensions': 'W/(m K)', 'type': 'EntrySet'}

membrane_thermal_frame_dict = \
    {'title': 'Membrane Thermal Properties', 'show_title': False,
     'font': 'Arial 10 bold', 'sticky': 'WEN', #'row': 2,
     'size_label': 'l', 'size_unit': 'l',
     'widget_dicts': [
                      in_plane_label_3,
                      through_plane_label_3,
                      thermal_conductivity_mem,
                      ]}
     #'highlightbackground': 'grey', 'highlightthickness': 1}

membrane_frame_dict = \
    {'title': 'Membrane Settings', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'sub_frame_dicts': [
                         membrane_model_frame_dict,
                         membrane_thermal_frame_dict,
                        ]}
     #'highlightbackground': 'grey', 'highlightthickness': 1}

anode_label_2 = {'label': 'Anode', 'row': 1, 'column': 1,
                 'pady': 0, 'type': 'Label', 'sticky': 'WENS'}
cathode_label_2 = {'label': 'Cathode', 'row': 1, 'column': 2,
                   'pady': 0, 'type': 'Label', 'sticky': 'WENS'}

exchange_current_density = \
    {'label': 'Exchange Current Density:', 'value': [5.0e8, 8.0e5],
     #'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'vol_ex_cd'], ['cathode', 'vol_ex_cd']],
     'dtype': 'float', 'dimensions': 'A/m³', 'type': 'EntrySet'}

reactant_gdl_diffusion_coefficient = \
    {'label': 'Effective GDL Diffusion Coefficient:', 'value': [10e-6, 6e-6],
     #'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'diff_coeff_gdl'], ['cathode', 'diff_coeff_gdl']],
     'dtype': 'float', 'dimensions': 'm²/s', 'type': 'EntrySet'}

reactant_catalyst_diffusion_coefficient = \
    {'label': 'Effective Catalyst Diffusion Coefficient:',
     'value': [1e-7, 1e-7],
     #'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'diff_coeff_cl'], ['cathode', 'diff_coeff_cl']],
     'dtype': 'float', 'dimensions': 'm²/s', 'type': 'EntrySet'}

catalyst_proton_conductivity = \
    {'label': 'Effective Catalyst Proton Conductivity:', 'value': [1.5, 1.5],
     # 'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'prot_con_cl'], ['cathode', 'prot_con_cl']],
     'dtype': 'float', 'dimensions': '1/(Ohm-m)', 'type': 'EntrySet'}

tafel_slope = \
    {'label': 'Tafel Slope:', 'value': [0.035, 0.035],
     # 'width': 11, 'columnspan': [1, 2, 2, 1],
     'sim_name': [['anode', 'tafel_slope'], ['cathode', 'tafel_slope']],
     'dtype': 'float', 'dimensions': 'V', 'type': 'EntrySet'}

open_circuit_voltage = \
    {'label': 'Open Circuit Voltage:', 'value': 1.0,
     'width': 21, 'columnspan': [1, 2, 1],
     'sim_name': ['cell', 'open_circuit_voltage'],
     'dtype': 'float', 'dimensions': 'V', 'type': 'EntrySet'}

electrochemistry_frame_dict = \
    {'title': 'Electrochemistry Model', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'size_label': 'l', 'size_unit': 'xl',
     'widget_dicts': [anode_label_2, cathode_label_2,
                      exchange_current_density,
                      reactant_gdl_diffusion_coefficient,
                      reactant_catalyst_diffusion_coefficient,
                      catalyst_proton_conductivity,
                      tafel_slope,
                      open_circuit_voltage]}
     #'highlightbackground': 'grey', 'highlightthickness': 1}

coolant_density = \
    {'label': 'Density:', 'value': 981.0,
     # 'width': 21, 'columnspan': [1, 2, 1],
     'sim_name': ['coolant_channel', 'fluid', 'density'],
     'dtype': 'float', 'dimensions': 'kg/m³', 'type': 'EntrySet'}

coolant_viscosity = \
    {'label': 'Viscosity:', 'value': 9.81e-3,
     # 'width': 21, 'columnspan': [1, 2, 1],
     'sim_name': ['coolant_channel', 'fluid', 'viscosity'],
     'dtype': 'float', 'dimensions': 'Pa-s', 'type': 'EntrySet'}

coolant_specific_heat = \
    {'label': 'Specific Heat Capacity:', 'value': 2010.0,
     # 'width': 21, 'columnspan': [1, 2, 1],
     'sim_name': ['coolant_channel', 'fluid', 'specific_heat'],
     'dtype': 'float', 'dimensions': 'J/(kg-K)', 'type': 'EntrySet'}

coolant_thermal_conductivity = \
    {'label': 'Thermal Conductivity:', 'value': 0.159,
     # 'width': 21, 'columnspan': [1, 2, 1],
     'sim_name': ['coolant_channel', 'fluid', 'thermal_conductivity'],
     'dtype': 'float', 'dimensions': 'W/(m-K)', 'type': 'EntrySet'}

coolant_properties_frame_dict = \
    {'title': 'Coolant Properties', 'show_title': True,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'size_label': 'xl', 'size_unit': 'l',
     'widget_dicts': [coolant_density,
                      coolant_viscosity,
                      coolant_specific_heat,
                      coolant_thermal_conductivity]}
     #'highlightbackground': 'grey', 'highlightthickness': 1}

frame_dict = \
    {'title': 'Physical Properties', 'show_title': False,
     'font': 'Arial 10 bold', 'sticky': 'WEN',
     'sub_frame_dicts': [porous_frame_dict,
                         membrane_frame_dict,
                         electrochemistry_frame_dict,
                         coolant_properties_frame_dict],
     'highlightbackground': 'grey', 'highlightthickness': 1}

tab_dict = {'title': 'Physical Properties', 'show_title': False,
            'sub_frame_dicts': [frame_dict]}
