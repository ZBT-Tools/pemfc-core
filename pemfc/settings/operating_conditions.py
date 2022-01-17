""" Operating conditions """
"""Electrochemistry settings"""
# if current_control = True, provided current_density value or array is used
# as operating point
# if current_control = False, provided average_cell_voltage value or array is
# used as operating point (not working as well, so current_control should be
# used at the moment)
current_control = True

# target current density [A/m^2]

# current_density = (100.0, 500.0, 1000.0, 1500.0, 2000.0, 5000.0, 10000.0,
#                    15000.0, 20000.0, 25000.0, 30000.0)

current_density = 20000.0

# average cell voltage [V] (not used, when current_control = True)
average_cell_voltage = 0.50

# open circuit voltage [V]
open_circuit_voltage = 1.0

# cathode stoichiometry
stoichiometry_cathode = 1.5

# anode stoichiometry
stoichiometry_anode = 1.5

# reaction educt stoichiometry
cathode_reaction_stoich = [-1.0, 0.0, 2.0]
anode_reaction_stoich = [-2.0, 0.0, 0.0]

# reaction electron stoichiometry
cathode_electron_number = 4.0
anode_electron_number = 4.0

"""Fluid settings"""
# components dictionary
# main fuel species (component for stoichiometry calculation)
# must be at first position
cathode_composition = \
    {'O2': {'state': 'gas', 'molar_fraction': 0.21},
     'N2': {'state': 'gas', 'molar_fraction': 0.79},
     'H2O': {'state': 'gas-liquid', 'molar_fraction': 0.0}}
anode_composition = \
    {'H2': {'state': 'gas', 'molar_fraction': 1.0},
     'N2': {'state': 'gas', 'molar_fraction': 0.0},
     'H2O': {'state': 'gas-liquid', 'molar_fraction': 0.0}}
# anode_species = {'H2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}
# anode_species = {'H2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'}

# relative humidity of inlet gases
# (overwrites H2O mole_fraction component above)
cathode_humidity = 0.5
anode_humidity = 0.5

# inlet composition (molar fractions)
# cathode_inlet_composition = [0.185, 0.697, 0.118]
# anode_inlet_composition = [1.0, 0.0, 0.0]

""""Thermal Settings"""
# temperature value used for most boundary conditions (see below and edit
# as required) [K]
temperature = 343.15

# # set the total mass flow of the coolant for the stack [kg/s]
# coolant_mass_flow = 1.e-1
# if coolant_mass_flow is not provided, set the desired coolant temperature
# difference [K] (mass flow overrides the temperature difference, if provided)
coolant_temperature_difference = 5.0

# coolant pinch point temperature difference [K]
temp_pinch = 0.0

# reactant inlet temperatures [K]
temp_cathode_in = temperature  # - coolant_temperature_difference - temp_pinch
# 443.15
temp_anode_in = temperature  # - coolant_temperature_difference - temp_pinch
# 443.15

# coolant inlet temperature [K]
temp_coolant_in = temperature \
    - coolant_temperature_difference * 0.5  # - temp_pinch

# environment temperature [K]
temp_environment = 293.15

# heat flux boundary conditions at the endplates [W/m²]
heat_flux_endplate = 0.0

# convection coefficient between the stack walls and the environment [W/(m^2K)]
convection_coefficient_environment = 10.0

# initial temperature [K]
temp_initial = temperature  # 443.15


"""Pressure settings"""
# pressure at the outlet of the manifolds [Pa]
p_manifold_cathode_out = 101325.0
# pressure at the outlet of the anode manifold [Pa]
p_manifold_anode_out = 101325.0
# pressure at the outlet of the cooling manifold [Pa]
p_manifold_coolant_out = 101325.0
