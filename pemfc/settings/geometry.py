""" Geometry settings """

"""Stack Settings"""
# number of cells in the stack
cell_number = 10

""""Cell Geometry """
# length of the cell, a side of the active area [m]
cell_length = 0.5
# height of the cell, b side of the active area [m]
cell_width = 0.1
# thickness of the bipolar plate [m]
cathode_bipolar_plate_thickness = 2.0e-3
anode_bipolar_plate_thickness = 2.0e-3
# thickness of the membrane [m]
membrane_thickness = 15.e-6
# thickness of the catalyst layer [m]
cathode_catalyst_layer_thickness = 10.e-6
anode_catalyst_layer_thickness = 10.e-6
# catalyst layer porosity (ionomer considered as porous volume) [-]
cathode_catalyst_layer_porosity = 0.5
anode_catalyst_layer_porosity = 0.5
# thickness of the gas diffusion layer [m]
cathode_gdl_thickness = 200.e-6
anode_gdl_thickness = 200.e-6
# gas diffusion layer porosity [-]
cathode_gdl_porosity = 0.8
anode_gdl_porosity = 0.8

"""Flow Field Geometry"""
# channel length [m]
cathode_channel_length = 0.5
anode_channel_length = 0.5
# channel width [m]
cathode_channel_width = 1.e-3
anode_channel_width = 1.e-3
# rib width [m]
cathode_channel_rib_width = cathode_channel_width
anode_channel_rib_width = anode_channel_width
# channel height [m]
cathode_channel_height = 1.e-3
anode_channel_height = 1.e-3
# number of channels
cathode_channel_number = 50
anode_channel_number = 50
# shape of channel (rectangular, trapezoidal, triangular)
cathode_channel_shape = 'rectangular'
anode_channel_shape = 'rectangular'
# width of channel at its base (only for trapezoidal) [m]
cathode_channel_base_width = 0.8 * cathode_channel_width
anode_channel_base_width = 0.8 * anode_channel_width
# channel bends [n]
cathode_channel_bends = 0
anode_channel_bends = 0
# bend pressure loss coefficient of the channel bends
bend_pressure_loss_coefficient = 0.1
# flow direction in channel along x-axis
cathode_flow_direction = 1
anode_flow_direction = -1

"""Coolant Settings"""
# set boolean to calculate coolant flow
coolant_circuit = True
# set boolean for coolant flow between the edge cells and endplates
cooling_bc = True
# channel length [m]
coolant_channel_length = 0.5
# height of the coolant channel [m]
coolant_channel_height = 1.e-3
# width of the coolant channel [m]
coolant_channel_width = 2.e-3
# coolant channel shape
cool_channel_shape = 'rectangular'
cool_channel_base_width = 0.8 * coolant_channel_width
# number of coolant channels per cell
coolant_channel_number = 20
# channel bends [n]
coolant_channel_bends = 0.0
# bend pressure loss coefficient of the channel bends
coolant_bend_pressure_loss_coefficient = 0.1

"""Manifold Geometry"""
# set boolean for calculation of the flow distribution
calc_cathode_distribution = True
calc_anode_distribution = True
calc_coolant_distribution = True

# Configuration: U- or Z-shape
anode_manifold_configuration = 'U'
cathode_manifold_configuration = 'U'
coolant_manifold_configuration = 'U'

# manifold cross-sectional shape ('circular', 'rectangular')
cathode_manifold_cross_shape = 'circular'
anode_manifold_cross_shape = 'circular'
coolant_manifold_cross_shape = 'circular'

# manifold diameters [m] (for circular shape)
cathode_in_manifold_diameter = 40e-3
cathode_out_manifold_diameter = 40e-3
anode_in_manifold_diameter = 40e-3
anode_out_manifold_diameter = 40e-3
coolant_in_manifold_diameter = 40e-3
coolant_out_manifold_diameter = 40e-3

# manifold height [m] (for rectangular shape)
cathode_in_manifold_height = 10e-3
cathode_out_manifold_height = 10e-3
anode_in_manifold_height = 10e-3
anode_out_manifold_height = 10e-3
coolant_in_manifold_height = 10e-3
coolant_out_manifold_height = 10e-3

# manifold width [m] (for rectangular shape)
cathode_in_manifold_width = 10e-3
cathode_out_manifold_width = 10e-3
anode_in_manifold_width = 10e-3
anode_out_manifold_width = 10e-3
coolant_in_manifold_width = 10e-3
coolant_out_manifold_width = 10e-3

# geometrical pressure loss coefficient of the manifolds
general_loss_coefficient = 0.4
cathode_in_manifold_pressure_loss_coefficient = general_loss_coefficient
cathode_out_manifold_pressure_loss_coefficient = general_loss_coefficient
anode_in_manifold_pressure_loss_coefficient = general_loss_coefficient
anode_out_manifold_pressure_loss_coefficient = general_loss_coefficient
coolant_in_manifold_pressure_loss_coefficient = general_loss_coefficient
coolant_out_manifold_pressure_loss_coefficient = general_loss_coefficient

# Model type
model_name = 'UpdatedKoh'
anode_manifold_model = model_name
cathode_manifold_model = model_name
coolant_manifold_model = model_name
# anode_manifold_model = 'Koh'
# cathode_manifold_model = 'Koh'
# coolant_manifold_model = 'Koh'
