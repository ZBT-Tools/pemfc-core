""" Simulation Settings"""

"""Numerical Settings"""
# discretization of the flow channel along the x-axis
elements = 10
# apply channel-land-discretization
# (higher accuracy, slower calculation)
channel_land_discretization = False
# convergence criteria of the simulation
convergence_criteria = 1.e-7
convergence_criteria_flow = 1e-7
# maximum number of iterations
minimum_iteration_number = 3
minimum_iteration_number_flow = 3
maximum_iteration_number = 100
maximum_iteration_number_flow = 100
# underrelaxation factor for current density updates (0.0 - 1.0)
# lower: faster convergence, higher: better stability
underrelaxation_factor = 0.7
underrelaxation_factor_flow = 0.3
# minimum number of iterations
# numerical concentration value for determining critical current density for
# linearization/regularization of Kulikovsky model at limiting current densities
c_eps = 0.1
# difference of current densities used for linear gradient determination
# for linearization/regularization at limiting current densities
delta_i = 5.0
# calculate the PEMFC stack temperatures
calc_temperature = True
# calculate the current density distribution
calc_current_density = True
# calculate the flow distribution (setting moved to individual flow circuits)
# calc_flow_distribution = True
# calculate the activation voltage losses
calc_activation_loss = True
# calculate the membrane voltage losses
calc_membrane_loss = True
# calculate gdl diffusion voltage losses
calc_gdl_loss = True
# calculate catalyst layer diffusion voltage losses
calc_cl_loss = True


