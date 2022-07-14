import os
import json
import inspect
from pemfc.main_app import main
import numpy as np
import matplotlib.pyplot as plt

# Variation of membrane thickness
membrane_thickness_array = np.linspace(10e-6, 30e-6, 3)

# Current density array (for simulation of polarization curve)
current_density_array = np.linspace(1000, 20000, 10)

# Load settings file
base_dir = \
    os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda: 0)))
with open(os.path.join(base_dir, '../pemfc/settings', 'settings.json')) as file:
    settings = json.load(file)

# Disable writing output to files to save simulation time
settings['output']['save_csv'] = False
settings['output']['save_plot'] = False
# Adjust stack parameters
settings['stack']['cell_number'] = 1

# Initialize result stores
stack_power_list = []
average_voltage_list = []
global_data = None

# Loop for parameter variation simulation
for membrane_thickness in membrane_thickness_array:
    # Adjust membrane thickness in settings dictionary
    settings['membrane']['thickness'] = membrane_thickness
    # Set current density array for polarization curve simulation
    settings['simulation']['current_density'] = current_density_array

    # Run simulation
    global_data, local_data, sim = main(settings)
    power_curve = []
    voltage_curve = []
    for operating_point in global_data:
        power_curve.append(operating_point['Stack Power']['value'])
        voltage_curve.append(operating_point['Average Cell Voltage']['value'])
    stack_power_list.append(power_curve)
    average_voltage_list.append(voltage_curve)
power_unit = global_data[0]['Stack Power']['units']
voltage_unit = global_data[0]['Average Cell Voltage']['units']

# Plot data with matplotlib
fig, ax = plt.subplots()
colors = ['k', 'r', 'b']
linestyles = ['solid', 'dashed']
labels = ['{:.2f} µm'.format(item * 1e6) for item in membrane_thickness_array]
for i in range(len(average_voltage_list)):
    ax.plot(current_density_array / 1e4, average_voltage_list[i],
            label=labels[i], color=colors[i], linestyle=linestyles[0])
ax2 = ax.twinx()
for i in range(len(stack_power_list)):
    ax2.plot(current_density_array / 1e4, stack_power_list[i],
             label=labels[i], color=colors[i], linestyle=linestyles[1])
ax.set_xlabel('Current Density / A/cm²')
ax.set_ylabel('Average Cell Voltage / ' + voltage_unit)
ax2.set_ylabel('Stack Power / ' + power_unit)
ax.legend(loc='lower center')
plt.show()
