import os
import json
import inspect
import sys
import numpy as np
import matplotlib.pyplot as plt
try:
    from import_pemfc import main
except:
    from .import_pemfc import main


# Variation of HFR
hfr_factor = [1.0, 1.5/10.0]
# membrane_thickness_array = np.linspace(10e-6, 30e-6, 3)
# membrane_conductivity_array = np.asarray([15.0, 15.0/6.67])

# Current density array (for simulation of polarization curve)
current_density_array = np.linspace(100, 20000, 40)

# Load settings file
base_dir = \
    os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda: 0)))
with open(os.path.join(base_dir, '../pemfc/settings', 'settings.json'), "r") as file:
    settings = json.load(file)

# Disable writing output to files to save simulation time
settings['output']['save_csv'] = False
settings['output']['save_plot'] = False
# Adjust stack parameters
settings['stack']['cell_number'] = 1
settings['membrane']['type'] = 'Constant'


# Initialize result stores
stack_power_list = []
average_voltage_list = []
global_data = None
membrane_conductivity = settings['membrane']['ionic_conductivity']
cat_bpp_conductivity = settings['cathode']["bpp"]['electrical_conductivity']
cat_gde_conductivity = settings['cathode']["gde"]['electrical_conductivity']
ano_bpp_conductivity = settings['anode']["bpp"]['electrical_conductivity']
ano_gde_conductivity = settings['anode']["gde"]['electrical_conductivity']

# Loop for parameter variation simulation
for i in range(len(hfr_factor)):
    # Adjust parameters
    settings['membrane']['ionic_conductivity'] = \
        membrane_conductivity * hfr_factor[i]
    settings['cathode']['electrical_conductivity_bpp'] = \
        [item * hfr_factor[i] for item in cat_bpp_conductivity]
    settings['cathode']['electrical_conductivity_gde'] = \
        [item * hfr_factor[i] for item in cat_gde_conductivity]
    settings['anode']['electrical_conductivity_bpp'] = \
        [item * hfr_factor[i] for item in ano_bpp_conductivity]
    settings['anode']['electrical_conductivity_gde'] = \
        [item * hfr_factor[i] for item in ano_gde_conductivity]
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
labels = ['HFR Factor: {:.2f}'.format(1.0/item) for item in hfr_factor]
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
