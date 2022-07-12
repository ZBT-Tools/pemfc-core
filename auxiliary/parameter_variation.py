import os
import json
import inspect
from pemfc.main_app import main
import numpy as np
import matplotlib.pyplot as plt

# Variation of membrane thickness
membrane_thickness_array = np.linspace(10e-6, 50e-6, 10)

# Load settings file
base_dir = \
    os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda: 0)))
with open(os.path.join(base_dir, '../pemfc/settings', 'settings.json')) as file:
    settings = json.load(file)
# Disable writing output to files to save simulation time
settings['output']['save_csv'] = False
settings['output']['save_plot'] = False

stack_power_array = []
global_data = None
# Loop for parameter variation simulation
for membrane_thickness in membrane_thickness_array:
    settings['membrane']['thickness'] = membrane_thickness
    global_data, local_data, sim = main(settings)
    stack_power_value = global_data[0]['Stack Power']['value']
    stack_power_array.append(stack_power_value)
power_unit = global_data[0]['Stack Power']['units']

# Plot data with matplotlib
fig, ax = plt.subplots()
ax.plot(membrane_thickness_array, stack_power_array)
ax.set_xlabel('Membrane Thickness / m')
ax.set_ylabel('Stack Power / ' + power_unit)
plt.show()
