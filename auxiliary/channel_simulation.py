# general imports
import numpy as np
import sys
import copy
import matplotlib.pyplot as plt

# local modul imports
from pemfc.src import channel as chl
from pemfc.src import fluid as fluids
from pemfc.src import species as species
from pemfc.src import interpolation as ip

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


nodes = 50
mass_flow_hydrogen = 0.0001    # kg/s
mass_flow_air = 0.0000124612
mass_flow_water = 0.0000124612
wall_temp = None  # 380
# wall_temp_1 = 380.0
# wall_temp_2 = 420.0
heat_flux = 0.0    # W/m²
inlet_temperature = 293.15
outlet_pressure = 101331.0

length = 0.65501
width = 0.004
diameter = 0.0016
height = 0.001

hydrogen_dict = {
    'name': 'Hydrogen',
    'fluid_components': {'H2': 'gas'},
    'inlet_composition': 1.0,
    'temp_init': inlet_temperature,
    'press_init': outlet_pressure,
    'nodes': nodes
}

air_dict = {
    'name': 'Air',
    'fluid_components': {'O2': 'gas', 'N2': 'gas'},
    'inlet_composition': [0.21, 0.79],
    'temp_init': inlet_temperature,
    'press_init': outlet_pressure,
    'nodes': nodes
}

water_dict = {
    'name': 'Water',
    'fluid_components': None,
    'inlet_composition': None,
    'liquid_props':
        species.ConstantProperties('Air',
                                   specific_heat=1007.0,
                                   density=1.20433,
                                   viscosity=1.824e-05,
                                   thermal_conductivity=0.0257),
    'temp_init': inlet_temperature,
    'press_init': outlet_pressure,
    'nodes': nodes
    }


channel_dict = {
    'name': 'Channel',
    'length': length,
    'cross_sectional_shape': 'rectangular',
    'p_out': outlet_pressure,
    'temp_in': inlet_temperature,
    'flow_direction': 1,
    'width': width,
    'height': height,
    'diameter': diameter,
    'bend_number': 0,
    'bend_friction_factor': 500.0,
    'constant_friction_fractor': 0.0
    }

hydrogen = fluids.factory(hydrogen_dict)
air = fluids.factory(air_dict)
water = fluids.factory(water_dict)
fluids = [hydrogen, air, water]

channel_dicts = [copy.deepcopy(channel_dict) for i in range(3)]

channels = [chl.Channel(channel_dicts[i], fluids[i]) for i in range(3)]


error = 1e5
iter_max = 50
mass_flows = [mass_flow_hydrogen, mass_flow_air, mass_flow_water]

for j, channel in enumerate(channels):
    channel.update(mass_flow_in=mass_flows[j])

x = ip.interpolate_1d(channels[0].x)
for channel in channels:
    plt.plot(x, ip.interpolate_1d(channel.pressure),
             label='Fluid Pressure - ' + channel.fluid.name)
plt.xlabel('Channel Location [m]')
plt.ylabel('Pressure [Pa]')
plt.legend()
plt.show()


print('Pressure:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ', channel.pressure)

print('Pressure drop:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ',
          channel.pressure[channel.id_in] - channel.pressure[channel.id_out])

print('Flow Velocity:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ', channel.velocity)

