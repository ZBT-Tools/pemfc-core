# general imports
import numpy as np
import sys
import copy
import matplotlib.pyplot as plt

# local modul imports
from pemfc import channel
from pemfc import fluid
from pemfc import species
from pemfc import interpolation as ip

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


nodes = 20
mass_flow_hydrogen = 0.0001    # kg/s
mass_flow_air = 0.0001
mass_flow_water = 0.002
wall_temp = None  # 380
# wall_temp_1 = 380.0
# wall_temp_2 = 420.0
heat_flux = 1000.0  # W/m²

inlet_temperature = 293.15
outlet_pressure = 101325.0

length = 0.1
width = 0.001
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
        species.ConstantProperties('Liquid',
                                   specific_heat=4000.0,
                                   density=1000.0,
                                   viscosity=1e-3,
                                   thermal_conductivity=0.6),
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
    'bend_number': 0,
    'bend_friction_factor': 500,
    'additional_friction_fractor': 0.01
    }

hydrogen = fluid.factory(hydrogen_dict)
air = fluid.factory(air_dict)
water = fluid.factory(water_dict)
fluids = [hydrogen, air, water]

channel_dicts = [copy.deepcopy(channel_dict) for i in range(3)]

channels = [channel.Channel(channel_dicts[i], fluids[i]) for i in range(3)]


error = 1e5
iter_max = 50
temp_old = np.asarray([channel.temperature for channel in channels])
mass_flows = [mass_flow_hydrogen, mass_flow_air, mass_flow_water]
delta_temp = 30.0
for i in range(iter_max):
    if error < 1e-4:
        break
    error = 0.0
    for j, channel in enumerate(channels):
        channel.update(mass_flow_in=mass_flows[j], heat_flux=heat_flux,
                       wall_temp=wall_temp)
        mass_flows[j] = np.sum(channel.heat) \
            / (np.average(channel.fluid.specific_heat) * delta_temp)
        error += np.sum(np.abs(((temp_old[j] - channel.temperature) / channel.temperature)))
        temp_old[j, :] = channel.temperature

x = ip.interpolate_1d(channels[0].x)
for channel in channels:
    plt.plot(x, channel.temp_ele,
             label='Fluid Temperature - ' + channel.fluid.name)
    plt.plot(x, channel.wall_temp,
             label='Wall Temperature - ' + channel.fluid.name)
plt.xlabel('Channel Location [m]')
plt.ylabel('Temperature [K]')
plt.legend()
plt.show()
for channel in channels:
    plt.plot(x, channel.heat,
             label='Heat Transfer - ' + channel.fluid.name)
plt.xlabel('Channel Location [m]')
plt.ylabel('Heat [W]')
plt.legend()
plt.show()

print('Pumping Power:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ',
          np.average(channel.vol_flow) * (channel.pressure[channel.id_in]
                                          - channel.pressure[channel.id_out]))
print('Mass Flows:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ', np.average(channel.mass_flow_total))

print('Pressure drop:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ', channel.pressure[channel.id_in]
          - channel.pressure[channel.id_out])

print('Flow Velocity:')
for i, channel in enumerate(channels):
    print(channel.fluid.name + ': ', channel.velocity)

