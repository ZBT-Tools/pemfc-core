# general imports
import numpy as np
import cProfile
import copy
import sys
import matplotlib.pyplot as plt
import matplotlib

# local module imports
from pemfc import channel as chl
from pemfc import fluid as fluid
from pemfc import flow_circuit as flow_circuit
from pemfc import interpolation as ip

matplotlib.use('TkAgg')

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


def do_c_profile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats('cumtime')
    return profiled_func


n_chl = 50
n_subchl = 1

temperature = 293.15
pressure = 101325.0
nodes = 10

channel_dict = {
    'name': 'Channel',
    'length': 0.400,
    'cross_sectional_shape': 'circular',
    'width': 2.4e-3,
    'height': 2e-3,
    'diameter': 0.003,
    'p_out': pressure,
    'temp_in': temperature,
    'flow_direction': 1,
    'bend_number': 0,
    'bend_friction_factor': 0.1,
    'constant_friction_factor': 0.2
    }


fluid_dict = {
    'name': 'Cathode Gas',
    'fluid_components': {'O2': 'gas', 'N2': 'gas'},
    'inlet_composition': [0.21, 0.79],
    # 'temp_init': temperature,
    # 'press_init': pressure,
    'nodes': nodes
}

constant_air_dict = {
    'name': 'Air',
    'fluid_components': None,
    'inlet_composition': None,
    'specific_heat': 1007.0,
    'density': 1.20433,
    'viscosity': 1.824e-05,
    'thermal_conductivity': 0.0257,
    'temp_init': temperature,
    'press_init': pressure,
    'nodes': nodes
    }
constant_water_dict = {
    'name': 'Air',
    'fluid_components': None,
    'inlet_composition': None,
    'specific_heat': 4180,
    'density': 997.13,
    'viscosity': 0.000891,
    'thermal_conductivity': 0.0257,
    'temp_init': temperature,
    'press_init': pressure,
    'nodes': nodes
    }

in_manifold_dict = {
    'name': 'Inlet Manifold',
    'length': 0.14,
    'p_out': channel_dict['p_out'],
    'temp_in': 293.15,
    'flow_direction': 1,
    'width': 25e-3,
    'height': 25e-3,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'constant_friction_factor': 0.4,
    'flow_split_factor': 0.0,
    'wall_friction': True
}

out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'

flow_circuit_dict = {
    'name': 'Flow Circuit',
    'type': 'UpdatedKoh',
    'shape': 'U'
    }

channels = [chl.Channel(channel_dict, fluid.create(constant_water_dict))
            for i in range(n_chl)]

flow_model = \
    flow_circuit.create(flow_circuit_dict, in_manifold_dict,
                        out_manifold_dict, channels,
                        channel_multiplier=n_subchl)


x = (ip.interpolate_1d(flow_model.manifolds[0].x)
     - flow_model.manifolds[0].dx * 0.5) \
    / (flow_model.manifolds[0].length - flow_model.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length

vol_flow_liter_per_min = 3.2
vol_flow = vol_flow_liter_per_min / 1000.0 / 60.0
mass_flow = vol_flow * constant_water_dict['density']

flow_model.update(inlet_mass_flow=mass_flow)
q = flow_model.normalized_flow_distribution  # * 100.0
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f} Model, Koh-Model'.format(reynolds), color='k')

in_manifold_dict['constant_friction_factor'] = 0.0
out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'
flow_circuit_dict['type'] = 'VariableResistance'

channels = [chl.Channel(channel_dict, fluid.create(constant_water_dict))
            for i in range(n_chl)]

flow_model_2 = \
    flow_circuit.create(flow_circuit_dict, in_manifold_dict,
                        out_manifold_dict, channels,
                        channel_multiplier=n_subchl)


x = (ip.interpolate_1d(flow_model_2.manifolds[0].x)
     - flow_model_2.manifolds[0].dx * 0.5) \
    / (flow_model_2.manifolds[0].length - flow_model_2.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length

vol_flow_liter_per_min = 3.2
vol_flow = vol_flow_liter_per_min / 1000.0 / 60.0
mass_flow = vol_flow * constant_water_dict['density']

flow_model_2.update(inlet_mass_flow=mass_flow)
q = flow_model_2.normalized_flow_distribution  # * 100.0
reynolds = flow_model_2.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f} Model, VariableResistance'.format(reynolds), color='b')

# ref_data_dir = r'D:\Software\Python\PycharmProjects\CFDManifoldAnalyzer\data\Manifold_Distribution_Dummy-Stack_OPIAT.csv'
# ref_data_raw = np.loadtxt(ref_data_dir, delimiter=';').transpose()
# x = ref_data_raw[0] / np.max(ref_data_raw[0])
# q = (ref_data_raw[1] - 1.0) * 100.0
# plt.plot(x, q, label='Re={0:.2f} CFD Reference'.format(reynolds), color='r')

plt.legend()
plt.show()

# flow_model.update(inlet_mass_flow=0.00059425)
# q = (flow_model.normalized_flow_distribution - 1.0) * 100.0
# reynolds = flow_model.manifolds[0].reynolds[0]
# plt.plot(x, q, label='Re={0:.2f}'.format(reynolds), color='b')
#
# flow_model.update(inlet_mass_flow=0.000297125)
# q = (flow_model.normalized_flow_distribution - 1.0) * 100.0
# reynolds = flow_model.manifolds[0].reynolds[0]
# plt.plot(x, q, label='Re={0:.2f}'.format(reynolds), color='r')

# print('Normalized Flow Distribution: ',
#       flow_model.normalized_flow_distribution)
# np.savetxt('output/flow_distribution.txt',
#            (flow_model.normalized_flow_distribution - 1.0) * 100.0)
m_in = flow_model.manifolds[0]
# plt.plot(m_in.pressure - 101325.0, color='b')
# plt.show()
m_out = flow_model.manifolds[1]
# plt.plot(m_out.pressure - 101325.0, color='r')
# plt.show()
m_in_2 = flow_model_2.manifolds[0]
# plt.plot(m_in.pressure - 101325.0, color='b')
# plt.show()
m_out_2 = flow_model_2.manifolds[1]
plt.plot(m_in.pressure, color='k', linestyle='solid', label='Inlet Manifold - Koh')
plt.plot(m_out.pressure, color='k', linestyle='dashed', label='Outlet Manifold - Koh')
plt.plot(m_in_2.pressure, color='r', linestyle='solid', label='Inlet Manifold - Variable Resistance')
plt.plot(m_out_2.pressure, color='r', linestyle='dashed', label='Outlet Manifold - Variable Resistance')
plt.show()
