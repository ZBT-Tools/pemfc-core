"""
Comparison with paper:
Huang, Fuxiang, Diankai Qiu, Zhutian Xu, Linfa Peng, and Xinmin Lai.
“Analysis and Improvement of Flow Distribution in Manifold for Proton Exchange
Membrane Fuel Cell Stacks.” Energy 226 (July 2021): 120427.
https://doi.org/10.1016/j.energy.2021.120427.
"""

# general imports
import numpy as np
import copy
import sys
import matplotlib.pyplot as plt
import matplotlib

# local module imports
from pemfc import channel as chl
from pemfc import fluid as fluid
from pemfc import flow_circuit as flow_circuit
from pemfc import interpolation as ip
from pemfc import constants
from pemfc import half_cell

matplotlib.use('TkAgg')

np.set_printoptions(threshold=sys.maxsize, linewidth=10000,
                    precision=9, suppress=True)
np.seterr(all='raise')


# Operating conditions
temperature = 348
pressure = 150000.0
current_density = 1.4 * 10000
stoichiometry = 2.0

# Channel discretization
nodes = 20

# Channel/cell configuration
active_area = 250 / 10000.0
channel_dict = {
    'name': 'Channel',
    'length': 0.269,
    'cross_sectional_shape': 'rectangular',
    'width': 8.5e-4,
    'height': 4e-4,
    'p_out': pressure,
    'temp_in': temperature,
    'flow_direction': 1,
    'bend_number': 0,
    # 'bend_friction_factor': 0.1,
    # 'constant_friction_factor': 0.2
    }

fluid_dict = {
    "name": "Cathode Gas Mixture",
    "components": {
        "O2": {
            "state": "gas",
            "molar_fraction": 0.21
        },
        "N2": {
            "state": "gas",
            "molar_fraction": 0.79
        },
        "H2O": {
            "state": "gas-liquid",
            "molar_fraction": 0.0
        }
    },
    "humidity": 1.0,
    'temp_init': temperature,
    'press_init': pressure,
    'nodes': nodes
}

# Stack configuration
n_chl = 200
n_subchl = 40
in_manifold_dict = {
    'name': 'Inlet Manifold',
    'length': 1.435e-3 * n_chl,
    'cross_sectional_shape': 'circular',
    'p_out': channel_dict['p_out'],
    'temp_in': temperature,
    'flow_direction': 1,
    'diameter': 23e-3,
    'bend_number': 0,
    'bend_friction_factor': 0.0,
    'flow_split_factor': 0.0,
    'wall_friction': True
}

# Create channels instances
channels = [chl.Channel(channel_dict, fluid.create(fluid_dict))
            for i in range(n_chl)]

# Calculate mass flow
species_mass_flow, _ = half_cell.HalfCell.calc_faraday_flow(
    channels[0].fluid, current_density * active_area, stoichiometry,
    reaction_stoichiometry=-1, charge_number=4, reactant_index=0)
mass_flow = np.sum(species_mass_flow)
mass_flow = 0.00025275611864800167

# First flow circuit model configuration
flow_circuit_dict = {
    'name': 'Flow Circuit',
    'type': 'UpdatedKoh',
    'shape': 'U'
    }

in_manifold_dict['flow_resistances'] = [{'type': 'Constant', 'value': 0.4}]
out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'

flow_model = \
    flow_circuit.create(flow_circuit_dict, in_manifold_dict,
                        out_manifold_dict, channels,
                        channel_multiplier=n_subchl)

x = (ip.interpolate_1d(flow_model.manifolds[0].x)
     - flow_model.manifolds[0].dx * 0.5) \
    / (flow_model.manifolds[0].length - flow_model.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length


flow_model.update(inlet_mass_flow=mass_flow)
q = flow_model.normalized_flow_distribution  # * 100.0
reynolds = flow_model.manifolds[0].reynolds[0]
plt.plot(x, q, label='Re={0:.2f} Model, Koh-Model'.format(reynolds), color='k')

# Another flow circuit model configuration
in_manifold_dict['flow_resistances'] = []
out_manifold_dict = copy.deepcopy(in_manifold_dict)
out_manifold_dict['name'] = 'Outlet Manifold'
flow_circuit_dict['type'] = 'VariableResistance'


flow_model_2 = \
    flow_circuit.create(flow_circuit_dict, in_manifold_dict,
                        out_manifold_dict, channels,
                        channel_multiplier=n_subchl)

x = (ip.interpolate_1d(flow_model_2.manifolds[0].x)
     - flow_model_2.manifolds[0].dx * 0.5) \
    / (flow_model_2.manifolds[0].length - flow_model_2.manifolds[0].dx[0])
# x = flow_model.manifolds[0].x / flow_model.manifolds[0].length

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
