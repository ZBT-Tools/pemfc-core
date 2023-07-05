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
pressure = 250000.0
current_density = 1.0 * 10000
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
n_cells = 200
n_subchl = 40
in_manifold_dict = {
    'name': 'Inlet Manifold',
    'length': 1.435e-3 * n_cells,
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
            for i in range(n_cells)]

# Calculate mass flow
species_mass_flow, _ = half_cell.HalfCell.calc_faraday_flow(
    channels[0].fluid, current_density * active_area * n_cells, stoichiometry,
    reaction_stoichiometry=-1, charge_number=4, reactant_index=0)
mass_flow = np.sum(species_mass_flow)
# mass_flow = 0.00025275611864800167

# First flow circuit model configuration
flow_circuit_dict = {
    'name': 'Flow Circuit',
    'type': 'VariableResistance',
    'shape': 'U'
    }

flow_resistance_list = [{'type': 'Constant', 'value': 0.004},
                        # {'type': 'RennelsTeeMain', 'branch_diameter': 0.005},
                        {'type': 'BassettTeeMain', 'branch_diameter': 0.005},
                        {'type': 'IdelchikTeeMain', 'branch_diameter': 0.005},
                        {'type': 'HuangTeeMain', 'branch_diameter': 0.005}
                        ]

labels = [(item['type'].strip('TeeMain') + ' ' + str(item.get('value', ''))).strip()
          for item in flow_resistance_list]
colors = ['k', 'r', 'b', 'g']

fig, ax = plt.subplots()
flow_models = []
for i, flow_res_dict in enumerate(flow_resistance_list):
    in_manifold_dict['flow_resistances'] = [flow_res_dict]
    out_manifold_dict = copy.deepcopy(in_manifold_dict)
    out_manifold_dict['name'] = 'Outlet Manifold'

    flow_models.append(
        flow_circuit.create(flow_circuit_dict, in_manifold_dict,
                            out_manifold_dict, channels,
                            channel_multiplier=n_subchl))
    flow_model = flow_models[i]

    flow_model.update(inlet_mass_flow=mass_flow)
    q = flow_model.normalized_flow_distribution  # * 100.0
    reynolds = flow_model.manifolds[0].reynolds[0]
    ax.plot(q, label=labels[i], color=colors[i])
ax.set_xlabel('Cell Number / -')
ax.set_ylabel('Normalized Flow Distribution / -')

plt.legend()
plt.show()

fig, ax = plt.subplots()
for i, flow_model in enumerate(flow_models):
    m_in = flow_model.manifolds[0]
    m_out = flow_model.manifolds[1]
    ax.plot(m_in.pressure, color=colors[i], linestyle='solid',
            label='Inlet Manifold - ' + labels[i])
    ax.plot(m_out.pressure, color=colors[i], linestyle='dashed',
            label='Outlet Manifold - ' + labels[i])
ax.set_xlabel('Cell Number / -')
ax.set_ylabel('Pressure / Pa')
plt.legend()
plt.show()
