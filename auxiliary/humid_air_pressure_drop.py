from pemfc.src.fluid import fluid as fl
import pemfc.src.channel as chl

channel_dict = {
    "name": "Cathode Channel",
    "fluid": {
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
                "molar_fraction": 0.00
            }
        },
        'nodes': 21
        # "humidity": 0.75
    },
    "cross_sectional_shape": "rectangular",
    "base_width": 0.0008,
    "length": 0.5,
    "p_out": 101325.0,
    "temp_in": 343.15,
    "flow_direction": 1,
    "width": 0.001,
    "height": 0.001,
    "bend_number": 0,
    "bend_friction_factor": 0.1
}

humid_air = fl.create(channel_dict['fluid'])
channel = chl.Channel(channel_dict, humid_air)
mole_source_h2o = 1e-2
mass_source_h2o = mole_source_h2o * humid_air.species_mw[2]
if isinstance(channel, chl.TwoPhaseMixtureChannel):
    channel.mass_source[2][:] = mass_source_h2o
    # channel.mole_source[2][:] = mole_source_h2o
for i in range(5):
    channel.update(mass_flow_in=5e-5, update_mass=True, update_fluid=True,
                   update_flow=True)
print('test')