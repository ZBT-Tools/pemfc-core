import numpy as np
import pemfc.src.fluid.fluid as fluid
import pemfc.src.fluid.diffusion_model as diffusion_model
import time

fluid_dict = \
    {
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
        "humidity": 0.5,
        "temperature": 333.15,
        "pressure": 101325.0,
        "nodes": (1)
    }

humid_air = fluid.factory(fluid_dict, backend='pemfc')
humid_air.update()
surface_tension = humid_air.phase_change_species.calc_surface_tension(
    humid_air.temperature)
diff_model = diffusion_model.MixtureAveragedDiffusionModel(humid_air.gas)
diff_model.updated = False
# diff_model.update(humid_air.temperature, humid_air.pressure,
#                   humid_air.mole_fraction)
d_eff = np.asarray([diff_model.calc_diffusion_coeff(name, humid_air.temperature,
                                                    humid_air.pressure,
                                                    humid_air.mole_fraction)
                    for name in humid_air.species_names])
diff_model.update(humid_air.temperature, humid_air.pressure,
                  humid_air.mole_fraction, update_names=['O2', 'H2O'])
print(d_eff)
print(diff_model.d_eff)
print(diff_model.d_ij[:, :, 0])

