import numpy as np
import pemfc.src.fluid.fluid as fluid
import cantera as ct

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
        "temperature": 343.15,
        "pressure": 101325.0,
        "nodes": (10, 5)
    }

humid_air_pemfc = fluid.factory(fluid_dict, backend='pemfc')
humid_air_pemfc.update(,

humid_air_ct = fluid.factory(fluid_dict, backend='cantera')
humid_air_ct.update(,
# print(humid_air_pemfc.mole_fraction)
# print(humid_air_pemfc.gas.concentration)

water_ct = ct.Water()
water_ct.TP = 373.15, 101325.0


def calc_humid_composition(humidity, temperature, pressure,
                           dry_molar_composition, id_pc=-1):
    if humidity > 1.0:
        raise ValueError('relative humidity must not exceed 1.0')
    water_ct.TP = temperature, pressure
    molar_fraction_water = humidity * water_ct.P_sat / pressure
    humid_composition = np.asarray(dry_molar_composition)
    humid_composition[id_pc] = 0.0
    humid_composition /= np.sum(humid_composition, axis=0)
    humid_composition *= (1.0 - molar_fraction_water)
    humid_composition[id_pc] = molar_fraction_water
    return humid_composition


humid_air_obj_ct = ct.Solution('gri30.yaml')

species_names = fluid_dict['components'].keys()

dry_fractions = [v['molar_fraction'] for k, v
                 in fluid_dict['components'].items()]
humid_fractions = calc_humid_composition(fluid_dict['humidity'],
                                         fluid_dict['temperature'],
                                         fluid_dict['pressure'],
                                         dry_fractions, -1)

humid_air_obj_ct.X = dict(zip(species_names, humid_fractions))
all_species_names = humid_air_obj_ct.species_names
species_ids = [all_species_names.index(item) for item in species_names]
molar_composition = np.asarray([item * np.ones(fluid_dict['nodes'])
                                for item in humid_fractions])

n_iter = 1000
start_time_pemfc = time.time()
for i in range(n_iter):
    humid_air_pemfc.update(343.15, 101325,,
end_time_pemfc = time.time()

start_time_ct = time.time()
temp = 343.15
for i in range(n_iter):
    # temp += 10.0
    humid_air_ct.update(temp, 101325,,
end_time_ct = time.time()

vap_enthalpy_pemfc = humid_air_pemfc.calc_vaporization_enthalpy()
vap_enthalpy_ct = humid_air_ct.calc_vaporization_enthalpy()

print(vap_enthalpy_pemfc)
print(vap_enthalpy_ct)

print('PEMFC Time: ', end_time_pemfc-start_time_pemfc)
print('Cantera Time: ', end_time_ct-start_time_ct)

