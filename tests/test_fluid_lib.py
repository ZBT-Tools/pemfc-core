import numpy as np
import pemfc.src.fluid as fluid
import cantera as ct
# import CoolProp
# import CoolProp.CoolProp as CP

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
        "nodes": 46
    }

humid_air_pemfc = fluid.factory(fluid_dict, backend='pemfc')
humid_air_pemfc.update()

humid_air_ct = fluid.factory(fluid_dict, backend='cantera')
humid_air_ct.update()
# print(humid_air_pemfc.mole_fraction)
# print(humid_air_pemfc.gas.concentration)

water_ct = ct.Water()


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


# humid_air_obj_ct.TP = fluid_dict['temperature'], fluid_dict['pressure']
# humid_air_ct = ct.SolutionArray(humid_air_obj_ct, (fluid_dict['nodes']))
# humid_air_ct.TP = 383.15, 101325
# humid_air_ct.X[:, species_ids] = molar_composition.transpose()

print(humid_air_obj_ct())

n_iter = 1000
start_time_pemfc = time.time()
for i in range(n_iter):
    humid_air_pemfc.update(343.15, 101325)

end_time_pemfc = time.time()

start_time_ct = time.time()
temp = 343.15
for i in range(n_iter):
    # temp += 10.0
    humid_air_ct.update(temp, 101325)

    print('test')
end_time_ct = time.time()

# import numpy as np
# import itertools
# from multiprocessing import Pool
# CP.set_config_bool(CP.DONT_CHECK_PROPERTY_LIMITS, True)
#
# def generate_values(TR,P=101325):
#     """ Starting with T,R as inputs, generate all other values """
#     T,R = TR
#     psi_w = CP.HAPropsSI('psi_w','T',T,'R',R,'P',P)
#     other_output_keys = ['T_wb','T_dp','Hda','Sda','Vda','Omega']
#     outputs = {'psi_w':psi_w,'T':T,'P':P,'R':R}
#     for k in other_output_keys:
#         outputs[k] = CP.HAPropsSI(k,'T',T,'R',R,'P',P)
#     return outputs
#
# def get_supported_input_pairs():
#     """ Determine which input pairs are supported """
#     good_ones = []
#     inputs = generate_values((300, 0.5))
#     for k1, k2 in itertools.product(inputs.keys(), inputs.keys()):
#         if 'P' in [k1,k2] or k1==k2:
#             continue
#         args = ('psi_w', k1, inputs[k1], k2, inputs[k2], 'P', inputs['P'])
#         try:
#             psi_w_new = CP.HAPropsSI(*args)
#             if not np.isfinite(psi_w_new):
#                 raise ValueError('Returned NaN; not ok')
#             good_ones.append((k1,k2))
#         except BaseException as BE:
#             pass
#             if 'currently at least one of' in str(BE) or 'cannot provide two inputs' in str(BE):
#                 pass
#             else:
#                 print(BE)
#                 good_ones.append((k1,k2))
#     return good_ones
# supported_pairs = get_supported_input_pairs()
#
# def calculate(inputs):
#     """ For a given input, try all possible input pairs """
#     errors = []
#     for k1, k2 in supported_pairs:
#         psi_w_input = inputs['psi_w']
#         args = 'psi_w',k1,inputs[k1],k2,inputs[k2],'P',inputs['P']
#         try:
#             psi_w_new = CP.HAPropsSI(*args)
#             if not np.isfinite(psi_w_new):
#                 raise ValueError('Returned NaN; not ok')
#         except BaseException as BE:
#             errors.append((str(BE),args, inputs))
#     return errors
#
# if __name__ == '__main__':
#     print(CoolProp.__version__)
#     TR = itertools.product(np.linspace(240, 345, 11), np.linspace(0, 1, 11))
#     with Pool(processes=2) as pool:
#         input_values = pool.map(generate_values, TR)
#         errors = pool.map(calculate, input_values)
#         for err in itertools.chain.from_iterable(errors):
#             print(err)

print('PEMFC Time: ', end_time_pemfc-start_time_pemfc)
print('Cantera Time: ', end_time_ct-start_time_ct)