import numpy as np


p_saturation_param = \
    (-4.66691122e-18, 2.19146750e-14, -4.56208833e-11, 5.54957241e-08,
     -4.37186346e-05, 2.33207549e-02, -8.53414571e+00, 2.11600925e+03,
     -3.40245294e+05, 3.20415279e+07, -1.34211567e+09)

vaporization_enthalpy_param =\
    (-2.01107201e-18, 8.92669752e-15, -1.76751771e-11,
     2.05547260e-08, -1.55445645e-05,  7.98692642e-03,
     -2.82333561e+00,  6.77951176e+02, -1.05826022e+05,
     9.69666280e+06, -3.95999179e+08)


class TwoPhaseSpecies:

    def __init__(self, p_sat_param, h_vap_param):
        self.p_sat_param = p_sat_param
        self.h_vap_param = h_vap_param

    def calc_p_sat(self, temperature):
        return np.polyval(self.p_sat_param, temperature)

    def calc_h_vap(self, temperature):
        return np.polyval(self.h_vap_param, temperature)

    def calc_humid_composition(self, humidity, temperature, pressure,
                               dry_molar_composition, id_pc):
        if humidity > 1.0:
            raise ValueError('relative humidity must not exceed 1.0')
        molar_fraction_phase_change_species = \
            humidity * self.calc_p_sat(temperature) / pressure
        humid_composition = np.asarray(dry_molar_composition)
        humid_composition[id_pc] = 0.0
        humid_composition /= np.sum(humid_composition, axis=0)
        humid_composition *= (1.0 - molar_fraction_phase_change_species)
        humid_composition[id_pc] = molar_fraction_phase_change_species
        return humid_composition


water = TwoPhaseSpecies(p_saturation_param, vaporization_enthalpy_param)
