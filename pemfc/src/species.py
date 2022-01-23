# general imports
import numpy as np
from numpy.polynomial.polynomial import polyval
from abc import ABC, abstractmethod

# local module imports
from data import material_properties as mat_prop
from . import global_functions as g_func


class FluidProperties(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def calc_property(self, property_name, temperature, pressure=101325.0):
        pass


class ConstantProperties(FluidProperties):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']
    ATTRIBUTE_NAMES = [name.replace(' ', '_').lower()
                       for name in PROPERTY_NAMES]

    def __init__(self, name, **kwargs):
        super().__init__()
        self.name = name
        self.property = {}
        const_coeffs = mat_prop.constant
        for i, attr_name in enumerate(self.ATTRIBUTE_NAMES):
            prop = kwargs.get(attr_name, None)
            if prop is None:
                prop = kwargs.get(self.PROPERTY_NAMES[i], None)
            if self.name in const_coeffs[self.PROPERTY_NAMES[i]] \
                    and prop is None:
                prop = const_coeffs[self.PROPERTY_NAMES[i]][self.name]
                setattr(self, attr_name, prop)
            elif prop is not None:
                setattr(self, attr_name, prop)
            else:
                raise ValueError('Value for {} must be in database or '
                                 'directly provided'.format(attr_name))
            self.property[self.PROPERTY_NAMES[i]] = prop

    def calc_property(self, property_name, temperature, pressure=101325.0):
        return self.property[property_name]


class PolynomialProperties(FluidProperties, ABC):

    def __init__(self, species_list, property_names, poly_coeffs):
        super().__init__()
        self.names = g_func.ensure_list(species_list)
        self.coeff_dict_dict = dict()
        self.property_names = property_names
        for prop_name in self.property_names:
            self.coeff_dict_dict[prop_name] = dict()
            for species_name in self.names:
                self.coeff_dict_dict[prop_name][species_name] = \
                    poly_coeffs[prop_name][species_name]
        self.coeff_dict_arr = dict()
        for prop_name in self.property_names:
            self.coeff_dict_arr[prop_name] = \
                np.stack([np.flip(self.coeff_dict_dict[prop_name][item],
                                  axis=-1)
                          for item in self.coeff_dict_dict[prop_name]], axis=-1)

    def calc_property(self, property_name, temperature, pressure=101325.0):
        if property_name in self.property_names:
            return polyval(temperature, self.coeff_dict_arr[property_name])
        else:
            raise ValueError('property_name {} not valid'.format(property_name))


class IncompressibleProperties(PolynomialProperties):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']

    def __init__(self, species_list):
        poly_coeffs = mat_prop.incompressible_polynomials
        super().__init__(species_list, self.PROPERTY_NAMES, poly_coeffs)
        self.name = self.names[0]

    def calc_specific_heat(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Specific Heat'])

    def calc_viscosity(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Viscosity'])

    def calc_thermal_conductivity(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Thermal Conductivity'])

    def calc_density(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Density'])


class GasProperties(PolynomialProperties):

    PROPERTY_NAMES = ['Specific Heat', 'Viscosity', 'Thermal Conductivity']

    def __init__(self, species_list):
        poly_coeffs = mat_prop.gas_polynomials
        super().__init__(species_list, self.PROPERTY_NAMES, poly_coeffs)

        mole_weights = mat_prop.molecular_weight
        self.coeff_dict_dict2 = dict()
        self.mw = list()
        for species_name in self.names:
            if species_name not in mole_weights:
                raise AttributeError('No data available for '
                                     'provided specie: ' + species_name)
            else:
                self.mw.append(mole_weights[species_name])
                self.coeff_dict_dict2[species_name] = dict()
                for prop_name in self.PROPERTY_NAMES:
                    self.coeff_dict_dict2[species_name][prop_name] = \
                        poly_coeffs[prop_name][species_name]
        self.mw = np.asarray(self.mw)

    def calc_specific_heat(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Specific Heat'])

    def calc_viscosity(self, temperature):
        return polyval(temperature, self.coeff_dict_arr['Viscosity'])

    def calc_thermal_conductivity(self, temperature, pressure):
        lambda_1_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['Thermal Conductivity'][:][0])
        lambda_10_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['Thermal Conductivity'][:][1])
        result = lambda_1_bar \
            + (pressure - 1.e5) / 9.e5 * (lambda_10_bar - lambda_1_bar)
        result *= 10.0
        return result

    def calc_property(self, property_name, temperature, pressure=101325.0):
        if property_name == 'Thermal Conductivity':
            return self.calc_thermal_conductivity(temperature, pressure)
        else:
            return super().calc_property(property_name, temperature, pressure)


class PhaseChangeProperties(PolynomialProperties):

    PROPERTY_NAMES = ['Saturation Pressure', 'Vaporization Enthalpy']

    def __init__(self, liquids_dict):
        # print("Constructor of Two Phase Species")
        poly_coeffs = mat_prop.phase_change_polynomials

        if not isinstance(liquids_dict, dict):
            raise TypeError('Input data must be provided as dictionary '
                            'with species names as keys and objects'
                            'type FluidProperties as values')
        self.names = list(liquids_dict.keys())

        for name in self.names:
            if name not in poly_coeffs[self.PROPERTY_NAMES[0]]:
                raise NotImplementedError('No phase change data available for '
                                          'provided liquid specie: ' + name)

        super().__init__(self.names, self.PROPERTY_NAMES, poly_coeffs)

        # self.gas = GasProperties(self.names)
        # if len(liquids_dict) == 1:
        #     self.liquid = next(iter(liquids_dict.values()))
        # else:
        #     density = [value.density for key, value in liquids_dict]
        #     specific_heat = \
        #       [value.specific_heat for key, value in liquids_dict]
        #     viscosity = [value.viscosity for key, value in liquids_dict]
        #     thermal_conductivity = \
        #         [value.thermal_conductivity for key, value in liquids_dict]
        #     self.liquid = \
        #         FluidProperties('liquids', density=np.asarray(density),
        #                         viscosity=np.asarray(viscosity),
        #                         specific_heat=np.asarray(specific_heat),
        #                         thermal_conductivity=
        #                         np.asarray(thermal_conductivity))

    def calc_saturation_pressure(self, temperature):
        return polyval(temperature,
                       self.coeff_dict_arr['Saturation Pressure'])

    def calc_vaporization_enthalpy(self, temperature):
        return polyval(temperature,
                       self.coeff_dict_arr['Vaporization Enthalpy'])

    def calc_property(self, property_name, temperature, **kwargs):
        if property_name == 'Saturation Pressure':
            return self.calc_saturation_pressure(temperature)
        elif property_name == 'Vaporization Enthalpy':
            return self.calc_vaporization_enthalpy(temperature)
        else:
            raise ValueError('property_name {} not valid'.format(
                             property_name))
