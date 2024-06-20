# general imports
import sys
import numpy as np
from numpy.polynomial.polynomial import polyval
from abc import ABC, abstractmethod

# local module imports
if 'main_app.py' in sys.argv[0]:
    from data import material_properties as mat_prop
else:
    from pemfc.data import material_properties as mat_prop
from .. import global_functions as g_func


class FluidProperties(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def calc_property(self, property_name, temperature, pressure=101325.0):
        pass


class ConstantProperties(FluidProperties):

    PROPERTY_NAMES = ['density', 'specific_heat', 'viscosity',
                      'thermal_conductivity']
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

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      tensor=False):
        if property_name in self.property_names:
            return polyval(temperature, self.coeff_dict_arr[property_name],
                           tensor=tensor)
        else:
            raise ValueError('property_name {} not valid'.format(property_name))


class IncompressibleProperties(PolynomialProperties):

    PROPERTY_NAMES = ['density', 'specific_heat', 'viscosity',
                      'thermal_conductivity']

    def __init__(self, species_list):
        poly_coeffs = mat_prop.incompressible_polynomials
        super().__init__(species_list, self.PROPERTY_NAMES, poly_coeffs)
        self.name = self.names[0]

    def calc_specific_heat(self, temperature, tensor=False):
        return polyval(temperature, self.coeff_dict_arr['specific_heat'],
                       tensor=tensor)

    def calc_viscosity(self, temperature, tensor=False):
        return polyval(temperature, self.coeff_dict_arr['viscosity'],
                       tensor=tensor)

    def calc_thermal_conductivity(self, temperature, tensor=False):
        return polyval(temperature, self.coeff_dict_arr[
            'thermal_conductivity'], tensor=tensor)

    def calc_density(self, temperature, tensor=False):
        return polyval(temperature, self.coeff_dict_arr['density'],
                       tensor=tensor)


class GasProperties(PolynomialProperties):

    PROPERTY_NAMES = ['specific_heat', 'viscosity', 'thermal_conductivity']

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

    def calc_specific_heat(self, temperature, tensor=True):
        """
        Specific heat capacity at constant pressure in J/(kg-K)
        """
        return polyval(temperature, self.coeff_dict_arr['specific_heat'],
                       tensor=tensor)

    def calc_viscosity(self, temperature, tensor=True):
        """
        Dynamic viscosity in Pa-s
        """
        return polyval(temperature, self.coeff_dict_arr['viscosity'],
                       tensor=tensor)

    def calc_thermal_conductivity(self, temperature, pressure, tensor=True):
        """
        Thermal conductivity in W/(m-K)
        """
        lambda_1_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['thermal_conductivity'][:][0],
                    tensor=tensor)
        lambda_10_bar = \
            polyval(temperature,
                    self.coeff_dict_arr['thermal_conductivity'][:][1],
                    tensor=tensor)
        result = lambda_1_bar \
            + (pressure - 1.e5) / 9.e5 * (lambda_10_bar - lambda_1_bar)
        result *= 10.0
        return result

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      tensor=True):
        if property_name == 'thermal_conductivity':
            return self.calc_thermal_conductivity(temperature, pressure,
                                                  tensor=tensor)
        else:
            return super().calc_property(property_name, temperature,
                                         pressure, tensor=tensor)


class PhaseChangeProperties(PolynomialProperties):

    PROPERTY_NAMES = ['saturation_pressure', 'vaporization_enthalpy',
                      'surface_tension']

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

    def calc_saturation_pressure(self, temperature, tensor=False):
        """
        Saturation pressure in Pa
        """
        return polyval(temperature,
                       self.coeff_dict_arr['saturation_pressure'],
                       tensor=tensor)

    def calc_vaporization_enthalpy(self, temperature, tensor=False):
        """
        Vaporization enthalpy in J/mol
        """
        return polyval(temperature,
                       self.coeff_dict_arr['vaporization_enthalpy'],
                       tensor=tensor)

    def calc_surface_tension(self, temperature, tensor=False):
        """
        Surface tension in N/m
        """
        return polyval(temperature,
                       self.coeff_dict_arr['surface_tension'], tensor=tensor)

    def calc_humid_composition(self, humidity, temperature, pressure,
                               dry_molar_composition, id_pc):
        if humidity > 1.0:
            raise ValueError('relative humidity must not exceed 1.0')
        molar_fraction_phase_change_species = \
            humidity * self.calc_saturation_pressure(temperature) / pressure
        humid_composition = np.asarray(dry_molar_composition)
        humid_composition[id_pc] = 0.0
        humid_composition /= np.sum(humid_composition, axis=0)
        humid_composition *= (1.0 - molar_fraction_phase_change_species)
        humid_composition[id_pc] = molar_fraction_phase_change_species
        return humid_composition

    def calc_property(self, property_name, temperature, **kwargs):
        if property_name == 'saturation_pressure':
            return self.calc_saturation_pressure(temperature)
        elif property_name == 'vaporization_enthalpy':
            return self.calc_vaporization_enthalpy(temperature)
        else:
            raise ValueError('property_name {} not valid'.format(
                             property_name))
