# general imports
import sys
import numpy as np
from abc import ABC, abstractmethod

# local module imports
from .output_object import OutputObject
from . import constants, global_functions as g_func, species
from collections import OrderedDict

if 'main_app.py' in sys.argv[0]:
    from data import water_properties as water_props
else:
    from pemfc.data import water_properties as water_props


class OneDimensionalFluid(ABC, OutputObject):

    PROPERTY_NAMES = ['Density', 'Specific Heat', 'Viscosity',
                      'Thermal Conductivity']
    TYPE_NAME = 'Base Fluid'

    def __init__(self, nx, name, temperature=298.15, pressure=101325.0,
                 **kwargs):
        # print("__init__ for Fluid")
        super().__init__(name)
        self.nodes = nx
        self.print_variables = \
            {
                'names': ['temperature', 'pressure'],
                'units': ['K', 'Pa'],
                # 'sub_names': ['None', 'None']
            }
        self.combine_print_variables(self.print_variables,
                                     kwargs.get('print_variables', None))
        temperature = np.asarray(temperature)
        pressure = np.asarray(pressure)
        self._temperature = np.zeros(self.nodes)
        self._pressure = np.zeros(self.nodes)
        self.write_input_to_array(temperature, self._temperature)
        self.write_input_to_array(pressure, self._pressure)

        self.property = dict()
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(self.nodes)

    @staticmethod
    def write_input_to_array(value, array):
        if np.ndim(value) == 0 \
                or value.shape == array.shape:
            array[:] = value
        else:
            raise ValueError('argument value must be scalar or '
                             'have the same shape as 1D array')

    @abstractmethod
    def update(self, temperature, pressure, *args, **kwargs):
        if np.any(temperature < 200.0):
            raise ValueError('temperature too low, check boundary conditions')
        if np.any(temperature > 1000.0):
            raise ValueError('temperature too high, check boundary conditions')
        self.temperature = temperature
        self.pressure = pressure

    @abstractmethod
    def calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass

    def calc_fraction(self, composition, axis=0, recursion_count=0):
        """
        Calculates mixture fractions based on a multi-dimensional
        array with different species along the provided axis.
        """
        np.seterr(all='raise')
        try:
            comp_sum = np.sum(composition, axis)
            # result = np.divide(composition, comp_sum, out=composition,
            #                    where=comp_sum != 0.0)
            result = composition / comp_sum
            return result
        except FloatingPointError:
            if axis == 0:
                composition = g_func.fill_zero_sum(composition, axis=-1)
            else:
                composition = g_func.fill_zero_sum(composition, axis=0)
            recursion_count += 1
            if recursion_count > 1:
                raise ValueError('Something wrong in '
                                 '{}.calc_fraction()'.format(self))
            return self.calc_fraction(composition, axis, recursion_count)

    @property
    def temperature(self):
        return self._temperature

    @property
    def pressure(self):
        return self._pressure

    @temperature.setter
    def temperature(self, value):
        self.write_input_to_array(value, self._temperature)

    @pressure.setter
    def pressure(self, value):
        self.write_input_to_array(value, self._pressure)

    @property
    def density(self):
        return self.property['Density']

    @property
    def viscosity(self):
        return self.property['Viscosity']

    @property
    def thermal_conductivity(self):
        return self.property['Thermal Conductivity']

    @property
    def specific_heat(self):
        return self.property['Specific Heat']

    # @density.setter
    # def density(self, value):
    #     self.property['Density'] = value
    #
    # @viscosity.setter
    # def viscosity(self, value):
    #     self.property['Viscosity'] = value
    #
    # @thermal_conductivity.setter
    # def thermal_conductivity(self, value):
    #     self.property['Thermal Conductivity'] = value
    #
    # @specific_heat.setter
    # def specific_heat(self, value):
    #     self.property['Specific Heat'] = value

    def rescale(self, new_nx):
        """
        linearly rescales all numpy arrays with the 1D discretization as first
        dimension
        :param new_nx: new discretization along first dimension
        :return: only modification of attributes
        """
        if new_nx != self.nodes:
            attr_list = [a for a in dir(self) if not a.startswith('__')]
            for name in attr_list:
                attr = getattr(self, name)
                if name == 'property':
                    for key in attr.keys():
                        rescaled = self.rescale_attribute(attr[key], new_nx)
                        if rescaled is not None:
                            attr[key] = rescaled
                else:
                    type_attr = getattr(type(self), name, None)
                    if not isinstance(type_attr, property):
                        rescaled = self.rescale_attribute(attr, new_nx)

                        if rescaled is not None:
                            setattr(self, name, rescaled)
        self.add_print_variables(self.print_variables)

    def rescale_attribute(self, attribute, new_nx):
        if isinstance(attribute, np.ndarray):
            return self.rescale_array(attribute, new_nx)
        else:
            return None

    def rescale_array(self, array, new_nx):
        try:
            return self.linear_rescale_1d(array, new_nx)
        except ValueError:
            return self.linear_rescale_2d(array, new_nx)

    def linear_rescale_1d(self, array, new_nodes):
        if new_nodes != self.nodes:
            if len(array.shape) == 1:
                first = array[0]
                last = array[-1]
                return np.linspace(first, last, new_nodes)
            else:
                raise ValueError('argument array must be one-dimensional')
        else:
            raise TypeError('argument array must be of type numpy.ndarray')

    def linear_rescale_2d(self, array, new_nodes):
        if new_nodes != self.nodes:
            if isinstance(array, np.ndarray):
                if len(array.shape) == 2:
                    first = array[0, :]
                    last = array[-1, :]
                    linear_rescaled = \
                        np.array([np.linspace(first[i], last[i], new_nodes)
                                  for i in range(len(first))])
                    return linear_rescaled.transpose()
                else:
                    raise ValueError('argument array must be one-dimensional')
            else:
                raise TypeError('argument array must be of type numpy.ndarray')


class ConstantFluid(OneDimensionalFluid):

    TYPE_NAME = 'Constant Fluid'

    def __init__(self, nx, name, fluid_props, temperature=298.15,
                 pressure=101325.0, **kwargs):
        # print("__init__ for IncompressibleFluid")
        super().__init__(nx, name, temperature, pressure, **kwargs)
        self.name = name
        if not isinstance(fluid_props, species.ConstantProperties):
            raise TypeError('Argument fluid_props must be of type '
                            'ConstantProperties')
        else:
            self.property['Density'][:] = fluid_props.density
            self.property['Specific Heat'][:] = fluid_props.specific_heat
            self.property['Viscosity'][:] = fluid_props.viscosity
            self.property['Thermal Conductivity'][:] = \
                fluid_props.thermal_conductivity
        self.add_print_variables(self.print_variables)

    def update(self, temperature, pressure, *args, **kwargs):
        super().update(temperature, pressure)

    def calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass


class IncompressibleFluid(OneDimensionalFluid):

    TYPE_NAME = 'Incompressible Fluid'

    def __init__(self, nx, name, fluid_props, temperature=298.15,
                 pressure=101325.0, **kwargs):
        # print("__init__ for IncompressibleFluid")
        if isinstance(fluid_props, species.IncompressibleProperties):
            self.properties = fluid_props
        else:
            raise TypeError('Argument fluid_props must be of type '
                            'IncompressibleProperties')
        super().__init__(nx, name, temperature, pressure, **kwargs)
        self.add_print_variables(self.print_variables)

    def update(self, temperature, pressure=101325.0, *args, **kwargs):
        super().update(temperature, pressure)
        self.calc_properties(self._temperature, self._pressure)

    def calc_property(self, property_name, temperature):
        """
        Wrapper function for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :return: the calculated 1D array of the specific property
        """
        return self.properties.calc_property(property_name, temperature)

    def calc_properties(self, temperature, pressure=101325.0, **kwargs):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: float or 1D pressure array
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self.calc_property(prop, temperature)


class GasMixture(OneDimensionalFluid):

    TYPE_NAME = 'Gas Mixture'

    def __init__(self, nx, name, species_dict, mole_fractions,
                 temperature=298.15, pressure=101325.0, **kwargs):
        # print("__init__ for Gas Mixture")
        print_variables = \
            {
                'names': ['mole_fraction'],
                'units': ['-'],
                'sub_names': ['self.species.names']
            }
        super().__init__(nx, name, temperature, pressure,
                         print_variables=print_variables, **kwargs)
        self.name = name
        species_names = list(species_dict.keys())
        # self.n_species = len(species_names)
        self.gas_constant = constants.GAS_CONSTANT
        self.species = species.GasProperties(species_names)
        self.n_species = len(self.species.names)
        self.species_viscosity = \
            self.species.calc_viscosity(self._temperature).transpose()

        mole_fractions = g_func.ensure_list(mole_fractions)
        mole_fractions = np.asarray(mole_fractions)
        if len(mole_fractions) != self.n_species:
            raise ValueError('Initial mole fractions must be provided '
                             'for all species')
        mole_fraction_sum = np.sum(mole_fractions)
        if mole_fraction_sum != 1.0:
            if np.round(mole_fraction_sum, 3) != 1.0:
                raise ValueError('Initial mole fractions must add up to unity')
            else:
                mole_fractions[-1] = 1.0 - np.sum(mole_fractions[:-1])

        self.array_shape_2d = (nx, self.n_species)
        self._mole_fraction = \
            np.ones(self.array_shape_2d) * mole_fractions
        self.mw = np.sum(self._mole_fraction * self.species.mw, axis=-1)
        mass_fraction_init = \
            mole_fractions * self.species.mw / \
            np.sum(mole_fractions * self.species.mw)
        self._mass_fraction = \
            np.ones(self.array_shape_2d) * mass_fraction_init
        self._concentration = np.zeros(self.array_shape_2d)
        if isinstance(type(self), GasMixture):
            self.calc_properties(self._temperature, self._pressure,
                                 method='ideal')
            self._concentration[:] = self.calc_concentration().transpose()

        # self.add_print_data(self.mole_fraction, 'Mole Fraction',
        #                     sub_names=self.species.names)
        self.add_print_variables(self.print_variables)

    @property
    def mole_fraction(self):
        return self._mole_fraction.transpose()

    @property
    def mass_fraction(self):
        return self._mass_fraction.transpose()

    @property
    def concentration(self):
        return self._concentration.transpose()

    @concentration.setter
    def concentration(self, value):
        self._concentration[:] = value.transpose()

    def update(self, temperature, pressure, mole_composition=None,
               method='ideal', *args, **kwargs):
        super().update(temperature, pressure)
        if mole_composition is None:
            mole_composition = self.mole_fraction
        elif np.shape(mole_composition)[0] != self.n_species:
            raise ValueError('First dimension of composition must be equal to '
                             'number of species')
        if np.max(mole_composition) > constants.SMALL:
            mole_composition[mole_composition < constants.SMALL] = 0.0
        self._mole_fraction[:] = self.calc_mole_fraction(mole_composition)
        self.calc_molar_mass()
        self._mass_fraction[:] = \
            self.calc_mass_fraction(self._mole_fraction)
        self.calc_properties(self._temperature, self._pressure, method)
        self._concentration[:] = \
            self.calc_concentration().transpose()

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self._pressure \
            / (self.gas_constant * self._temperature)
        return self._mole_fraction.transpose() * total_mol_conc

    def calc_mole_fraction(self, mole_composition):
        if np.min(mole_composition) < 0.0:
            raise ValueError('mole_composition must not be smaller zero')
        mole_fraction = self.calc_fraction(mole_composition).transpose()
        return mole_fraction

    def calc_molar_mass(self, mole_fraction=None):
        if mole_fraction is None:
            self.mw[:] = np.sum(self._mole_fraction * self.species.mw, axis=-1)
        else:
            return np.sum(mole_fraction.transpose() * self.species.mw, axis=-1)

    def calc_mass_fraction(self, mole_fraction, molar_mass=None):
        if molar_mass is None:
            # mw = np.where(self.mw == 0.0, 1.0, self.mw)
            return np.outer(1.0 / self.mw, self.species.mw) * mole_fraction
        else:
            # molar_mass = np.where(molar_mass == 0.0, 1.0, molar_mass)
            return np.outer(1.0 / molar_mass, self.species.mw) * mole_fraction

    def calc_specific_heat(self, temperature, mass_fraction=None):
        species_cp = self.species.calc_specific_heat(temperature).transpose()
        if mass_fraction is None:
            return np.sum(self._mass_fraction * species_cp, axis=-1)
        else:
            return np.sum(mass_fraction.transpose() * species_cp, axis=-1)

    def calc_viscosity(self, temperature):
        """
        Calculates the mixture viscosity of a
        gas according to Herning and Zipperer.
        """
        self.species_viscosity[:] = \
            self.species.calc_viscosity(temperature).transpose()
        x_sqrt_mw = self._mole_fraction * np.sqrt(self.species.mw)
        x_sqrt_mw = np.where(x_sqrt_mw == 0.0, 1.0, x_sqrt_mw)
        return np.sum(self.species_viscosity * x_sqrt_mw, axis=-1) \
            / np.sum(x_sqrt_mw, axis=-1)

    def calc_wilke_coefficients(self):
        """
        Calculates the wilke coefficients for
        each species combination of a gas.
        """
        visc = self.species_viscosity.transpose()
        mw = self.species.mw
        alpha = []
        for i in range(self.n_species):
            beta = []
            for j in range(self.n_species):
                a = np.power((1. + np.power((visc[i] / visc[j]), 0.5)
                             * np.power(mw[j] / mw[i], 0.25)),  2.0)
                b = np.power(np.sqrt(8.) * (1. + mw[j] / mw[i]), -0.5)
                beta.append(a/b)
            alpha.append(beta)
        return np.asarray(alpha)

    def calc_thermal_conductivity(self, temperature, pressure):
        """
        Calculates the heat conductivity of a gas mixture,
        according to Wilkes equation.
        """
        self._mole_fraction[1:] = np.maximum(1e-16, self._mole_fraction[1:])
        wilke_coeffs = self.calc_wilke_coefficients()
        lambda_species = \
            self.species.calc_thermal_conductivity(temperature, pressure)
        lambda_mix = np.zeros(len(temperature))
        for i in range(self.n_species):
            a = self._mole_fraction[:, i] * lambda_species[i]
            b = np.sum(self._mole_fraction.transpose() * wilke_coeffs[i],
                       axis=0)
            b += 1e-16
            lambda_mix += a / b
        return lambda_mix

    def calc_density(self, temperature, pressure, method="ideal"):
        """
        Calculate gas mixture density
        :param temperature: temperature array
        :param pressure: pressure array
        :param method: string indicating the calculation method
        :return: density of the mixture
        """
        if method == "ideal":
            return pressure * self.mw / (self.gas_constant * temperature)
        else:
            raise NotImplementedError('Method {} to calculate '
                                      'gas mixture density has not '
                                      'been implemented'.format(method))

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      method='ideal'):
        """
        Wrapper function for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        if property_name == 'Specific Heat':
            return self.calc_specific_heat(temperature)
        elif property_name == 'Viscosity':
            return self.calc_viscosity(temperature)
        elif property_name == 'Thermal Conductivity':
            return self.calc_thermal_conductivity(temperature, pressure)
        elif property_name == 'Density':
            return self.calc_density(temperature, pressure, method)
        else:
            raise ValueError('property_name '
                             '{} not valid'.format(property_name))

    def calc_properties(self, temperature, pressure=101325.0, method='ideal'):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self.calc_property(prop, temperature, pressure, method)


class TwoPhaseMixture(OneDimensionalFluid):

    TYPE_NAME = 'Two-Phase Mixture'

    def __init__(self, nx, name, species_dict, mole_fractions,
                 liquid_props=None, temperature=298.15, pressure=101325.0,
                 **kwargs):
        # print("__init__ for TwoPhaseMixture")
        print_variables = \
            {
                'names': ['humidity'],
                'units': ['-'],
                'sub_names': ['None']
            }
        super().__init__(nx, name, temperature, pressure,
                         print_variables=print_variables, **kwargs)
        if not isinstance(species_dict, dict):
            raise TypeError('Argument species_dict must be a dict '
                            'containing all species names and their '
                            'expected aggregate states in terms of "gas", '
                            '"gas-liquid", or "liquid"')
        gas_species_dict = {k: v for k, v in species_dict.items() if 'gas' in v}
        phase_change_species_names = \
            [key for key in species_dict
             if 'gas' and 'liquid' in species_dict[key]]
        if liquid_props is None:
            liquid_props = \
                species.IncompressibleProperties(phase_change_species_names[0])
            self.phase_change_species = \
                species.PhaseChangeProperties({liquid_props.name: liquid_props})
        elif isinstance(liquid_props, dict):
            self.phase_change_species = \
                species.PhaseChangeProperties(liquid_props)
        else:
            raise TypeError('Data for PhaseChangeSpecies object '
                            'can only be provided as dictionary with species '
                            'name as key and FluidProperties object as value '
                            'for the liquid properties')
        self.liquid = IncompressibleFluid(nx, name + ': Liquid Phase',
                                          fluid_props=liquid_props,
                                          temperature=self._temperature,
                                          pressure=self._pressure)
        self.gas = GasMixture(nx, name + ': Gas Phase',
                              species_dict=gas_species_dict,
                              mole_fractions=mole_fractions,
                              temperature=self._temperature,
                              pressure=self._pressure)
        ids_pc = [self.species.names.index(name) for name in
                  phase_change_species_names]
        if len(ids_pc) > 1:
            raise NotImplementedError('at the moment only one species '
                                      'undergoing phase change is allowed')
        self.id_pc = ids_pc[0]
        all_ids = np.array(list(range(len(self.species.names))), dtype='int32')
        ids_no_pc = np.delete(all_ids, ids_pc)
        self.ids_no_pc = list(ids_no_pc)

        # update mole_fractions if humidity is provided
        humidity = kwargs.get('humidity', None)
        if humidity is not None:
            mole_fractions = \
                water_props.water.calc_humid_composition(humidity, temperature,
                                                         pressure,
                                                         mole_fractions,
                                                         self.id_pc)

        # Total properties (both phases)
        self.n_species = len(species_dict)
        self.array_shape_2d = self.gas.array_shape_2d
        self._mole_fraction = \
            np.ones(self.array_shape_2d) * mole_fractions
        mass_fraction_init = \
            mole_fractions * self.gas.species.mw / \
            np.sum(mole_fractions * self.gas.species.mw)
        self._mass_fraction = \
            np.ones(self.array_shape_2d) * mass_fraction_init
        self.mw = np.zeros(nx)

        self.liquid_mass_fraction = np.zeros(nx)
        self.liquid_mole_fraction = np.zeros(nx)
        self.humidity = np.zeros(nx)
        self.saturation_pressure = np.zeros(nx)

        # Print data
        self.add_print_variables(self.print_variables)
        # self.add_print_data(self.humidity, 'Humidity')
        # print(self.print_data)
        # print(self.gas.print_data)

    def rescale(self, new_nx):
        self.gas.rescale(new_nx)
        self.liquid.rescale(new_nx)
        super().rescale(new_nx)
    # @property
    # def mole_fraction_gas(self):
    #     return self._mole_fraction_gas.transpose()
    #
    # @property
    # def mass_fraction_gas(self):
    #     return self._mass_fraction_gas.transpose()

    @property
    def temperature(self):
        return self._temperature

    @property
    def pressure(self):
        return self._pressure

    @temperature.setter
    def temperature(self, value):
        self.write_input_to_array(value, self._temperature)
        self.gas._temperature = self._temperature
        self.liquid._temperature = self._temperature

    @pressure.setter
    def pressure(self, value):
        self.write_input_to_array(value, self._pressure)
        self.gas._pressure = self._pressure
        self.liquid._pressure = self._pressure

    @property
    def mole_fraction(self):
        return self._mole_fraction.transpose()

    @property
    def mass_fraction(self):
        return self._mass_fraction.transpose()

    @property
    def species(self):
        return self.gas.species

    def _set_name(self, name):
        super()._set_name(name)
        self.gas.name = self.name + ': Gas Phase'
        self.liquid.name = self.name + ': Liquid Phase'

    def update(self, temperature, pressure, mole_composition=None,
               gas_mole_composition=None, method='ideal', *args, **kwargs):
        super().update(temperature, pressure)
        if mole_composition is not None:
            if np.max(mole_composition) > constants.SMALL:
                mole_composition[mole_composition < constants.SMALL] = 0.0
                self._mole_fraction[:] = \
                    self.gas.calc_mole_fraction(mole_composition)

        self.mw[:] = self.gas.calc_molar_mass(self.mole_fraction)

        self._mass_fraction[:] = \
            self.gas.calc_mass_fraction(self._mole_fraction, self.mw)
        self.saturation_pressure[:] = \
            self.phase_change_species.calc_saturation_pressure(temperature)
        if gas_mole_composition is None \
                or np.max(gas_mole_composition) < constants.SMALL:
            gas_mole_composition = self.calc_concentration()
        self.gas.update(temperature, pressure, gas_mole_composition, method)
        self.liquid.update(temperature, pressure)
        self.liquid_mole_fraction[:] = \
            np.sum(mole_composition, axis=0) \
            - np.sum(gas_mole_composition, axis=0)
        self.liquid_mole_fraction[self.liquid_mole_fraction < 0.0] = 0.0
        self.liquid_mass_fraction = \
            np.sum(mole_composition * self.mw, axis=0) \
            - np.sum(gas_mole_composition * self.gas.mw, axis=0)
        self.liquid_mass_fraction[self.liquid_mass_fraction < 0.0] = 0.0

        # dry_conc = np.copy(gas_conc)
        # dry_conc[self.id_pc] = 0.0
        # total_conc = dry_conc
        # total_conc[self.id_pc] = np.sum(dry_conc, axis=0) \
        #     * self.mole_fraction[self.id_pc] \
        #     / (1.0 - self.mole_fraction[self.id_pc])

        # try:
        #     self.liquid_mole_fraction[:] = \
        #         1.0 - np.sum(gas_conc, axis=0) / np.sum(total_conc, axis=0)
        # except FloatingPointError:
        #     raise FloatingPointError('error in calculation of total '
        #                              'concentration, should not be zero')
        #
        # self.liquid_mass_fraction[:] = \
        #     1.0 - np.sum(gas_conc * self.gas.mw, axis=0) \
        #     / np.sum(total_conc * self.mw, axis=0)

        self.calc_properties(temperature, pressure, method)
        self.calc_humidity()

    def calc_property(self, property_name, temperature, pressure=101325.0,
                      method='ideal'):
        """
        Calculate a single property listed in PROPERTY_NAMES. Wrapper function
        for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        if property_name == 'Specific Heat':
            return self.calc_specific_heat()
        else:
            return self.gas.property[property_name]

    def calc_properties(self, temperature, pressure=101325.0, method='ideal'):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self.calc_property(prop, temperature, pressure, method)

    def calc_specific_heat(self):
        return self.liquid_mass_fraction \
            * self.liquid.specific_heat \
            + (1.0 - self.liquid_mass_fraction) * self.gas.specific_heat

    def calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        :return: gas phase concentration
        (n_species x n_nodes)
        """
        r_t = self.gas.gas_constant * self._temperature
        total_gas_conc = self._pressure / r_t
        conc = self.mole_fraction * total_gas_conc
        # all_conc = np.copy(conc)
        sat_conc = self.saturation_pressure / r_t
        dry_mole_fraction = np.copy(self.mole_fraction)
        dry_mole_fraction[self.id_pc] = 0.0
        dry_mole_fraction = self.calc_fraction(dry_mole_fraction)
        for i in range(self.n_species):
            if i == self.id_pc:
                conc[self.id_pc] = np.where(conc[self.id_pc] > sat_conc,
                                            sat_conc, conc[self.id_pc])
            else:
                try:
                    conc[i] = \
                        np.where(conc[self.id_pc] > sat_conc,
                                 (total_gas_conc - sat_conc)
                                 * dry_mole_fraction[i],
                                 conc[i])
                except FloatingPointError:
                    raise FloatingPointError
        return np.maximum(conc, 0.0)

    # def calc_concentration(self):
    #     """
    #     Calculates the gas phase molar concentrations.
    #     :return: gas phase concentration and total concentration arrays
    #     (n_species x n_nodes)
    #     """
    #     r_t = self.gas.gas_constant * self.temperature
    #     total_gas_conc = self.pressure / r_t
    #     conc = self.mole_fraction * total_gas_conc
    #
    #     all_gas_conc = np.copy(conc)
    #     sat_conc = self.saturation_pressure / r_t
    #     no_pc_total_gas_conc = total_gas_conc - sat_conc
    #     conc[self.id_pc] = \
    #         np.where(conc[self.id_pc] > sat_conc, sat_conc, conc[self.id_pc])
    #     liquid_conc = all_gas_conc - conc
    #     gas_mole_fraction = self.calc_fraction(conc)
    #     conc = gas_mole_fraction * total_gas_conc
    #     print(conc)
    #     print(all_gas_conc)
    #     return np.maximum(conc, 0.0), all_gas_conc + liquid_conc

    def calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        p_sat = self.saturation_pressure
        self.humidity[:] = \
            self._mole_fraction[:, self.id_pc] * self._pressure / p_sat


def liquid_factory(nx, name, liquid_props, temperature, pressure):
    if isinstance(liquid_props, species.ConstantProperties):
        return ConstantFluid(nx, name, liquid_props, temperature, pressure)
    elif isinstance(liquid_props, species.IncompressibleProperties):
        return IncompressibleFluid(nx, name, liquid_props,
                                   temperature, pressure)
    else:
        raise TypeError('argument liquid_props must be of type '
                        'ConstantProperties or IncompressibleProperties')


def arg_factory(nx, name, liquid_props=None, species_dict=None,
                mole_fractions=None, temperature=293.15,
                pressure=101325.0, **kwargs):
    if species_dict is None:
        return liquid_factory(nx, name, liquid_props, temperature, pressure)
    else:
        species_types = list(species_dict.values())
        species_types_str = ' '.join(species_types)
        if 'gas' in species_types_str and 'liquid' not in species_types_str:
            return GasMixture(nx, name, species_dict, mole_fractions,
                              temperature, pressure)
        elif 'gas' in species_types_str and 'liquid' in species_types_str:
            return TwoPhaseMixture(nx, name, species_dict, mole_fractions,
                                   liquid_props, temperature, pressure,
                                   **kwargs)
        elif 'liquid' in species_types_str and 'gas' \
                not in species_types_str:
            return liquid_factory(nx, name, liquid_props, temperature, pressure)
        else:
            raise NotImplementedError


def factory(fluid_dict):
    nx = fluid_dict['nodes']
    name = fluid_dict['name']
    liquid_props = fluid_dict.get('liquid_props', None)

    components = fluid_dict.get('components', None)
    # species_dict = fluid_dict.get('fluid_components', None)
    # mole_fractions = fluid_dict.get('inlet_composition', None)
    if components is None:
        species_dict = None
        mole_fractions = None
    else:
        if isinstance(components, dict):
            components = OrderedDict(((k, v) for k, v in components.items()))
        else:
            raise TypeError('input variable "components" '
                            'must be a python dictionary')
        species_dict = OrderedDict(((k, v['state'])
                                    for k, v in components.items()))
        mole_fractions = [v['molar_fraction'] for k, v in components.items()]
    temperature = fluid_dict.get('temperature', 293.15)
    pressure = fluid_dict.get('pressure', 101325.0)
    humidity = fluid_dict.get('humidity', None)

    prop_names = \
        ('density', 'viscosity', 'thermal_conductivity', 'specific_heat')
    if any(key in fluid_dict for key in prop_names):
        props_dict = {key: fluid_dict[key] for key in prop_names}
        fluid_props = species.ConstantProperties(name, **props_dict)
        return ConstantFluid(nx, name, fluid_props, temperature, pressure)
    else:
        return arg_factory(nx, name, liquid_props=liquid_props,
                           species_dict=species_dict,
                           mole_fractions=mole_fractions,
                           temperature=temperature,
                           pressure=pressure, humidity=humidity)

# test_species = species.GasSpecies(['O2', 'N2', 'H2'])

# temp = np.array([[300.0, 400.0], [300.0, 400.0]])
# temp = np.array([300.0, 400.0])
# press = np.array([[100000.0, 100000.0], [500000.0, 500000.0]])
# press = 101325.0
# print(species.coeff_dict_dict)
# print(species.coeff_dict_dict2)

# print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
# test = species.calc_thermal_conductivity(temp, press)
# test = species.calc_specific_heat(temp)
# test = species.calc_viscosity(temp)


# print(species.coeff_dict_arr['Thermal Conductivity'][0][0])
# temp = np.linspace(300, 400, 10)
# press = np.linspace(100000, 100000, 10)
# air = Fluid(10, 'air', {'O2': 'gas', 'N2': 'gas'},
#             mole_fractions_init=[0.21, 0.79])
# print(temp)
# print(gas.calc_viscosity(temp))
# print(gas.calc_thermal_conductivity(temp, press))
# print(air.calc_specific_heat(temp))
# print(gas.mw)
# print(gas.calc_density(temp, press))
# liquid_water_props = species.FluidProperties('H2O')
# water = species.PhaseChangeSpecies({'H2O': liquid_water_props})

# print(water.gas_props.calc_specific_heat(temp))
# print(water.calc_saturation_pressure(temp))
# print(water.calc_vaporization_enthalpy(temp))
# print(water.names)
# liquid_water = Fluid(10, 'liquid water', fluid_props=liquid_water_props)
# print(type(liquid_water))
#
# wet_air = Fluid(10, 'wet air', {'O2': 'gas', 'N2': 'gas', 'H2O': 'gas-liquid'},
#                 mole_fractions_init=[0.205, 0.785, 0.01],
#                 liquid_props={'H2O': liquid_water_props})
#
# wet_air.update(temp, press, (1.5, 2.3, 4.0))
#
# print(wet_air.mole_fraction)
