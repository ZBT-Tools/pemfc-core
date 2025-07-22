# general imports
import sys
import numpy as np
from abc import ABC, abstractmethod

# local module imports
from ..output_object import OutputObject
from .. import constants, global_functions as gf
from . import species
from ..global_functions import move_axis
from collections import OrderedDict


try:
    import cantera as ct
    FOUND_CANTERA = True
except ImportError:
    FOUND_CANTERA = False


# class Fluid(ABC, OutputObject):
#     """
#     Abstract factory class to retrieve concrete implementations of the fluid
#     classes.
#     """
#     def __new__(cls, fluid_dict, **kwargs):
#         _get_fluid(fluid_dict, **kwargs)
#
#     def __init__(self, fluid_dict, **kwargs):
#         super().__init__()


class DiscreteFluid(ABC, OutputObject):

    PROPERTY_NAMES = ['density', 'specific_heat', 'viscosity',
                      'thermal_conductivity']
    TYPE_NAME = 'Base Fluid'

    def __init__(self, array_shape, name, temperature=298.15, pressure=101325.0,
                 **kwargs):
        # print("__init__ for Fluid")
        super().__init__(name)

        if isinstance(array_shape, int):
            array_shape = (array_shape,)

        self.array_shape = array_shape
        self.print_variables = \
            {
                'names': ['temperature', 'pressure'],
                'units': ['K', 'Pa'],
                # 'sub_names': ['None', 'None']
            }
        self.print_variables = self.combine_print_variables(
            self.print_variables, kwargs.get('print_variables', None))
        temperature = np.asarray(temperature)
        pressure = np.asarray(pressure)
        self._temperature = np.zeros(self.array_shape)
        self._pressure = np.zeros(self.array_shape)
        self._write_input_to_array(temperature, self._temperature)
        self._write_input_to_array(pressure, self._pressure)

        self.property = dict()
        for name in self.PROPERTY_NAMES:
            self.property[name] = np.zeros(self.array_shape)

    @staticmethod
    def _write_input_to_array(value, array):
        if np.ndim(value) == 0 \
                or value.shape == array.shape:
            array[:] = value
        else:
            raise ValueError('argument value must be scalar or '
                             'have the same shape as 1D array')

    @abstractmethod
    def update(self, temperature, pressure, *args, **kwargs):
        """
        Update local fluid properties depending on temperature and pressure input
        :param temperature: one-dimensional temperature array
        :param pressure: one-dimensional pressure array
        :param args: additional args for concrete implementations
        :param kwargs: additional keyword arguments for concrete implementations
        """
        if temperature is None:
            temperature = self._temperature
        if pressure is None:
            pressure = self._pressure
        if np.any(temperature < 200.0):
            raise ValueError('temperature too low, check boundary conditions')
        if np.any(temperature > 1000.0):
            raise ValueError('temperature too high, check boundary conditions')
        self.temperature = temperature
        self.pressure = pressure

    @abstractmethod
    def _calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass

    def _calc_fraction(self, composition, axis=0, recursion_count=0):
        """
        Calculates mass or mole fractions based on a multi-dimensional
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
                composition = gf.fill_zero_sum(composition, axis=-1)
            else:
                composition = gf.fill_zero_sum(composition, axis=0)
            recursion_count += 1
            if recursion_count > 1:
                raise ValueError('Something wrong in '
                                 '{}.calc_fraction()'.format(self))
            return self._calc_fraction(composition, axis, recursion_count)

    @property
    def temperature(self):
        return self._temperature

    @property
    def pressure(self):
        return self._pressure

    @temperature.setter
    def temperature(self, value):
        self._write_input_to_array(value, self._temperature)

    @pressure.setter
    def pressure(self, value):
        self._write_input_to_array(value, self._pressure)

    @property
    def density(self):
        return self.property['density']

    @property
    def viscosity(self):
        return self.property['viscosity']

    @property
    def thermal_conductivity(self):
        return self.property['thermal_conductivity']

    @property
    def specific_heat(self):
        return self.property['specific_heat']

    # @density.setter
    # def density(self, value):
    #     self.property['density'] = value
    #
    # @viscosity.setter
    # def viscosity(self, value):
    #     self.property['viscosity'] = value
    #
    # @thermal_conductivity.setter
    # def thermal_conductivity(self, value):
    #     self.property['thermal_conductivity'] = value
    #
    # @specific_heat.setter
    # def specific_heat(self, value):
    #     self.property['specific_heat'] = value

    def rescale(self, new_array_shape):
        """
        linearly rescales all numpy arrays with the 1D discretization as first
        dimension
        :param new_array_shape: new discretization along first dimension
        :return: only modification of attributes
        """
        if new_array_shape != self.array_shape:
            attr_list = [a for a in dir(self) if not a.startswith('__')]
            for name in attr_list:
                attr = getattr(self, name)
                if name == 'property':
                    for key in attr.keys():
                        rescaled = self._rescale_attribute(attr[key],
                                                           new_array_shape)
                        if rescaled is not None:
                            attr[key] = rescaled
                else:
                    type_attr = getattr(type(self), name, None)
                    if not isinstance(type_attr, property):
                        rescaled = self._rescale_attribute(attr, new_array_shape)

                        if rescaled is not None:
                            setattr(self, name, rescaled)
            self.array_shape = new_array_shape
        self.add_print_variables(self.print_variables)

    def _rescale_attribute(self, attribute, new_nx):
        if isinstance(attribute, np.ndarray):
            return self._rescale_array(attribute, new_nx)
        else:
            return None

    def _rescale_array(self, array, new_array_shape):
        if new_array_shape != self.array_shape:
            if isinstance(array, np.ndarray):
                if len(array.shape) == 1:
                    first = array[0]
                    last = array[-1]
                    return np.linspace(first, last, new_array_shape)
                elif len(array.shape) == 2:
                    first = array[:, 0]
                    last = array[:, -1]
                    linear_rescaled = \
                        np.array([np.linspace(first[i], last[i],
                                              new_array_shape)
                                  for i in range(len(first))])
                    return linear_rescaled
                else:
                    raise ValueError(
                        'argument array must be one- or two-dimensional')
            else:
                raise TypeError(
                    'argument array must be of type numpy.ndarray')


class ConstantFluid(DiscreteFluid):

    TYPE_NAME = 'Constant Fluid'

    def __init__(self, array_shape, name, fluid_props, temperature=298.15,
                 pressure=101325.0, **kwargs):
        # print("__init__ for IncompressibleFluid")
        super().__init__(array_shape, name, temperature, pressure, **kwargs)
        self.name = name
        if not isinstance(fluid_props, species.ConstantProperties):
            raise TypeError('Argument fluid_props must be of type '
                            'ConstantProperties')
        else:
            self.property['density'][:] = fluid_props.density
            self.property['specific_heat'][:] = fluid_props.specific_heat
            self.property['viscosity'][:] = fluid_props.viscosity
            self.property['thermal_conductivity'][:] = \
                fluid_props.thermal_conductivity
        self.add_print_variables(self.print_variables)
        self.update(self._temperature, self._pressure)

    def update(self, temperature=None, pressure=None, *args, **kwargs):
        super().update(temperature, pressure)

    def _calc_properties(self, temperature, pressure=101325.0, **kwargs):
        pass


class IncompressibleFluid(DiscreteFluid):

    TYPE_NAME = 'Incompressible Fluid'

    def __init__(self, array_shape, name, fluid_props, temperature=298.15,
                 pressure=101325.0, **kwargs):
        # print("__init__ for IncompressibleFluid")
        if isinstance(fluid_props, species.IncompressibleProperties):
            self.properties = fluid_props
        else:
            raise TypeError('Argument fluid_props must be of type '
                            'IncompressibleProperties')
        super().__init__(array_shape, name, temperature, pressure, **kwargs)
        self.add_print_variables(self.print_variables)
        self.update(self._temperature, self._pressure)

    def update(self, temperature=None, pressure=None, *args, **kwargs):
        super().update(temperature, pressure)
        self._calc_properties(self._temperature, self._pressure)

    def calc_property(self, property_name, temperature):
        """
        Wrapper function for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :return: the calculated 1D array of the specific property
        """
        return self.properties.calc_property(property_name, temperature)

    def _calc_properties(self, temperature, pressure=101325.0, **kwargs):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: float or 1D pressure array
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self.calc_property(prop, temperature)


class GasMixture(DiscreteFluid):

    TYPE_NAME = 'Gas Mixture'

    def __init__(self, array_shape, name, species_dict, mole_fractions,
                 temperature=298.15, pressure=101325.0, **kwargs):
        # print("__init__ for Gas Mixture")
        print_variables = \
            {
                'names': ['mole_fraction'],
                'units': ['-'],
                'sub_names': ['self.species.names']
            }
        super().__init__(array_shape, name, temperature, pressure,
                         print_variables=print_variables, **kwargs)
        self.name = name
        species_names = list(species_dict.keys())
        # self.n_species = len(species_names)
        self.gas_constant = constants.GAS_CONSTANT
        self.species = species.GasProperties(species_names)
        self.n_species = len(self.species.names)
        self.species_id = \
            dict(zip(self.species.names, list(range(self.n_species))))
        self.species_viscosity = self.species.calc_viscosity(self._temperature)

        mole_fractions = gf.ensure_list(mole_fractions)
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

        self.array_shape_2d = (self.n_species, *self.array_shape)
        self._mole_fraction = \
            gf.array_vector_multiply(np.ones(self.array_shape_2d),
                                     mole_fractions)
        self.mw = \
            np.sum(gf.array_vector_multiply(self._mole_fraction,
                                            self.species.mw),
                   axis=0)
        mass_fraction_init = \
            mole_fractions * self.species.mw / \
            np.sum(mole_fractions * self.species.mw)
        self._mass_fraction = \
            gf.array_vector_multiply(np.ones(self.array_shape_2d),
                                     mass_fraction_init)
        self._concentration = np.zeros(self.array_shape_2d)

        self.update(self._temperature, self._pressure,)

        print_data = kwargs.get('print_data', True)
        if print_data:
            self.add_print_data(self.mole_fraction, 'Mole Fraction',
                                sub_names=self.species.names)
            self.add_print_variables(self.print_variables)

    @property
    def mole_fraction(self):
        return self._mole_fraction

    @property
    def mass_fraction(self):
        return self._mass_fraction

    @property
    def concentration(self):
        return self._concentration

    @concentration.setter
    def concentration(self, value):
        self._concentration[:] = value

    @property
    def species_mw(self):
        return self.species.mw

    @property
    def species_names(self):
        return self.species.names

    def update(self, temperature=None, pressure=None, mole_composition=None,
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
        self._calc_properties(self._temperature, self._pressure, method)
        self._concentration[:] = self._calc_concentration()

    def _calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        """
        total_mol_conc = self._pressure \
            / (self.gas_constant * self._temperature)
        return self._mole_fraction * total_mol_conc

    def calc_mole_fraction(self, mole_composition):
        if np.min(mole_composition) < 0.0:
            raise ValueError('mole_composition must not be smaller zero')
        mole_fraction = self._calc_fraction(mole_composition)
        return mole_fraction

    def calc_molar_mass(self, mole_fraction=None):
        if mole_fraction is None:
            self.mw[:] = np.sum(gf.array_vector_multiply(self._mole_fraction,
                                                         self.species.mw),
                                axis=0)
            return self.mw
        else:
            return np.sum(gf.array_vector_multiply(mole_fraction,
                                                   self.species.mw),
                          axis=0)

    def calc_mass_fraction(self, mole_fraction, molar_mass=None):
        if molar_mass is None:
            # mw = np.where(self.mw == 0.0, 1.0, self.mw)
            return np.multiply.outer(self.species.mw, 1.0 / self.mw) \
                * mole_fraction
        else:
            # a = np.multiply.outer(self.species.mw, 1.0 / molar_mass)
            return np.multiply.outer(self.species.mw, 1.0 / molar_mass) \
                * mole_fraction

    def _calc_specific_heat(self, temperature, mass_fraction=None):
        species_cp = self.species.calc_specific_heat(temperature)
        if mass_fraction is None:
            return np.sum(self._mass_fraction * species_cp, axis=0)
        else:
            return np.sum(mass_fraction * species_cp, axis=0)

    def _calc_viscosity(self, temperature):
        """
        Calculates the mixture viscosity of a
        gas according to Herning and Zipperer.
        """
        self.species_viscosity[:] = \
            self.species.calc_viscosity(temperature)
        x_sqrt_mw = gf.array_vector_multiply(
            self._mole_fraction, np.sqrt(self.species.mw))
        x_sqrt_mw = np.where(x_sqrt_mw == 0.0, 1.0, x_sqrt_mw)
        return np.sum(self.species_viscosity * x_sqrt_mw, axis=0) \
            / np.sum(x_sqrt_mw, axis=0)

    def _calc_wilke_coefficients(self):
        """
        Calculates the wilke coefficients for
        each species combination of a gas.
        """
        visc = self.species_viscosity
        mw = self.species.mw
        alpha = []
        for i in range(self.n_species):
            beta = []
            for j in range(self.n_species):
                a = (1. + np.sqrt((visc[i] / visc[j]))
                     * (mw[j] / mw[i]) ** 0.25) ** 2.0
                b = (np.sqrt(8.) * (1. + mw[j] / mw[i])) ** -0.5
                beta.append(a / b)
            alpha.append(beta)
        return np.asarray(alpha)

    def _calc_thermal_conductivity(self, temperature, pressure):
        """
        Calculates the heat conductivity of a gas mixture,
        according to Wilkes equation.
        """
        self._mole_fraction[1:] = np.maximum(1e-16, self._mole_fraction[1:])
        wilke_coeffs = self._calc_wilke_coefficients()
        lambda_species = \
            self.species.calc_thermal_conductivity(temperature, pressure)
        lambda_mix = np.zeros(temperature.shape)
        for i in range(self.n_species):
            a = self._mole_fraction[i] * lambda_species[i]
            b = np.sum(self._mole_fraction * wilke_coeffs[i],
                       axis=0)
            b += 1e-16
            lambda_mix += a / b
        return lambda_mix

    def _calc_density(self, temperature, pressure, method="ideal"):
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

    def _calc_property(self, property_name, temperature, pressure=101325.0,
                       method='ideal'):
        """
        Wrapper function for the native property functions
        :param property_name: name of self.PROPERTY_NAMES to calculate
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        if property_name == 'specific_heat':
            return self._calc_specific_heat(temperature)
        elif property_name == 'viscosity':
            return self._calc_viscosity(temperature)
        elif property_name == 'thermal_conductivity':
            return self._calc_thermal_conductivity(temperature, pressure)
        elif property_name == 'density':
            return self._calc_density(temperature, pressure, method)
        else:
            raise ValueError('property_name '
                             '{} not valid'.format(property_name))

    def _calc_properties(self, temperature, pressure=101325.0, method='ideal'):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self._calc_property(prop, temperature, pressure, method)


class CanteraGasMixture(DiscreteFluid):
    """
    Wrapper for discrete fluid properties calculated by Cantera library
    """
    def __init__(self, array_shape, name, species_dict, mole_fractions,
                 temperature=298.15, pressure=101325.0, **kwargs):
        super().__init__(array_shape, name, temperature, pressure, **kwargs)

        self.dict = kwargs.get('fluid_dict')
        self.species_names = list(species_dict.keys())
        self.n_species = len(self.species_names)

        # Property name reference to Cantera names
        self.prop_name_ref = \
            {'density': 'density',
             'viscosity': 'viscosity',
             'thermal_conductivity': 'thermal_conductivity',
             'specific_heat': 'cp'}
        # Create Cantera gas solution
        self.solution = ct.Solution('gri30.yaml')
        # Set initial species and their mole fractions
        self.solution.X = dict(zip(self.species_names, mole_fractions))
        # Initialize state
        self.solution.TP = temperature, pressure
        # Create an array based solution
        self.solution_array = ct.SolutionArray(self.solution, array_shape)

        self.all_species_names = self.solution.species_names
        # Get ids for used species from complete list
        self.species_ids = [self.all_species_names.index(item) for item
                            in self.species_names]
        self.species_mw = \
            self.solution.molecular_weights[self.species_ids] * 1e-3
        self._mole_fraction = self._get_composition()
        self._mass_fraction = self._get_composition(mole_fractions=False)
        self.update(temperature=self._temperature, pressure=self._pressure)

    def update(self, temperature=None, pressure=None, mole_composition=None,
               *args, **kwargs):
        super().update(temperature, pressure)
        self.solution_array.TP = self._temperature, self._pressure
        if mole_composition is not None:
            self.solution_array.X = self._set_composition(mole_composition)
        self._calc_properties(self._temperature, self._pressure)
        self._mass_fraction = self._get_composition(mole_fractions=False)

    def _set_composition(self, mole_composition):
        if np.shape(mole_composition)[0] != self.n_species:
            raise ValueError('First dimension of composition must be '
                             'equal to number of species')
        self.solution_array.X[:, self.species_ids] = move_axis(
            mole_composition, 0, -1)

    def _get_composition(self, mole_fractions=True):
        if mole_fractions:
            return move_axis(self.solution_array.X, -1, 0)[self.species_ids]
        else:
            return move_axis(self.solution_array.Y, -1, 0)[self.species_ids]

    def rescale(self, new_array_shape):
        super().rescale(new_array_shape)
        self.solution_array(self.solution, new_array_shape)
        self.update(self._temperature, self._pressure, self._mole_fraction)

    def _calc_properties(self, temperature, pressure=101325.0, **kwargs):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: float or 1D pressure array
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop] = \
               getattr(self.solution_array, self.prop_name_ref[prop])

    @property
    def mole_fraction(self):
        return self._mole_fraction

    @property
    def mass_fraction(self):
        return self._mass_fraction

    @property
    def concentration(self):
        return self.solution_array.concentrations[:, self.species_ids] * 1e3

    @concentration.setter
    def concentration(self, value):
        self.solution_array.concentrations[:, self.species_ids] = value * 1000.0


class TwoPhaseMixture(DiscreteFluid):

    TYPE_NAME = 'Two-Phase Mixture'

    def __init__(self, array_shape, name, species_dict, mole_fractions,
                 liquid_props=None, temperature=298.15, pressure=101325.0,
                 **kwargs):
        # print("__init__ for TwoPhaseMixture")

        super().__init__(array_shape, name, temperature, pressure, **kwargs)
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
        self.liquid = IncompressibleFluid(
            self.array_shape, name + ': Liquid Phase', fluid_props=liquid_props,
            temperature=self._temperature, pressure=self._pressure)
        self.gas = GasMixture(
            self.array_shape, name + ': Gas Phase',
            species_dict=gas_species_dict, mole_fractions=mole_fractions,
            temperature=self._temperature, pressure=self._pressure,
            print_data=False)
        ids_pc = [self.species.names.index(name) for name in
                  phase_change_species_names]
        if len(ids_pc) > 1:
            raise NotImplementedError('at the moment only one species '
                                      'undergoing phase change is allowed')
        self.id_pc = ids_pc[0]
        all_ids = np.array(list(range(len(self.species.names))), dtype='int32')
        ids_no_pc = np.delete(all_ids, ids_pc)
        self.ids_no_pc = list(ids_no_pc)
        self.n_species = len(species_dict)

        # update mole_fractions if humidity is provided
        humidity = kwargs.get('humidity', None)
        if humidity is not None:
            mole_fractions = \
                self.phase_change_species.calc_humid_composition(
                    humidity, temperature, pressure, mole_fractions, self.id_pc)

        # Total properties (both phases)
        self.array_shape_2d = self.gas.array_shape_2d
        self._mole_fraction = \
            gf.array_vector_multiply(np.ones(self.array_shape_2d),
                                     mole_fractions)
        mass_fraction_init = \
            mole_fractions * self.gas.species.mw / \
            np.sum(mole_fractions * self.gas.species.mw)
        self._mass_fraction = \
            gf.array_vector_multiply(np.ones(self.array_shape_2d),
                                     mass_fraction_init)
        self.mw = np.zeros(self.array_shape)

        self.liquid_mass_fraction = np.zeros(self.array_shape)
        self.liquid_mole_fraction = np.zeros(self.array_shape)
        self.humidity = np.zeros(self.array_shape)
        self.saturation_pressure = np.zeros(self.array_shape)
        self.surface_tension = np.zeros(self.array_shape)

        # Print data
        print_variables = \
            {
                'names': ['humidity', 'mole_fraction'],
                'units': ['-', '-'],
                'sub_names': ['None', 'self.species.names']
            }
        self.add_print_variables(print_variables)
        # self.add_print_data(self.humidity, 'Humidity')
        # print(self.print_data)
        # print(self.gas.print_data)
        self.update(temperature=self._temperature, pressure=self._pressure)

    def update(self, temperature=None, pressure=None,
               mole_composition=None, gas_mole_composition=None,
               method='ideal', *args, **kwargs):
        super().update(temperature, pressure)
        if mole_composition is not None:
            if np.max(mole_composition) > constants.SMALL:
                mole_composition[mole_composition < constants.SMALL] = 0.0
                self._mole_fraction[:] = \
                    self.gas.calc_mole_fraction(mole_composition)

        self.mw[:] = self.gas.calc_molar_mass(self.mole_fraction)
        test = self.mw[0]
        self._mass_fraction[:] = \
            self.gas.calc_mass_fraction(self._mole_fraction, self.mw)
        self.saturation_pressure[:] = \
            self.phase_change_species.calc_saturation_pressure(self._temperature)
        if gas_mole_composition is None \
                or np.max(gas_mole_composition) < constants.SMALL:
            gas_mole_composition = self._calc_concentration()
        if mole_composition is None:
            mole_composition = gas_mole_composition
        self.gas.update(self._temperature,
                        self._pressure, gas_mole_composition, method)
        self.liquid.update(self._temperature, self._pressure)
        self.liquid_mole_fraction[:] = \
            np.sum(mole_composition, axis=0) \
            - np.sum(gas_mole_composition, axis=0)
        self.liquid_mole_fraction[self.liquid_mole_fraction < 0.0] = 0.0
        self.liquid_mass_fraction = \
            np.sum(mole_composition * self.mw, axis=0) \
            - np.sum(gas_mole_composition * self.gas.mw, axis=0)
        self.liquid_mass_fraction[self.liquid_mass_fraction < 0.0] = 0.0
        self.surface_tension[:] = self.calc_surface_tension()

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

        self._calc_properties(temperature, pressure, method)
        self._calc_humidity()

    def rescale(self, new_array_shape):
        self.gas.rescale(new_array_shape)
        self.liquid.rescale(new_array_shape)
        super().rescale(new_array_shape)
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
        self._write_input_to_array(value, self._temperature)
        self.gas._temperature = self._temperature
        self.liquid._temperature = self._temperature

    @pressure.setter
    def pressure(self, value):
        self._write_input_to_array(value, self._pressure)
        self.gas._pressure = self._pressure
        self.liquid._pressure = self._pressure

    @property
    def mole_fraction(self):
        return self._mole_fraction

    @property
    def mass_fraction(self):
        return self._mass_fraction

    @property
    def species(self):
        return self.gas.species

    @property
    def species_id(self):
        return self.gas.species_id

    @property
    def species_mw(self):
        return self.gas.species_mw

    @property
    def species_names(self):
        return self.gas.species_names

    def _set_name(self, name):
        super()._set_name(name)
        self.gas.name = self.name + ': Gas Phase'
        self.liquid.name = self.name + ': Liquid Phase'

    def calc_surface_tension(self, temperature=None):
        if temperature is None:
            temperature = self.temperature
        return self.phase_change_species.calc_surface_tension(temperature)

    def _calc_property(self, property_name, temperature, pressure=101325.0,
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
        if property_name == 'specific_heat':
            return self._calc_specific_heat()
        else:
            return self.gas.property[property_name]

    def _calc_properties(self, temperature, pressure=101325.0, method='ideal'):
        """
        Wrapper function to calculate the classes properties
        :param temperature: 1D temperature array
        :param pressure: 1D pressure array
        :param method: method density calculation (at the moment only ideal)
        :return: the calculated 1D array of the specific property
        """
        for prop in self.PROPERTY_NAMES:
            self.property[prop][:] = \
                self._calc_property(prop, temperature, pressure, method)

    def _calc_specific_heat(self):
        return self.liquid_mass_fraction \
            * self.liquid.specific_heat \
            + (1.0 - self.liquid_mass_fraction) * self.gas.specific_heat

    def _calc_concentration(self):
        """
        Calculates the gas phase molar concentrations.
        :return: gas phase concentration
        (n_species x n_nodes)
        """
        r_t = self.gas.gas_constant * self._temperature
        total_gas_conc = self._pressure / r_t
        conc = self._mole_fraction * total_gas_conc
        # all_conc = np.copy(conc)
        sat_conc = self.saturation_pressure / r_t
        dry_mole_fraction = np.copy(self.mole_fraction)
        dry_mole_fraction[self.id_pc] = 0.0
        dry_mole_fraction = self._calc_fraction(dry_mole_fraction)
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

    def _calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        p_sat = self.saturation_pressure
        self.humidity[:] = \
            self._mole_fraction[self.id_pc] * self._pressure / p_sat

        # Simplification of limiting humidity to 100%
        self.humidity[self.humidity > 1.0] = 1.0

    def calc_vaporization_enthalpy(self, temperature=None):
        if temperature is None:
            temperature = self._temperature
        return self.phase_change_species.calc_vaporization_enthalpy(temperature)


class CanteraTwoPhaseMixture(CanteraGasMixture):
    """
    Wrapper for discrete fluid properties calculated by Cantera library
    """
    def __init__(self, array_shape, name, species_dict, mole_fractions,
                 temperature=298.15, pressure=101325.0, **kwargs):

        phase_change_species_names = \
            [key for key in species_dict
             if 'gas' and 'liquid' in species_dict[key]]
        species_names = list(species_dict.keys())

        ids_pc = [species_names.index(name) for name in
                  phase_change_species_names]
        self.id_pc = ids_pc[0]
        if len(ids_pc) > 1:
            raise NotImplementedError('at the moment only one species '
                                      'undergoing phase change is allowed')

        self.water_base = ct.Water()
        # if isinstance(temperature, (np.ndarray, list, tuple)):
        #     self.water.TP = temperature[0], pressure[0]
        # else:
        #     self.water.TP = temperature
        # self.water = ct.SolutionArray(self.water_base, shape=array_shape)
        # update mole_fractions if humidity is provided
        humidity = kwargs.get('humidity', None)

        if humidity is not None:
            self.humidity = np.ones(array_shape) * humidity
            mole_fractions = \
                self.calc_humid_composition(
                    humidity, temperature, pressure, mole_fractions, self.id_pc)
        else:
            self.humidity = np.zeros(array_shape)
        self.saturation_pressure = np.zeros(array_shape)
        super().__init__(array_shape, name, species_dict, mole_fractions,
                         temperature, pressure, **kwargs)
        self.update(self._temperature, self._pressure, mole_fractions)

    def update(self, temperature=None, pressure=None,
               mole_composition=None, gas_mole_composition=None,
               *args, **kwargs):
        super().update(temperature, pressure, *args, **kwargs)
        # self.water.TP = self._temperature, self._pressure
        temperature = self._temperature.ravel()
        pressure = self._pressure.ravel()
        sat_pressure = self.saturation_pressure.ravel()
        for i in range(len(sat_pressure)):
            self.water_base.TP = temperature[i], pressure[i]
            sat_pressure[i] = self.water_base.P_sat
        # self.saturation_pressure[:] = sat_pressure.reshape(self.array_shape)
        self._calc_humidity()
        # self.water.TP = self._temperature, self._pressure

    def calc_humid_composition(self, humidity, temperature, pressure,
                               dry_molar_composition, id_pc=-1):
        if humidity > 1.0:
            raise ValueError('relative humidity must not exceed 1.0')
        self.water_base.TP = temperature, pressure
        molar_fraction_water = humidity * self.water_base.P_sat / pressure
        humid_composition = np.asarray(dry_molar_composition)
        humid_composition[id_pc] = 0.0
        humid_composition /= np.sum(humid_composition, axis=0)
        humid_composition *= (1.0 - molar_fraction_water)
        humid_composition[id_pc] = molar_fraction_water
        return humid_composition

    def _calc_humidity(self):
        """
        Calculates the relative humidity of the fluid.
        """
        self.humidity[:] = self._mole_fraction[self.id_pc] \
            * self._pressure / self.saturation_pressure

    def calc_vaporization_enthalpy(self, temperature=None):
        """
        Returns molar enthalpy of vaporization (J/mol)
        """
        if temperature is None:
            temperature = self._temperature

        vap_enthalpy = np.zeros(self.array_shape)
        vap_enth_ravel = vap_enthalpy.ravel()
        temperature = np.asarray(gf.ensure_list(temperature)).ravel()
        for i in range(len(vap_enth_ravel)):
            self.water_base.TQ = temperature[i], 0.0
            enthalpy_sat_liq = self.water_base.enthalpy_mole
            self.water_base.TQ = temperature[i], 1.0
            enthalpy_sat_vap = self.water_base.enthalpy_mole
            vap_enth_ravel[i] = np.abs(enthalpy_sat_vap - enthalpy_sat_liq)
        # self.water_base.TP = self._temperature[0], self._pressure[0]
        return vap_enthalpy * 1e-3


def liquid_factory(nx, name, liquid_props, temperature, pressure):
    if isinstance(liquid_props, species.ConstantProperties):
        return ConstantFluid(nx, name, liquid_props, temperature, pressure)
    elif isinstance(liquid_props, species.IncompressibleProperties):
        return IncompressibleFluid(nx, name, liquid_props,
                                   temperature, pressure)
    else:
        raise TypeError('argument liquid_props must be of type '
                        'ConstantProperties or IncompressibleProperties')


def _arg_factory(array_shape, name, liquid_props=None, species_dict=None,
                 mole_fractions=None, temperature=293.15,
                 pressure=101325.0, backend='pemfc', **kwargs):
    if species_dict is None:
        return liquid_factory(array_shape, name, liquid_props,
                              temperature, pressure)
    else:
        species_types = list(species_dict.values())
        species_types_str = ' '.join(species_types)
        if 'gas' in species_types_str and 'liquid' not in species_types_str:
            if backend == 'pemfc':
                return GasMixture(
                    array_shape, name, species_dict, mole_fractions,
                    temperature, pressure, **kwargs)
            elif backend == 'cantera':
                return CanteraGasMixture(
                    array_shape, name, species_dict, mole_fractions,
                    temperature, pressure, **kwargs)
            else:
                raise NotImplementedError
        elif 'gas' in species_types_str and 'liquid' in species_types_str:
            if backend == 'pemfc':
                return TwoPhaseMixture(
                    array_shape, name, species_dict, mole_fractions,
                    liquid_props, temperature, pressure, **kwargs)
            elif backend == 'cantera':
                return CanteraTwoPhaseMixture(
                    array_shape, name, species_dict, mole_fractions,
                    temperature, pressure, **kwargs)
        elif 'liquid' in species_types_str and 'gas' \
                not in species_types_str:
            return liquid_factory(array_shape, name, liquid_props,
                                  temperature, pressure)
        else:
            raise NotImplementedError


def create(fluid_dict, **kwargs):

    backend = kwargs.pop('backend', 'pemfc')
    if not FOUND_CANTERA and backend == 'cantera':
        print('Warning: Cantera module not found, '
              'using internal gas property library')
        backend = 'pemfc'

    array_shape = fluid_dict['nodes']
    if isinstance(array_shape, int):
        array_shape = (array_shape,)
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
        return ConstantFluid(array_shape, name, fluid_props,
                             temperature, pressure)
    else:
        return _arg_factory(
            array_shape, name, liquid_props=liquid_props,
            species_dict=species_dict,
            mole_fractions=mole_fractions, temperature=temperature,
            pressure=pressure, humidity=humidity, backend=backend,
            fluid_dict=fluid_dict)
