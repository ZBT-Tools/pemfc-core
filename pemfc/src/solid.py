# general imports
import numpy as np
from abc import ABC, ABCMeta, abstractmethod
from collections.abc import Callable

# local modul imports
from pemfc.src import output_object as oo


class MetaAbstractSolid(ABCMeta):
    """
    Definition to avoid metaclass conflicts
    """
    pass


class AbstractSolid(ABC, oo.OutputObject, metaclass=MetaAbstractSolid):
    """
    Abstract Base Class for describing a basic solid layer as a cuboid (width, length,
    thickness) with physical properties (electrical and thermal conductivity).
    Porous materials are considered with effective properties if the 'porosity'
    (default: 0.0) and optionally a Bruggeman exponent (default: 1.5) is provided.
    """

    def __init__(self, layer_dict: dict):
        name = layer_dict.get('name', 'unnamed')
        super().__init__(name)
        self.layer_dict = layer_dict
        self.thickness = layer_dict['thickness']
        self.width = layer_dict['width']
        self.length = layer_dict['length']
        self.porosity = layer_dict.get('porosity', 0.0)
        self.bruggeman_exponent = layer_dict.get('bruggeman_exponent', 1.5)
        self.thermal_conductivity = \
            layer_dict.get('thermal_conductivity', 0.0)
        self.electrical_conductivity = \
            layer_dict.get('electrical_conductivity', 0.0)
        self.area = self.width * self.length

    @abstractmethod
    def calc_conductance(self, conductivity, effective=True):
        pass


class Solid1D(AbstractSolid):
    """
    Class for describing a basic solid layer as a cuboid (width, length,
    thickness) with physical properties (electrical and thermal conductivity).
    Porous materials are considered with effective properties if the 'porosity'
    (default: 0.0) and optionally a Bruggeman exponent (default: 1.5) is provided.
    """
    def __init__(self, layer_dict: dict, discretization: int):
        super().__init__(layer_dict)
        self.dx = discretization
        self.area_dx = self.width * self.dx
        self.thermal_conductance = \
            self.calc_conductance(self.thermal_conductivity)
        self.electrical_conductance = \
            self.calc_conductance(self.electrical_conductivity)

    def calc_conductance(self, conductivity, effective=True):
        if np.ndim(conductivity) == 0:
            conductance_x = \
                self.width * self.thickness * conductivity / self.dx
            conductance_z = self.area_dx * conductivity / self.thickness
        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 2:
            conductance_x = \
                self.width * self.thickness * conductivity[1] / self.dx
            conductance_z = self.area_dx * conductivity[0] / self.thickness
        else:
            raise ValueError('conductivity must be either single scalar or '
                             'an iterable with two entries')
        if effective:
            conductance_z *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_x *= (1.0 - self.porosity) ** self.bruggeman_exponent
        return np.asarray([conductance_z, conductance_x])

    def calc_voltage_loss(self, current_density, area=None, **kwargs):
        if area is None:
            current = current_density * self.area_dx
        else:
            current = current_density * area
        return current / self.electrical_conductance[0]


class Solid2D(AbstractSolid):
    """
    Class for describing a basic solid layer as a cuboid (width, length,
    thickness) with physical properties (electrical and thermal conductivity).
    Porous materials are considered with effective properties if the 'porosity'
    (default: 0.0) and optionally a Bruggeman exponent (default: 1.5) is provided.
    """
    def __init__(self, layer_dict: dict, discretization: tuple):
        super().__init__(layer_dict)
        self.dx = np.asarray(discretization)
        self.area_dx = self.width * self.dx
        self.thermal_conductance = \
            self.calc_conductance(self.thermal_conductivity)
        self.electrical_conductance = \
            self.calc_conductance(self.electrical_conductivity)

    def calc_conductance(self, conductivity, effective=True):
        if np.ndim(conductivity) == 0:
            conductance_x = \
                self.width * self.thickness * conductivity / self.dx
            conductance_z = self.area_dx * conductivity / self.thickness
        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 2:
            conductance_x = \
                self.width * self.thickness * conductivity[1] / self.dx
            conductance_z = self.area_dx * conductivity[0] / self.thickness
        else:
            raise ValueError('conductivity must be either single scalar or '
                             'an iterable with two entries')
        if effective:
            conductance_z *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_x *= (1.0 - self.porosity) ** self.bruggeman_exponent
        return np.asarray([conductance_z, conductance_x])

    def calc_voltage_loss(self, current_density, area=None, **kwargs):
        if area is None:
            current = current_density * self.area_dx
        else:
            current = current_density * area
        return current / self.electrical_conductance[0]


def create_wrapper(object1d: Callable, object2d: Callable, *args, **kwargs):
    if len(args) == 2:
        discretization = args[1]
        if isinstance(discretization, int):
            return object1d(*args, **kwargs)
        elif (isinstance(discretization, (tuple, list))
              and len(discretization) == 2):
            return object2d(*args, **kwargs)
        else:
            raise NotImplementedError("parameter 'discretization' must be "
                                      "either integer or tuple of length 2")


def create(layer_dict: dict, discretization: int | tuple):
    return create_wrapper(Solid1D, Solid2D, layer_dict, discretization)


test_layer_dict = {'name': 'test_solid_2D', 'width': 3.0, 'length': 3.0,
                   'thickness': 0.1, 'electrical_conductivity': 1.0,
                   'thermal_conductivity': 1.0}
solid2d = create(test_layer_dict, (20, 1))
# print('test')