# general imports
import numpy as np
from abc import ABC, ABCMeta, abstractmethod
from collections.abc import Callable
# local modul imports
from pemfc.src import output_object as oo
from pemfc.src import discretization as dsct


class SolidLayer(oo.OutputObject):
    """
    Class for describing a basic solid layer as a cuboid (width,
    length, thickness) with physical properties (electrical and thermal
    conductivity). Porous materials are considered with effective properties if
    the 'porosity' (default: 0.0) and optionally a Bruggeman exponent
    (default: 1.5) is provided.
    """

    def __init__(self, layer_dict: dict, discretization: dsct.Discretization2D):
        # Initialize super class
        name = layer_dict.get('name', 'unnamed')

        super().__init__(name=name)
        self.dsct = discretization
        self.dict = layer_dict
        self.thickness = layer_dict['thickness']

        self.porosity = layer_dict.get('porosity', 0.0)
        self.bruggeman_exponent = layer_dict.get('bruggeman_exponent', 1.5)
        self.thermal_conductivity = \
            layer_dict.get('thermal_conductivity', 0.0)
        self.electrical_conductivity = \
            layer_dict.get('electrical_conductivity', 0.0)
        self.thermal_conductance = \
            self.calc_conductance(self.thermal_conductivity)
        self.electrical_conductance = \
            self.calc_conductance(self.electrical_conductivity)

    def calc_conductance(self, conductivity, effective=True):
        if np.ndim(conductivity) == 0:
            conductance_z = \
                self.dsct.d_area * conductivity / self.thickness
            conductance_x = \
                self.dsct.dx[1] / self.dsct.dx[0] \
                * (self.thickness * conductivity)
            conductance_y = \
                self.dsct.dx[0] / self.dsct.dx[1] \
                * self.thickness * conductivity

        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 2:
            conductance_z = self.dsct.d_area * conductivity[0] \
                            / self.thickness
            conductance_x = \
                self.dsct.dx[1] / self.dsct.dx[0] \
                * self.thickness * conductivity[1]
            conductance_y = \
                self.dsct.dx[0] / self.dsct.dx[1] * self.thickness * conductivity[1]
        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 3:
            conductance_z = self.dsct.d_area * conductivity[0] / self.thickness
            conductance_x = \
                self.dsct.dx[1] / self.dsct.dx[0] * self.thickness * conductivity[1]
            conductance_y = \
                self.dsct.dx[0] / self.dsct.dx[1] * self.thickness * conductivity[2]
        else:
            raise ValueError('conductivity must be either single scalar or '
                             'an iterable with two entries')
        if effective:
            conductance_z *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_x *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_y *= (1.0 - self.porosity) ** self.bruggeman_exponent

        return np.asarray([conductance_z, conductance_x, conductance_y])

    def calc_voltage_loss(self, current_density, area=None, axis=0, **kwargs):
        if area is None:
            current = current_density * self.dsct.d_area
        else:
            current = current_density * area
        return current / self.electrical_conductance[axis]


# solid_dict = {
#     'width': 1.0,
#     'length': 10.0,
#     'thickness': 0.1,
#     'electrical_conductivity': 1.0,
#     'thermal_conductivity': 1.0,
#     'ionic_conductivity': 15.0,
#     'type': 'Constant',
# }
#
# test_solid = SolidLayer(solid_dict, (10, 2))
# print(test_solid.dict)
