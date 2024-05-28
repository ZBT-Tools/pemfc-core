# general imports
import numpy as np
from abc import ABC, ABCMeta, abstractmethod
from collections.abc import Callable
from collections.abc import Iterable
# local modul imports
from pemfc.src import output_object as oo
from pemfc.src import discretization as dsct
from typing import Self


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
        self.conductance = {'electrical': self.electrical_conductance,
                            'thermal': self.thermal_conductance}

    def calc_conductance(self, conductivity, effective=False):
        if np.ndim(conductivity) == 0:
            conductance_x = \
                self.dsct.d_area * conductivity / self.thickness
            conductance_y = \
                self.dsct.dx[1] / self.dsct.dx[0] \
                * (self.thickness * conductivity)
            conductance_z = \
                self.dsct.dx[0] / self.dsct.dx[1] \
                * self.thickness * conductivity

        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 2:
            conductance_x = self.dsct.d_area * conductivity[0] \
                            / self.thickness
            conductance_y = \
                self.dsct.dx[1] / self.dsct.dx[0] \
                * self.thickness * conductivity[1]
            conductance_z = \
                self.dsct.dx[0] / self.dsct.dx[1] * self.thickness * conductivity[1]
        elif np.ndim(conductivity) == 1 and np.shape(conductivity)[0] == 3:
            conductance_x = self.dsct.d_area * conductivity[0] / self.thickness
            conductance_y = \
                self.dsct.dx[1] / self.dsct.dx[0] * self.thickness * conductivity[1]
            conductance_z = \
                self.dsct.dx[0] / self.dsct.dx[1] * self.thickness * conductivity[2]
        else:
            raise ValueError('conductivity must be either single scalar or '
                             'an iterable with two entries')
        if effective:
            conductance_x *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_y *= (1.0 - self.porosity) ** self.bruggeman_exponent
            conductance_z *= (1.0 - self.porosity) ** self.bruggeman_exponent

        return np.asarray([conductance_x, conductance_y, conductance_z])

    def calc_voltage_loss(self, current_density, area=None, axis=0, **kwargs):
        if area is None:
            current = current_density * self.dsct.d_area
        else:
            current = current_density * area
        return current / self.electrical_conductance[axis]

    def reduce_conductance(self, factor, indices, axis=0):
        if axis == 0:
            for conductance in self.conductance:
                conductance[indices, :, :] *= factor
        elif axis == 1:
            for conductance in self.conductance:
                conductance[:, indices, :] *= factor
        elif axis in (2, -1):
            for conductance in self.conductance:
                conductance[:, :, indices] *= factor
        else:
            raise ValueError('only three-dimensional arrays supported')

    @classmethod
    def connect(cls, conductance_a: Iterable[Self],
                conductance_b: Iterable[Self], mode):
        conductance_a = np.asarray(conductance_a)
        conductance_b = np.asarray(conductance_b)
        if mode == 'serial':
            return 1.0 / (1.0 / (conductance_a + conductance_b))
        elif mode == 'parallel':
            return conductance_a + conductance_b
        else:
            raise ValueError("argument 'mode' must either be 'serial' or 'parallel'")

    @classmethod
    def calc_inter_node_conductance(cls, conductance: Iterable[Self],
                                    axis=0, mode='serial'):
        conductance_by_2 = 2.0 * np.asarray(conductance)
        if axis == 0:
            return cls.connect(conductance_by_2[:-1, :, :],
                               conductance_by_2[1:, :, :],
                               mode)
        elif axis == 1:
            return cls.connect(conductance_by_2[:, :-1, :],
                               conductance_by_2[:, 1:, :],
                               mode)
        elif axis in (2, -1):
            return cls.connect(conductance_by_2[:, :, :-1],
                               conductance_by_2[:, :, 1:],
                               mode)
        else:
            raise ValueError('axis argument must be either 0, 1, 2 (-1)')


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
