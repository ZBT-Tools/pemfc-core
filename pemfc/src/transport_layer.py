# General imports
import numpy as np
from abc import ABC, abstractmethod
# from collections.abc import Callable
from collections.abc import Iterable
# Local modul imports
from . import output_object as oo
from . import discretization as dsct
from . import matrix_functions as mf
from typing import Self


class TransportLayer(oo.OutputObject2D, ABC):
    """
    Abstract class as an interface for describing a basic solid or porous layer
    as a cuboid (width, length, thickness) with physical properties
    (e.g. electrical conductivity, thermal conductivity, diffusion
    coefficient).
    Porous materials are considered with effective properties if the
    'porosity' (default: 0.0) and optionally a Bruggeman exponent
    (default: 1.5) is provided.
    """

    def __init__(self, input_dict: dict, transport_properties: dict,
                 discretization: dsct.Discretization, *args, **kwargs):
        # Initialize super class
        name = input_dict.get('name', 'unnamed')

        super().__init__(name=name)
        self.types = [key for key in transport_properties]
        self.shift_axis = input_dict.get('shift_axis', None)
        self.discretization = discretization
        self.dict = input_dict
        self.thickness: float = 0.0
        self.transport_properties = transport_properties
        self.volume_fraction = input_dict.get('volume_fraction', 1.0)
        self.initial_volume_fraction = np.copy(self.volume_fraction)
        self.bruggeman_exponent = input_dict.get('bruggeman_exponent', 1.5)
        self.effective = input_dict.get('effective', False)
        self.geometric_factors: np.ndarray = np.asarray([0.0, 0.0, 0.0])
        self.conductance: dict = {}
        self.d_volume = self.discretization.d_volume
        self.dx = self.discretization.dx
        self.first_update = True

    @classmethod
    def create(cls, input_dict: dict, transport_properties: dict,
               discretization: dsct.Discretization, *args, **kwargs):
        if isinstance(discretization, dsct.Discretization2D):
            return TransportLayer2D(input_dict, transport_properties,
                                    discretization)
        elif isinstance(discretization, dsct.Discretization3D):
            return TransportLayer3D(input_dict, transport_properties,
                                    discretization)
        else:
            raise NotImplementedError("argument 'discretization' must be of "
                                      "type 'Discretization'")

    @abstractmethod
    def calc_geometric_factors(self, effective=True, *args, **kwargs):
        pass

    def modify_geometric_factors(self, result, effective=True, *args, **kwargs):
        if effective:
            result *= self.volume_fraction ** self.bruggeman_exponent
        if self.shift_axis is not None:
            result = mf.shift_nodes(result, axis=self.shift_axis,
                                    include_axis=True, inverse=True)
            self.d_volume = mf.shift_nodes(self.d_volume, axis=self.shift_axis,
                                           include_axis=True,
                                           inverse=False,
                                           except_first_axis=False)
            self.dx = mf.shift_nodes(self.dx, axis=self.shift_axis,
                                     include_axis=True, inverse=False)
        return np.asarray(result)

    def calc_conductance(self, conductivity, effective=False):
        conductivity = np.asarray(conductivity)
        if np.ndim(conductivity) == 0:
            conductance = self.geometric_factors * conductivity
        elif np.ndim(conductivity) == 1:
            if len(conductivity) == 2:
                conductivity = np.asarray(
                    [conductivity[0], conductivity[1], conductivity[1]])
                conductance = (
                    conductivity * self.geometric_factors.transpose()).transpose()
            elif len(conductivity) == 3:
                conductance = (
                    conductivity * self.geometric_factors.transpose()).transpose()
            else:
                raise ValueError('conductivity array is limited to three '
                                 'entries in first dimension')
        elif (np.ndim(conductivity) == np.ndim(self.geometric_factors) and
                conductivity.shape == self.geometric_factors.shape):
            conductance = conductivity * self.geometric_factors
        elif conductivity.shape == self.geometric_factors[0].shape:
            conductivity = np.asarray([conductivity for i in range(len(
                self.geometric_factors))])
            conductance = conductivity * self.geometric_factors
        else:
            raise ValueError('conductivity must be either single scalar,'
                             'an iterable with two (tp, ip) or three (x, '
                             'y, z) entries, or a full array with the same '
                             'shape as the attribute "geometric_factors"')
        if effective:
            conductance *= self.volume_fraction ** self.bruggeman_exponent
        return conductance

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
    def connect(cls, conductance_a: (list, np.ndarray),
                conductance_b: (list, np.ndarray), mode):
        conductance_a = np.asarray(conductance_a)
        conductance_b = np.asarray(conductance_b)
        if mode == 'serial':
            inf_array = np.ones(conductance_a.shape) * 1e16
            zero_array = np.zeros(conductance_a.shape)
            resistance_a = np.divide(
                1.0, conductance_a, out=inf_array, where=conductance_a != 0.0)
            resistance_b = np.divide(
                1.0, conductance_b, out=inf_array, where=conductance_a != 0.0)
            resistance_sum = resistance_a + resistance_b
            return np.divide(
                1.0, resistance_sum, out=zero_array,
                where=resistance_sum != 0.0)
        elif mode == 'parallel':
            return conductance_a + conductance_b
        else:
            raise ValueError("argument 'mode' must either be "
                             "'serial' or 'parallel'")

    @classmethod
    def calc_inter_node_conductance(
            cls, conductance: (list, np.ndarray), axis,
            mode='serial', shifted=False):
        if shifted:
            conductance_by_2 = np.copy(conductance)
            np.moveaxis(conductance_by_2, axis, 0)[1:-1] *= 2.0
        else:
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

    @staticmethod
    def calc_average_inter_node_value(values: (list, np.ndarray), axis):
        return mf.transform_nodes(values, axis, mode='average')

    def update(self, transport_properties: dict,
               volume_fraction: np.ndarray = None, *args, **kwargs):
        if volume_fraction is not None and self.first_update:
            self.geometric_factors = self.calc_geometric_factors(
                effective=False)
            self.effective = True
            self.first_update = False
        for key in self.conductance:
            self.transport_properties.update(transport_properties)
            self.conductance[key][:] = self.calc_conductance(
                transport_properties[key], effective=self.effective)


class TransportLayer2D(TransportLayer):
    """
    Class to describe a basic solid or porous layer
    as a cuboid (width, length, thickness) with physical properties
    (e.g. electrical conductivity, thermal conductivity, diffusion
    coefficient).
    Porous materials are considered with effective properties if the
    'porosity' (default: 0.0) and optionally a Bruggeman exponent
    (default: 1.5) is provided.
    """

    def __init__(self, input_dict: dict, transport_properties: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        # Initialize super class
        super().__init__(input_dict, transport_properties, discretization,
                         *args, **kwargs)
        self.thickness = input_dict['thickness']
        self.geometric_factors = self.calc_geometric_factors(
            effective=self.effective)
        self.conductance = {key: self.calc_conductance(value) for key, value
                            in transport_properties.items()}

    def calc_geometric_factors(self, effective=True, *args, **kwargs):
        result = np.asarray([
            self.discretization.d_area / self.thickness,
            self.discretization.dx[1] / self.discretization.dx[0] * self.thickness,
            self.discretization.dx[0] / self.discretization.dx[1] * self.thickness
            ])
        return self.modify_geometric_factors(result, effective=effective)


class TransportLayer3D(TransportLayer):
    """
    Class to describe a basic solid or porous layer
    as a cuboid (depth, length, width) with physical properties
    (e.g. electrical conductivity, thermal conductivity, diffusion
    coefficient).
    Porous materials are considered with effective properties if the
    'porosity' (default: 0.0) and optionally a Bruggeman exponent
    (default: 1.5) is provided.
    """

    def __init__(self, input_dict: dict, transport_properties: dict,
                 discretization: dsct.Discretization3D, *args, **kwargs):
        # Initialize super class
        super().__init__(input_dict, transport_properties, discretization,
                         *args, **kwargs)
        self.geometric_factors = self.calc_geometric_factors(
            effective=self.effective)
        self.conductance = {key: self.calc_conductance(value) for key, value
                            in transport_properties.items()}
        print('test')

    def calc_geometric_factors(self, effective=True, *args, **kwargs):
        result = np.asarray([
            self.discretization.d_area[0] / self.discretization.dx[0],
            self.discretization.d_area[1] / self.discretization.dx[1],
            self.discretization.d_area[2] / self.discretization.dx[2]
            ])
        return self.modify_geometric_factors(result, effective=effective,
                                             *args, **kwargs)




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
