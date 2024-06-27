# General imports
import numpy as np
from abc import ABC, abstractmethod
# from collections.abc import Callable
from collections.abc import Iterable
# Local modul imports
from . import output_object as oo
from . import discretization as dsct
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

    # def __new__(cls, input_dict: dict, transport_properties: dict,
    #             discretization: dsct.Discretization, *args, **kwargs):
    #     if isinstance(discretization, dsct.Discretization2D):
    #         return super(TransportLayer, cls).__new__(TransportLayer2D)
    #     elif isinstance(discretization, dsct.Discretization3D):
    #         return super(TransportLayer, cls).__new__(TransportLayer3D)
    #     else:
    #         raise NotImplementedError("argument 'discretization' must be of "
    #                                   "type 'Discretization'")

    def __init__(self, input_dict: dict, transport_properties: dict,
                 discretization: dsct.Discretization, *args, **kwargs):
        # Initialize super class
        name = input_dict.get('name', 'unnamed')

        super().__init__(name=name)
        self.discretization = discretization
        self.dict = input_dict
        self.thickness: float = 0.0
        self.transport_properties = transport_properties
        self.porosity = input_dict.get('porosity', 0.0)
        self.bruggeman_exponent = input_dict.get('bruggeman_exponent', 1.5)
        self.effective = input_dict.get('effective', False)
        self.geometric_factors: np.ndarray = np.asarray([0.0, 0.0, 0.0])
        self.conductance: dict = {}

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
    def calc_geometric_factors(self):
        pass

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
        else:
            raise ValueError('conductivity must be either single scalar,'
                             'an iterable with two (tp, ip) or three (x, '
                             'y, z) entries, or a full array with the same '
                             'shape as the attribute "geometric_factors"')
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
    def connect(cls, conductance_a: Iterable[Self],
                conductance_b: Iterable[Self], mode):
        conductance_a = np.asarray(conductance_a)
        conductance_b = np.asarray(conductance_b)
        if mode == 'serial':
            try:
                return 1.0 / (1.0 / (conductance_a + conductance_b))
            except FloatingPointError:
                raise
        elif mode == 'parallel':
            return conductance_a + conductance_b
        else:
            raise ValueError("argument 'mode' must either be "
                             "'serial' or 'parallel'")

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

    def update(self, transport_properties, *args, **kwargs):
        pass


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

    # def __new__(cls, input_dict: dict, transport_properties: dict,
    #             discretization: dsct.Discretization2D, *args, **kwargs):
    #     instance = super().__new__(cls, input_dict, transport_properties,
    #                                discretization, *args, **kwargs)
    #     return instance

    def __init__(self, input_dict: dict, transport_properties: dict,
                 discretization: dsct.Discretization2D, *args, **kwargs):
        # Initialize super class
        name = input_dict.get('name', 'unnamed')
        self.thickness = input_dict['thickness']
        super().__init__(input_dict, transport_properties, discretization)
        self.geometric_factors = self.calc_geometric_factors()
        self.conductance = {key: self.calc_conductance(value) for key, value
                            in transport_properties.items()}

    def calc_geometric_factors(self):
        result = np.asarray([
            self.discretization.d_area / self.thickness,
            self.discretization.dx[1] / self.discretization.dx[0] * self.thickness,
            self.discretization.dx[0] / self.discretization.dx[1] * self.thickness
            ])
        if self.effective:
            result *= (
                    (1.0 - self.porosity) ** self.bruggeman_exponent)
        return result


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
        name = input_dict.get('name', 'unnamed')

        super().__init__(input_dict, transport_properties, discretization)
        self.geometric_factors = self.calc_geometric_factors()
        self.conductance = {key: self.calc_conductance(value) for key, value
                            in transport_properties.items()}

    def calc_geometric_factors(self):
        result = np.asarray([
            self.discretization.d_area[0] / self.discretization.dx[0],
            self.discretization.d_area[1] / self.discretization.dx[1],
            self.discretization.d_area[2] / self.discretization.dx[2]
            ])
        if self.effective:
            result *= (
                    (1.0 - self.porosity) ** self.bruggeman_exponent)
        return result

    def update(self, transport_properties, *args, **kwargs):
        for key in self.conductance:
            self.conductance[key][:] = self.calc_conductance(
                transport_properties[key])

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
