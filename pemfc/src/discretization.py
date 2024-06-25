import numpy as np
from functools import reduce
from abc import ABC, abstractmethod
import pemfc.src.global_functions as gf


class Discretization(ABC):
    """
    Abstract class to define the general interface for the specific
    implementations of discretization classes in cartesian coordinates
    """
    def __new__(cls, discretization_dict, **kwargs):

        # Derive discretization shape
        shape = discretization_dict['shape']
        exception_msg = \
            "parameter 'discretization' must be either integer or tuple " \
            "of length 2 or less"
        if isinstance(shape, int):
            shape = (shape, 1)
            discretization_dict['shape'] = shape
            return super(Discretization, cls).__new__(Discretization2D)
        elif isinstance(shape, (tuple, list)):
            if len(shape) == 1:
                shape = shape + (1, )
                discretization_dict['shape'] = shape
                return super(Discretization, cls).__new__(Discretization2D)
            elif len(shape) == 2:
                shape = shape
                discretization_dict['shape'] = shape
                return super(Discretization, cls).__new__(Discretization2D)
            elif len(shape) == 3:
                shape = shape
                discretization_dict['shape'] = shape
                return super(Discretization, cls).__new__(Discretization3D)
            else:
                raise NotImplementedError(exception_msg)
        else:
            raise NotImplementedError(exception_msg)

    def __init__(self, discretization_dict, **kwargs):
        # Shape parameter defines the dimensionality and the division of each
        # dimension
        self.shape: tuple = discretization_dict['shape']
        # Vector defining the complete length of the discretization shape
        self.length: np.ndarray
        # Vector containing the element lengths in each direction
        self.dx: np.ndarray
        # Vector containing the local evolution of each axis according to the
        # discretization
        self.x: np.ndarray
        # Array containing the total surface area of each plane
        # (i.e. x-y, y-z, z-x for a three-dimensional discretization)
        self.area: np.ndarray
        # Size ratio between adjacent element lengths in each dimension
        self.ratio: tuple
        # Direction of each axis (limited to (-1, 1))
        self.direction: tuple

    @staticmethod
    def calculate_spacing(elements, ratio, direction=1):
        if ratio == 1.0:
            return np.linspace(0, 1.0, elements + 1)
        else:
            array = np.ones(elements + 1)
            # iterator = range(elements) if direction is -1 else reversed(range(elements))
            ratio = ratio if direction == 1 else 1.0 / ratio
            for i in range(elements):
                array[i+1] = array[i] * ratio
            array -= np.min(array)
            if ratio < 1.0:
                array -= np.max(array)
                array = np.abs(array)
            array *= 1.0 / np.max(array)
            # if ratio < 1.0:
            #     array = np.flip(array)
            return array


class Discretization2D(Discretization):
    """
    Class to precalculate discretization of 2D domains
    """
    def __init__(self, discretization_dict, **kwargs):
        super().__init__(discretization_dict, **kwargs)

        # Length of each axis (x, y)
        self.length = np.asarray([discretization_dict['length'],
                                  discretization_dict['width']])

        self.area = np.asarray([self.length[0] * self.length[1], ])
        self.ratio = discretization_dict.get('ratio', (1.0, 1.0))
        self.direction = \
            discretization_dict.get('direction', (1, 1))

        x = [self.length[i] * self.calculate_spacing(
            self.shape[i], self.ratio[i], self.direction[i]) for i in range(
            len(self.length))]
        # x = [np.linspace(0, self.length, self.shape[0] + 1),
        #      np.linspace(0, self.width, self.shape[1] + 1)]
        self.x = np.asarray(
            [np.asarray([x[0] for j in range(x[1].shape[0])]).transpose(),
             np.asarray([x[1] for j in range(x[0].shape[0])])])
        # for i in range(len(x)):
        #     self.x.append(np.asarray([x[i] for j in range(x[1 - i].shape[0])])
        #                   .reshape(x_shape, order='C'))
        dx = [np.abs(np.diff(self.x[0], axis=0)),
              np.abs(np.diff(self.x[1], axis=1))]
        self.dx = np.asarray([dx[0][:, :self.shape[1]], dx[1][:self.shape[0], :]])
        # Temporary reduce discretization to single dx
        # self.dx = self.dx.mean(axis=(1, 2))

        # dx = np.diff(self.x)
        # dy = np.diff(self.y)

        # self.dx = np.asarray([dx for i in range(self.shape[1])]).transpose()
        # self.dy = np.asarray([dy for i in range(self.shape[0])])
        self.d_area = self.dx[0] * self.dx[1]


class Discretization3D(Discretization):
    """
    Class to precalculate discretization of 3D domains
    """
    def __init__(self, discretization_dict, **kwargs):
        super().__init__(discretization_dict, **kwargs)

        # Length of each axis (x, y, z)
        self.length = np.asarray([discretization_dict['depth'],
                                  discretization_dict['length'],
                                  discretization_dict['width'],
                                  ])

        self.area = np.asarray([self.length[1] * self.length[2],
                                self.length[0] * self.length[2],
                                self.length[0] * self.length[1]])
        self.volume = np.prod(self.length)

        self.ratio = discretization_dict.get('ratio', (1.0, 1.0, 1.0))
        self.direction = \
            discretization_dict.get('direction', (1, 1, 1))

        x = [self.length[i] * self.calculate_spacing(
            self.shape[i], self.ratio[i], self.direction[i]) for i in range(
            len(self.length))]
        # x = [np.linspace(0, self.length, self.shape[0] + 1),
        #      np.linspace(0, self.width, self.shape[1] + 1)]
        self.x = np.asarray([
            np.moveaxis(np.asarray([[x[0] for i in range(x[1].shape[0])]
                                    for j in range(x[2].shape[0])]),
                        (0, 1, 2), (2, 1, 0)),

            np.moveaxis(np.asarray([[x[1] for i in range(x[0].shape[0])]
                                    for j in range(x[2].shape[0])]),
                        (0, 1, 2), (2, 0, 1)),

            np.moveaxis(np.asarray([[x[2] for i in range(x[0].shape[0])]
                                    for j in range(x[1].shape[0])]),
                        (0, 1, 2), (1, 0, 2))])

        # for i in range(len(x)):
        #     self.x.append(np.asarray([x[i] for j in range(x[1 - i].shape[0])])
        #                   .reshape(x_shape, order='C'))

        dx = [np.abs(np.diff(self.x[0], axis=0)),
              np.abs(np.diff(self.x[1], axis=1)),
              np.abs(np.diff(self.x[2], axis=2))]
        self.dx = np.asarray([dx[0][:, :self.shape[1], :self.shape[2]],
                              dx[1][:self.shape[0], :, :self.shape[2]],
                              dx[2][:self.shape[0], :self.shape[1], :]])
        # Temporary reduce discretization to single dx
        # self.dx = self.dx.mean(axis=(1, 2))

        # dx = np.diff(self.x)
        # dy = np.diff(self.y)

        # self.dx = np.asarray([dx for i in range(self.shape[1])]).transpose()
        # self.dy = np.asarray([dy for i in range(self.shape[0])])
        self.d_area = np.asarray([self.dx[1] * self.dx[2],
                                  self.dx[0] * self.dx[2],
                                  self.dx[0] * self.dx[1]])
        self.d_volume = self.dx[0] * self.dx[1] * self.dx[2]


# disc_dict = {
#     'length': 0.5,
#     'width': 0.1,
#     'depth': 0.2,
#     'shape': (2, 10, 3),
#     'ratio': (1, 0.9, 0.6),
#     'direction': (1, 1, 1),
# }
#
# discretization_object = Discretization(disc_dict)
