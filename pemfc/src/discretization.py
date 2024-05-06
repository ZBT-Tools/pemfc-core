import numpy as np
from functools import reduce


class Discretization2D:
    """
    Abstract base class to precalculate discretization of domains
    """
    def __init__(self, discretization_dict, **kwargs):
        super().__init__(**kwargs)

        # Derive discretization shape
        shape = discretization_dict['shape']
        exception_msg = \
            "parameter 'discretization' must be either integer or tuple " \
            "of length 2 or less"
        if isinstance(shape, int):
            self.shape = (shape, 1)
        elif isinstance(shape, (tuple, list)):
            if len(shape) == 1:
                self.shape = shape + (1, )
            elif len(shape) == 2:
                self.shape = shape
            else:
                raise NotImplementedError(exception_msg)
        else:
            raise NotImplementedError(exception_msg)

        # Length: dimension in x-coordinate
        self.length = discretization_dict['length']
        # Width: dimension in y-coordinate
        self.width = discretization_dict['width']

        self.area = self.width * self.length

        self.ratio = discretization_dict.get('ratio', (1.0, 1.0))
        self.direction = \
            discretization_dict.get('direction', (1, 1))

        x_length = [self.length, self.width]
        x = [x_length[i] * self.calculate_spacing(self.shape[i], self.ratio[i],
             self.direction[i]) for i in range(2)]

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


# disc_dict = {
#     'length': 0.5,
#     'width': 0.1,
#     'shape': (10, 2),
#     'ratio': (1.0, 0.5),
#     'direction': (1, -1),
# }
#
# discretization_object = Discretization2D(disc_dict)
