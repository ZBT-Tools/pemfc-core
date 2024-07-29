import numpy as np


def interpolate_1d(array, add_edge_points=False, weights=None):
    """
    Linear interpolation in between the given array data. If
    add_edge_points is True, the neighbouring value from the settings array is
    used at the edges and the returned array will larger than the settings array.
    """
    array = np.asarray(array)
    if weights is not None:
        weights = np.asarray(weights)
        if len(weights) != 2 or weights.ndim != array.ndim + 1:
            raise ValueError('weights must have one dimension more than array '
                             'and the highest dimension must equal 2')
        if array.ndim > 1:
            raise ValueError('only 1d array supported for interpolation with '
                             'weights provided')
        interpolated = (array[:-1] * weights[0] + array[1:] * weights[1])

    else:
        interpolated = (array[:-1] + array[1:]) * 0.5

        weights = np.ones((2, len(array) - 1)) * 0.5

    if add_edge_points:
        first = np.asarray([array[0] * weights[0, 0]])
        last = np.asarray([array[-1] * weights[1, -1]])
        return np.concatenate((first, interpolated, last), axis=0)
    else:
        return interpolated


def interpolate_along_axis(array, axis, add_edge_points=False, weights=None):
    """
    Linear interpolation in between the given array data along the given
    axis.
    If add_edge_points is True, the neighbouring value from the settings array is
    used at the edges and the returned array will be larger than the settings
    array.
    """
    if axis == 0:
        return interpolate_1d(array, add_edge_points, weights=weights)
    else:
        a = interpolate_1d(np.moveaxis(array, axis, 0), add_edge_points,
                           weights=weights)
        return np.moveaxis(a, 0, axis)


def to_nodes_1d(elements, endpoints=None):
    """
    Calculates an node 1-d-array from an element 1-d-array,
    uses the [:, 1], [:, -2] entries of the calculated node 1-d-array
    to fill the first als last row of the node 1-d-array.
    """
    nodes = interpolate_1d(elements)
    return np.hstack([elements[0], nodes, elements[-1]])


def elements(elements, axis=0):
    """
    Calculates an node 2-d-array from an element 2-d-array,
    uses the [:, 1], [:, -2] entries of the calculated node 2-d-array
    to fill the first als last row of the node 2-d-array.
    """
    nodes = np.asarray((elements[:, :-1] + elements[:, 1:])) * .5
    return np.hstack([elements[:, [0]], nodes, elements[:, [-1]]])