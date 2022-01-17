import numpy as np


def interpolate_1d(array, add_edge_points=False):
    """
    Linear interpolation in between the given array data. If
    add_edge_points is True, the neighbouring value from the settings array is
    used at the edges and the returned array will larger than the settings array.
    """

    interpolated = np.asarray(array[:-1] + array[1:]) * .5
    if add_edge_points:
        first = np.asarray([array[0]])
        last = np.asarray([array[-1]])
        return np.concatenate((first, interpolated, last), axis=0)
    else:
        return interpolated


def interpolate_along_axis(array, axis, add_edge_points=False):
    """
    Linear interpolation in between the given array data along the given
    axis.
    If add_edge_points is True, the neighbouring value from the settings array is
    used at the edges and the returned array will be larger than the settings
    array.
    """
    if axis == 0:
        return interpolate_1d(array, add_edge_points)
    else:
        a = interpolate_1d(np.moveaxis(array, axis, 0), add_edge_points)
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