from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from scipy import linalg as sp_la
# import scipy as sp
# from . import cell
# from . import global_functions as g_func

if TYPE_CHECKING:
    from pemfc.src.cell import Cell


def tile_add_overlap(array, n, m=1):
    if n <= 1:
        return array
    else:
        add_on = np.copy(array[:-m])
        add_on[:m] += array[-m:]
        last = np.copy(array)
        last[:m] += array[-m:]
        return np.hstack((array[:-m], np.tile(add_on, n-2), last))


def overlapping_vector(vector, reps, overlap_size):
    n_vector = len(vector)
    n_result = reps * n_vector - (reps - 1) * overlap_size
    # n_result = reps * (n_vector - overlap_size) + overlap_size
    result = np.zeros(n_result)
    non_overlap_size = int(n_vector - overlap_size)
    for i in range(reps):
        start_id = i * non_overlap_size
        end_id = start_id + n_vector
        result[start_id:end_id] += vector
    return result


def block_diag_overlap(block_list, overlap):
    m_sblks = [block.shape[0] for block in block_list]
    n_sblks = [block.shape[1] for block in block_list]
    n_blocks = len(block_list)
    m_ovl = overlap[0]
    n_ovl = overlap[1]
    m_final = np.sum(np.asarray(m_sblks)) - (n_blocks - 1) * m_ovl
    n_final = np.sum(np.asarray(n_sblks)) - (n_blocks - 1) * n_ovl
    block_array = np.zeros((m_final, n_final))
    for i, block in enumerate(block_list):
        x_id_start = i * m_sblks[i - 1] - i * m_ovl
        x_id_end = (i + 1) * m_sblks[i] - i * m_ovl
        y_id_start = i * n_sblks[i - 1] - i * n_ovl
        y_id_end = (i + 1) * n_sblks[i] - i * n_ovl
        block_array[x_id_start:x_id_end, y_id_start:y_id_end] += block[:, :]
    return block_array


def build_1d_conductance_matrix(cond_vector, offset=1):
    n_layer = len(cond_vector) + 1
    center_diag = overlapping_vector(cond_vector, 2, n_layer-2)
    center_diag *= -1.0
    off_diag = np.copy(cond_vector)
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-offset) \
        + np.diag(off_diag, k=offset)


def build_x_cell_conductance_matrix(cond_vector):
    list_mat = []
    for j in range(cond_vector.shape[2]):
        for i in range(cond_vector.shape[1]):
            list_mat.append(build_1d_conductance_matrix(cond_vector[:, i, j]))
    return sp_la.block_diag(*list_mat)


def build_y_cell_conductance_matrix(cond_vector, axis, n_layer=None):
    n_ele = cond_vector.shape[axis]
    if not n_ele > 1:
        raise ValueError('y-conductance matrix can only be built for n_ele > 1')
    if n_layer is None:
        n_layer = len(cond_vector)
    # center_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele)])
    # flatten in column-major order
    center_diag = cond_vector.flatten(order='F')
    off_diag = np.copy(center_diag[:-n_layer])
    center_diag[n_layer:-n_layer] *= 2.0
    center_diag *= -1.0
    # off_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele-1)])
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-n_layer) \
        + np.diag(off_diag, k=n_layer)


def build_z_cell_conductance_matrix(cond_vector, axis):
    if axis == -1:
        axis = len(cond_vector.shape) - 1
    offset = 1
    for i in range(axis):
        offset *= cond_vector.shape[i]

    # center_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele)])
    # flatten in column-major order
    center_diag = cond_vector.flatten(order='F') * 0.5
    # cond_vector_new = cond_vector[]
    # for i in range(cond_vector.shape[axis]):

    off_diag = np.copy(center_diag[:-offset])
    center_diag[offset:-offset] *= 2.0
    center_diag *= -1.0
    # off_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele-1)])
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-offset) \
        + np.diag(off_diag, k=offset)


def calculate_center_diagonal(array, axis):
    if not array.ndim == 3:
        raise ValueError('only three-dimensional arrays supported at the moment')
    if axis not in (0, 1, 2, -1):
        raise ValueError('axis must be 0, 1, 2, or -1')
    shape = list(array.shape)
    shape[axis] += 1
    result = np.zeros(shape)
    if axis == 0:
        result[1:-1, :, :] = array[:-1, :, :] + array[1:, :, :]
        result[0, :, :] = array[0, :, :]
        result[-1, :, :] = array[-1, :, :]
    elif axis == 1:
        result[:, 1:-1, :] = array[:, :-1, :] + array[:, 1:, :]
        result[:, 0, :] = array[:, 0, :]
        result[:, -1, :] = array[:, -1, :]
    elif axis == 2 or axis == -1:
        result[:, :, 1:-1] = array[:, :, :-1] + array[:, :, 1:]
        result[:, :, 0] = array[:, :, 0]
        result[:, :, -1] = array[:, :, -1]
    else:
        raise ValueError('only three-dimensional arrays supported')
    return result.flatten(order='F')


def calculate_off_diagonal(array, axis):
    if not array.ndim == 3:
        raise ValueError('only three-dimensional arrays supported at the moment')
    if axis not in (0, 1, 2, -1):
        raise ValueError('axis must be 0, 1, 2, or -1')
    shape = list(array.shape)
    shape[axis] += 1
    result = np.zeros(shape)
    if axis == 0:
        result[:-1, :, :] = array
    elif axis == 1:
        result[:, :-1, :] = array
    elif axis == 2:
        result[:, :, :-1] = array
    else:
        raise ValueError('only three-dimensional arrays supported at the moment')
    return result.flatten(order='F')


def create_one_dimensional_conductance_matrix(conductance_array, axis):
    """
    Calculate a tri-diagonal matrix for conduction in one dimension along the
    given axis, while the order of arranging from inner to outer axis in the
    matrix structure is from lowest to highest axis index 0, 1, 2 i.e. x, y, z
    The conductance array for the provided axis should already give the
    inter-nodal conductances so that the array is one element smaller compared
    to the final solution vector. This should not be true for the remaining
    axis.

    Args:
        conductance_array: k, l, m - sized conductance array while the axis index
                           indicates the connected axis with a reduced dimension
        axis: axis for one-dimensional conduction

    Returns: tri-diagonal matrix for conduction with off-diagonal position
             according to the conduction axis
    """
    if axis == -1:
        axis = len(conductance_array.shape) - 1
    offset = 1
    for i in range(axis):
        offset *= conductance_array.shape[i]
    center_diag = calculate_center_diagonal(conductance_array, axis=axis)
    off_diag = calculate_off_diagonal(conductance_array, axis)[:-offset]
    center_diag *= -1.0
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-offset) \
        + np.diag(off_diag, k=offset)


def build_cell_conductance_matrix(conductance_list: list[np.ndarray]):
    cond_mtx_list = []
    for i, conductance_array in enumerate(conductance_list):
        if (isinstance(conductance_array, np.ndarray)
                and conductance_array.shape[i] > 0):
            cond_mtx_list.append(create_one_dimensional_conductance_matrix(
                conductance_array, axis=i))
    return np.sum(cond_mtx_list, axis=0)


def transform_nodes(node_values: (list, np.ndarray), axis,
                    include_axis=False, inverse=False,
                    mode='shift', except_first_axis=True, **kwargs):
    if len(node_values) == node_values[0].ndim and except_first_axis:
        node_values = [np.asarray(item) for item in node_values]
        if include_axis:
            transform_axes = tuple(range(len(node_values)))
        else:
            transform_axes = tuple(i for i in range(len(node_values))
                                   if i != axis)
        if axis == -1:
            axis = node_values[0].ndim - 1
        for i in transform_axes:
            if inverse and i == axis:
                inverse = True
            else:
                inverse = False
            if mode == 'shift':
                new_values = shift_array_nodes(node_values[i], axis,
                                               inverse=inverse)
            elif mode == 'average':
                new_values = average_array_nodes(node_values[i], axis)
            else:
                raise ValueError('only mode "shift" or "average" is allowed')
            node_values[i] = new_values
        return node_values
    else:
        if mode == 'shift':
            result = shift_array_nodes(np.asarray(node_values), axis,
                                       inverse=inverse)
        elif mode == 'average':
            result = average_array_nodes(np.asarray(node_values), axis)
        else:
            raise ValueError('only mode "shift" or "average" is allowed')
        return result


def average_array_nodes(array: np.ndarray, axis: int):
    old_shape = array.shape
    new_shape = tuple(size - 1 if j == axis else size
                      for j, size in enumerate(old_shape))
    new_array = np.zeros(new_shape)
    # factor = 2.0 if inverse_factor else 0.5
    half_values = 0.5 * np.moveaxis(array, axis, 0)
    np.moveaxis(new_array, axis, 0)[:] += half_values[:-1]
    np.moveaxis(new_array, axis, 0)[:] += half_values[1:]
    return new_array


def shift_array_nodes(array: np.ndarray, axis: int, inverse=False):
    old_shape = array.shape
    new_shape = tuple(size + 1 if j == axis else size
                      for j, size in enumerate(old_shape))
    new_array = np.zeros(new_shape)
    # factor = 2.0 if inverse_factor else 0.5
    if inverse:
        inf_array = np.ones(array.shape) * 1e16
        inverse_array = np.divide(1.0, array, out=inf_array,
                                  where=array != 0.0)
        inverse_half_values = 0.5 * inverse_array

        inverse_result = np.zeros(np.moveaxis(new_array, axis, 0).shape)
        # inf_array = np.ones(half_values.shape) * 1e16
        # inverse_half_values = np.divide(1.0, half_values, out=inf_array,
        #                                 where=half_values != 0.0)
        inverse_result[:-1] += inverse_half_values
        inverse_result[1:] += inverse_half_values
        result = 1.0 / inverse_result
        np.moveaxis(new_array, axis, 0)[:] = result
    else:
        half_values = 0.5 * np.moveaxis(array, axis, 0)
        np.moveaxis(new_array, axis, 0)[:-1] += half_values
        np.moveaxis(new_array, axis, 0)[1:] += half_values
    return new_array


def shift_nodes(node_values: (list, np.ndarray), axis,
                include_axis=False, inverse=False,
                except_first_axis=True, **kwargs):
    return transform_nodes(node_values, axis, include_axis=include_axis,
                           inverse=inverse, mode='shift',
                           except_first_axis=except_first_axis, **kwargs)


def connect_cells(matrix, cell_ids, layer_ids, values, mtx_ids,
                  replace=False):
    if np.isscalar(values):
        values = np.full(len(cell_ids), values)
    if not len(cell_ids) == len(layer_ids):
        raise ValueError('cell and layer index lists must have equal length')
    for i in range(len(cell_ids)):
        mtx_id_0 = mtx_ids[cell_ids[i, 0]][:][layer_ids[i, 0]]
        mtx_id_1 = mtx_ids[cell_ids[i, 1]][:][layer_ids[i, 1]]
        if replace:
            matrix[mtx_id_0, mtx_id_1] = values[i].flatten('F')
            matrix[mtx_id_0, mtx_id_0] = -values[i].flatten('F')
            matrix[mtx_id_1, mtx_id_1] = -values[i].flatten('F')
            matrix[mtx_id_1, mtx_id_0] = values[i].flatten('F')
        else:
            matrix[mtx_id_0, mtx_id_1] += values[i].flatten('F')
            matrix[mtx_id_0, mtx_id_0] += -values[i].flatten('F')
            matrix[mtx_id_1, mtx_id_1] += -values[i].flatten('F')
            matrix[mtx_id_1, mtx_id_0] += values[i].flatten('F')


def set_single_unit_entry(matrix, index):
    row = np.zeros(matrix.shape[0])
    row[index] = 1.0
    matrix[index, :] = row
    return matrix


def get_single_axis_values(array, axis: int,
                           indices: (int, tuple, np.ndarray) = None):
    if indices is None:
        indices = np.arange(array.shape[axis])
    return np.take(array, indices, axis=axis)


def get_axis_values(array, axes: tuple, indices: tuple):
    if isinstance(axes, int):
        axes = (axes,)
    if isinstance(indices, int):
        indices = (indices,)
    for i in reversed(range(len(axes))):
        array = get_single_axis_values(array, axes[i], indices[i])
    return array


def set_axis_values(array: np.ndarray, values: np.ndarray,
                    axes: tuple, indices: tuple):
    if isinstance(axes, int):
        axes = (axes,)
    if isinstance(indices, int):
        indices = (indices,)
    shape = array.shape
    ids_flat = np.arange(int(np.prod(shape)))
    indices = get_axis_values(np.reshape(ids_flat, shape, order='F'), axes,
                              indices).ravel(order='F')
    array = array.ravel(order='F')
    values = values.ravel(order='F')
    np.put(array, indices, values)
    result = array.reshape(shape, order='F')
    return result


def spai(matrix: np.ndarray, m: int):
    """Perform m step of the SPAI iteration."""
    from scipy.sparse import identity
    from scipy.sparse import diags
    from scipy.sparse.linalg import onenormest

    n = matrix.shape[0]

    ident = identity(n, format='csr')
    alpha = 2 / onenormest(matrix @ matrix.T)
    M = alpha * matrix

    for index in range(m):
        C = matrix @ M
        G = ident - C
        AG = matrix @ G
        trace = (G.T @ AG).diagonal().sum()
        alpha = trace / np.linalg.norm(AG.data) ** 2
        M = M + alpha * G

    return M
