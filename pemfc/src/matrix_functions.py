from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from scipy import linalg as sp_la
import scipy as sp
from . import cell
from . import global_functions as g_func

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


def calculate_center_diagonal(array, axis, weights=(0.5, 1.0, 0.5)):
    if not array.ndim == 3:
        raise ValueError('only three-dimensional arrays supported at the moment')
    result = np.zeros(array.shape)
    if axis == 0:
        result[1:-1, :, :] = \
            sp.ndimage.convolve1d(array, weights, mode='constant',
                                  axis=axis, cval=0.0)[1:-1, :, :]
        result[0, :, :] = (array[0, :, :] + array[1, :, :]) * 0.5
        result[-1, :, :] = (array[-2, :, :] + array[-1, :, :]) * 0.5

    elif axis == 1:
        result[:, 1:-1, :] = \
            sp.ndimage.convolve1d(array, weights, mode='constant',
                                  axis=axis, cval=0.0)[:, 1:-1, :]
        result[:, 0, :] = (array[:, 0, :] + array[:, 1, :]) * 0.5
        result[:, -1, :] = (array[:, -2, :] + array[:, -1, :]) * 0.5

    elif axis == 2 or axis == -1:
        result[:, :, 1:-1] = \
            sp.ndimage.convolve1d(array, weights, mode='constant',
                                  axis=axis, cval=0.0)[:, :, 1:-1]
        result[:, :, 0] = (array[:, :, 0] + array[:, :, 1]) * 0.5
        result[:, :, -1] = (array[:, :, -2] + array[:, :, -1]) * 0.5
    else:
        raise ValueError('only three-dimensional arrays supported')
    return result.flatten(order='F')


def calculate_off_diagonal(array, axis, weights=(0.5, 0.5)):
    off_coeff_matrix = sp.ndimage.convolve1d(
        array, weights, mode='constant', axis=axis, cval=0.0)
    if axis == 0:
        off_coeff_matrix[-1, :, :] = 0.0
    elif axis == 1:
        off_coeff_matrix[:, -1, :] = 0.0
    elif axis == 2:
        off_coeff_matrix[:, :, -1] = 0.0
    else:
        raise ValueError('only three-dimensional arrays supported at the moment')

    return off_coeff_matrix.flatten(order='F')


def build_one_dimensional_conductance_matrix(conductance_array, axis,
                                             center_weights=(0.5, 1.0, 0.5),
                                             off_weights=(0.5, 0.5)):
    if axis == -1:
        axis = len(conductance_array.shape) - 1
    offset = 1
    for i in range(axis):
        offset *= conductance_array.shape[i]
    center_diag = calculate_center_diagonal(
        conductance_array, axis=axis, weights=center_weights)
    off_diag = calculate_off_diagonal(
        conductance_array, axis, weights=off_weights)[:-offset]
    center_diag *= -1.0
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-offset) \
        + np.diag(off_diag, k=offset)


def calculate_center_diagonal_2(array, axis):
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


def calculate_off_diagonal_2(array, axis):
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
    center_diag = calculate_center_diagonal_2(conductance_array, axis=axis)
    off_diag = calculate_off_diagonal_2(conductance_array, axis)[:-offset]
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


def add_explicit_layer_source(rhs_vector, source_term, index_array,
                              layer_id=None, replace=False):
    if isinstance(source_term, np.ndarray):
        if rhs_vector.ndim != source_term.ndim:
            raise ValueError('source_term must have same dimensions as '
                             'rhs_vector')
    if layer_id is None:
        if np.isscalar(source_term):
            source_vector = np.full_like(rhs_vector, -source_term)
        else:
            source_vector = np.asarray(-source_term)
    else:
        if replace is True:
            source_vector = np.copy(rhs_vector)
            np.put(source_vector, index_array[layer_id], -source_term)
        else:
            source_vector = np.zeros(rhs_vector.shape)
            np.put(source_vector, index_array[layer_id], -source_term)
    if replace is True:
        rhs_vector[:] = source_vector
    else:
        rhs_vector += source_vector
    return rhs_vector, source_vector


def add_implicit_layer_source(matrix, coefficients, index_array,
                              layer_id=None, replace=False):
    # matrix_size = matrix.shape[0]
    diag_vector = np.diagonal(matrix).copy()
    # if layer_id is None:
    #     if np.isscalar(coefficients):
    #         source_vector = np.full_like(diag_vector, coefficients)
    #     else:
    #         source_vector = np.asarray(coefficients)
    # else:
    #     if replace is True:
    #         source_vector = np.copy(diag_vector)
    #         np.put(source_vector, index_array[layer_id], coefficients)
    #     else:
    #         source_vector = np.zeros(matrix_size)
    #         np.put(source_vector, index_array[layer_id], coefficients)
    # if replace is True:
    #     diag_vector = source_vector
    # else:
    #     diag_vector += source_vector
    new_diag_vector, source_vector = add_explicit_layer_source(
        diag_vector, -coefficients.flatten(order='F'), index_array,
        layer_id=layer_id, replace=replace)
    np.fill_diagonal(matrix, new_diag_vector)
    return matrix, source_vector


def set_implicit_layer_fixed(matrix, index_array, layer_id):
    row_length = matrix.shape[0]
    for row_id in index_array[layer_id]:
        row = np.zeros(row_length)
        row[row_id] = 1.0
        matrix[row_id] = row
    return matrix


def set_single_unit_entry(matrix, index):
    row = np.zeros(matrix.shape[0])
    row[index] = 1.0
    matrix[index, :] = row
    return matrix
