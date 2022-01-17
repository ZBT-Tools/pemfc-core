import numpy as np
from scipy import linalg as sp_la


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
    n_result = reps*n_vector - (reps-1)*overlap_size
    # n_result = reps * (n_vector - overlap_size) + overlap_size
    result = np.zeros(n_result)
    non_overlap_size = int(n_vector-overlap_size)
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


def build_z_cell_conductance_matrix(cond_vector, n_ele):
    list_mat = []
    for i in range(n_ele):
        list_mat.append(build_1d_conductance_matrix(cond_vector[:, i]))
    return sp_la.block_diag(*list_mat)


def build_x_cell_conductance_matrix(cond_vector, n_ele, n_layer=None):
    if not n_ele > 1:
        raise ValueError('x-conductance matrix can only be built for n_ele > 1')
    if n_layer is None:
        n_layer = len(cond_vector)
    center_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele)])
    center_diag[n_layer:-n_layer] *= 2.0
    center_diag *= -1.0
    off_diag = np.concatenate([cond_vector[:, i] for i in range(n_ele-1)])
    return np.diag(center_diag, k=0) \
        + np.diag(off_diag, k=-n_layer) \
        + np.diag(off_diag, k=n_layer)


def build_cell_conductance_matrix(x_cond_vector, z_cond_vector, n_ele):
    z_cond_mtx = build_z_cell_conductance_matrix(z_cond_vector, n_ele)
    if n_ele > 1:
        x_cond_mtx = build_x_cell_conductance_matrix(x_cond_vector, n_ele)
    else:
        x_cond_mtx = 0.0
    return x_cond_mtx + z_cond_mtx


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
            matrix[mtx_id_0, mtx_id_1] = values[i]
            matrix[mtx_id_0, mtx_id_0] = -values[i]
            matrix[mtx_id_1, mtx_id_1] = -values[i]
            matrix[mtx_id_1, mtx_id_0] = values[i]
        else:
            matrix[mtx_id_0, mtx_id_1] += values[i]
            matrix[mtx_id_0, mtx_id_0] += -values[i]
            matrix[mtx_id_1, mtx_id_1] += -values[i]
            matrix[mtx_id_1, mtx_id_0] += values[i]


def create_index_lists(cells):
    n_cells = len(cells)
    index_list = []
    layer_ids = [[] for _ in range(cells[-1].n_layer)]
    for i, cell in enumerate(cells):
        index_array = \
            (cells[i-1].n_ele * cells[i-1].n_layer) * i \
            + cell.index_array
        index_list.append(index_array.tolist())

    for i in range(n_cells):
        for j in range(cells[i].n_layer):
            layer_ids[j].append(index_list[i][j])
    layer_index_list = []
    for sub_list in layer_ids:
        layer_index_list.append(np.hstack(sub_list))
    return index_list, layer_index_list
