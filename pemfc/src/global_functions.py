# general imports
import numpy as np
# from scipy import ndimage

# local module imports
from . import constants


def ensure_list(variable, length=1):
    if isinstance(variable, (list, tuple, np.ndarray)):
        return variable
    else:
        return [variable for i in range(length)]


def dim(a):
    if not isinstance(a, (list, tuple)):
        return []
    return [len(a)] + dim(a[0])


def full_like(array):
    """faster than native numpy version"""
    result = np.zeros(array.shape)
    result[:] = array
    return result


def full(shape, value):
    """faster than native numpy version"""
    result = np.zeros(shape)
    result[:] = value
    return result


def zeros_like(array):
    """faster than native numpy version"""
    if isinstance(array, np.ndarray):
        return np.zeros(array.shape)
    else:
        return 0.0


def fill_transposed(in_array, shape):
    transposed_array = np.zeros(shape).transpose()
    transposed_array[:] = in_array
    return transposed_array.transpose()


def move_axis(array, source=0, target=-1):
    if array.ndim > 2:
        return np.moveaxis(array, source, target)
    else:
        return array.transpose()


def array_vector_multiply(array, vector):
    return (array.transpose() * vector).transpose()


def add_source(var, source, direction=1, tri_mtx=None):
    """
    Add discrete 1d source of length n-1 to var of length n
    :param var: 1d array of quantity variable
    :param source: 1d array of source to add to var
    :param direction: flow direction (1: along array counter, -1: opposite to
    array counter)
    :param tri_mtx: if triangle matrix (2D array, nxn) is not provided,
    it will be created temporarily
    :return:
    """
    n = len(var) - 1
    if len(source) != n:
        raise ValueError('parameter source must be of length (var-1)')
    if direction == 1:
        if tri_mtx is None:
            ones = np.zeros((n, n))
            ones.fill(1.0)
            fwd_mat = np.tril(ones)
        else:
            fwd_mat = tri_mtx
        var[1:] += np.matmul(fwd_mat, source)
    elif direction == -1:
        if tri_mtx is None:
            ones = np.zeros((n, n))
            ones.fill(1.0)
            bwd_mat = np.triu(ones)
        else:
            bwd_mat = tri_mtx
        var[:-1] += np.matmul(bwd_mat, source)
    else:
        raise ValueError('parameter direction must be either 1 or -1')
    return var


def np_log(array):
    """
    Use numpy log function and set zero values where solution is not defined
    """
    return np.log(array, out=zeros_like(array), where=(array != 0))


def np_div(array1, array2):
    """
    Use numpy divide function and set zero values where solution is not defined
    """
    if isinstance(array1, np.ndarray):
        return np.divide(array1, array2, out=zeros_like(array1),
                         where=(array2 != 0))
    elif isinstance(array2, np.ndarray):
        return np.divide(array1, array2, out=zeros_like(array2),
                         where=(array2 != 0))
    else:
        if array2 == 0.0:
            return 0.0
        else:
            return array1 / array2


def exponential_distribution(y_avg, nx, a=1.0, b=0.0):
    n_nodes = nx + 1
    x = np.linspace(0, 1, n_nodes)
    dx = np.diff(x)
    m = a * (y_avg - b) / (1.0 - np.exp(-a))
    return m / (a * dx) * (np.exp(-a * x[:-1]) - np.exp(-a * x[1:])) + b


def fill_last_zeros(array, axis=-1, axis_sum=None):
    if axis == 0:
        array_t = array.transpose()
        return fill_last_zeros(array_t, axis=-1).transpose()
    if axis_sum is None:
        axis_sum = np.abs(np.sum(array, axis=0))
    shape = array.shape
    prev = np.arange(shape[-1])
    prev[axis_sum < constants.SMALL] = 0
    prev = np.maximum.accumulate(prev)
    return array[:, prev]


def fill_first_zeros(array, axis=-1, axis_sum=None):
    array = np.flip(array, axis)
    return np.flip(fill_last_zeros(array, axis), axis, axis_sum)


def fill_zero_sum(array, axis=-1, axis_sum=None):
    if axis == 0:
        array_t = array.transpose()
        return fill_zero_sum(array_t, axis=-1).transpose()
    elif axis in (-1, 1):
        if axis_sum is None:
            try:
                axis_sum = np.abs(np.sum(array, axis=0))
            except RecursionError:
                print('RecursionError in fill_zero_sum reached')
        else:
            axis_sum = np.abs(axis_sum)
    else:
        raise ValueError('axis must be 0, 1 or -1. only 2D-arrays allowed.')
    nonzero = np.nonzero(axis_sum)[0]
    try:
        if nonzero[-1] != array.shape[axis] - 1:
            array = fill_last_zeros(array, axis_sum=axis_sum)
        if nonzero[0] != 0:
            array = fill_first_zeros(array, axis_sum=axis_sum)
    except IndexError:
        raise IndexError
    return array


# def fill_surrounding_average_1d(array, axis=0):
#     footprint = np.zeros((3, 3))
#     weights = np.array([1.0, 0.0, 1.0])
#     if axis == 0:
#         footprint[1, :] = weights
#     elif axis in (-1, 1):
#         footprint[:, 1] = weights
#     else:
#         raise ValueError('argument axis can only be 0, 1 or -1')
#
#     mask_array = np.sum(np.abs(array), axis) * np.ones(array.shape)
#     averaged = ndimage.generic_filter(array, np.nanmean, footprint=footprint,
#                                       mode='constant', cval=np.NaN)
#     return np.where(mask_array < constants.SMALL, averaged, array)


def construct_empty_stack_array(cell_array, n_cells):
    """
    Construct zeroed stack array from one- or multidimensional cell variable
    array and number of cells
    :param cell_array: array of variable discretized for a unit cell
    :param n_cells: number of cells in stack
    :return:
    """
    if isinstance(cell_array, np.ndarray):
        cell_array_shape = cell_array.shape
    else:
        cell_array_shape = (len(cell_array))
    stack_array_shape = (n_cells,) + cell_array_shape
    return np.zeros(stack_array_shape)


def calc_temp_heat_transfer(wall_temp, fluid_temp, capacity_rate, heat_coeff,
                            flow_direction):
    wall_temp = np.asarray(wall_temp)
    fluid_temp = np.asarray(fluid_temp)
    fluid_temp_old = np.copy(fluid_temp)
    capacity_rate = np.asarray(capacity_rate)
    heat_coeff = np.asarray(heat_coeff)
    assert capacity_rate.shape == wall_temp.shape
    assert heat_coeff.shape == wall_temp.shape
    fluid_temp_avg = np.asarray(fluid_temp[:-1] + fluid_temp[1:]) * .5
    # fluid_temp_avg = np.zeros(wall_temp.shape)
    # if flow_direction == 1:
    #     fluid_temp_avg[:] = fluid_temp[0]
    # else:
    #     fluid_temp_avg[:] = fluid_temp[-1]
    id_range = range(len(wall_temp))
    if flow_direction == -1:
        id_range = reversed(id_range)
    for i in id_range:
        fluid_avg = fluid_temp_avg[i]
        fluid_out_old = 5e5
        error = 1e3
        iter = 0
        itermax = 10
        while error > 1e-4 and iter <= itermax:
            if flow_direction == -1:
                fluid_in = fluid_temp[i + 1]
            else:
                fluid_in = fluid_temp[i]
            if wall_temp[i] == fluid_in:
                delta_temp = wall_temp[i] - fluid_avg
            else:
                temp_diff_ratio = (wall_temp[i] - fluid_avg) \
                    / (wall_temp[i] - fluid_in)
                if temp_diff_ratio > 0.0:
                    delta_temp = wall_temp[i] - fluid_avg
                else:
                    delta_temp = wall_temp[i] - fluid_in
            q = heat_coeff[i] * delta_temp
            fluid_out = fluid_in + q / capacity_rate[i]
            if fluid_in < wall_temp[i]:
                fluid_out = np.minimum(wall_temp[i] - 1e-3, fluid_out)
            else:
                fluid_out = np.maximum(wall_temp[i] + 1e-3, fluid_out)
            fluid_avg = (fluid_in + fluid_out) * 0.5
            error = np.abs(fluid_out_old - fluid_out) / fluid_out
            fluid_out_old = np.copy(fluid_out)
            iter += 1
        if flow_direction == -1:
            fluid_temp[i] = fluid_out
        else:
            fluid_temp[i + 1] = fluid_out
    if np.any(fluid_temp < 200.0):
        raise ValueError
    fluid_temp_avg = np.asarray(fluid_temp[:-1] + fluid_temp[1:]) * .5
    heat = heat_coeff * (wall_temp - fluid_temp_avg)
    return fluid_temp, heat


def calc_diff(vec):
    """
    Calculates the difference between the i+1 and i position of an 1-d-array.
    """
    return vec[:-1] - vec[1:]


def sub(char):
    """
    return subscript string of inputs string
    :param char: characters to be subscript
    :return: characters as subscript
    """
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
             "abcdefghijklmnopqrstuvwxyz0123456789+-=()"
    sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZ" \
            "ₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
    res = char.maketrans(''.join(normal), ''.join(sub_s))
    return char.translate(res)


def sup(char):
    """
    return subscript string of inputs string
    :param char: characters to be subscript
    :return: characters as subscript
    """
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" \
             "abcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻ" \
              "ᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = char.maketrans(''.join(normal), ''.join(super_s))
    return char.translate(res)



