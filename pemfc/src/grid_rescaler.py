"""Precomputed axis-wise grid rescaler for fixed source/target shapes.

Equivalent in numerical behaviour to repeated :func:`global_functions.rescale`
calls with the same shapes, but all Python-level loops and linspace
constructions are done once in ``__init__``. The :meth:`apply` method then
uses pure-numpy indexing to produce the rescaled array without any Python
iteration over 1D slices.

This module is intended as a drop-in replacement for the
``gf.rescale``-based ``DiffusionTransport.rescale_input`` path. The original
path remains available via the ``USE_VECTORIZED_RESCALE`` flag in
:mod:`diffusion_transport` — set it to ``False`` to fall back to
``gf.rescale`` without touching any other code.
"""
import numpy as np


class GridRescaler:
    """Rescale arrays of a fixed source shape to a fixed target shape.

    Parameters
    ----------
    source_shape, target_shape : tuple of int
        Shape of the arrays passed to :meth:`apply` and the shape returned.
    axes : tuple of int, optional
        Axes to rescale. Defaults to all axes of ``target_shape``.
    method : str or sequence of str
        Interpolation method per axis. Accepted values: ``'linear'``
        (cell-centered linear interpolation, same grid convention as
        :func:`global_functions.linear_rescale_1d`), ``'constant_segments'``
        (replicated segments, same as :func:`global_functions.segment_upscale`).
    """

    def __init__(self, source_shape, target_shape, axes=None,
                 method='linear'):
        self.source_shape = tuple(source_shape)
        self.target_shape = tuple(target_shape)
        if axes is None:
            axes = tuple(range(len(self.target_shape)))
        self.axes = tuple(axes)
        methods = ([method] * len(self.axes)
                   if isinstance(method, str)
                   else list(method))

        # Per-axis application functions. None entries = axis unchanged.
        self._steps = []
        current_shape = list(self.source_shape)
        for ax, m in zip(self.axes, methods):
            src_len = current_shape[ax]
            tgt_len = self.target_shape[ax]
            if src_len == tgt_len:
                self._steps.append(None)
            else:
                self._steps.append(
                    _build_axis_rescaler(src_len, tgt_len, ax, m))
            current_shape[ax] = tgt_len

    def apply(self, array: np.ndarray) -> np.ndarray:
        if array.shape != self.source_shape:
            raise ValueError(
                f"input shape {array.shape} does not match expected "
                f"source shape {self.source_shape}")
        for step in self._steps:
            if step is not None:
                array = step(array)
        return array


# ---------------------------------------------------------------------------
# Per-axis rescaler factories
# ---------------------------------------------------------------------------

def _cell_centers(n: int) -> np.ndarray:
    """Mid-points of n equal intervals in [0, 1].

    Matches the grid convention used by
    :func:`global_functions.linear_rescale_1d`.
    """
    edges = np.linspace(0.0, 1.0, n + 1)
    return (edges[:-1] + edges[1:]) * 0.5


def _build_axis_rescaler(src_len: int, tgt_len: int, axis: int, method: str):
    """Return a callable that rescales a single axis from ``src_len`` to
    ``tgt_len`` using the given ``method``."""
    if method == 'linear':
        return _build_linear_interp(src_len, tgt_len, axis)
    elif method == 'constant_segments':
        return _build_segment_upscale(src_len, tgt_len, axis)
    else:
        raise NotImplementedError(
            f"method {method!r} not supported — use 'linear' or "
            f"'constant_segments'")


def _build_linear_interp(src_len: int, tgt_len: int, axis: int):
    """Vectorised equivalent of ``apply_along_axis(linear_rescale_1d, ...)``.

    Pre-computes the low/high source-index pair and the interpolation weight
    for each target index, then applies them as broadcast numpy operations —
    no Python loop over 1D slices.
    """
    x_old = _cell_centers(src_len)
    x_new = _cell_centers(tgt_len)
    idx_high = np.clip(np.searchsorted(x_old, x_new), 1, src_len - 1)
    idx_low = idx_high - 1
    denom = x_old[idx_high] - x_old[idx_low]
    alpha = np.divide(x_new - x_old[idx_low], denom,
                      out=np.zeros_like(x_new), where=denom > 0.0)
    # Clamp to [0, 1] so out-of-range x_new values map to the nearest
    # source value, matching np.interp's edge behaviour.
    alpha = np.clip(alpha, 0.0, 1.0)

    def apply(arr: np.ndarray) -> np.ndarray:
        low = np.take(arr, idx_low, axis=axis)
        high = np.take(arr, idx_high, axis=axis)
        weight_shape = [1] * arr.ndim
        weight_shape[axis] = -1
        a = alpha.reshape(weight_shape)
        return (1.0 - a) * low + a * high

    return apply


def _build_segment_upscale(src_len: int, tgt_len: int, axis: int):
    """Vectorised equivalent of :func:`global_functions.segment_upscale`."""
    segment_size = tgt_len // src_len
    if segment_size == 0:
        raise ValueError(
            f"target length {tgt_len} must be at least source length "
            f"{src_len} for 'constant_segments' method")
    indices = np.arange(tgt_len) // segment_size
    # segment_upscale pads any remainder with the last source value
    indices = np.minimum(indices, src_len - 1)

    def apply(arr: np.ndarray) -> np.ndarray:
        return np.take(arr, indices, axis=axis)

    return apply
