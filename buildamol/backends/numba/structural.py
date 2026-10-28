"""
Numba backend for structural operations.

This module is auto-loaded by the backend dispatcher in buildamol.backends.api.
It provides Numba JIT-compiled implementations of key geometric functions.

Functions in this module override the NumPy defaults when backend="numba" is selected.
"""

from buildamol.backends.numba import (
    rotate_coords,
    rotation_matrix,
    euclidean_distances,
    IC_to_xyz,
)

__all__ = [
    "rotate_coords",
    "rotation_matrix",
    "euclidean_distances",
    "IC_to_xyz",
]
