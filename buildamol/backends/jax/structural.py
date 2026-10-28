"""
JAX backend for structural operations.

This module is auto-loaded by the backend dispatcher in buildamol.backends.api.
It provides JAX implementations of key geometric functions.

Functions in this module override the NumPy defaults when backend="jax" is selected.
"""

from buildamol.backends.jax import (
    rotate_coords,
    rotation_matrix,
    euclidean_distances,
)

__all__ = [
    "rotate_coords",
    "rotation_matrix",
    "euclidean_distances",
]
