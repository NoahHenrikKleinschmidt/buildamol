"""
JAX backend for structural and optimization operations.

This module provides JAX implementations of key geometric and algorithmic functions.
JAX enables functional transformations like jit, vmap, and grad.

Functions in this module override the NumPy defaults when backend="jax" is selected.
"""

import buildamol.utils.auxiliary as aux
from . import algorithms as _algorithms


__all__ = [
    "rotate_coords",
    "rotation_matrix",
    "euclidean_distances",
    "swarm_optimize",
    "anneal_optimize",
    "genetic_optimize",
]

# Re-export algorithms
swarm_optimize = _algorithms.swarm_optimize
anneal_optimize = _algorithms.anneal_optimize
genetic_optimize = _algorithms.genetic_optimize


def rotate_coords(coords, angle, axis):
    """
    JAX-based coordinate rotation.

    Parameters
    ----------
    coords : jax.Array
        Coordinates to rotate (N, 3)
    angle : float
        Rotation angle in radians
    axis : jax.Array
        Rotation axis (3,)

    Returns
    -------
    jax.Array
        Rotated coordinates (N, 3)
    """
    jnp = aux.get_jax_numpy()

    x, y, z = axis
    c = jnp.cos(angle)
    s = jnp.sin(angle)
    t = 1.0 - c

    rot = jnp.array(
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * y * x + s * z, t * y * y + c, t * y * z - s * x],
            [t * z * x - s * y, t * z * y + s * x, t * z * z + c],
        ],
        dtype=coords.dtype,
    )

    rot = jnp.transpose(rot)
    return jnp.dot(coords, rot)


def rotation_matrix(axis, angle):
    """
    JAX-based rotation matrix computation.

    Parameters
    ----------
    axis : jax.Array
        Unit rotation axis (3,)
    angle : float
        Rotation angle in radians

    Returns
    -------
    jax.Array
        3x3 rotation matrix
    """
    jnp = aux.get_jax_numpy()

    x, y, z = axis
    c = jnp.cos(angle)
    s = jnp.sin(angle)
    t = 1.0 - c

    return jnp.array(
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * y * x + s * z, t * y * y + c, t * y * z - s * x],
            [t * z * x - s * y, t * z * y + s * x, t * z * z + c],
        ],
        dtype=axis.dtype,
    )


def euclidean_distances(X, Y):
    """
    JAX-based pairwise Euclidean distances.

    Parameters
    ----------
    X : jax.Array
        (N, 3) coordinate array
    Y : jax.Array
        (M, 3) coordinate array

    Returns
    -------
    jax.Array
        (N, M) distance matrix
    """
    jnp = aux.get_jax_numpy()

    X_sq = jnp.sum(X**2, axis=1, keepdims=True)
    Y_sq = jnp.sum(Y**2, axis=1, keepdims=True).T
    XY = jnp.dot(X, Y.T)

    dist_sq = X_sq + Y_sq - 2 * XY
    dist_sq = jnp.maximum(dist_sq, 0.0)
    return jnp.sqrt(dist_sq)
