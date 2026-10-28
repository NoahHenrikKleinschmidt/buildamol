"""
Numba-accelerated backend implementations for structural and optimization operations.

This module provides Numba JIT-compiled versions of key geometric,
distance computation, and optimization functions.
"""

import numpy as np
import buildamol.utils.auxiliary as aux
from . import algorithms as _algorithms


__all__ = [
    "rotate_coords",
    "rotation_matrix",
    "euclidean_distances",
    "IC_to_xyz",
    "swarm_optimize",
    "anneal_optimize",
    "genetic_optimize",
]

# Re-export algorithms
swarm_optimize = _algorithms.swarm_optimize
anneal_optimize = _algorithms.anneal_optimize
genetic_optimize = _algorithms.genetic_optimize


@aux.njit
def _numba_rotation_matrix(axis, angle):
    """
    Compute a 3x3 rotation matrix using Rodrigues' rotation formula.

    Parameters
    ----------
    axis : np.ndarray
        Unit rotation axis (3,)
    angle : float
        Rotation angle in radians

    Returns
    -------
    np.ndarray
        3x3 rotation matrix
    """
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1.0 - c

    return np.array(
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * y * x + s * z, t * y * y + c, t * y * z - s * x],
            [t * z * x - s * y, t * z * y + s * x, t * z * z + c],
        ]
    )


@aux.njit
def _numba_rotate_coords(coords, angle, axis):
    """
    Rotate coordinates around an axis by an angle.

    Parameters
    ----------
    coords : np.ndarray
        Coordinates to rotate (N, 3)
    angle : float
        Rotation angle in radians
    axis : np.ndarray
        Rotation axis (3,)

    Returns
    -------
    np.ndarray
        Rotated coordinates (N, 3)
    """
    rot = _numba_rotation_matrix(axis, angle)
    rot = np.transpose(rot)
    return np.dot(coords, rot)


def rotate_coords(coords, angle, axis):
    """
    Numba-accelerated coordinate rotation.

    Parameters
    ----------
    coords : np.ndarray
        Coordinates to rotate (N, 3)
    angle : float
        Rotation angle in radians
    axis : np.ndarray
        Rotation axis (3,)

    Returns
    -------
    np.ndarray
        Rotated coordinates (N, 3)
    """
    return _numba_rotate_coords(coords, angle, axis)


def rotation_matrix(axis, angle):
    """
    Numba-accelerated rotation matrix computation.

    Parameters
    ----------
    axis : np.ndarray
        Unit rotation axis (3,)
    angle : float
        Rotation angle in radians

    Returns
    -------
    np.ndarray
        3x3 rotation matrix
    """
    return _numba_rotation_matrix(axis, angle)


@aux.njit
def _numba_euclidean_distances(X, Y):
    """
    Compute pairwise Euclidean distances.

    Parameters
    ----------
    X : np.ndarray
        (N, 3) coordinate array
    Y : np.ndarray
        (M, 3) coordinate array

    Returns
    -------
    np.ndarray
        (N, M) distance matrix
    """
    n, m = X.shape[0], Y.shape[0]
    dists = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            dx = X[i, 0] - Y[j, 0]
            dy = X[i, 1] - Y[j, 1]
            dz = X[i, 2] - Y[j, 2]
            dists[i, j] = np.sqrt(dx * dx + dy * dy + dz * dz)
    return dists


def euclidean_distances(X, Y):
    """
    Numba-accelerated pairwise Euclidean distances.

    Parameters
    ----------
    X : np.ndarray
        (N, 3) coordinate array
    Y : np.ndarray
        (M, 3) coordinate array

    Returns
    -------
    np.ndarray
        (N, M) distance matrix
    """
    return _numba_euclidean_distances(X, Y)


@aux.njit
def _numba_IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    Numba-accelerated internal coordinate to Cartesian coordinate conversion.

    Parameters
    ----------
    a, b, c : np.ndarray
        Reference atoms (3,)
    anchor : np.ndarray
        Anchor atom (3,)
    r : float
        Bond length
    theta : float
        Bond angle
    dihedral : float
        Dihedral angle

    Returns
    -------
    np.ndarray
        New atom coordinate (3,)
    """
    ab = b - a
    bc = c - b

    plane_ABC = np.cross(ab, bc)
    plane_ABC = plane_ABC / np.linalg.norm(plane_ABC)

    plane_bcd = np.cross(bc, np.array([0.0, 0.0, 1.0]))
    plane_bcd = plane_bcd / np.linalg.norm(plane_bcd)

    bc_normalized = bc / np.linalg.norm(bc)

    _rot = _numba_rotation_matrix(plane_bcd, theta)
    _rot = np.transpose(_rot)
    direction = np.dot(np.array([r, 0.0, 0.0]), _rot)

    _rot2 = _numba_rotation_matrix(bc_normalized, dihedral)
    _rot2 = np.transpose(_rot2)
    direction = np.dot(direction, _rot2)

    return anchor + direction


def IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    Numba-accelerated internal coordinate to Cartesian conversion.

    Parameters
    ----------
    a, b, c : np.ndarray
        Reference atoms (3,)
    anchor : np.ndarray
        Anchor atom (3,)
    r : float
        Bond length
    theta : float
        Bond angle
    dihedral : float
        Dihedral angle

    Returns
    -------
    np.ndarray
        New atom coordinate (3,)
    """
    return _numba_IC_to_xyz(a, b, c, anchor, r, theta, dihedral)
