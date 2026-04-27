"""
Tests for backend selection and dispatch in structural operations.
"""

import numpy as np
import pytest
import buildamol as bam


# ==================================================================================
# Backend Selection Tests
# ==================================================================================


def test_backend_resolution_numpy():
    """Test backend resolution defaults to NumPy."""
    resolved = bam.utils.auxiliary.resolve_compute_backend(None)
    assert resolved == "numpy"


def test_backend_resolution_explicit():
    """Test explicit backend resolution."""
    assert bam.utils.auxiliary.resolve_compute_backend("numpy") == "numpy"
    assert bam.utils.auxiliary.resolve_compute_backend("NUMPY") == "numpy"


def test_backend_resolution_invalid():
    """Test that invalid backends raise ValueError."""
    with pytest.raises(ValueError):
        bam.utils.auxiliary.resolve_compute_backend("invalid")


def test_backend_resolution_numba_optional():
    """Test that unavailable backends fall back to NumPy."""
    # This should not raise, but return numpy if numba isn't available
    resolved = bam.utils.auxiliary.resolve_compute_backend("numba")
    assert resolved in {"numba", "numpy"}


def test_backend_resolution_jax_optional():
    """Test that unavailable JAX falls back to NumPy."""
    resolved = bam.utils.auxiliary.resolve_compute_backend("jax")
    assert resolved in {"jax", "numpy"}


# ==================================================================================
# Compose Rotation Tests
# ==================================================================================


def test_rotate_coords_numpy_backend():
    """Test rotate_coords with explicit NumPy backend."""
    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 2  # 90 degrees

    rotated = bam.structural.rotate_coords(coords, angle, axis, backend="numpy")
    assert rotated.shape == coords.shape

    # First coordinate should rotate from [1, 0, 0] toward [0, 1, 0]
    assert np.isclose(rotated[0, 0], 0.0, atol=1e-10)
    assert np.isclose(rotated[0, 1], 1.0, atol=1e-10)


def test_rotate_coords_default_backend():
    """Test rotate_coords without specifying backend (uses default)."""
    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 4

    rotated = bam.structural.rotate_coords(coords, angle, axis)
    assert rotated.shape == coords.shape
    assert not np.allclose(rotated, coords)  # Should be rotated


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_NUMBA, reason="Numba not available")
def test_rotate_coords_numba_backend():
    """Test rotate_coords with Numba backend if available."""
    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 4

    # NumPy reference
    rotated_numpy = bam.structural.rotate_coords(coords, angle, axis, backend="numpy")

    # Numba version
    rotated_numba = bam.structural.rotate_coords(coords, angle, axis, backend="numba")

    # Results should match
    assert np.allclose(rotated_numpy, rotated_numba, atol=1e-10)


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_rotate_coords_jax_backend():
    """Test rotate_coords with JAX backend if available."""
    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )

    axis = np.array([0.0, 0.0, 1.0])
    angle = np.pi / 4

    # NumPy reference
    rotated_numpy = bam.structural.rotate_coords(coords, angle, axis, backend="numpy")

    # JAX version
    rotated_jax = bam.structural.rotate_coords(coords, angle, axis, backend="jax")

    # Results should match
    assert np.allclose(rotated_numpy, rotated_jax, atol=1e-9)


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_rotate_coords_consistency_jax():
    """Test JAX backend consistency against NumPy."""
    np.random.seed(42)
    coords = np.random.randn(10, 3)
    axis = np.array([1.0, 1.0, 1.0])
    axis = axis / np.linalg.norm(axis)
    angle = 0.5

    rotated_numpy = bam.structural.rotate_coords(coords, angle, axis, backend="numpy")
    rotated_jax = bam.structural.rotate_coords(coords, angle, axis, backend="jax")
    assert np.allclose(
        rotated_numpy, rotated_jax, atol=1e-8
    ), "Backend jax produces different results"


def test_rotate_coords_consistency_numba():
    """Test Numba backend consistency against NumPy when available."""
    if not bam.utils.auxiliary.HAS_NUMBA:
        pytest.xfail("Numba not available on this machine")

    np.random.seed(42)
    coords = np.random.randn(10, 3)
    axis = np.array([1.0, 1.0, 1.0])
    axis = axis / np.linalg.norm(axis)
    angle = 0.5

    rotated_numpy = bam.structural.rotate_coords(coords, angle, axis, backend="numpy")
    rotated_numba = bam.structural.rotate_coords(coords, angle, axis, backend="numba")
    assert np.allclose(
        rotated_numpy, rotated_numba, atol=1e-8
    ), "Backend numba produces different results"


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_direct_jax_backend_returns_jax_array():
    """Direct JAX backend outputs should be JAX arrays and match NumPy results."""
    import jax
    from buildamol.backends import jax as jax_backend
    import buildamol.structural.base as structural_base

    jnp = bam.utils.auxiliary.get_jax_numpy()
    coords = jnp.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    axis = jnp.array([0.0, 0.0, 1.0])
    angle = 0.25

    rotated = jax_backend.rotate_coords(coords, angle, axis)
    dists = jax_backend.euclidean_distances(coords, coords)
    rot = jax_backend.rotation_matrix(axis, angle)

    assert isinstance(rotated, jax.Array)
    assert isinstance(dists, jax.Array)
    assert isinstance(rot, jax.Array)

    coords_np = np.asarray(coords)
    axis_np = np.asarray(axis)

    rotated_np = bam.structural.rotate_coords(
        coords_np, angle, axis_np, backend="numpy"
    )
    dists_np = structural_base._euclidean_distances(coords_np, coords_np)
    rot_np = structural_base._rotation_matrix(axis_np, angle)

    assert isinstance(rotated_np, np.ndarray)
    assert isinstance(dists_np, np.ndarray)
    assert isinstance(rot_np, np.ndarray)

    assert np.allclose(rotated_np, np.asarray(rotated), atol=1e-9)
    assert np.allclose(dists_np, np.asarray(dists), atol=1e-9)
    assert np.allclose(rot_np, np.asarray(rot), atol=1e-9)


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_dispatched_jax_backend_returns_jax_array():
    """Public API dispatch to JAX backend should preserve JAX array outputs."""
    import jax

    jnp = bam.utils.auxiliary.get_jax_numpy()
    coords = jnp.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    axis = jnp.array([0.0, 0.0, 1.0])
    angle = 0.5

    rotated = bam.structural.rotate_coords(coords, angle, axis, backend="jax")
    assert isinstance(rotated, jax.Array)


def test_direct_numba_backend_returns_numpy_array():
    """Direct import from Numba backend should return NumPy arrays."""
    if not bam.utils.auxiliary.HAS_NUMBA:
        pytest.xfail("Numba not available on this machine")

    from buildamol.backends import numba as numba_backend

    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    axis = np.array([0.0, 0.0, 1.0])
    angle = 0.25

    rotated = numba_backend.rotate_coords(coords, angle, axis)
    dists = numba_backend.euclidean_distances(coords, coords)
    rot = numba_backend.rotation_matrix(axis, angle)

    assert isinstance(rotated, np.ndarray)
    assert isinstance(dists, np.ndarray)
    assert isinstance(rot, np.ndarray)


# ==================================================================================
# Global Flag Compatibility Tests
# ==================================================================================


def test_use_all_numba_sets_backend():
    """Test that use_all_numba() affects backend resolution."""
    original_backend = bam.utils.auxiliary.resolve_compute_backend(None)

    if bam.utils.auxiliary.HAS_NUMBA:
        bam.use_all_numba()
        resolved = bam.utils.auxiliary.resolve_compute_backend(None)
        assert resolved == "numba" or resolved == "numpy"

        # Reset
        bam.dont_use_numba()


def test_use_jax_sets_backend():
    """Test that use_jax() affects backend resolution."""
    if bam.utils.auxiliary.HAS_JAX:
        bam.use_jax()
        resolved = bam.utils.auxiliary.resolve_compute_backend(None)
        assert resolved == "jax"

        # Reset
        bam.dont_use_jax()


def test_dont_use_numba_resets_backend():
    """Test that dont_use_numba() resets backend."""
    if bam.utils.auxiliary.HAS_NUMBA:
        bam.use_all_numba()
        bam.dont_use_numba()
        resolved = bam.utils.auxiliary.resolve_compute_backend(None)
        assert resolved == "numpy"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
