"""
Backend implementations for optimized compute paths.

This package provides backend-specific implementations of key functions
across BuildAMol. Backends are selectable at runtime via global toggles
or per-call via the `backend=` parameter on dispatched functions.

Supported backends:
  - numpy: Default CPU backend (always available)
  - numba: JIT compilation via Numba (optional, requires numba package)
  - jax: JAX functional arrays (optional, requires jax package)

Example
-------
>>> import buildamol as bam
>>> from buildamol.backends.jax import rotate_coords
>>> # Use JAX-specific implementation directly
>>> coords_rotated = rotate_coords(coords, angle, axis)

Or use the public API with backend selection:
>>> coords_rotated = bam.structural.rotate_coords(coords, angle, axis, backend="jax")
"""

__all__ = ["api"]
