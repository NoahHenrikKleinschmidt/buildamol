"""Compatibility shim for PSF functionality.

The implementation lives in ``buildamol.extensions.charmm.psf``.
This module re-exports the same public API to avoid breaking existing imports.
"""

from buildamol.extensions.charmm.psf import PSFMaker, write_psf

__all__ = ["PSFMaker", "write_psf"]
