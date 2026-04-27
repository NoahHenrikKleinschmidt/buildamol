"""
An extension to facilitate solvating molecules in a box of water.
"""

import importlib
import warnings

AVAILABLE_SOLVATION_BACKENDS = ["pdbfixer", "biobb"]

SOLVATION_BACKEND = AVAILABLE_SOLVATION_BACKENDS[0]
SOLVATION_IMPORT_ERROR = None


def _missing_backend_solvate(*args, **kwargs):
    """
    Placeholder used when the currently selected backend is not installed.
    """
    raise ImportError(
        "The selected solvation backend is not available. Install optional "
        "dependencies for pdbfixer or biobb_amber, then call "
        "set_solvation_backend(...)."
    )


def available_solvation_backends():
    """
    Get the available solvation backends

    Returns
    -------
    list of str
        The available solvation backends
    """
    return AVAILABLE_SOLVATION_BACKENDS


def set_solvation_backend(backend: str):
    """
    Set the default solvation backend to use

    Parameters
    ----------
    backend : str
        The backend to use, either 'easydock'
    """
    global SOLVATION_BACKEND
    if backend not in AVAILABLE_SOLVATION_BACKENDS:
        raise ValueError(f"Invalid solvation backend '{backend}'")
    SOLVATION_BACKEND = backend

    module_name = f"buildamol.extensions.molecular_dynamics.solvate.backend_{backend}"
    try:
        module = importlib.import_module(module_name)
        globals()["solvate"] = module.solvate
        globals()["solvate"].__doc__ = module.solvate.__doc__
        globals()["SOLVATION_IMPORT_ERROR"] = None
    except Exception as exc:  # pragma: no cover - defensive fallback for optional deps
        globals()["solvate"] = _missing_backend_solvate
        globals()["SOLVATION_IMPORT_ERROR"] = exc
        warnings.warn(
            f"Solvation backend '{backend}' is unavailable: {exc}",
            RuntimeWarning,
            stacklevel=2,
        )


def get_solvation_backend():
    """
    Get the default solvation backend

    Returns
    -------
    str
        The default solvation backend
    """
    return SOLVATION_BACKEND


set_solvation_backend(SOLVATION_BACKEND)


__all__ = [
    "AVAILABLE_SOLVATION_BACKENDS",
    "SOLVATION_BACKEND",
    "SOLVATION_IMPORT_ERROR",
    "available_solvation_backends",
    "set_solvation_backend",
    "get_solvation_backend",
    "solvate",
]
