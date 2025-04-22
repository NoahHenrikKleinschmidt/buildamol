"""
An extension to facilitate solvating molecules in a box of water.
"""

AVAILABLE_SOLVATION_BACKENDS = ["pdbfixer"]

SOLVATION_BACKEND = AVAILABLE_SOLVATION_BACKENDS[0]


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

    import importlib

    # do from backend import solvate as solvate
    globals()["solvate"] = importlib.import_module(
        f"buildamol.extensions.molecular_dynamics.solvate.backend_{backend}"
    ).solvate
    # make sure that the docstring is updated
    globals()["solvate"].__doc__ = importlib.import_module(
        f"buildamol.extensions.molecular_dynamics.solvate.backend_{backend}"
    ).solvate.__doc__


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
