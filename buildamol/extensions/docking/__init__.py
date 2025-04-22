"""
An extension for facilitating docking of ligand molecules to proteins.
"""

AVAILABLE_DOCKING_BACKENDS = [
    "easydock", "dockstring"
]

DOCKING_BACKEND = AVAILABLE_DOCKING_BACKENDS[0]

def available_docking_backends():
    """
    Get the available docking backends

    Returns
    -------
    list of str
        The available docking backends
    """
    return AVAILABLE_DOCKING_BACKENDS

def set_docking_backend(backend: str):
    """
    Set the default docking backend to use

    Parameters
    ----------
    backend : str
        The backend to use, either 'easydock'
    """
    global DOCKING_BACKEND
    if backend not in AVAILABLE_DOCKING_BACKENDS:
        raise ValueError(f"Invalid docking backend '{backend}'")
    DOCKING_BACKEND = backend

    import importlib
    # do from backend import __dock__ as dock
    globals()["dock"] = importlib.import_module(f"buildamol.extensions.docking.{backend}").__dock__
    # make sure that the docstring is updated
    globals()["dock"].__doc__ = importlib.import_module(f"buildamol.extensions.docking.{backend}").dock.__doc__

def get_docking_backend():
    """
    Get the default docking backend

    Returns
    -------
    str
        The default docking backend
    """
    return DOCKING_BACKEND

set_docking_backend(DOCKING_BACKEND)