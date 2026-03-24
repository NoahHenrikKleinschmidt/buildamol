"""
Backend dispatch API and utilities.

Provides the decorator and resolver logic for routing function calls
to backend-specific implementations.
"""

import functools
import importlib
import warnings
from typing import Optional, Callable, Any


def backend_dispatched(module_name: str, func_name: Optional[str] = None):
    """
    Decorator for dispatching function calls to backend-specific implementations.

    Automatically adds a `backend=None` keyword argument and routes to the
    appropriate backend function if available, falling back to the decorated
    (NumPy) implementation.

    Parameters
    ----------
    module_name : str
        Backend module name (e.g., "structural", "optimizers")
    func_name : str, optional
        Backend function name. If None, uses the decorated function's name.

    Returns
    -------
    Callable
        Decorated function with automatic backend dispatch.

    Example
    -------
    >>> from buildamol.backends.api import backend_dispatched

    >>> @backend_dispatched("structural", "rotate_coords")
    ... def rotate_coords(coords, angle, axis):
    ...     # NumPy implementation
    ...     return rotated_coords

    >>> # Users can now call:
    >>> rotate_coords(coords, angle, axis)  # Uses NumPy
    >>> rotate_coords(coords, angle, axis, backend="jax")  # Uses JAX if available
    """

    def decorator(func: Callable) -> Callable:
        actual_func_name = func_name or func.__name__

        @functools.wraps(func)
        def wrapper(*args, backend=None, **kwargs):
            import buildamol.utils.auxiliary as aux

            backend = aux.resolve_compute_backend(backend)

            if backend != "numpy":
                try:
                    backend_fn = get_backend_function(
                        module_name, actual_func_name, backend
                    )
                    if backend_fn is not None:
                        return backend_fn(*args, **kwargs)
                except (ImportError, AttributeError, ValueError):
                    pass  # Fall back to default NumPy implementation

            return func(*args, **kwargs)

        return wrapper

    return decorator


def get_backend_function(
    module_name: str, func_name: str, backend: str
) -> Optional[Callable]:
    """
    Fetch a backend-specific function implementation.

    Parameters
    ----------
    module_name : str
        Backend module name (e.g., "structural", "optimizers")
    func_name : str
        Function name to retrieve
    backend : str
        Backend name (numpy, numba, jax)

    Returns
    -------
    Callable or None
        The backend-specific function, or None if using NumPy (default)

    Raises
    ------
    ImportError
        If the backend module cannot be imported
    AttributeError
        If the function is not found in the backend module
    ValueError
        If the backend name is not recognized
    """
    if backend == "numpy":
        return None

    backend_module = importlib.import_module(
        f"buildamol.backends.{backend}.{module_name}"
    )
    return getattr(backend_module, func_name)
