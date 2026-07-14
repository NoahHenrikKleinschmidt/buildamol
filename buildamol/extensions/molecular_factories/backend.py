_BACKEND = "buildamol-native"
_VALID_BACKENDS = {"buildamol-native", "rdkit-accelerated"}


def set_backend(name: str) -> None:
    """Set the active molecular-factories backend globally."""
    global _BACKEND
    if name not in _VALID_BACKENDS:
        raise ValueError(
            f"Unknown backend {name!r}. Choose from {sorted(_VALID_BACKENDS)}."
        )
    _BACKEND = name


def get_backend() -> str:
    """Return the currently active backend name."""
    return _BACKEND


class backend:
    """Context manager for transient backend switches.

    Example
    -------
    with mf.backend("rdkit-accelerated"):
        mol = pipeline().molecule
    """

    def __init__(self, name: str) -> None:
        self._name = name
        self._prev: str = _BACKEND

    def __enter__(self) -> "backend":
        self._prev = _BACKEND
        set_backend(self._name)
        return self

    def __exit__(self, *_) -> None:
        set_backend(self._prev)
