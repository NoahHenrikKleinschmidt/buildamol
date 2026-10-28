"""
Utility and auxiliary functions and constants used by buildamol.
Most of these functions are not relevant to the user, but are used internally.
"""

from importlib import import_module


_MODULE_ALIASES = {
    "constants": "buildamol.utils.constants",
    "defaults": "buildamol.utils.defaults",
    "abstract": "buildamol.utils.abstract",
    "visual": "buildamol.utils.visual",
    "convert": "buildamol.utils.convert",
    "pdb": "buildamol.utils.pdb",
    "cif": "buildamol.utils.cif",
    "xml": "buildamol.utils.xml",
    "json": "buildamol.utils.json",
    "sdmol": "buildamol.utils.sdmol",
    "ic": "buildamol.utils.ic",
    "auxiliary": "buildamol.utils.auxiliary",
    "pdbqt": "buildamol.utils.pdbqt",
    "xyz": "buildamol.utils.xyz",
}


_STAR_EXPORT_MODULES = (
    "buildamol.utils.auxiliary",
    "buildamol.utils.defaults",
)


def __getattr__(name):
    module_path = _MODULE_ALIASES.get(name)
    if module_path is not None:
        module = import_module(module_path)
        globals()[name] = module
        return module

    for module_name in _STAR_EXPORT_MODULES:
        module = import_module(module_name)
        if hasattr(module, name):
            value = getattr(module, name)
            globals()[name] = value
            return value

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


def __dir__():
    return sorted(set(globals()) | set(_MODULE_ALIASES))
