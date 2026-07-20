import buildamol.utils as utils
from importlib import import_module

_utils = utils

import buildamol.structural as structural
import buildamol.resources as resources
import buildamol.graphs as graphs
import buildamol.optimizers as optimizers

from buildamol.core import *
from buildamol.resources import *
from buildamol.utils.auxiliary import (
    use_numba,
    use_all_numba,
    dont_use_numba,
    use_jax,
    dont_use_jax,
    use_ic,
    dont_use_ic,
)
from buildamol.optimizers import *

from buildamol.utils.info import __version__, __author__

from buildamol.structural.stitch_algorithm import (
    stitching_algorithm,
    set_stitching_algorithm,
    get_stitching_algorithm,
    grid_optimize,
)

# Keep CHARMM APIs in `buildamol.resources.charmm`, but avoid exposing them at top-level.
for _name in getattr(resources.charmm, "__all__", []):
    globals().pop(_name, None)

# a little hack to make sure the utils module is not the optimizers.utils...
utils = _utils
del _utils


def __getattr__(name):
    if name == "visual":
        module = import_module("buildamol.utils.visual")
        globals()[name] = module
        return module

    if name in ("MoleculeViewer3D", "Py3DmolViewer"):
        module = import_module("buildamol.utils.visual")
        value = getattr(module, name)
        globals()[name] = value
        return value

    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
