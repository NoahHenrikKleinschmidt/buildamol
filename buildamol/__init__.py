import buildamol.utils as utils

_utils = utils

import buildamol.utils.visual as visual
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
    use_ic,
    dont_use_ic,
)
from buildamol.utils.visual import MoleculeViewer3D, Py3DmolViewer
from buildamol.optimizers import *

from buildamol.utils.info import __version__, __author__

# a little hack to make sure the utils module is not the optimizers.utils...
utils = _utils
del _utils
