"""
The `structural` module contains classes and functions to work with molecular structures.
It is at the heart of BuildAMol functionality and provides most of the useful features.

Almost all functions and classes in this module are integrated into the main BuildAMol API
through methods of the `Molecule` and `Scaffold` classes - so it is usually not necessary to 
use this module directly. However there are some useful features that are not directly integrated
into the API, and in some cases users may want to use this module directly to access them.
"""

from buildamol.structural.infer import *
from buildamol.structural.smiles import *
from buildamol.structural.base import *
from buildamol.structural.patch import (
    Patcher,
    patch,
    __default_keep_keep_patcher__,
)

from buildamol.structural.stitch import (
    Stitcher,
    stitch,
    __default_keep_keep_stitcher__,
)

from buildamol.structural.neighbors import (
    AtomNeighborhood,
    ResidueNeighborhood,
    compute_quartets,
    compute_triplets,
    generate_triplets,
    generate_quartets,
    constraints,
)

import buildamol.structural.geometry as geometry
import buildamol.structural.groups as groups
