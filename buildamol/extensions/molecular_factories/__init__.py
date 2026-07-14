"""
This package contains classes that can be used to generate molecular structures automatically.
"""

from .derivator import Derivator
from .assembler import Assembler
from .generator import Generator
from .fragment_library import (
    load_fragment_library,
    make_brics_fragments,
    get_chembl_fragments,
)

# Pipeline building blocks
from .base import ChainableBlock, ChainedBlock, Context, MultiContext
from .sources import Choice, Compound
from .modifiers import (
    Modify,
    Random,
    SetLinkerAtoms,
    SetAttachResidue,
    FindLinkerAtoms,
    Connect,
    Forge,
)
from .backend import set_backend, get_backend, backend
from .optimizer import Optimizer
