"""Unified hub for CHARMM-related functionality.

This module provides a single import surface that re-exports CHARMM-specific
naming, atom typing, topology, and helper labeling utilities without moving
implementation code.

Examples
--------
>>> from buildamol.extensions import charmm
>>> typer = charmm.CHARMMTyper.from_file("top_all36_carb.rtf")
>>> top = charmm.read_topology("topology.pkl")
"""

from buildamol.extensions.naming.charmm import (
    CHARMMAtomNameGraphEngine,
    CHARMMResidueNameGraphEngine,
    CHARMMResidueNameLookupEngine,
    charmm_to_pdb,
    # get_note,
    # is_ambiguous,
    pdb_to_charmm,
    # pdb_to_charmm_mapping,
    # print_summary,
    translate_resname,
)
from buildamol.extensions.molecular_dynamics.atom_typing.charmm_typer import (
    CHARMMTyper,
    rename_for_charmm_typing,
    type_with_charmm,
)
from buildamol.extensions.charmm.psf import PSFMaker, write_psf

from buildamol.resources.charmm import (
    CHARMMTopology,
    add_linkage,
    add_patch,
    available_linkages,
    available_patches,
    export_topology,
    get_default_topology,
    get_linkage,
    get_patch,
    has_linkage,
    has_patch,
    read_topology,
    restore_default_topology,
    save_topology,
    set_default_topology,
)

__all__ = [
    "CHARMMTopology",
    "CHARMMTyper",
    "CHARMMAtomNameGraphEngine",
    "CHARMMResidueNameGraphEngine",
    "CHARMMResidueNameLookupEngine",
    "type_with_charmm",
    "rename_for_charmm_typing",
    "PSFMaker",
    "write_psf",
    "translate_resname",
    "pdb_to_charmm",
    # "pdb_to_charmm_mapping",
    "charmm_to_pdb",
    # "get_note",
    # "is_ambiguous",
    # "print_summary",
    "read_topology",
    "save_topology",
    "export_topology",
    "get_default_topology",
    "set_default_topology",
    "restore_default_topology",
    "has_patch",
    "available_patches",
    "available_linkages",
    "add_patch",
    "add_linkage",
    "has_linkage",
    "get_patch",
    "get_linkage",
]
