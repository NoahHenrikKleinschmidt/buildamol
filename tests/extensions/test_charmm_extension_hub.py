from buildamol.extensions import charmm as hub
from buildamol.extensions.charmm.psf import PSFMaker, write_psf

from buildamol.extensions.molecular_dynamics.atom_typing.charmm_typer import (
    CHARMMTyper,
    rename_for_charmm_typing,
    type_with_charmm,
)
from buildamol.extensions.molecular_dynamics.psf import (
    PSFMaker as LegacyPSFMaker,
    write_psf as legacy_write_psf,
)
from buildamol.extensions.naming.charmm import (
    CHARMMAtomNameGraphEngine,
    CHARMMResidueNameGraphEngine,
    CHARMMResidueNameLookupEngine,
    charmm_to_pdb,
    pdb_to_charmm,
    translate_resname,
)
from buildamol.resources.charmm import CHARMMTopology, read_topology


def test_charmm_extension_hub_reexports_primary_symbols():
    assert hub.CHARMMTyper is CHARMMTyper
    assert hub.rename_for_charmm_typing is rename_for_charmm_typing
    assert hub.type_with_charmm is type_with_charmm

    assert hub.CHARMMAtomNameGraphEngine is CHARMMAtomNameGraphEngine
    assert hub.CHARMMResidueNameGraphEngine is CHARMMResidueNameGraphEngine
    assert hub.CHARMMResidueNameLookupEngine is CHARMMResidueNameLookupEngine

    assert hub.translate_resname is translate_resname
    assert hub.pdb_to_charmm is pdb_to_charmm
    assert hub.charmm_to_pdb is charmm_to_pdb

    assert hub.CHARMMTopology is CHARMMTopology
    assert hub.read_topology is read_topology
    assert hub.PSFMaker is PSFMaker
    assert hub.write_psf is write_psf
    assert LegacyPSFMaker is PSFMaker
    assert legacy_write_psf is write_psf


def test_charmm_extension_hub_exports_documented_symbols():
    expected = {
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
        "charmm_to_pdb",
        "read_topology",
        "save_topology",
        "export_topology",
        "get_default_topology",
        "set_default_topology",
        "restore_default_topology",
    }
    assert expected.issubset(set(hub.__all__))
