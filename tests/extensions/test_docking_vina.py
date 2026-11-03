"""
Tests for the in-process AutoDock Vina docking backend.
"""

import warnings

import numpy as np
import pytest

import buildamol as bam
import buildamol.utils.auxiliary as aux
import tests.base as base

pytestmark = pytest.mark.skipif(
    not (aux.has_package("vina") and aux.has_package("meeko")),
    reason="requires vina and meeko",
)

POCKET = base.FILES / "vina_pocket.pdbqt"
CENTER = (17.61, -18.43, 32.75)
BOX = (15, 15, 15)


@pytest.fixture(scope="module")
def docker():
    from buildamol.extensions.docking.backend_vina import Docker

    with Docker(seed=1) as d:
        d.setup(str(POCKET), CENTER, BOX)
        yield d


def _ligand(smiles="CC(=O)Nc1ccc(O)cc1", id="LIG"):
    lig = bam.Molecule.from_smiles(smiles, id=id)
    lig.move(np.array(CENTER) - lig.center_of_geometry)
    return lig


def _bond_lengths(mol):
    return np.array(
        [np.linalg.norm(b.atom1.coord - b.atom2.coord) for b in mol.get_bonds()]
    )


def test_vina_is_default_backend():
    import buildamol.extensions.docking as docking

    assert docking.get_docking_backend() == "vina"


def test_score_does_not_move_ligand(docker):
    lig = _ligand()
    scored = docker.score(lig)
    assert scored is not lig
    assert len(scored.models) == 1
    assert isinstance(scored.docking_score, float)
    assert np.allclose(scored.get_coords(), lig.get_coords(), atol=1e-3)
    assert [a.id for a in scored.get_atoms()] == [a.id for a in lig.get_atoms()]


def test_optimize_improves_score(docker):
    lig = _ligand()
    before = docker.score(lig).docking_score
    optimized = docker.optimize(lig)
    assert optimized.docking_score <= before
    assert len(optimized.models) == 1
    assert np.allclose(_bond_lengths(optimized), _bond_lengths(lig), atol=0.1)


def test_dock_returns_sorted_poses(docker):
    lig = _ligand()
    docked = docker.dock(lig, exhaustiveness=2, n_poses=3)
    assert 1 <= len(docked.models) <= 3
    scores = [m.docking_score for m in docked.models]
    assert scores == sorted(scores)
    assert docked.docking_score == scores[0]
    assert [a.id for a in docked.get_atoms()] == [a.id for a in lig.get_atoms()]
    # the poses are mapped back onto the right atoms
    assert np.allclose(_bond_lengths(docked), _bond_lengths(lig), atol=0.1)
    # the input ligand is not modified
    assert len(lig.models) == 1


def test_dock_smiles(docker):
    docked = docker.dock("CCOc1ccccc1C(=O)N", exhaustiveness=2, n_poses=1)
    assert isinstance(docked, bam.Molecule)
    assert docked.docking_score < 0


def test_score_outside_box_raises(docker):
    lig = _ligand()
    lig.move(np.array([50.0, 0, 0]))
    with pytest.raises(ValueError, match="outside the docking box"):
        docker.score(lig)


def test_unsupported_atom_types_raise(docker):
    with pytest.raises(ValueError, match="not supported by Vina"):
        docker.score(_ligand("C[Si](C)(C)c1ccccc1"))


def test_xs_types_depend_on_bonding():
    from buildamol.extensions.docking import backend_vina as bv

    # both oxygens are AutoDock 'OA' but Vina needs separate maps for them
    assert bv._pdbqt_xs_types(bv._smiles_to_pdbqt("CCO", CENTER)) == {"C_H", "C_P", "O_DA"}
    assert bv._pdbqt_xs_types(bv._smiles_to_pdbqt("CC(=O)O", CENTER)) == {
        "C_H", "C_P", "O_A", "O_DA"
    }  # fmt: skip


def test_maps_extended_for_new_xs_types(monkeypatch):
    from buildamol.extensions.docking import backend_vina as bv

    monkeypatch.setattr(bv, "_REFERENCE_SMILES", ("CC",))
    session = bv._VinaSession(str(POCKET), CENTER, BOX, "vina", 0, 1)
    try:
        assert session._map_types == {"C_H"}
        # without extending the maps first, Vina would terminate the process here
        for smiles, recompute in (("CCO", True), ("CC(=O)O", True), ("CO", False)):
            n = len(session._map_ligands)
            energies, _ = session.run(bv._smiles_to_pdbqt(smiles, CENTER), "score")
            assert np.isfinite(energies[0])
            assert (len(session._map_ligands) > n) == recompute
        assert session._map_types == {"C_H", "C_P", "O_A", "O_DA"}
    finally:
        session.close()


def test_dock_batch_in_parallel(docker):
    ligands = [
        _ligand("CC(=O)Nc1ccc(O)cc1", id="A"),
        _ligand("c1ccc2[nH]ccc2c1", id="B"),
        _ligand("C[Si](C)(C)c1ccccc1", id="C"),
    ]
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        docked = docker.dock(ligands, exhaustiveness=2, n_poses=1, ncpu=2)
    assert len(docked) == 3
    assert docked[0].docking_score < 0 and docked[1].docking_score < 0
    assert docked[2] is None
    assert any("not supported by Vina" in str(i.message) for i in w)


def test_dock_function():
    from buildamol.extensions.docking.backend_vina import dock

    docked = dock(str(POCKET), _ligand(), CENTER, BOX, mode="optimize")
    assert isinstance(docked.docking_score, float)


def test_pose_mapping_keeps_atom_order():
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from buildamol.extensions.docking import backend_vina as bv

    # this atom order differs from the canonical SMILES order meeko writes into the PDBQT
    rdmol = Chem.AddHs(Chem.MolFromSmiles("OCC(N)C(=O)OCC"))
    AllChem.EmbedMolecule(rdmol, randomSeed=1)
    ref = rdmol.GetConformer().GetPositions()
    pdbqt = bv._rdkit_to_pdbqt(rdmol)
    # the input PDBQT is a valid "pose": mapping it back must reproduce the input (incl. merged hydrogens)
    pose = bv._pdbqt_pose_coords(pdbqt)[0]
    mapped = bv._map_pose(rdmol, ref, pose, bv._pdbqt_index_map(pdbqt))
    assert np.allclose(mapped, ref, atol=5e-3)


def test_optimize_keeps_structure_of_noncanonical_ligand(docker):
    # a ligand whose atom order does not follow the canonical SMILES order
    lig = _ligand("OCC(N)C(=O)NCc1ccccc1O")
    optimized = docker.optimize(lig)
    assert np.allclose(_bond_lengths(optimized), _bond_lengths(lig), atol=0.05)
    docked = docker.dock(lig, exhaustiveness=2, n_poses=2)
    for i in range(len(docked.models)):
        docked.set_model(i)
        assert np.allclose(_bond_lengths(docked), _bond_lengths(lig), atol=0.05)


def test_box_errors_are_informative(docker):
    lig = _ligand()
    lig.move(np.array([9.0, 0, 0]))
    with pytest.raises(ValueError, match="ligand atoms are outside the docking box.*box_around"):
        docker.optimize(lig)


def test_box_around():
    from buildamol.extensions.docking.backend_vina import box_around

    lig = _ligand()
    center, size = box_around(lig, padding=3)
    xyz = lig.get_coords()
    assert np.allclose(center, (xyz.min(0) + xyz.max(0)) / 2)
    assert np.allclose(size, np.ptp(xyz, axis=0) + 6)
