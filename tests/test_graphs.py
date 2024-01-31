"""
Tests to check the behaviour of the bam.AtomGraph and bam.ResidueGraph object
"""

import pytest

from copy import deepcopy
import random
import numpy as np
import buildamol as bam
import tests.base as base
import timeit

bam.load_sugars()

# =================================================================
# AtomGraph Tests
# =================================================================


def test_atom_graph_from_molecule():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol.lock_all()
    graph = bam.graphs.AtomGraph.from_molecule(mol)

    assert graph is not None, "No molecule is made"

    _received = len(list(graph.bonds))
    _expected = len(list(mol.bonds))
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    assert len(graph._locked_edges) == 0, "Molecule is not locked"

    graph = bam.graphs.AtomGraph.from_molecule(mol, locked=True)

    assert len(graph._locked_edges) == len(mol.locked_bonds), "Molecule is not locked"


def test_atom_graph_pdb_one_residue_is_non_empty():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = bam.graphs.AtomGraph.from_molecule(mol)
    _received = len(list(mol.bonds))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _received = len(list(mol.nodes))
    _expected = 24
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_pdb_multi_residue_is_non_empty():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds()
    mol = bam.graphs.AtomGraph.from_molecule(mol)
    _received = len(list(mol.bonds))
    assert _received > 0, f"Expected to find bonds, got {_received}"

    _received = len(list(mol.atoms))
    _expected = 246
    assert _received == _expected, f"Expected {_expected} atoms, got {_received}"


def test_atom_graph_one_residue_get_neighbors():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    graph = bam.graphs.AtomGraph.from_molecule(mol)

    try:
        neigs = graph.get_neighbors("C1")
    except Exception as e:
        pass
    else:
        raise AssertionError("Expected an exception, got none")

    neigs = graph.get_neighbors(mol.get_atom("C1"))
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_multi_residue_get_neighbors():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds()
    graph = bam.graphs.AtomGraph.from_molecule(mol)

    try:
        neigs = graph.get_neighbors("C1")
    except Exception as e:
        pass
    else:
        raise AssertionError("Expected an exception, got none")

    neigs = graph.get_neighbors(mol.get_atoms("C1", by="id"))
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 11
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_atom_graph_get_descendants():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = bam.graphs.AtomGraph.from_molecule(mol)

    o3 = next(i for i in mol.structure.get_atoms() if i.id == "O3")
    ho3 = next(i for i in mol.structure.get_atoms() if i.id == "HO3")

    _received = mol.get_descendants(ho3, o3)
    _received = len(_received)
    _expected = 22
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = mol.get_descendants(o3, ho3)
    _received = len(_received)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    # if there is only one atom, then no directionality is available and no descendants
    # can be found!
    pytest.raises(KeyError, mol.get_descendants, o3, o3)

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    _received = mol.get_descendants(c5, c6)
    _received = len(_received)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"


def test_atom_graph_rotate_descendants_only():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = bam.graphs.AtomGraph.from_molecule(mol)

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    descendants = mol.get_descendants(c5, c6)
    others = set(i for i in mol.atoms if i not in descendants)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])

    mol.rotate_around_edge(c5, c6, np.radians(35), descendants_only=True)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])

    assert np.all(current_others == new_others), "Other atoms have also moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"


def test_atom_graph_rotate_all():
    mol = bam.Molecule.from_pdb(base.MANNOSE)
    mol.infer_bonds()
    mol = bam.graphs.AtomGraph.from_molecule(mol)

    c6 = next(i for i in mol.structure.get_atoms() if i.id == "C6")
    c5 = next(i for i in mol.structure.get_atoms() if i.id == "C5")

    descendants = mol.get_descendants(c5, c6)
    others = set(i for i in mol.atoms if i not in descendants)
    others.remove(c6)
    others.remove(c5)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])
    current_ref = np.array((c5.coord, c6.coord))

    mol.rotate_around_edge(c5, c6, np.radians(35), descendants_only=False)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])
    new_ref = np.array((c5.coord, c6.coord))

    assert not np.all(current_others == new_others), "Other atoms have not moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"
    assert np.allclose(current_ref, new_ref), "Reference atoms have moved"


# =================================================================
# ResidueGraph tests
# =================================================================


def test_residue_graph_from_molecule():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.reindex()
    mol.infer_bonds(restrict_residues=False)
    mol.get_residue_connections()
    mol.lock_all()
    graph_simple = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    v = mol.draw()
    for edge in mol.get_residue_connections():
        v.draw_vector(
            f"""{edge[0].full_id[3:]} ---> {edge[1].full_id[3:]}""",
            edge[0].coord,
            edge[1].coord,
            color="magenta",
            elongate=1.2,
        )
    v.show()

    assert graph_simple is not None, "No molecule is made"

    _received = len(list(graph_simple.bonds))
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    assert len(graph_simple._locked_edges) == 0, "Molecule is not locked"

    v = graph_simple.draw()
    v.show()

    graph_detailed = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)
    graph_detailed.make_detailed(False)
    assert graph_detailed is not None, "No molecule is made"

    graph_detailed.lock_centers()
    v = graph_detailed.draw()
    v.draw_edges(*graph_detailed.get_locked_edges(), color="magenta")
    v.draw_edges(*graph_detailed.get_unlocked_edges(), color="limegreen")
    v.show()

    _received = len(list(graph_detailed.bonds))
    _expected = 40
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"


def test_residue_graph_multi_residue_get_neighbors():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)

    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)
    graph.make_detailed(False)
    graph.show()

    pytest.raises(KeyError, graph.get_neighbors, "C1")

    # get the C1 of the BMA residue (seqid=4)
    neigs = graph.get_neighbors(mol.get_atom("C1", residue=4))
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 2
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    pytest.raises(KeyError, graph.get_neighbors, "MAN")

    neigs = graph.get_neighbors(mol.get_residues("MAN", by="name"))
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    bma = mol.get_residue("BMA", by="name")
    neigs = graph.get_neighbors(bma)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = graph.get_neighbors(bma, 2)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 6
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = graph.get_neighbors(bma, 2, "at")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"


def test_residue_graph_residue_order():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol = mol.make_residue_graph(detailed=False)

    mol2 = bam.Molecule.from_pdb(base.MANNOSE9)
    mol2.infer_bonds(restrict_residues=False)
    mol2 = mol2.make_residue_graph(detailed=False)

    mol3 = bam.Molecule.from_pdb(base.MANNOSE9)
    mol3.infer_bonds(restrict_residues=False)
    mol3 = mol3.make_residue_graph(detailed=False)

    mol4 = bam.Molecule.from_pdb(base.MANNOSE9)
    mol4.infer_bonds(restrict_residues=False)
    mol4 = mol4.make_residue_graph(detailed=False)

    for i in range(10):
        rdx = random.randint(0, len(mol.residues) - 1)
        assert (
            mol.residues[rdx].serial_number == mol2.residues[rdx].serial_number
        ), "Residue order is not preserved"
        assert (
            mol.residues[rdx].serial_number == mol3.residues[rdx].serial_number
        ), "Residue order is not preserved"
        assert (
            mol.residues[rdx].serial_number == mol4.residues[rdx].serial_number
        ), "Residue order is not preserved"


def test_residue_graph_get_descendants():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    v = graph.draw()

    nag2 = mol.residues[0]
    nag3 = mol.residues[1]

    v.draw_point("nag2", nag2.coord, color="red")
    v.draw_point("nag3", nag3.coord, color="blue")
    v.show()

    _received = graph.get_descendants(nag2, nag3)
    _received = len(_received)
    _expected = 9
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = graph.get_descendants(nag3, nag2)
    _received = len(_received)
    _expected = 0
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    bma = mol.residues[2]

    _received = graph.get_descendants(nag3, bma)
    _received = len(_received)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"

    _received = graph.get_descendants(bma, nag3)
    _received = len(_received)
    _expected = 1
    assert _received == _expected, f"Expected {_expected} descendants, got {_received}"


def test_residue_graph_rotate_descendants_only():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    v = graph.draw()

    nag3 = mol.residues[1]
    bma = mol.residues[2]

    v.draw_point("nag3", nag3.coord, color="red")
    v.draw_point("bma", bma.coord, color="blue")

    descendants = graph.get_descendants(nag3, bma)
    others = set(i for i in graph.residues if i not in descendants)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])

    for i in range(5):
        graph.rotate_around_edge(nag3, bma, np.radians(10), descendants_only=True)
        v.draw_edges(*graph.edges, color="magenta")

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])

    assert np.all(current_others == new_others), "Other residues have also moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"

    v.show()


def test_residue_graph_rotate_descendants_only_detailed():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=True)

    v = graph.draw()

    t1 = timeit.timeit()
    cons = mol.get_residue_connections()
    for i in range(20):
        atom1, atom2 = cons.pop()

        v.draw_point(str(i) + " atom1", atom1.coord, color="red")
        v.draw_point(str(i) + " atom2", atom2.coord, color="blue")

        descendants = graph.get_descendants(atom1, atom2)
        others = set(i for i in graph.residues if i not in descendants)

        current_descendants = np.array([i.coord for i in descendants])
        current_others = np.array([i.coord for i in others])

        for j in range(5):
            graph.rotate_around_edge(
                atom1, atom2, np.radians(10), descendants_only=True
            )
            v.draw_edges(*graph.edges, color="magenta")

        new_descendants = np.array([i.coord for i in descendants])
        new_others = np.array([i.coord for i in others])

        assert np.all(current_others == new_others), "Other residues have also moved!"
        assert not np.allclose(
            current_descendants, new_descendants
        ), "Descendants have not moved"

    t2 = timeit.timeit()
    v.show()
    print("Total: ", t2 - t1)


def test_residue_graph_rotate_all():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    nag3 = mol.residues[1]
    bma = mol.residues[2]

    descendants = graph.get_descendants(nag3, bma)
    others = set(i for i in graph.residues if i not in descendants)
    others.remove(nag3)
    others.remove(bma)

    current_descendants = np.array([i.coord for i in descendants])
    current_others = np.array([i.coord for i in others])
    current_ref = np.array((nag3.coord, bma.coord))

    graph.rotate_around_edge(nag3, bma, np.radians(35), descendants_only=False)

    new_descendants = np.array([i.coord for i in descendants])
    new_others = np.array([i.coord for i in others])
    new_ref = np.array((nag3.coord, bma.coord))

    assert not np.all(current_others == new_others), "Other residues have not moved!"
    assert not np.allclose(
        current_descendants, new_descendants
    ), "Descendants have not moved"
    assert np.allclose(current_ref, new_ref), "Reference residues have moved"


def test_residue_graph_detailed():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    mol = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)

    assert len(mol.nodes) == 11, "Wrong number of nodes"
    assert len(mol.bonds) == 10, "Wrong number of bonds"

    mol.make_detailed(False)

    assert len(mol.edges) == 40, "Wrong number of edges"
    assert len(mol.nodes) == 41, "Wrong number of nodes"


def test_residue_graph_detailed_get_neighbors():
    mol = bam.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    graph = bam.graphs.ResidueGraph.from_molecule(mol, detailed=False)
    graph.make_detailed(False)

    bma = mol.get_residue(3)
    neigs = graph.get_neighbors(bma)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"
    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    ids = [i.id for i in neigs]
    assert "C6" in ids, "Wrong neighbors: " + str(ids)
    assert "C1" in ids, "Wrong neighbors: " + str(ids)

    has_residues = False
    has_atoms = False
    for node in graph.nodes:
        if node.__class__.__name__ == "Residue":
            has_residues = True
        if node.__class__.__name__ == "Atom":
            has_atoms = True

    assert has_residues, "Residues not found"
    assert has_atoms, "Atoms not found"


def test_atom_graph_lock():
    glc = bam.Molecule.from_compound("GLC")
    glc.repeat(5, "14bb")
    glc.lock_all()

    g = glc._AtomGraph

    assert len(g.nodes) == len(glc.atoms)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) > 0

    glc.unlock_all()
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 0

    glc.lock_bond(4, 5)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 1

    old = len(glc.bonds)
    old_edges = len(g.edges)

    glc.add_bond(5, 12)
    assert len(glc.bonds) == old + 1
    assert len(g.edges) == old_edges + 1
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 1

    glc.lock_bond(5, 12)
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 2

    glc.remove_bond(5, 12)
    assert len(glc.bonds) == old
    assert len(g.edges) == old_edges
    assert len(set(g._locked_edges)) == len(set(glc.locked_bonds)) == 2
