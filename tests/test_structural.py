"""
Tests for the auxiliary structure module
"""

from copy import deepcopy
import numpy as np
import pytest
import buildamol as bam
import tests.base as base
import Bio.PDB as bio
import re

MARGIN = 1.5 * 1e-2
MANNOSE = bio.PDBParser().get_structure("MAN", base.MANNOSE)
bam.load_small_molecules()
bam.load_sugars()


# def test_missing_proper_1():
#     to_deletes = {"C1", "C2", "C3", "C4", "C5", "O5"}

#     for to_delete in to_deletes:
#         _man = MANNOSE.copy()
#         _man = next(_man.get_residues())
#         true_coords = _man.child_dict.get(to_delete)
#         assert true_coords is not None, f"No atom {to_delete} found!"

#         true_coords = true_coords.coord

#         _man.detach_child(to_delete)
#         assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

#         bam.structural.fill_missing_atoms(_man)

#         assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

#         new_coords = _man.child_dict.get(to_delete).coord

#         _diff = np.sum(np.abs(new_coords - true_coords))
#         assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


# def test_missing_proper_4():
#     to_deletes = {"H1", "HO1", "HO2", "HO3", "HO4", "HO6", "O1"}

#     for to_delete in to_deletes:
#         _man = MANNOSE.copy()
#         _man = next(_man.get_residues())

#         true_coords = _man.child_dict.get(to_delete)
#         assert true_coords is not None, f"No atom {to_delete} found!"

#         true_coords = true_coords.coord

#         _man.detach_child(to_delete)
#         assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

#         bam.structural.fill_missing_atoms(_man)

#         assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

#         new_coords = _man.child_dict.get(to_delete).coord

#         _diff = np.sum(np.abs(new_coords - true_coords))
#         assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


# def test_missing_improper_1():
#     to_deletes = {"C6", "O3", "O4", "O6"}

#     for to_delete in to_deletes:
#         _man = MANNOSE.copy()
#         _man = next(_man.get_residues())

#         top = deepcopy(bam.resources.get_default_topology())
#         abstract = top.get_residue(_man.resname)
#         for idx, i in enumerate(abstract.internal_coordinates):
#             if i.is_proper:
#                 del abstract.internal_coordinates[idx]

#         true_coords = _man.child_dict.get(to_delete)
#         assert true_coords is not None, f"No atom {to_delete} found!"

#         true_coords = true_coords.coord

#         _man.detach_child(to_delete)
#         assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

#         bam.structural.fill_missing_atoms(_man, top)

#         assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

#         new_coords = _man.child_dict.get(to_delete).coord

#         _diff = np.sum(np.abs(new_coords - true_coords))
#         assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


# def test_missing_multi():
#     man = bam.Molecule.from_compound("MAN")
#     man.repeat(2, "14bb")
#     to_delete = ("O5", "C3", "O2", "O3")
#     man.remove_atoms(*to_delete)
#     for atom in man.get_atoms():
#         atom.__mol__ = man
#     for res in man.residues:
#         bam.structural.fill_missing_atoms(res)


# def test_missing_improper_4():
#     to_deletes = {"H2", "H3", "H4", "H5", "H61", "H62"}

#     for to_delete in to_deletes:
#         _man = MANNOSE.copy()
#         _man = next(_man.get_residues())
#         top = deepcopy(bam.resources.get_default_topology())
#         abstract = top.get_residue(_man.resname)
#         for idx, i in enumerate(abstract.internal_coordinates):
#             if i.is_proper:
#                 del abstract.internal_coordinates[idx]

#         true_coords = _man.child_dict.get(to_delete)
#         assert true_coords is not None, f"No atom {to_delete} found!"

#         true_coords = true_coords.coord

#         _man.detach_child(to_delete)
#         assert _man.child_dict.get(to_delete) is None, "Atom was not deleted!"

#         bam.structural.fill_missing_atoms(_man, top)

#         assert _man.child_dict.get(to_delete) is not None, "Atom was not added again!"

#         new_coords = _man.child_dict.get(to_delete).coord

#         _diff = np.sum(np.abs(new_coords - true_coords))
#         assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


# def test_missing_one_random_atom():
#     _man = MANNOSE.copy()
#     _man = next(_man.get_residues())

#     to_delete = np.random.choice(list(_man.child_dict.keys()), 1)
#     to_delete = to_delete[0]
#     true_coords = _man.child_dict.get(to_delete)
#     assert true_coords is not None, f"No atom {to_delete} found!"

#     true_coords = true_coords.coord

#     _man.detach_child(to_delete)

#     bam.structural.fill_missing_atoms(_man)

#     assert (
#         _man.child_dict.get(to_delete) is not None
#     ), f"Atom {to_delete} was not added again!"

#     new_coords = _man.child_dict.get(to_delete).coord

#     _diff = np.sum(np.abs(new_coords - true_coords))
#     assert _diff < MARGIN, f"[{to_delete}] Difference in coordinates is {_diff=}"


# def test_missing_multiple_random_atoms():
#     _man = MANNOSE.copy()
#     _man = next(_man.get_residues())

#     to_delete = np.random.choice(list(_man.child_dict.keys()), 5, replace=False)

#     true_coords = {i: _man.child_dict.get(i).coord for i in to_delete}
#     for i in to_delete:
#         _man.detach_child(i)

#     bam.structural.fill_missing_atoms(_man)

#     for i in to_delete:
#         assert _man.child_dict.get(i) is not None, f"Atom {i} was not added again!"

#         new_coords = _man.child_dict.get(i).coord

#         _diff = np.sum(np.abs(new_coords - true_coords[i]))
#         assert _diff < MARGIN, f"[{i}] Difference in coordinates is {_diff=}"


# def test_missing_multiple_random_atoms_galactose():
#     GALACTOSE = bio.PDBParser().get_structure("GAL", base.GALACTOSE)

#     _gal = GALACTOSE.copy()
#     _gal = next(_gal.get_residues())

#     to_delete = np.random.choice(list(_gal.child_dict.keys()), 5, replace=False)

#     true_coords = {i: _gal.child_dict.get(i).coord for i in to_delete}
#     for i in to_delete:
#         _gal.detach_child(i)

#     bam.structural.fill_missing_atoms(_gal)

#     for i in to_delete:
#         assert _gal.child_dict.get(i) is not None, f"Atom {i} was not added again!"

#         new_coords = _gal.child_dict.get(i).coord

#         _diff = np.sum(np.abs(new_coords - true_coords[i]))
#         assert _diff < MARGIN, f"[{i}] Difference in coordinates is {_diff=}"


# def test_missing_multiple_random_atoms_mannose9():
#     _man = bio.PDBParser().get_structure("MAN9", base.MANNOSE9)

#     atoms = list(_man.get_atoms())
#     to_delete = np.random.choice(atoms, 15, replace=False)

#     true_coords = [None] * len(to_delete)
#     parents = [i.get_parent() for i in to_delete]
#     for idx, i in enumerate(to_delete):
#         parent = i.get_parent()
#         true_coords[idx] = i.coord
#         parent.detach_child(i.id)

#     bam.structural.fill_missing_atoms(_man)

#     for i, true_coord, parent in zip(to_delete, true_coords, parents):
#         assert parent.child_dict.get(i.id) is not None, f"Atom {i} was not added again!"

#         new_coords = i.coord

#         _diff = np.sum(np.abs(new_coords - true_coord))
#         assert _diff < MARGIN, f"[{i}] Difference in coordinates is {_diff=}"


def test_apply_standard_bonds():
    bam.load_sugars()

    bonds = bam.structural.apply_reference_bonds(MANNOSE)

    _recieved = len(bonds)
    _expected = 24
    _what = "bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    bonds = [set((i.id, j.id)) for i, j, order in bonds]

    _bond = set(("C5", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C6", "O6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C4", "C3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "O3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"
    bam.unload_sugars()


def test_apply_standard_bonds_one_atom():
    bam.load_sugars()
    atom = {i.id: i for i in MANNOSE.get_atoms()}
    atom = atom.get("C1")

    bonds = bam.structural.apply_reference_bonds(atom)
    bonds = [set((i.id, j.id)) for i, j, order in bonds]

    _recieved = len(bonds)
    _expected = 4
    _what = "bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "O1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"
    bam.unload_sugars()


def test_infer_bonds():
    bonds = bam.structural.infer_bonds(MANNOSE)

    _recieved = len(bonds)
    _expected = 24
    _what = "bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    bonds = [set((i.id, j.id)) for i, j in bonds]

    _bond = set(("C5", "O5"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C5", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C6", "O6"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C4", "C3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C1", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C4"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "O3"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("C3", "C2"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_infer_residue_connections():
    _man9 = bio.PDBParser().get_structure("MANNOSE9", base.MANNOSE9)
    bonds = bam.structural.infer_residue_connections(_man9)

    connections = [
        set((i.get_parent()._id[1], j.get_parent()._id[1])) for i, j in bonds
    ]
    bonds = [set((i.id, j.id)) for i, j in bonds]

    _received = len(bonds)
    _expected = 10
    assert _received == _expected, f"Expected {_expected} bonds, got {_received}"

    _bond = set(("O4", "C1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set(("O3", "C1"))
    _recieved = _bond in bonds
    _expected = True
    _what = f"for {_bond} in bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((2, 3))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((5, 6))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _bond = set((5, 8))
    _recieved = _bond in connections
    _expected = True
    _what = f"for {_bond} in connections"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_infer_residue_connections_triplet():
    _man9 = bio.PDBParser().get_structure("MANNOSE9", base.MANNOSE9)
    bonds = bam.structural.infer_residue_connections(_man9, triplet=True)
    _no_triplets = bam.structural.infer_residue_connections(_man9)

    assert len(bonds) == 2 * len(_no_triplets), "Not all triplets are found!"


def test_atom_neighborhood_basic():
    man = bam.molecule(MANNOSE)
    man.infer_bonds()
    mannose = man._AtomGraph

    _recieved = len(man.bonds)
    _expected = 24
    _what = "bonds"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    neighborhood = bam.structural.AtomNeighborhood(mannose)
    assert neighborhood is not None, "No neighborhood object is made..."

    _recieved = len(neighborhood.atoms)
    _expected = 24
    _what = "atoms"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    # Neigborhood is not designed for this anymore...
    # a = neighborhood.get_atom(1)
    # assert a is not None

    # a = neighborhood.get_atom("C1")
    # assert a is not None
    # assert not isinstance(a, list), "we get a list of C1s although there is only 1"


def test_atom_neighborhood_get():
    man = bam.molecule(MANNOSE)
    man.infer_bonds()
    mannose = man._AtomGraph
    neighborhood = bam.structural.AtomNeighborhood(mannose)

    c1 = next(atom for atom in man.get_atoms() if atom.id == "C1")

    _recieved = set(i.id for i in neighborhood.get_neighbors(c1))
    _expected = {"H1", "C2", "O1", "O5"}
    _what = "as n=1 neighbors of C1"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _recieved = set(i.id for i in neighborhood.get_neighbors(c1, 2))
    _n2 = {"HO1", "H2", "O2", "C3", "C5"}
    _expected.update(_n2)
    _what = "as n<=2 neighbors of C1"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    _recieved = set(i.id for i in neighborhood.get_neighbors(c1, 2, mode="at"))
    _expected = _n2
    _what = "as n==2 neighbors of C1"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_residue_neighborhood_basic():
    mannose = bam.Molecule.from_pdb(base.MANNOSE9)
    mannose.infer_bonds(restrict_residues=False)
    graph = mannose.make_residue_graph()

    neighborhood = bam.structural.ResidueNeighborhood(graph)
    assert neighborhood is not None, "No neighborhood object is made..."

    _recieved = len(neighborhood.residues)
    _expected = 11
    _what = "residues"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"

    # # Graphs and Neighborhoods are not designed for this anymore...
    # # nag2 = mannose.get_residue(2)  # first residue is labelled as NAG (resseq=2)
    # # a = neighborhood.get_residue(nag2)
    # # assert a is not None

    # man = mannose.get_residues("MAN", by="name")
    # a = neighborhood.get_residue(man)
    # assert a is not None
    # assert isinstance(a, list), "we expect a list of MANs!"

    # bma = mannose.get_residue("BMA")
    # a = neighborhood.get_residue(bma)
    # assert a is not None
    # assert not isinstance(
    #     a, list
    # ), "we expect a single residue of BMA since there is only one!"


def test_residue_neighborhood_get():
    mannose = bam.Molecule.from_pdb(base.MANNOSE9)
    mannose.infer_bonds(restrict_residues=False)
    graph = mannose.make_residue_graph()

    neighborhood = bam.structural.ResidueNeighborhood(graph)
    assert neighborhood is not None, "No neighborhood object is made..."

    # because the graph is not detailed there should be no
    # node and neighbors for C1 from BMA residue
    try:
        c1 = mannose.get_atom("C1", residue=4)
        neigs = neighborhood.get_neighbors(c1)
    except KeyError:
        pass
    else:
        raise AssertionError(
            "We should not be able to get neighbors of C1 from BMA residue!"
        )

    man = mannose.get_residues("MAN")
    neigs = neighborhood.get_neighbors(man)
    assert isinstance(neigs, list), f"Expected a list but received {type(neigs)}"

    _received = len(neigs)
    _expected = 8
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    bma = mannose.get_residue("BMA")
    neigs = neighborhood.get_neighbors(bma)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = neighborhood.get_neighbors(bma, 2)
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 7
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    neigs = neighborhood.get_neighbors(bma, 2, "at")
    assert isinstance(neigs, set), f"Expected a set but received {type(neigs)}"

    _received = len(neigs)
    _expected = 4
    assert _received == _expected, f"Expected {_expected} neighbors, got {_received}"

    _recieved = set(i.id[1] for i in neighborhood.get_neighbors(bma, 2, "at"))
    _expected = {2, 8, 6, 11}
    _what = "as n=2 neighbors of BMA"
    assert (
        _recieved == _expected
    ), f"Recieved {_recieved} {_what}, expected {_expected} {_what}!"


def test_compute_angle():
    mannose = bam.molecule(MANNOSE)
    mannose.infer_bonds()

    for triplet, angle in mannose.compute_angles().items():
        assert 90 < angle < 120, f"Angle {angle} is not in range 90-120°!"

    # top = bam.resources.get_default_topology()
    # man = top.get_residue("MAN")

    # _atom = "O5"  # some ref atom to get ICs for
    # ic, *_ = man.get_internal_coordinates(_atom, None, None, None, mode="partial")

    # missing = man.get_missing_atoms(mannose)
    # assert missing == [], f"There are missing atoms: {missing}"

    # refs = ic.get_reference_atoms(mannose)
    # assert len(refs) == 4, f"We got weird reference atoms: {refs}"

    # _true_angle = ic.bond_angle_123
    # _recieved = bam.structural.compute_angle(*refs[:-1])
    # _what = "° between 1-2-3"
    # assert _recieved == pytest.approx(
    #     _true_angle, 1e-3
    # ), f"Recieved {_recieved} {_what}, expected {_true_angle} {_what}!"

    # _true_angle = ic.bond_angle_234
    # _recieved = bam.structural.compute_angle(*refs[1:])
    # _what = "° between 2-3-4"
    # assert _recieved == pytest.approx(
    #     _true_angle, 1e-3
    # ), f"Recieved {_recieved} {_what}, expected {_true_angle} {_what}!"


def test_compute_dihedral():
    mannose = bam.molecule(MANNOSE)

    for quartet, dihedral in mannose.compute_dihedrals().items():
        assert -120 < dihedral < 120, f"Dihedral {dihedral} is not in range -120-120°!"

    # mannose = bam.utils.defaults.__bioPDBParser__.get_structure("MAN", base.MANNOSE)
    # mannose = next(mannose.get_residues())

    # top = bam.resources.get_default_topology()
    # man = top.get_residue("MAN")

    # _atom = "O5"  # some ref atom to get ICs for
    # ic, *_ = man.get_internal_coordinates(_atom, None, None, None, mode="partial")

    # missing = man.get_missing_atoms(mannose)
    # assert missing == [], f"There are missing atoms: {missing}"

    # refs = ic.get_reference_atoms(mannose)
    # assert len(refs) == 4, f"We got weird reference atoms: {refs}"

    # _true_dihedral = ic.dihedral
    # _recieved = bam.structural.compute_dihedral(*refs)
    # _what = "° between 1-2-3-4"
    # assert _recieved == pytest.approx(
    #     _true_dihedral, 1e-3
    # ), f"Recieved {_recieved} {_what}, expected {_true_dihedral} {_what}!"


def test_compute_triplets():
    bonds = [(1, 2), (1, 3), (2, 4), (3, 5)]
    triplets = bam.structural.compute_triplets(bonds, unique=False)
    _expected = set(((2, 1, 3), (3, 1, 2), (1, 2, 4), (4, 2, 1), (1, 3, 5), (5, 3, 1)))
    assert (
        set(triplets) == _expected
    ), f"Expected {len(_expected)} triplets, got {len(triplets)}"
    triplets = bam.structural.compute_triplets(bonds, unique=True)
    assert (
        len(set(triplets).intersection(_expected)) == 3
    ), "Unique triplets are not unique!"


def test_quartet_class():
    a = bam.structural.neighbors.Quartet(1, 2, 3, 4, False)
    b = bam.structural.neighbors.Quartet(1, 2, 3, 4, False)
    c = bam.structural.neighbors.Quartet(5, 3, 4, 6, True)

    assert a == b, "Quartets are not equal!"
    assert a != c, "Quartets are equal!"

    assert (1, 2, 3, 4) == a
    assert (1, 2, 3, 4, False) == a
    assert (1, 2, 3, 4, True) != a

    assert a[0] == 1


def test_compute_quartets():
    bonds = bonds = [(1, 2), (2, 3), (2, 4), (3, 5)]
    quartets = bam.structural.compute_quartets(bonds)

    _received = len(quartets)
    _expected = 3
    assert _received == _expected, f"Expected {_expected} quartets, got {_received}"
    assert sum(1 for i in quartets if i.improper) == 1, "Expected 1 improper quartet"
    assert sum(1 for i in quartets if not i.improper) == 2, "Expected 2 proper quartets"

    bonds = [(1, 2), (2, 3), (2, 4), (3, 5), (4, 6), (5, 7)]
    quartets = bam.structural.compute_quartets(bonds)

    _received = len(quartets)
    _expected = 6
    assert _received == _expected, f"Expected {_expected} quartets, got {_received}"

    # Quartet = bam.structural.neighbors.Quartet
    # assert Quartet(1, 2, 4, 6, False) in quartets
    # assert Quartet(1, 4, 2, 3, True) in quartets


def test_patcher_anchors():
    man1 = bam.Molecule.from_pdb(base.MANNOSE)
    man1.infer_bonds()
    man2 = deepcopy(man1)

    top = bam.get_default_topology()
    patch = top.get_patch("12aa")

    p = bam.structural.Patcher()
    p.target = man1
    p.source = man2
    p.patch = patch
    anchors = p.get_anchors()

    assert anchors[0] and anchors[1]
    assert anchors[0].id == "O2"
    assert anchors[1].id == "C1"


def test_patcher_anchors_2():
    bam.load_sugars()
    glc = bam.Molecule.from_compound("GLC")

    top = bam.get_default_topology()
    patch = top.get_patch("14bb")

    p = bam.structural.Patcher(True, True)
    p.target = glc
    p.source = glc
    p.patch = patch
    anchors = p.get_anchors()

    assert len(anchors) == 2
    assert anchors[0].id == "O4"
    assert anchors[1].id == "C1"
    bam.unload_sugars()


def test_patcher_two_man():
    man1 = bam.Molecule.from_pdb(base.MANNOSE)
    man1.infer_bonds()
    man2 = deepcopy(man1)

    man1.lock_all()
    man2.lock_all()

    top = bam.get_default_topology()
    patches = ("12aa", "12ab", "14bb")
    p = bam.structural.Patcher(copy_target=True, copy_source=True)
    for patch in patches:
        patch = top.get_patch(patch)

        assert len(man1.atoms) == len(man2.atoms) == 24
        assert len(man1.bonds) == len(man2.bonds) == 24

        _man1, _man2 = p.apply(patch, man1, man2)
        new = p.merge()

        v = new.draw()
        v.draw_edges(*new.locked_bonds, color="red")
        v.draw_edges(*new.bonds, color="cyan", linewidth=2)
        v.show()

        assert new is not man1 and new is not man2
        assert len(new.residues) == 2

        # check that the right amount of atoms are available
        # since the different patches remove different numbers of atoms
        # we simply check that there are now more but not quite double the atoms
        assert 1.5 * len(man1.atoms) < len(new.atoms) < 2 * len(man1.atoms)
        assert 1.5 * len(man1.bonds) < len(new.bonds) < 2 * len(man1.bonds)
        assert len(new.locked_bonds) < 2 * len(
            new.bonds
        )  # there is a new connection that should not be locked
        # however, we need x2 because the locked bonds are in "either way"
        # mode which means we get twice the bonds once as A-B, and once as B-A

        connection = new.infer_residue_connections(triplet=False)
        assert len(connection) == 1
        assert len(new.get_bonds(*connection[0])) == 1
        assert new.is_locked(*connection[0]) is False

        assert len(new.infer_residue_connections(triplet=True)) == 2

        # check that the atoms have been properly renumbered
        _seen_serials = set()
        for atom in new.atoms:
            assert atom.serial_number not in _seen_serials
            _seen_serials.add(atom.serial_number)

        # check that there are no major clashes
        # e.g. check that all atoms are at least 1Angstrom apart from each other
        for atom in new.atoms:
            for other in new.atoms:
                if atom == other:
                    continue
                assert np.sum(np.abs(atom.coord - other.coord)) > 1

        # check that there are now super weird angles
        # e.g. angles that are less than 100 or more than 130 degrees
        for angle in new.compute_angles().values():
            assert 100 < angle < 130


def test_patcher_multiple_man():
    bam.load_sugars()

    man1 = bam.Molecule.from_compound("MAN")
    man1.lock_all()

    man2 = man1.copy()
    man3 = man1.copy()
    man4 = man1.copy()

    top = bam.get_default_topology()

    orig_residues = len(man1.residues)
    orig_atoms = len(man1.atoms)

    p = bam.structural.Patcher(False, False)

    p.apply(top.get_patch("14bb"), man3, man4)
    man_34 = p.merge()

    p.apply(top.get_patch("14bb"), man_34, man2)
    man_23 = p.merge()

    p.apply(top.get_patch("12aa"), man1, man_23)
    new = p.merge()

    assert new is man1
    assert len(new.residues) == len(man1.residues) == 4 * orig_residues
    assert 3.5 * orig_atoms < len(man1.atoms) < 4 * orig_atoms

    # check that the atoms have been properly renumbered
    _seen_serials = set()
    for atom in man1.atoms:
        assert atom.serial_number not in _seen_serials
        _seen_serials.add(atom.serial_number)

    # check that there are no major clashes
    # e.g. check that all atoms are at least 1Angstrom apart from each other
    for atom in man1.atoms:
        for other in man1.atoms:
            if atom == other:
                continue
            assert np.sum(np.abs(atom.coord - other.coord)) > 1

    # check that there are now super weird angles
    # e.g. angles that are less than 100 or more than 130 degrees
    for angle in man1.compute_angles().values():
        assert 100 < angle < 130

    v = man1.draw()
    res_con = man1.infer_residue_connections(triplet=True)
    v.draw_edges(*res_con, color="limegreen", linewidth=3)
    v.show()

    g = man1.make_residue_graph()
    g.unlock_all()
    g.lock_centers()
    v = g.draw()
    v.draw_edges(*g.edges, color="limegreen", linewidth=2)
    v.draw_edges(*g._locked_edges, color="red", linewidth=2)
    v.show()


def test_keep_copy_patcher():
    bam.load_sugars()

    glc = bam.Molecule.from_compound("GLC")

    patcher = bam.structural.Patcher(copy_target=True, copy_source=True)
    patch = bam.get_default_topology().get_patch("12aa")

    patcher.apply(patch, glc, glc)
    new = patcher.merge()

    assert new is not glc
    assert len(new.atoms) == len(glc.atoms) * 2 - 3
    assert len(glc.atoms) == 24


def test_stitcher_two_glucose():
    bam.load_sugars()

    glc = bam.Molecule.from_compound("GLC")
    glc2 = bam.Molecule.from_compound("GLC")

    glc2.rotate_around_bond(6, 5, 68)
    glc2.rotate_around_bond(3, 4, 41)

    s = bam.structural.Stitcher(True, True)

    at_glc = "C1"
    at_glc2 = "O4"
    remove_on_glc = ("O1", "HO1")
    remove_on_glc2 = ("HO4",)

    old_glc_coords = np.array([at.coord for at in glc.atoms])
    old_glc2_coords = np.array([at.coord for at in glc2.atoms])

    new_glc, new_glc2 = s.apply(
        target=glc,
        source=glc2,
        target_removals=remove_on_glc,
        source_removals=remove_on_glc2,
        target_atom=at_glc,
        source_atom=at_glc2,
        optimization_steps=1e6,
    )

    new_glc_coords = np.array([at.coord for at in new_glc.atoms])
    new_glc2_coords = np.array([at.coord for at in new_glc2.atoms])

    assert new_glc is not glc
    assert new_glc2 is not glc2
    assert len(new_glc.atoms) == len(glc.atoms) - len(remove_on_glc)
    assert len(new_glc2.atoms) == len(glc2.atoms) - len(remove_on_glc2)

    r_glc_1 = glc.get_atom(remove_on_glc[0])
    r_glc_2 = glc.get_atom(remove_on_glc[1])
    r_glc = r_glc_1 if r_glc_1.serial_number < r_glc_2.serial_number else r_glc_2
    r_glc = r_glc.serial_number

    r_glc2_1 = glc2.get_atom(remove_on_glc2[0])
    r_glc2 = r_glc2_1.serial_number

    assert np.allclose(old_glc_coords[:r_glc, :], new_glc_coords[:r_glc, :])
    assert not np.allclose(old_glc2_coords[:r_glc2, :], new_glc2_coords[:r_glc2, :])

    final = s.merge()
    final.show()

    assert len(final.atoms) == len(glc.atoms) + len(glc2.atoms) - len(
        remove_on_glc
    ) - len(remove_on_glc2)
    assert len(final.residues) == 2
    assert len(final.bonds) == 46

    for angle in final.compute_angles().values():
        assert 90 < angle < 130

    for angle in final.compute_dihedrals().values():
        assert -180 < angle < 180

    _seen_indices = set()
    for atom in final.atoms:
        assert atom.serial_number not in _seen_indices
        _seen_indices.add(atom.serial_number)


def test_stitcher_three_glucose():
    bam.load_sugars()

    glc = bam.Molecule.from_compound("GLC")
    glc2 = bam.Molecule.from_compound("GLC")

    glc2.rotate_around_bond(6, 5, 68)
    glc2.rotate_around_bond(3, 4, 41)

    s = bam.structural.Stitcher(True, True)

    at_glc = "C1"
    at_glc2 = "O4"
    remove_on_glc = ("O1", "HO1")
    remove_on_glc2 = ("HO4",)

    s.apply(
        target=glc,
        source=glc2,
        target_removals=remove_on_glc,
        source_removals=remove_on_glc2,
        target_atom=at_glc,
        source_atom=at_glc2,
        optimization_steps=10,
    )

    new_glc = s.merge()

    s.apply(
        target=new_glc,
        source=glc,
        target_removals=remove_on_glc,
        source_removals=remove_on_glc2,
        target_atom=at_glc,
        source_atom=at_glc2,
        target_residue=2,
        optimization_steps=10,
    )
    final = s.merge()

    assert final is not glc and final is not new_glc
    assert len(final.residues) == 3
    assert len(final.atoms) == 3 * len(glc.atoms) - 2 * (
        len(remove_on_glc) + len(remove_on_glc2)
    )
    assert len(final.bonds) == 68

    for angle in final.compute_angles().values():
        assert 90 < angle < 130

    for angle in final.compute_dihedrals().values():
        assert -180 < angle < 180

    _seen_indices = set()
    for atom in final.atoms:
        assert atom.serial_number not in _seen_indices
        _seen_indices.add(atom.serial_number)

    final.show()


# def test_stitcher_two_glucose_root_atoms():
#     glc = bam.Molecule.from_compound("GLC")
#     glc2 = bam.Molecule.from_compound("GLC")

#     glc2.rotate_around_bond(6, 5, 68)
#     glc2.rotate_around_bond(3, 4, 41)

#     s = bam.structural.Stitcher(True, True)

#     at_glc = "C1"
#     at_glc2 = "O4"
#     remove_on_glc = ("O1", "HO1")
#     remove_on_glc2 = ("HO4",)

#     old_glc_coords = np.array([at.coord for at in glc.atoms])
#     old_glc2_coords = np.array([at.coord for at in glc2.atoms])

#     glc.set_root(at_glc)
#     glc2.set_root(at_glc2)

#     new_glc, new_glc2 = s.apply(
#         target=glc,
#         source=glc2,
#         target_removals=remove_on_glc,
#         source_removals=remove_on_glc2,
#         target_atom=None,
#         source_atom=None,
#         optimization_steps=1e6,
#     )

#     new_glc_coords = np.array([at.coord for at in new_glc.atoms])
#     new_glc2_coords = np.array([at.coord for at in new_glc2.atoms])

#     assert new_glc is not glc
#     assert new_glc2 is not glc2
#     assert len(new_glc.atoms) == len(glc.atoms) - len(remove_on_glc)
#     assert len(new_glc2.atoms) == len(glc2.atoms) - len(remove_on_glc2)

#     r_glc_1 = glc.get_atom(remove_on_glc[0])
#     r_glc_2 = glc.get_atom(remove_on_glc[1])
#     r_glc = r_glc_1 if r_glc_1.serial_number < r_glc_2.serial_number else r_glc_2
#     r_glc = r_glc.serial_number

#     r_glc2_1 = glc2.get_atom(remove_on_glc2[0])
#     r_glc2 = r_glc2_1.serial_number

#     assert np.allclose(old_glc_coords[:r_glc, :], new_glc_coords[:r_glc, :])
#     assert not np.allclose(old_glc2_coords[:r_glc2, :], new_glc2_coords[:r_glc2, :])

#     final = s.merge()
#     assert len(final.atoms) == len(glc.atoms) + len(glc2.atoms) - len(
#         remove_on_glc
#     ) - len(remove_on_glc2)
#     assert len(final.residues) == 2
#     assert len(final.bonds) == 46

#     for angle in final.compute_angles().values():
#         assert 90 < angle < 130

#     for angle in final.compute_dihedrals().values():
#         assert -180 < angle < 180

#     _seen_indices = set()
#     for atom in final.atoms:
#         assert atom.serial_number not in _seen_indices
#         _seen_indices.add(atom.serial_number)

#     final.show()


def test_patch_and_stich():
    bam.load_sugars()

    glc = bam.Molecule.from_compound("GLC")
    man = bam.Molecule.from_compound("MAN")

    # ------------------------------------------
    # using the built-in patcher within molecule
    # ------------------------------------------
    # make a cellulose
    glc.repeat(4, "14bb")

    # make a mannose chain
    man.repeat(3, "16ab")

    # ------------------------------------------
    # now stitch them together
    # ------------------------------------------
    stitcher = bam.structural.Stitcher()
    stitcher.apply(
        target=glc,
        source=man,
        target_removals=("HO3",),
        source_removals=("O1", "HO1"),
        target_atom="O3",
        source_atom="C1",
        target_residue=2,
        source_residue=1,
    )
    final = stitcher.merge()
    final.show()

    # ------------------------------------------

    assert final is glc
    assert len(final.residues) == 7

    for angle in final.compute_angles().values():
        assert 90 < angle < 130

    for angle in final.compute_dihedrals().values():
        assert -180 < angle < 180

    _seen_indices = set()
    for atom in final.atoms:
        assert atom.serial_number not in _seen_indices
        _seen_indices.add(atom.serial_number)


def test_relabel_Hydrogens():
    bam.load_sugars()

    ref = bam.Molecule.from_compound("GLC")
    mol = bam.Molecule.from_pubchem("D-glucose")

    bam.structural.relabel_hydrogens(mol)

    hydrogens = (a for a in mol.get_atoms() if a.element == "H")
    for h in hydrogens:
        assert h.id.startswith("H")
        assert h.id[1].isdigit() or h.id[1] in ("C", "O", "N", "S")
        assert h.id[-1].isdigit()


def test_autolabel():
    bam.load_sugars()

    ref = bam.Molecule.from_compound("GLC")
    mol = bam.Molecule.from_pubchem("GLC")
    mol.autolabel()
    refs = set((i.id) for i in ref.get_atoms())
    mols = set((i.id) for i in mol.get_atoms())
    assert refs == mols

    for atom in mol.get_atoms():
        assert atom.id == atom.full_id[-1][0]


def test_autolabel2():
    bam.unload_all_compounds()
    bam.load_small_molecules()
    bam.load_sugars()
    assert len(bam.get_default_compounds()) == 3180  # small + sugars
    assert bam.has_compound("CH3") == True, "CH3 is not a compound!"
    mol = bam.molecule("CH3")
    assert mol is not None
    mol = bam.Molecule.from_compound("CH4", by="formula")
    assert isinstance(mol, bam.Molecule)
    mol.autolabel()
    assert set(i.id for i in mol.get_atoms()) == set(("C1", "H11", "H12", "H13", "H14"))
