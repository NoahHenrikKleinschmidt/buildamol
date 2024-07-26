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
MANNOSE = bio.PDBParser().get_structure("MAN", base.MANPDB)
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
#     _man = bio.PDBParser().get_structure("MAN9", base.MAN9PDB)

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
    _man9 = bio.PDBParser().get_structure("MANNOSE9", base.MAN9PDB)
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
    _man9 = bio.PDBParser().get_structure("MANNOSE9", base.MAN9PDB)
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
    mannose = bam.Molecule.from_pdb(base.MAN9PDB)
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
    mannose = bam.Molecule.from_pdb(base.MAN9PDB)
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

    # mannose = bam.utils.defaults.__bioPDBParser__.get_structure("MAN", base.MANPDB)
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
    man1 = bam.Molecule.from_pdb(base.MANPDB)
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
    man1 = bam.Molecule.from_pdb(base.MANPDB)
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

        if base.ALLOW_VISUAL:
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

    if base.ALLOW_VISUAL:
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
    if base.ALLOW_VISUAL:
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

    if base.ALLOW_VISUAL:
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
    if base.ALLOW_VISUAL:
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


def test_rotate_molecule():
    mol = bam.Molecule.from_compound("GLC")

    old_coords = np.array([a.coord for a in mol.get_atoms()])

    d = mol.draw()
    bam.structural.rotate_molecule(mol, 90, np.array([1, 0, 0]))
    d.draw_edges(*mol.get_bonds(), color="red", linewidth=2)

    bam.structural.rotate_molecule(mol, 30, np.array([0, 1, 0.2]))
    d.draw_edges(*mol.get_bonds(), color="green", linewidth=2)

    bam.structural.rotate_molecule(mol, 45, np.array([4.5, 2.3, 1.4]))
    d.draw_edges(*mol.get_bonds(), color="blue", linewidth=2)

    new_coords = np.array([a.coord for a in mol.get_atoms()])

    assert not np.allclose(old_coords, new_coords)

    d.show()


def test_flip():
    mol = bam.Molecule.from_compound("GLC")
    old_coords = np.array([a.coord for a in mol.get_atoms()])
    d = mol.draw()
    bam.structural.flip_molecule(mol, [0, 0, 1], center=mol.center_of_geometry)
    d.draw_vector("axis", [0, 0, 0], [0, 0, 1], color="orange")
    d.draw_edges(*mol.get_bonds(), color="red", linewidth=2)
    bam.structural.flip_molecule(mol, [0, 1, 0])
    d.draw_edges(*mol.get_bonds(), color="green", linewidth=2)
    bam.structural.flip_molecule(mol, [1, 0, 0])
    d.draw_edges(*mol.get_bonds(), color="blue", linewidth=2)
    new_coords = np.array([a.coord for a in mol.get_atoms()])
    assert not np.allclose(old_coords, new_coords)
    d.show()


def test_infer_hydrogens_glucose():
    mol = bam.Molecule.from_compound("GLC")
    ref = mol.copy()

    mol.remove_atoms("H1", "H61", "H62", "HO1")

    # v = mol.draw()

    hydrogenator = bam.structural.infer.Hydrogenator()
    hydrogenator.infer_hydrogens(mol)

    # v.draw_points(
    #     mol.get_coords("H1", "H61", "H62", "HO1"), colors="orange"
    #     )

    # v.show()

    assert len(mol.get_atoms("H1")) == 1
    assert len(mol.get_atoms("H61")) == 1
    assert len(mol.get_atoms("H62")) == 1
    assert len(mol.get_atoms("HO1")) == 1

    assert mol.get_atom("HO1") in mol.get_neighbors(mol.get_atom("O1"))
    assert mol.get_atom("H1") in mol.get_neighbors(mol.get_atom("C1"))

    new_coords = mol.get_coords("H1", "H61", "H62", "HO1")
    ref_coords = ref.get_coords("H1", "H61", "H62", "HO1")

    # evaluate that each of the inferred atoms is close
    # by one of the references. (this is primarily because
    # H61/H62 may have swapped labels when inferred)
    d = bam.structural.cdist(new_coords, ref_coords)

    assert d[0].min() < 0.5
    assert d[1].min() < 0.5
    assert d[2].min() < 0.5
    assert (
        d[3].min() < 1.1
    )  # this is the HO1 which can freely rotate around the C-O axis and therefore is not necessarily in the same position


def test_infer_hydrogens_glucose_all():
    mol = bam.Molecule.from_compound("GLC")
    ref = mol.copy()

    mol.remove_atoms(*mol.get_atoms("H", by="element"))

    # v = mol.draw()

    hydrogenator = bam.structural.infer.Hydrogenator()
    hydrogenator.infer_hydrogens(mol, bond_length=1.05)

    # v.draw_points(mol.get_coords("H", by="element"), colors="orange")

    # v.show()

    assert len(mol.get_atoms("H1")) == 1
    assert len(mol.get_atoms("H61")) == 1
    assert len(mol.get_atoms("H62")) == 1
    assert len(mol.get_atoms("HO1")) == 1

    assert mol.get_atom("HO1") in mol.get_neighbors(mol.get_atom("O1"))
    assert mol.get_atom("H1") in mol.get_neighbors(mol.get_atom("C1"))

    new_coords = mol.get_coords("H1", "H61", "H62", "HO1")
    ref_coords = ref.get_coords("H1", "H61", "H62", "HO1")

    # evaluate that each of the inferred atoms is close
    # by one of the references. (this is primarily because
    # H61/H62 may have swapped labels when inferred)
    d = bam.structural.cdist(new_coords, ref_coords)

    assert d[0].min() < 0.1
    assert d[1].min() < 0.1
    assert d[2].min() < 0.1
    assert (
        d[3].min() < 2.1
    )  # this is the HO1 which can freely rotate around the C-O axis and therefore is not necessarily in the same position


def test_infer_hydrogens_tyrosine_all():
    bam.load_amino_acids()
    mol = bam.Molecule.from_compound("TYR")
    ref = mol.copy()

    mol.remove_atoms(*mol.get_atoms("H", by="element"))

    # v = mol.draw()

    hydrogenator = bam.structural.infer.Hydrogenator()
    hydrogenator.infer_hydrogens(mol, 1.05)

    # v.draw_points(mol.get_coords(
    #     "H", by="element"
    # ), colors="orange")

    # v.show()

    ref_coords = ref.get_coords("HD2", "HB2", "HB3", "HXT")
    new_coords = mol.get_coords("HD2", "HB1", "HB2", "HOXT")

    d = bam.structural.cdist(new_coords, ref_coords)
    assert d[0].min() < 0.1
    assert d[1].min() < 0.1
    assert d[2].min() < 0.1
    assert (
        d[3].min() < 2.1
    )  # the OH may freely rotate and be therefore not in the same position


def test_geometry_tetrahedral():
    bam.load_small_molecules()
    mol = bam.get_compound("CH4")
    assert mol is not None

    tetrahedron = bam.structural.geometry.Tetrahedral()
    assert tetrahedron is not None

    C = mol.get_atom("C", by="element")
    Hs = mol.get_atoms("H", by="element")

    # make just from one central atom and one H
    coords = tetrahedron.make_coords(C, Hs[0])

    v = mol.draw()
    v.draw_points(coords, colors="orange")
    assert np.all(coords[0] == C.coord)
    assert np.all(coords[1] == Hs[0].coord)
    theta = bam.structural.angle_between(coords[1], coords[0], coords[4])
    assert 109 < theta < 110
    theta = bam.structural.angle_between(coords[2], coords[0], coords[4])
    assert 109 < theta < 110
    theta = bam.structural.angle_between(coords[3], coords[0], coords[4])
    assert 109 < theta < 110

    # make from central atom and 2 neighbors
    coords = tetrahedron.make_coords(C, *Hs[:2])

    v.draw_points(coords, colors="red")
    if base.ALLOW_VISUAL:
        v.show()

    assert np.all(coords[0] == C.coord)
    assert np.all(coords[1] == Hs[0].coord)
    assert np.all(coords[2] == Hs[1].coord)

    theta = bam.structural.angle_between(coords[1], coords[0], coords[2])
    assert 95 < theta < 130
    theta = bam.structural.angle_between(coords[1], coords[0], coords[3])
    assert 95 < theta < 130
    theta = bam.structural.angle_between(coords[1], coords[0], coords[4])
    assert 95 < theta < 130
    theta = bam.structural.angle_between(coords[2], coords[0], coords[4])
    assert 95 < theta < 130
    theta = bam.structural.angle_between(coords[3], coords[0], coords[4])
    assert 95 < theta < 130


def test_tetrahedral_applied_glucose():
    bam.load_sugars()
    mol = bam.get_compound("GLC")
    assert mol is not None

    C = mol.get_atom("C1")
    C2 = mol.get_atom("C2")
    O5 = mol.get_atom("O5")
    O1 = mol.get_atom("O1")
    H1 = mol.get_atom("H1")

    tetrahedron = bam.structural.geometry.Tetrahedral()

    # reconstruct coordinates
    coords = tetrahedron.make_coords(C, C2, O1)

    coords[3] = bam.structural.adjust_distance(coords[0], coords[3], C - O5)
    coords[4] = bam.structural.adjust_distance(coords[0], coords[4], C - H1)

    if base.ALLOW_VISUAL:
        v = mol.draw()
        v.draw_points(coords, colors="orange")
        v.show()

    assert np.all(coords[0] == C.coord)
    assert np.all(coords[1] == C2.coord)
    assert np.all(coords[2] == O1.coord)

    assert np.abs(coords[3] - O5.coord).sum() < 0.01
    assert np.abs(coords[4] - H1.coord).sum() < 0.01


def test_planar():
    mol = bam.molecule("ethene").autolabel()
    assert mol is not None

    planar = bam.structural.geometry.TriangularPlanar()

    C1 = mol.get_atom("C1")
    C2 = mol.get_atom("C2")
    H1 = mol.get_atom("H11")
    H2 = mol.get_atom("H12")

    coords = planar.make_coords(C1, C2, H1)

    assert np.all(coords[0] == C1.coord)
    assert np.all(coords[1] == C2.coord)
    assert np.all(coords[2] == H1.coord)
    coords[3] = bam.structural.adjust_distance(coords[0], coords[3], C1 - H2)
    assert np.abs(coords[3] - H2.coord).sum() < 0.01

    mol.transpose(
        [0, 2, 12],
        45.6,
        np.array([0.2, 1, 0.3]) / np.linalg.norm(np.array([0.2, 1, 0.3])),
    )

    coords = planar.make_coords(C1, C2, H1)
    assert np.all(coords[0] == C1.coord)
    assert np.all(coords[1] == C2.coord)
    assert np.all(coords[2] == H1.coord)
    coords[3] = bam.structural.adjust_distance(coords[0], coords[3], C1 - H2)
    assert np.abs(coords[3] - H2.coord).sum() < 0.01

    # v = mol.draw()
    # v.draw_points(coords, colors="orange")
    # v.show()

    coords = planar.make_coords(C1, C2)
    assert np.all(coords[0] == C1.coord)
    assert np.all(coords[1] == C2.coord)
    assert bam.structural.angle_between(coords[0], coords[1], coords[2]) - 120 < 0.5
    assert bam.structural.angle_between(coords[0], coords[1], coords[3]) - 120 < 0.5
    assert bam.structural.angle_between(coords[2], coords[0], coords[3]) - 120 < 0.5


def test_linear():
    mol = bam.molecule("C#C")[0]
    mol.autolabel()
    assert mol is not None

    linear = bam.structural.geometry.Linear()

    C1 = mol.get_atom("C1")
    C2 = mol.get_atom("C2")
    H1 = mol.get_atom("H1")
    H2 = mol.get_atom("H2")

    coords = linear.make_coords(C1, C2)

    assert np.all(coords[0] == C1.coord)
    assert np.all(coords[1] == C2.coord)
    assert bam.structural.angle_between(coords[0], coords[1], coords[2]) - 180 < 0.5

    coords[2] = bam.structural.adjust_distance(coords[0], coords[2], C1 - H1)
    assert np.abs(coords[2] - H1.coord).sum() < 0.01

    mol.transpose(
        [0, 2, 12],
        45.6,
        np.array([0.2, 1, 0.3]) / np.linalg.norm(np.array([0.2, 1, 0.3])),
    )

    coords = linear.make_coords(C1, C2)

    assert np.all(coords[0] == C1.coord)
    assert np.all(coords[1] == C2.coord)
    # assert bam.structural.angle_between(coords[0], coords[1], coords[2]) - 180 < 0.5

    # coords[2] = bam.structural.adjust_distance(coords[0], coords[2], C1 - H1)
    # assert np.abs(coords[2] - H1.coord).sum() < 0.01


def test_bipyramid():

    mol = bam.Molecule.empty()
    res = bam.Residue("PCL5")
    mol.add_residues(res)

    bipyramid = bam.structural.geometry.TrigonalBipyramidal()

    P = bam.Atom.new("P")
    C1 = bam.Atom.new("Cl", [0, 0, 1])

    # atoms without coordinates yet
    C2 = bam.Atom.new("Cl")
    C3 = bam.Atom.new("Cl")
    C4 = bam.Atom.new("Cl")
    C5 = bam.Atom.new("Cl")

    mol.add_atoms(P, C1, C2, C3, C4, C5)

    coords = bipyramid.make_coords(P, C1, direction="planar")
    assert np.all(coords[0] == P.coord)
    assert np.all(coords[1] == C1.coord)

    assert bam.structural.angle_between(coords[1], coords[0], coords[2]) - 120 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[3]) - 120 < 0.5

    assert bam.structural.angle_between(coords[1], coords[0], coords[4]) - 90 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[5]) - 90 < 0.5

    C2.coord = coords[2]
    C3.coord = coords[3]
    C4.coord = coords[4]
    C5.coord = coords[5]

    v = mol.draw()
    coords = bipyramid.make_coords(P, C2, C4, direction="mixed")

    v.draw_points(coords, colors="orange")
    if base.ALLOW_VISUAL:
        v.show()

    assert np.all(coords[0] == P.coord)
    assert np.all(coords[1] == C2.coord)
    assert np.all(coords[2] == C4.coord)

    assert bam.structural.angle_between(coords[1], coords[0], coords[4]) - 120 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[3]) - 120 < 0.5

    assert bam.structural.angle_between(coords[1], coords[0], coords[2]) - 90 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[5]) - 90 < 0.5

    # now test with non-specified directionality to see if it could infer them
    direction = bam.structural.geometry._infer_point_relations(
        [i.coord for i in (P, C4, C5)], bipyramid.angle
    )
    assert direction == "axial"
    direction = bam.structural.geometry._infer_point_relations(
        [i.coord for i in (P, C2, C3)], bipyramid.angle
    )
    assert direction == "planar"
    direction = bam.structural.geometry._infer_point_relations(
        [i.coord for i in (P, C2, C4)], bipyramid.angle
    )
    assert direction == "mixed"
    direction = bam.structural.geometry._infer_point_relations(
        [i.coord for i in (P, C4, C2)], bipyramid.angle
    )
    assert direction == "mixed"
    direction = bam.structural.geometry._infer_point_relations(
        [i.coord for i in (P, C2, C5)], bipyramid.angle
    )
    assert direction == "mixed"


def test_octahedral():
    mol = bam.Molecule.new(None, "PF6")

    octahedral = bam.structural.geometry.Octahedral()

    P = bam.Atom.new("P")
    F1 = bam.Atom.new("F", [1, 0, 0])
    F2 = bam.Atom.new("F", [0, 1, 0])
    F3 = bam.Atom.new("F")
    F4 = bam.Atom.new("F")
    F5 = bam.Atom.new("F")
    F6 = bam.Atom.new("F")

    mol.add_atoms(P, F1, F2, F3, F4, F5, F6)

    coords = octahedral.make_coords(P, F1, F2, direction="planar", length=1)
    assert np.all(coords[0] == P.coord)
    assert np.all(coords[1] == F1.coord)
    assert np.all(coords[2] == F2.coord)

    assert bam.structural.angle_between(coords[1], coords[0], coords[2]) - 90 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[3]) - 180 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[4]) - 90 < 0.5
    assert bam.structural.angle_between(coords[1], coords[0], coords[5]) - 90 < 0.5

    assert bam.structural.angle_between(coords[2], coords[0], coords[3]) - 90 < 0.5
    assert bam.structural.angle_between(coords[2], coords[0], coords[4]) - 180 < 0.5
    assert bam.structural.angle_between(coords[2], coords[0], coords[5]) - 90 < 0.5

    F3.coord = coords[3]
    F4.coord = coords[4]
    F5.coord = coords[5]
    F6.coord = coords[6]

    # v = mol.draw()

    mol.rotate(23, "x").rotate(45, "y").move([1, 2, 3])

    coords = octahedral.make_coords(P, F1, F5, direction="mixed", length=1)

    # v.draw_points(coords, colors="orange")

    # v.show()

    assert np.all(coords[0] == P.coord)
    assert np.all(coords[1] == F1.coord)
    assert np.all(coords[2] == F5.coord)

    assert np.abs(F2.coord - coords[3]).sum() < 0.1
    assert np.abs(F3.coord - coords[4]).sum() < 0.1
    assert np.abs(F4.coord - coords[5]).sum() < 0.1
    assert np.abs(F6.coord - coords[6]).sum() < 0.1


def test_geometries_make_from_one():
    for geometry in (
        bam.structural.geometry.Linear(),
        bam.structural.geometry.TrigonalPlanar(),
        bam.structural.geometry.TrigonalBipyramidal(),
        bam.structural.geometry.Octahedral(),
        bam.structural.geometry.Tetrahedral(),
    ):
        atom = bam.Atom.new("C")
        new_coords = geometry.make_coords(atom)
        assert not np.isnan(
            new_coords
        ).any(), f"geometry {geometry.__class__.__name__} failed"
        assert not np.isinf(
            new_coords
        ).any(), f"geometry {geometry.__class__.__name__} failed"


def test_tetrahedral_apply():
    bam.load_small_molecules()
    mol = bam.get_compound("CH4").autolabel()

    # v = mol.draw()

    tetra = bam.structural.geometry.Tetrahedral()

    old = mol.get_coords()

    atoms = [mol.get_atom("C1"), *mol.get_atoms("H", by="element")]
    tetra.apply(atoms)
    new = mol.get_coords()

    # v.draw_points(new, ids=[i.id for i in atoms], colors="green")
    # v.show()

    assert ((old - new) ** 2).sum() < 0.1


def test_geometry_fill_hydrogens_from_one():
    for geometry in (
        bam.structural.geometry.Linear(),
        bam.structural.geometry.TrigonalPlanar(),
        bam.structural.geometry.TrigonalBipyramidal(),
        bam.structural.geometry.Octahedral(),
        bam.structural.geometry.Tetrahedral(),
    ):
        atom = bam.Atom.new("C")

        atoms, bonds = geometry.fill_hydrogens(atom)
        assert len(atoms) == geometry.size

        mol = bam.Molecule.new(None, "NEW")
        mol.add_atoms(*atoms)
        mol.infer_bonds()

        if base.ALLOW_VISUAL:
            mol.show()


def test_geometry_fill_hydrogens_from_two():
    for geometry in (
        bam.structural.geometry.Linear(),
        bam.structural.geometry.TrigonalPlanar(),
        bam.structural.geometry.TrigonalBipyramidal(),
        bam.structural.geometry.Octahedral(),
        bam.structural.geometry.Tetrahedral(),
    ):
        atom = bam.Atom.new("C")
        atom2 = bam.Atom.new("C", [1, 0, 0])

        atoms, bonds = geometry.fill_hydrogens(atom, atom2, direction="planar")
        assert len(atoms) == geometry.size

        mol = bam.Molecule.new(None, "NEW")
        mol.add_atoms(atoms)
        mol.add_bonds(bonds)

        if base.ALLOW_VISUAL:
            mol.show()


def test_change_element_add_hydrogens():
    mol = bam.Molecule.from_compound("GLC")

    old_neighbors = mol.get_neighbors("O1")
    old_atom = mol.get_atom("O1")

    bam.structural.infer.change_element(mol.get_atom("O1"), "N", mol)

    new = mol.get_atom("N1")
    assert old_atom is new

    assert mol.get_atom("N1").element == "N"
    assert mol.get_atom("N1").id == "N1"
    assert len(mol.get_neighbors("N1")) == len(old_neighbors) + 1

    old_neighbors = mol.get_neighbors("O5")

    old = mol.get_atom("O5")
    bam.structural.infer.change_element(mol.get_atom("O5"), "N", mol)
    new = mol.get_atom("N5")
    assert old is new

    assert mol.get_atom("N5").element == "N"
    assert mol.get_atom("N5").id == "N5"
    assert len(mol.get_neighbors("N5")) == len(old_neighbors) + 1

    if base.ALLOW_VISUAL:
        mol.show()


def test_change_element_remove_hydrogens():
    mol = bam.Molecule.from_compound("GLC")

    old_neighbors = mol.get_neighbors("O1")

    old = mol.get_atom("O1")
    bam.structural.infer.change_element(mol.get_atom("O1"), "N", mol)
    new = mol.get_atom("N1")
    assert old is new

    assert mol.get_atom("N1").element == "N"
    assert mol.get_atom("N1").id == "N1"
    assert len(mol.get_neighbors("N1")) == len(old_neighbors) + 1

    old_neighbors = mol.get_neighbors("N1")
    old = mol.get_atom("N1")
    # now change back again
    bam.structural.change_element(mol.get_atom("N1"), "O", mol)
    new = mol.get_atom("O1")
    assert old is new

    assert mol.get_atom("O1").element == "O"
    assert mol.get_atom("O1").id == "O1"
    assert len(mol.get_neighbors("O1")) == len(old_neighbors) - 1

    if base.ALLOW_VISUAL:
        mol.show()


def test_change_bond_order_add_hydrogens():
    mol = bam.Molecule.from_smiles("CC=NC")
    mol.autolabel()

    C1 = mol.get_atom("C1")
    N3 = mol.get_atom("N3")

    old_C_neighbors = mol.get_neighbors(C1)
    old_N_neighbors = mol.get_neighbors(N3)

    # change the double bond to a single bond
    bam.structural.change_bond_order(mol, C1, N3, 1)

    assert len(mol.get_neighbors(C1)) == len(old_C_neighbors) + 1
    assert len(mol.get_neighbors(N3)) == len(old_N_neighbors) + 1

    if base.ALLOW_VISUAL:
        mol.show()


def test_change_bond_order_remove_hydrogens():
    mol = bam.Molecule.from_smiles("CCNC")
    mol.autolabel()

    C1 = mol.get_atom("C1")
    N3 = mol.get_atom("N3")

    old_C_neighbors = mol.get_neighbors(C1)
    old_N_neighbors = mol.get_neighbors(N3)

    # change the double bond to a single bond
    bam.structural.change_bond_order(mol, C1, N3, 2)

    assert len(mol.get_neighbors(C1)) == len(old_C_neighbors) - 1
    assert len(mol.get_neighbors(N3)) == len(old_N_neighbors) - 1

    if base.ALLOW_VISUAL:
        mol.show()
    mol.to_pdb("mol.pdb")


def test_find_equatorial_hydrogens():
    bam.load_small_molecules()
    mol = bam.molecule("cyclohexane")

    v = mol.draw()
    H = bam.structural.infer.find_equatorial_hydrogens(mol)
    v.draw_atoms(*H, colors="orange")

    for C in mol.get_atoms("C", by="element"):
        H = bam.structural.infer.get_equatorial_hydrogen_neighbor(mol, C)
        v.draw_atom(H, color="purple")

    if base.ALLOW_VISUAL:
        v.show()


def test_find_axial_hydrogens():
    bam.load_small_molecules()
    mol = bam.molecule("cyclohexane")

    v = mol.draw()
    H = bam.structural.infer.find_axial_hydrogens(mol)
    v.draw_atoms(*H, colors="orange")

    for C in mol.get_atoms("C", by="element"):
        H = bam.structural.infer.get_axial_hydrogen_neighbor(mol, C)
        v.draw_atom(H, color="purple")

    if base.ALLOW_VISUAL:
        v.show()


def test_find_equatorial_hydrogens2():
    bam.load_sugars()
    mol = bam.molecule("GLC")

    v = mol.draw()

    H = bam.structural.infer.find_equatorial_hydrogens(mol)

    v.draw_atoms(*H, colors="orange")

    for carbon in mol.get_atoms("C", by="element"):
        H = bam.structural.infer.get_equatorial_hydrogen_neighbor(mol, carbon)
        if H:
            v.draw_atom(H, color="purple")

    if base.ALLOW_VISUAL:
        v.show()


def test_find_axial_hydrogens2():
    bam.load_sugars()
    mol = bam.molecule("GLC")

    v = mol.draw()

    H = bam.structural.infer.find_axial_hydrogens(mol)

    v.draw_atoms(*H, colors="orange")

    for carbon in mol.get_atoms("C", by="element"):
        H = bam.structural.infer.get_axial_hydrogen_neighbor(mol, carbon)
        if H:
            v.draw_atom(H, color="purple")

    if base.ALLOW_VISUAL:
        v.show()


def test_get_left_and_right_hydrogen():
    mol = bam.read_smiles("CCC(=O)O")
    mol.autolabel()

    center = mol.get_atom("C2")
    left_hydrogen = bam.structural.infer.get_left_hydrogen(mol, center)
    assert left_hydrogen is not None

    right_hydrogen = bam.structural.infer.get_right_hydrogen(mol, center)
    assert right_hydrogen is not None

    if base.ALLOW_VISUAL:
        v = mol.draw()
        v.draw_atom(left_hydrogen, color="orange")
        v.draw_atom(right_hydrogen, color="pink")
        v.show()


def test_superimpose_bonds():
    mol = bam.Molecule.from_compound("GLC")

    mol2 = mol.copy()
    mol2.transpose([0, 2, 12], 45, [1, 0, 0])

    v = mol.draw() + mol2.draw(atoms=False, line_color="red")

    a, b = "C2", "C3"
    bond1 = tuple(i.coord for i in mol.get_bond(a, b))
    bond2 = tuple(i.coord for i in mol2.get_bond(a, b))

    new_coords = bam.structural.superimpose_points(mol2.get_coords(), bond2, bond1)

    for idx, atom in enumerate(mol2.get_atoms()):
        atom.coord = new_coords[idx]

    v += mol2.draw(atoms=False, line_color="green")

    if base.ALLOW_VISUAL:
        v.show()


def test_superimpose_triplet():
    mol = bam.Molecule.from_compound("GLC")

    mol2 = mol.copy()
    mol2.transpose([0, 2, 12], 45, [1, 0, 0])

    v = mol.draw() + mol2.draw(atoms=False, line_color="red")

    points = "C2", "C3", "C1"
    points1 = tuple(i.coord for i in (mol.get_atom(i) for i in points))
    points2 = tuple(i.coord for i in (mol2.get_atom(i) for i in points))

    new_coords = bam.structural.superimpose_points(mol2.get_coords(), points2, points1)

    for idx, atom in enumerate(mol2.get_atoms()):
        atom.coord = new_coords[idx]

    v += mol2.draw(atoms=False, line_color="green")

    if base.ALLOW_VISUAL:
        v.show()


def test_infer_reactivity_atoms_from_groups():
    bam.load_amino_acids()
    mol = bam.Molecule.from_compound("TYR")

    carboxyl = bam.structural.groups.carboxyl
    amine = bam.structural.groups.amine

    link = bam.Linkage.from_functional_groups(mol, carboxyl, mol, amine)
    out = mol % link + mol
    out.show()


def test_infer_bond_orders():
    mols = [
        "2-butene",
        "NAG",
        "TYR",
        "acetylbenzene",
    ]
    for mol in mols:
        mol = bam.molecule(mol)
        for bond in mol.get_bonds():
            mol.set_bond_order(*bond, 1)

        bam.structural.infer_bond_orders(mol)
        assert any(bond.is_double() for bond in mol.get_bonds()), (
            "No double bonds found in molecule" + mol.id
        )
        if base.ALLOW_VISUAL:
            mol.show()


def test_match_aromatic():
    aromatic = bam.structural.groups.aromatic

    bam.load_small_molecules()
    bam.load_amino_acids()

    for mol in ["TYR", "BNZ", "acetophenone"]:
        mol = bam.molecule(mol)

        for bond in mol.bonds:
            mol.set_bond_order(*bond, 1)

        atoms = mol.atoms

        matches = aromatic.find_matches(mol, atoms)
        if len(matches) == 0:
            raise ValueError("No matches found")

        aromatic.apply_connectivity(mol, atoms)

        if base.ALLOW_VISUAL:
            mol.show()


def test_plane_from_points():
    from buildamol.extensions import polymers

    c = polymers.cyclic_alkane(10)

    vec = bam.structural.plane_of_points(c.get_coords())

    v = c.draw()
    v.draw_vector("vec", c.center_of_geometry, c.center_of_geometry + vec, color="red")
    if base.ALLOW_VISUAL:
        v.show()


def test_adjust_protonation_hydroxyl():
    mol = bam.read_smiles("OC")
    mol.autolabel()

    n_hydrogens = len(mol.get_atoms("H", by="element"))

    O = mol.get_atom("O", by="element")

    # add one proton to the oxygen
    bam.structural.adjust_protonation(mol, O, 1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens + 1
    mol.to_rdkit()

    # remove one proton from the oxygen
    bam.structural.adjust_protonation(mol, O, 0)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens
    mol.to_rdkit()

    # remove another proton from the oxygen
    bam.structural.adjust_protonation(mol, O, -1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens - 1
    mol.to_rdkit()


def test_adjust_protonation_carbonyl():
    mol = bam.read_smiles("C(=O)C")
    mol.autolabel()

    n_hydrogens = len(mol.get_atoms("H", by="element"))

    O = mol.get_atom("O", by="element")

    # add one proton to the oxygen
    bam.structural.adjust_protonation(mol, O, 1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens + 1
    mol.to_rdkit()

    # remove one proton from the oxygen
    bam.structural.adjust_protonation(mol, O, 0)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens
    mol.to_rdkit()

    # remove another proton from the oxygen
    try:
        bam.structural.adjust_protonation(mol, O, -1)
    except ValueError:
        pass
    else:
        raise RuntimeError("Should not be able to remove a proton from a carbonyl")


def test_adjust_protonation_amino():
    mol = bam.read_smiles("C(=O)N")
    mol.autolabel()

    n_hydrogens = len(mol.get_atoms("H", by="element"))

    N = mol.get_atom("N1")

    # add one proton to the nitrogen
    bam.structural.adjust_protonation(mol, N, 1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens + 1
    mol.to_rdkit()

    # remove one proton from the nitrogen
    bam.structural.adjust_protonation(mol, N, 0)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens
    mol.to_rdkit()

    # remove another proton from the nitrogen
    bam.structural.adjust_protonation(mol, N, -1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens - 1
    mol.to_rdkit()


def test_adjust_protonation_imine():
    mol = bam.read_smiles("CC=NC")
    mol.autolabel()

    n_hydrogens = len(mol.get_atoms("H", by="element"))

    N = mol.get_atom("N", by="element")

    # add one proton to the nitrogen
    bam.structural.adjust_protonation(mol, N, 1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens + 1
    mol.to_rdkit()

    # remove one proton from the nitrogen
    bam.structural.adjust_protonation(mol, N, 0)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens
    mol.to_rdkit()

    # remove another proton from the nitrogen
    try:
        bam.structural.adjust_protonation(mol, N, -1)
    except ValueError:
        pass
    else:
        raise RuntimeError("Should not be able to remove a proton from an imine")
