"""
Tests to check the behaviour of the bam.Molecule object
"""

import os
from copy import deepcopy
import numpy as np
import buildamol as bam
import Bio.PDB as bio

bam.load_sugars()
bam.load_amino_acids()

import tests.base as base

# set ALLOW_VISUAL to True to allow visual tests


def test_molecule_basic():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 0

    _mol = bam.utils.defaults.__bioPDBParser__.get_structure("MAN", base.MANPDB)
    mol = bam.Molecule(_mol)
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 0

    mol = bam.Molecule.from_smiles("OCC1OC(O)C(C(C1O)O)O", add_hydrogens=False)
    assert mol is not None

    assert len(mol.atoms) == 12
    assert len(mol.bonds) == 12

    mol = bam.Molecule.from_smiles("OCC1OC(O)C(C(C1O)O)O")
    assert mol is not None

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 24

    mol = bam.Molecule(_mol)

    a = mol.get_atom(1)
    assert a is not None

    b = mol.get_atom("C1")
    assert b is not None

    a.id = "HIHIHI"
    assert b.id == "HIHIHI"
    assert a is b

    a = mol.get_residue(1)
    assert a is not None

    b = mol.get_residue("MAN", by="name")
    assert b is not None

    assert a is b


def test_can_read_xml():
    mol = bam.Molecule.from_xml(base.GLCXML)
    assert mol.count_atoms() > 0
    assert mol.count_bonds() > 0
    if base.ALLOW_VISUAL:
        mol.show()


def test_can_write_pdb():
    import os

    glc = bam.Molecule.from_compound("GLC")
    assert glc is not None

    glc.to_pdb("GLC.pdb")

    glc2 = bam.Molecule.from_pdb("GLC.pdb")
    assert glc2 is not None

    assert len(glc.atoms) == len(glc2.atoms)
    assert len(glc.bonds) == len(glc2.bonds)

    assert glc.to_biopython() == glc2.to_biopython()

    os.remove("GLC.pdb")


def test_can_write_cif():
    import os

    glc = bam.Molecule.from_compound("GLC")
    assert glc is not None

    try:
        glc.to_cif("test.cif")
        with open("test.cif", "r") as f:
            lines = f.readlines()
        assert len(lines) != 0
    except Exception as e:
        raise e

    os.remove("test.cif")


def test_molecule_from_compound():
    glc = bam.Molecule.from_compound("GLC")
    assert glc is not None
    assert len(glc.atoms) == 24
    assert len(glc.bonds) == 24

    a = glc.atoms[0]
    assert a.full_id[0] == "GLC"
    assert a.full_id[1] == 0
    assert a.full_id[2] == "A"
    assert a.full_id[3] == ("H_GLC", 1, " ")
    assert a.full_id[4] == ("H1", " ")


def test_atomgraph_sync():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol2 = bam.Molecule.from_pdb(base.MANPDB)

    mol.apply_standard_bonds()
    mol2.apply_standard_bonds()

    assert mol is not mol2

    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 24

    assert len(mol.bonds) == len(mol2.bonds)

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph"

    for atom in mol2.get_atoms():
        assert atom in mol2._AtomGraph.nodes, f"Atom {atom} not in graph"

    # make the molecule larger
    mol % "14bb"
    mol += mol2

    assert len(mol.atoms) == 45

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph"

    mol.reindex(2, 100, 5)
    assert len(mol.atoms) == 45

    assert mol.chains[0].id == "B"
    assert mol.residues[0].id[1] == 100

    for atom in mol.get_atoms():
        assert atom in mol._AtomGraph.nodes, f"Atom {atom} not in graph after reindex"


def test_molecule_bonds():
    mol = bam.Molecule.from_pdb(base.MANPDB)

    assert len(mol.bonds) == 0

    mol.apply_standard_bonds()
    assert len(mol.bonds) == 24

    mol.bonds = []
    assert len(mol.bonds) == 0

    mol.infer_bonds()
    assert len(mol.bonds) == 24


def test_molecule_get_bonds():
    glc = bam.Molecule.from_compound("GLC")
    glc.repeat(2, "14bb")

    assert len(glc.bonds) == 46
    v = glc.draw()

    b = glc.get_bonds("O4")
    v.draw_edges(*b, color="red")
    assert len(b) == 4

    c1 = glc.get_atom("C1", residue=2)
    b = glc.get_bonds("O4", c1)
    v.draw_edges(*b, color="green")
    assert len(b) == 1

    if base.ALLOW_VISUAL:
        v.show()


def test_angles():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    # top = bam.resources.get_default_topology()
    # abstract = top.get_residue("MAN")

    for triplet, angle in mol.compute_angles().items():
        # triplet = [i.id for i in triplet]
        assert 90 < angle < 120, f"Angle {angle} is not in the expected range"

        # ics = abstract.get_internal_coordinates(*triplet, None, mode="partial")
        # _angle = "bond_angle_123"
        # if len(ics) == 0:
        #     ics = abstract.get_internal_coordinates(None, *triplet, mode="partial")
        #     _angle = "bond_angle_234"
        #     if len(ics) == 0:
        #         continue
        # _angle = getattr(ics[0], _angle)
        # assert (
        #     np.abs(_angle - angle) < 0.01
        # ), f"Angle {angle} does not match reference {_angle}"


def test_dihedrals():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    # top = bam.resources.get_default_topology()
    # abstract = top.get_residue("MAN")

    for quartet, dihedral in mol.compute_dihedrals().items():
        assert (
            -180 < dihedral < 180
        ), f"Dihedral {dihedral} is not in the expected range"
    #     quartet = [i.id for i in quartet]
    #     ics = abstract.get_internal_coordinates(*quartet)
    #     if len(ics) == 0:
    #         ics = abstract.get_internal_coordinates(*quartet[::-1])
    #         if len(ics) == 0:
    #             continue
    #     _dihedral = ics[0].dihedral

    #     assert (
    #         np.abs(_dihedral - dihedral) < 0.01


def test_add_atoms():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    pre = len(mol.atoms)

    new = bam.Atom("C99", np.array((0.5, 1.23, -0.5)))
    mol.add_atoms(new)

    assert len(mol.atoms) == pre + 1

    _new = mol.get_atom("C99")
    assert _new is not None
    assert _new.serial_number != 1

    mol.remove_atoms(_new)
    assert len(mol.atoms) == pre

    assert "C99" not in [i.id for i in mol.atoms]


def test_remove_atoms():
    glc = bam.Molecule.from_compound("GLC")
    glc.remove_atoms("C1", "O4")
    assert len(glc.atoms) == 22
    assert len(glc.bonds) == 18


def test_remove_residues():
    glc = bam.Molecule.from_compound("GLC")
    glc = glc % "14bb" * 6
    b = glc.get_bonds(glc.get_residue(3))
    b_before = glc.count_bonds()
    before = glc.count_atoms()
    r = glc.remove_residues(3)
    assert glc.count_residues() == 5
    assert glc.count_atoms() == before - len(r[0].child_list)
    assert glc.count_bonds() == b_before - len(b) - 2


def test_add_residues():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    other = bam.Molecule.from_compound("GLC")

    residues_pre = len(mol.residues)
    atoms_pre = len(mol.atoms)

    new = bam.Residue("NEW", " ", 2)
    mol.add_residues(new)

    assert len(mol.residues) == residues_pre + 1

    _new = mol.get_residue("NEW")
    assert _new is not None
    assert _new.id[1] != 1

    mol.remove_residues(_new)
    assert len(mol.residues) == residues_pre

    assert "NEW" not in [i.resname for i in mol.residues]

    atoms_pre = len(mol.atoms)
    mol.add_residues(*other.residues, _copy=True)
    assert len(other.residues) != 0
    assert len(mol.residues) == residues_pre + len(other.residues)
    assert len(mol.atoms) == atoms_pre + len(other.atoms)

    _seen_serials = set()
    for atom in mol.atoms:
        assert atom.serial_number not in _seen_serials
        _seen_serials.add(atom.serial_number)

    mol.remove_residues(2)
    assert len(mol.residues) == residues_pre
    assert len(mol.atoms) == atoms_pre

    new = bam.Molecule.from_compound("GLC")
    new = new.residues[0]
    mol.add_residues(new)
    assert mol.residues[-1].resname == "GLC"


def test_get_descendants():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    descendants = mol.get_descendants("C5", "C6")
    _received = len(descendants)
    _expected = 4
    assert (
        _received == _expected
    ), f"Expected {_expected} descendants, received {_received}"

    descendants = mol.get_descendants("O3", "C3")
    _received = len(descendants)
    _expected = 21
    assert (
        _received == _expected
    ), f"Expected {_expected} descendants, received {_received}"


def test_rotate_all():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = np.radians(45)

    current_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )
    current_refs = np.array((first.coord, second.coord))
    mol.rotate_around_bond(first, second, angle)

    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )
    new_refs = np.array((first.coord, second.coord))

    assert not np.allclose(current_coords, new_coords)
    assert np.allclose(current_refs, new_refs)

    # and rotate back to revert...
    mol.rotate_around_bond(first, second, -angle)
    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second]
    )

    assert np.allclose(current_coords, new_coords)


def test_rotate_some():
    mol = bam.Molecule.from_compound("GLC")
    mol.apply_standard_bonds()

    v = mol.draw()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = 10  # 10 degrees

    v.draw_point("first (O3)", first.coord, color="magenta")
    v.draw_point("second (C3)", second.coord, color="magenta")

    # HO3 should not change in dicrection O3---C3
    anchor = mol.get_atom("HO3")

    for i in range(5):
        current_coords = np.array(
            [
                i.coord
                for i in mol.get_atoms()
                if i != first and i != second and i != anchor
            ]
        )
        current_refs = np.array((first.coord, second.coord))
        current_anchor = anchor.coord

        mol.rotate_around_bond(first, second, angle, descendants_only=True)

        new_coords = np.array(
            [
                i.coord
                for i in mol.get_atoms()
                if i != first and i != second and i != anchor
            ]
        )
        new_refs = np.array((first.coord, second.coord))
        new_anchor = anchor.coord

        assert not np.allclose(current_coords, new_coords)
        assert np.allclose(current_refs, new_refs)
        assert np.allclose(current_anchor, new_anchor)

        v.draw_edges(*mol.bonds)
        for atom in mol.atoms:
            v.draw_point(atom.id, atom.coord, opacity=0.4, showlegend=False)

    if base.ALLOW_VISUAL:
        v.show()


def test_rotate_some_inverse():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.apply_standard_bonds()

    first = mol.get_atom("O3")
    second = mol.get_atom("C3")
    angle = np.radians(45)

    # HO3 should be the only change in dicrection C3---O3
    anchor = mol.get_atom("HO3")

    current_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    current_refs = np.array((first.coord, second.coord))
    current_anchor = anchor.coord

    mol.rotate_around_bond(second, first, angle, descendants_only=True)

    new_coords = np.array(
        [i.coord for i in mol.get_atoms() if i != first and i != second and i != anchor]
    )
    new_refs = np.array((first.coord, second.coord))
    new_anchor = anchor.coord

    assert np.allclose(current_coords, new_coords)
    assert np.allclose(current_refs, new_refs)
    assert not np.allclose(current_anchor, new_anchor)


def test_adjust_indexing():
    mol = bam.Molecule.from_compound("MAN")
    other = mol.copy()

    mol.adjust_indexing(other)

    assert len(mol.atoms) == len(other.atoms)
    assert len(mol.residues) == len(other.residues)

    assert mol.residues[0].id[1] == 1
    assert other.residues[0].id[1] == 2

    assert mol.get_atom(1) is not None
    assert (
        other.get_atom(1) is None and other.get_atom(mol.count_atoms() + 1) is not None
    )


def test_adjust_indexing_with_add_residues():
    mol = bam.Molecule.from_compound("MAN")
    other = mol.copy()

    mol.adjust_indexing(other)

    before = len(mol.atoms)
    mol.add_residues(*other.residues)

    assert len(mol.atoms) == before * 2

    _seen_serials = set()
    for atom in mol.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    assert mol.residues[0].id[1] == 1
    assert len(mol.residues) == 2


def test_set_linkage():
    mol = bam.Molecule.from_compound("GLC")
    mol.set_linkage("14bb")

    assert mol._linkage is not None
    assert not isinstance(mol._linkage, str)

    mol.set_linkage()
    assert mol._linkage is None

    # set the patch with fancy dunder methods
    mol % "14bb"

    assert mol._linkage is not None
    assert not isinstance(mol._linkage, str)

    mol % None
    assert mol._linkage is None

    mol %= "14bb"

    assert mol is not None
    assert mol._linkage is not None
    assert not isinstance(mol._linkage, str)


def test_attach_with_patch():
    glc = bam.Molecule.from_compound("GLC")
    glc.set_linkage("14bb")

    glc2 = deepcopy(glc)

    _current_residues = len(glc.residues)
    glc.attach(glc2)

    assert len(glc.residues) == _current_residues * 2

    # glc2 should not have been affected since it was a copy
    assert len(glc2.residues) == _current_residues

    # attach with fancy dunder methods
    glc = deepcopy(glc2)

    new = glc + glc2
    assert new is not glc
    assert len(new.residues) == _current_residues * 2
    assert len(glc.residues) == _current_residues
    assert len(glc2.residues) == _current_residues

    glc += glc2
    assert len(glc.residues) == _current_residues * 2
    assert len(glc2.residues) == _current_residues


def test_infer_bonds():
    mol = bam.molecule("GLC")
    assert mol.count_bonds() == 24
    mol.bonds = None
    assert mol.count_bonds() == 0

    mol.infer_bonds()
    assert mol.count_bonds() == 24

    b = mol.bonds[0]
    b = mol._AtomGraph[b.atom1][b.atom2]
    assert b.get("bond_order") == 1
    assert b.get("bond_obj", None) is not None


def test_infer_residue_connections():
    glc = bam.Molecule.from_compound("GLC")
    glc.repeat(10, "14bb")
    glc.bonds = []

    connections = glc.infer_residue_connections(triplet=False)
    assert len(connections) == 9


def test_find_clashes():
    mol = bam.molecule("GLC")
    mol = mol.repeat(10, "14bb")
    clashes = mol.find_clashes()
    assert len(clashes) == 0

    edges = mol.get_residue_connections()
    for edge in edges:
        mol.rotate_around_bond(
            *edge, np.random.uniform(-180, 180), descendants_only=True
        )

    clashes = mol.find_clashes()
    assert len(clashes) != 0


def test_attach_with_recipe():
    recipe = bam.Linkage()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    glc = bam.Molecule.from_compound("GLC")
    glc.set_linkage(recipe)

    glc2 = deepcopy(glc)

    _current_residues = len(glc.residues)
    glc.attach(glc2)

    assert len(glc.residues) == _current_residues * 2

    # glc2 should not have been affected since it was a copy
    assert len(glc2.residues) == _current_residues

    # attach with fancy dunder methods
    glc = deepcopy(glc2)

    new = glc + glc2
    assert new is not glc
    assert len(new.residues) == _current_residues * 2
    assert len(glc.residues) == _current_residues
    assert len(glc2.residues) == _current_residues

    glc += glc2
    assert len(glc.residues) == _current_residues * 2
    assert len(glc2.residues) == _current_residues


def test_multiply_with_patch():
    man = bam.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)

    # set 14bb patch
    n = 10
    man % "14bb"
    man = man * n

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.compute_angles().values():
        assert 100 < angle < 130

    if base.ALLOW_VISUAL:
        man.show()


def test_multiply_with_recipe():
    man = bam.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man.locked_bonds)

    # make recipe for 14bb
    recipe = bam.Linkage()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    # set 14bb recipe
    n = 10
    man % recipe
    man = man * n

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man.locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert pre_locked < new_locked

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.compute_angles().values():
        assert 100 < angle < 130

    if base.ALLOW_VISUAL:
        man.show()


def test_repeat():
    man = bam.Molecule.from_compound("GLC")
    man.lock_all()

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man.locked_bonds)

    n = 10
    man = man.repeat(n, "14bb")

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.compute_angles().values():
        assert 100 < angle < 130

    if base.ALLOW_VISUAL:
        man.show()


def test_repeat_with_recipe():
    man = bam.Molecule.from_compound("GLC")
    man.lock_all()

    # make recipe for 14bb
    recipe = bam.Linkage()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    pre_residues = len(man.residues)
    pre_atoms = len(man.atoms)
    pre_bonds = len(man.bonds)
    pre_locked = len(man.locked_bonds)

    n = 10
    man = man.repeat(n, recipe)

    new_residues = len(man.residues)
    new_atoms = len(man.atoms)
    new_bonds = len(man.bonds)
    new_locked = len(man.locked_bonds)

    assert new_residues == pre_residues * n
    assert n * 0.75 * pre_atoms < new_atoms < pre_atoms * n
    assert n * 0.75 * pre_bonds < new_bonds < pre_bonds * n
    assert pre_locked < new_locked

    # test that the new molecule has no weird bond lengths
    for atom1, atom2 in man.bonds:
        dist = np.linalg.norm(atom1.coord - atom2.coord)
        assert 0.95 < dist < 1.8

    # test that the new molecule has no weird bond angles
    for angle in man.compute_angles().values():
        assert 100 < angle < 130

    if base.ALLOW_VISUAL:
        man.show()


def test_make_mannose8():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN

    """
    bma = bam.Molecule.from_compound("BMA")
    nag = bam.Molecule.from_compound("NAG")
    man = bam.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    nag % "14bb"
    nag *= 2
    man8 = nag + bma

    # now we attach the 13ab MAN to the BMA
    man8 % "13ab"
    man8 += man

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man % "16ab"
    man_branch = man * 2

    man_branch % "13ab"
    man_branch @ 1  # attach to residue 1 (not the default last one)
    man_branch += man

    man_branch @ None  # reset the attachment point
    man_branch % "12aa"
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    man8 @ -2
    man_branch @ 1

    man8 % "16ab"
    man8 += man_branch

    v = man8.draw()
    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        if not length:
            v.draw_edges([bond], color="magenta", linewidth=3)

    for triplet, angle in man8.compute_angles().items():
        if not 100 < angle < 150:
            v.draw_vector(
                None,
                triplet[0].coord,
                triplet[1].coord,
                color="red",
            )
            v.draw_vector(
                None,
                triplet[2].coord,
                triplet[1].coord,
                color="red",
            )

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    v.draw_edges(*man8.get_residue_connections(triplet=True), color="red", linewidth=3)

    if base.ALLOW_VISUAL:
        v.show()

    g = man8.make_residue_graph(detailed=False)
    g2 = man8.make_residue_graph(detailed=True)

    assert g is not g2
    assert len(g.nodes) < len(g2.nodes)

    v = g.draw()
    if base.ALLOW_VISUAL:
        v.show()

    v = g2.draw()
    if base.ALLOW_VISUAL:
        v.show()

    try:
        man8.to_pdb("man8.pdb")
    except Exception as e:
        raise e
    else:
        os.remove("man8.pdb")


def test_make_mannose8_2():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN

    """

    bma = bam.Molecule.from_compound("BMA")
    nag = bam.Molecule.from_compound("NAG")
    man = bam.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    man8 = nag % "14bb" * 2 + bma

    # now we attach the 13ab MAN to the BMA
    man8 % "13ab"
    man8 += man

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man_branch = man % "16ab" * 2

    man_branch @ 1 % "13ab"  # attach to residue 1 (not the default last one)
    man_branch += man

    man_branch % "12aa" @ None  # reset the attachment point
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    man8 @ -2 % "16ab"
    man_branch @ 1

    _man8 = man8.copy()
    man8 += man_branch

    # just checkin if the one line syntax is the same as the two line syntax
    _man8 = _man8 @ -2 % "16ab" + man_branch @ 1
    assert len(man8.atoms) == len(_man8.atoms)

    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        assert length, "Bond length is not in range 0.95 - 1.8"

    for angle in man8.compute_angles().values():
        assert 100 < angle < 130

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    if base.ALLOW_VISUAL:
        v = man8.make_residue_graph(detailed=False).draw()
        v.show()


def test_make_mannose8_3():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN
    ```
    """

    bma = bam.Molecule.from_compound("BMA")
    nag = bam.Molecule.from_compound("NAG")
    man = bam.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    nag.set_linkage("14bb")
    nag2 = nag + nag

    man8 = nag2.attach(bma)

    # now we attach the 13ab MAN to the BMA
    man8.attach(man, "13ab")

    # now we make the mannose branch
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man.set_linkage("16ab")
    man_branch = man * 2

    man_branch.set_attach_residue(1)
    man_branch.set_linkage("13ab")
    man_branch.attach(man)

    man_branch.set_linkage("12aa")
    man_branch.set_attach_residue()
    man_branch += man

    # and now we attach the man branch to the NAG-NAG--BMA---MAN
    # but at the second last residue (BMA), not the last one
    # man8.set_attach_residue(-2)
    # man_branch.set_attach_residue(1)
    # man8.set_linkage("16ab")
    man8.attach(man_branch, "16ab", at_residue=-2, other_residue=1)

    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        assert length, "Bond length is not in range 0.95 - 1.8"

    for angle in man8.compute_angles().values():
        assert 100 < angle < 130

    _seen_serials = set()
    for atom in man8.atoms:
        assert atom.get_serial_number() not in _seen_serials
        _seen_serials.add(atom.get_serial_number())

    if base.ALLOW_VISUAL:

        v = man8.draw()
        colors = [
            "red",
            "green",
            "blue",
            "magenta",
            "cyan",
            "orange",
            "purple",
            "pink",
            "brown",
            "grey",
            "black",
        ]
        idx = 0
        for residue in man8.residues:
            for bond in man8.bonds:
                if bond[0].get_parent() == residue and bond[1].get_parent() == residue:
                    v.draw_edges(bond, color=colors[idx], linewidth=3)
            idx += 1

        v += man8.optimize(inplace=True).move([20, 0, 0]).draw(atoms=False)
        v.show()


def test_make_mannose8_with_recipe():
    """
    Structure to build:

    ```
                               MAN
                                |
                              (16ab)
                                |
    ~ --- NAG                  MAN -(13ab)- MAN -(12aa)- MAN
            \\                  /
            (14bb)          (16ab)
              \\              /
              NAG -(14bb)- BMA
                            \\
                            (13ab)
                               \\
                               MAN

    """
    bam.load_sugars()

    bma = bam.Molecule.from_compound("BMA")
    nag = bam.Molecule.from_compound("NAG")
    man = bam.Molecule.from_compound("MAN")

    # make the NAG-NAG--BMA (we can always use the 14bb patch)
    nag % "14bb"
    nag.repeat(2)
    man8 = nag + bma

    man8 % "13ab"
    man8 += man

    # now we start adding the mannoses
    # MAN --- MAN
    #  \
    #  MAN --- MAN
    man % "16ab"
    man_branch = man * 2
    man_branch @ 1 % "13ab"
    man_branch += man

    man_branch.attach(man, "12aa", at_residue=-1, use_patch=False)

    man8 @ -2
    man_branch @ 1

    recipe_16ab = bam.Linkage("16ab")
    recipe_16ab.add_bond(("O6", "C1"))
    recipe_16ab.add_delete("HO6", "target")
    recipe_16ab.add_delete("HO1", "source")
    recipe_16ab.add_delete("O1", "source")

    man8 % recipe_16ab
    man8 += man_branch

    for bond in man8.bonds:
        assert bond[0] in man8.atoms
        assert bond[1] in man8.atoms
        length = 0.95 < np.linalg.norm(bond[1].coord - bond[0].coord) < 1.8
        assert length, "Bond length is not in range 0.95 - 1.8"

    for angle in man8.compute_angles().values():
        assert 100 < angle < 130 if angle != 0 else True

    all_serials = [atom.serial_number for atom in man8.get_atoms()]
    assert len(set(all_serials)) == len(all_serials)

    if base.ALLOW_VISUAL:
        v = man8.draw()
        colors = [
            "red",
            "green",
            "blue",
            "magenta",
            "cyan",
            "orange",
            "purple",
            "pink",
            "brown",
            "grey",
            "black",
        ]
        idx = 0
        for residue in man8.residues:
            for bond in man8.bonds:
                if bond[0].get_parent() == residue and bond[1].get_parent() == residue:
                    v.draw_edges(bond, color=colors[idx], linewidth=3)
            idx += 1

        v.show()


def test_relabel():
    scrambled = bam.Molecule.from_compound("BMA")

    # relabel to elementwise order without specific connectivity
    counts = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}

    for atom in scrambled.atoms:
        counts[atom.element] += 1
        atom.id = atom.element + str(counts[atom.element])

    # randomly rotate the molecule
    random_vector = np.random.rand(3) * 10

    # choose some rotation axis
    axis = np.random.randint(0, len(scrambled.bonds) - 1)
    axis = scrambled.bonds[axis]

    for atom in scrambled.atoms:
        atom.coord += random_vector

    scrambled.rotate_around_bond(*axis, np.random.randint(0, 180))

    old_scrambled = set(i.id for i in scrambled.atoms)
    old_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    # now we relabel the atoms
    scrambled.autolabel()

    new_scrambled = set(i.id for i in scrambled.atoms)
    new_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    # now we check if the relabeling worked

    assert scrambled.get_bonds("C1", "C2") != []
    assert scrambled.get_bonds("C6", "H61") != []
    assert old_scrambled.difference(new_scrambled) != set()
    assert np.allclose(old_scrambled_coords, new_scrambled_coords)


def test_rotate_descendants_2():
    glc = bam.Molecule.from_compound("GLC")
    glc.lock_all()
    glc.repeat(4, "14bb")
    connections = glc.get_residue_connections()

    v = glc.draw()

    v.draw_point(
        "root",
        glc.get_atom(1).coord,
        color="red",
        showlegend=True,
    )

    cdx = 1
    for c in connections:
        glc.unlock_bond(*c)

        v.draw_vector(
            f"""[{cdx}]    {c[0].full_id[3:]} -> {c[1].full_id[3:]}""",
            c[0].coord,
            c[1].coord,
            color="cyan",
            elongate=1.2,
            showlegend=True,
        )
        cdx += 1

    v.draw_edges(
        *glc.get_bonds(glc.residues[0]),
        color="limegreen",
        linewidth=5,
        opacity=1,
    )

    opacities = np.linspace(0.2, 0.8, len(connections))
    cdx = 0
    for c in connections:
        assert c not in glc.locked_bonds

        descendants = glc.get_descendants(*c)
        ancestors = glc.get_ancestors(*c)

        old_coords_descendants = np.array([i.coord for i in descendants])
        old_coords_ancestors = np.array([i.coord for i in ancestors])

        glc.rotate_around_bond(*c, 5, descendants_only=True)

        new_coords_descendants = np.array([i.coord for i in descendants])
        new_coords_ancestors = np.array([i.coord for i in ancestors])

        assert not np.allclose(old_coords_descendants, new_coords_descendants)
        assert np.allclose(old_coords_ancestors, new_coords_ancestors)

        v.draw_edges(*glc.bonds, color="lightgreen", opacity=opacities[cdx])

        v.draw_edges(
            *glc.get_bonds(glc.residues[0]),
            color="green",
            linewidth=3,
            opacity=1,  # opacities[cdx],
        )

        cdx += 1
    # for atom in glc.atoms:
    #     v.draw_point(atom.id, atom.coord, color="red", opacity=0.3, showlegend=False)

    # v.draw_edges(glc.bonds, color="red")

    v.draw_edges(*glc.get_bonds(glc.residues[0]), color="teal", linewidth=5, opacity=1)
    if base.ALLOW_VISUAL:
        v.show()


def test_from_rdkit():
    from rdkit import Chem

    rdkit_mol = Chem.MolFromPDBFile(str(base.GLCPDB), removeHs=False)
    assert sum(1 for i in rdkit_mol.GetAtoms()) == 24

    mol = bam.Molecule.from_rdkit(rdkit_mol, "myman")
    assert len(mol.atoms) == 24
    assert len(mol.bonds) == 24
    assert len(mol.residues) == 1
    assert len(mol.chains) == 1
    assert mol.atoms[0].coord.sum() != 0  # check if coords are set
    assert mol.atoms[5].coord.sum() != 0
    assert mol.id == "myman"


def test_to_rdkit():
    glc = bam.Molecule.from_compound("GLC")
    rdkit_mol = glc.to_rdkit()
    assert sum(1 for i in rdkit_mol.GetAtoms()) == 24
    assert sum(1 for i in rdkit_mol.GetBonds()) == 24


def test_to_and_from_xml():
    glc = bam.Molecule.from_compound("GLC")

    glc.to_xml("glc.xml")
    glc2 = bam.Molecule.from_xml("glc.xml")
    assert (glc2.get_coords()).sum() != 0
    assert glc.count_atoms() == glc2.count_atoms()
    assert glc.count_bonds() == glc2.count_bonds()
    assert glc.to_biopython() == glc2.to_biopython()

    glc.add_model(0)
    glc.set_model(1)
    glc.move([10, 0, 0])
    glc.add_model(1)
    glc.set_model(2)
    glc.move([0, 10, 0])
    glc.to_xml("glc.xml")
    glc2 = bam.Molecule.from_xml("glc.xml")

    assert glc.count_models() == glc2.count_models()
    assert (glc2.get_coords()).sum() != 0
    assert glc.count_atoms() == glc2.count_atoms()
    assert glc.count_bonds() == glc2.count_bonds()
    assert glc.to_biopython() == glc2.to_biopython()

    os.remove("glc.xml")
    if base.ALLOW_VISUAL:

        v = glc2.draw()
        for model in glc2.get_models():
            glc2.set_model(model)
            v.draw_edges(*glc2.bonds, color="red", opacity=0.5)
        v.show()


def test_work_with_pubchem():
    phprop = bam.Molecule.from_pubchem("2-[4-(2-methylpropyl)phenyl]propanal")
    phprop.autolabel()
    phprop.residues[0].resname = "PRO"

    l = bam.linkage(
        "C8",
        "C8",
        ["H8"],
        ["H8"],
    )
    phprop % l
    phprop = phprop + phprop
    if base.ALLOW_VISUAL:
        phprop.show()


def test_chlorine():
    benz = bam.Molecule.from_pubchem("1,2,4-trichloro-5-methylbenzene")
    benz.autolabel()
    benz.residues[0].resname = "CBZ"
    if base.ALLOW_VISUAL:
        benz.show()


# def test_infer_missing():
#     man = bam.molecule("MAN")
#     to_delete = ("HO4", "O4", "HO3")
#     man.remove_atoms(*to_delete)

#     v = man.draw()

#     for i in to_delete:
#         assert man.get_atom(i) is None

#     assert bam.has_residue("MAN")

#     man.infer_missing_atoms()

#     for i in to_delete:
#         assert man.get_atom(i) is not None

#     v.draw_edges(man.bonds, color="red")
#     v.show()


def test_peptide_link():
    bam.load_amino_acids()
    his = bam.Molecule.from_compound("HIS")
    ser = bam.Molecule.from_compound("SER")

    peptide_link = bam.linkage("C", "N", ["OXT", "HXT"], ["H"])
    his % peptide_link
    peptide = his.attach(ser, inplace=False)
    peptide.attach(his)
    assert peptide.count_residues() == 3
    if base.ALLOW_VISUAL:
        peptide.show()


def test_peptide_link_multiple():
    bam.load_amino_acids()
    met = bam.Molecule.from_compound("MET")
    his = bam.Molecule.from_compound("HIS")
    ser = bam.Molecule.from_compound("SER")

    peptide_link = bam.linkage("C", "N", ["OXT", "HXT"], ["H"])
    peptide = bam.connect(met, his, peptide_link)
    peptide = bam.connect(peptide, ser, peptide_link)
    peptide = bam.connect(peptide, his, peptide_link)
    peptide = bam.connect(peptide, ser, peptide_link)
    assert peptide.count_residues() == 5

    if base.ALLOW_VISUAL:
        peptide.show()


def test_ferrocynyl():
    bam.load_small_molecules()
    core = bam.read_smiles("[H]P1([H])=NP([H])([H])=NP([H])([H])=N1")
    core.rename_residue(1, "CRE")

    phenol = bam.molecule("phenol")
    phenol.autolabel().rename_residue(1, "PHE")

    linker = bam.read_smiles("C=NNC")
    linker.autolabel().rename_residue(1, "LNK")

    phos = bam.read_smiles("[P](=S)([Cl])([Cl])[Cl]")
    phos.rename_residue(1, "PHO")

    link2 = bam.linkage("N2", "P1", None, ["CL1"])
    link3 = bam.linkage("C4", "C1")

    per1 = phenol % link3 + (linker % link2 + phos) @ 1

    link4 = bam.linkage("O1", "P1", None, ["CL1"])

    per2 = per1 % link4 + phos

    link5 = bam.linkage("P1", "O1", ["CL2"], None)

    per2 = per2 % link5 + per1

    link6 = bam.linkage("P1", "N2", ["CL3"])
    per3 = per2.attach(linker, link6, at_residue=4, inplace=False)
    if base.ALLOW_VISUAL:
        per3.show()


def test_keeps_pdb_enumeration():
    glc3 = bam.molecule("GLC") % "14bb" * 3

    glc3.reindex(start_chainid=3, start_atomid=8)

    glc3.get_residue(1).serial_number = 99
    glc3.get_residue(2).serial_number = 120
    glc3.get_residue(3).serial_number = 150

    glc3.atoms[10].serial_number = 999

    glc3.to_pdb("test.pdb")

    _new = bam.Molecule.from_pdb("test.pdb")

    assert _new.residues[0].serial_number == 99
    assert _new.residues[1].serial_number == 120
    assert _new.residues[2].serial_number == 150

    assert min(i.serial_number for i in _new.atoms) == 8
    assert max(i.serial_number for i in _new.atoms) == 999

    os.remove("test.pdb")


def test_molecule_move():
    mol = bam.Molecule.from_compound("GLC")
    old_coords = np.array([i.coord for i in mol.atoms])

    mol.move((1, 0, 0))

    new_coords = np.array([i.coord for i in mol.atoms])

    assert new_coords.sum() != old_coords.sum()
    assert np.all(new_coords[:, 0] == old_coords[:, 0] + 1)
    assert np.all(new_coords[:, 1] == old_coords[:, 1])
    assert np.all(new_coords[:, 2] == old_coords[:, 2])


def test_molecule_rotate():
    mol = bam.Molecule.from_compound("GLC")
    d = mol.draw()
    mol.rotate(90, (1, 0, 0))
    d.draw_edges(*mol.bonds, color="red")

    mol.rotate(90, (0, 1, 0), center=mol.get_atom("HO6").coord)
    d.draw_edges(*mol.bonds, color="blue")

    if base.ALLOW_VISUAL:
        d.show()


def test_molecule_transpose():
    mol = bam.Molecule.from_compound("GLC")

    d = mol.draw()

    mol.transpose((0, 0, 0), 90, (1, 0, 0))
    d.draw_edges(*mol.get_bonds(), color="red", showlegend=False)

    mol.transpose((1, 2, 0), 90, (0, 1, 0))

    d.draw_edges(*mol.get_bonds(), color="blue", showlegend=False)

    d.draw_atom(mol.get_atom("HO6"), color="pink")

    mol.transpose((0, 0, 0), 90, (0, 0, 1), center=mol.get_atom("HO6").coord)

    d.draw_edges(*mol.get_bonds(), color="green", showlegend=False)

    if base.ALLOW_VISUAL:
        d.show()


def test_base_classes_move():
    mol = bam.Molecule.from_compound("GLC")

    old_coords = np.array([i.coord for i in mol.atoms])

    d = mol.draw()

    atom = mol.atoms[0]
    residue = mol.residues[0]
    chain = mol.chains[0]
    model = mol.model
    structure = mol.structure

    residue.move((1, 0, 0))
    chain.move((1, 0, 0))
    model.move((1, 0, 0))
    structure.move((1, 0, 0))
    atom.move((1, 0, 0))

    new_coords = np.array([i.coord for i in mol.atoms])
    assert new_coords.sum() != old_coords.sum()

    d.draw_edges(*mol.bonds, color="red")

    residue.rotate(10, (1, 0, 0), axis_is_absolute=False)
    chain.rotate(10, (1, 0, 0), axis_is_absolute=False)
    model.rotate(10, (1, 0, 0), axis_is_absolute=False)
    structure.rotate(10, (1, 0, 0), axis_is_absolute=False)
    atom.rotate(10, (1, 0, 0), axis_is_absolute=False)

    d.draw_edges(*mol.bonds, color="blue")

    residue.rotate(10, (1, 0, 0), axis_is_absolute=True)
    chain.rotate(10, (1, 0, 0), axis_is_absolute=True)
    model.rotate(10, (1, 0, 0), axis_is_absolute=True)
    structure.rotate(10, (1, 0, 0), axis_is_absolute=True)
    atom.rotate(10, (1, 0, 0), axis_is_absolute=True)

    d.draw_edges(*mol.bonds, color="orange")

    # the final output will have one bond that is really
    # bad because we move and rotate one individual atom
    # to an impossible position, so that's fine...
    if base.ALLOW_VISUAL:
        d.show()


def test_linkage_apply_deletes():
    mol = bam.Molecule.from_compound("GLC")

    # d = mol.draw()

    link = bam.linkage("O3", "C1", ["HO3"], ["H1"])

    link.apply_deletes(mol)

    assert mol.get_atom("HO3") is None
    assert mol.get_atom("H1") is not None

    link.apply_deletes(None, mol)

    assert mol.get_atom("HO3") is None
    assert mol.get_atom("H1") is None

    # d.draw_edges(*mol.bonds, color="red", linewidth=3)

    # d.show()


def test_linkage_apply_bond():
    mol = bam.Molecule.from_compound("GLC")

    d = mol.draw()

    link = bam.linkage("O3", "C1", ["HO3"], ["H1"])

    link.apply_deletes(mol, mol)
    link.apply_bond(mol, mol)

    assert mol.get_bond("O3", "C1") is not None

    d.draw_edges(*mol.bonds, color="red", linewidth=3)
    if base.ALLOW_VISUAL:
        d.show()


def test_squash():
    mol = bam.Molecule.from_compound("GLC")

    mol = mol % "14bb" * 3

    assert len(mol.residues) == 3
    atoms = len(mol.atoms)
    bonds = len(mol.bonds)

    mol.squash()
    assert len(mol.residues) == 1
    assert len(mol.atoms) == atoms
    assert len(mol.bonds) == bonds

    # mol.show()


def test_merge():
    mol = bam.Molecule.from_compound("GLC")

    mol2 = mol.copy()

    mol2 = mol2 % "14bb" * 3

    mol2.move((10, 10, 0))

    mol.merge(mol2)

    assert len(mol.chains) == 2
    assert len(mol.residues) == 4

    # mol.show()


def test_atom_new():
    a = bam.Atom.new("C")
    assert a.element == "C"
    assert a.id == "C1"
    assert a.coord.sum() == 0
    assert a.serial_number == 1

    a = bam.Atom.new("C ")
    assert a.element == "C"
    assert a.id == "C"

    a = bam.Atom.new("C2", coord=(1, 2, 3))
    assert a.element == "C"
    assert a.id == "C2"
    assert a.coord.sum() == 6

    a = bam.Atom.new("CA")
    assert a.element == "CA"
    assert a.id == "CA1"

    a = bam.Atom.new(" CA")
    assert a.element == "C"
    assert a.id == "CA"

    a = bam.Atom.new("CA ")
    assert a.element == "CA"
    assert a.id == "CA"

    a = bam.Atom.new(" CA  ")
    assert a.element == "C"
    assert a.id == "CA"

    a = bam.Atom.new("CA_D")
    assert a.element == "CA"
    assert a.id == "CAD"

    a = bam.Atom.new("5C")
    assert a.element == "C"
    assert a.id == "5C"

    a = bam.Atom.from_element("Fe")
    assert a.element == "FE"
    assert a.id == "FE1"


def test_base_matches_and_equal():
    mol = bam.Molecule.from_compound("GLC")
    mol2 = mol.copy()

    for a1, a2 in zip(mol.atoms, mol2.atoms):
        assert a1.matches(a2)
        assert a2.matches(a1)
        assert a1.equals(a2)
        assert a2.equals(a1)

    for r1, r2 in zip(mol.residues, mol2.residues):
        assert r1.matches(r2)
        assert r2.matches(r1)
        assert r1.equals(r2)
        assert r2.equals(r1)

    for c1, c2 in zip(mol.chains, mol2.chains):
        assert c1.matches(c2)
        assert c2.matches(c1)
        assert c1.equals(c2)
        assert c2.equals(c1)

    for m1, m2 in zip(mol.structure.get_models(), mol2.structure.get_models()):
        assert m1.matches(m2)
        assert m2.matches(m1)
        assert m1.equals(m2)
        assert m2.equals(m1)

    assert mol.structure.matches(mol2.structure)
    assert mol2.structure.matches(mol.structure)
    assert mol.structure.equals(mol2.structure)
    assert mol2.structure.equals(mol.structure)


def test_base_matches_and_equal_2():
    mol = bam.Molecule.from_compound("GLC")
    mol2 = bam.Molecule.from_compound("GAL")
    mol_1 = (mol % "14bb" * 3) % "12ab" + mol2
    mol_2 = (mol % "14bb" * 3) % "12ab" + mol2
    mol = mol_1
    mol2 = mol_2

    for a1, a2 in zip(mol.atoms, mol2.atoms):
        assert a1.matches(a2)
        assert a2.matches(a1)
        assert a1.equals(a2)
        assert a2.equals(a1)
        assert not a1 == a2

    for r1, r2 in zip(mol.residues, mol2.residues):
        assert r1.matches(r2)
        assert r2.matches(r1)
        assert r1.equals(r2)
        assert r2.equals(r1)
        assert not r1 == r2

    for c1, c2 in zip(mol.chains, mol2.chains):
        assert c1.matches(c2)
        assert c2.matches(c1)
        assert c1.equals(c2)
        assert c2.equals(c1)
        assert not c1 == c2

    for m1, m2 in zip(mol.structure.get_models(), mol2.structure.get_models()):
        assert m1.matches(m2)
        assert m2.matches(m1)
        assert m1.equals(m2)
        assert m2.equals(m1)
        assert not m1 == m2

    assert mol.structure.matches(mol2.structure)
    assert mol2.structure.matches(mol.structure)
    assert mol.structure.equals(mol2.structure)
    assert mol2.structure.equals(mol.structure)
    assert not mol.structure == mol2.structure


def test_base_matches_with_reindex():
    mol = bam.Molecule.from_compound("GLC")
    mol2 = mol.copy()
    mol2.reindex(1, 2, len(mol2.atoms))

    for a1, a2 in zip(mol.atoms, mol2.atoms):
        assert a1.matches(a2)
        assert a2.matches(a1)
        assert not a1.equals(a2)
        assert not a2.equals(a1)

    for r1, r2 in zip(mol.residues, mol2.residues):
        assert r1.matches(r2)
        assert r2.matches(r1)
        assert r1.equals(r2)
        assert r2.equals(r1)

    for c1, c2 in zip(mol.chains, mol2.chains):
        assert c1.matches(c2)
        assert c2.matches(c1)
        assert c1.equals(c2)
        assert c2.equals(c1)

    for m1, m2 in zip(mol.structure.get_models(), mol2.structure.get_models()):
        assert m1.matches(m2)
        assert m2.matches(m1)
        assert m1.equals(m2)
        assert m2.equals(m1)

    assert mol.structure.matches(mol2.structure)
    assert mol2.structure.matches(mol.structure)
    assert mol.structure.equals(mol2.structure)
    assert mol2.structure.equals(mol.structure)


def test_base_no_matches():
    mol = bam.Molecule.from_compound("GLC")
    mol2 = bam.Molecule.from_compound("GAL")

    for a1, a2 in zip(mol.atoms, mol2.atoms):
        assert not a1.matches(a2)
        assert not a2.matches(a1)
        assert not a1.equals(a2)  # on atom level all should still be equal
        assert not a2.equals(a1)

    for r1, r2 in zip(mol.residues, mol2.residues):
        assert not r1.matches(r2)
        assert not r2.matches(r1)
        assert not r1.equals(r2)
        assert not r2.equals(r1)

    for c1, c2 in zip(mol.chains, mol2.chains):
        assert not c1.matches(c2)
        assert not c2.matches(c1)
        assert not c1.equals(c2)
        assert not c2.equals(c1)

    for m1, m2 in zip(mol.structure.get_models(), mol2.structure.get_models()):
        assert not m1.matches(m2)
        assert not m2.matches(m1)
        assert not m1.equals(m2)
        assert not m2.equals(m1)

    assert not mol.structure.matches(mol2.structure)
    assert not mol2.structure.matches(mol.structure)
    assert not mol.structure.equals(mol2.structure)
    assert not mol2.structure.equals(mol.structure)


def test_polyphenylene():
    bam.load_small_molecules()
    benzene = bam.get_compound("benzene")
    # don't worry about the skewed visuals
    # (that's because there are miniscule difference in the z-axis, < 1e-2A)
    # benzene.show()

    # set up the linkage instructions
    link = bam.linkage("C1", "C1")

    # start with the centrla benzene ring
    periphery = benzene.copy()

    for carbon in range(1, 6):

        # update the linkage to the next carbon
        link.atom1 = f"C{carbon}"

        # attach one benzene to the central one
        periphery.attach(benzene, link, at_residue=1)

    # periphery.show()

    # setup a new linkage
    link2 = bam.linkage("C1", "C4")

    # set a new attach residue for the periphery (now it's the second residue)
    periphery.set_attach_residue(2)

    # now make the central core
    core = benzene.copy()

    # and attach the periphery to the core
    for carbon in bam.utils.element_range("C", 6):
        link2.atom1 = carbon

        core.attach(periphery, link2, at_residue=1)

    if base.ALLOW_VISUAL:
        core.show()


def test_add_hydrogens_glucose():
    mol = bam.Molecule.from_compound("GLC")
    ref_n = mol.count_atoms()

    mol.remove_atoms(*mol.get_atoms("H", by="element"))
    assert mol.count_atoms() < ref_n

    mol.add_hydrogens()

    assert len(mol._atoms) == ref_n

    if base.ALLOW_VISUAL:
        mol.show()


def test_add_hydrogens_mannose9():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False, infer_bond_orders=True)
    mol.autolabel()

    ref_n = mol.count_atoms()
    old_hydrogens = mol.get_atoms("H", by="element")
    n_hydrogens = len(old_hydrogens)
    mol.remove_atoms(old_hydrogens)

    assert mol.count_atoms() == ref_n - n_hydrogens

    v = mol.draw()
    mol.add_hydrogens()

    new_hydrogens = mol.get_atoms("H", by="element")

    # we should have the same number of hydrogens as before
    assert len(new_hydrogens) == n_hydrogens

    v.draw_points(mol.get_coords(*new_hydrogens), colors="limegreen")
    v.draw_points(bam.utils.coord_array(*old_hydrogens), colors="purple")

    if base.ALLOW_VISUAL:
        v.show()

    assert len(mol._atoms) == ref_n


def test_optimize():
    mol = bam.Molecule.from_json(base.ROTSCAFJSON)
    mol.optimize(residue_graph=False)

    for b in mol.get_bonds():
        assert 0.9 < b.compute_length() < 2
    # mol.show()


def test_phosphorylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.phosphorylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_hydroxylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.hydroxylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_methylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.methylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_acetylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.acetylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_amidate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.amidate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_carboxylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.carboxylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_benzylate():
    mol = bam.Molecule.from_smiles("CC")
    n = mol.count_atoms()
    bam.benzylate(mol, "C1")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_phenolate():
    mol = bam.Molecule.from_smiles("CCC")
    n = mol.count_atoms()
    bam.phenolate(mol, "C1", how="ortho")
    bam.phenolate(mol, "C2", how="meta")
    bam.phenolate(mol, "C3", how="para")
    m = mol.count_atoms()
    assert m > n, "No atoms seem to have been added"
    if base.ALLOW_VISUAL:
        mol.show()


def test_carboxylate_multiple():
    mol = bam.Molecule.from_smiles("CCCCC")
    carbons = mol.get_atoms("C3", "C2", "C4")
    hydrogens = [mol.get_left_hydrogen(i) for i in carbons]
    hydrogens[0] = None
    bam.carboxylate(mol, carbons, hydrogens)
    if base.ALLOW_VISUAL:
        mol.show()


def test_get_left_right_hydrogens():
    mol = bam.Molecule.from_smiles("OC(=O)CC=O").autolabel()
    center = mol.get_atom("C2")
    left = mol.get_left_hydrogen(center)
    right = mol.get_right_hydrogen(center)
    assert left.element == "H"
    assert right.element == "H"
    if base.ALLOW_VISUAL:
        v = mol.draw()
        v.draw_point("center", center.coord, color="red")
        v.draw_point("left", left.coord, color="blue")
        v.draw_point("right", right.coord, color="green")
        v.show()


def test_get_atoms_filter_only():
    mol = bam.Molecule.from_compound("GLC")
    carbons = mol.get_atoms(filter=lambda x: x.element == "C")
    carbons = list(carbons)
    assert len(carbons) == 6


def test_get_atoms_residue_only():
    mol = bam.Molecule.from_compound("GLC")
    carbons = mol.get_atoms(residue=1)
    carbons = list(carbons)
    assert carbons == list(mol.get_residue(1).get_atoms())
    assert len(carbons) == 24


def test_superimpose():
    mol = bam.Molecule.from_compound("GLC")

    mol2 = mol.copy()
    mol2.transpose([0, 2, 12], 45, [1, 0, 0])
    mol3 = mol2.copy()
    mol4 = mol2.copy()

    v = mol.draw() + mol2.draw(atoms=False, line_color="red")

    mol2.superimpose_to_atom(
        mol2.get_atom("C1"),
        mol.get_atom("C1"),
    )

    v += mol2.draw(atoms=False, line_color="blue")

    mol3.superimpose_to_bond(
        mol3.get_bond("C1", "C2"),
        mol.get_bond("C1", "C2"),
    )

    v += mol3.draw(atoms=False, line_color="green")

    mol4.superimpose_to_triplet(
        ("C1", "O1", "HO1"),
        mol.get_atoms("C1", "O1", "HO1"),
    )

    v += mol4.draw(atoms=False, line_color="purple")
    if base.ALLOW_VISUAL:
        v.show()


def test_react_with():
    mol = bam.Molecule.from_compound("TYR")

    out = mol.react_with(
        mol, bam.structural.groups.carboxyl, bam.structural.groups.amine
    )
    assert out.count_residues() == 2
    if base.ALLOW_VISUAL:
        out.show()


def test_react_with2():
    mol = bam.Molecule.from_compound("TYR")

    out = mol.react_with(
        mol,
        bam.structural.groups.carboxyl,
        bam.structural.groups.hydroxyl,
        inplace=False,
    ).react_with(
        mol, bam.structural.groups.carboxyl, bam.structural.groups.amine, inplace=False
    )
    assert out.count_residues() == 3
    if base.ALLOW_VISUAL:
        out.show()


def test_from_geometry():
    P = bam.Atom.new("P")

    mol1 = bam.Molecule.from_geometry(bam.structural.geometry.trigonal_bipyramidal, [P])
    assert len(mol1.atoms) == 6
    assert len(mol1.bonds) == 5

    Fs = [bam.Atom.new("F") for _ in range(5)]
    mol2 = bam.Molecule.from_geometry(
        bam.structural.geometry.trigonal_bipyramidal, [P, *Fs]
    )

    assert len(mol2.atoms) == 6
    assert len(mol2.bonds) == 5

    mol3 = bam.Molecule.from_geometry(
        bam.structural.geometry.trigonal_bipyramidal, [P, *Fs[:3]]
    )

    if base.ALLOW_VISUAL:
        mol1.show()
        mol2.show()
        mol3.show()


def test_add_with_functional_groups():
    mol = bam.Molecule.from_compound("TYR")
    mol1, mol2 = mol.copy(2)

    out = mol1 % bam.structural.groups.carboxyl + mol2 % bam.structural.groups.amine
    assert out.count_residues() == 2
    if base.ALLOW_VISUAL:
        out.show()


def test_adjust_length_with_descendants():
    mol = bam.molecule("CC")
    bam.acetylate(mol, "C1")

    mol.adjust_bond_length("C", "C1", 4, True)
    assert 3.9 < mol.get_bond("C", "C1").compute_length() < 4.1
    if base.ALLOW_VISUAL:
        mol.show()


def test_get_atoms_generator_implementation():
    mol = bam.read_pdb(base.PADPDB)
    out = mol.get_atoms(
        "C",
        by="element",
        residue=1,
        filter=lambda x: int(x.id[1:]) % 2 == 0,
        keep_order=True,
    )
    assert len(out) != 0


# def test_sampler():
#     import numpy as np

#     bam.load_small_molecules()
#     N = 8
#     fragments = []
#     for id in bam.get_default_compounds().ids:
#         info, pdb = bam.get_compound(id, return_type="dict")
#         n_atoms = sum(1 for i in pdb["atoms"]["elements"] if i != "H")
#         if n_atoms <= N:
#             fragments.append(bam.get_compound(id, return_type="molecule"))
#             a = []
#             for adx, atom in enumerate(fragments[-1].get_atoms()):
#                 if fragments[-1].get_hydrogen(atom):
#                     a.append(adx)
#             if len(a) == 0:
#                 fragments.pop()
#                 continue

#     class MoleculeSampler:
#         def __init__(self, fragments):
#             attachment_points = []
#             atom_ids = []

#             to_drop = []
#             for fdx, fragment in enumerate(fragments):
#                 n_atoms = sum(1 for i in fragment.get_atoms() if i.element != "H")
#                 a = []
#                 for adx, atom in enumerate(fragment.get_atoms()):
#                     if atom.element == "H":
#                         continue
#                     if fragment.get_hydrogen(atom):
#                         a.append(adx)
#                 if len(a) == 0:
#                     to_drop.append(fdx)
#                     continue
#                 attachment_points.append(a)
#                 atom_ids.append([atom.id for atom in fragment.get_atoms()])

#             for fragment in to_drop:
#                 del fragments[fragment]

#             self.fragments = fragments
#             self.attachment_points = attachment_points
#             self.atom_ids = atom_ids
#             self._others = {
#                 i: np.array([j for j in range(len(fragments)) if j != i])
#                 for i in range(len(fragments))
#             }
#             self._n = 1

#         def random(self, n_fragments: int):
#             """
#             Make a random matrix encoding for a molecule assembled from fragments.
#             """
#             self._n = n_fragments
#             molecule_vec = np.full((n_fragments, 4), -1, dtype=int)
#             molecule_vec[0, 0] = np.random.choice(len(self.fragments))

#             _used_atoms = {i: set() for i in range(n_fragments)}
#             for i in range(1, n_fragments):
#                 molecule_vec[i, 0] = np.random.choice(len(self.fragments))

#                 while molecule_vec[i, 1] == -1:
#                     atom = np.random.choice(self.attachment_points[molecule_vec[i, 0]])
#                     if atom not in _used_atoms[i]:
#                         molecule_vec[i, 1] = atom
#                         _used_atoms[i].add(atom)

#                 while molecule_vec[i, 2] == -1:
#                     target = i - 1  # np.random.choice(range(i))
#                     available = [
#                         i
#                         for i in self.attachment_points[molecule_vec[target, 0]]
#                         if i not in _used_atoms[target]
#                     ]
#                     if len(available) > 0:
#                         molecule_vec[i, 2] = target
#                         molecule_vec[i, 3] = np.random.choice(available)
#                         _used_atoms[target].add(molecule_vec[i, 3])

#             return molecule_vec

#         def make_molecule(self, matrix):
#             """
#             Assemble a molecule based on an instruction matrix
#             """
#             _used_atoms = {i: set() for i in range(len(matrix))}

#             mol = self.fragments[matrix[0, 0]].copy()

#             for i in range(1, len(matrix)):
#                 source, source_atom, target, target_atom = matrix[i, :]

#                 if target_atom in _used_atoms[target]:
#                     raise ValueError("Target atom already used")
#                 if source_atom in _used_atoms[i]:
#                     raise ValueError("Source atom already used")

#                 link = bam.linkage(
#                     self.atom_ids[matrix[target, 0]][target_atom],
#                     self.atom_ids[source][source_atom],
#                 )
#                 mol.stitch_attach(
#                     self.fragments[source], link, at_residue=int(target + 1)
#                 )

#                 _used_atoms[target].add(target_atom)
#                 _used_atoms[i].add(source_atom)

#             return mol

#     sampler = MoleculeSampler(fragments)
#     for i in range(3):
#         m = sampler.random(4)
#         mol = sampler.make_molecule(m)
#         v = mol.draw(atoms=False)
#         v += (
#             mol.optimize(inplace=True)
#             .move([10, 0, 0])
#             .draw(atoms=False, line_color="blue")
#         )
#         if base.ALLOW_VISUAL:
#             v.show()


def test_search_by_constraints():
    mol = bam.read_smiles("CC(O)=CC")
    constraints = bam.structural.neighbors.constraints
    out = mol.search_by_constraints(
        [
            constraints.multi_constraint(
                constraints.has_element("C"),
                constraints.has_neighbor_hist({"H": 3}),
            ),
        ]
    )
    assert len(out) == 2


def test_search_by_constraints2():
    mol = bam.read_smiles("CC(O)=CC")
    constraints = bam.structural.neighbors.constraints
    out = mol.search_by_constraints(
        [
            constraints.multi_constraint(
                constraints.has_element("C"),
                constraints.has_neighbor_hist({"O": 1, "C": 2}),
            ),
            constraints.has_element("O"),
        ]
    )
    assert len(out) == 1


def test_add_hydrogens_specific():
    glc = bam.Molecule.from_compound("GLC")

    glc.remove_hydrogens()
    glc.add_hydrogens("O1")

    assert len(glc.get_atoms("H", by="element")) == 1

    glc.add_hydrogens("C3", "C4", "C5")

    assert glc.get_hydrogen("C3") is not None


def test_remove_hydrogens_specific():
    glc = bam.Molecule.from_compound("GLC")

    current = len(glc.get_atoms("H", by="element"))
    glc.remove_hydrogens("O1")
    assert len(glc.get_atoms("H", by="element")) == current - 1

    glc.remove_hydrogens("C3", "C4", "C5")
    assert len(glc.get_atoms("H", by="element")) == current - 4
    assert glc.get_hydrogen("C3") is None


def test_set_charge_with_protonation():
    mol = bam.Molecule.from_compound("GLC")

    n_hydrogens = len(mol.get_atoms("H", by="element"))
    mol.set_charge("O1", 1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens + 1

    mol.set_charge("O1", 0)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens

    mol.set_charge("O1", -1)
    assert len(mol.get_atoms("H", by="element")) == n_hydrogens - 1

    try:
        mol.set_charge("O1", -2)
    except ValueError:
        pass
    else:
        assert False, "Should have raised an error"


def test_isomers():
    # GAL and GIV are L/D isomers of galactose
    glc = bam.Molecule.from_compound("GLC")
    gal = bam.Molecule.from_compound("GAL")
    giv = bam.Molecule.from_compound("GIV")

    bam.dont_use_ic()

    A = glc % "14bb" + gal
    B = glc % "14bb" + giv

    A.superimpose_to_residue(1, B.get_residue(1))

    b_dihe = B.compute_dihedral(
        B.get_atom("C5", residue=1),
        B.get_atom("C4", residue=1),
        B.get_atom("O4", residue=1),
        B.get_atom("C1", residue=2),
    )

    a_dihe = A.compute_dihedral(
        A.get_atom("C5", residue=1),
        A.get_atom("C4", residue=1),
        A.get_atom("O4", residue=1),
        A.get_atom("C1", residue=2),
    )

    A.rotate_descendants(
        A.get_atom("C4", residue=1),
        A.get_atom("O4", residue=1),
        180 - a_dihe + b_dihe,
    )

    assert np.allclose(A.get_coords(residue=1), B.get_coords(residue=1))
    assert not np.allclose(A.get_coords(residue=2), B.get_coords(residue=2))

    if base.ALLOW_VISUAL:
        v = A.draw()
        B.move([10, 0, 0])
        v += B.draw(atoms=True, line_color="red")
        v.show()
        pass
