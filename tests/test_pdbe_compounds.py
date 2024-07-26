"""
Tests for the PDBe compounds class.
"""

import tests.base as base
import Bio.PDB as bio
import buildamol as bam
from buildamol.resources import pdbe_compounds
import numpy as np

bam.load_sugars()


def test_from_cif():
    current = len(bam.get_default_compounds())
    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    assert comps is not None, "Could not load the PDBe compounds from a CIF file."
    assert len(comps.ids) != 0, "The number of compounds is not correct."
    assert (
        len(bam.get_default_compounds()) == current
    ), "The compounds were added to the default compounds!"


def test_getting_compounds():
    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man = comps.get("D-Mannose", by="name")
    assert man is not None

    man = comps.get("C6 H12 O6", by="formula")
    assert man is not None

    man = comps.get("OC1C(O)C(OC(O)C1O)CO", by="smiles")
    assert man is not None

    man = comps.get("MAN")
    assert man is not None

    assert len(man.atoms) == 24
    assert len(man.bonds) == 24


def test_compound_is_same():
    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man = comps.get("D-Mannose", by="name")

    ref = bam.Molecule.from_pdb(base.MANPDB)
    ref.infer_bonds()

    man_atoms = [(i.id, i.serial_number) for i in man.atoms]
    ref_atoms = [(i.id, i.serial_number) for i in ref.atoms]
    assert man_atoms == ref_atoms

    man_bonds = set((i.id, j.id) for i, j in man.bonds)
    ref_bonds = set((i.id, j.id) for i, j in ref.bonds)

    assert man_bonds == ref_bonds


def test_compound_getting_types():
    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    man_mol = comps.get("D-Mannose", by="name")
    assert isinstance(man_mol, bam.Molecule)

    man_dict = comps.get("D-Mannose", by="name", return_type="dict")
    assert isinstance(man_dict, dict)

    man_struct = comps.get("D-Mannose", by="name", return_type="structure")
    assert isinstance(man_struct, bio.Structure.Structure)
    assert len(list(man_struct.get_atoms())) == 24

    man_res = comps.get("D-Mannose", by="name", return_type="residue")
    assert isinstance(man_res, bio.Residue.Residue)
    assert len(list(man_res.get_atoms())) == 24

    glc_mol = comps.get("Glucose", by="name")
    assert isinstance(glc_mol, list)
    assert len(glc_mol) == 2
    assert isinstance(glc_mol[0], bam.Molecule)

    glc_dict = comps.get("Glucose", by="name", return_type="dict")
    assert isinstance(glc_dict, list)
    assert len(glc_dict) == 2
    assert isinstance(glc_dict[0], dict)

    glc_struct = comps.get("Glucose", by="name", return_type="structure")
    assert isinstance(glc_struct, list)
    assert len(glc_struct) == 2
    assert isinstance(glc_struct[0], bio.Structure.Structure)
    assert len(list(glc_struct[0].get_atoms())) == 24


def test_relabel():
    bam.load_sugars()
    scrambled = bam.Molecule.from_compound("MAN")

    # randomly rotate the molecule
    random_vector = np.random.rand(3) * 10

    # choose some rotation axis
    axis = np.random.randint(0, len(scrambled.bonds) - 1)
    axis = scrambled.bonds[axis]

    # relabel to elementwise order without specific connectivity
    counts = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}

    for atom in scrambled.atoms:
        counts[atom.element] += 1
        atom.id = atom.element + str(counts[atom.element])
        atom.coord += random_vector

    scrambled.rotate_around_bond(*axis, np.random.randint(0, 180))

    # now relabel the atoms to match the original
    old_scrambled = [i.id for i in scrambled.atoms]
    old_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    comps = bam.resources.get_default_compounds()
    comps.relabel_atoms(scrambled)

    new_scrambled = [i.id for i in scrambled.atoms]
    new_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    assert "H12" in old_scrambled
    assert "H61" not in old_scrambled
    assert "H12" not in new_scrambled
    assert "H61" in new_scrambled
    assert np.allclose(old_scrambled_coords, new_scrambled_coords)

    # v = bam.utils.visual.MoleculeViewer3D(scrambled)
    # scrambled.show()
    bam.unload_sugars()


def test_relabel_2():
    bam.unload_all_compounds()
    comps = bam.resources.get_default_compounds()
    assert (
        len(comps) == 0
    ), f"The number of compounds is not correct: {len(comps)} should be empty"
    bam.load_sugars()
    assert (
        len(comps) == 1068
    ), f"The number of compounds is not correct: {len(comps)} should be filled with 1068 sugar compounds"

    for i in ("MAN", "GLC", "BMA", "FUC"):
        scrambled = bam.Molecule.from_compound(i)

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

        # now relabel the atoms to match the original
        old_scrambled = set(i.id for i in scrambled.atoms)
        old_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

        comps.relabel_atoms(scrambled)

        new_scrambled = set(i.id for i in scrambled.atoms)
        new_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

        assert scrambled.get_bonds("C1", "C2") != []
        assert scrambled.get_bonds("C6", "H61") != []
        assert old_scrambled.difference(new_scrambled) != set()
        assert np.allclose(old_scrambled_coords, new_scrambled_coords)

    bam.unload_sugars()
    # v = bam.utils.visual.MoleculeViewer3D(scrambled)
    # v.show()


def test_get_2FJ():
    assert bam.has_compound("2FJ") == False, "The compound was already added!"

    comps = bam.read_compounds(str(base.PDBE_TEST_FILE), set_default=False)
    _dict = comps.get("2FJ", return_type="dict")
    assert isinstance(_dict, dict)
    mol = comps.get("2FJ", return_type="molecule")
    assert isinstance(mol, bam.Molecule)
    assert (
        comps.has_residue("2FJ") == True
    ), "The compound was not added to the compounds!"
    assert (
        bam.has_compound("2FJ") == False
    ), "The compound was added to the default compounds!"
    assert (
        bam.get_default_compounds().has_residue("2FJ") == False
    ), "The compound was added to the default compounds!"


def test_get_all_molecule():
    bam.unload_all_compounds()
    bam.load_sugars()
    comps = bam.get_default_compounds()
    assert len(comps) != 0, "No compounds were loaded!"
    assert len(comps) == 1068, "The number of compounds is not correct!"
    for comp, d_data, d_pdb in comps:
        try:
            assert isinstance(comp, str)
            assert isinstance(d_data, dict)
            assert isinstance(d_pdb, dict)
            mol = comps.get(comp)
            assert isinstance(mol, bam.Molecule)
            assert sum(1 for i in mol.get_atoms()) == len(d_pdb["atoms"]["ids"])
            assert sum(1 for i in mol.get_bonds()) == len(d_pdb["bonds"]["bonds"])
        except StopIteration as e:
            w = Warning(f"Failed for {comp}: {e}")
            print(w)
    bam.unload_all_compounds()


def test_get_compounds_preserve_coords():
    bam.load_sugars()
    A = bam.get_compound("GLC")
    B = bam.get_compound("GLC")

    A.move((10, 0, 0))
    assert not np.allclose(A.get_atom("C1").coord, B.get_atom("C1").coord)
    C = bam.get_compound("GLC")
    assert np.allclose(C.get_atom("C1").coord, B.get_atom("C1").coord)
