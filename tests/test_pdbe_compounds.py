"""
Tests for the PDBe compounds class.
"""

import tests.base as base
import Bio.PDB as bio
import biobuild as bb
from biobuild.resources import pdbe_compounds
import numpy as np

bb.load_sugars()


def test_from_cif():
    comps = pdbe_compounds.PDBECompounds.from_file(base.PDBE_TEST_FILE)

    assert comps is not None, "Could not load the PDBe compounds from a CIF file."
    assert len(comps.ids) != 0, "The number of compounds is not correct."


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

    ref = bb.Molecule.from_pdb(base.MANNOSE)
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
    assert isinstance(man_mol, bb.Molecule)

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
    assert isinstance(glc_mol[0], bb.Molecule)

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
    scrambled = bb.Molecule.from_compound("MAN")

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

    comps = bb.resources.get_default_compounds()
    comps.relabel_atoms(scrambled)

    new_scrambled = [i.id for i in scrambled.atoms]
    new_scrambled_coords = np.array([i.coord for i in scrambled.atoms])

    assert "H12" in old_scrambled
    assert "H61" not in old_scrambled
    assert "H12" not in new_scrambled
    assert "H61" in new_scrambled
    assert np.allclose(old_scrambled_coords, new_scrambled_coords)

    # v = bb.utils.visual.MoleculeViewer3D(scrambled)
    scrambled.show()


def test_relabel_2():
    comps = bb.resources.get_default_compounds()

    for i in ("MAN", "GLC", "BMA", "FUC"):
        scrambled = bb.Molecule.from_compound(i)

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

        # v = bb.utils.visual.MoleculeViewer3D(scrambled)
        # v.show()


def test_get_3BU():
    comps = bb.read_compounds(base.PDBE_TEST_FILE)
    _dict = comps.get("3BU", return_type="dict")
    assert isinstance(_dict, dict)
    mol = comps.get("3BU", return_type="molecule")
    assert isinstance(mol, bb.Molecule)


def test_get_all_molecule():
    comps = bb.get_default_compounds()
    for comp, d_data, d_pdb in comps:
        try:
            assert isinstance(comp, str)
            assert isinstance(d_data, dict)
            assert isinstance(d_pdb, dict)
            mol = comps.get(comp)
            assert isinstance(mol, bb.Molecule)
            assert sum(1 for i in mol.get_atoms()) == len(d_pdb["atoms"]["ids"])
            assert sum(1 for i in mol.get_bonds()) == len(d_pdb["bonds"])
        except StopIteration as e:
            w = Warning(f"Failed for {comp}: {e}")
            print(w)
