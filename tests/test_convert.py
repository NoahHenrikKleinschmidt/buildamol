"""
Tests to check the conversion of pybel to biopython
"""

import numpy as np

from openbabel import pybel
import Bio.PDB as bio

import buildamol as bam
import tests.base as base


def test_biopython():
    mol = bam.read_smiles("C1=CC=CC=C1")
    assert mol is not None

    _biopython = mol.to_biopython()
    assert _biopython is not None

    assert sum(1 for i in _biopython.get_atoms()) == mol.count_atoms()

    reverse = bam.Molecule(_biopython)
    assert reverse is not None

    assert reverse.count_atoms() == 12
    assert reverse.count_bonds() == 0


def test_openbabel():
    mol = bam.read_smiles("C=CC=CC=C")
    assert mol is not None

    _pybel = mol.to_pybel()
    assert _pybel is not None

    assert _pybel.OBMol.NumAtoms() == mol.count_atoms()
    assert _pybel.OBMol.NumBonds() == mol.count_bonds()

    reverse = bam.Molecule.from_pybel(_pybel)
    assert reverse is not None

    assert reverse.count_atoms() == 14
    assert reverse.count_bonds() == 13
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 3


def test_rdkit():
    mol = bam.read_smiles("C1=CC=CC=C1")
    assert mol is not None

    rdmol = mol.to_rdkit()
    assert rdmol is not None

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()

    reverse = bam.Molecule.from_rdkit(rdmol)
    assert reverse is not None

    assert reverse.count_atoms() == 12
    assert reverse.count_bonds() == 12
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 3


def test_rdkit2():
    mol = bam.read_pdb(base.GLCPDB)
    mol.infer_bonds()
    assert mol is not None

    rdmol = mol.to_rdkit()
    assert rdmol is not None

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()

    reverse = bam.Molecule.from_rdkit(rdmol)
    assert reverse is not None

    assert reverse.count_atoms() == mol.count_atoms()
    assert reverse.count_bonds() == mol.count_bonds()

    mol.remove_atoms("H1", "H2")
    mol.get_bond("C1", "C2").double()

    rdmol = mol.to_rdkit()

    assert sum(1 for bond in rdmol.GetBonds() if bond.GetBondTypeAsDouble() == 2) == 1

    reverse = bam.Molecule.from_rdkit(rdmol)

    assert reverse.count_atoms() == mol.count_atoms()
    assert reverse.count_bonds() == mol.count_bonds()
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 1


def test_rdkit3():
    mol = bam.Molecule.from_pubchem("benzoic acid")
    mol.autolabel()
    assert mol is not None

    rdmol = mol.to_rdkit()
    assert rdmol is not None

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()

    reverse = bam.Molecule.from_rdkit(rdmol)
    assert reverse is not None

    assert reverse.count_atoms() == mol.count_atoms()
    assert reverse.count_bonds() == mol.count_bonds()
    # 3 aromatic + C=O
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 4

    conv = bam.utils.convert.RDKITBiopythonConverter()

    conv.molecule_to_pdbio(mol)
    classical_rdkit = conv._pdbio_to_rdkit()
    assert classical_rdkit is not None

    conv._rdkit_to_pdbio(rdmol)
    classical = bam.Molecule.from_pdb(conv.__fileio__, id="classical")
    assert classical is not None

    classical.move((0, 0, 5))

    if base.ALLOW_VISUAL:
        v = classical.draw()
        for bond in reverse.get_bonds():
            v.draw_vector(
                None, bond[0].coord, bond[1].coord, color="red", linewidth=bond.order**3
            )
        v.show()


def test_rdkit4():
    mol = bam.Molecule.from_pubchem("GlcNAc")
    mol.autolabel()
    assert mol is not None

    rdmol = mol.to_rdkit()
    assert rdmol is not None

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()

    reverse = bam.Molecule.from_rdkit(rdmol)
    assert reverse is not None

    assert reverse.count_atoms() == mol.count_atoms()
    assert reverse.count_bonds() == mol.count_bonds()
    # C=O
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 1

    conv = bam.utils.convert.RDKITBiopythonConverter()

    conv.molecule_to_pdbio(mol)
    classical_rdkit = conv._pdbio_to_rdkit()
    assert classical_rdkit is not None

    conv._rdkit_to_pdbio(rdmol)
    classical = bam.Molecule.from_pdb(conv.__fileio__, id="classical")
    assert classical is not None

    classical.move((0, 0, 5))

    if base.ALLOW_VISUAL:
        v = classical.draw()
        for bond in reverse.get_bonds():
            v.draw_vector(
                None, bond[0].coord, bond[1].coord, color="red", linewidth=bond.order**3
            )
        v.show()


def test_to_rdkit5():

    mol = bam.Molecule.from_pdb(base.EX8PDB)

    rdmol = mol.to_rdkit()

    ref = bam.utils.Chem.MolFromPDBFile(
        base.EX8PDB,
        proximityBonding=False,
        removeHs=False,
    )

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()
    assert ref.GetNumAtoms() == mol.count_atoms()
    assert ref.GetNumBonds() == mol.count_bonds()

    for atom, atom_ref in zip(rdmol.GetAtoms(), ref.GetAtoms()):
        assert atom.GetAtomicNum() == atom_ref.GetAtomicNum()
        assert atom.GetFormalCharge() == atom_ref.GetFormalCharge()
        assert atom.GetIsAromatic() == atom_ref.GetIsAromatic()
        assert atom.GetHybridization() == atom_ref.GetHybridization()
        assert atom.GetTotalNumHs() == atom_ref.GetTotalNumHs()
        assert atom.GetTotalValence() == atom_ref.GetTotalValence()
        assert atom.GetIdx() == atom_ref.GetIdx()
        assert atom.GetDegree() == atom_ref.GetDegree()

    for bond, bond_ref in zip(rdmol.GetBonds(), ref.GetBonds()):
        assert bond.GetBondType() == bond_ref.GetBondType()
        assert bond.GetBeginAtomIdx() == bond_ref.GetBeginAtomIdx()
        assert bond.GetEndAtomIdx() == bond_ref.GetEndAtomIdx()
        assert bond.GetIdx() == bond_ref.GetIdx()
        assert bond.GetIsAromatic() == bond_ref.GetIsAromatic()
        assert bond.GetIsConjugated() == bond_ref.GetIsConjugated()
        assert bond.GetStereo() == bond_ref.GetStereo()


def test_stk():
    mol = bam.Molecule.from_pdb(base.GLCPDB)
    mol.infer_bonds()

    stk = mol.to_stk()
    assert stk is not None

    reverse = bam.Molecule.from_stk(stk)
    assert reverse is not None

    assert reverse.count_atoms() == mol.count_atoms()
    assert reverse.count_bonds() == mol.count_bonds()
