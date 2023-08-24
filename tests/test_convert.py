"""
Tests to check the conversion of pybel to biopython
"""

import numpy as np

from openbabel import pybel
import Bio.PDB as bio

import biobuild as bb
import tests.base as bas


def test_biopython():
    mol = bb.read_smiles("C1=CC=CC=C1")
    assert mol is not None

    _biopython = mol.to_biopython()
    assert _biopython is not None

    assert sum(1 for i in _biopython.get_atoms()) == mol.count_atoms()

    reverse = bb.Molecule(_biopython)
    assert reverse is not None

    assert reverse.count_atoms() == 12
    assert reverse.count_bonds() == 0


def test_openbabel():
    mol = bb.read_smiles("C=CC=CC=C")
    assert mol is not None

    _pybel = mol.to_pybel()
    assert _pybel is not None

    assert _pybel.OBMol.NumAtoms() == mol.count_atoms()
    assert _pybel.OBMol.NumBonds() == mol.count_bonds()

    reverse = bb.Molecule.from_pybel(_pybel)
    assert reverse is not None

    assert reverse.count_atoms() == 14
    assert reverse.count_bonds() == 13
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 3


def test_rdkit():
    mol = bb.read_smiles("C1=CC=CC=C1")
    assert mol is not None

    rdmol = mol.to_rdkit()
    assert rdmol is not None

    assert rdmol.GetNumAtoms() == mol.count_atoms()
    assert rdmol.GetNumBonds() == mol.count_bonds()

    reverse = bb.Molecule.from_rdkit(rdmol)
    assert reverse is not None

    assert reverse.count_atoms() == 12
    assert reverse.count_bonds() == 12
    assert sum(1 for bond in reverse.get_bonds() if bond.order == 2) == 3
