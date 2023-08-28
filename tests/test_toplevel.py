"""
Tests to check the behaviour of the bb.Molecule object
"""

import os
from copy import deepcopy
import numpy as np
import biobuild as bb
import Bio.PDB as bio

import tests.base as base


def test_molecule():
    bb.load_sugars()
    mol = bb.molecule("GLC")
    assert isinstance(mol, bb.Molecule)

    mol = bb.molecule("2-acetamido-2-deoxy-beta-D-glucopyranose")
    assert isinstance(mol, bb.Molecule)

    mol = bb.molecule("CC(=O)NC1C(CC(OC1C(C(CO)O)O)(C(=O)O)O)O")
    assert isinstance(mol, list)
    assert len(mol) == 2
    assert isinstance(mol[0], bb.Molecule)
    assert isinstance(mol[1], bb.Molecule)
    bb.unload_sugars()

    assert bb.utils.auxiliary.HAS_RDKIT == True, "RDKit is not installed!"
    mol2 = bb.read_smiles("CC(=O)NC1C(CC(OC1C(C(CO)O)O)(C(=O)O)O)O")

    assert mol2 is not None
    assert mol2.count_atoms() == mol[0].count_atoms()
    assert mol2.count_bonds() == mol[0].count_bonds()


def test_polymerize():
    bb.load_sugars()
    glc = bb.molecule("GLC")
    glc5 = bb.polymerize(glc, 5, "14bb")
    assert glc5 is not glc
    assert len(glc5.residues) == 5
    bb.unload_sugars()


def test_connect():
    bb.load_sugars()
    glc1 = bb.molecule("GLC")
    glc2 = bb.molecule("GLC")
    out = bb.connect(glc1, glc2, "14bb", 1, 1)
    assert out is not glc1
    assert out is not glc2
    assert len(out.residues) == 2
    bb.unload_sugars()
