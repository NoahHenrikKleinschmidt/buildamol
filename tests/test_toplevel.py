"""
Tests to check the behaviour of the bb.Molecule object
"""

import os
from copy import deepcopy
import numpy as np
import biobuild as bb
import Bio.PDB as bio

import base


def test_molecule():
    mol = bb.molecule("GLC")
    assert isinstance(mol, bb.Molecule)

    mol = bb.molecule("2-acetamido-2-deoxy-beta-D-glucopyranose")
    assert isinstance(mol, bb.Molecule)

    mol = bb.molecule("CC(=O)NC1C(CC(OC1C(C(CO)O)O)(C(=O)O)O)O")
    assert isinstance(mol, bb.Molecule)


def test_polymerize():
    glc = bb.molecule("GLC")
    glc5 = bb.polymerize(glc, 5, "14bb")
    assert glc5 is not glc
    assert len(glc5.residues) == 5


def test_connect():
    glc1 = bb.molecule("GLC")
    glc2 = bb.molecule("GLC")
    out = bb.connect(glc1, glc2, "14bb", 1, 1)
    assert out is not glc1
    assert out is not glc2
    assert len(out.residues) == 2
