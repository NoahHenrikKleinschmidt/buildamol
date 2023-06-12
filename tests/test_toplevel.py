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
