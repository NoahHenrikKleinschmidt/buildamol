"""
Tests to check the behaviour of the bam.Molecule object
"""

import os
from copy import deepcopy
import numpy as np
import buildamol as bam
import Bio.PDB as bio

import tests.base as base


def test_molecule():
    bam.load_sugars()
    mol = bam.molecule("GLC")
    assert isinstance(mol, bam.Molecule)

    mol = bam.molecule("2-acetamido-2-deoxy-beta-D-glucopyranose")
    assert isinstance(mol, bam.Molecule)

    mol = bam.molecule("CC(=O)NC1C(CC(OC1C(C(CO)O)O)(C(=O)O)O)O")
    assert isinstance(mol, list)
    assert len(mol) == 2
    assert isinstance(mol[0], bam.Molecule)
    assert isinstance(mol[1], bam.Molecule)
    bam.unload_sugars()

    assert bam.utils.auxiliary.HAS_RDKIT == True, "RDKit is not installed!"
    mol2 = bam.read_smiles("CC(=O)NC1C(CC(OC1C(C(CO)O)O)(C(=O)O)O)O")

    assert mol2 is not None
    assert mol2.count_atoms() == mol[0].count_atoms()
    assert mol2.count_bonds() == mol[0].count_bonds()


def test_polymerize():
    bam.load_sugars()
    glc = bam.molecule("GLC")
    glc5 = bam.polymerize(glc, 5, "14bb")
    assert glc5 is not glc
    assert len(glc5.residues) == 5
    bam.unload_sugars()


def test_connect():
    bam.load_sugars()
    glc1 = bam.molecule("GLC")
    glc2 = bam.molecule("GLC")
    out = bam.connect(glc1, glc2, "14bb", 1, 1)
    assert out is not glc1
    assert out is not glc2
    assert len(out.residues) == 2
    bam.unload_sugars()


def test_reaction_can_apply_with_single_target():
    A = bam.molecule("OCCO")
    B = bam.molecule("c1ccccc1C")

    R = bam.Reaction(
        atom1=lambda mol: mol.get_atoms(
            "C",
            by="element",
            filter=bam.structural.constraints_v2.has_single_bond_with("O"),
        ).pop(),
        atom2=lambda mol: 1,
        delete_in_target=lambda atom: atom.get_neighbors(
            filter=bam.structural.constraints_v2.has_element("O")
        ),
        delete_in_source=lambda atom: atom.get_hydrogens().pop(),
    )
    assert R.can_apply(A, B)

    out1 = R(A, B)

    atom1 = A.get_atoms(
        "C",
        by="element",
        filter=bam.structural.constraints_v2.has_single_bond_with("O"),
    ).pop()
    atom2 = B.get_atom(1)

    R = bam.Reaction(
        atom1=atom1,
        atom2=atom2,
        delete_in_target=atom1.get_neighbors(
            filter=bam.structural.constraints_v2.has_element("O")
        ).pop(),
        delete_in_source=atom2.get_hydrogens().pop(),
    )
    assert R.can_apply(A, B)
    out2 = R(A, B)
    assert out1.to_smiles() == out2.to_smiles()


def test_reaction_can_apply_with_multiple_targets():
    A = bam.molecule("OCCO")
    B = bam.molecule("c1ccccc1C")

    R = bam.Reaction(
        atom1=lambda mol: mol.get_atoms(
            "C",
            by="element",
            filter=bam.structural.constraints_v2.has_single_bond_with("O"),
        ),
        atom2=lambda mol: 1,
        delete_in_target=lambda atom: atom.get_neighbors(
            filter=bam.structural.constraints_v2.has_element("O")
        ),
    )
    assert R.can_apply(A, B)

    out1 = R(A, B)

    atom1 = A.get_atoms(
        "C",
        by="element",
        filter=bam.structural.constraints_v2.has_single_bond_with("O"),
    )
    atom2 = B.get_atom(1)

    R = bam.Reaction(
        atom1=atom1,
        atom2=atom2,
        delete_in_target=lambda atom: atom.get_neighbors(
            filter=bam.structural.constraints_v2.has_element("O")
        ).pop(),
        delete_in_source=atom2.get_hydrogens().pop(),
    )
    assert R.can_apply(A, B)
    out2 = R(A, B)
    assert out1.to_smiles() == out2.to_smiles()
