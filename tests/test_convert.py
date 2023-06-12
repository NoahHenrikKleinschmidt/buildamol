"""
Tests to check the conversion of pybel to biopython
"""

import numpy as np

from openbabel import pybel
import Bio.PDB as bio

import biobuild as bb
import base


def test_class():
    converter = bb.utils.convert.PybelBioPythonConverter()
    assert converter is not None


def test_convert_molecule():
    _pybel = pybel.readfile("pdb", base.MANNOSE)
    _pybel = next(_pybel)

    _biopython = bio.PDBParser().get_structure("MAN", base.MANNOSE)

    converter = bb.utils.convert.PybelBioPythonConverter()
    _converted = converter.convert(_pybel)
    assert _converted is not None

    assert len(list(_converted.get_atoms())) == len(list(_biopython.get_atoms()))

    _ref_atoms = [(i.get_parent().id[1], i.coord) for i in _biopython.get_atoms()]
    for residue in _converted.get_residues():
        for atom in residue.get_atoms():
            for r in _ref_atoms:
                if r[0] == residue.id[1]:
                    if np.abs(np.sum((atom.coord + r[1]) - 2 * atom.coord)) < 0.01:
                        break
            else:
                raise AssertionError("Atom not found")
