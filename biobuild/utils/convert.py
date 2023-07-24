"""
Functions to convert different data formats
"""

# raise DeprecationWarning("This module is currently dropped...")

import tempfile
import os
import numpy as np

from periodictable import elements

import Bio.PDB as bio

try:
    from openbabel import pybel
except ImportError:
    pybel = None

try:
    from rdkit import Chem
except ImportError:
    Chem = None

import biobuild.utils.auxiliary as aux
import biobuild.utils.defaults as defaults


class PybelBioPythonConverter:
    """
    Convert Pybel (openbabel) data structures to Biopython
    """

    def __init__(self) -> None:
        self.__fileio__ = tempfile.mktemp(suffix=".pdb")

    __element_counts__ = {}
    _current_residue = None

    def convert(self, obj, full_structure: bool = False):
        """
        Convert between Biopython and Pybel (openbabel) data structures

        Parameters
        ----------
        obj : object
            The object to convert. The following conversions are supported:
            - A pybel `Molecule` -> biopython `Structure`
            - A pybel `Residue` -> biopython `Residue`
            - A pybel `Atom` -> biopython `Atom`
            - A pybel `Residue` -> biopython `Structure` (if `full_structure=True`)
            - A pybel `Atom` -> biopython `Structure` (if `full_structure=True`)

        full_structure : bool, optional
            If `True` a new full structure object is made instead of simply
            the hierarchical class equivalent.

        Returns
        -------
        object
            The converted object
        """
        if isinstance(obj, pybel.Molecule):
            return self.pybel_molecule_to_biopython(obj)

        elif isinstance(obj, pybel.Residue):
            residue = self.pybel_residue_to_biopython(obj)
            if full_structure:
                struct = self._new_biopython_structure()
                struct[0]["A"].add(residue)
                for a in residue.get_atoms():
                    residue.detach_child(a.id)
                    residue.add(a)
                return struct
            return residue

        elif isinstance(obj, pybel.Atom):
            atom = self.pybel_atom_to_biopython(obj)
            if full_structure:
                struct = self._new_biopython_structure()
                residue = bio.Residue.Residue(("H_UNL", 1, " "), "UNL", " ")
                residue.add(atom)
                struct[0]["A"].add(atom.get_parent())
                return struct
            return atom

        else:
            raise ValueError(f"Cannot convert object of type {type(obj)}")

    @staticmethod
    def _new_biopython_structure():
        """
        Create a new biopython structure object with a Model and Chain
        """
        struct = bio.Structure.Structure("New Molecule")
        struct.add(bio.Model.Model(0))
        struct[0].add(bio.Chain.Chain("A"))
        return struct

    def pybel_atom_to_biopython(self, obj):
        """
        Convert a pybel `Atom` to a biopython `Atom`

        Parameters
        ----------
        obj : pybel.Atom
            The pybel `Atom` to convert

        Returns
        -------
        bio.Atom.Atom
            The converted biopython `Atom`
        """

        # just for testing purposes...
        if obj.atomicnum == 67:
            element = "H"
        else:
            element = str(elements[obj.atomicnum])

        # # due to the issue of Pybel reading the CHARMM force-field
        # # designation of hydrogens by ids relating to their connectivity
        # # e.g. HO2 = Hydrogen at Oxygen2, the element is falsely inferred
        # # as Holmium (Ho) instead of H, so we need to manually correct for
        # # this - thereby, of course, destroying biobuild capacity to work
        # # with Holmium, should anybody actually wish to use it...
        # if element == "Ho":
        #     element = "H"

        if not element in self.__element_counts__:
            self.__element_counts__[element] = {None: 0}

        self.__element_counts__[element].setdefault(self._current_residue, 0)
        self.__element_counts__[element][self._current_residue] += 1
        name = f"{element}{self.__element_counts__[element][self._current_residue]}"

        new = bio.Atom.Atom(
            name=name,
            coord=np.asarray(obj.coords),
            serial_number=obj.idx,
            bfactor=0.0,
            occupancy=0.0,
            altloc=" ",
            fullname=name,
            element=element,
        )
        return new

    def pybel_residue_to_biopython(self, obj):
        """
        Convert a pybel `Residue` to a biopython `Residue`

        Parameters
        ----------
        obj : pybel.Residue
            The pybel `Residue` to convert

        Returns
        -------
        bio.Residue.Residue
            The converted biopython `Residue`
        """
        self._current_residue = (obj.idx, obj.name)
        new = bio.Residue.Residue(
            id=(f"H_{obj.name}", obj.idx + 1, " "), resname=obj.name, segid=" "
        )
        for atom in obj.atoms:
            atom = self.pybel_atom_to_biopython(atom)
            new.add(atom)

        self._current_residue = None
        return new

    def pybel_molecule_to_biopython(self, obj):
        """
        Convert a pybel `Molecule` to a biopython `Structure`

        Parameters
        ----------
        obj : pybel.Molecule
            The pybel `Molecule` to convert

        Returns
        -------
        bio.Structure.Structure
            The converted biopython `Structure`
        """
        if os.path.exists(obj.title):
            id = aux.filename_to_id(obj.title)
        else:
            id = obj.title

        new = bio.Structure.Structure(id)
        new.add(bio.Model.Model(0))
        new[0].add(bio.Chain.Chain("A"))

        if len(obj.residues) == 0:
            residue = bio.Residue.Residue(("H_UNL", 1, " "), "UNL", " ")
            new[0]["A"].add(residue)
            for atom in obj.atoms:
                residue.add(self.pybel_atom_to_biopython(atom))
        else:
            for residue in obj.residues:
                residue = self.pybel_residue_to_biopython(residue)
                new[0]["A"].add(residue)

        return new

    def biobuild_to_pybel(self, obj):
        """
        Convert a biobuild object to a pybel object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if type(obj).__name__ == "Molecule":
            obj.to_pdb(self.__fileio__)
        elif type(obj).__name__ == "Structure":
            io = bio.PDBIO()
            io.set_structure(obj.to_biopython())
            io.save(self.__fileio__)
        else:
            while not type(obj).__name__ in ["Molecule", "Structure"]:
                obj = obj.get_parent()
            return self.biobuild_to_pybel(obj)

        # read the file back in with pybel
        pybel_obj = next(pybel.readfile("pdb", self.__fileio__), None)
        if not pybel_obj:
            raise ValueError("Could not convert to pybel object")
        return pybel_obj


class RDKITBiopythonConverter:
    """
    A class to convert between RDKit and biopython objects
    This is done simply by using a temporary PDB pseudo-file
    """

    def __init__(self):
        self.__fileio__ = tempfile.mktemp(suffix=".pdb")

    def rdkit_to_biopython(self, obj):
        """
        Convert an RDKit object to a biopython object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not is_rdkit(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        mol = self.rdkit_to_pdbio(obj)
        new = self.pdbio_to_biopython(mol)
        os.remove(self.__fileio__)
        return new

    def biopython_to_rdkit(self, obj):
        """
        Convert a biopython object to an RDKit object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not is_biopython(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        self.biopython_to_pdbio(obj)
        new = self.pdbio_to_rdkit()
        os.remove(self.__fileio__)
        return new

    def rdkit_to_pdbio(self, obj):
        """
        Convert an RDKit object to a PDB file

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        str
            The PDB file
        """
        if not is_rdkit(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        Chem.MolToPDBFile(obj, self.__fileio__)

    def pdbio_to_rdkit(self) -> "Chem.rdchem.Mol":
        """
        Convert the internal FileIO to an RDKit Mol object

        Returns
        -------
        Mol
            The converted object
        """
        mol = Chem.MolFromPDBFile(self.__fileio__, proximityBonding=False)
        if mol is None:
            raise ValueError("Could not convert PDB file to RDKit Mol")
        return mol

    def biopython_to_pdbio(self, obj):
        """
        Convert a biopython object to a PDB file

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        str
            The PDB file
        """
        if not is_biopython(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        io = bio.PDBIO()
        io.set_structure(obj)
        io.save(self.__fileio__)

    def pdbio_to_biopython(self, id: str = "temp") -> "bio.Structure.Structure":
        """
        Convert the internal FileIO to a biopython structure object

        Parameters
        ----------
        id : str, optional
            The ID of the biopython object, by default "temp"

        Returns
        -------
        Structure
            The converted object
        """
        return defaults.__bioPDBParser__.get_structure(id, self.__fileio__)

    def molecule_to_pdbio(self, mol):
        """
        Convert an biobuild molecule to a PDB file

        Parameters
        ----------
        mol : object
            The molecule to convert

        Returns
        -------
        str
            The PDB file
        """
        mol.to_pdb(self.__fileio__, symmetric=False)


def is_biopython(obj):
    """
    Check if an object is a biopython object

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
        `True` if the object is a biopython object, `False` otherwise
    """
    _valids = {
        bio.Structure.Structure,
        bio.Model.Model,
        bio.Chain.Chain,
        bio.Residue.Residue,
        bio.Atom.Atom,
    }
    return any([isinstance(obj, _valid) for _valid in _valids])


def is_pybel(obj):
    """
    Check if an object is a pybel object

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
        `True` if the object is a pybel object, `False` otherwise
    """
    _valids = {
        pybel.Molecule,
        pybel.Residue,
        pybel.Atom,
    }
    return any([isinstance(obj, _valid) for _valid in _valids])


def is_rdkit(obj):
    """
    Check if an object is an rdkit object

    Parameters
    ----------
    obj : object
        The object to check

    Returns
    -------
    bool
        `True` if the object is an rdkit object, `False` otherwise
    """
    _valids = {
        Chem.rdchem.Mol,
        Chem.rdchem.Atom,
    }
    return any([isinstance(obj, _valid) for _valid in _valids])


if __name__ == "__main__":
    # MANNOSE = "support/examples/MAN.pdb"
    # _pybel = pybel.readfile("pdb", MANNOSE)
    # _pybel = next(_pybel)

    # converter = PybelBioPythonConverter()

    # _biopython = converter.pybel_molecule_to_biopython(_pybel)
    # print(_biopython)

    # mol = pybel.readstring("smi", "OCC1OC(O)C(C(C1O)O)O")
    # mol.addh()
    # mol.make3D()

    # converter = PybelBioPythonConverter()
    # _biopython = converter.pybel_molecule_to_biopython(mol)
    # print(_biopython)

    rdkitconv = RDKITBiopythonConverter()
    mol = Chem.MolFromPDBFile("/Users/noahhk/GIT/biobuild/support/examples/GAL.pdb")
    bp = rdkitconv.rdkit_to_biopython(mol)
    mo2 = rdkitconv.biopython_to_rdkit(bp)

    pass
