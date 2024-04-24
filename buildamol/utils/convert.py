"""
Functions to convert different data formats
"""

# raise DeprecationWarning("This module is currently dropped...")

import tempfile
import os
import numpy as np

from periodictable import elements

import Bio.PDB as bio
import buildamol.utils.auxiliary as aux
import buildamol.utils.defaults as defaults

has_pybel = aux.HAS_PYBEL
has_rdkit = aux.HAS_RDKIT
has_openmm = aux.HAS_OPENMM


def mol_to_numpy_array(mol, include_bonds: bool = False):
    """
    Convert a molecule to a numpy array

    Parameters
    ----------
    mol : object
        The molecule to convert
    include_bonds : bool, optional
        Whether to include the bonds in the array
        If False, the bond array will be empty

    Returns
    -------
    tuple,
        An array of atomic numbers and atom coordinates
        and an array of bonds between atoms as serial numbers of a,b and bond order
    """
    coords = np.array(
        [
            (atom.serial_number, atom.atomic_number, atom.coords)
            for atom in mol.get_atoms()
        ]
    )
    if include_bonds:
        bonds = np.zeros((coords.shape[0], 3))
        for idx, b in enumerate(mol.get_bonds()):
            bonds[idx, :] = [b.atom1.serial_number, b.atom2.serial_number, b.order]
    else:
        bonds = np.empty((0, 1))
    return coords, bonds


class PDBIO:
    """
    The base class for intermediary PDB file-based conversions
    """

    def __init__(self):
        self.__fileio__ = tempfile.mktemp(suffix=".pdb")

    def biopython_to_pdbio(self, obj):
        """
        Store a biopython object in a PDB file
        """
        if not is_biopython(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        io = bio.PDBIO()
        io.set_structure(obj)
        io.save(self.__fileio__)

    def molecule_to_pdbio(self, obj):
        """
        Store a BuildAMol Molecule in a PDB file
        """
        obj.to_pdb(self.__fileio__, symmetric=False)

    def pdbio_to_biopython(self, id: str = "temp") -> "bio.Structure.Structure":
        """
        Convert the internal FileIO to a biopython structure object
        """
        return defaults.__bioPDBParser__.get_structure(id, self.__fileio__)

    def cleanup(self):
        """
        Remove the temporary PDB file
        """
        os.remove(self.__fileio__)


class STKBuildAMolConverter(PDBIO):
    """
    Convert STK to buildamol and vice versa
    """

    def __init__(self):
        self.__fileio__ = tempfile.mktemp(suffix=".mol")

    def buildamol_to_stk(self, obj):
        """
        Convert a BuildAMol object to an STK object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not aux.HAS_STK:
            raise ImportError("Could not import stk")
        obj.to_molfile(self.__fileio__)
        return aux.stk.BuildingBlock.init_from_file(self.__fileio__)

    def stk_to_pdbio(self, obj):
        """
        Convert an STK object to a BuildAMol object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not aux.HAS_STK:
            raise ImportError("Could not import stk")
        aux.stk.PdbWriter().write(obj, self.__fileio__)


class OpenMMBioPythonConverter(PDBIO):
    """
    Convert OpenMM data structures to Biopython
    """

    def biopython_to_openmm(self, obj):
        """
        Convert a biopython object to an OpenMM object

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
        new = self._pdbio_to_openmm()
        self.cleanup()
        return new

    def buildamol_to_openmm(self, obj):
        """
        Convert a BuildAMol Molecule to an OpenMM object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        self.molecule_to_pdbio(obj)
        new = self._pdbio_to_openmm()
        self.cleanup()
        return new

    def openmm_to_biopython(self, topology, positions):
        """
        Convert an OpenMM object to a biopython object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not has_openmm:
            raise ImportError("Could not import OpenMM")
        self._openmm_to_pdbio(topology, positions)
        new = self.pdbio_to_biopython()
        self.cleanup()
        return new

    def _pdbio_to_openmm(self):
        """
        Convert the internal FileIO to an OpenMM object

        Returns
        -------
        object
            The converted object
        """
        if not has_openmm:
            raise ImportError("Could not import OpenMM")
        pdb = aux.openmm.PDBFile(self.__fileio__)
        return pdb

    def _openmm_to_pdbio(self, topology, positions):
        """
        Convert an OpenMM object to a PDB file

        Parameters
        ----------
        topology : object
            The topology to convert
        positions : object
            The positions to convert

        Returns
        -------
        str
            The PDB file
        """
        if not has_openmm:
            raise ImportError("Could not import OpenMM")
        aux.openmm.PDBFile.writeFile(topology, positions, open(self.__fileio__, "w"))


class PybelBioPythonConverter(PDBIO):
    """
    Convert Pybel (openbabel) data structures to Biopython
    """

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
        if isinstance(obj, aux.pybel.Molecule):
            return self.pybel_molecule_to_biopython(obj)

        elif isinstance(obj, aux.pybel.Residue):
            residue = self.pybel_residue_to_biopython(obj)
            if full_structure:
                struct = self._new_biopython_structure()
                struct[0]["A"].add(residue)
                for a in residue.get_atoms():
                    residue.detach_child(a.id)
                    residue.add(a)
                return struct
            return residue

        elif isinstance(obj, aux.pybel.Atom):
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
        # # this - thereby, of course, destroying BuildAMol capacity to work
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

    def buildamol_to_pybel(self, obj):
        """
        Convert a BuildAMol object to a pybel object

        Parameters
        ----------
        obj : object
            The object to convert

        Returns
        -------
        object
            The converted object
        """
        if not has_pybel:
            raise ImportError("Could not import pybel")

        if type(obj).__name__ == "Molecule":
            obj.to_pdb(self.__fileio__)
        elif type(obj).__name__ == "Structure":
            io = bio.PDBIO()
            io.set_structure(obj.to_biopython())
            io.save(self.__fileio__)
        else:
            while not type(obj).__name__ in ["Molecule", "Structure"]:
                obj = obj.get_parent()
            return self.buildamol_to_pybel(obj)

        # read the file back in with pybel
        pybel_obj = next(aux.pybel.readfile("pdb", self.__fileio__), None)
        if not pybel_obj:
            raise ValueError("Could not convert to pybel object")
        return pybel_obj


class RDKITBiopythonConverter(PDBIO):
    """
    A class to convert between RDKit and biopython objects
    This is done simply by using a temporary PDB pseudo-file
    """

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
        if not has_rdkit:
            raise ImportError("Could not import RDKit")

        if not is_rdkit(obj):
            raise ValueError(f"Cannot convert object of type {type(obj)}")
        mol = self._rdkit_to_pdbio(obj)
        new = self.pdbio_to_biopython(mol)
        self.cleanup()
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
        new = self._pdbio_to_rdkit()
        self.cleanup()
        return new

    def _rdkit_to_pdbio(self, obj):
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
        aux.Chem.MolToPDBFile(obj, self.__fileio__)

    def _pdbio_to_rdkit(self) -> "Chem.rdchem.Mol":
        """
        Convert the internal FileIO to an RDKit Mol object

        Returns
        -------
        Mol
            The converted object
        """
        try:
            mol = aux.Chem.MolFromPDBFile(
                self.__fileio__, proximityBonding=False, removeHs=False
            )
        except:
            mol = aux.Chem.MolFromPDBFile(
                self.__fileio__, proximityBonding=False, removeHs=True
            )

        if mol is None:
            raise ValueError("Could not convert PDB file to RDKit Mol")
        return mol


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
        aux.pybel.Molecule,
        aux.pybel.Residue,
        aux.pybel.Atom,
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
        aux.Chem.rdchem.Mol,
        aux.Chem.rdchem.Atom,
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
    mol = aux.Chem.MolFromPDBFile("/Users/noahhk/GIT/biobuild/support/examples/GAL.pdb")
    bp = rdkitconv.rdkit_to_biopython(mol)
    mo2 = rdkitconv.biopython_to_rdkit(bp)

    pass
