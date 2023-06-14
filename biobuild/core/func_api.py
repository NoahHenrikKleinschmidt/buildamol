"""
Toplevel API functions
"""
import os
from typing import Union
import Bio.PDB as bio

import biobuild.structural as structural
import biobuild.utils as utils
import biobuild.core.Molecule as Molecule
import biobuild.core.Linkage as _linkage
import biobuild.resources as resources


def read_pdb(filename: str, id: str = None) -> "Molecule.Molecule":
    """
    Read a PDB file and return a molecule.

    Parameters
    ----------
    filename : str
        The path to the PDB file
    id : str
        The id of the molecule

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_pdb(filename, id=id)


def write_pdb(mol: "Molecule.Molecule", filename: str) -> None:
    """
    Write a molecule to a PDB file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The path to the PDB file
    """
    mol.to_pdb(filename)


def read_cif(filename: str, id: str = None) -> "Molecule.Molecule":
    """
    Read a CIF file and return a molecule.

    Parameters
    ----------
    filename : str
        The path to the CIF file
    id : str
        The id of the molecule

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_cif(filename, id=id)


def write_cif(mol: "Molecule.Molecule", filename: str) -> None:
    """
    Write a molecule to a CIF file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The path to the CIF file
    """
    mol.to_cif(filename)


def read_smiles(smiles: str, id: str = None) -> "Molecule.Molecule":
    """
    Read a SMILES string and return a molecule.

    Parameters
    ----------
    smiles : str
        The SMILES string

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_smiles(smiles, id=id)


def molecule(mol: str) -> "Molecule":
    """
    Generate a molecule from an input string. This string can be a PDB id, filename, SMILES or InChI string, IUPAC name or abbreviation.
    This function will try its best to automatically generate the molecule with minimal user effort. However, using a dedicated class method is
    recommended for more efficient and predictable results.

    Parameters
    ----------
    mol : str
        The input string

    Returns
    -------
    molecule : Molecule
        The generated molecule
    """
    if isinstance(mol, bio.Structure.Structure):
        return Molecule(mol)

    if not isinstance(mol, str):
        raise ValueError("input must be a string")

    # ------------------
    # mol may be a PDB id
    # ------------------
    if resources.has_compound(mol):
        return resources.get_compound(mol)

    if os.path.isfile(mol):
        if mol.endswith(".pdb"):
            return Molecule.from_pdb(mol)
        elif mol.endswith(".cif"):
            return Molecule.from_cif(mol)
        elif mol.endswith(".pkl"):
            return Molecule.load(mol)

    try:
        return Molecule.from_pubchem(mol)
    except:
        pass

    try:
        return Molecule.from_smiles(mol)
    except:
        pass

    raise ValueError(f"Could not generate molecule from input: {mol}")


def polymerize(
    mol: "Molecule.Molecule",
    n: int,
    link: Union[str, "_linkage.Linkage"] = None,
    inplace: bool = False,
) -> "Molecule.Molecule":
    """
    Polymerize a molecule

    Parameters
    ----------
    mol : Molecule
        The molecule to polymerize
    n : int
        The number of monomers to add
    link : str or Linkage
        The linkage to use for polymerization. If None, the default linkage of the molecule is used.
    inplace : bool
        Whether to polymerize the molecule in place or return a new molecule

    Returns
    -------
    Molecule
        The polymerized molecule
    """
    if link is None and mol._linkage is None:
        raise ValueError(
            "No patch or recipe provided and no default is set on the molecule"
        )
    return mol.repeat(n, link, inplace=inplace)


def connect(
    mol_a: "Molecule.Molecule",
    mol_b: "Molecule.Molecule",
    link: Union[str, "_linkage.Linkage"],
    at_residue_a: Union[int, "bio.Residue.Residue"] = None,
    at_residue_b: Union[int, "bio.Residue.Residue"] = None,
    copy_a: bool = True,
    copy_b: bool = True,
    _topology=None,
) -> "Molecule.Molecule":
    """
    Connect two molecules together

    Parameters
    ----------
    mol_a : Molecule
        The first molecule
    mol_b : Molecule
        The second molecule
    link : Linkage or str
        The linkage to use for connection. This can be either an instance of the Linkage class or a string identifier
        of a pre-defined patch in the (currently loaded default or specified) CHARMMTopology.
    at_residue_a : int or bio.PDB.Residue
        The residue of the first molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    at_residue_b : int or bio.PDB.Residue
        The residue of the second molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    copy_a : bool
        Whether to copy the first molecule before connecting
    copy_b : bool
        Whether to copy the second molecule before connecting.
        If False, all atoms of the second molecule will be added to the first molecule.
    _topology : CHARMMTopology
        A specific topology to use in case a pre-existing patch is used as link and only the string identifier
        is supplied.

    Returns
    -------
    Molecule
        The connected molecule
    """
    if copy_b:
        mol_b = mol_b.copy()
    new = mol_a.attach(
        mol_b,
        link,
        at_residue=at_residue_a,
        other_residue=at_residue_b,
        inplace=not copy_a,
        _topology=_topology,
    )
    return new


def phosphorylate(
    mol: "Molecule.Molecule",
    at_atom: Union[int, str, bio.Atom.Atom],
    delete: Union[int, str, bio.Atom.Atom],
    inplace: bool = True,
) -> "Molecule.Molecule":
    """
    Phosphorylate a molecule at a specific atom

    Parameters
    ----------
    mol : Molecule
        The molecule to phosphorylate
    at_atom : int or str or Atom
        The atom to phosphorylate. If an integer is provided, the atom seqid must be used, starting at 1.
    delete : int or str or Atom
        The atom to delete. If an integer is provided, the atom seqid must be used, starting at 1.
    inplace : bool
        Whether to phosphorylate the molecule in place or return a new molecule

    Returns
    -------
    Molecule
        The phosphorylated molecule
    """
    phos = Molecule.from_compound("PO4")
    at_atom = mol.get_atom(at_atom)
    delete = mol.get_atom(delete)
    at_residue = at_atom.get_parent()

    l = _linkage.Linkage("phosphorylation")
    l.add_bond(at_atom.id, "P")
    l.add_delete(delete.id, "target")
    l.add_delete("O2", "source")

    if not inplace:
        mol = mol.copy()

    mol.attach(phos, l, at_residue=at_residue)
    return mol


__all__ = [
    "read_pdb",
    "write_pdb",
    "read_cif",
    "write_cif",
    "read_smiles",
    "molecule",
    "polymerize",
    "connect",
    "phosphorylate",
]
