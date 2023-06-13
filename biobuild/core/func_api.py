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
    return Molecule.Molecule.from_pdb(filename, id=id)


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
    return Molecule.Molecule.from_cif(filename, id=id)


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
    return Molecule.Molecule.from_smiles(smiles, id=id)


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
        return Molecule.Molecule(mol)

    if not isinstance(mol, str):
        raise ValueError("input must be a string")

    # ------------------
    # mol may be a PDB id
    # ------------------
    if resources.has_compound(mol):
        return resources.get_compound(mol)

    if os.path.isfile(mol):
        if mol.endswith(".pdb"):
            return Molecule.Molecule.from_pdb(mol)
        elif mol.endswith(".cif"):
            return Molecule.Molecule.from_cif(mol)
        elif mol.endswith(".pkl"):
            return Molecule.Molecule.load(mol)

    try:
        return Molecule.Molecule.from_pubchem(mol)
    except:
        pass

    try:
        return Molecule.Molecule.from_smiles(mol)
    except:
        pass

    raise ValueError(f"Could not generate molecule from input: {mol}")


def polymerize(
    mol: "Molecule.Molecule", n: int, link=None, inplace: bool = False
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
    link: "_linkage.Linkage",
    at_residue_a,
    at_residue_b,
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


__all__ = [
    "read_pdb",
    "read_cif",
    "read_smiles",
    "molecule",
    "polymerize",
    "connect",
]

if __name__ == "__main__":
    # glycans = [
    #     ("NAG@1", "NAG@2", "14bb"),
    #     ("NAG@2", "BMA@1", "14bb"),
    #     ("BMA@1", "MAN@1", "13ab"),
    #     ("BMA@1", "MAN@2", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    # glycans = [
    #     ("NAG@1", "NAG@2", "14bb"),
    #     ("NAG@2", "BMA@3", "14bb"),
    #     ("BMA@3", "MAN@4", "13ab"),
    #     ("BMA@3", "MAN@5", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    all_carb_file = "/Users/noahhk/Downloads/charmm.carbs.36.all.txt"
    alltop = resources.CHARMMTopology.from_file(all_carb_file)
    utils.defaults.set_default_topology(alltop)

    iupac = "Man(b1-6)[Man(b1-3)]b-Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    mol = glycan("test", iupac)
    mol.show()

    # glycans = [
    #     ("NAG@A", "NAG@B", "14bb"),
    #     ("NAG@B", "BMA@C", "14bb"),
    #     ("BMA@C", "MAN@D", "13ab"),
    #     ("BMA@C", "MAN@E", "16ab"),
    # ]

    # mol = _parse_iupac_graph("test", glycans)
    # mol.show()

    pass
