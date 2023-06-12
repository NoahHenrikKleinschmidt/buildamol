"""
Some basic functions to interact with the core of biobuild.
"""
import os
from typing import Union
import Bio.PDB as bio

import biobuild.structural as structural
import biobuild.utils as utils
import biobuild.core.Molecule as Molecule
import biobuild.core.Linkage as _linkage
import biobuild.resources as resources



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
    molecule: Molecule, n: int, link=None, inplace: bool = False
) -> Molecule:
    """
    Polymerize a molecule

    Parameters
    ----------
    molecule : Molecule
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
    if link is None and molecule._linkage is None:
        raise ValueError(
            "No patch or recipe provided and no default is set on the molecule"
        )
    return molecule.repeat(n, link, inplace=inplace)


__all__ = [
    "molecule",
    "polymerize",
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
