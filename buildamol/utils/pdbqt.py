"""
Write and read pdbqt files.
"""

import buildamol.utils.auxiliary as aux
has_meeko = aux.has_package("meeko")


def encode_pdbqt(molecule: "Molecule"):
    """
    Encode a molecule to a pdbqt string

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    molecule : Molecule
        The molecule to encode

    Returns
    -------
    str
        The pdbqt string
    """
    if not has_meeko:
        raise ImportError("PDBQT encoding requires Meeko")

    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy

    rdmol = molecule.to_rdkit()
    molprep = MoleculePreparation()
    pdbqt = molprep.prepare(rdmol)
    if len(pdbqt) == 0:
        raise ValueError("PDBQT encoding failed! Check the input molecule...")
    pdbqt = PDBQTWriterLegacy.write_string(pdbqt[0])
    return pdbqt[0]

def write_pdbqt(molecule: "Molecule", filename: str):
    """
    Write a molecule to a pdbqt file

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    molecule : Molecule
        The molecule to write
    filename : str
        The filename of the pdbqt file
    """
    pdbqt = encode_pdbqt(molecule)
    with open(filename, "w") as f:
        f.write(pdbqt)


def read_pdbqt(filename: str):
    """
    Read a pdbqt file into an array of atoms

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    filename : str
        The filename of the pdbqt file
    
    Returns
    -------
    list
        The list of atoms. Each entry is a tuple of form:
        ('idx', 'serial', 'name/element', 'resid', 'resname', 'chain', 'xyz/coord', 'partial_charges', 'atom_type')
    """
    if not has_meeko:
        raise ImportError("PDBQT reading requires Meeko")

    from meeko import PDBQTMolecule

    pdbqt_string = open(filename, "r").read()
    atoms = PDBQTMolecule(pdbqt_string=pdbqt_string).atoms()
    return atoms