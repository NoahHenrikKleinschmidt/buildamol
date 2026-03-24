"""
Functions for handling SMILES strings
"""

# import buildamol.utils.convert as convert
import buildamol.utils.auxiliary as aux

Chem = aux.Chem
AllChem = aux.AllChem
aux.RDLogger.DisableLog("rdApp.*")


def read_smiles(smiles: str, add_hydrogens: bool = True):
    """
    Read a SMILES string using RDKit

    Parameters
    ----------
    smiles : str
        The SMILES string to read
    add_hydrogens : bool
        Whether to add hydrogens to the structure

    Returns
    -------
    Chem.Mol
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES string {smiles}")

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    # some molecules fail to embed, this may fix it
    AllChem.EmbedMolecule(mol, useRandomCoords=True)
    AllChem.UFFOptimizeMolecule(mol)

    return mol


def make_smiles(
    molecule: "Molecule", isomeric: bool = True, add_hydrogens: bool = False
) -> str:
    """
    Generate a SMILES string from a molecule

    Parameters
    ----------
    molecule : Molecule
        The molecule to convert
    isomeric : bool
        Whether to include isomeric information
    add_hydrogens : bool
        Whether to add hydrogens to the SMILES string
    Returns
    -------
    smiles : str
        The SMILES string
    """
    rdmol = molecule.to_rdkit()
    if not add_hydrogens:
        rdmol = Chem.RemoveHs(rdmol)
    return Chem.MolToSmiles(rdmol, isomericSmiles=isomeric)


__all__ = [
    "read_smiles",
    "make_smiles",
]
