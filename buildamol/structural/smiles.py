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

    # Use ETKDGv3 with a fixed seed for deterministic, stereo-correct conformers.
    # useRandomCoords fallback handles edge cases where the standard embedding fails.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    ret = AllChem.EmbedMolecule(mol, params)
    if ret == -1:
        params.useRandomCoords = True
        AllChem.EmbedMolecule(mol, params)
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
        try:
            rdmol = Chem.RemoveHs(rdmol)
        except Exception:
            rdmol = Chem.RemoveHs(rdmol, sanitize=False)
    return Chem.MolToSmiles(rdmol, isomericSmiles=isomeric)


__all__ = [
    "read_smiles",
    "make_smiles",
]
