"""
Functions for handling SMILES strings
"""

from openbabel import pybel
import biobuild.utils.convert as convert

def read_smiles(smiles: str, add_hydrogens: bool = True):
    """
    Read a SMILES string into a Bio.PDB.Structure

    Parameters
    ----------
    smiles : str
        The SMILES string to read
    add_hydrogens : bool
        Whether to add hydrogens to the structure

    Returns
    -------
    structure : Bio.PDB.Structure
        The structure
    """
    mol = pybel.readstring("smi", smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES string {smiles}")
    mol.make3D()

    if not add_hydrogens:
        mol.removeh()

    converter = convert.PybelBioPythonConverter()
    mol = converter.pybel_molecule_to_biopython(mol)

    return mol
