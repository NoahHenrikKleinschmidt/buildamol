"""
Functions for handling SMILES strings
"""

try:
    from openbabel import pybel

    use_openbabel = True
except ImportError:
    use_openbabel = False

try:
    from rdkit import Chem

    use_rdkit = True
except ImportError:
    use_rdkit = False

import biobuild.utils.convert as convert

if not use_openbabel and not use_rdkit:

    def read_smiles(smiles: str, add_hydrogens: bool = True):
        raise ImportError("Could not import either OpenBabel or RDKit")

elif use_rdkit:

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
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not parse SMILES string {smiles}")

        if add_hydrogens:
            mol = Chem.AddHs(mol)

        converter = convert.RDKITBiopythonConverter()
        mol = converter.rdkit_to_biopython(mol)

        return mol

elif use_openbabel:

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


__all__ = ["read_smiles"]
