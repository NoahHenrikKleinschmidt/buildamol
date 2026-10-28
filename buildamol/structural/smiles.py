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

    for serial_number, atom in enumerate(mol.GetAtoms(), start=1):
        info = Chem.AtomPDBResidueInfo()
        info.SetSerialNumber(serial_number)
        info.SetName(f"{atom.GetSymbol()}{serial_number}".rjust(4))
        info.SetResidueName("UNL")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        atom.SetMonomerInfo(info)

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
    molecule: "Molecule",
    isomeric: bool = True,
    add_hydrogens: bool = False,
    assign_stereo: bool = False,
) -> str:
    """
    Generate a SMILES string from a molecule

    Parameters
    ----------
    molecule : Molecule
        The molecule to convert
    isomeric : bool
        Whether to include isomeric information (E/Z, @/@@)
    add_hydrogens : bool
        Whether to include hydrogens in the SMILES string
    assign_stereo : bool
        Whether to assign stereochemistry from the 3D coordinates before
        generating SMILES. When True, chiral centres are perceived from the
        3D structure and encoded as ``[C@H]``/``[C@@H]`` etc. in the output.
        Implies ``isomeric=True``. Default is False to preserve the previous
        behaviour (no stereo assignment).

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

    if assign_stereo:
        isomeric = True
        # Derive chiral tags from 3-D coordinates, then assign CIP descriptors.
        # AssignAtomChiralTagsFromStructure sets @/@@ on each atom from the
        # conformer geometry. AssignStereochemistry then propagates those tags
        # into the CIP R/S labels used by MolToSmiles.
        # AssignStereochemistryFromStructure (newer RDKit) does both in one call.
        if hasattr(Chem, "AssignStereochemistryFromStructure"):
            Chem.AssignStereochemistryFromStructure(rdmol)
        else:
            Chem.AssignAtomChiralTagsFromStructure(rdmol)
            Chem.AssignStereochemistry(rdmol, cleanIt=True, force=True)

    return Chem.MolToSmiles(rdmol, isomericSmiles=isomeric)


__all__ = [
    "read_smiles",
    "make_smiles",
]
