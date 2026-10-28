nucleic_acid_names = [
    "DA",
    "DC",
    "DG",
    "DT",
    "A",
    "C",
    "G",
    "U",
]

from buildamol import core


def get_nucleic_acid_residues(mol: core.Molecule) -> list:
    """
    Get the nucleic acid residues from a molecule.

    Parameters
    ----------
    mol : Molecule
        The molecule to get the nucleic acid residues from.

    Returns
    -------
    list
        A list of nucleic acid residues.
    """
    return [i for i in mol.get_residues() if i.name.upper() in nucleic_acid_names]


def split_nucleic_acid_and_others(
    mol: core.Molecule,
) -> tuple[core.Molecule, core.Molecule]:
    """
    Split a molecule into two parts: the nucleic acid part and the rest.
    This function will create two copies of the input molecule:
    one containing only the nucleic acid residues and the other containing any other residues.

    Parameters
    ----------
    mol : Molecule
        The molecule to split.

    Returns
    -------
    tuple[Molecule, Molecule]
        A tuple containing the nucleic acid part and the rest of the molecule.
    """
    nucleic_acid_residues = set(get_nucleic_acid_residues(mol))
    if not nucleic_acid_residues:
        raise ValueError("No nucleic acid residues found in the molecule.")

    residues_to_remove = set(mol.get_residues()) - nucleic_acid_residues
    nucleic_acid_mol = mol.copy()
    nucleic_acid_mol.drop_residues(residues_to_remove)

    other_mol = mol.copy()
    other_mol.drop_residues(nucleic_acid_residues)

    return nucleic_acid_mol, other_mol


def select_nucleic_acid(
    mol: core.Molecule, other_residues: list = None
) -> core.Molecule:
    """
    Select the nucleic acid residues from a molecule, removing all other residues. This is an in-place operation that modifies the input molecule.

    Parameters
    ----------
    mol : Molecule
        The molecule to select the nucleic acid from.
    other_residues : list, optional
        A list of residue names to keep in addition to standard nucleotides.

    Returns
    -------
    Molecule
        The selected nucleic acid.
    """
    acceptable_residues = list(nucleic_acid_names)
    if other_residues is not None:
        acceptable_residues.extend(other_residues)
    residues_to_keep = set(
        i for i in mol.get_residues() if i.name in acceptable_residues
    )
    if not residues_to_keep:
        raise ValueError("No nucleic acid residues found in the molecule.")

    residues_to_remove = set(mol.get_residues()) - residues_to_keep
    mol.drop_residues(residues_to_remove)
    return mol


def is_nucleic_acid(mol: core.Molecule, allow_non_nucleic_acid: bool = False) -> bool:
    """
    Check if a molecule is a nucleic acid.

    Parameters
    ----------
    mol : Molecule
        The molecule to check.
    allow_non_nucleic_acid : bool, optional
        If True, the function will return True if the molecule contains any nucleic acid residues,
        even if it also contains other residues. If False, the function will return True only if all
        residues in the molecule are nucleic acids. Default is False.

    Returns
    -------
    bool
        True if the molecule is a nucleic acid, False otherwise.
    """
    if not allow_non_nucleic_acid:
        return all(res.name.upper() in nucleic_acid_names for res in mol.get_residues())
    return any(res.name.upper() in nucleic_acid_names for res in mol.get_residues())


__all__ = [
    "get_nucleic_acid_residues",
    "nucleic_acid_names",
    "select_nucleic_acid",
    "is_nucleic_acid",
    "split_nucleic_acid_and_others",
]
