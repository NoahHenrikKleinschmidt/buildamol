"""
Other utilities for working with proteins and peptides.
"""

from .peptides import amino_acid_names_3letter
import buildamol.core as core


def select_protein(mol: core.Molecule, other_residues: list = None) -> core.Molecule:
    """
    Select and keep the protein part of a molecule.
    Note that this is an in-place operation that modifies the input molecule.
    It will remove all non-protein residues, keeping only those that are standard amino acids.

    Parameters
    ----------
    mol : Molecule
        The molecule to select the protein from.
    other_residues : list, optional
        A list of residue names to keep in addition to standard amino acids.

    Returns
    -------
    Molecule
        The selected protein.
    """
    acceptable_residues = list(amino_acid_names_3letter)
    if other_residues is not None:
        acceptable_residues.extend(other_residues)
    residues_to_keep = set(
        i for i in mol.get_residues() if i.name in acceptable_residues
    )
    if not residues_to_keep:
        raise ValueError("No protein residues found in the molecule.")

    residues_to_remove = set(mol.get_residues()) - residues_to_keep
    mol.remove_residues(residues_to_remove)
    return mol


def split_protein_and_others(mol: core.Molecule) -> tuple[core.Molecule, core.Molecule]:
    """
    Split a molecule into two parts: the protein part and the rest.
    This function will create two copies of the input molecule:
    one containing only the protein residues and the other containing any other residues.

    Parameters
    ----------
    mol : Molecule
        The molecule to split.

    Returns
    -------
    tuple[Molecule, Molecule]
        A tuple containing the protein part and the rest of the molecule.
    """
    protein_part = select_protein(mol.copy())
    other_part = mol.copy()
    other_part.remove_residues([i.serial_number for i in protein_part.get_residues()])
    return protein_part, other_part


def is_protein(mol: core.Molecule, allow_non_protein: bool = True) -> bool:
    """
    Check if a molecule is a protein.

    Parameters
    ----------
    mol : Molecule
        The molecule to check.
    allow_non_protein : bool, optional
        If True, allows non-protein residues in the molecule.

    Returns
    -------
    bool
        True if the molecule is a protein, False otherwise.
    """
    if not allow_non_protein:
        return all(res.name in amino_acid_names_3letter for res in mol.get_residues())
    return any(res.name in amino_acid_names_3letter for res in mol.get_residues())


if __name__ == "__main__":
    import buildamol.core as core
    from buildamol.extensions.bio.proteins import is_protein, split_protein_and_others

    test_prot = core.read_pdb("/Users/noahhk/Downloads/4g0d.pdb")
    print("Is protein:", is_protein(test_prot))
    protein_part, other_part = split_protein_and_others(test_prot)
    print("Protein part residues:", len(protein_part.residues))
    print("Other part residues:", len(other_part.residues))
