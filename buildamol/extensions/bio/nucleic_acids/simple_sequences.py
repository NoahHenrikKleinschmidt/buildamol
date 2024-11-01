"""
Functions to work with simple DNA and RNA molecules
"""

import buildamol.core as core
import buildamol.resources as resources

resources.load_nucleotides()

__all__ = ["dna", "rna", "nucleic_acid", "get_5prime", "get_3prime"]

nucleotide_linkage = core.linkage(
    "C3'", "OP3", delete_in_target=["O3'", "HO3'"], id="phosphodiester"
)
"""
The phosphodiester linkage between nucleotides
"""
resources.add_linkage(nucleotide_linkage)


def get_5prime(mol: core.Molecule) -> core.Residue:
    """
    Get the 5' residue of a nucleic acid

    Parameters
    ----------
    mol : Molecule
        The nucleic acid molecule

    Returns
    -------
    Residue
        The 5' residue
    """
    hop3 = mol.get_atom("HOP3", by="id")
    if not hop3:
        raise ValueError("No 5' residue found based on HOP3 atom")
    return hop3.parent


def get_3prime(mol: core.Molecule) -> core.Residue:
    """
    Get the 3' residue of a nucleic acid

    Parameters
    ----------
    mol : Molecule
        The nucleic acid molecule

    Returns
    -------
    Residue
        The 3' residue
    """
    ho3 = mol.get_atom("HO3'", by="id")
    if not ho3:
        raise ValueError("No 3' residue found based on HO3'atom")
    return ho3.parent


def dna(sequence: str) -> core.Molecule:
    """
    Create a DNA molecule from a sequence

    Parameters
    ----------
    sequence : str
        The DNA sequence

    Returns
    -------
    Molecule
        The DNA molecule
    """
    sequence = sequence.upper()
    if not all(c in "ACGT" for c in sequence):
        raise ValueError("Invalid DNA sequence")
    mol = _construct_from_seq(sequence)
    mol.id = sequence
    return mol


def rna(sequence: str) -> core.Molecule:
    """
    Create an RNA molecule from a sequence

    Parameters
    ----------
    sequence : str
        The RNA sequence

    Returns
    -------
    Molecule
        The RNA molecule
    """
    sequence = sequence.upper()
    if not all(c in "ACGU" for c in sequence):
        raise ValueError("Invalid RNA sequence")
    mol = _construct_from_seq(sequence)
    mol.id = sequence
    return mol


def nucleic_acid(sequence: str) -> core.Molecule:
    """
    Create a generic nucleic acid molecule from a sequence (DNA or RNA)

    Parameters
    ----------
    sequence : str
        The nucleic acid sequence

    Returns
    -------
    Molecule
        The nucleic acid molecule
    """
    sequence = sequence.upper()
    if not all(c in "ACGTU" for c in sequence):
        raise ValueError("Invalid nucleic acid sequence")
    mol = _construct_from_seq(sequence)
    mol.id = sequence
    return mol


def _construct_from_seq(sequence: str) -> core.Molecule:
    mol = resources.get_compound(sequence[0])
    mol.set_linkage(nucleotide_linkage)
    for base in sequence[1:]:
        incoming = resources.get_compound(base)
        mol.attach(incoming)
    return mol
