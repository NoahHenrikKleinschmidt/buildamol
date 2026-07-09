"""
Functions to work with simple DNA and RNA molecules
"""

import buildamol.core as core
import buildamol.resources as resources

resources.load_nucleotides()

__all__ = [
    "dna",
    "rna",
    "nucleic_acid",
    "get_5prime",
    "get_3prime",
    "apply_phosphodiester_bonds",
]


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


def apply_phosphodiester_bonds(mol: core.Molecule) -> int:
    """
    Apply phosphodiester backbone bonds (O3'→P) between consecutive nucleotide residues.

    This is a fast, rule-based alternative to `infer_residue_connections` for nucleic
    acids — no distance search. Bonds that already exist are silently skipped.

    Parameters
    ----------
    mol : Molecule
        The nucleic acid (or mixed) molecule.

    Returns
    -------
    int
        Number of new bonds added.
    """
    added = 0
    for chain in mol.get_chains():
        residues = sorted(chain.get_residues(), key=lambda r: r.serial_number)
        for i in range(len(residues) - 1):
            res1, res2 = residues[i], residues[i + 1]
            if res2.serial_number - res1.serial_number > 1:
                continue  # gap in sequence numbering — not a real phosphodiester bond
            O3 = next((a for a in res1.get_atoms() if a.id == "O3'"), None)
            P = next((a for a in res2.get_atoms() if a.id == "P"), None)
            if O3 is not None and P is not None:
                mol.set_bond(O3, P)
                added += 1
    return added


def _construct_from_seq(sequence: str) -> core.Molecule:
    mol = resources.get_compound(sequence[0])
    mol.set_linkage(nucleotide_linkage)
    for base in sequence[1:]:
        incoming = resources.get_compound(base)
        mol.attach(incoming)
    return mol
