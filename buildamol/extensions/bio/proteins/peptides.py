import buildamol.core as core
import buildamol.resources as resources
import buildamol.structural as structural
import numpy as np
from typing import Union

__all__ = ["peptide", "phi", "psi", "omega"]

__1to3 = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def peptide(seq: str) -> core.Molecule:
    """
    Create a peptide from a sequence

    Parameters
    ----------
    seq : str
        The sequence of the peptide in one-letter code

    Returns
    -------
    Molecule
        The peptide
    """
    resources.load_amino_acids()
    amino_acids = {
        aa: (resources.get_compound(__1to3[aa]) if aa in __1to3 else None) for aa in seq
    }
    for aa in amino_acids:
        if amino_acids[aa] is None:
            raise ValueError(f"Unknown amino acid: '{aa}'")

    mol: core.Molecule = amino_acids[seq[0]].copy()
    mol.set_linkage("LINK")
    for aa in seq[1:]:
        mol.attach(amino_acids[aa], use_patch=False)

    if mol.count_clashes():
        mol.optimize()
    return mol


def phi(
    mol: core.Molecule, res: Union[int, core.Residue] = None
) -> Union[float, np.ndarray]:
    """
    Compute the phi angle of a residue in a protein

    Parameters
    ----------
    mol : Molecule
        The protein
    res : int
        The residue number of the residue having the alpha carbon.
        If not provided, all residues are considered.

    Returns
    -------
    float or ndarray
        The phi angle(s) in degrees
    """
    if res is None:
        res = range(1, mol.count_residues() + 1)
        return np.array([phi(mol, r) for r in res])

    res = mol.get_residue(res)
    if res is None:
        raise ValueError(f"Residue {res} not found")

    _prev = mol.get_residue(res.serial_number - 1)
    if _prev is None:
        return np.nan

    # get the atoms
    N = res.get_atom("N")
    CA = res.get_atom("CA")
    C = res.get_atom("C")
    C_prev = _prev.get_atom("C")

    return structural.compute_dihedral(C_prev, N, CA, C)


def psi(
    mol: core.Molecule, res: Union[int, core.Residue] = None
) -> Union[float, np.ndarray]:
    """
    Compute the psi angle of a residue in a protein

    Parameters
    ----------
    mol : Molecule
        The protein
    res : int
        The residue number of the residue having the alpha carbon.
        If not provided, all residues are considered.

    Returns
    -------
    float or ndarray
        The psi angle(s) in degrees
    """
    if res is None:
        res = range(1, mol.count_residues() + 1)
        return np.array([psi(mol, r) for r in res])

    res = mol.get_residue(res)
    if res is None:
        raise ValueError(f"Residue {res} not found")

    _next = mol.get_residue(res.serial_number + 1)
    if _next is None:
        return np.nan

    # get the atoms
    N = res.get_atom("N")
    CA = res.get_atom("CA")
    C = res.get_atom("C")
    N_next = _next.get_atom("N")

    return structural.compute_dihedral(N, CA, C, N_next)


def omega(
    mol: core.Molecule, res: Union[int, core.Residue] = None
) -> Union[float, np.ndarray]:
    """
    Compute the omega angle of a residue in a protein

    Parameters
    ----------
    mol : Molecule
        The protein
    res : int
        The residue number of the residue having the carboxyl carbon.
        If not provided, all residues are considered.

    Returns
    -------
    float or ndarray
        The omega angle(s) in degrees
    """
    if res is None:
        res = range(1, mol.count_residues() + 1)
        return np.array([omega(mol, r) for r in res])

    res = mol.get_residue(res)
    if res is None:
        raise ValueError(f"Residue {res} not found")

    _next = mol.get_residue(res.serial_number + 1)
    if _next is None:
        return np.nan

    # get the atoms
    CA = res.get_atom("CA")
    C = res.get_atom("C")
    N_next = _next.get_atom("N")
    CA_next = _next.get_atom("CA")

    return structural.compute_dihedral(CA, C, N_next, CA_next)


if __name__ == "__main__":
    p = peptide("ACDEFGHIKLMNPQRSTVWY")
    p.show()
