import buildamol.core as core
import buildamol.resources as resources
import buildamol.structural as structural
import numpy as np
from typing import Union

__all__ = [
    "peptide",
    "phi",
    "psi",
    "omega",
    "sequence",
    "sequence_to_3letter",
    "sequence_to_1letter",
    "amino_acids",
    "amino_acid_names_3letter",
    "amino_acid_names_1letter",
]

_1to3 = {
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

_3to1 = {v: k for k, v in _1to3.items()}

# for compatibility with older versions
__1to3 = _1to3
__3to1 = _3to1

amino_acid_names_3letter = set(_3to1.keys())
"""
The 3-letter codes of standard amino acids.
This includes only the names of the standard amino acids; does not load any molecule mobjects. 
Use the `amino_acids` object to access the actual molecules.
"""

amino_acid_names_1letter = set(_1to3.keys())
"""
The 1-letter codes of standard amino acids.
This includes only the names of the standard amino acids; does not load any molecule mobjects.
Use the `amino_acids` object to access the actual molecules.
"""


class _amino_acids_generator:
    """
    The standard amino acids obtainable via attribute access.
    """

    def __getattr__(self, name: str) -> core.Molecule:
        if len(name) == 1:
            name = _1to3.get(name.upper(), None)
            if name is None:
                raise ValueError(f"Unknown amino acid: '{name}'")
        elif len(name) == 3:
            name = name.upper()
            if name not in amino_acid_names_3letter:
                raise ValueError(f"Unknown amino acid: '{name}'")
        resources.load_amino_acids()
        mol = resources.get_compound(name)
        if mol is None:
            raise ValueError(f"Unknown amino acid: '{name}'")
        return mol

    def __dir__(self):
        return sorted(amino_acid_names_3letter)

    def __repr__(self):
        return f"amino_acids({', '.join(sorted(amino_acid_names_3letter))})"


amino_acids = _amino_acids_generator()
"""
Access point for standard amino acids. This supports 1-letter, 3-letter codes, and full names.

Example
-------
>>> from buildamol.extensions.bio.proteins.peptides import amino_acids
>>> amino_acids.ALA  # Access by 3-letter code
Molecule(ALA)  
>>> amino_acids.D  # Access by 1-letter code
Molecule(ASP)
>>> amino_acids.proline  # Access by full name
Molecule(PRO)

Each amino acid is a new and unique `Molecule` object!
>>> amino_acids.arginine == amino_acids.arginine
False # Each access returns a new Molecule object
"""


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
        aa: (resources.get_compound(_1to3[aa]) if aa in _1to3 else None) for aa in seq
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


def sequence(mol: core.Molecule, unknown: Union[str, callable] = "X") -> str:
    """
    Get the 1-letter code sequence of a peptide. This also works for proteins with multiple chains. Chains are separated by a colon.

    Parameters
    ----------
    mol : Molecule
        The peptide
    unknown : str or callable, optional
        The character to use for unknown residues (default: 'X')
        This can also be set to a function that takes a the molecule and residue object
        and returns a string. Set to None to ignore unknown residues.

    Returns
    -------
    str
        The sequence of the peptide in one-letter code
    """
    if not isinstance(mol, core.Molecule):
        raise TypeError("Expected a Molecule object")
    if unknown is not None and not callable(unknown):
        _unknown = lambda mol, res: unknown
    elif unknown is None:
        _unknown = lambda mol, res: ""
    elif callable(unknown):
        _unknown = unknown
    else:
        raise TypeError(f"Unknown must be a string or a callable, got {type(unknown)}")

    total_seq = []
    for chain in mol.get_chains():
        chain_seq = [None] * len(chain.child_list)
        for r, res in enumerate(chain.child_list):
            name = _3to1.get(res.name, None)
            if name is None:
                name = _unknown(mol, res)
            chain_seq[r] = name
        total_seq.append("".join(chain_seq))

    return ":".join(total_seq)


def sequence_to_3letter(seq: str, sep=" ") -> str:
    """
    Convert a sequence in one-letter code to three-letter code.

    Parameters
    ----------
    seq : str
        The sequence in one-letter code.

    Returns
    -------
    str
        The sequence in three-letter code.
    sep : str, optional
        The separator to use between the three-letter codes (default: space)

    Exampl
    -------
    >>> sequence_to_3letter("ACDE")
    'ALA CYS ASP GLU'
    """
    return sep.join(_1to3[aa] if aa in _1to3 else aa for aa in seq)


def sequence_to_1letter(seq: str, sep=" ") -> str:
    """
    Convert a sequence in three-letter code to one-letter code.

    Parameters
    ----------
    seq : str
        The sequence in three-letter code.

    Returns
    -------
    str
        The sequence in one-letter code.
    sep : str, optional
        The separator to use between the three-letter codes (default: space)
    Example
    -------
    >>> sequence_to_1letter("ALA CYS ASP GLU")
    'ACDE'
    """
    return "".join(_3to1[aa] if aa in _3to1 else aa for aa in seq.split(sep))


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
    # p.show()

    print(sequence(p))
    print(amino_acids.ALA)
