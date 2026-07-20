"""
Other utilities for working with proteins and peptides.
"""

from .peptides import amino_acid_names_3letter, amino_acids
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


def get_protein_residues(mol: core.Molecule) -> list:
    """
    Get the protein residues from a molecule.

    Parameters
    ----------
    mol : Molecule
        The molecule to get the protein residues from.

    Returns
    -------
    list
        A list of protein residues.
    """
    return [i for i in mol.get_residues() if i.name.upper() in amino_acid_names_3letter]


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
        return all(
            res.name.upper() in amino_acid_names_3letter for res in mol.get_residues()
        )
    return any(
        res.name.upper() in amino_acid_names_3letter for res in mol.get_residues()
    )


# from pdbfixer
__substitutions = {
    "2AS": "ASP",
    "3AH": "HIS",
    "5HP": "GLU",
    "5OW": "LYS",
    "ACL": "ARG",
    "AGM": "ARG",
    "AIB": "ALA",
    "ALM": "ALA",
    "ALO": "THR",
    "ALY": "LYS",
    "ARM": "ARG",
    "ASA": "ASP",
    "ASB": "ASP",
    "ASK": "ASP",
    "ASL": "ASP",
    "ASQ": "ASP",
    "AYA": "ALA",
    "BCS": "CYS",
    "BHD": "ASP",
    "BMT": "THR",
    "BNN": "ALA",
    "BUC": "CYS",
    "BUG": "LEU",
    "C5C": "CYS",
    "C6C": "CYS",
    "CAS": "CYS",
    "CCS": "CYS",
    "CEA": "CYS",
    "CGU": "GLU",
    "CHG": "ALA",
    "CLE": "LEU",
    "CME": "CYS",
    "CSD": "ALA",
    "CSO": "CYS",
    "CSP": "CYS",
    "CSS": "CYS",
    "CSW": "CYS",
    "CSX": "CYS",
    "CXM": "MET",
    "CY1": "CYS",
    "CY3": "CYS",
    "CYG": "CYS",
    "CYM": "CYS",
    "CYQ": "CYS",
    "DAH": "PHE",
    "DAL": "ALA",
    "DAR": "ARG",
    "DAS": "ASP",
    "DCY": "CYS",
    "DGL": "GLU",
    "DGN": "GLN",
    "DHA": "ALA",
    "DHI": "HIS",
    "DIL": "ILE",
    "DIV": "VAL",
    "DLE": "LEU",
    "DLY": "LYS",
    "DNP": "ALA",
    "DPN": "PHE",
    "DPR": "PRO",
    "DSN": "SER",
    "DSP": "ASP",
    "DTH": "THR",
    "DTR": "TRP",
    "DTY": "TYR",
    "DVA": "VAL",
    "EFC": "CYS",
    "FLA": "ALA",
    "FME": "MET",
    "GGL": "GLU",
    "GL3": "GLY",
    "GLZ": "GLY",
    "GMA": "GLU",
    "GSC": "GLY",
    "HAC": "ALA",
    "HAR": "ARG",
    "HIC": "HIS",
    "HIP": "HIS",
    "HMR": "ARG",
    "HPQ": "PHE",
    "HTR": "TRP",
    "HYP": "PRO",
    "IAS": "ASP",
    "IIL": "ILE",
    "IYR": "TYR",
    "KCX": "LYS",
    "LLP": "LYS",
    "LLY": "LYS",
    "LTR": "TRP",
    "LYM": "LYS",
    "LYZ": "LYS",
    "MAA": "ALA",
    "MEN": "ASN",
    "MHS": "HIS",
    "MIS": "SER",
    "MK8": "LEU",
    "MLE": "LEU",
    "MPQ": "GLY",
    "MSA": "GLY",
    "MSE": "MET",
    "MVA": "VAL",
    "NEM": "HIS",
    "NEP": "HIS",
    "NLE": "LEU",
    "NLN": "LEU",
    "NLP": "LEU",
    "NMC": "GLY",
    "OAS": "SER",
    "OCS": "CYS",
    "OMT": "MET",
    "PAQ": "TYR",
    "PCA": "GLU",
    "PEC": "CYS",
    "PHI": "PHE",
    "PHL": "PHE",
    "PR3": "CYS",
    "PRR": "ALA",
    "PTR": "TYR",
    "PYX": "CYS",
    "SAC": "SER",
    "SAR": "GLY",
    "SCH": "CYS",
    "SCS": "CYS",
    "SCY": "CYS",
    "SEL": "SER",
    "SEP": "SER",
    "SET": "SER",
    "SHC": "CYS",
    "SHR": "LYS",
    "SMC": "CYS",
    "SOC": "CYS",
    "STY": "TYR",
    "SVA": "SER",
    "TIH": "ALA",
    "TPL": "TRP",
    "TPO": "THR",
    "TPQ": "ALA",
    "TRG": "LYS",
    "TRO": "TRP",
    "TYB": "TYR",
    "TYI": "TYR",
    "TYQ": "TYR",
    "TYS": "TYR",
    "TYY": "TYR",
}


def replace_non_standard_amino_acids(mol: core.Molecule) -> core.Molecule:
    """
    Replace non-standard amino acids in a molecule with their standard counterparts.
    This function will modify the input molecule in place.

    Parameters
    ----------
    mol : Molecule
        The molecule to modify.

    Returns
    -------
    Molecule
        The modified molecule with non-standard amino acids replaced.
    """
    standard_atom_names = {
        name: set(i.name for i in amino_acids[name].get_atoms())
        for name in amino_acid_names_3letter
    }
    atoms_to_drop = set()
    for res in mol.get_residues():
        if res.name in __substitutions:
            standard_name = __substitutions[res.name]
            res.name = standard_name
            for atom in res.get_atoms():
                if atom.name not in standard_atom_names[standard_name]:
                    atoms_to_drop.add(atom.serial_number)
    mol.remove_atoms(atoms_to_drop)
    return mol


if __name__ == "__main__":
    import buildamol.core as core
    from buildamol.extensions.bio.proteins import is_protein, split_protein_and_others

    test_prot = core.read_pdb("/Users/noahhk/Downloads/4g0d.pdb")
    print("Is protein:", is_protein(test_prot))
    protein_part, other_part = split_protein_and_others(test_prot)
    print("Protein part residues:", len(protein_part.residues))
    print("Other part residues:", len(other_part.residues))
