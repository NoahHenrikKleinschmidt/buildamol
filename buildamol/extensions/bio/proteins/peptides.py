import buildamol.core as core
import buildamol.resources as resources

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


if __name__ == "__main__":
    p = peptide("ACDEFGHIKLMNPQRSTVWY")
    p.show()
