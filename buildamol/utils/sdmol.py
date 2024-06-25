"""
Functions to deal with SD and Mol files
"""

import os

import buildamol.utils.auxiliary as aux

has_rdkit = aux.HAS_RDKIT
Chem = aux.Chem


# import buildamol.core.Molecule as Molecule


def read_mol(filename: str):
    """
    Read a mol file into a Molecule

    Parameters
    ----------
    filename : str
        The filename of the mol file

    Returns
    -------
    Chem.Mol
        An RDKit molecule
    """
    if not has_rdkit:
        raise ImportError("Molfile reading requires RDKit")

    if not os.path.exists(filename):
        raise FileNotFoundError(f"Could not find file {filename}")

    mol = Chem.MolFromMolFile(filename, removeHs=False, sanitize=False)
    return mol


def write_mol(mol: "Molecule", filename: str):
    """
    Write a Molecule to a mol file

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The filename of the mol file
    """
    if not has_rdkit:
        raise ImportError("Molfile writing requires RDKit")

    mol = mol.to_rdkit()
    Chem.MolToMolFile(mol, filename)


#     ct = encode_mol(mol)
#     ct.write(open(filename, "r"), "ctfile")


# def write_mol(mol: Molecule, filename: str):
#     """
#     Write a Molecule to a mol file

#     Parameters
#     ----------
#     mol : Molecule
#         The molecule to write
#     filename : str
#         The filename of the mol file
#     """

#     ct = encode_mol(mol)
#     ct.write(open(filename, "r"), "ctfile")


# def encode_mol(mol: Molecule):
#     """
#     Encode a Molecule as a MolFile
#     """
#     atoms = {
#         str(atom.serial_number): ctfile.Atom(
#             atom_number=str(atom.serial_number),
#             atom_symbol=atom.element,
#             x=str(atom.coord[0]),
#             y=str(atom.coord[1]),
#             z=str(atom.coord[2]),
#             mass_difference="0",
#             charge=str(atom.pqr_charge or 0),
#             atom_stereo_parity="0",
#             hydrogen_count="0",
#             stereo_care_box="0",
#             valence="0",
#             h0designator="0",
#             not_used1="0",
#             not_used2="0",
#             atom_atom_mapping_number="0",
#             inversion_retention_flag="0",
#             exact_change_flag="0",
#         )
#         for atom in mol.get_atoms()
#     }
#     bond_counts = {}
#     for a, b in mol.get_bonds():
#         _tuple = (str(a.serial_number), str(b.serial_number))
#         bond_counts.setdefault(_tuple, 0)
#         bond_counts[_tuple] += 1

#     bonds = [
#         ctfile.Bond(
#             first_atom=atoms[a],
#             second_atom=atoms[b],
#             bond_type=str(c),
#             bond_stereo="0",
#             not_used1="0",
#             bond_topology="0",
#             reacting_center_status="0",
#         )
#         for (a, b), c in bond_counts.items()
#     ]
#     atoms = [i for i in atoms.values()]
#     molfile = ctfile.Molefile()
#     molfile["Ctab"].atoms = atoms
#     molfile["Ctab"].bonds = bonds
#     molfile["Ctab"].n_atoms = len(atoms)
#     molfile["Ctab"].n_bonds = len(bonds)


__all__ = ["read_mol", "write_mol"]

if __name__ == "__main__":
    import buildamol as bam

    bam.load_amino_acids()

    out = read_mol("/Users/noahhk/GIT/biobuild/biobuild/resources/rdser.mol")
    mol = bam.Molecule.from_compound("SER")
    write_mol(mol, "testser.mol")
    pass
