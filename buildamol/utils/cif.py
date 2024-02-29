"""
Functions to work with CIF files.
"""

from tabulate import tabulate

_bond_order_map = {
    1: "SING",
    2: "DOUB",
    3: "TRIP",
}
"""
Bond order mappings to strings
"""

_rev_bond_order_map = {
    "SING": 1,
    "DOUB": 2,
    "TRIP": 3,
}
"""
Reverse bond order mappings to strings
"""

_categories = (
    "comp_id",
    "atom_id_1",
    "atom_id_2",
    "label_atom_id_1",
    "label_atom_id_2",
    "value_order",
    "pdbx_ordinal",
)
"""
The categories of the bond table
"""


def write_bond_table(mol, filename, fmt_header: str = "_bond"):
    """
    Write a bond table to a CIF file.

    Parameters
    ----------
    mol : bam.Molecule
        The molecule to generate the connectivity for.
    filename : str
        The filename to write the CIF file to.
    fmt_header : str, optional
        The format string to use for the header of the table, by default "_bond"
    """
    with open(filename, "r") as f:
        c = f.read()
    with open(filename, "w") as f:
        f.write(c + "\n")
        f.write(make_bond_table(mol, fmt_header))


def make_bond_table(mol, fmt_header: str = "_bond"):
    """
    Generate a table of bond information for a molecule.

    Parameters
    ----------
    mol : bam.Molecule
        The molecule to generate the connectivity for.
    fmt_header : str, optional
        The format string to use for the header of the table, by default "_bond"

    Returns
    -------
    table : str
        The table of bond information.
    """
    lines = make_bond_lines(mol)
    table = tabulate(lines, tablefmt="plain")
    header = make_header(fmt_header)
    return header + table + "\n#   #\n"


def make_header(fmt_header: str = "_bond"):
    """
    Generate a header for a CIF file.

    Parameters
    ----------
    fmt_header : str, optional
        The format string to use for the header of the table, by default "_bond"

    Returns
    -------
    header : str
        The header of the table.
    """
    header = "#   #\nloop_\n"
    for category in _categories:
        header += fmt_header + "." + category + "\n"
    return header


def make_bond_lines(mol):
    """
    Generate a list of lines for a CIF file that describe the connectivity
    of the molecule.

    Parameters
    ----------
    mol : bam.Molecule
        The molecule to generate the connectivity for.

    Returns
    -------
    lines : list
        A list of lines that describe the connectivity of the molecule.
    """
    # bond_counts = {}
    # for bond in mol.bonds:
    #     if bond not in bond_counts:
    #         bond_counts[bond] = 1
    #     else:
    #         bond_counts[bond] += 1

    lines = []
    comp_id = mol.id
    for bdx, bond in enumerate(mol.get_bonds()):
        atom_a = bond[0]
        atom_b = bond[1]
        order = _bond_order_map[bond.order]
        lines.append(
            (
                comp_id,
                atom_a.serial_number,
                atom_b.serial_number,
                atom_a.id,
                atom_b.id,
                order,
                bdx + 1,
            )
        )
    return lines


def parse_bond_table(filename, fmt_header: str = "_bond"):
    """
    Parse the bond table from a BuildAMol generated CIF file.
    """
    bonds = {}
    at_bonds_table = False
    starter = fmt_header + "." + _categories[-1]  # the last header line
    with open(filename, "r") as f:
        for line in f:
            if line.startswith(starter):
                at_bonds_table = True
                continue
            if not at_bonds_table:
                continue
            if line[0] == "#":
                continue
            line = line.strip().split()
            if len(line) == 0:
                break
            a, b, order = int(line[1]), int(line[2]), line[-2]
            order = _rev_bond_order_map[order]
            bonds.setdefault((a, b), 0)
            bonds[(a, b)] += order

    return [(a, b, order) for (a, b), order in bonds.items()]


if __name__ == "__main__":
    from buildamol import molecule

    mol = molecule("GLC")
    print(make_bond_table(mol, fmt_header="_bond"))
