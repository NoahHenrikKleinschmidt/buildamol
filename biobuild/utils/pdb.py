"""
Auxiliary tools for PDB files.
"""

from tabulate import tabulate


def write_pdb(mol, filename, symmetric: bool = True):
    """
    Write a molecule to a PDB file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write.
    filename : str
        The filename to write to.
    symmetric : bool, optional
        Whether to write the molecule in a symmetric way, by default True.
    """
    with open(filename, "w") as f:
        f.write(make_atoms_table(mol))
        f.write("\n")
        f.write(make_connect_table(mol, symmetric))
        f.write("\nEND\n")


def write_connect_lines(mol, filename):
    """
    Write "CONECT" lines to a PDB file.
    This is necessary since Biopython by default does not do that...

    Parameters
    ----------
    mol : bb.Molecule
        The molecule to generate the connectivity for.
    filename : str
        The filename to write to.
    """
    with open(filename, "r") as f:
        c = f.read()
    with open(filename, "w") as f:
        f.write(c.replace("END", "").rstrip())
        f.write("\n")
        f.write(make_connect_table(mol))
        f.write("\nEND\n")


def parse_connect_lines(filename):
    """
    Parse "CONECT" lines from a PDB file.
    This is necessary since Biopython by default does not do that...

    Parameters
    ----------
    filename : str
        The filename to parse.

    Returns
    -------
    bonds: list
        A list of tuples of atom serial numbers that are bonded.
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    bonds = []
    known_bonds = set()
    for line in lines:
        if line.startswith("CONECT"):
            tokens = line.split()[1:]
            atom_a = int(tokens[0])
            for token in tokens[1:]:
                b = (atom_a, int(token))
                # make sure we don't add the same bond twice
                if b[::-1] in known_bonds:
                    continue
                bonds.append(b)
                known_bonds.add(b)
    return bonds


def make_connect_table(mol, symmetric=True):
    """
    Make a "CONECT" table for a PDB file.
    This is necessary since Biopython by default does not do that...

    Parameters
    ----------
    mol : Molecule
        The molecule to generate the connectivity for.
    symmetric : bool, optional
        Whether to generate symmetric bonds (i.e. if A is bonded to B, then
        B is bonded to A as well). Default is True. And both are written to the file.

    Returns
    -------
    connect_lines : str
        The lines to add to the PDB file.
    """
    connectivity = {}
    for atom_a, atom_b in mol.bonds:
        a = atom_a.serial_number
        b = atom_b.serial_number
        if a not in connectivity:
            connectivity[a] = [b]
        else:
            connectivity[a].append(b)

        if symmetric:
            if b not in connectivity:
                connectivity[b] = [a]
            else:
                connectivity[b].append(a)

    lines = []
    for atom in connectivity:
        lines.append(("CONECT", atom, *connectivity[atom]))
    table = tabulate(lines, tablefmt="plain")
    table = table.replace("  ", "   ")
    return table


def make_atoms_table(mol):
    """
    Make a PDB atom table

    Parameters
    ----------
    mol : bb.Molecule
        The molecule to generate the table for.

    Returns
    -------
    str
        The table
    """
    lines = []
    for atom in mol.get_atoms():
        c1 = atom.coord[0]
        if len(str(c1)) > 6:
            c1 = round(c1, 2)
        c2 = atom.coord[1]
        if len(str(c2)) > 6:
            c2 = round(c2, 2)
        c3 = atom.coord[2]
        if len(str(c3)) > 6:
            c3 = round(c3, 2)
        lines.append(
            (
                "HETATM#",
                atom.serial_number,
                atom.id + "@",
                atom.get_parent().resname
                + "#"
                + atom.get_parent().get_parent().id
                + "#",
                atom.get_parent().id[1],
                "##",
                c1,
                c2,
                c3,
                atom.occupancy,
                atom.bfactor,
                "#######",
                atom.element,
            )
        )
    table = tabulate(
        lines,
        tablefmt="plain",
        floatfmt=(
            "g",
            "g",
            "g",
            "g",
            "g",
            "g",
            ".3f",
            ".3f",
            ".3f",
            ".2f",
            ".2f",
            "g",
            "g",
        ),
    )
    table = table.replace("#", " ").replace("@ ", "")
    return table


if __name__ == "__main__":
    import biobuild as bb

    glc = bb.molecule("GLC")  # ("/Users/noahhk/GIT/iupac_labeller/data/myglc3.pdb")
    b = write_pdb(glc, "ser.pdb")
    print(b)
    glc.show()
