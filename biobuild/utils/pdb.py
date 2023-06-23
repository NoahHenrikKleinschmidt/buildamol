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
    # table = table.replace("  ", "   ")
    return table


atom_line = "HETATM{serial}{neg_adj}{element}{id}{altloc}{residue} {chain}{res_serial}{icode}    {x}{y}{z}{occ}{temp}       {seg}{element}{charge}"


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
        if len(atom.id) > 3:
            neg_adj = ""
        else:
            neg_adj = " "
        if atom.pqr_charge is None:
            charge = ""
        else:
            charge = str(int(atom.pqr_charge)) + ("-" if atom.pqr_charge < 0 else "+")
        lines.append(
            atom_line.format(
                serial=left_adjust(str(atom.serial_number), 5),
                neg_adj=neg_adj,
                id=right_adjust(
                    atom.id.replace(atom.element.upper(), "").replace(
                        atom.element.title(), ""
                    ),
                    2,
                ),
                altloc=atom.altloc,
                residue=left_adjust(atom.get_parent().resname, 3),
                chain=atom.get_parent().get_parent().id,
                # resseq=left_adjust(" ", 3),
                res_serial=left_adjust(str(atom.get_parent().id[1]), 3),
                icode=" ",  # atom.get_parent().id[2],
                x=left_adjust(f"{atom.coord[0]:.3f}", 8),
                y=left_adjust(f"{atom.coord[1]:.3f}", 8),
                z=left_adjust(f"{atom.coord[2]:.3f}", 8),
                occ=left_adjust(f"{atom.occupancy:.2f}", 6),
                temp=left_adjust(f"{atom.bfactor:.2f}", 6),
                seg=right_adjust("", 3),
                element=left_adjust(atom.element.upper(), 2),
                charge=left_adjust(charge, 2),
            )
        )
    return "\n".join(lines)


def right_adjust(s, n):
    """
    Right adjust a string to a certain length.

    Parameters
    ----------
    s : str
        The string to adjust.
    n : int
        The length to adjust to.

    Returns
    -------
    str
        The adjusted string.
    """
    return s + " " * (n - len(s))


def left_adjust(s, n):
    """
    Left adjust a string to a certain length.

    Parameters
    ----------
    s : str
        The string to adjust.
    n : int
        The length to adjust to.

    Returns
    -------
    str
        The adjusted string.
    """
    return " " * (n - len(s)) + s


# def make_atoms_table(mol):
#     """
#     Make a PDB atom table

#     Parameters
#     ----------
#     mol : bb.Molecule
#         The molecule to generate the table for.

#     Returns
#     -------
#     str
#         The table
#     """
#     lines = []
#     min_coord = min(atom.coord.min() for atom in mol.get_atoms())
#     for atom in mol.get_atoms():
#         c1 = _adjust_coord(atom.coord[0] - min_coord)
#         c2 = _adjust_coord(atom.coord[1] - min_coord)
#         c3 = _adjust_coord(atom.coord[2] - min_coord)
#         lines.append(
#             (
#                 "HETATM#",
#                 atom.serial_number,
#                 atom.id + "@",
#                 atom.get_parent().resname + "#" + atom.get_parent().get_parent().id,
#                 atom.get_parent().id[1],
#                 "##",
#                 c1,
#                 c2,
#                 c3,
#                 atom.occupancy,
#                 atom.bfactor,
#                 "#######",
#                 atom.element,
#             )
#         )
#     table = tabulate(
#         lines,
#         tablefmt="plain",
#         floatfmt=(
#             "g",
#             "g",
#             "g",
#             "g",
#             "g",
#             "g",
#             ".3f",
#             ".3f",
#             ".3f",
#             ".2f",
#             ".2f",
#             "g",
#             "g",
#         ),
#         colalign=(
#             "right",
#             "right",
#             "left",
#             "left",
#             "right",
#             "right",
#             "right",
#             "left",
#             "left",
#             "right",
#             "right",
#             "right",
#             "right",
#         ),
#     )
#     table = table.replace("#", " ").replace("@ ", "")
#     return table


# def _adjust_coord(coord):
#     """
#     Adjust a coordinate entry to conform to the 6 character limit of PDB files.

#     Parameters
#     ----------
#     coord : float
#         The coordinate to adjust.

#     Returns
#     -------
#     str
#         The adjusted coordinate.
#     """
#     coord = str(coord)
#     if len(coord) > 6:
#         coord = coord[:7]
#         n = len(coord[coord.find(".") + 1 :]) - 1
#         coord = str(round(float(coord), n))
#     return coord


if __name__ == "__main__":
    import biobuild as bb

    glc = bb.molecule(
        "/Users/noahhk/GIT/biobuild/docs/_tutorials/large.pkl"
    )  # ("/Users/noahhk/GIT/iupac_labeller/data/myglc3.pdb")
    glc.autolabel()
    b = write_pdb(glc, "ser.pdb")
    print(b)
    # glc.show()
