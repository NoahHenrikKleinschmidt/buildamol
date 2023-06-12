"""
Auxiliary tools for PDB files.
"""

from tabulate import tabulate


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
        f.write(_molecule_connect_lines(mol))
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


def _molecule_connect_lines(mol):
    """
    Generate a list of lines for a PDB file that describe the connectivity
    of the molecule.

    Parameters
    ----------
    mol : bb.Molecule
        The molecule to generate the connectivity for.

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


if __name__ == "__main__":
    import biobuild as bb

    glc = bb.molecule("GLC")  # ("/Users/noahhk/GIT/iupac_labeller/data/myglc3.pdb")
    b = _molecule_connect_lines(glc)
    print(b)
    glc.show()
