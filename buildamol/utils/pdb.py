"""
Auxiliary tools for PDB files.
"""

__amino_acids = set(
    (
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "CSE",  # selenocysteines
        "SEC",
    )
)


def find_models(filename):
    """
    Get the models from a PDB file.

    Parameters
    ----------
    filename : str
        The filename to read.

    Returns
    -------
    list
        A list of models.
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    models = []
    for line in lines:
        if line.startswith("MODEL"):
            models.append(line.replace("MODEL", "").strip())
    return models


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
        if len(mol.models) > 1:
            for model in mol.get_models():
                f.write(f"MODEL {model.id}\n")
                mol.set_model(model)
                f.write(make_atoms_table(mol))
                f.write("\nENDMDL\n")
            f.write(make_connect_table(mol, symmetric))
            f.write("\nEND\n")
        else:
            f.write(make_atoms_table(mol))
            f.write("\n")
            f.write(make_connect_table(mol, symmetric))
            f.write("\nEND\n")


def encode_pdb(mol, symmetric: bool = True) -> str:
    """
    Encode a molecule to a PDB file without actually writing it to disk.

    Parameters
    ----------
    mol : Molecule
        The molecule to write.
    symmetric : bool, optional
        Whether to write the molecule in a symmetric way, by default True.

    Returns
    -------
    str
        The PDB file contents.
    """
    lines = []
    if len(mol.models) > 1:
        for model in mol.get_models():
            lines.append(f"MODEL {model.id}")
            mol.set_model(model)
            lines.append(make_atoms_table(mol))
            lines.append("ENDMDL")
        lines.append(make_connect_table(mol, symmetric))
        lines.append("END")
    else:
        lines.append(make_atoms_table(mol))
        lines.append(make_connect_table(mol, symmetric))
        lines.append("END")
    return "\n".join(lines)


def decode_pdb(pdb: str):
    """
    Decode the contents of a PDB file.

    Parameters
    ----------
    pdb : str
        The PDB file contents.

    Returns
    -------
    atoms_table: dict
        The per-model atom tables
    connect_table: list
        A list of tuples of atom serial numbers that are bonded.
    """
    lines = pdb.split("\n")
    atoms_table = _parse_atom_lines(lines)
    connect_table = _parse_connect_lines(lines)
    return atoms_table, connect_table


def write_connect_lines(mol, filename):
    """
    Write "CONECT" lines to a PDB file.
    This is necessary since Biopython by default does not do that...

    Parameters
    ----------
    mol : bam.Molecule
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
        return _parse_connect_lines(f)


def _parse_connect_lines(lines):
    bonds = {}
    known_bonds = set()
    for line in lines:
        if line.startswith("CONECT"):
            # split the line into tokens of length 5
            line = line[6:]
            tokens = [
                line[i : i + 5].strip()
                for i in range(0, len(line), 5)
                if len(line[i : i + 5].strip()) > 0
            ]

            atom_a = int(tokens[0])
            for token in tokens[1:]:
                b = (atom_a, int(token))
                # make sure we don't add the same bond twice
                if b[::-1] in known_bonds:
                    continue
                bonds.setdefault(b, 0)
                bonds[b] += 1
                known_bonds.add(b)
    return [(*k, v) for k, v in bonds.items()]


def parse_atom_lines(filename, model=None):
    with open(filename, "r") as f:
        return _parse_atom_lines(f, model)


def _parse_atom_lines(lines, model=None):
    atoms = {-1: []}
    _model = -1
    _skip_lines = False
    for line in lines:
        if line.startswith("MODEL"):
            _model = int(line.split()[-1])
            if model is not None and not _model == model:
                _skip_lines = True
            else:
                _skip_lines = False
                atoms[_model] = []
            continue
        if _skip_lines:
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atoms[_model].append(_split_atom_line(line))
    return atoms


def _split_atom_line(line) -> tuple:
    info = {
        "serial": int(line[6:11].strip()),
        "id": line[12:16].strip(),
        "alt_loc": line[16].strip(),
        "residue": line[17:20].strip(),
        "chain": line[21].strip(),
        "res_seq": int(line[22:26].strip()),
        "icode": line[26].strip(),
        "x": line[30:38].strip(),
        "y": line[38:46].strip(),
        "z": line[46:54].strip(),
        "occ": float(line[54:60].strip()),
        "temp": float(line[60:66].strip()),
        "element": line[76:78].strip(),
        "charge": eval(line[78:80] + "0"),
    }
    return info


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
    for bond in mol.get_bonds():
        a = bond.atom1.serial_number
        b = bond.atom2.serial_number
        if a not in connectivity:
            connectivity[a] = [b] * bond.order
        else:
            connectivity[a].extend([b] * bond.order)

        if symmetric:
            if b not in connectivity:
                connectivity[b] = [a] * bond.order
            else:
                connectivity[b].extend([a] * bond.order)

    lines = []
    for atom in connectivity:
        line = "CONECT" + left_adjust(str(atom), 5)
        for c in connectivity[atom]:
            line += left_adjust(str(c), 5)
            if len(line) > 70:
                lines.append(line)
                line = "CONECT" + left_adjust(str(atom), 5)
        lines.append(line)
    return "\n".join(lines)


# atom_line = "{prefix}{serial}{neg_adj}{element}{id}{altloc}{residue} {chain}{res_serial}{icode}    {x}{y}{z}{occ}{temp}       {seg}{element}{charge}"
atom_line = "{prefix}{serial}{neg_adj}{id}{altloc}{residue} {chain}{res_serial}{icode}    {x}{y}{z}{occ}{temp}       {seg}{element}{charge}"


def make_atoms_table(mol):
    """
    Make a PDB atom table

    Parameters
    ----------
    mol : bam.Molecule
        The molecule to generate the table for.

    Returns
    -------
    str
        The table
    """
    lines = []
    for atom in mol.get_atoms():
        new_line = encode_atom(atom)
        lines.append(new_line)
    return "\n".join(lines)


def encode_atom(atom) -> str:
    """
    Make an ATOM line for a PDB file.
    """
    neg_adj = " "
    # if len(atom.id) > 3:
    #     neg_adj = ""
    #     altloc_len = 2
    # else:
    #     neg_adj = " "
    # altloc_len = 1
    if atom.pqr_charge is None or atom.pqr_charge == 0:
        charge = ""
    else:
        charge = str(abs(int(atom.pqr_charge))) + ("-" if atom.pqr_charge < 0 else "+")

    if atom.get_parent().resname in __amino_acids:
        prefix = "ATOM  "
    else:
        prefix = "HETATM"

    new_line = atom_line.format(
        prefix=prefix,
        serial=left_adjust(str(atom.serial_number), 5),
        neg_adj=neg_adj,
        id=right_adjust(atom.id.upper()[:4], 4),
        # id=right_adjust(
        #     atom.id.replace(atom.element.upper(), "").replace(
        #         atom.element.title(), ""
        #     ),
        #     2,
        # ),
        altloc=right_adjust(atom.altloc, 1),
        residue=left_adjust(atom.get_parent().resname, 3),
        chain=atom.get_parent().get_parent().id or " ",
        # resseq=left_adjust(" ", 3),
        res_serial=left_adjust(str(atom.get_parent().serial_number), 4),
        icode="",  # atom.get_parent().id[2],
        x=left_adjust(f"{atom.coord[0]:.3f}", 8),
        y=left_adjust(f"{atom.coord[1]:.3f}", 8),
        z=left_adjust(f"{atom.coord[2]:.3f}", 8),
        # occ=left_adjust("1.00", 6),
        occ=left_adjust(f"{atom.occupancy:.2f}", 6),
        # temp=left_adjust("0.00", 6),
        temp=left_adjust(f"{atom.bfactor:.2f}", 6),
        seg=right_adjust("", 3),
        element=right_adjust(atom.element.upper(), 2),
        charge=left_adjust(charge, 2),
    )

    return new_line


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


if __name__ == "__main__":
    import buildamol as bam

    out = parse_connect_lines(
        "/Users/noahhk/GIT/glycosylator/__projects__/solf2/solf2_man5_glycosylated_raw.pdb"
    )
    out
    pass
