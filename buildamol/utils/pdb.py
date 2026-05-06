"""
Auxiliary tools for PDB files.
"""

_BASE36_ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_EXTENDED_PREFIXES = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def _to_base36(value: int) -> str:
    if value == 0:
        return "0"
    encoded = []
    while value > 0:
        value, remainder = divmod(value, 36)
        encoded.append(_BASE36_ALPHABET[remainder])
    return "".join(reversed(encoded))


def _parse_extended_int(value: str) -> int:
    token = value.strip()
    if len(token) == 0:
        return 0

    sign = 1
    if token[0] == "-":
        sign = -1
        token = token[1:]
    elif token[0] == "+":
        token = token[1:]

    if len(token) == 0:
        return 0

    if token.isdigit():
        return sign * int(token)

    lead = token[0]
    if lead not in _EXTENDED_PREFIXES:
        return 0

    suffix = token[1:].upper()
    if any(char not in _BASE36_ALPHABET for char in suffix):
        return 0

    width = len(token)
    decimal_limit = 10**width - 1
    suffix_value = 0
    for char in suffix:
        suffix_value = suffix_value * 36 + _BASE36_ALPHABET.index(char)

    prefix_index = _EXTENDED_PREFIXES.index(lead)
    value = decimal_limit + 1 + prefix_index * (36 ** (width - 1)) + suffix_value
    return sign * value


def _format_extended_int(value: int, width: int) -> str:
    value = int(value)
    if value >= 0 and len(str(value)) <= width:
        return f"{value:>{width}d}"

    if value < 0:
        raise ValueError(
            f"Cannot encode negative integer {value} into an extended PDB field."
        )

    decimal_limit = 10**width - 1
    offset = value - (decimal_limit + 1)
    if offset < 0:
        return f"{value:>{width}d}"

    span = 36 ** (width - 1)
    prefix_index, suffix_value = divmod(offset, span)
    if prefix_index >= len(_EXTENDED_PREFIXES):
        raise ValueError(
            f"Cannot encode integer {value} into a width-{width} PDB field."
        )

    prefix = _EXTENDED_PREFIXES[prefix_index]
    suffix = _to_base36(suffix_value).rjust(width - 1, "0")
    return f"{prefix}{suffix}"


def _parse_float_field(value, default=0.0):
    token = value.strip()
    if len(token) == 0:
        return default
    try:
        return float(token)
    except ValueError:
        return default


def _format_pdb_text_field(value, width: int, align: str = "<", safe: bool = True) -> str:
    text = "" if value is None else str(value)
    if safe:
        text = text[:width]

    if align == ">":
        return f"{text:>{width}}"
    return f"{text:<{width}}"


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


def write_pdb(mol, filename, symmetric: bool = True, safe: bool = True):
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
    safe : bool, optional
        Whether to truncate fixed-width text fields to PDB-compatible widths,
        by default True.
    """
    with open(filename, "w") as f:
        if len(mol.models) > 1:
            for model in mol.get_models():
                f.write(f"MODEL {model.id}\n")
                mol.set_model(model)
                f.write(make_atoms_table(mol, safe=safe))
                f.write("\nENDMDL\n")
            f.write(make_connect_table(mol, symmetric))
            f.write("\nEND\n")
        else:
            f.write(make_atoms_table(mol, safe=safe))
            f.write("\n")
            f.write(make_connect_table(mol, symmetric))
            f.write("\nEND\n")


def encode_pdb(mol, symmetric: bool = True, safe: bool = True) -> str:
    """
    Encode a molecule to a PDB file without actually writing it to disk.

    Parameters
    ----------
    mol : Molecule
        The molecule to write.
    symmetric : bool, optional
        Whether to write the molecule in a symmetric way, by default True.
    safe : bool, optional
        Whether to truncate fixed-width text fields to PDB-compatible widths,
        by default True.

    Returns
    -------
    str
        The PDB file contents.
    """
    lines = []
    n_models = len(getattr(mol, "models", []))
    if n_models > 1:
        for model in mol.get_models():
            lines.append(f"MODEL {model.id}")
            mol.set_model(model)
            lines.append(make_atoms_table(mol, safe=safe))
            lines.append("ENDMDL")
        lines.append(make_connect_table(mol, symmetric))
        lines.append("END")
    else:
        lines.append(make_atoms_table(mol, safe=safe))
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

            atom_a = _parse_extended_int(tokens[0])
            for token in tokens[1:]:
                b = (atom_a, _parse_extended_int(token))
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
            _model = _parse_extended_int(line.split()[-1])
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
        "serial": _parse_extended_int(line[6:11]),
        "id": line[12:16].strip(),
        "alt_loc": line[16].strip(),
        "residue": line[17:20].strip(),
        "chain": line[21].strip(),
        "res_seq": _parse_extended_int(line[22:26]),
        "icode": line[26].strip(),
        "x": _parse_float_field(line[30:38]),
        "y": _parse_float_field(line[38:46]),
        "z": _parse_float_field(line[46:54]),
        "occ": _parse_float_field(line[54:60], default=1.0),
        "temp": _parse_float_field(line[60:66]),
        "element": line[76:78].strip(),
        "charge": _parse_charge(line),
    }
    return info


_plusminus = ("+", "-")


def _parse_charge(line):
    charge = line[78:80].strip()
    if len(charge) == 0:
        return 0
    elif charge[0] in _plusminus:
        try:
            return int(charge)
        except ValueError:
            return 0
    elif charge[-1] in _plusminus:
        try:
            return int(charge[::-1])
        except ValueError:
            return 0
    return 0



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
        line = "CONECT" + _format_extended_int(atom, 5)
        for c in connectivity[atom]:
            line += _format_extended_int(c, 5)
            if len(line) > 70:
                lines.append(line)
                line = "CONECT" + _format_extended_int(atom, 5)
        lines.append(line)
    return "\n".join(lines)


# atom_line = "{prefix}{serial}{neg_adj}{element}{id}{altloc}{residue} {chain}{res_serial}{icode}    {x}{y}{z}{occ}{temp}       {seg}{element}{charge}"
atom_line = "{prefix}{serial}{neg_adj}{id}{altloc}{residue} {chain}{res_serial}{icode}    {x}{y}{z}{occ}{temp}       {seg}{element}{charge}"


def make_atoms_table(mol, safe: bool = True):
    """
    Make a PDB atom table

    Parameters
    ----------
    mol : bam.Molecule
        The molecule to generate the table for.
    safe : bool, optional
        Whether to truncate fixed-width text fields to PDB-compatible widths,
        by default True.

    Returns
    -------
    str
        The table
    """
    lines = []
    for atom in mol.get_atoms():
        new_line = encode_atom(atom, safe=safe)
        lines.append(new_line)
    return "\n".join(lines)


def encode_atom(atom, safe: bool = True) -> str:
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

    occupancy = atom.occupancy if atom.occupancy is not None else 1.0
    bfactor = atom.bfactor if atom.bfactor is not None else 0.0

    new_line = atom_line.format(
        prefix=prefix,
        serial=_format_extended_int(atom.serial_number, 5),
        neg_adj=neg_adj,
        id=_format_pdb_text_field(atom.id.upper(), 4, safe=safe),
        altloc=_format_pdb_text_field(atom.altloc, 1, safe=safe),
        residue=_format_pdb_text_field(
            atom.get_parent().resname, 3, align=">", safe=safe
        ),
        chain=_format_pdb_text_field(
            atom.get_parent().get_parent().id or " ", 1, safe=safe
        ),
        res_serial=_format_extended_int(atom.get_parent().serial_number, 4),
        icode="",  # atom.get_parent().id[2],
        x=f"{atom.coord[0]:>8.3f}",
        y=f"{atom.coord[1]:>8.3f}",
        z=f"{atom.coord[2]:>8.3f}",
        occ=f"{occupancy:>6.2f}",
        temp=f"{bfactor:>6.2f}",
        seg=f"{'':<3}",
        element=_format_pdb_text_field(atom.element.upper(), 2, align=">", safe=safe),
        charge=f"{charge:>2}",
    )

    return new_line


if __name__ == "__main__":
    import buildamol as bam

    out = parse_connect_lines(
        "/Users/noahhk/GIT/glycosylator/__projects__/solf2/solf2_man5_glycosylated_raw.pdb"
    )
    out
    pass
