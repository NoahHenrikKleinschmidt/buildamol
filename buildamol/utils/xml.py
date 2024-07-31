import re


class XMLEntry:
    def __init__(self, name: str):
        self.name = name
        self.children = []
        self.attributes = {}
        self.indent = 0
        self.has_children = False
        self.has_parent = False

    @classmethod
    def from_string(cls, string: str):
        if not string.startswith("<"):
            raise ValueError("Invalid XML string")
        string = string[1:-1]

        has_children = True
        if string[-1] == "/":
            has_children = False
            string = string[:-1]

        if "=" not in string:
            name = string
            entry = cls(name)
        else:
            idx = string.find(" ")
            name = string[:idx]
            entry = cls(name)

            attrs = string[idx:]
            attrs = re.findall(r'(\w+)=("[^"]*")', attrs)
            for key, value in attrs:
                entry.attributes[key] = value.strip('"')

        entry.has_children = has_children
        return entry

    @property
    def opener(self):
        indent = "\t" * self.indent
        end = ">" if self.children else "/>"
        if self.attributes:
            return f"""{indent}<{self.name} {' '.join(f'{k}="{v}"' for k, v in self.attributes.items())}{end}"""
        return f"{indent}<{self.name}{end}"

    @property
    def closer(self):
        if not self.children:
            return ""
        indent = "\t" * self.indent
        return f"{indent}</{self.name}>"

    def adjust_indent(self):
        if not self.has_parent:
            self.indent = 0
        for child in self.children:
            child.indent = self.indent + 1
            child.adjust_indent()

    def add_child(self, child):
        child.indent = self.indent + 1
        child.has_parent = True
        self.has_children = True
        self.children.append(child)

    def add_children(self, children):
        for child in children:
            self.add_child(child)

    def add_attribute(self, key, value):
        self.attributes[key] = value

    def get_attribute(self, key, dtype=None):
        v = self.attributes.get(key, None)
        if v and dtype:
            return dtype(v)
        return v

    def get_children(self):
        return self.children

    def get_child(self, name):
        for child in self.children:
            if child.name == name:
                return child
        return None

    def to_string(self):
        if not self.children:
            return f"{self.opener}{self.closer}"
        return (
            self.opener + "\n" + "\n".join(map(str, self.children)) + "\n" + self.closer
        )

    def __iter__(self):
        for child in self.children:
            yield child

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.children[key]
        value = self.attributes.get(key, None)
        if not value:
            value = self.get_child(key)
        return value

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return f"XMLEntry('{self.name}', n_attributes={len(self.attributes)}, n_children={len(self.children)})"


def encode_molecule(mol: "Molecule", atom_attributes: list = None) -> XMLEntry:
    """
    Encode a Molecule as an XML tree

    Parameters
    ----------
    mol : Molecule
        The molecule to encode
    atom_attributes : list, optional
        A list of atom attributes to include in the XML, the following attributes are always included:
            - serial_number
            - id
            - element

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    if atom_attributes is None:
        atom_attributes = []

    root = XMLEntry("molecule")
    root.attributes["id"] = mol.id
    structure_entry = XMLEntry("structure")
    structure_entry.attributes["atoms"] = mol.count_atoms()
    root.add_child(structure_entry)

    for chain in mol.get_chains():
        chain_entry = XMLEntry("chain")
        chain_entry.attributes["id"] = chain.id
        structure_entry.add_child(chain_entry)

        for residue in chain.child_list:
            residue_entry = XMLEntry("residue")
            residue_entry.attributes["serial"] = residue.serial_number
            residue_entry.attributes["name"] = residue.resname
            chain_entry.add_child(residue_entry)

            for atom in residue.child_list:
                atom_entry = XMLEntry("atom")
                atom_entry.attributes["serial"] = atom.serial_number
                atom_entry.attributes["id"] = atom.id
                atom_entry.attributes["element"] = atom.element
                for attr in atom_attributes:
                    atom_entry.attributes[attr] = getattr(atom, attr, "none")
                residue_entry.add_child(atom_entry)

    bonds_entry = XMLEntry("connectivity")
    bonds_entry.attributes["bonds"] = mol.count_bonds()
    root.add_child(bonds_entry)
    for bond in mol.get_bonds():
        bond_entry = XMLEntry("bond")
        bond_entry.attributes["atom1"] = bond[0].serial_number
        bond_entry.attributes["atom2"] = bond[1].serial_number
        bond_entry.attributes["order"] = bond.order
        bonds_entry.add_child(bond_entry)

    coords_entry = XMLEntry("coordinates")
    coords_entry.attributes["models"] = len(mol.structure.child_list)
    root.add_child(coords_entry)
    for model in mol.get_models():
        model_entry = XMLEntry("model")
        model_entry.attributes["id"] = model.id
        coords_entry.add_child(model_entry)

        for atom in model.get_atoms():
            atom_entry = XMLEntry("atomcoord")
            atom_entry.attributes["serial"] = atom.serial_number
            atom_entry.attributes["x"] = round(atom.coord[0], 4)
            atom_entry.attributes["y"] = round(atom.coord[1], 4)
            atom_entry.attributes["z"] = round(atom.coord[2], 4)
            model_entry.add_child(atom_entry)
    return root


def encode_linkage(linkage: "Linkage") -> XMLEntry:
    """
    Encode a Linkage as an XML tree

    Parameters
    ----------
    linkage : Linkage
        The linkage to encode

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    root = XMLEntry("linkage")
    root.attributes["id"] = linkage.id
    root.attributes["descr"] = linkage.description
    root.attributes["atom1"] = linkage._stitch_ref_atoms[0]
    root.attributes["atom2"] = linkage._stitch_ref_atoms[1]
    root.attributes["deletes1"] = linkage.deletes[0]
    root.attributes["deletes2"] = linkage.deletes[1]
    ic = XMLEntry("ics")
    root.add_child(ic)
    for i in linkage.internal_coordinates:
        ic.add_child(encode_ic(i))
    return root


def encode_ic(ic: "InternalCoordinates") -> XMLEntry:
    """
    Encode an InternalCoordinates object as an XML tree

    Parameters
    ----------
    ic : InternalCoordinates
        The internal coordinates to encode

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    root = XMLEntry("ic")
    root.attributes["atom1"] = ic.atom1
    root.attributes["atom2"] = ic.atom2
    root.attributes["atom3"] = ic.atom3
    root.attributes["atom4"] = ic.atom4
    root.attributes["improper"] = ic.improper
    root.attributes["length12"] = ic.bond_length_12
    root.attributes["length13"] = ic.bond_length_13
    root.attributes["length34"] = ic.bond_length_34
    root.attributes["angle123"] = ic.bond_angle_123
    root.attributes["angle234"] = ic.bond_angle_234
    root.attributes["dihedral"] = ic.dihedral
    return root


def encode_topology(topology: "CHARMMTopology") -> XMLEntry:
    """
    Encode a CHARMMTopology object as an XML tree

    Parameters
    ----------
    topology : CHARMMTopology
        The topology to encode

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    root = XMLEntry("topology")
    root.attributes["id"] = topology.id
    patches = XMLEntry("patches")
    root.add_child(patches)
    for patch in topology.patches:
        patches.add_child(encode_linkage(patch))
    root.adjust_indent()
    return root


def encode_pdbe_compounds(compounds: "PDBECompounds") -> XMLEntry:
    """
    Encode a PDBECompounds object as an XML tree

    Parameters
    ----------
    compounds : PDBECompounds
        The compounds to encode

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    root = XMLEntry("pdbe_compounds")
    root.attributes["id"] = compounds.id
    metadata = XMLEntry("metadata")
    root.add_child(metadata)
    for compound in compounds._compounds:
        c = compounds._compounds[compound]
        c["id"] = compound
        metadata.add_child(_encode_meta_data_compound(c))

    pdbdata = XMLEntry("pdbdata")
    root.add_child(pdbdata)
    for compound in compounds._pdb:
        c = compounds._pdb[compound]
        c["id"] = compound
        pdbdata.add_child(_encode_pdb_data_compound(c))

    root.adjust_indent()
    return root


def _encode_pdb_data_compound(_dict):
    root = XMLEntry("compound")
    root.attributes["id"] = _dict["id"]
    residues = XMLEntry("residues")
    root.add_child(residues)
    for i in range(len(_dict["residues"]["serials"])):
        residue = XMLEntry("residue")
        residue.attributes["serial"] = _dict["residues"]["serials"][i]
        residue.attributes["name"] = _dict["residues"]["names"][i]
        residues.add_child(residue)
    atoms = XMLEntry("atoms")
    root.add_child(atoms)
    for i in range(len(_dict["atoms"]["serials"])):
        atom = XMLEntry("atom")
        atom.attributes["serial"] = _dict["atoms"]["serials"][i]
        atom.attributes["id"] = _dict["atoms"]["ids"][i]
        atom.attributes["element"] = _dict["atoms"]["elements"][i]
        atom.attributes["residue"] = _dict["atoms"]["residue"][i]
        atom.attributes["charge"] = _dict["atoms"]["charges"][i] or 0
        atom.attributes["x"] = _dict["atoms"]["coords"][i][0]
        atom.attributes["y"] = _dict["atoms"]["coords"][i][1]
        atom.attributes["z"] = _dict["atoms"]["coords"][i][2]
        atoms.add_child(atom)
    bonds = XMLEntry("bonds")
    root.add_child(bonds)
    for i in range(len(_dict["bonds"]["bonds"])):
        bond = XMLEntry("bond")
        bond.attributes["atom1"] = _dict["bonds"]["bonds"][i][0]
        bond.attributes["residue1"] = _dict["bonds"]["parents"][i][0]
        bond.attributes["atom2"] = _dict["bonds"]["bonds"][i][1]
        bond.attributes["residue2"] = _dict["bonds"]["parents"][i][1]
        bond.attributes["order"] = _dict["bonds"]["orders"][i]
        bonds.add_child(bond)
    return root


def _encode_meta_data_compound(_dict):
    root = XMLEntry("compound")
    root.attributes["id"] = _dict["id"]
    root.attributes["name"] = _dict["name"]
    root.attributes["formula"] = _dict["formula"]
    root.attributes["type"] = _dict["type"]
    root.attributes["oneletter"] = _dict["one_letter_code"]
    root.attributes["threeletter"] = _dict["three_letter_code"]
    names = XMLEntry("names")
    root.add_child(names)
    for name in _dict["names"]:
        _name = XMLEntry("value")
        _name.attributes["value"] = name
        names.add_child(_name)
    descriptors = XMLEntry("descriptors")
    root.add_child(descriptors)
    for desc in _dict["descriptors"]:
        _desc = XMLEntry("value")
        _desc.attributes["value"] = desc
        descriptors.add_child(_desc)
    return root


def decode_xml(xml: str) -> XMLEntry:
    """
    Decode an XML string into an XML tree

    Parameters
    ----------
    xml : str
        The XML string

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    xml = xml.strip().split("\n")
    root = XMLEntry.from_string(xml[0])
    parent_stack = [root]
    for line in xml[1:]:
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith("</"):
            parent_stack.pop()
            continue
        entry = XMLEntry.from_string(line)
        parent_stack[-1].add_child(entry)
        if entry.has_children:
            parent_stack.append(entry)
    return root


def write_xml(filename: str, entry: XMLEntry):
    """
    Write an XML tree to a file

    Parameters
    ----------
    filename : str
        The filename of the XML file
    entry : XMLEntry
        The root of the XML tree
    """
    with open(filename, "w") as f:
        f.write(str(entry))


def read_xml(filename: str) -> XMLEntry:
    """
    Read an XML file into an XML tree

    Parameters
    ----------
    filename : str
        The filename of the XML file

    Returns
    -------
    XMLEntry
        The root of the XML tree
    """
    f = open(filename)
    xml = f.read()
    f.close()
    return decode_xml(xml)


if __name__ == "__main__":
    import buildamol as bam

    bam.load_sugars()

    out = encode_topology(bam.get_default_topology())
    write_xml("test.xml", out)
    rev = decode_xml(open("test.xml").read())
    print(rev)
