class XMLEntry:
    def __init__(self, name: str):
        self.name = name
        self.children = []
        self.attributes = {}
        self.indent = 0
        self.has_children = False

    @classmethod
    def from_string(cls, string: str):
        if not string.startswith("<"):
            raise ValueError("Invalid XML string")
        string = string[1:-1]
        has_children = True
        if string[-1] == "/":
            has_children = False
            string = string[:-1]
        name, *attrs = string.split(" ")
        entry = cls(name)
        for attr in attrs:
            key, value = attr.split("=")
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

    def add_child(self, child):
        child.indent = self.indent + 1
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
        return self.get_child(key)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return f"XMLEntry('{self.name}', n_attributes={len(self.attributes)}, n_children={len(self.children)})"


def encode_xml(mol: "Molecule", atom_attributes: list = None) -> XMLEntry:
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
        if line.startswith("</"):
            parent_stack.pop()
            continue
        entry = XMLEntry.from_string(line)
        parent_stack[-1].add_child(entry)
        if entry.has_children:
            parent_stack.append(entry)
    return root


def write_xml(mol: "Molecule", filename: str, atom_attributes: list = None):
    """
    Write a Molecule to an XML file

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The filename of the XML file
    atom_attributes : list, optional
        A list of atom attributes to include in the XML, the following attributes are always included:
            - serial_number
            - id
            - element
    """
    with open(filename, "w") as f:
        f.write(str(encode_xml(mol, atom_attributes=atom_attributes)))


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
    xml = open(filename).read()
    return decode_xml(xml)


if __name__ == "__main__":
    import buildamol as bam

    mol = bam.molecule("dimethylbenzene")
    mol.stack("x", 10, 3)

    from time import time

    t = time()
    write_xml(mol, "dimethylbenzene.xml", atom_attributes=["type"])
    out = read_xml("dimethylbenzene.xml")
    print(time() - t)
