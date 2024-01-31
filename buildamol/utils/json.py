"""
This module defines a JSON encoding for biobuild Molecules.
"""
from copy import deepcopy
import json
import numpy as np
import buildamol.utils.auxiliary as aux
import buildamol.utils.info as info

DICTS = {
    "molecule": {
        "version": info.__version__,
        "id": None,
        "type": None,
        "names": [],
        "identifiers": [],
        "formula": None,
        "one_letter_code": None,
        "three_letter_code": None,
        "structure": {
            "atoms": {
                "serial": [],
                "id": [],
                "element": [],
                "parent": [],
                "coords_3d": [],
            },
            "bonds": {
                "serial": [],
                "order": [],
            },
            "residues": {
                "serial": [],
                "name": [],
                "parent": [],
            },
            "chains": {
                "serial": [],
                "id": [],
            },
            "model": {
                "id": None,
            },
        },
    },
    "linkage": {
        "version": info.__version__,
        "id": None,
        "description": None,
        "bond": {
            "target": None,
            "source": None,
        },
        "to_delete": {
            "target": [],
            "source": [],
        },
        "ics": [],
    },
    "ic": {
        "version": info.__version__,
        "atoms": [],
        "improper": False,
        "bond_angle_123": None,
        "bond_angle_234": None,
        "dihedral": None,
        "bond_length_12": None,
        "bond_length_34": None,
        "bond_length_13": None,
    },
    "charmm_topology": {
        "version": info.__version__,
        "id": None,
        "patches": [],
    },
    "pdbe_compounds": {
        "version": info.__version__,
        "id": None,
        "compounds": [],
    },
}


def read(filename: str) -> dict:
    """
    Read an object's dictionary from a JSON file

    Parameters
    ----------
    filename : str
        The name of the file to read

    Returns
    -------
    dict
        The object's encoded dictionary form
    """
    return json.load(open(filename, "r"))


def get_dict(key: str) -> dict:
    """
    Get the basic dictionary structure for a given key
    """
    return deepcopy(DICTS[key])


class NpEncoder(json.JSONEncoder):
    """
    An encoder for numpy data

    From: https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
    """

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, set):
            return list(obj)
        return super(NpEncoder, self).default(obj)


def encode_linkage(link: "Linkage"):
    """
    Encode a Linkage into a dictionary
    """
    _dict = get_dict("linkage")
    _dict["id"] = link.id
    _dict["description"] = link.description
    a, b = link.bonds[0]
    a = getattr(a, "id", a)
    b = getattr(b, "id", b)
    if a.startswith("2") and b.startswith("1"):
        a, b = b, a
    _dict["bond"]["target"] = a
    _dict["bond"]["source"] = b

    for i in link.deletes[0]:
        i = getattr(i, "id", i)
        _dict["to_delete"]["target"].append(i)
    for i in link.deletes[1]:
        i = getattr(i, "id", i)
        _dict["to_delete"]["source"].append(i)

    if link.has_IC:
        for ic in link.internal_coordinates:
            _dict["ics"].append(encode_ic(ic))

    return _dict


def encode_ic(ic: "InternalCoordinates"):
    """
    Encode an InternalCoordinates object into a dictionary
    """
    _dict = get_dict("ic")
    _dict["atoms"] = [getattr(i, "id", i) for i in ic.atoms]
    _dict["improper"] = ic.improper
    _dict["bond_angle_123"] = ic.bond_angle_123
    _dict["bond_angle_234"] = ic.bond_angle_234
    _dict["dihedral"] = ic.dihedral
    _dict["bond_length_12"] = ic.bond_length_12
    _dict["bond_length_34"] = ic.bond_length_34
    _dict["bond_length_13"] = ic.bond_length_13
    return _dict


def encode_molecule(
    mol: "Molecule",
    type: str = None,
    names: list = None,
    identifiers: list = None,
    formula: str = None,
    one_letter_code: str = None,
    three_letter_code: str = None,
) -> dict:
    """
    Encode a Molecule as into a dictionary

    Parameters
    ----------
    mol : Molecule
        The molecule to encode
    type : str, optional
        The type of the molecule (e.g. protein, nucleic acid, etc.).
    names : list, optional
        A list of additional names of the molecules.
    identifiers : list, optional
        A list of additional identifiers of the molecules (e.g. InChI).
    formula : str, optional
        The molecular formula of the molecule. By default this is inferred from the molecule's atoms.
    one_letter_code : str, optional
        The one-letter code of the molecule.
    three_letter_code : str, optional
        The three-letter code of the molecule.

    Returns
    -------
    dict
        The encoded molecule
    """
    _dict = get_dict("molecule")
    _dict["id"] = mol.id
    _dict["names"] = names or []
    _dict["identifiers"] = identifiers or []
    _dict["formula"] = formula or aux.make_formula(mol.get_atoms())
    _dict["one_letter_code"] = one_letter_code
    _dict["three_letter_code"] = three_letter_code
    _dict["type"] = type

    # --------------------
    # encode the structure
    # --------------------

    # ----------------------------------------------------------------------------------------
    _dict["structure"]["atoms"]["serial"] = [a.serial_number for a in mol.get_atoms()]
    _dict["structure"]["atoms"]["id"] = [a.id for a in mol.get_atoms()]
    _dict["structure"]["atoms"]["element"] = [
        a.element.title() for a in mol.get_atoms()
    ]
    _dict["structure"]["atoms"]["parent"] = [
        a.get_parent().id[1] for a in mol.get_atoms()
    ]
    _dict["structure"]["atoms"]["coords_3d"] = [
        a.get_coord().tolist() for a in mol.get_atoms()
    ]
    # ----------------------------------------------------------------------------------------

    bonds = {
        (bond[0].serial_number, bond[1].serial_number): bond.order
        for bond in mol.get_bonds()
    }
    # for a, b in mol.get_bonds():
    #     bond = (a.serial_number, b.serial_number)
    #     if bond in bonds:
    #         bonds[bond] += 1
    #     else:
    #         bonds[bond] = 1

    _dict["structure"]["bonds"]["serial"] = list(bonds.keys())
    _dict["structure"]["bonds"]["order"] = list(bonds.values())

    # ----------------------------------------------------------------------------------------

    _dict["structure"]["residues"]["serial"] = [r.id[1] for r in mol.get_residues()]
    _dict["structure"]["residues"]["name"] = [r.resname for r in mol.get_residues()]
    _dict["structure"]["residues"]["parent"] = [
        r.get_parent().id for r in mol.get_residues()
    ]

    # ----------------------------------------------------------------------------------------

    _dict["structure"]["chains"]["serial"] = [
        idx for idx, c in enumerate(mol.get_chains())
    ]
    _dict["structure"]["chains"]["id"] = [c.id for c in mol.get_chains()]

    # ----------------------------------------------------------------------------------------

    _dict["structure"]["model"]["id"] = mol._model.id

    return _dict


def encode_charmm_topology(top: "CHARMMTopology"):
    """
    Encode a CHARMMTopology into a dictionary

    Parameters
    ----------
    top : CHARMMTopology
        The topology to encode
    """
    _dict = get_dict("charmm_topology")
    _dict["id"] = top.id
    for link in top.patches:
        _dict["patches"].append(encode_linkage(link))
    return _dict


def encode_pdbe_compounds(compounds: "PDBECompounds"):
    """
    Encode a PDBECompounds object into a dictionary

    Parameters
    ----------
    compounds : PDBECompounds
        The compounds to encode
    """
    _dict = get_dict("pdbe_compounds")
    _dict["id"] = getattr(compounds, "id", None)
    for id in compounds.ids:
        comp_dict = compounds._compounds[id]
        mol = compounds.get(id)
        _dict["compounds"].append(
            encode_molecule(
                mol,
                type=comp_dict["type"],
                names=comp_dict["names"],
                identifiers=comp_dict["descriptors"],
                formula=comp_dict["formula"],
                one_letter_code=comp_dict["one_letter_code"],
                three_letter_code=comp_dict["three_letter_code"],
            )
        )
    return _dict


def write_molecule(
    mol: "Molecule",
    filename: str,
    type: str = None,
    names: list = None,
    identifiers: list = None,
    formula: str = None,
    one_letter_code: str = None,
    three_letter_code: str = None,
):
    """
    Write a molecule to a JSON file

    Parameters
    ----------
    mol : Molecule
        The molecule to encode
    filename : str
        The name of the file to write
    type : str, optional
        The type of the molecule (e.g. protein, nucleic acid, etc.).
    names : list, optional
        A list of additional names of the molecules.
    identifiers : list, optional
        A list of additional identifiers of the molecules (e.g. InChI).
    formula : str, optional
        The molecular formula of the molecule. By default this is inferred from the molecule's atoms.
    one_letter_code : str, optional
        The one letter code of the molecule.
    three_letter_code : str, optional
        The three letter code of the molecule.
    """
    json.dump(
        encode_molecule(
            mol, type, names, identifiers, formula, one_letter_code, three_letter_code
        ),
        open(filename, "w"),
        indent=4,
        cls=NpEncoder,
    )


def write_linkage(link: "Linkage", filename: str):
    """
    Write a linkage to a JSON file

    Parameters
    ----------
    link : Linkage
        The linkage to encode
    filename : str
        The name of the file to write
    """
    json.dump(
        encode_linkage(link),
        open(filename, "w"),
        indent=4,
        cls=NpEncoder,
    )


def write_charmm_topology(top, filename: str):
    """
    Write a CHARMM topology to a JSON file

    Parameters
    ----------
    top : Topology
        The topology to encode
    filename : str
        The name of the file to write
    """
    json.dump(
        encode_charmm_topology(top),
        open(filename, "w"),
        indent=4,
        cls=NpEncoder,
    )


def write_pdbe_compounds(compounds, filename: str):
    """
    Write a PDBECompounds object to a JSON file

    Parameters
    ----------
    compounds : PDBECompounds
        The PDBECompounds object to encode
    filename : str
        The name of the file to write
    """
    json.dump(
        encode_pdbe_compounds(compounds),
        open(filename, "w"),
        indent=4,
        cls=NpEncoder,
    )


if __name__ == "__main__":
    import buildamol as bam

    top = bam.get_default_compounds()
    for i in range(5, len(top.ids)):
        top.remove(top.ids[-1])

    write_pdbe_compounds(top, "comps.json")
