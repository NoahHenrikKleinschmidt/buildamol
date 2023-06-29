"""
This module defines a JSON encoding for biobuild Molecules.
"""

import json
import numpy as np
import biobuild.utils.auxiliary as aux

DICTS = {
    "molecule": {
        "id": None,
        "names": [],
        "identifiers": [],
        "formula": None,
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
        "id": None,
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
        "atoms": [],
        "improper": False,
        "bond_angle_123": None,
        "bond_angle_234": None,
        "dihedral": None,
        "bond_length_12": None,
        "bond_length_34": None,
        "bond_length_13": None,
    },
}


def get_dict(key: str) -> dict:
    """
    Get the basic dictionary structure for a given key
    """
    return dict(DICTS[key])


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
        return super(NpEncoder, self).default(obj)


def encode_linkage(link: "Linkage"):
    """
    Encode a Linkage into a dictionary
    """
    _dict = get_dict("linkage")
    _dict["id"] = link.id
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
    mol: "Molecule", names: list = None, identifiers: list = None, formula: str = None
) -> dict:
    """
    Encode a Molecule as into a dictionary

    Parameters
    ----------
    mol : Molecule
        The molecule to encode
    names : list, optional
        A list of additional names of the molecules.
    identifiers : list, optional
        A list of additional identifiers of the molecules (e.g. InChI).
    formula : str, optional
        The molecular formula of the molecule. By default this is inferred from the molecule's atoms.

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

    bonds = {}
    for a, b in mol.get_bonds():
        bond = (a.serial_number, b.serial_number)
        if bond in bonds:
            bonds[bond] += 1
        else:
            bonds[bond] = 1

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


def encode(obj, **kwargs):
    """
    Encode an object into a dictionary

    Parameters
    ----------
    obj : object
        The object to encode
    kwargs : dict
        Additional keyword arguments to pass to the encoder

    Returns
    -------
    dict
        The encoded object
    """
    if type(obj).__name__ == "Molecule":
        return encode_molecule(obj, **kwargs)
    elif type(obj).__name__ == "Linkage":
        return encode_linkage(obj)
    elif type(obj).__name__ == "InternalCoordinates":
        return encode_ic(obj)
    else:
        raise TypeError(f"Cannot encode object of type {type(obj)} into a dictionary")


def read(filename: str) -> dict:
    """
    Read a molecule from a JSON file

    Parameters
    ----------
    filename : str
        The name of the file to read

    Returns
    -------
    dict
        The molecule's encoded dictionary form
    """
    return json.load(open(filename, "r"))


def write(
    obj,
    filename: str,
    **kwargs,
):
    """
    Write a molecule to a JSON file

    Parameters
    ----------
    obj : object
        The object to encode
    filename : str
        The name of the file to write
    kwargs : dict
        Additional keyword arguments to pass to the encoder
    """
    json.dump(
        encode(obj, **kwargs),
        open(filename, "w"),
        indent=4,
        cls=NpEncoder,
    )


def write_molecule(
    mol: "Molecule",
    filename: str,
    names: list = None,
    identifiers: list = None,
    formula: str = None,
):
    """
    Write a molecule to a JSON file

    Parameters
    ----------
    mol : Molecule
        The molecule to encode
    filename : str
        The name of the file to write
    names : list, optional
        A list of additional names of the molecules.
    identifiers : list, optional
        A list of additional identifiers of the molecules (e.g. InChI).
    formula : str, optional
        The molecular formula of the molecule. By default this is inferred from the molecule's atoms.
    """
    json.dump(
        encode_molecule(mol, names, identifiers, formula),
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


if __name__ == "__main__":
    import biobuild as bb

    man = bb.molecule("MAN")
    man.to_json("man.json")

    link = bb.get_linkage("14bb")
    write_linkage(link, "link.json")
