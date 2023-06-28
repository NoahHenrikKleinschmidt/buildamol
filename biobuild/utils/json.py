"""
This module defines a JSON encoding for biobuild Molecules.
"""

import json
import numpy as np
import biobuild.utils.auxiliary as aux


def base_dict() -> dict:
    """
    The basic dictionary structure
    """
    _dict = {
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
    }
    return _dict


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


def encode(
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
    _dict = base_dict()
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
        The molecule to write
    filename : str
        The name of the file to write to
    names : list, optional
        A list of additional names of the molecules.
    identifiers : list, optional
        A list of additional identifiers of the molecules (e.g. InChI).
    """
    json.dump(
        encode(mol, names, identifiers), open(filename, "w"), indent=4, cls=NpEncoder
    )
