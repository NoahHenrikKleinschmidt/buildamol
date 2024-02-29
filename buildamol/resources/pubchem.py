"""
This module contains functions for interacting with PubChem.

Remotely accessing PubChem
--------------------------

The `query` function can be used to query PubChem for a compound. It returns the 2D and 3D
representations of the compound as `pubchempy.Compound` objects. The `pubchempy` package is
used to perform the query. The `query` function takes a query string and a query type as
input. The query type can be one of "name", "cid", "smiles", "inchi", "inchikey".
The returned outputs can be used to create ``Molecule`` objects using a function ``buildamol.Molecule._molecule_from_pubchem`` which is by default integrated
into the ``buildamol.Molecule.Molecule.from_pubchem`` class method.

Converting PubChem data to CIF
------------------------------
PubChem allows data downloads in the form of JSON and SDF files. The JSON file contains the descriptive data of a chemical compound while the SDF contains
the 3D conformer of the compound. The `pubchem_to_cif` function can be used to convert these files into a CIF file which can be used to load the compound
into BuildAMol using the ``PDBECompounds`` class. The function takes the paths to the JSON and SDF files as input and optionally an ID for the compound
and a path to the CIF file to write. If no ID is specified, the function will try to infer an ID from the JSON file.

.. code-block:: python

    from buildamol.resources import pubchem

    # ... download JSON and SDF files from PubChem for some compound
    json_file = "some_compound.json"
    sdf_file = "some_compound.sdf"

    # convert to CIF
    pubchem.pubchem_to_cif(json_file, sdf_file, id="SOMECOMPOUND", cif_file="some_compound.cif")
"""

import json
import os
import re
import pubchempy as pcp
from tabulate import tabulate

_grammar_dict = {
    1: "st",
    2: "nd",
    3: "rd",
}
"""
A mapping for the first three numbers to their ordinal suffixes.
"""

_sdf_coord_pattern = re.compile(
    r"^ {3,}(-?\d+\.\d*) {3,}(-?\d+\.\d*) {3,}(-?\d+\.\d*) +(\w+) +(\d)"
)
"""
Regex pattern for extracting atom coordinates from SDF file.
The pattern has 5 groups:
    1. x coordinate
    2. y coordinate
    3. z coordinate
    4. atom symbol
    5. charge
"""

_sdf_bond_pattern = re.compile(r"^ +(\d+) +(\d+) +(\d) {2}\d")
"""
Regex pattern for extracting bond information from SDF file.
The pattern has 3 groups:
    1. atom 1
    2. atom 2
    3. bond type
"""

_sdf_bond_type_map = {
    "1": "SING",
    "2": "DOUB",
    "3": "TRIP",
    "4": "DOUB",
}
"""
Mapping of SDF bond types to BuildAMol bond types.
"""

_sdf_charge_map = {
    "0": 0,
    "3": 1,
    "2": 2,
    "1": 3,
    "5": -1,
    "6": -2,
    "7": -3,
    "4": 0,  # would be radical, but we don't support that
}
"""
Mapping of SDF charge types to proper values.
"""

_cif_template_file = os.path.join(os.path.dirname(__file__), "pdbe_template.cif")
"""
Path to the CIF template file.
"""

_json_headers = ("TOCHeading", "Name", "String")
"""
Id providing headers in the JSON data structure
"""

# =========================================================================== #
#                               Public API                                   #
# =========================================================================== #


def query(_query: str, by: str = "name", idx: int = 0):
    """
    Query PubChem for a compound.

    Parameters
    ----------
    _query : str
        The query string
    by : str, optional
        The type of query to perform. One of "name", "cid", "smiles", "inchi", "inchikey".
        Defaults to "name".
    idx : int, optional
        The index of the compound to return if multiple compounds are found. Defaults to 0.
        Only one compound will be returned at a time.

    Returns
    -------
    tuple
        The 2D and 3D representations of the compound
    """
    _2d_comp = pcp.get_compounds(_query, by)
    _3d_comp = pcp.get_compounds(_query, by, record_type="3d")
    if len(_2d_comp) == 0 or len(_3d_comp) == 0:
        raise ValueError(f"Could not find molecule {_query} in PubChem")
    elif len(_2d_comp) > 1 or len(_3d_comp) > 1:
        Warning(
            f"Found multiple molecules for {_query} in PubChem, using the {idx+1}{_grammar_dict.get(idx+1, 'th')} one!"
        )
    return _2d_comp[idx], _3d_comp[idx]


def pubchem_to_cif(json_file: str, sdf_file: str, id: str = None, cif_file: str = None):
    """
    Convert a downloaded PubChem data entry to a CIF file. This is useful for preparing
    a PubChem component for use in BuildAMol - i.e. to relabel atoms as required by the
    CHARMM force field. Simply download the PubChem data for a molecule in form of the JSON
    data and the 3D conformer SDF file and pass the file names to this function to generate
    a CIF file from them which can be loaded into BuildAMol using the `PDBECompounds` class.

    Parameters
    ----------
    json_file : str
        The path to the JSON file containing the PubChem data of the molecule.
    sdf_file : str
        The path to the SDF file containing the 3D conformer of the molecule.
    id : str, optional
        The ID of the molecule. If not specified, an id will be inferred from the JSON file.
    cif_file : str, optional
        The path to the CIF file to write. If not specified, the CIF file will be written to
        the same directory as the JSON file with the same name (but adjusted suffix).
    """
    atoms, bonds = load_sdf(sdf_file)
    name, data = load_json(json_file)
    if id is None:
        id = make_id(name, data)

    # make the atom and bond tables
    atoms, atom_table = make_atom_table(id, atoms)
    _, bond_table = make_bond_table(id, bonds, atoms)

    # now make the CIF file
    cif = make_cif(
        id=id,
        name=name,
        type="Ligand",
        atoms=atom_table,
        bonds=bond_table,
        create_date=get_create_date(data),
        modify_date=get_modify_date(data),
        formal_charge=get_formal_charge(data),
        mol_weight=get_mol_weight(data),
        formula=get_formula(data),
        synonyms=get_synonyms(id, name, data)[1],
        descriptors=get_descriptors(id, data)[1],
    )

    # now write the CIF file
    if cif_file is None:
        # just in case the json file has a different suffix
        cif_file = json_file.replace(".json", "") + ".cif"
    with open(cif_file, "w") as f:
        f.write(cif)


__all__ = ["query", "pubchem_to_cif"]

# =========================================================================== #
#                               Private API                                  #
# =========================================================================== #


def load_json(json_file):
    """
    Load a PubChem JSON file and return the data as a dictionary.

    Parameters
    ----------
    json_file : str
        The path to the JSON file.

    Returns
    -------
    tuple
        The name of the molecule and the data as a dictionary.
    """
    with open(json_file, "r") as f:
        data = json.load(f)["Record"]
        name = data["RecordTitle"]
        data = data["Section"]

    # now remove the list structures
    data = {i[_json_headers[0]]: i for i in data}
    data = list_to_dict(data)

    return name, data


def load_sdf(sdf_file):
    """
    Load the SDF file and return the data a list of atom data and bond data.

    Parameters
    ----------
    sdf_file : str
        The path to the SDF file.

    Returns
    -------
    tuple
        The atom data and bond data.
    """
    atom_data = []
    bond_data = []
    with open(sdf_file, "r") as f:
        for line in f:
            if _sdf_coord_pattern.match(line):
                m = _sdf_coord_pattern.match(line)
                atom_data.append(m.groups())
            elif _sdf_bond_pattern.match(line):
                m = _sdf_bond_pattern.match(line)
                bond_data.append(m.groups())

    # Convert the data to a better format
    atom_data = [
        (
            idx + 1,  # atom index
            i[-2],  # atom element
            _sdf_charge_map[i[-1]],  # charge
            i[0],
            i[1],
            i[2],  # coordinatex xyz
        )
        for idx, i in enumerate(atom_data)
    ]
    bond_data = [(int(i[0]), int(i[1]), _sdf_bond_type_map[i[2]]) for i in bond_data]
    return atom_data, bond_data


def list_to_dict(_dict):
    """
    Flatten list-substructures into dictionaries.

    For some reason the PubChem JSON data structure contains lists of dictionaries
    where pure dictionaries would be sufficient. This function removes the lists
    and replaces them with dictionaries.
    """
    if isinstance(_dict, list):
        if len(_dict) == 1:
            if isinstance(_dict[0], dict):
                return list_to_dict(_dict[0])
            return _dict[0]
        id_header = match_header(_dict)
        if id_header is None:
            return _dict
        new = {i[id_header]: i for i in _dict}
        new = list_to_dict(new)
        return new
    elif isinstance(_dict, dict):
        for k in _dict:
            _dict[k] = list_to_dict(_dict[k])
    return _dict


def successive_get(_dict, *keys):
    """
    Get successively deeper keys from a dictionary.

    Parameters
    ----------
    _dict : dict
        The dictionary to get the keys from.
    *keys
        The keys to get from the dictionary.

    Returns
    -------
    object
        The object at the end of the key chain.
        If any of the keys is not found, None is returned.
    """
    for key in keys:
        _dict = _dict.get(key, None)
        if _dict is None:
            return None
    return _dict


def bracketize(s):
    """
    Ensure that a string is bracketized if it contains any special characters.
    """
    if any(i in s for i in " /:[](),@=+"):
        return f'"{s}"'
    return s


def match_header(_list):
    """
    Find an ID header in a list of dictionaries.
    """
    for id_header in _json_headers:
        for j in _list:
            if isinstance(j, dict) and id_header in j:
                return id_header


def make_cif(
    id: str,
    name: str,
    type: str,
    create_date: str,
    modify_date: str,
    mol_weight: str,
    formula: str,
    formal_charge: str,
    atoms: str,
    bonds: str,
    synonyms: str,
    descriptors: str,
) -> str:
    with open(_cif_template_file, "r") as f:
        cif = f.read()
    cif = cif.replace("{ID}", id)
    cif = cif.replace("{FULLNAME}", name)
    cif = cif.replace("{TYPE}", type)
    cif = cif.replace("{CREATE_DATE}", create_date)
    cif = cif.replace("{MODIFY_DATE}", modify_date)
    cif = cif.replace("{MOL_WEIGHT}", mol_weight)
    cif = cif.replace("{FORMULA}", formula)
    cif = cif.replace("{FORMAL_CHARGE}", formal_charge)
    cif = cif.replace("{ATOMS_TABLE}", atoms)
    cif = cif.replace("{BONDS_TABLE}", bonds)
    cif = cif.replace("{SYNONYMS_TABLE}", synonyms)
    cif = cif.replace("{DESCRIPTIONS_TABLE}", descriptors)
    return cif


def make_id(name, _dict):
    """
    Make an ID for the molecule from the PubChem data.

    Parameters
    ----------
    name : str
        The name of the molecule.

    _dict : dict
        The PubChem data.

    Returns
    -------
    str
        The ID.
    """
    id = get_iupac_condensed_id(_dict)
    if id is None:
        _names = get_synonyms(_dict)
        if len(_names) > 0:
            id = sorted(_names, key=len)[0]
    if id is None:
        id = name.replace(" ", "_").replace("-", "")[:3]
    else:
        id = id.replace(" ", "").replace("-", "")
    return id


def get_iupac_condensed_id(_dict):
    """
    Get the IUPAC condensed ID from the PubChem data.
    If it is available.
    """
    id = successive_get(
        _dict,
        "Biologic Description",
        "Section",
        "Biologic Line Notation",
        "Information",
        "IUPAC Condensed",
        "Value",
        "StringWithMarkup",
        "String",
    )
    return id


def get_synonyms(id, name, _dict):
    """
    Generate a list of synonyms of the molecule from the PubChem data.

    Parameters
    ----------
    id : str
        The ID of the molecule.
    name : str
        The name of the molecule.
    _dict : dict
        The PubChem data dictionary.

    Returns
    -------
    list
        A list of synonyms.
    str
        A cif-formatted table of the synonyms.
    """
    _names = _get_synonyms(_dict)
    _names[name] = name
    iupac_id = get_iupac_condensed_id(_dict)
    if iupac_id:
        _names[iupac_id] = iupac_id
    _synonyms = []
    for idx, name in enumerate(_names):
        _synonyms.append(
            (
                idx + 1,  # ordinal
                id,  # comp_id
                bracketize(name),  # name
            )
        )
    return _synonyms, tabulate(_synonyms, tablefmt="plain")


def make_atom_table(id, atoms):
    """
    Reformat the SDF atom data into a CIF-ready atom table.

    Parameters
    ----------
    id : str
        The ID of the molecule.
    atoms : list
        The SDF-parsed atom data.

    Returns
    -------
    list
        A reformatted list of atoms (required for `make_bond_table`).
    str
        A cif-formatted table of the atoms.
    """
    _atoms = []
    _atom_counts = {}
    for adx, element, charge, x, y, z in atoms:
        if element not in _atom_counts:
            _atom_counts[element] = 0
        _atom_counts[element] += 1
        a_id = f"{element}{_atom_counts[element]}"
        _atoms.append(
            (
                id,  # comp_id
                a_id,  # atom_id
                a_id,  # alt_atom_id
                element,  # type_symbol
                charge,  # charge
                x,  # x
                y,  # y
                z,  # z
                a_id,  # component_atom_id
                id,  # component_comp_id
                adx,  # ordinal
            )
        )
    return _atoms, tabulate(_atoms, tablefmt="plain")


def make_bond_table(id, bonds, atom_table):
    """
    Reformat the SDF bond data into a CIF-ready bond table.

    Parameters
    ----------
    id : str
        The ID of the molecule.
    bonds : list
        The SDF-parsed bond data.
    atom_table : list
        The atom table generated by `make_atom_table`.

    Returns
    -------
    list
        A reformatted bond list (not required by any other functions).
    str
        A cif-formatted table of the bonds.
    """
    _bonds = []
    for bdx, (a1, a2, bond_type) in enumerate(bonds):
        _bonds.append(
            (
                id,  # comp_id
                atom_table[a1 - 1][1],  # atom_id_1
                atom_table[a2 - 1][1],  # atom_id_2
                bond_type,  # value_order
                bdx + 1,  # ordinal
            )
        )
    return _bonds, tabulate(_bonds, tablefmt="plain")


def get_create_date(_dict):
    """
    Get the create date from the PubChem data.
    """
    date = successive_get(
        _dict, "Names and Identifiers", "Section", "Create Date", "Information", "Value"
    )
    if date:
        date = list(date.values())[0]
    else:
        date = "0000-00-00"
    return date


def get_modify_date(_dict):
    """
    Get the modify date from the PubChem data.
    """
    date = successive_get(
        _dict, "Names and Identifiers", "Section", "Modify Date", "Information", "Value"
    )
    if date:
        date = list(date.values())[0]
    else:
        date = get_create_date(_dict)
    return date


def get_mol_weight(_dict):
    """
    Get the molecular weight from the PubChem data.
    """
    weight = successive_get(
        _dict,
        "Chemical and Physical Properties",
        "Section",
        "Section",
        "Molecular Weight",
        "Information",
        "Value",
        "StringWithMarkup",
        "String",
    )
    if not weight:
        weight = "0.000"
    return weight


def get_formal_charge(_dict):
    """
    Get the formal charge from the PubChem data.
    """
    charge = successive_get(
        _dict,
        "Chemical and Physical Properties",
        "Section",
        "Section",
        "Formal Charge",
        "Information",
        "Value",
        "Number",
    )
    if not charge:
        charge = "0"
    return str(charge)


def get_formula(_dict):
    """
    Get the molecular formula for the compound.
    """
    formula = successive_get(
        _dict,
        "Names and Identifiers",
        "Section",
        "Molecular Formula",
        "Information",
        "Value",
        "StringWithMarkup",
        "String",
    )
    return formula


def get_descriptors(id, _dict):
    """
    Generate a list of descriptors of the molecule from the PubChem data.

    Parameters
    ----------
    id : str
        The ID of the molecule.
    _dict : dict
        The PubChem data dictionary.

    Returns
    -------
    list
        A list of descriptors.
    str
        A cif-formatted table of the descriptors.
    """
    _descriptors = _get_descriptors(_dict)
    _descriptions = []
    for key, value in _descriptors.items():
        _descriptions.append(
            (
                id,  # comp_id
                bracketize(key),  # type
                bracketize(value),  # descriptor
            )
        )
    return _descriptions, tabulate(_descriptions, tablefmt="plain")


def _get_synonyms(_dict):
    """
    Get a dictionary of names for the compound.
    """
    _names = successive_get(
        _dict,
        "Names and Identifiers",
        "Section",
        "Synonyms",
        "Section",
        "Information",
        "Value",
        "StringWithMarkup",
    )
    return _names


def _get_descriptors(_dict):
    """
    Get a dictionary of descriptors for the compound.
    """
    descriptors = successive_get(
        _dict, "Names and Identifiers", "Section", "Computed Descriptors", "Section"
    )
    _descriptors = {
        i: successive_get(j, "Information", "Value", "StringWithMarkup", "String")
        for i, j in descriptors.items()
    }
    return _descriptors


if __name__ == "__main__":
    pubchem_to_cif(
        "/Users/noahhk/Downloads/COMPOUND_CID_126961780-2.json",
        "/Users/noahhk/Downloads/Conformer3D_COMPOUND_CID_126961780-2.sdf",
        "Epi4",
        "Epi4.cif",
    )
