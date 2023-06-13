"""
This module defines a Parser to extract information from the PDBE compound library.

Toplevel functions
==================

For convenience, the module provides a number of toplevel functions to parse the PDBE compound library.

.. autofunction:: read_pdbe_compounds


Reading the PDBE Compound Library
=================================

The PDBE Compound Library is a database of small molecules that are found in PDB structures. It can be downloaded as an mmCIF file. 
`biobuild` provides a parser to extract information from the mmCIF file.

To parse a PDBE Compound Library file, use the `PDBECompounds` class:

.. code-block:: python

    from biobuild.resources import PDBECompounds

    compounds = PDBECompounds.from_file("path/to/pdbe-compounds.cif")

Because parsing may be an intensive operation, the `PDBECompounds` class implements a `save` method to save the parsed data to a pickle file.
For future sessions the pre-parsed object can be directly loaded using the `load` method:

.. code-block:: python

    # save an existing PDBECompounds object
    compounds.save("path/to/pdbe-compounds.pkl")
    
    # load a pre-parsed PDBECompounds object
    compounds = PDBECompounds.load("path/to/pdbe-compounds.pkl")

    
Working with PDBECompounds
==========================

The `PDBECompounds` class provides a dictionary-like interface to the compounds in the library. It supports a number of query methods unified within the `get`. 
Compounds can be obtained from the library from their 3-letter _PDB ID_, their _name_, _chemical formula_, or _SMILES_ or _InChI_ string.


.. code-block:: python

    # get a compound by its PDB ID
    glc = compounds.get("GLC")

    # get a compound by its name
    glc = compounds.get("alpha d-glucose", by="name")

    # get a compound by its chemical formula
    glc = compounds.get("C6 H12 O6", by="formula")

    # get a compound by its SMILES string
    glc = compounds.get("C([C@@H]1[C@H]([C@@H]([C@H](O[C@@H]1O)O)O)O)O", by="smiles")

    # get a compound by its InChI string (note the "InChI=" prefix)
    # THIS ALSO USES `by="smiles"`!
    glc = compounds.get("InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1", by="smiles")

    
By default the `get` method will create a `Molecule` object from the compound. Other output formats such as `biopython.Structure` or simple `dict` are also supported and can be specified using the `return_type` argument.
If multiple compounds match a query, they are returned as a list (unless `return_type` is a dictionary, in which they are kept as a dictionary).    


.. code-block:: python

    # get a compound as a `biopython.Structure` object
    glc = compounds.get("GLC", return_type="structure")

    # get a compound as a `dict`
    glc = compounds.get("GLC", return_type="dict")

    

Setting default PDBECompounds
=============================

_biobuild_  loads a default PDBECompounds object for convenience. The default instance can be accessed using the `get_default_compounds` function. A custom instance can be set as the default using the `set_default_compounds` function.

.. code-block:: python

    import biobuild as bb

    # get the default PDBECompounds instance
    compounds = bb.get_default_compounds()

    # ... do something with the compounds

    # set a custom PDBECompounds instance as the default
    bb.set_default_compounds(compounds)

.. warning::

    The `set_default_compounds` has an additional keyword argument `overwrite` which is set to `False` by default. If set to `True` the default instance
    is permanently overwritten by the new one and will be used automatically by all future sessions. This is not recommended as it may lead to unexpected behaviour.
    
"""

import os
import warnings
import numpy as np
import pdbecif.mmcif_tools as mmcif_tools
import pickle

import Bio.PDB as bio
from Bio import SVDSuperimposer

import biobuild.utils.defaults as defaults
import biobuild.utils.auxiliary as aux
import biobuild.core.Molecule as Molecule

__to_float__ = {
    "_chem_comp": set(
        (
            "formula_weight",
            "pdbx_formal_charge",
        )
    ),
}

__to_delete__ = {
    "_chem_comp": set(
        (
            "pdbx_type",
            "mon_nstd_parent_comp_id",
            "pdbx_synonyms",
            "pdbx_initial_date",
            "pdbx_modified_date",
            "pdbx_release_status",
            "pdbx_ambiguous_flag",
            "pdbx_replaced_by",
            "pdbx_replaces",
            "pdbx_processing_site",
            "pdbx_model_coordinates_details",
            "pdbx_model_coordinates_db_code",
            "pdbx_ideal_coordinates_details",
            "pdbx_subcomponent_list",
            "pdbx_model_coordinates_missing_flag",
            "pdbx_ideal_coordinates_missing_flag",
        )
    ),
}

# a list of categories to ignore while parsing
__to_ignore__ = [
    "_pdbx_chem_comp_audit",
    "_pdbx_chem_comp_feature",
]

__needs_to_have__ = set(
    (
        "_chem_comp",
        "_pdbx_chem_comp_descriptor",
        "_pdbx_chem_comp_identifier",
        "_chem_comp_atom",
        "_chem_comp_bond",
    )
)

__search_by__ = ("id", "name", "formula", "smiles")


_bond_order_rev_map = {
    "SING": 1,
    "DOUB": 2,
    "TRIP": 3,
}
"""
Reverse map of bond order names to numbers.
"""


def read_compounds(filename: str, set_default: bool = True) -> "PDBECompounds":
    """
    Reads a PDBECompounds object from a CIF file.

    Parameters
    ----------
    filename : str
        The path to the CIF file to read.
    set_default : bool, optional
        Whether to set the read PDBECompounds object as the default.
        Defaults to True.

    Returns
    -------
    PDBECompounds
        The PDBECompounds object parsed from the CIF file.
    """
    compounds = PDBECompounds.from_file(filename)
    if set_default:
        set_default_compounds(compounds)
    return compounds


def load_compounds(filename: str, set_default: bool = True) -> "PDBECompounds":
    """
    Loads a PDBECompounds object from a pickle file.

    Parameters
    ----------
    filename : str
        The path to the pickle file to load.
    set_default : bool, optional
        Whether to set the loaded PDBECompounds object as the default.
        Defaults to True.

    Returns
    -------
    PDBECompounds
        The PDBECompounds object loaded from the pickle file.
    """
    compounds = PDBECompounds.load(filename)
    if set_default:
        set_default_compounds(compounds)
    return compounds


def save_compounds(filename: str, compounds: "PDBECompounds" = None):
    """
    Saves a PDBECompounds object to a file.

    Parameters
    ----------
    filename : str
        The path to the pickle file to save to.

    compounds : PDBECompounds, optional
        The PDBECompounds object to save. If not provided, the default
        PDBECompounds object is used.
    """
    if compounds is None:
        compounds = get_default_compounds()
    compounds.save(filename)


def has_compound(compound: str, search_by: str = None) -> bool:
    """
    Checks if a given compound is available in the currently loaded default
    PDBECompounds instance.

    Parameters
    ----------
    compound : str
        The compound to check.

    search_by : str, optional
        The type of search to perform. One of `id`, `name`, `formula`, `smiles` (includes `inchi`).
        By default, all search types are used.
    """
    if search_by is not None:
        return get_default_compounds().has_residue(compound, by=search_by)

    for by in __search_by__:
        if get_default_compounds().has_residue(compound, by=by):
            return True
    return False


def get_compound(
    compound: str, search_by: str = None, return_type: str = "molecule"
) -> Molecule:
    """
    Get a compound from the currently loaded default PDBECompounds instance.

    Parameters
    ----------
    compound : str
        The compound to get.

    search_by : str, optional
        The type of search to perform. One of `id`, `name`, `formula`, `smiles` (includes `inchi`).
        By default, all search types are used.

    return_type : str, optional
        The type of object to return. One of `molecule` (a biobuild Molecule), `structure` (a biopython structure only, does not include connectivity), `dict` (the compound data and coordinate dictionaries).
        Defaults to `molecule`.

    Returns
    -------
    Molecule
        The molecule object.
    """
    if not has_compound(compound, search_by=search_by):
        raise ValueError(
            f"Compound {compound} not found in the default PDBECompounds instance."
        )
    comps = get_default_compounds()
    for by in __search_by__:
        if comps.has_residue(compound, by=by):
            comp = comps.get(compound, by=by, return_type=return_type)
            if return_type == "dict":
                comp = comp, comps._pdb[comp["id"]]
            return comp


def get_default_compounds() -> "PDBECompounds":
    """
    Get the currently loaded default PDBECompounds instance.

    Returns
    -------
    PDBECompounds
        The currently loaded default PDBECompounds instance.
    """
    if "PDBECompounds" not in defaults.__default_instances__:
        defaults.__default_instances__["PDBECompounds"] = PDBECompounds.load(
            defaults.DEFAULT_PDBE_COMPOUNDS_FILE
        )
    return defaults.__default_instances__["PDBECompounds"]


def set_default_compounds(obj, overwrite: bool = False):
    """
    Set a PDBECompounds object as the new default.

    Parameters
    ----------
    obj : PDBECompounds
        The PDBECompounds object to set as the new default
    overwrite : bool
        If set to `True`, the new object will be permanently saved as
        the default. Otherwise, the new object will only be used for
        the current session.
    """
    if not obj.__class__.__name__ == "PDBECompounds":
        raise TypeError("The object must be a PDBECompounds instance.")
    if overwrite:
        current = defaults.__default_instances__.get("PDBECompounds", None)
        if current:
            current.save(defaults.DEFAULT_PDBE_COMPOUNDS_FILE + ".bak")
    defaults.__default_instances__["PDBECompounds"] = obj
    if overwrite:
        obj.save(defaults.DEFAULT_PDBE_COMPOUNDS_FILE)


def restore_default_compounds():
    """
    Restore the default PDBECompounds object from the backup file
    """
    defaults.__default_instances__["PDBECompounds"] = pickle.load(
        open(defaults.DEFAULT_PDBE_COMPOUNDS_FILE + ".bak", "rb")
    )
    os.remove(defaults.DEFAULT_PDBE_COMPOUNDS_FILE + ".bak")


class PDBECompounds:
    """
    This class wraps a dictionary of PDBE compounds.
    It facilitates easy data access and searchability for compounds.

    Parameters
    ----------
    compounds : dict
        A dictionary of PDBE compounds.
    """

    def __init__(self, compounds: dict) -> None:
        self._compounds = {k: None for k in compounds.keys()}
        self._pdb = dict(self._compounds)
        self._setup_dictionaries(compounds)

        self._filename = None
        self._orig = compounds

    @classmethod
    def from_file(cls, filename: str) -> "PDBECompounds":
        """
        Create a PDBECompounds object from a `cif` file.

        Parameters
        ----------
        filename : str
            The path to the file.

        Returns
        -------
        PDBECompounds
            The PDBECompounds object.
        """
        _dict = mmcif_tools.MMCIF2Dict()
        _dict = _dict.parse(filename, ignoreCategories=__to_ignore__)

        # Make sure all required fields are present
        _keys = set(_dict.keys())
        for key in _keys:
            v = _dict[key]
            if not all(i in v.keys() for i in __needs_to_have__):
                del _dict[key]
                warnings.warn(
                    f"Compound '{key}' does not have all required fields. It will be ignored."
                )

        new = cls(_dict)
        new._filename = filename
        return new

    @classmethod
    def load(self, filename: str) -> "PDBECompounds":
        """
        Load a PDBECompounds object from a `pickle` file.

        Parameters
        ----------
        filename : str
            The path to the file.

        Returns
        -------
        PDBECompounds
            The PDBECompounds object.
        """
        with open(filename, "rb") as f:
            return pickle.load(f)

    @property
    def ids(self) -> list:
        """
        Get a list of all compound ids.

        Returns
        -------
        list
            A list of all compound ids.
        """
        return list(self._compounds.keys())

    @property
    def formulas(self) -> list:
        """
        Get a list of all compound formulas.

        Returns
        -------
        list
            A list of all compound formulas.
        """
        return [i["formula"] for i in self._compounds.values()]

    def save(self, filename: str = None) -> None:
        """
        Save the PDBECompounds object to a `pickle` file.

        Parameters
        ----------
        filename : str
            The path to the file. By default the same file used to load the object is used.
            If no file is available, a `ValueError` is raised.
        """
        if filename is None:
            if self._filename is None:
                raise ValueError("No filename specified.")
            filename = self._filename

        if not filename.endswith(".pkl"):
            filename = aux.change_suffix(filename, ".pkl")

        self._filename = filename

        with open(filename, "wb") as f:
            pickle.dump(self, f)

    def get(
        self,
        query: str,
        by: str = "id",
        return_type: str = "molecule",
    ):
        """
        Get a compound that matches the given criteria.

        Parameters
        ----------
        query : str
            The query to search for.

        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        return_type : str, optional
            The type of object to return, by default "molecule". Possible values are:
            - "molecule": return a biobuild `Molecule` object
            - "dict": return a dictionary of the compound data
            - "structure": return a biopython `Structure` object
            - "residue": return a biopython `Residue` object

        Returns
        -------
        object
            The object that matches the given criteria.
        """
        _dict = self._get(query, by)
        if len(_dict.keys()) == 0:
            raise ValueError(f"No compound found for query '{query}'.")
        elif len(_dict.keys()) > 1:
            return [self.get(i, "id", return_type) for i in _dict.keys()]

        _dict = list(_dict.values())[0]
        if return_type == "molecule":
            return self._molecule(_dict)
        elif return_type == "dict":
            return _dict
        elif return_type == "structure":
            return self._make_structure(_dict)
        elif return_type == "residue":
            res = self._make_residue(_dict)
            self._fill_residue(res, _dict)
            return res
        else:
            raise ValueError(f"Invalid return_type '{return_type}'.")

    def remove(self, id: str) -> None:
        """
        Remove a compound from the dictionary.

        Parameters
        ----------
        id : str
            The id of the compound to remove.
        """
        self._compounds.pop(id, None)
        self._pdb.pop(id, None)

    def has_residue(self, query: str, by: str = "id") -> bool:
        """
        Check if a compound has a residue definition.

        Parameters
        ----------
        query : str
            The query to search for.

        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        Returns
        -------
        bool
            True if the compound has a residue definition, False otherwise.
        """
        _dict = self._get(query, by)
        return len(_dict.keys()) > 0

    def translate_ids_3_to_1(self, ids: list) -> list:
        """
        Translate a list of 3-letter compound ids to 1-letter ids.
        Any ids that cannot be translated are returned as "X".

        Parameters
        ----------
        ids : list
            A list of 3-letter compound ids.

        Returns
        -------
        list
            A list of 1-letter compound ids.
        """
        new = []
        for i in ids:
            i = self.get(i, "id", "dict")
            i = i.get("one_letter_code")
            if i is None:
                i = "X"
            new.append(i)
        return new

    def translate_ids_1_to_3(self, ids: list) -> list:
        """
        Translate a list of 1-letter compound ids to 3-letter ids.
        Any unknown ids are replaced with "XXX".

        Parameters
        ----------
        ids : list
            A list of 1-letter compound ids.

        Returns
        -------
        list
            A list of 3-letter compound ids.
        """
        new = []
        for i in ids:
            i = self.get(i, "id", "dict")
            i = i.get("three_letter_code")
            if i is None:
                i = "XXX"
            new.append(i)
            return new

    def relabel_atoms(self, structure):
        """
        Relabel the atoms of a `Molecule` object to match the atom names of the given compound.

        Parameters
        ----------
        structure
            A `Molecule` object or biopython object holding at least one residue.

        Returns
        -------
        object
            The object with relabeled atoms.
        """
        imposer = SVDSuperimposer.SVDSuperimposer()
        for residue in structure.get_residues():
            # get a reference
            name = residue.resname
            ref = self.get(name, "id", "residue")
            if ref is None:
                warnings.warn(f"Could not find residue '{name}' in PDBECompounds.")
                continue

            residue_coords = np.asarray([atom.coord for atom in residue.get_atoms()])
            ref_coords = np.asarray([atom.coord for atom in ref.get_atoms()])
            imposer.set(ref_coords, residue_coords)
            imposer.run()

            # get the residue coordinates
            new_coords = imposer.get_transformed()

            ref_atom_ids = [atom.id for atom in ref.get_atoms()]
            for idx, atom in enumerate(list(residue.get_atoms())):
                _new = new_coords[idx]
                min_idx = np.argmin(np.linalg.norm(ref_coords - _new, axis=1))
                atom.id = ref_atom_ids[min_idx]

        return structure

    def _get(self, q: str, by: str = "id"):
        """
        Get a dictionary of compounds that match the given criteria.

        Parameters
        ----------
        q : str
            The query to search for.
        by : str, optional
            The type of query, by default "id". Possible values are:
            - "id": search by compound id
            - "name": search by compound name (must match any available synonym exactly)
            - "formula": search by compound formula
            - "smiles": search by compound smiles (this also works for inchi)

        Returns
        -------
        dict
            A dictionary of compounds that match the given criteria.
        """
        if by == "id":
            _q = self._compounds.get(q)
            if _q is None:
                return {}
            return {q: _q}
        elif by == "name":
            return {k: v for k, v in self._compounds.items() if q.lower() in v["names"]}
        elif by == "formula":
            return {
                k: v
                for k, v in self._compounds.items()
                if q.upper().replace(" ", "") == v["formula"]
            }
        elif by == "smiles":
            return {k: v for k, v in self._compounds.items() if q in v["descriptors"]}
        else:
            raise ValueError(f"Invalid search type: {by}")

    def _setup_dictionaries(self, data_dict):
        """
        Fill in the dictionaries with the appropriate data

        Parameters
        ----------
        data_dict : dict
            A dictionary of data from the cif file.
        """
        for key, value in data_dict.items():
            comp = value["_chem_comp"]

            for k in __to_delete__["_chem_comp"]:
                if k in comp:
                    comp.pop(k)
            for k in __to_float__["_chem_comp"]:
                comp[k] = float(comp[k])
            for k in comp:
                if comp[k] == "?":
                    comp[k] = None

            comp["descriptors"] = value["_pdbx_chem_comp_descriptor"]["descriptor"]

            comp["names"] = [comp["name"]]
            synonyms = value.get("_pdbx_chem_comp_synonyms", None)
            if synonyms:
                comp["names"].extend(synonyms["name"])
            identifiers = value.get("_pdbx_chem_comp_identifier", None)
            if identifiers:
                comp["names"].extend(identifiers["identifier"])
            comp["names"] = set(i.lower() for i in comp["names"])

            if comp["formula"] is not None:
                comp["formula"] = comp["formula"].replace(" ", "")

            self._compounds[key] = comp

            atoms = value["_chem_comp_atom"]
            bonds = value["_chem_comp_bond"]

            atoms["pdbx_model_Cartn_x_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_x_ideal"]
            ]
            atoms["pdbx_model_Cartn_y_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_y_ideal"]
            ]
            atoms["pdbx_model_Cartn_z_ideal"] = [
                i.replace("?", "NaN") for i in atoms["pdbx_model_Cartn_z_ideal"]
            ]
            atoms["charge"] = [i.replace("?", "NaN") for i in atoms["charge"]]

            pdb = {
                "atoms": {
                    "full_ids": [i.replace(",", "'") for i in atoms["atom_id"]],
                    "ids": [
                        i.replace(",", "'") for i in atoms["pdbx_component_atom_id"]
                    ],
                    "serials": np.array(atoms["pdbx_ordinal"], dtype=int),
                    "coords": np.array(
                        [
                            (float(i), float(j), float(k))
                            for i, j, k in zip(
                                atoms["pdbx_model_Cartn_x_ideal"],
                                atoms["pdbx_model_Cartn_y_ideal"],
                                atoms["pdbx_model_Cartn_z_ideal"],
                            )
                        ]
                    ),
                    "elements": atoms["type_symbol"],
                    "charges": np.array(atoms["charge"], dtype=float),
                    "residue": atoms["pdbx_component_comp_id"],
                },
                "bonds": {
                    "bonds": [
                        (a.replace(",", "'"), b.replace(",", "'"))
                        for a, b in zip(bonds["atom_id_1"], bonds["atom_id_2"])
                    ],
                    "orders": np.array(bonds["value_order"]),
                },
            }
            self._pdb[key] = pdb

    def _molecule(self, compound: dict) -> Molecule:
        """
        Make a biobuild Molecule from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        Molecule
            A biobuild Molecule.
        """
        struct = self._make_structure(compound)
        mol = Molecule(struct)
        pdb = self._pdb[compound["id"]]
        for bdx in range(len(pdb["bonds"]["bonds"])):
            bond = pdb["bonds"]["bonds"][bdx]
            for _ in range(_bond_order_rev_map.get(pdb["bonds"]["orders"][bdx])):
                mol.add_bond(*bond)

        return mol

    def _make_structure(self, compound: dict) -> "bio.PDB.Structure":
        """
        Make a structure from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        bio.PDB.Structure
            A structure.
        """
        if len(compound) == 0:
            return None

        struct = bio.Structure.Structure(compound["id"])
        model = bio.Model.Model(0)
        struct.add(model)
        chain = bio.Chain.Chain("A")
        model.add(chain)
        res = self._make_residue(compound)
        chain.add(res)
        # add atoms to residue (at the end to ensure that the residue's full id is propagated)
        self._fill_residue(res, compound)
        return struct

    def _fill_residue(self, residue, compound: dict):
        """
        Fill a residue with atoms from a compound.

        Parameters
        ----------
        residue : bio.PDB.Residue
            A residue.
        compound : dict
            A dictionary of a compound.
        """
        pdb = self._pdb.get(compound["id"], None)
        if pdb is None:
            raise ValueError("No pdb data for compound.")

        atoms = pdb["atoms"]
        for i in range(len(atoms["ids"])):
            atom = self._make_atom(
                # SWITCHED FULL_ID AND ID<
                atoms["full_ids"][i],
                atoms["serials"][i],
                atoms["coords"][i],
                atoms["elements"][i],
                atoms["ids"][i],
                atoms["charges"][i],
            )
            residue.add(atom)

    def _make_residue(self, compound: dict) -> "bio.PDB.Residue":
        """
        Make a residue from a compound.

        Parameters
        ----------
        compound : dict
            A dictionary of a compound.

        Returns
        -------
        bio.PDB.Residue
            A residue.
        """
        if len(compound) == 0:
            return None

        res = bio.Residue.Residue(("H_" + compound["id"], 1, " "), compound["id"], " ")

        return res

    def _make_atom(
        self,
        id: str,
        serial: int,
        coords: np.ndarray,
        element: str,
        full_id: str,
        charge: float,
    ) -> "bio.PDB.Atom":
        """
        Make an atom.

        Parameters
        ----------
        id : str
            The atom id.
        serial : int
            The atom serial.
        coords : np.ndarray
            The atom coordinates.
        element : str
            The atom element.
        full_id : str
            The atom full id.
        Returns
        -------
        bio.PDB.Atom
            An atom.
        """
        atom = bio.Atom.Atom(
            id,
            coord=coords,
            serial_number=serial,
            bfactor=0.0,
            occupancy=0.0,
            element=element,
            fullname=full_id,
            altloc=" ",
            pqr_charge=charge,
        )
        return atom

    def __getitem__(self, key):
        return self._compounds[key], self._pdb[key]

    def __len__(self):
        return len(self._compounds)

    def __iter__(self):
        for key in self._compounds:
            yield key, self._compounds[key], self._pdb[key]


__all__ = [
    "PDBECompounds",
    "set_default_compounds",
    "get_default_compounds",
    "restore_default_compounds",
    "get_compound",
    "has_compound",
    "read_compounds",
    "load_compounds",
    "save_compounds",
]


if __name__ == "__main__":
    f = "support/pdbe_compounds/test.cif"
    compounds = PDBECompounds.from_file(f)

    from biobuild.core import Molecule

    real = Molecule.from_compound("GLC")

    scrambled = Molecule.from_compound("GLC")

    counts = {"C": 0, "H": 0, "O": 0, "N": 0, "S": 0, "P": 0}
    for atom in scrambled.atoms:
        counts[atom.element] += 1
        atom.id = atom.element + str(counts[atom.element])

    print(scrambled.atoms)
    compounds.relabel_atoms(scrambled)
    print(scrambled.atoms)

    from biobuild.utils.visual import MoleculeViewer3D

    v = MoleculeViewer3D(scrambled)
    v.show()
