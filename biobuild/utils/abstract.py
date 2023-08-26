"""
Abstract classes for storing force field data from CHARMM topology and parameter files
"""

from typing import Union
import attr
import Bio.PDB as bio
import numpy as np


class AbstractEntity:
    """
    A representation of a single collection entity (Residue, Patch, etc.)
    """

    def __init__(self, id=None):
        self.id = id
        self._atoms = {}
        self.bonds = []

    @property
    def atoms(self):
        return list(self._atoms.values())

    def get_atom(self, id) -> "AbstractAtom":
        """
        Get an atom by its ID

        Parameters
        ----------
        id: str
            The ID or type of the atom

        Returns
        -------
        atom: AbstractAtom
            The atom with the given ID.
            If no atom is found, None is returned.
        """
        if isinstance(id, (list, tuple)):
            return [self.get_atom(i) for i in id]
        atom = self._atoms.get(id, None)
        return atom

    def get_atoms_by_type(self, _type):
        """
        Get a list of atoms by their type

        Parameters
        ----------
        _type: str
            The type of the atom

        Returns
        -------
        atoms: list
            A list of atoms with the given type.
        """
        atoms = [i for i in self._atoms.values() if i.type == _type]
        return atoms

    def has_atom(self, _atom):
        """
        Check if the entity has an atom

        Parameters
        ----------
        _atom: AbstractAtom
            The atom to check for

        Returns
        -------
        has_atom: bool
            True if the atom is in the entity, False otherwise
        """
        if isinstance(_atom, AbstractAtom):
            _atom = _atom.id
        return _atom in self._atoms

    def get_bond(self, id1, id2) -> "AbstractBond":
        """
        Get a bond by its atom IDs
        """
        if isinstance(id1, AbstractAtom):
            id1 = id1.id
        if isinstance(id2, AbstractAtom):
            id2 = id2.id

        # ----------------------------     FUTURE FIX    ----------------------------
        # The bond.atom1.id assumes that the atom objects actually have an `id` attribute.
        # If only strings were provided as atoms then this will fail.
        # ----------------------------     FUTURE FIX    ----------------------------
        for bond in self.bonds:
            if bond.atom1.id == id1 and bond.atom2.id == id2:
                return bond
            elif bond.atom1.id == id2 and bond.atom2.id == id1:
                return bond
        return None

    def add_atom(self, atom):
        """
        Add an atom to the residue
        """
        self._atoms[atom.id] = atom

    def add_bond(self, bond):
        """
        Add a bond to the residue
        """
        self.bonds.append(bond)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.id})"


class AbstractEntity_with_IC(AbstractEntity):
    def __init__(self, id=None):
        super().__init__(id)
        self.internal_coordinates = []

    @property
    def IC_atom_ids(self):
        """
        Returns a set of all atom IDs of all atoms for which the patch also stores
        internal coordinates.
        """
        ids = set()
        for ic in self.internal_coordinates:
            if isinstance(ic.atom1, str):
                ids.update(ic.atoms)
                continue
            ids.update(ic.ids)
        return ids

    @property
    def has_IC(self):
        return len(self.internal_coordinates) > 0

    def add_internal_coordinates(self, ic):
        """
        Add an internal coordinate to the residue
        """
        self.internal_coordinates.append(ic)

    def add_ic(self, ic):
        """
        Add an internal coordinate to the residue
        """
        self.add_internal_coordinates(ic)

    def get_internal_coordinates(self, *ids, mode: str = "exact"):
        """
        Get internal coordinates by their atom IDs

        Parameters
        ----------
        ids: str
            The IDs of the atoms in the internal coordinate (this can also be any data object that has an `id` attribute)
        mode: str
            The mode to use for matching the internal coordinates.
            Supported modes are:
            - `exact`: The internal coordinate must match the given atom IDs exactly (requires that four ids are given)
            - `partial`: The internal coordinate must match the given atom IDs where they are provided, but wildcards can be set to None (requires that four ids or None are given in order)
            - `multi_partial`: The internal coordinate must match any given atom IDs where they are provided, but wildcards can be set to None (requires that four ids or None are given in order)
            - `anywhere`: The internal coordinate must contain all of the given atom IDs in any order (requires at least one id)
            - `anywhere_partial`: The internal coordinate must contain any of the given atom IDs in any order (requires at least one id)

        Returns
        -------
        ics: list
            A list of internal coordinates
        """

        if len(ids) == 0:
            return self.internal_coordinates

        ids = tuple(i.id if hasattr(i, "id") else i for i in ids)

        if mode == "exact":
            if len(ids) != 4:
                raise ValueError(
                    "Exact mode requires that four ids are given to match the internal coordinates"
                )

            for ic in self.internal_coordinates:
                if ids == ic.ids or ids[::-1] == ic.ids:
                    return [ic]
            return []

        elif mode == "partial":
            if len(ids) != 4:
                raise ValueError(
                    "Partial mode requires that four ids or None are given to match the internal coordinates"
                )

            ids = np.array(ids)
            mask = ids != None

            ics = [
                ic
                for ic in self.internal_coordinates
                if np.all(ids[mask] == np.array(ic.ids)[mask])
            ]
            return ics

        elif mode == "anywhere":
            ids = set(ids)
            ics = [ic for ic in self.internal_coordinates if ids.issubset(set(ic.ids))]
            return ics

        elif mode == "anywhere_partial":
            ids = set(ids)
            ics = [
                ic
                for ic in self.internal_coordinates
                if len(ids.intersection(set(ic.ids))) != 0
            ]
            return ics

        elif mode == "multi_partial":
            if len(ids) != 4:
                raise ValueError(
                    "Multi partial mode requires that four ids or None are given to match the internal coordinates"
                )

            ids = np.array(ids)
            mask = ids != None

            ics = [
                ic
                for ic in self.internal_coordinates
                if np.any(ids[mask] == np.array(ic.ids)[mask])
            ]
            return ics

        else:
            raise ValueError(f"Unknown mode {mode}")


@attr.s(hash=True)
class AbstractAtom:
    """
    A representation of a single Atom
    """

    id = attr.ib(type=str, hash=True)
    type = attr.ib(default=None, type=str, hash=True)
    charge = attr.ib(default=None, type=float, repr=False)
    mass = attr.ib(default=None, type=float, repr=False)
    _element = attr.ib(default=None, type=str, repr=False)
    _parent = attr.ib(default=None, repr=False)
    is_wildcard = attr.ib(default=False, type=bool, repr=False)
    coord = attr.ib(default=None, type=np.ndarray, repr=False)

    @property
    def element(self):
        if self._element is None:
            return self.id[0]
        return self._element

    def to_biopython(self, serial_number: int = 1):
        """
        Returns a Bio.PDB.Atom object

        Parameters
        ----------
        serial_number: int
            The serial number of the atom
        """
        if self.coord:
            coord = self.coord
        else:
            coord = np.full(3, np.nan)
        new = bio.Atom.Atom(
            self.id,
            coord=coord,
            serial_number=serial_number,
            bfactor=0.0,
            occupancy=0.0,
            altloc="",
            fullname=self.id,
            element=self.element,
            pqr_charge=self.charge,
        )
        return new

    def get_parent(self):
        """
        Get the parent of the atom
        """
        return self._parent

    def set_parent(self, obj):
        """
        Set the parent of the atom
        """
        self._parent = obj

    def __repr__(self):
        return f"AbstractAtom({self.id})"


@attr.s(hash=True)
class AbstractBond:
    """
    A representation of a single Bond between two atoms (or atom types)
    optionally, a bond spring constant (K) and length can be supplied.
    """

    atom1 = attr.ib(type=Union[AbstractAtom, str], hash=True)
    atom2 = attr.ib(type=Union[AbstractAtom, str], hash=True)

    K = attr.ib(type=float, default=None)
    length = attr.ib(type=float, default=None)

    @property
    def atoms(self):
        return self.atom1, self.atom2

    def get_parent(self):
        """
        Get the parent of the bond (i.e. it's residue)
        """
        if isinstance(self.atom1, str):
            return None
        return self.atom1.get_parent()

    def find_atoms(self, obj):
        """
        Find the atoms of the bond in an object.
        This will return a tuple of identified atoms with the
        same id if they exist in the object, None for any atom
        that was not found.
        """
        if isinstance(self.atom1, str):
            atom1 = obj.get_atom(self.atom1)
        else:
            atom1 = obj.get_atom(self.atom1.id)

        if isinstance(self.atom2, str):
            atom2 = obj.get_atom(self.atom2)
        else:
            atom2 = obj.get_atom(self.atom2.id)

        return atom1, atom2

    def migrate_atoms(self, obj):
        """
        Migrate the atoms of the bond to a new object.
        This will update the atom1 and atom2 attributes of the bond
        to point to the atoms in the new object if they can be found.
        """
        atom1, atom2 = self.find_atoms(obj)
        self.atom1 = atom1 if atom1 else self.atom1
        self.atom2 = atom2 if atom2 else self.atom2

    def __getitem__(self, key):
        return self.atoms[key]

    def __iter__(self):
        return iter(self.atoms)

    def __contains__(self, item):
        return item in self.atoms

    def __repr__(self):
        id1 = getattr(self.atom1, "id", self.atom1)
        id2 = getattr(self.atom2, "id", self.atom2)
        return f"AbstractBond({id1}, {id2})"
