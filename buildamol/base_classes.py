"""
The base_classes are deriviatives of the original Biopython classes, but with
the change that they use a UUID4 as their identifier (full_id) instead of a hierarchical
tuple. This makes each object unique and allows for easy comparison where `a == b` is akin to `a is b`.
Consequently, the `__hash__` method is overwritten to use the UUID4 as the hash.

.. warning::
    
    Each class has its own `copy` method that returns a deep copy of the object with a new UUID4. So `a.copy() == a` is `False`, while a standard `deepcopy(a) == a` is `True` since the UUID4 will not have been updated automatically.

Converting to and from `biopython`
----------------------------------

Each BuildAMol class can be generated from a biopython class using the `from_biopython` class method. And each BuildAMol class has a `to_biopython` method that returns the pure-biopython equivalent.
It is important to note, that for most purposes, however, the BuildAMol classes should work fine as trop-in replacements for the original biopython classes. 

.. code-block:: python

    import Bio.PDB as bio
    from buildamol.base_classes import Atom

    bio_atom = bio.Atom("CA", (0, 0, 0))
    atom = Atom.from_biopython(bio_atom)

    assert atom == bio_atom # False since atom uses a UUID4 as its identifier
    assert atom.to_biopython() == bio_atom # True 

The conversion from and to biopython works hierarchically, so if an entire biopython structure is converted to BuildAMol
then all atoms, residues, chains and models will be converted to their BuildAMol equivalents.

.. code-block:: python

    import Bio.PDB as bio
    from buildamol.base_classes import Structure

    bio_structure = bio.PDBParser().get_structure("test", "test.pdb")
    structure = Structure.from_biopython(bio_structure)

    atoms = list(structure.get_atoms())
    bio_atoms = list(bio_structure.get_atoms())
    assert len(atoms) == len(bio_atoms) # True
    
"""

from copy import deepcopy
from typing import Union, List

# from uuid import uuid4
import Bio.PDB as bio
import numpy as np
import periodictable as pt


__all__ = ["Atom", "Residue", "Chain", "Model", "Structure", "Bond"]

__global_counters__ = {
    "A": 1,
    "R": 1,
    "C": 1,
    "M": 1,
    "S": 1,
}

__global_element_counters__ = {}


class ID:
    """
    The base class for BuildAMol's internal object identification.
    All classes that inheret from this class will be recorded as unique objects.
    """

    __global_idx__ = 0

    def __init__(self):
        self.__id = ID.__global_idx__ + 1
        ID.__global_idx__ += 1

    def copy(self):
        new = deepcopy(self)
        new._new_id()
        return new

    def get_id(self):
        return self.__id

    def has_id(self, id):
        return id in self.child_dict

    def _new_id(self):
        self.__id = ID.__global_idx__ + 1
        ID.__global_idx__ += 1

    def _adopt_id(self, id):
        self.__id = id

    def __hash__(self):
        return self.__id

    def __eq__(self, other):
        if not isinstance(other, ID):
            return False
        return self.__id == other.__id

    def __ne__(self, other):
        if not isinstance(other, ID):
            return True
        return self.__id != other.__id


class Atom(ID, bio.Atom.Atom):
    """
    An Atom object that inherits from Biopython's Atom class.

    Parameters
    ----------
    id : str
        The atom identifier
    coord : ndarray
        The atom coordinates
    serial_number : int, optional
        The atom serial number. The default is 1.
    bfactor : float, optional
        The atom bfactor. The default is 0.0.
    occupancy : float, optional
        The atom occupancy. The default is 1.0.
    fullname : str, optional
        The atom fullname. The default is None, in which case the id is used again.
    element : str, optional
        The atom element. The default is None, in which case it is inferred based on the id.
    altloc : str, optional
        The atom altloc. The default is " ".
    pqr_charge : float, optional
        The atom pqr_charge. The default is None.
    radius : float, optional
        The atom radius. The default is None.
    """

    __slots__ = (
        "id",
        "parent",
        "name",
        "fullname",
        "coord",
        "mass",
        "serial_number",
        "bfactor",
        "occupancy",
        "altloc",
        "element",
        "pqr_charge",
        "radius",
        "level",
        "disordered_flag",
        "anisou_array",
        "siguij_array",
        "sigatm_array",
        "xtra",
        "_sorting_keys",
    )

    def __init__(
        self,
        id: str,
        coord: "np.ndarray",
        serial_number: int = 1,
        bfactor: float = 0.0,
        occupancy: float = 1.0,
        fullname: str = None,
        element: str = None,
        altloc=" ",
        pqr_charge=None,
        radius=None,
    ):
        if not fullname:
            fullname = id
        ID.__init__(self)
        if element:
            element = element.upper()
        bio.Atom.Atom.__init__(
            self,
            id,
            np.asarray(coord, dtype=np.float64),
            bfactor,
            occupancy,
            altloc,
            fullname,
            serial_number,
            element,
            pqr_charge,
            radius,
        )
        self.level = "A"

    @classmethod
    def new(
        cls,
        element_or_id: str,
        coord: "np.ndarray" = None,
        generate_id: bool = True,
        **kwargs,
    ) -> "Atom":
        """
        Create a blank atom with a given element and coordinates.

        Parameters
        ----------
        element_or_id : str
            The atom element. If the element is not found in the periodic table, it will be used as the atom id.
            To ensure the correct element is assigned to the atom either pass it as it as keyword argument directly or
            to let BuildAMol correctly infer the element use one of the following patterns with your atom id:
                - <element+><number> -> C1, O2, CA2 (calcium id=CA2)
                - <number><element+> ->  1C, 2O, 1CA (calcium id=1CA)
                - <space><element><string|number+> -> " CA" (Carbon id=CA), " ND2  " (Nitrongen id=ND2), " OXT" (Oxygen id=OXT)
                - <element+><space+> -> "CA " (Calcium id=CA), "FE " (Iron id=FE)
                - <element+>_<string|number+> -> "CA_" (Calcium id=CA), "C_A" (Carbon id=CA)
            (the + indicates multiple characters)
        coord : ndarray, optional
            The atom coordinates. The default is None.
        generate_id : bool, optional
            Whether to automatically generate a new id for the atom to avoid identically named atoms. The default is True.

        **kwargs
            Additional keyword arguments to pass to the Atom initializer

        Returns
        -------
        Atom
            The blank atom.
        """
        if not generate_id:
            _id = element_or_id
        if coord is None:
            coord = np.zeros(3)
        if pt.elements.__dict__.get(element_or_id) or pt.elements.__dict__.get(
            element_or_id.title()
        ):
            element = element_or_id.upper()
            __global_element_counters__[element] = (
                __global_element_counters__.get(element, 0) + 1
            )
            id = f"{element}{__global_element_counters__[element]}"
        else:
            id, element = Atom._infer_element_from_id(element_or_id)

        element = kwargs.pop("element", element)
        if not generate_id:
            id = _id
        return cls(id=id, coord=coord, element=element, **kwargs)

    @classmethod
    def from_element(self, element: str, **kwargs):
        """
        Create a blank atom with a given element and coordinates.

        Parameters
        ----------
        element : str
            The atom element.
        **kwargs
            Additional keyword arguments to pass to the new Atom initializer.
        """
        return Atom.new(element, **kwargs)

    @property
    def charge(self):
        """
        The atom charge.
        """
        return self.pqr_charge

    @charge.setter
    def charge(self, value):
        self.pqr_charge = value

    @property
    def full_id(self):
        p = self.get_parent()
        if p:
            return (*p.get_full_id(), (self.name, self.altloc))
        else:
            return (None, None, None, None, (self.id, self.altloc))

    @full_id.setter
    def full_id(self, value):
        pass

    def matches(
        self, other, include_id: bool = True, include_coord: bool = False
    ) -> bool:
        """
        Check if the atom matches another atom.
        This will return True if the two atoms have the same element, id, parent-residue name, and altloc.
        """
        verdict = (
            self.element == other.element
            and self.parent.resname == other.parent.resname
            and self.altloc == other.altloc
        )
        if include_id:
            verdict = verdict and self.id == other.id

        if include_coord:
            verdict = verdict and np.allclose(self.coord, other.coord)
        return verdict

    def equals(self, other, include_coord: bool = False) -> bool:
        """
        Check if the atom is equal to another atom.
        This will return True if the two atoms match and have same the parent-serial number.
        """
        return self.matches(other, include_id=True, include_coord=include_coord) and (
            self.parent.id[1] == other.parent.id[1]
        )

    @classmethod
    def from_biopython(cls, atom) -> "Atom":
        """
        Convert a Biopython atom to an Atom object

        Parameters
        ----------
        atom
            The Biopython atom

        Returns
        -------
        Atom
            The Atom object
        """
        return cls(
            atom.id,
            atom.coord,
            atom.serial_number,
            atom.bfactor,
            atom.occupancy,
            atom.fullname,
            atom.element,
            atom.altloc,
            atom.pqr_charge,
            atom.radius,
        )

    def to_biopython(self):
        """
        Convert the Atom object to a Biopython atom

        Returns
        -------
        Atom
            The Biopython atom
        """
        return bio.Atom.Atom(
            self.id,
            self.coord,
            self.bfactor,
            self.occupancy,
            self.altloc,
            self.fullname,
            self.serial_number,
            self.element.upper(),
            self.pqr_charge,
            self.radius,
        )

    @property
    def atomic_number(self):
        """
        The atomic number of the atom's element.
        """
        return pt.elements.symbol(self.element).number

    def set_id(self, id):
        """
        Set the atom identifier.

        Parameters
        ----------
        id : str
            The identifier to set.
        """
        self.id = id
        return self

    def set_element(self, element, adjust_id: bool = True):
        """
        Set the atom element.

        Parameters
        ----------
        element : str
            The element to set.
        adjust_id : bool, optional
            Whether to adjust the atom id to the new element. The default is True.
        """
        current = self.element
        self.element = element.upper()
        self.mass = pt.elements.symbol(element).mass
        if adjust_id:
            count = __global_element_counters__.get(self.element, None)
            if not count:
                __global_element_counters__[self.element] = 1
                count = 0
            count += 1
            self.id = self.element + str(count)
        else:
            self.id = self.id.replace(current, self.element)
        return self

    def move(self, vector):
        """
        Move the atom by a vector.

        Parameters
        ----------
        vector : ndarray
            The vector to move the atom by.
        """
        self.coord += vector
        return self

    @staticmethod
    def _infer_element_from_id(id: str) -> str:
        """
        A method to infer the element from the atom id.

        Rules that are considered (+ indicates multiple characters):
        - <element+><number>, e.g. C1, O2, CA2 (calcium id=CA2)
        - <number><element+>, e.g. 1C, 2O, 1CA (calcium id=1CA)
        - <space><element><string|number+> -> " CA" (Carbon id=CA), " ND2  " (Nitrongen id=ND2), " OXT" (Oxygen id=OXT)
        - <element+><space+> -> "CA " (Calcium), "FE " (Iron)
        - <element+>_<string|number+> -> "CA_" (Calcium id=CA), "C_A" (Carbon id=CA)

        """
        id = id.upper()
        element = None

        # <number><element+>, e.g. 1C, 2O, 3N, 1H
        # or
        # <element+><number>, e.g. C1, O2, N3, H1
        if id[0].isdigit() or (id[0] != " " and id[-1].isdigit()):
            element = id.translate(str.maketrans("", "", "0123456789"))

        # <space><element><string+> -> " CA" (Carbon id=CA), " ND2" (Nitrongen id=ND2), " OXT" (Oxygen id=OXT)
        elif id[0] == " ":
            id = id.strip()
            element = id[0]

        # <element+>_<string+> -> "CA_" (Calcium id=CA), "C_A" (Carbon id=CA)
        elif "_" in id:
            element = id.split("_")[0]
            id = id.replace("_", "")

        # <element+><space+> -> "CA " (Calcium), "FE " (Iron)
        elif id[-1] == " ":
            element = id = id.strip()

        return id, element

    def __repr__(self):
        return f"Atom({self.id}, {self.serial_number})"

    def __lt__(self, other):
        return (self.serial_number < other.serial_number) or (
            pt.elements.symbol(self.element.title()).number
            < pt.elements.symbol(other.element.title()).number
        )

    def __gt__(self, other):
        return self.serial_number > other.serial_number or (
            pt.elements.symbol(self.element.title()).number
            > pt.elements.symbol(other.element.title()).number
        )

    def __le__(self, other):
        return self.serial_number <= other.serial_number or (
            pt.elements.symbol(self.element.title()).number
            <= pt.elements.symbol(other.element.title()).number
        )

    def __ge__(self, other):
        return self.serial_number >= other.serial_number or (
            pt.elements.symbol(self.element.title()).number
            >= pt.elements.symbol(other.element.title()).number
        )

    def __hash__(self):
        return ID.__hash__(self)

    # def __eq__(self, other):
    #     return self.serial_number == other.serial_number and (
    #         pt.elements.symbol(self.element.title()).number
    #         == pt.elements.symbol(other.element.title()).number
    #     )

    # def __ne__(self, other):
    #     return self.serial_number != other.serial_number or (
    #         pt.elements.symbol(self.element.title()).number
    #         != pt.elements.symbol(other.element.title()).number
    #     )


class Residue(ID, bio.Residue.Residue):
    """
    A Residue object that inherits from Biopython's Residue class.

    Parameters
    ----------
    resname : str
        The residue name
    segid : str
        The residue segid.
    icode : int
        The residue icode.
        This is the residue serial number.
    """

    __slots__ = (
        "level",
        "disordered",
        "resname",
        "segid",
        "internal_coord",
        "_id",
        "parent",
        "child_list",
        "child_dict",
        "xtra",
        "_coord",
    )

    def __init__(self, resname, segid=" ", icode=1):
        ID.__init__(self)
        bio.Residue.Residue.__init__(
            self, ("H_" + resname, icode, segid), resname, segid
        )
        self.level = "R"
        self.serial_number = icode
        self._coord = None

    @classmethod
    def new(cls, resname: str, segid: str = " ", icode: int = None) -> "Residue":
        """
        Create a blank residue with a given name and segid.

        Parameters
        ----------
        resname : str
            The residue name.
        segid : str, optional
            The residue segid. The default is " ".
        icode : int, optional
            The residue icode. The default is None.

        Returns
        -------
        Residue
            The blank residue.
        """
        if icode is None:
            icode = __global_counters__["R"]
            __global_counters__["R"] += 1
        return cls(resname, segid, icode)

    @property
    def id(self):
        return ("H_" + self.resname, self.serial_number, self.segid)

    @id.setter
    def id(self, value):
        pass

    @property
    def full_id(self):
        p = self.get_parent()
        if p:
            return (*p.get_full_id(), self.id)
        else:
            return (
                None,
                None,
                None,
                self.id,
            )

    @full_id.setter
    def full_id(self, value):
        pass

    @property
    def coord(self):
        if self._coord is None:
            return self.center_of_mass()
        else:
            return self._coord

    @coord.setter
    def coord(self, value):
        self._coord = value

    def get_coord(self) -> "np.ndarray":
        """
        Get the center of mass of the residue.
        """
        return self.coord

    def set_coord(self, value):
        """
        Set the center of mass of the residue.
        """
        self.coord = value

    def get_coords(self) -> "np.ndarray":
        """
        Get the coordinates of all atoms in the residue.
        """
        return np.array([atom.coord for atom in self.get_atoms()])

    def get_atom(self, atom: Union[str, int]) -> Atom:
        """
        Get an atom by its name or serial number.

        Parameters
        ----------
        atom : str or int
            The atom name or serial number.

        Returns
        -------
        Atom
            The atom.
        """
        if isinstance(atom, int):
            return next((i for i in self.child_list if i.serial_number == atom), None)
        elif isinstance(atom, str):
            return next((i for i in self.child_list if i.id == atom), None)
        else:
            raise TypeError(f"atom must be either a string or an integer, not {atom=}")

    def matches(self, other) -> bool:
        """
        Check if the residue matches another residue.
        This will return True if the two residues have the same resname, segid, and parent-chain id.
        """
        return (
            self.resname == other.resname
            and self.segid == other.segid
            and self.parent.id == other.parent.id
        )

    def equals(self, other, include_serial: bool = False) -> bool:
        """
        Check if the residue is equal to another residue.
        This will check if the two residues are in the same parent and if all atoms are matching.
        """
        verdict = self.matches(other) and all(
            i.matches(j) for i, j in zip(self.get_atoms(), other.get_atoms())
        )
        if include_serial:
            verdict = verdict and self.serial_number == other.serial_number
        return verdict

    def link(self, atom):
        """
        Softlink an atom into this residue's child_list without touching the atom's own parent references.
        """
        self.child_dict[atom.get_id()] = atom
        self.child_list.append(atom)

    def unlink(self, atom):
        """
        Unlink an atom from this residue's child_list without touching the atom's own parent references.
        """
        del self.child_dict[atom.get_id()]
        self.child_list.remove(atom)

    # def add(self, atom):
    #     if atom.get_id() not in self.child_dict:
    #         self.child_list.append(atom)
    #         self.child_dict[atom.id] = atom
    #         atom.set_parent(self)

    # def remove(self, atom):
    #     if atom.get_id() in self.child_dict:
    #         self.child_list.remove(atom)
    #         del self.child_dict[atom.id]
    #         atom.set_parent(None)

    @classmethod
    def from_biopython(cls, residue) -> "Residue":
        """
        Convert a BioPython Residue object to a Residue object.

        Parameters
        ----------
        residue : BioPython Residue object
            The residue to convert.

        Returns
        -------
        Residue
            The converted residue
        """
        new = cls(residue.id[0], residue.id[1], residue.id[-1])
        for atom in residue.get_atoms():
            new.add(Atom.from_biopython(atom))
        return new

    def to_biopython(self) -> bio.Residue.Residue:
        """
        Convert a Residue object to a pure BioPython Residue object.

        Returns
        -------
        bio.Residue.Residue
            The converted residue.
        """
        new = bio.Residue.Residue(self.id, self.resname, self.segid)
        for atom in self.get_atoms():
            new.add(atom.to_biopython())
        return new

    def add(self, atom):
        if not isinstance(atom, Atom):
            atom = Atom.from_biopython(atom)
        bio.Residue.Residue.add(self, atom)

    def move(self, vector):
        """
        Move the residue by a vector.

        Parameters
        ----------
        vector : ndarray
            The vector to move the residue by.
        """
        for atom in self.get_atoms():
            atom.move(vector)

    def copy(self) -> "Residue":
        """
        Return a deep copy of the residue with a new UUID4.

        Returns
        -------
        Residue
            The copied residue.
        """
        new = ID.copy(self)
        for atom in new.get_atoms():
            ID._new_id(atom)
        return new

    def __repr__(self):
        return f"Residue({self.resname}, {self.serial_number})"

    def __lt__(self, other):
        return self.serial_number < other.serial_number

    def __gt__(self, other):
        return self.serial_number > other.serial_number

    def __le__(self, other):
        return self.serial_number <= other.serial_number

    def __ge__(self, other):
        return self.serial_number >= other.serial_number

    # def __eq__(self, other):
    #     return (
    #         self.serial_number == other.serial_number
    #         and self.resname == other.resname
    #         and self.parent == other.parent
    #     )

    # def __ne__(self, other):
    #     return (
    #         self.serial_number != other.serial_number
    #         or self.resname != other.resname
    #         or self.parent != other.parent
    #     )


class Chain(ID, bio.Chain.Chain):
    """
    A Chain object that inherits from Biopython's Chain class.

    Parameters
    ----------
    id : str
        The chain identifier
    """

    __slots__ = (
        "level",
        "internal_coord",
        "_id",
        "parent",
        "child_list",
        "child_dict",
        "xtra",
    )

    def __init__(self, id):
        ID.__init__(self)
        super(bio.Chain.Chain, self).__init__(id)
        self.level = "C"

    @classmethod
    def new(cls, id: str) -> "Chain":
        """
        Create a blank chain with a given id.

        Parameters
        ----------
        id : str, optional
            The chain identifier. The default is None.

        Returns
        -------
        Chain
            The blank chain.
        """
        return cls(id)

    @property
    def full_id(self):
        p = self.get_parent()
        if p:
            return (*p.get_full_id(), self.id)
        else:
            return (None, None, self.id)

    @full_id.setter
    def full_id(self, value):
        pass

    def add(self, residue):
        if not isinstance(residue, Residue):
            residue = Residue.from_biopython(residue)
        bio.Chain.Chain.add(self, residue)

    def get_residue(self, residue: Union[str, int]) -> Residue:
        """
        Get a residue by its name or serial number.

        Note
        ----
        If there are multiple residues with the same name, the first one will be returned.

        Parameters
        ----------
        residue : str or int
            The residue name or serial number.

        Returns
        -------
        Residue
            The residue.
        """
        if isinstance(residue, int):
            return next(
                (i for i in self.child_list if i.serial_number == residue), None
            )
        elif isinstance(residue, str):
            return next((i for i in self.child_list if i.resname == residue), None)
        else:
            raise TypeError(
                f"residue must be either a string or an integer, not {residue=}"
            )

    def get_residues(self, residue: Union[str, int] = None) -> "List[Residue]":
        """
        Get all residues in the chain.

        Parameters
        ----------
        residue : str or int, optional
            The residue name or serial number to filter by. The default is None (normal behavior returns a generator of all residues)

        Returns
        -------
        List[Residue]
            The list of residues. If no residue argument is specified the default generator is returned.
        """
        if residue:
            if isinstance(residue, str):
                return [i for i in self.child_list if i.resname == residue]
            # makes no sense but for consistency
            elif isinstance(residue, int):
                return [i for i in self.child_list if i.serial_number == residue]
            else:
                raise TypeError(
                    f"residue must be either a string or an integer, not {residue=}"
                )
        return super().get_residues()

    def get_coords(self) -> "np.ndarray":
        """
        Get the coordinates of all atoms in the chain.
        """
        return np.array([atom.coord for atom in self.get_atoms()])

    def matches(self, other) -> bool:
        """
        Check if the chain matches another chain.
        This will return True if the two chains have matching residues.
        """
        return all(
            i.matches(j) for i, j in zip(self.get_residues(), other.get_residues())
        )

    def equals(self, other) -> bool:
        """
        Check if the chain is equal to another chain.
        This will check if the two chains have the same id, the same parent-model id, and have equal residues.
        """
        return (
            self.id == other.id
            and self.parent.id == other.parent.id
            and all(
                i.equals(j) for i, j in zip(self.get_residues(), other.get_residues())
            )
        )

    @classmethod
    def from_biopython(cls, chain) -> "Chain":
        """
        Convert a BioPython Chain object to a Chain object.

        Parameters
        ----------
        chain : BioPython Chain object
            The chain to convert.

        Returns
        -------
        Chain
            The converted chain.
        """
        new = cls(chain.id)
        for residue in chain.get_residues():
            new.add(Residue.from_biopython(residue))
        return new

    def to_biopython(self) -> bio.Chain.Chain:
        """
        Convert a Chain object to a pure BioPython Chain object.

        Parameters
        ----------
        with_children : bool, optional
            Whether to convert the residues of the chain as well. The default is True.

        Returns
        -------
        bio.Chain.Chain
            The converted chain.
        """
        new = bio.Chain.Chain(self.id)
        for residue in self.get_residues():
            new.add(residue.to_biopython())
        return new

    def move(self, vector):
        """
        Move the chain by a vector.

        Parameters
        ----------
        vector : ndarray
            The vector to move the chain by.
        """
        for residue in self.get_residues():
            residue.move(vector)

    def copy(self):
        """
        Return a deep copy of the chain with a new UUID4.

        Returns
        -------
        Chain
            The copied chain.
        """
        new = ID.copy(self)
        for residue in new.get_residues():
            ID._new_id(residue)
            for atom in residue.get_atoms():
                ID._new_id(atom)
        return new

    def link(self, residue):
        """
        Softlink a residue into this chain's child_list without touching the residue's own parent references.
        """
        self.child_dict[residue.get_id()] = residue
        self.child_list.append(residue)

    def unlink(self, residue):
        """
        Unlink a residue from this chain's child_list without touching the residue's own parent references.
        """
        del self.child_dict[residue.get_id()]
        self.child_list.remove(residue)

    def __repr__(self):
        return f"Chain({self._id})"

    def __lt__(self, other):
        return ord(self.id) < ord(other.id)

    def __gt__(self, other):
        return ord(self.id) > ord(other.id)

    def __le__(self, other):
        return ord(self.id) <= ord(other.id)

    def __ge__(self, other):
        return ord(self.id) >= ord(other.id)

    # def __eq__(self, other):
    #     return ord(self.id) == ord(other.id)

    # def __ne__(self, other):
    #     return ord(self.id) != ord(other.id)


class Model(bio.Model.Model, ID):
    """
    A Model object that inherits from Biopython's Model class.

    Parameters
    ----------
    id : int or str
        The model identifier
    """

    __slots__ = (
        "level",
        # "serial_num",
        "_id",
        "parent",
        "child_list",
        "child_dict",
        "xtra",
    )

    def __init__(self, id):
        ID.__init__(self)
        super(bio.Model.Model, self).__init__(id)
        self.level = "M"

    @classmethod
    def new(cls, id: int = None) -> "Model":
        """
        Create a blank model with a given id.

        Parameters
        ----------
        id : int
            The model identifier.

        Returns
        -------
        Model
            The blank model.
        """
        if id is None:
            id = __global_counters__["M"]
            __global_counters__["M"] += 1
        return cls(id)

    @property
    def serial_number(self):
        if isinstance(self.id, int):
            return self.id
        else:
            return ord(self.id)

    @serial_number.setter
    def serial_number(self, value):
        pass

    @property
    def serial_num(self):
        return self.serial_number

    @serial_num.setter
    def serial_num(self, value):
        pass

    @property
    def full_id(self):
        p = self.get_parent()
        if p:
            return (p.get_full_id(), self.id)
        else:
            return (None, self.id)

    @full_id.setter
    def full_id(self, value):
        pass

    def add(self, chain):
        if not isinstance(chain, Chain):
            chain = Chain.from_biopython(chain)
        bio.Model.Model.add(self, chain)

    def get_coords(self) -> "np.ndarray":
        """
        Get the coordinates of all atoms in the model.
        """
        return np.array([atom.coord for atom in self.get_atoms()])

    def move(self, vector):
        """
        Move the model by a vector.

        Parameters
        ----------
        vector : ndarray
            The vector to move the model by.
        """
        for chain in self.get_chains():
            chain.move(vector)

    def link(self, chain):
        """
        Softlink a chain into this model's child_list without touching the chain's own parent references.
        """
        self.child_dict[chain.get_id()] = chain
        self.child_list.append(chain)

    def unlink(self, chain):
        """
        Unlink a chain from this model's child_list without touching the chain's own parent references.
        """
        del self.child_dict[chain.get_id()]
        self.child_list.remove(chain)

    def matches(self, other) -> bool:
        """
        Check if the model matches another model.
        This will return True if the two models have matching chains.
        """
        return all(i.matches(j) for i, j in zip(self.get_chains(), other.get_chains()))

    def equals(self, other) -> bool:
        """
        Check if the model is equal to another model.
        This will return True if the two models have the same id, same parent-structure id, and have matching chains.
        """
        return (
            self.matches(other)
            and self.id == other.id
            and self.parent.id == other.parent.id
        )

    def copy(self):
        """
        Return a deep copy of the model with a new UUID4.

        Returns
        -------
        Model
            The copied model.
        """
        new = ID.copy(self)
        for chain in new.get_chains():
            ID._new_id(chain)
            for residue in chain.get_residues():
                ID._new_id(residue)
                for atom in residue.get_atoms():
                    ID._new_id(atom)
        return new

    @classmethod
    def from_biopython(cls, model):
        """
        Convert a BioPython Model object to a Model object.

        Parameters
        ----------
        model : BioPython Model object
            The model to convert.

        Returns
        -------
        Model
            The converted model.
        """
        new = cls(model.id)
        for chain in model.get_chains():
            new.add(Chain.from_biopython(chain, with_children=True))
        return new

    def to_biopython(self):
        """
        Convert a Model object to a pure BioPython Model object.

        Returns
        -------
        bio.Model.Model
            The converted model.
        """
        new = bio.Model.Model(self.id)
        for chain in self.get_chains():
            new.add(chain.to_biopython(with_children=True))
        return new

    def link(self, chain):
        """
        Softlink a chain into this model's child_list without touching the chain's own parent references.
        """
        self.child_dict[chain.get_id()] = chain
        self.child_list.append(chain)

    def __repr__(self):
        return f"Model({self._id})"

    # somehow the __eq__ was not inherited from ID...
    def __eq__(self, other):
        if not isinstance(other, Model):
            return False
        return self._ID__id == other._ID__id

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __le__(self, other):
        return self.id <= other.id

    def __ge__(self, other):
        return self.id >= other.id


class Structure(ID, bio.Structure.Structure):
    """
    A Structure object that inherits from Biopython's Structure class.

    Parameters
    ----------
    id : str
        The structure identifier
    """

    __slots__ = (
        "level",
        "_id",
        "parent",
        "child_list",
        "child_dict",
        "xtra",
    )

    def __init__(self, id):
        ID.__init__(self)
        super(bio.Structure.Structure, self).__init__(id)
        self.level = "S"

    @classmethod
    def new(cls, id: str) -> "Structure":
        """
        Create a blank structure with a given id.

        Parameters
        ----------
        id : str
            The structure identifier.

        Returns
        -------
        Structure
            The blank structure.
        """
        return cls(id)

    @property
    def full_id(self):
        return self.id

    @full_id.setter
    def full_id(self, value):
        pass

    def add(self, model):
        if not isinstance(model, Model):
            model = Model.from_biopython(model)
        bio.Structure.Structure.add(self, model)

    def get_coords(self) -> "np.ndarray":
        """
        Get the coordinates of all atoms in the structure.
        """
        return np.array([atom.coord for atom in self.get_atoms()])

    def matches(self, other) -> bool:
        """
        Check if the structure matches another structure.
        This will return True if the two structures have the same id.
        """
        return self.id == other.id

    def equals(self, other) -> bool:
        """
        Check if the structure is equal to another structure.
        This will return True if the two structures have the same id and have equal models.
        """
        return self.matches(other) and all(
            i.equals(j) for i, j in zip(self.get_models(), other.get_models())
        )

    def copy(self):
        """
        Return a deep copy of the structure with a new UUID4.

        Returns
        -------
        Structure
            The copied structure.
        """
        new = ID.copy(self)
        for model in new.get_models():
            ID._new_id(model)
            for chain in model.get_chains():
                ID._new_id(chain)
                for residue in chain.get_residues():
                    ID._new_id(residue)
                    for atom in residue.get_atoms():
                        ID._new_id(atom)
        return new

    @classmethod
    def from_biopython(cls, structure: "bio.Structure.Structure") -> "Structure":
        """
        Convert a BioPython Structure object to a Structure object.

        Parameters
        ----------
        structure : BioPython Structure object
            The structure to convert.

        Returns
        -------
        Structure
            The converted structure.
        """
        s = cls(structure.id)
        for model in structure.get_models():
            m = Model(model.id)
            for chain in model.get_chains():
                c = Chain(chain.id)
                for residue in chain.get_residues():
                    # ------------------------ NOTE -------------------------
                    # This is a little weird bugfix where I found
                    # that sometimes the segid was "     " instead of " ".
                    # This is a little hacky, but it works.
                    # It could be that there is a problem with the pdb module
                    # but that one has already seen enough modification so
                    # I don't want to tinker with it again...
                    # -------------------------------------------------------
                    segid = residue.segid
                    if len(segid) > 1:
                        segid = segid[0]
                    # -------------------------------------------------------
                    r = Residue(residue.resname, segid, residue.id[1])
                    for atom in residue.get_atoms():
                        a = Atom(
                            atom.id,
                            atom.coord,
                            atom.serial_number,
                            atom.bfactor,
                            atom.occupancy,
                            atom.fullname,
                            atom.element,
                            atom.altloc,
                            atom.pqr_charge,
                            atom.radius,
                        )
                        r.add(a)
                    c.add(r)
                m.add(c)
            s.add(m)
        return s

    def to_biopython(self) -> "bio.Structure.Structure":
        """
        Convert a Structure object to a pure BioPython Structure object.

        Returns
        -------
        bio.Structure.Structure
            The converted structure.
        """
        s = bio.Structure.Structure(self.id)
        for model in self.get_models():
            m = bio.Model.Model(model.id)
            for chain in model.get_chains():
                c = bio.Chain.Chain(chain.id)
                for residue in chain.get_residues():
                    r = bio.Residue.Residue(
                        (residue.resname, residue.serial_number, residue.segid),
                        residue.resname,
                        residue.segid,
                    )
                    for atom in residue.get_atoms():
                        a = bio.Atom.Atom(
                            atom.id,
                            atom.coord,
                            atom.bfactor,
                            atom.occupancy,
                            atom.altloc,
                            atom.fullname,
                            atom.serial_number,
                            atom.element,
                            atom.pqr_charge,
                            atom.radius,
                        )
                        r.add(a)
                    c.add(r)
                m.add(c)
            s.add(m)
        return s

    def move(self, vector):
        """
        Move the structure by a vector.

        Parameters
        ----------
        vector : ndarray
            The vector to move the structure by.
        """
        for model in self.get_models():
            model.move(vector)

    def link(self, model):
        """
        Softlink a model into this structure's child_list without touching the model's own parent references.
        """
        self.child_dict[model.get_id()] = model
        self.child_list.append(model)

    def __repr__(self):
        return f"Structure({self._id})"

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __le__(self, other):
        return self.id <= other.id

    def __ge__(self, other):
        return self.id >= other.id


class Bond:
    """
    A class representing a bond between two atoms.

    Attributes
    ----------
    atom1 : Atom
        The first atom in the bond.
    atom2 : Atom
        The second atom in the bond.
    """

    __linkers = {0: "<none>", 1: "--", 2: "==", 3: "#"}

    __slots__ = ("atom1", "atom2", "order")

    def __init__(self, *atoms) -> None:
        if len(atoms) == 1:
            self = Bond(*atoms[0])
        elif len(atoms) == 2:
            self.atom1 = atoms[0]
            self.atom2 = atoms[1]
            self.order = 1
        elif len(atoms) == 3:
            self.atom1 = atoms[0]
            self.atom2 = atoms[1]
            self.order = atoms[2]
        else:
            raise ValueError("Bond must be initialized with one tuple or two atoms")

    @property
    def length(self) -> float:
        return self.compute_length()

    def invert(self):
        """
        Invert the bond, i.e. swap the two atoms.
        """
        self.atom1, self.atom2 = self.atom2, self.atom1
        return self

    def single(self):
        """
        Make the bond a single bond.
        """
        self.order = 1
        return self

    def double(self):
        """
        Make the bond a double bond.
        """
        self.order = 2
        return self

    def triple(self):
        """
        Make the bond a triple bond.
        """
        self.order = 3
        return self

    def is_single(self) -> bool:
        """
        Check if the bond is a single bond.

        Returns
        -------
        bool
            True if the bond is a single bond, False otherwise.
        """
        return self.order == 1

    def is_double(self) -> bool:
        """
        Check if the bond is a double bond.

        Returns
        -------
        bool
            True if the bond is a double bond, False otherwise.
        """
        return self.order == 2

    def is_triple(self) -> bool:
        """
        Check if the bond is a triple bond.

        Returns
        -------
        bool
            True if the bond is a triple bond, False otherwise.
        """
        return self.order == 3

    def compute_length(self) -> float:
        """
        Compute the bond length.

        Returns
        -------
        float
            The bond length.
        """
        return self.atom1 - self.atom2

    def to_vector(self) -> "np.ndarray":
        """
        Convert the bond to a vector.

        Returns
        -------
        ndarray
            The bond as a vector from atom1 to atom2.
        """
        vec = np.zeros(3)
        vec[0] = self.atom2.coord[0] - self.atom1.coord[0]
        vec[1] = self.atom2.coord[1] - self.atom1.coord[1]
        vec[2] = self.atom2.coord[2] - self.atom1.coord[2]
        return vec

    def to_list(self) -> list:
        """
        Convert the bond to a list of
        atom1, atom2, bond_order.

        Returns
        -------
        list
            The bond as a list.
        """
        return [self.atom1, self.atom2, self.order]

    def to_tuple(self) -> tuple:
        """
        Convert the bond to a tuple of
        atom1, atom2, bond_order.

        Returns
        -------
        tuple
            The bond as a tuple.
        """
        return (self.atom1, self.atom2, self.order)

    def __iter__(self):
        yield self.atom1
        yield self.atom2

    def __getitem__(self, idx):
        if idx == 0:
            return self.atom1
        elif idx == 1:
            return self.atom2
        else:
            raise IndexError("Bond only has two atoms")

    def __repr__(self) -> str:
        return f"Bond({self.atom1}, {self.atom2})"

    def __str__(self) -> str:
        return f"({self.atom1} {self.__linkers.get(self.order, '?') } {self.atom2})"

    def __eq__(self, other):
        a = self.atom1 == other[0] and self.atom2 == other[1]
        b = self.atom1 == other[1] and self.atom2 == other[0]
        return a or b

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.atom1) + hash(self.atom2)

    def __len__(self):
        return 2

    def __contains__(self, item):
        return item == self.atom1 or item == self.atom2


if __name__ == "__main__":
    s = bio.PDBParser().get_structure(
        "test", "/Users/noahhk/GIT/biobuild/support/examples/GLC.pdb"
    )

    _s = Structure.from_biopython(s)
    __s = _s.to_biopython()
    print(_s)
