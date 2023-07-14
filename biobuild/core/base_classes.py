"""
The base_classes are deriviatives of the original Biopython classes, but with
the change that they use a UUID4 as their identifier (full_id) instead of a hierarchical
tuple. This makes each object unique and allows for easy comparison where `a == b` is akin to `a is b`.
Consequently, the `__hash__` method is overwritten to use the UUID4 as the hash.

.. warning::
    
    Each class has its own `copy` method that returns a deep copy of the object with a new UUID4. So `a.copy() == a` is `False`, while a standard `deepcopy(a) == a` is `True` since the UUID4 will not have been updated automatically.

Converting to and from `biopython`
----------------------------------

Each biobuild class can be generated from a biopython class using the `from_biopython` class method. And each biobuild class has a `to_biopython` method that returns the pure-biopython equivalent.
It is important to note, that for most purposes, however, the biobuild classes should work fine as trop-in replacements for the original biopython classes. 

.. code-block:: python

    import Bio.PDB as bio
    from biobuild.core.base_classes import Atom

    bio_atom = bio.Atom("CA", (0, 0, 0))
    atom = Atom.from_biopython(bio_atom)

    assert atom == bio_atom # False since atom uses a UUID4 as its identifier
    assert atom.to_biopython() == bio_atom # True 

The conversion from and to biopython works hierarchically, so if an entire biopython structure is converted to biobuild
then all atoms, residues, chains and models will be converted to their biobuild equivalents.

.. code-block:: python

    import Bio.PDB as bio
    from biobuild.core.base_classes import Structure

    bio_structure = bio.PDBParser().get_structure("test", "test.pdb")
    structure = Structure.from_biopython(bio_structure)

    atoms = list(structure.get_atoms())
    bio_atoms = list(bio_structure.get_atoms())
    assert len(atoms) == len(bio_atoms) # True
    
"""
from copy import deepcopy
from uuid import uuid4
import Bio.PDB as bio
import periodictable as pt


class ID:
    def __init__(self):
        self.__id = str(uuid4())

    def copy(self):
        new = deepcopy(self)
        new._new_id()
        return new

    def get_id(self):
        return self.__id

    def has_id(self, id):
        return id in self.child_dict

    def _new_id(self):
        self.__id = str(uuid4())

    def __hash__(self):
        return hash(self.__id)

    def __eq__(self, other):
        if not isinstance(other, ID):
            return False
        return self.__id == other.__id

    def __ne__(self, other):
        if not isinstance(other, ID):
            return True
        return self.__id != other.__id


class Atom(ID, bio.Atom.Atom):
    def __init__(
        self,
        id: str,
        coord: "ndarray",
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
            coord,
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
    def __init__(self, resname, segid, icode):
        ID.__init__(self)
        bio.Residue.Residue.__init__(
            self, ("H_" + resname, icode, segid), resname, segid
        )
        self.level = "R"
        self.serial_number = icode

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
    def __init__(self, id):
        ID.__init__(self)
        super(bio.Chain.Chain, self).__init__(id)
        self.level = "C"

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
    def __init__(self, id):
        ID.__init__(self)
        super(bio.Model.Model, self).__init__(id)
        self.level = "M"

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

    def __repr__(self):
        return f"Model({self._id})"

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __le__(self, other):
        return self.id <= other.id

    def __ge__(self, other):
        return self.id >= other.id


class Structure(ID, bio.Structure.Structure):
    def __init__(self, id):
        ID.__init__(self)
        super(bio.Structure.Structure, self).__init__(id)
        self.level = "S"

    @property
    def full_id(self):
        return self.id

    @full_id.setter
    def full_id(self, value):
        pass

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
        rdx = 1
        adx = 1
        for model in structure.get_models():
            m = Model(model.id)
            for chain in model.get_chains():
                c = Chain(chain.id)
                for residue in chain.get_residues():
                    r = Residue(residue.resname, residue.segid, rdx)
                    rdx += 1
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


if __name__ == "__main__":
    s = bio.PDBParser().get_structure(
        "test", "/Users/noahhk/GIT/biobuild/support/examples/GLC.pdb"
    )

    _s = Structure.from_biopython(s)
    __s = _s.to_biopython()
    print(_s)
