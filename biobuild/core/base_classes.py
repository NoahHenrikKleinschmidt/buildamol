"""
Wrapper classes for the Biopython classes to make them more convenient to use.
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
        return self.__id == other.__id

    def __ne__(self, other):
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
        self.icode = icode

    @property
    def id(self):
        return ("H_" + self.resname, self.icode, self.segid)

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

    def __repr__(self):
        return f"Residue({self.resname}, {self.icode}, chain={self.parent.id if self.parent else None})"

    def __lt__(self, other):
        return self.icode < other.icode

    def __gt__(self, other):
        return self.icode > other.icode

    def __le__(self, other):
        return self.icode <= other.icode

    def __ge__(self, other):
        return self.icode >= other.icode

    # def __eq__(self, other):
    #     return (
    #         self.icode == other.icode
    #         and self.resname == other.resname
    #         and self.parent == other.parent
    #     )

    # def __ne__(self, other):
    #     return (
    #         self.icode != other.icode
    #         or self.resname != other.resname
    #         or self.parent != other.parent
    #     )


class Chain(ID, bio.Chain.Chain):
    def __init__(self, id):
        ID.__init__(self)
        super(bio.Chain.Chain, self).__init__(id)

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

    @property
    def full_id(self):
        return self.id

    @full_id.setter
    def full_id(self, value):
        pass

    @classmethod
    def from_biopython(cls, biopython_structure: "bio.Structure.Structure"):
        s = cls(biopython_structure.id)
        rdx = 1
        adx = 1
        for model in biopython_structure.get_models():
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

    def to_biopython(self):
        s = bio.Structure.Structure(self.id)
        for model in self.get_models():
            m = bio.Model.Model(model.id)
            for chain in model.get_chains():
                c = bio.Chain.Chain(chain.id)
                for residue in chain.get_residues():
                    r = bio.Residue.Residue(
                        residue.resname, residue.segid, residue.icode
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
