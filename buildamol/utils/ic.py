"""
Abstract classes for storing force field data from CHARMM topology and parameter files
"""

import attr


@attr.s
class InternalCoordinates:
    """
    A representation of a single Internal Coordinate

    Parameters
    ----------
    atom1
        The first atom in the internal coordinate
    atom2
        The second atom in the internal coordinate
    atom3
        The third atom in the internal coordinate
    atom4
        The fourth atom in the internal coordinate
    bond_length_12
        The bond length between atoms 1 and 2
    bond_length_34
        The bond length between atoms 3 and 4
    bond_angle_123
        The bond angle between atoms 1, 2 and 3
    bond_angle_234
        The bond angle between atoms 2, 3 and 4
    dihedral
        The dihedral angle between atoms 1, 2, 3 and 4
    bond_length_13
        The bond length between atoms 1 and 3 (optional, for impropers)
    improper
        Whether the internal coordinate is an improper (optional)
    """

    # With this set of data, the position of atom 1 may be determined based on the
    # positions of atoms 2-4, and the position of atom 4 may be determined from the
    # positions of atoms 1-3, allowing the recursive generation of coordinates for
    # all atoms in the structure based on a three-atom seed.

    atom1 = attr.ib(type=str)
    atom2 = attr.ib(type=str)
    atom3 = attr.ib(type=str)
    atom4 = attr.ib(type=str)
    bond_length_12 = attr.ib(type=float, repr=False)
    bond_length_34 = attr.ib(type=float, repr=False)
    bond_angle_123 = attr.ib(type=float, repr=False)
    bond_angle_234 = attr.ib(type=float, repr=False)
    dihedral = attr.ib(type=float, repr=False)
    bond_length_13 = attr.ib(default=None, type=float, repr=False)
    improper = attr.ib(default=False, type=bool, repr=False)

    @classmethod
    def _from_dict(cls, _dict):
        """
        Make an InternalCoordinate from a JSON dictionary
        """
        new = cls(
            *_dict["atoms"],
            bond_length_12=_dict["bond_length_12"],
            bond_length_34=_dict["bond_length_34"],
            bond_angle_123=_dict["bond_angle_123"],
            bond_angle_234=_dict["bond_angle_234"],
            dihedral=_dict["dihedral"],
            bond_length_13=_dict["bond_length_13"],
            improper=_dict["improper"],
        )
        return new

    @classmethod
    def _from_xml(cls, xml):
        """
        Make an InternalCoordinate from an XML element
        """
        new = cls(
            xml.attributes["atom1"],
            xml.attributes["atom2"],
            xml.attributes["atom3"],
            xml.attributes["atom4"],
            bond_length_12=eval(xml.attributes["length12"]),
            bond_length_34=float(xml.attributes["length34"]),
            bond_angle_123=float(xml.attributes["angle123"]),
            bond_angle_234=float(xml.attributes["angle234"]),
            dihedral=float(xml.attributes["dihedral"]),
            bond_length_13=eval(xml.attributes["length13"]),
            improper=eval(xml.attributes["improper"]),
        )
        return new

    @classmethod
    def from_quartet(cls, quartet):
        """
        Make an InternalCoordinate from a Quartet
        """
        new = cls(
            quartet.atom1,
            quartet.atom2,
            quartet.atom3,
            quartet.atom4,
            bond_length_12=quartet.dist_12,
            bond_length_34=quartet.dist_34,
            bond_angle_123=quartet.angle_123,
            bond_angle_234=quartet.angle_234,
            dihedral=quartet.dihedral,
            bond_length_13=quartet.dist_13,
            improper=quartet.improper,
        )
        return new

    @property
    def angles(self):
        """
        Returns the bond angles that constitute the internal coordinate
        """
        return (self.bond_angle_123, self.bond_angle_234)

    @property
    def lengths(self):
        """
        Returns the bond lengths that constitute the internal coordinate
        """
        if self.improper:
            return (self.bond_length_13, self.bond_length_34)
        else:
            return (self.bond_length_12, self.bond_length_34)

    @property
    def atoms(self):
        """
        Returns the atoms that constitute the internal coordinate
        """
        return (self.atom1, self.atom2, self.atom3, self.atom4)

    @property
    def ids(self):
        """
        Returns the ids of the atoms that constitute the internal coordinate
        """
        if hasattr(self.atom1, "id"):
            return (self.atom1.id, self.atom2.id, self.atom3.id, self.atom4.id)
        else:
            return (self.atom1, self.atom2, self.atom3, self.atom4)

    @property
    def is_improper(self):
        """
        Returns True if the internal coordinate is improper
        """
        return self.improper

    @property
    def is_proper(self):
        """
        Returns True if the internal coordinate is proper
        """
        return not self.improper

    def get_reference_atoms(self, _src):
        """
        Get reference atoms with available coordinates from a source.
        Note, that this method does not check if reference coordinates
        are actually available, so make sure to only provide source objects
        that have cartesian coordinates contained.

        Parameters
        ----------
        _src : object
            The source object to get the reference atom from. This can be:
            - a list or tuple of Atoms
            - an object with a `get_atoms` method (e.g. a biopython Residue)

        Returns
        -------
        list
            A list of reference atoms with available coordinates
        """

        if hasattr(_src, "get_atoms"):
            _src = _src.get_atoms()
        elif not isinstance(_src, (list, tuple, set)):
            raise ValueError("Invalid source object")

        _src = {a.id: a for a in _src}
        _set_src = set(_src.keys())

        _intersection = set(self.ids).intersection(_set_src)
        if len(_intersection) == 0:
            return []

        _ref_atoms = []
        for _id in self.ids:
            if _id in _intersection:
                _ref_atoms.append(_src[_id])

        return _ref_atoms

    def __repr__(self):
        return f"InternalCoordinates({self.ids})"
