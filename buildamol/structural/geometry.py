"""
Geometric units for elemental molecular geometry
"""

import numpy as np

import buildamol.structural.base as base


class Geometry:
    """
    The base class for geometries
    """

    max_points = -1
    size = -1
    angle = -1

    def make_and_apply(
        self, atoms_to_make_from: list, atoms_to_apply_to: list, **kwargs
    ):
        """
        Automatically make and apply the geometry to a list of atom objects

        Parameters
        ----------
        atoms_to_make_from : list
            The list of atoms to make the geometry from.
            Note that there must not be more atoms than the geometry can generate coordinates for.
            Note that in each case the central atom must be the first atom in the list.
            Other atoms (planar or axial, where applicable) must come after. Be sure to pass
            a 'direction' argument to specify the inference direction if the points are ambiguous.
        atoms_to_apply_to : list
            The list of atoms to apply the geometry to.
            Note that there must not be more atoms than the geometry can generate coordinates for.
            Note that in each case the central atom must be the first atom in the list.
            Other atoms (planar or axial, where applicable) must come after.
        **kwargs
            Additional keyword arguments to pass to the coordinate generation function.
        """
        coords = self.make_coords(*atoms_to_make_from, **kwargs)
        if len(atoms_to_apply_to) > self.size:
            raise ValueError(
                f"Too many atoms for this geometry. {self.__class__.__name__} generates {self.size} coordinates but the input has {len(atoms_to_apply_to)} atoms."
            )
        for i, atom in enumerate(atoms_to_apply_to):
            atom.coord = coords[i]

    def apply(self, atoms: list, bonds: list = None, **kwargs):
        """
        Automatically apply coordinate generation to a list of atom objects

        Parameters
        ----------
        atoms : list
            The list of atoms to apply the geometry to.
            Note that there must not be more atoms than the geometry can generate coordinates for.
            Note that in each case the central atom must be the first atom in the list.
            Other atoms (planar or axial, where applicable) must come after. Be sure to pass
            a 'direction' argument to specify the inference direction if the points are ambiguous.
        bonds : list
            A list of tuples of the same atoms as in 'atoms' that are bonded. This can be used to adjust the bond lengths after the geometry is applied.

        **kwargs
            Additional keyword arguments to pass to the coordinate generation function.

        Returns
        -------
        list
            The list of atoms with the new coordinates applied.
        """
        if hasattr(atoms[0], "get_coord"):
            coord_getter = lambda x: x.get_coord()
        elif hasattr(atoms[0], "coord"):
            coord_getter = lambda x: x.coord
        else:
            raise ValueError(
                "Invalid input for 'apply'"
                + str(atoms)
                + ". If you have coordinates rather than a molecule or atoms use 'make_coords' instead."
            )

        old = [coord_getter(i) for i in atoms]

        if len(old) > self.size:
            raise ValueError(
                f"Too many atoms for this geometry. {self.__class__.__name__} generates {self.size} coordinates but the input has {len(old)} atoms."
            )

        new = self.make_coords(*old[: self.max_points], **kwargs)

        if hasattr(atoms[0], "set_coord"):
            coord_setter = lambda x, y: x.set_coord(y)
        else:
            coord_setter = lambda x, y: setattr(x, "coord", y)

        if bonds is not None:
            old_lengths = {
                b: np.linalg.norm(coord_getter(a) - coord_getter(b)) for a, b in bonds
            }

        for i, atom in enumerate(atoms):
            coord_setter(atom, new[i])

        if bonds is not None:
            for bond in bonds:
                base.adjust_distance(*bond, new_length=old_lengths[bond])

        return atoms

    def fill_hydrogens(self, *atoms, make_bonds: bool = True, **kwargs):
        """
        Fill the geometry with hydrogen atoms where empty coordinates are found

        Parameters
        ----------
        atoms : list
            The list of atoms to apply the geometry to.
            Note that there must not be more atoms than the geometry can generate coordinates for.
            Note that in each case the central atom must be the first atom in the list.
            Other atoms (planar or axial, where applicable) must come after. Be sure to pass a 'direction' argument to specify the inference direction if the points are ambiguous.
        make_bonds : bool
            If True, bond objects will also be returned connecting the central atom to the other atoms (provided or newly made).
        **kwargs
            Additional keyword arguments to pass to the coordinate generation function.

        Returns
        -------
        atoms : list
            The list of atoms with the new coordinates applied and hydrogens added where necessary.
        bonds : list
            A list of Bonds of the same atoms as in 'atoms' that are bonded. Or None if 'make_bonds' is False.
        """
        if len(atoms) > self.max_points:
            _atoms = atoms[: self.max_points]
        else:
            _atoms = atoms

        coords = self.make_coords(*(a.coord for a in _atoms), **kwargs)

        # # set the coordinates of the atoms
        # for i, atom in enumerate(atoms):
        #     if atom.coord is None:
        #         atom.coord = coords[i]

        # now fill the rest with hydrogens
        import buildamol.base_classes as core

        atoms = list(atoms)
        for coord in coords[len(atoms) :]:
            atom = core.Atom.new("H", coord=coord)
            atoms.append(atom)

        if make_bonds:
            bonds = [core.Bond(atoms[0], a) for a in atoms[1:]]
            return atoms, bonds
        else:
            return atoms, None

    def make_coords(*coords, **kwargs):
        """
        Make the coordinates of the geometry
        """
        raise NotImplementedError


class Tetrahedral(Geometry):
    """
    Tetrahedral geometry for Sp3 hybridized atoms
    (or atoms with 4 substituents)
    """

    max_points = 3
    size = 5
    angle = np.radians(120)
    dihedral = np.radians(109.5)

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, **kwargs):
        """
        Make the coordinates of a tetrahedron

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the tetrahedron
            The first atom is the center of the tetrahedron and the other 4 atoms follow.
        float : float
            The bond length to use for the tetrahedron.
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).

        Returns
        -------
        array-like
            The coordinates of the tetrahedron with center at the 0th index and the other 4 atoms following.
            Provided points always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            return self.make_coords_from_two(*coords, length=length)
        elif len(coords) == 3:
            return self.make_coords_from_three(*coords, length=length)
        # elif len(coords) == 4:
        #     return self.make_coords_from_four(*coords, length=length)
        else:
            raise ValueError("Invalid number of atoms")

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a tetrahedron

        Parameters
        ----------
        center : Atom or array-like
            The center of the tetrahedron
        float : float
            The bond length to use for the tetrahedron.
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the tetrahedron with center at the 0th index and the other 4 atoms following.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two(center, other, length=length)

    def make_coords_from_two(self, center, other, length: float = None):
        """
        Get the coordinates of a tetrahedron

        Parameters
        ----------
        center : Atom or array-like
            The center of the tetrahedron
        other : Atom or array-like
            The other atom to define one axis of the tetrahedron
        float : float
            The bond length to use for the tetrahedron.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the tetrahedron with center at the 0th index and other the 1st index, and remaing 3 atoms following.
            Provided points always precede the ones that were generated. E.g. other will be at index 1.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")

        axis = center - other
        axis /= np.linalg.norm(axis)
        a = np.abs(axis)
        if a[0] < 1e-6 and a[1] < 1e-6:
            axis2 = base.x_axis
        elif a[0] < 1e-6 and a[2] < 1e-6:
            axis2 = base.z_axis
        else:
            axis2 = base.y_axis

        coords = np.zeros((5, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-axis * length, self.dihedral, np.cross(axis, axis2))
            + center
        )

        _new_axis = center - coords[2]
        _new_axis /= np.linalg.norm(_new_axis)

        coords[3] = base.rotate_coords(-_new_axis * length, self.angle, axis) + center
        coords[4] = base.rotate_coords(-_new_axis * length, -self.angle, axis) + center

        return coords

    def make_coords_from_three(self, center, other1, other2, length: float = None):
        """
        Get the coordinates of a tetrahedron

        Parameters
        ----------
        center : Atom or array-like
            The center of the tetrahedron
        other1 : Atom or array-like
            The first atom to define one axis of the tetrahedron
        other2 : Atom or array-like
            The second atom to define one axis of the tetrahedron
        float : float
            The bond length to use for the tetrahedron.
            If not provided the distance between center and other1 is used.

        Returns
        -------
        array-like
            The coordinates of the tetrahedron with center at the 0th index and the other 4 atoms following.
            Provided points always precede the ones that were generated. E.g. other1 will be at index 1 and other2 at index 2.
        """
        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other1)
        length2 = np.linalg.norm(center - other2)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1

        axis1 = center - other1
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other2
        axis2 /= np.linalg.norm(axis2)

        coords = np.zeros((5, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = base.rotate_coords(-axis2 * length, self.angle, axis1) + center
        coords[4] = base.rotate_coords(-axis2 * length, 2 * self.angle, axis1) + center

        return coords

    # def make_coords_from_four(
    #     self, center, other1, other2, other3, length: float = None
    # ):
    #     """
    #     Get the coordinates of a tetrahedron

    #     Parameters
    #     ----------
    #     center : Atom or array-like
    #         The center of the tetrahedron
    #     other1 : Atom or array-like
    #         The first atom to define one axis of the tetrahedron
    #     other2 : Atom or array-like
    #         The second atom to define one axis of the tetrahedron
    #     other3 : Atom or array-like
    #         The third atom to define one axis of the tetrahedron
    #     float : float
    #         The bond length to use for the tetrahedron.
    #         If not provided the distance between center and other1 is used.

    #     Returns
    #     -------
    #     array-like
    #         The coordinates of the tetrahedron with center at the 0th index and the other 4 atoms following.
    #     """
    #     center = getattr(center, "coord", center)
    #     other1 = getattr(other1, "coord", other1)
    #     other2 = getattr(other2, "coord", other2)
    #     other3 = getattr(other3, "coord", other3)
    #     length1 = np.linalg.norm(center - other1)
    #     length2 = np.linalg.norm(center - other2)
    #     length3 = np.linalg.norm(center - other3)
    #     if length1 == 0 or length2 == 0 or length3 == 0:
    #         raise ValueError("The two atoms are at the same position")
    #     if length is None:
    #         length = length1

    #     axis1 = center - other1
    #     axis1 /= np.linalg.norm(axis1)
    #     # axis2 = center - other2
    #     # axis2 /= np.linalg.norm(axis2)
    #     axis3 = center - other3
    #     axis3 /= np.linalg.norm(axis3)

    #     coords = np.zeros((5, 3))
    #     coords[0] = center
    #     coords[1] = other1
    #     coords[2] = other2
    #     coords[3] = other3
    #     coords[4] = base.rotate_coords(-axis3 * length, self.angle, axis1) + center

    #     return coords


class TrigonalPlanar(Geometry):
    """
    Planar geometry for Sp2 hybridized atoms
    (or atoms with 3 substituents)
    """

    max_points = 3
    size = 4
    angle = np.radians(120)

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, **kwargs):
        """
        Make the coordinates of a planar triangle

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the planar triangle
            The first atom is the center of the planar triangle  and the other 3 atoms follow.
        float : float
            The bond length to use for the planar triangle .
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).

        Returns
        -------
        array-like
            The coordinates of the planar triangle  with center at the 0th index and the other 3 atoms following.
            Provided points always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            return self.make_coords_from_two(*coords, length=length)
        elif len(coords) == 3:
            return self.make_coords_from_three(*coords, length=length)
        else:
            raise ValueError("Invalid number of atoms")

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a planar triangle

        Parameters
        ----------
        center : Atom or array-like
            The center of the planar triangle
        float : float
            The bond length to use for the planar triangle .
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the planar triangle  with center at the 0th index and the other 3 atoms following.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two(center, other, length=length)

    def make_coords_from_two(self, center, other, length: float = None):
        """
        Get the coordinates of a planar triangle

        Parameters
        ----------
        center : Atom or array-like
            The center of the planar triangle
        other : Atom or array-like
            The other atom to define one axis of the planar triangle
        float : float
            The bond length to use for the planar triangle .
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the planar triangle  with center at the 0th index and other the 1st index, and remaing 2 atoms following.
            Provided points always precede the ones that were generated. E.g. other will be at index 1.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)

        coords = np.zeros((4, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-axis * length, self.angle, np.cross(axis, base.z_axis))
            + center
        )
        coords[3] = (
            base.rotate_coords(
                -axis * length, 2 * self.angle, np.cross(axis, base.z_axis)
            )
            + center
        )

        return coords

    def make_coords_from_three(self, center, other1, other2, length: float = None):
        """
        Get the coordinates of a planar triangle

        Parameters
        ----------
        center : Atom or array-like
            The center of the planar triangle
        other1 : Atom or array-like
            The first atom to define one axis of the planar triangle
        other2 : Atom or array-like
            The second atom to define one axis of the planar triangle
        float : float
            The bond length to use for the planar triangle .
            If not provided the distance between center and other1 is used.

        Returns
        -------
        array-like
            The coordinates of the planar triangle  with center at the 0th index and the other 3 atoms following.
            Provided points always precede the ones that were generated. E.g. other1 will be at index 1 and other2 at index 2.
        """
        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other1)
        length2 = np.linalg.norm(center - other2)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other1
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other2
        axis2 /= np.linalg.norm(axis2)

        coords = np.zeros((4, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2

        coords[3] = base.rotate_coords(-axis2 * length, np.pi, axis1) + center

        return coords


class SquarePlanar(Geometry):
    """
    Square Planar geometry (or partial octahedral without axial atoms)
    (has 4 substituents)
    """

    max_points = 3
    size = 5
    angle = np.pi / 2

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, **kwargs):
        """
        Make the coordinates of a square planar

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the square planar
            The first atom is the center of the square planar and two more atoms may follow.
        float : float
            The bond length to use for the square planar.
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).

        Returns
        -------
        array-like
            The coordinates of the square planar with center at the 0th index and the other 4 atoms following.
            Provided points always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            return self.make_coords_from_two(*coords, length=length)
        elif len(coords) == 3:
            return self.make_coords_from_three(*coords, length=length)

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a square planar

        Parameters
        ----------
        center : Atom or array-like
            The center of the square planar
        float : float
            The bond length to use for the square planar.
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the square planar with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two(center, other, length=length)

    def make_coords_from_two(self, center, other, length: float = None):
        """
        Make the coordinates of a square planar given the central node and another node


        Parameters
        ----------
        center : Atom or array-like
            The center of the square planar
        other : Atom or array-like
            The planar node
        float : float
            The bond length to use for the square planar.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the square planar with center at the 0th index and the other 4 atoms following.
            Provided points always precede the ones that were generated. E.g. other will be at the 1st index.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)
        perpendicular_axis = np.cross(axis, base.z_axis)

        coords = np.zeros((5, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-axis * length, self.angle, perpendicular_axis) + center
        )
        coords[3] = (
            base.rotate_coords(-axis * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-axis * length, 3 * self.angle, perpendicular_axis)
            + center
        )

        return coords

    def make_coords_from_three(self, center, other1, other2, length: float = None):
        """
        Make the coordinates of a square planar given the two other nodes

        Parameters
        ----------
        other1 : Atom or array-like
            The first planar node
        other2 : Atom or array-like
            The second planar node
        float : float
            The bond length to use for the square planar.
            If not provided the distance between other1 and other2 is used.

        Returns
        -------
        array-like
            The coordinates of the square planar with center at the 0th index and the other 4 atoms following.
            Provided points always precede the ones that were generated. E.g. other1 and other2 will be at the 1st and 2nd index.
        """

        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other1)
        length2 = np.linalg.norm(center - other2)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other1
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other2
        axis2 /= np.linalg.norm(axis2)

        if np.allclose(np.abs(axis1), np.abs(axis2)):
            a = np.abs(axis1)
            if a[0] < 1e-6 and a[1] < 1e-6:
                axis2 = base.y_axis
            elif a[0] < 1e-6 and a[2] < 1e-6:
                axis2 = base.z_axis
            else:
                axis2 = base.x_axis
        perpendicular_axis = np.cross(axis1, axis2)

        coords = np.zeros((5, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = (
            base.rotate_coords(-axis1 * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-axis1 * length, 3 * self.angle, perpendicular_axis)
            + center
        )

        return coords


class Linear(Geometry):
    """
    Linear geometry for Sp1 hybridized atoms
    (or atoms with 2 substituents)
    """

    max_points = 2
    size = 3
    angle = np.pi

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, **kwargs):
        """
        Make the coordinates of a line

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the line
            The first atom is the center of the line and the other 2 atoms follow.
        float : float
            The bond length to use for the line.
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).

        Returns
        -------
        array-like
            The coordinates of the line with center at the 0th index and the other 2 atoms following.
            Provided points always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            return self.make_coords_from_two(*coords, length=length)
        else:
            raise ValueError("Invalid number of atoms")

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a line

        Parameters
        ----------
        center : Atom or array-like
            The center of the line
        float : float
            The bond length to use for the line.
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the line with center at the 0th index and the other 2 atoms following.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two(center, other, length=length)

    def make_coords_from_two(self, center, other, length: float = None):
        """
        Get the coordinates of a line

        Parameters
        ----------
        center : Atom or array-like
            The center of the line
        other : Atom or array-like
            The other atom to define one axis of the line
        float : float
            The bond length to use for the line.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the line with center at the 0th index and other the 1st index, and
            the remaining atom following. Provided points always precede the ones that were generated.
            E.g. other will be at index 1.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)

        coords = np.zeros((3, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = center + axis * length

        return coords


class TrigonalBipyramidal(Geometry):
    """
    Trigonal Bipyramidal geometry for Sp3d hybridized atoms
    (or atoms with 5 substituents)
    """

    max_points = 3
    size = 6
    angle = np.radians(120)

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, direction: str = None):
        """
        Make the coordinates of a trigonal bipyramidal

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the trigonal bipyramidal
            The first atom is the center of the trigonal bipyramidal and two more atoms may follow.
        float : float
            The bond length to use for the trigonal bipyramidal.
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).
        direction : str
            If the points are ambiguous with respect to the axial and planar atoms, the direction can be specified.
            This can be either "axial", "planar", or "mixed".

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index
            and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Points that were provided always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            if direction == "axial":
                return self.make_coords_from_two_axial(*coords, length=length)
            elif direction == "planar":
                return self.make_coords_from_two_planar(*coords, length=length)
            else:
                raise ValueError(
                    "No direction specified or invalid direction for two points"
                )
        elif len(coords) == 3:
            if direction is None:
                direction = _infer_point_relations(coords, self.angle)
            if direction == "axial":
                return self.make_coords_from_two_axial(*coords[:2], length=length)
            elif direction == "planar":
                return self.make_coords_from_three_planar(*coords, length=length)
            elif direction == "mixed":
                return self.make_coords_from_three_mixed(*coords, length=length)
            else:
                raise ValueError("Invalid direction for three points")
        elif len(coords) > 3:
            raise ValueError(
                "Too many atoms for this geometry, provide at most 3 atoms"
            )

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a trigonal bipyramidal

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two_planar(center, other, length=length)

    def make_coords_from_two_axial(self, center, other, length: float = None):
        """
        Make the coordinates of a trigonal bipyramidal given the central node and an axial node
        (node that is in the linear axis with the central node)

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other : Atom or array-like
            The axial node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other will be at index 1.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)
        perpendicular_axis = np.cross(axis, base.z_axis)

        coords = np.zeros((6, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-perpendicular_axis * length, self.angle, axis) + center
        )
        coords[3] = (
            base.rotate_coords(-perpendicular_axis * length, 2 * self.angle, axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-perpendicular_axis * length, 3 * self.angle, axis)
            + center
        )
        coords[5] = center + axis * length

        return coords

    def make_coords_from_two_planar(self, center, other, length: float = None):
        """
        Make the coordinates of a trigonal bipyramidal given the central node and a planar node
        (node that is in the plane with the central node)

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other : Atom or array-like
            The planar node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other will be at index 1.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)

        a = np.abs(axis)
        if a[0] < 1e-6 and a[1] < 1e-6:
            axis2 = base.y_axis
        elif a[0] < 1e-6 and a[2] < 1e-6:
            axis2 = base.z_axis
        else:
            axis2 = base.x_axis

        perpendicular_axis = np.cross(axis2, base.z_axis)

        coords = np.zeros((6, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-axis * length, self.angle, perpendicular_axis) + center
        )
        coords[3] = (
            base.rotate_coords(-axis * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = center + perpendicular_axis * length
        coords[5] = center - perpendicular_axis * length

        return coords

    # def make_coords_from_three_axial(
    #     self,
    #     center,
    #     other1,
    #     other2,
    #     length: float = None,
    # ):
    #     """
    #     Make the coordinates of a trigonal bipyramidal given the axial nodes

    #     Parameters
    #     ----------
    #     other1 : Atom or array-like
    #         The first axial node
    #     other2 : Atom or array-like
    #         The second axial node
    #     float : float
    #         The bond length to use for the trigonal bipyramidal.
    #         If not provided the distance between other1 and other2 is used.

    #     Returns
    #     -------
    #     array-like
    #         The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
    #     """

    #     center = getattr(center, "coord", center)
    #     other1 = getattr(other1, "coord", other1)
    #     other2 = getattr(other2, "coord", other2)
    #     length1 = np.linalg.norm(center - other1)
    #     length2 = np.linalg.norm(center - other2)
    #     if length1 == 0 or length2 == 0:
    #         raise ValueError("The two atoms are at the same position")
    #     if length is None:
    #         length = length1
    #     axis1 = center - other1
    #     axis1 /= np.linalg.norm(axis1)
    #     axis2 = center - other2
    #     axis2 /= np.linalg.norm(axis2)

    #     if np.allclose(np.abs(axis1), np.abs(axis2)):
    #         a = np.abs(axis1)
    #         if a[0] < 1e-6 and a[1] < 1e-6:
    #             axis2 = base.y_axis
    #         elif a[0] < 1e-6 and a[2] < 1e-6:
    #             axis2 = base.z_axis
    #         else:
    #             axis2 = base.x_axis
    #     perpendicular_axis = np.cross(axis1, axis2)

    #     coords = np.zeros((6, 3))
    #     coords[0] = center
    #     coords[1] = (
    #         base.rotate_coords(-perpendicular_axis * length, self.angle, axis1) + center
    #     )
    #     coords[2] = (
    #         base.rotate_coords(-perpendicular_axis * length, 2 * self.angle, axis1)
    #         + center
    #     )
    #     coords[3] = (
    #         base.rotate_coords(-perpendicular_axis * length, 3 * self.angle, axis1)
    #         + center
    #     )
    #     coords[4] = other1
    #     coords[5] = other2

    #     return coords

    def make_coords_from_three_planar(
        self, center, other1, other2, length: float = None
    ):
        """
        Make the coordinates of a trigonal bipyramidal given the planar nodes

        Parameters
        ----------
        other1 : Atom or array-like
            The first planar node
        other2 : Atom or array-like
            The second planar node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between other1 and other2 is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other1 will be at index 1 and other2 at index 2.
        """

        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other1)
        length2 = np.linalg.norm(center - other2)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other1
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other2
        axis2 /= np.linalg.norm(axis2)

        if np.allclose(np.abs(axis1), np.abs(axis2)):
            a = np.abs(axis1)
            if a[0] < 1e-6 and a[1] < 1e-6:
                axis2 = base.y_axis
            elif a[0] < 1e-6 and a[2] < 1e-6:
                axis2 = base.z_axis
            else:
                axis2 = base.x_axis
        perpendicular_axis = np.cross(axis1, axis2)

        coords = np.zeros((6, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = (
            base.rotate_coords(-axis1 * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = center - perpendicular_axis * length
        coords[5] = center + perpendicular_axis * length

        return coords

    def make_coords_from_three_mixed(
        self, center, other1, other2, length: float = None
    ):
        """
        Make the coordinates of a trigonal bipyramidal given mixed nodes

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other1 : Atom or array-like
            The first node (planar)
        other2 : Atom or array-like
            The second node (axial)
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between other1 and other2 is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other1 will be at index 1 and other2 at index 2.
        """

        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other2)
        length2 = np.linalg.norm(center - other1)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other2
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other1
        axis2 /= np.linalg.norm(axis2)

        coords = np.zeros((6, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = base.rotate_coords(-axis2 * length, self.angle, axis1) + center
        coords[4] = base.rotate_coords(-axis2 * length, 2 * self.angle, axis1) + center
        coords[5] = center + axis1 * length

        return coords


def _infer_point_relations(points, planar_angle, mixed_angle=np.pi / 2):
    """
    Infer the directionality relation between points

    Parameters
    ----------
    points : array-like
        The coordinates of the points
    planar_angle : float
        The angle that defines the planar relation
    mixed_angle : float
        The angle that defines the planar-axial mixed relation

    Returns
    -------
    str
        The directionality relation between the points
    """
    if len(points) < 3:
        raise ValueError("At least 3 points are required")
    if len(points) == 3:

        if hasattr(points[0], "coord"):
            points = [p.coord for p in points]

        a = points[0] - points[1]
        b = points[0] - points[2]

        theta = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

        if abs(theta - mixed_angle) < 1e-3:
            return "mixed"

        elif abs(theta - planar_angle) < 1e-3:
            return "planar"

        elif abs(theta - np.pi) < 1e-3:
            return "axial"

        else:
            raise ValueError(
                "Cannot infer point relationships. Specify 'direction' manually."
            )

    raise ValueError("Invalid number of points")


class Octahedral(Geometry):
    """
    Octahedral geometry for Sp3d2 hybridized atoms
    (or atoms with 6 substituents)
    """

    max_points = 3
    size = 7
    angle = np.pi / 2

    def __init__(self, bond_length=1.2):
        self.bond_length = bond_length

    def make_coords(self, *coords, length: float = None, direction: str = None):
        """
        Make the coordinates of a trigonal bipyramidal

        Parameters
        ----------
        coords : array-like
            The coordinates of the atoms that define the trigonal bipyramidal
            The first atom is the center of the trigonal bipyramidal and two more atoms may follow.
        float : float
            The bond length to use for the trigonal bipyramidal.
            By default the bond length is either the default value or the distance between center and first other point (unless specified using this argument).
        direction : str
            If the points are ambiguous with respect to the axial and planar atoms, the direction can be specified.
            This can be either "axial", "planar", or "mixed".

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index
            and the other 6 atoms following, where first come the planar ones, then the axial ones.
            Points that were provided always precede the ones that were generated.
        """
        if len(coords) == 1:
            return self.make_coords_from_one(coords[0], length=length)
        elif len(coords) == 2:
            if direction == "axial":
                return self.make_coords_from_two_axial(*coords, length=length)
            elif direction == "planar":
                return self.make_coords_from_two_planar(*coords, length=length)
            else:
                raise ValueError(
                    "No direction specified or invalid direction for two points"
                )
        elif len(coords) == 3:
            if direction is None:
                direction = _infer_point_relations(coords, self.angle)
            if direction == "axial":
                return self.make_coords_from_two_axial(*coords[:2], length=length)
            elif direction == "planar":
                return self.make_coords_from_three_planar(*coords, length=length)
            elif direction == "mixed":
                return self.make_coords_from_three_mixed(*coords, length=length)
            else:
                raise ValueError("Invalid direction for three points")
        elif len(coords) > 3:
            raise ValueError(
                "Too many atoms for this geometry, provide at most 3 atoms"
            )

    def make_coords_from_one(self, center, length: float = None):
        """
        Get the coordinates of a trigonal bipyramidal

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the default bond length is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
        """
        length = length or self.bond_length
        center = getattr(center, "coord", center)
        other = center + base.x_axis * length
        return self.make_coords_from_two_planar(center, other, length=length)

    def make_coords_from_two_axial(self, center, other, length: float = None):
        """
        Make the coordinates of a trigonal bipyramidal given the central node and an axial node
        (node that is in the linear axis with the central node)

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other : Atom or array-like
            The axial node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 6 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. center and other will be at the 0th and 1st index.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)
        perpendicular_axis = np.cross(axis, base.z_axis)

        coords = np.zeros((7, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-perpendicular_axis * length, self.angle, axis) + center
        )
        coords[3] = (
            base.rotate_coords(-perpendicular_axis * length, 2 * self.angle, axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-perpendicular_axis * length, 3 * self.angle, axis)
            + center
        )
        coords[5] = (
            base.rotate_coords(-perpendicular_axis * length, 4 * self.angle, axis)
            + center
        )

        coords[6] = center + axis * length

        return coords

    def make_coords_from_two_planar(self, center, other, length: float = None):
        """
        Make the coordinates of a trigonal bipyramidal given the central node and a planar node
        (node that is in the plane with the central node)

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other : Atom or array-like
            The planar node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between center and other is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other will be at the 1st index.
        """
        center = getattr(center, "coord", center)
        other = getattr(other, "coord", other)
        if length is None:
            length = np.linalg.norm(center - other)
        if length == 0:
            raise ValueError("The two atoms are at the same position")
        axis = center - other
        axis /= np.linalg.norm(axis)
        perpendicular_axis = np.cross(axis, base.z_axis)

        coords = np.zeros((7, 3))
        coords[0] = center
        coords[1] = other
        coords[2] = (
            base.rotate_coords(-axis * length, self.angle, perpendicular_axis) + center
        )
        coords[3] = (
            base.rotate_coords(-axis * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-axis * length, 3 * self.angle, perpendicular_axis)
            + center
        )
        coords[5] = center + perpendicular_axis * length
        coords[6] = center - perpendicular_axis * length

        return coords

    # def make_coords_from_three_axial(
    #     self,
    #     center,
    #     other1,
    #     other2,
    #     length: float = None,
    # ):
    #     """
    #     Make the coordinates of a trigonal bipyramidal given the axial nodes

    #     Parameters
    #     ----------
    #     other1 : Atom or array-like
    #         The first axial node
    #     other2 : Atom or array-like
    #         The second axial node
    #     float : float
    #         The bond length to use for the trigonal bipyramidal.
    #         If not provided the distance between other1 and other2 is used.

    #     Returns
    #     -------
    #     array-like
    #         The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
    #     """

    #     center = getattr(center, "coord", center)
    #     other1 = getattr(other1, "coord", other1)
    #     other2 = getattr(other2, "coord", other2)
    #     length1 = np.linalg.norm(center - other1)
    #     length2 = np.linalg.norm(center - other2)
    #     if length1 == 0 or length2 == 0:
    #         raise ValueError("The two atoms are at the same position")
    #     if length is None:
    #         length = length1
    #     axis1 = center - other1
    #     axis1 /= np.linalg.norm(axis1)
    #     axis2 = center - other2
    #     axis2 /= np.linalg.norm(axis2)

    #     if np.allclose(np.abs(axis1), np.abs(axis2)):
    #         a = np.abs(axis1)
    #         if a[0] < 1e-6 and a[1] < 1e-6:
    #             axis2 = base.y_axis
    #         elif a[0] < 1e-6 and a[2] < 1e-6:
    #             axis2 = base.z_axis
    #         else:
    #             axis2 = base.x_axis
    #     perpendicular_axis = np.cross(axis1, axis2)

    #     coords = np.zeros((7, 3))
    #     coords[0] = center
    #     coords[1] = (
    #         base.rotate_coords(-perpendicular_axis * length, self.angle, axis1) + center
    #     )
    #     coords[2] = (
    #         base.rotate_coords(-perpendicular_axis * length, 2 * self.angle, axis1)
    #         + center
    #     )
    #     coords[3] = (
    #         base.rotate_coords(-perpendicular_axis * length, 3 * self.angle, axis1)
    #         + center
    #     )
    #     coords[4] = (
    #         base.rotate_coords(-perpendicular_axis * length, 4 * self.angle, axis1)
    #         + center
    #     )
    #     coords[5] = other1
    #     coords[6] = other2

    #     return coords

    def make_coords_from_three_planar(
        self, center, other1, other2, length: float = None
    ):
        """
        Make the coordinates of a trigonal bipyramidal given the planar nodes

        Parameters
        ----------
        other1 : Atom or array-like
            The first planar node
        other2 : Atom or array-like
            The second planar node
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between other1 and other2 is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other1 and other2 will be at the 1st and 2nd index.
        """

        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other1)
        length2 = np.linalg.norm(center - other2)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other1
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other2
        axis2 /= np.linalg.norm(axis2)

        if np.allclose(np.abs(axis1), np.abs(axis2)):
            a = np.abs(axis1)
            if a[0] < 1e-6 and a[1] < 1e-6:
                axis2 = base.y_axis
            elif a[0] < 1e-6 and a[2] < 1e-6:
                axis2 = base.z_axis
            else:
                axis2 = base.x_axis
        perpendicular_axis = np.cross(axis1, axis2)

        coords = np.zeros((7, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = (
            base.rotate_coords(-axis1 * length, 2 * self.angle, perpendicular_axis)
            + center
        )
        coords[4] = (
            base.rotate_coords(-axis1 * length, 3 * self.angle, perpendicular_axis)
            + center
        )
        coords[5] = center - perpendicular_axis * length
        coords[6] = center + perpendicular_axis * length

        return coords

    def make_coords_from_three_mixed(
        self, center, other1, other2, length: float = None
    ):
        """
        Make the coordinates of a trigonal bipyramidal given mixed nodes

        Parameters
        ----------
        center : Atom or array-like
            The center of the trigonal bipyramidal
        other1 : Atom or array-like
            The first node (planar)
        other2 : Atom or array-like
            The second node (axial)
        float : float
            The bond length to use for the trigonal bipyramidal.
            If not provided the distance between other1 and other2 is used.

        Returns
        -------
        array-like
            The coordinates of the trigonal bipyramidal with center at the 0th index and the other 5 atoms following, where first come the planar ones, then the axial ones.
            Provided points always precede the ones that were generated. E.g. other1 will be at index 1 and other2 will be at index 2.
        """

        center = getattr(center, "coord", center)
        other1 = getattr(other1, "coord", other1)
        other2 = getattr(other2, "coord", other2)
        length1 = np.linalg.norm(center - other2)
        length2 = np.linalg.norm(center - other1)
        if length1 == 0 or length2 == 0:
            raise ValueError("The two atoms are at the same position")
        if length is None:
            length = length1
        axis1 = center - other2
        axis1 /= np.linalg.norm(axis1)
        axis2 = center - other1
        axis2 /= np.linalg.norm(axis2)

        coords = np.zeros((7, 3))
        coords[0] = center
        coords[1] = other1
        coords[2] = other2
        coords[3] = base.rotate_coords(-axis2 * length, self.angle, axis1) + center
        coords[4] = base.rotate_coords(-axis2 * length, 2 * self.angle, axis1) + center
        coords[5] = base.rotate_coords(-axis2 * length, 3 * self.angle, axis1) + center
        coords[6] = center + axis1 * length

        return coords


tetrahedral = Tetrahedral()
"""
The default tetrahedral geometry
"""
trigonal_planar = TrigonalPlanar()
"""
The default trigonal planar geometry
"""

square_planar = SquarePlanar()
"""
The default square planar geometry
"""

linear = Linear()
"""
The default linear geometry
"""

trigonal_bipyramidal = TrigonalBipyramidal()
"""
The default trigonal bipyramidal geometry
"""

octahedral = Octahedral()
"""
The default octahedral geometry
"""


if __name__ == "__main__":
    import buildamol as bam

    base.origin = np.array([0, 0, 0], dtype=np.float64)

    bam.load_small_molecules()
    # mol = bam.get_compound("CH4")

    # mol.rename_atom("C", "A").rename_atom("HC1", "B").rename_atom(
    #     "HC2", "C"
    # ).rename_atom("HC3", "D").rename_atom("HC4", "E")

    # ics = infer.compute_internal_coordinates(mol.bonds)

    # c = mol.get_atom("A").coord
    # d = base.norm_vector(*mol.bonds[0])

    # # coords = compute_tetrahedron_coords(c, d)

    # tet = Sp3()
    # coords = tet.make_coords_from_two(mol.get_atom("A").coord, mol.get_atom("B").coord)

    # v = mol.draw()
    # v.draw_points(coords, colors=["orange"] * 5)

    # v.draw_points(
    #     tet.make_coords_from_one(mol.get_atom("A").coord),
    #     colors=["blue"] * 5,
    # )

    # out3 = tet.make_coords_from_three(
    #     mol.get_atom("A").coord, mol.get_atom("B").coord, mol.get_atom("C").coord
    # )
    # v.draw_points(
    #     out3,
    #     colors=["green"] * 5,
    # )
    # for i in range(1, 5):
    #     v.draw_vector(
    #         str(i),
    #         out3[0],
    #         out3[i],
    #         color="limegreen",
    #     )

    # v.show()

    # mol = bam.get_compound("formaldehyde")[0]

    # mol.rename_atom("C1", "A").rename_atom("O1", "B").rename_atom(
    #     "H11", "C"
    # ).rename_atom("H12", "D")

    # sp2 = Sp2()
    # coords = sp2.make_coords_from_two(mol.get_atom("A").coord, mol.get_atom("B").coord)

    # v = mol.draw()
    # v.draw_points(coords, colors=["orange"] * 4)

    # coords = sp2.make_coords_from_one(mol.get_atom("A").coord)
    # v.draw_points(
    #     coords,
    #     colors=["blue"] * 4,
    # )

    # coords = sp2.make_coords_from_three(
    #     mol.get_atom("A").coord,
    #     mol.get_atom("B").coord,
    #     mol.get_atom("C").coord,
    #     length=1.8,
    # )
    # v.draw_points(
    #     coords,
    #     colors=["green"] * 4,
    # )

    # v.show()

    # mol = bam.get_compound("C2H2")[0].autolabel()
    # mol.rename_atom("C1", "A").rename_atom("C2", "B").rename_atom(
    #     "H1", "C"
    # ).rename_atom("H2", "D")

    # sp1 = Linear()
    # coords = sp1.make_coords_from_two(mol.get_atom("A").coord, mol.get_atom("B").coord)

    # v = mol.draw()
    # v.draw_points(coords, colors=["orange"] * 3)

    # coords = sp1.make_coords_from_one(mol.get_atom("A").coord)
    # v.draw_points(
    #     coords,
    #     colors=["blue"] * 3,
    # )

    # v.show()
    bipl = TrigonalBipyramidal()

    v = bam.MoleculeViewer3D()
    # coords = bipl.make_coords_from_one(np.array([0, 0, 0], dtype=np.float64))

    # v.draw_points(coords, colors=["orange"] * 6)

    # coords = bipl.make_coords_from_two_axial(
    #     np.array([0, 0, 0], dtype=np.float64), np.array([0, 1, 0], dtype=np.float64)
    # )
    # v.draw_points(coords, colors=["blue"] * 6)

    # coords = bipl.make_coords_from_two_planar(
    #     np.array([0, 0, 0], dtype=np.float64), np.array([0, 1, 0], dtype=np.float64)
    # )
    # v.draw_points(coords, colors=["green"] * 6)

    # coords = bipl.make_coords_from_three_axial(
    #     np.array([0, 0, 0], dtype=np.float64),
    #     np.array([0, 1, 0], dtype=np.float64),
    #     np.array([0, -1, 0], dtype=np.float64),
    # )
    # v.draw_points(coords, colors=["red"] * 6)

    # p2 = base.rotate_coords(
    #     np.array([0, 1, 0], dtype=np.float64), np.radians(120), np.array([1, 0, 0])
    # )
    # coords = bipl.make_coords_from_three_planar(
    #     np.array([0, 0, 0], dtype=np.float64), np.array([0, 1, 0], dtype=np.float64), p2
    # )
    # v.draw_points(coords, colors=["purple"] * 6)

    # coords = bipl.make_coords_from_three_mixed(coords[0], coords[5], coords[1])
    # v.draw_points(coords, colors=["orange"] * 6)

    octa = TrigonalBipyramidal()
    coords = octa.make_coords(
        np.array([0, 0, 0], dtype=np.float64),
        # np.array([0.5, 0, 0], dtype=np.float64),
        # np.array([1, 0, 0], dtype=np.float64),
        # direction="axial",
    )

    v.draw_points(coords, colors=["orange"] * 7)

    coords = octa.make_coords(
        np.array([0, 0, 0], dtype=np.float64),
        np.array([0.5, 0, 0], dtype=np.float64),
        # np.array([1, 0, 0], dtype=np.float64),
        direction="planar",
    )

    v.draw_points(coords, colors=["blue"] * 7)

    coords = octa.make_coords(
        np.array([0, 0, 0], dtype=np.float64),
        np.array([1, 0, 0], dtype=np.float64),
        np.array([0, 1, 0], dtype=np.float64),
        direction="mixed",
    )

    v.draw_points(coords, colors=["green"] * 7)

    v.show()
    pass
