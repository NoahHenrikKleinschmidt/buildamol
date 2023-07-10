"""
Linkage definitions
===================

A linkage is a connection between two _molecules_. At its core each linkage simply defines two atoms that should be connected,
and what atoms to remove in the process. It is a "pseudo" chemical reaction, so to speak. 

Building on the CHARMM force field, biobuild distinguishes two kinds of linkages: patches and recipies.

A **patch** is a linkage that can be applied purely geometrically and does not require numeric optimization. This is because a 
patch includes geometric data in form of _internal coordinates_ of the atoms in the immediate vicinity of the newly formed bond.
Using this data, biobuild is able to attach molecule to one another through simple matrix transformations. Conesquently, patches
are the most efficient way to connect molecules and are preferable to **recipes** - the other type of linkage.

A **recipe** on the other hand is a linkage that requires numeric optimization. This is because a recipe does not include any
geometric data, but only the atoms that should be connected. The numeric optimization is then used to find the optimal (or at least suitable) conformation.
This is useful for most users who wish to define their own linkage types, but who will likely not wish to painstakingly define the detailed geometry of angles and dihedrals of the atom neighborhood.

The distinction between patches and recipies is purely nominal, as both are represented by the `Linkage` class. However, there are functional
wrappers available to create either a patch or recipe, respectively, which require different arguments (to make sure they are not forgotten and to make the code more readable).

.. code-block:: python

    from biobuild import recipe

    # Create a custom recipe
    my_link = recipe(
        atom1 = "C1",
        atom2 = "O4",
        delete_in_target = ["O1", "HO1"],
        delete_in_source = ["HO4"],
        id = "my_link"
    )


Pre-defined patches
-------------------
biobuild comes with a number of pre-defined patches from the CHARMM force field. These can be accessed through the `resources` module:

.. code-block:: python

    from biobuild import resources

    # Get a list of all pre-defined patches
    patches = resources.available_patches()

    # Check for a specific patch
    resources.has_patch("some_patch")

    # Get a specific patch
    my_patch = resources.get_patch("some_patch")

A custom linkage can be added to the list of pre-defined patches by using the `add_patch` function:

.. code-block:: python

    # add the above defined my_link to the list of pre-defined patches
    resources.add_patch(my_link)

.. note::

    Despite the use of "patch" in the function nomenclature, there is no difference between a patch and a recipe in terms of
    how they are used. Patches and Recipies are represented by the same data class and thus behave identically.
    Hence, there are also functional wrappers with the "linkage" available
    that can be used instead (if a user feels more comfortable with this) - they perform the same function.

    .. code-block:: python

        resources.add_linkage(my_link)
        # performs the same as 
        resources.add_patch(my_link)

        # check for a specific linkage
        resources.has_linkage("my_link")
        # performs the same as
        resources.has_patch("my_link")

        # etc.


Pre-defined patches can be accessed directly by their `id` and need not be obtained first through the `resources` module. They can be directly passed
to the ``Molecule``'s ``attach`` method or any other function that requires a linkage:

.. code-block:: python

    import biobuild as bb 

    mol1 = bb.read_pdb("my_molecule.pdb")
    mol2 = bb.read_pdb("my_other_molecule.pdb")

    # Attach mol2 to mol1 using the pre-defined patch "some_patch"
    mol1.attach(mol2, "some_patch")

    # works the same as doing
    some_patch = bb.get_patch("some_patch")
    mol1.attach(mol2, some_patch)
    

"""

import biobuild.utils as utils
import biobuild.structural.neighbors as neighbors


def patch(
    atom1,
    atom2,
    delete_in_target,
    delete_in_source,
    internal_coordinates: dict,
    id: str = None,
    description: str = None,
) -> "Linkage":
    """
    Make a new `Linkage` instance that describes a "patch" between two molecules.
    A patch is a linkage that can be applied purely geometrically and does not require numeric optimization.
    As such, it requires the internal coordinates of the atoms in the immediate vicinity of the newly formed bond.

    Parameters
    ----------
    atom1 : str or tuple of str
        The atom in the first (target) molecule to connect.
    atom2 : str or tuple of str
        The atom in the second (source) molecule to connect.
    delete_in_target : str or tuple of str
        The atom(s) in the first molecule to delete.
    delete_in_source : str or tuple of str
        The atom(s) in the second molecule to delete.
    internal_coordinates : dict, optional
        The internal coordinates of the atoms in the immediate vicinity of the newly formed bond.
        If provided, the link can be applied purely geometrically and will not require numeric optimization.
        If provided, this must be a dictionary where keys are tuples of four atoms ids and values tuples containing (in order):

            - the bond length between the first and second atom (first and third in case of an improper)
            - the bond length between the third and fourth atom
            - the bond angle between the first, second and third atom
            - the bond angle between the second, third and fourth atom
            - the dihedral angle between the first, second, third and fourth atom
            - True if the internal coordinate is improper, False otherwise

    id : str, optional
        The id of the linkage.
    description : str, optional
        A description of the linkage.

    Returns
    -------
    Linkage
        The new linkage.
    """
    return linkage(
        atom1,
        atom2,
        delete_in_target=delete_in_target,
        delete_in_source=delete_in_source,
        internal_coordinates=internal_coordinates,
        id=id,
        description=description,
    )


def recipe(
    atom1,
    atom2,
    delete_in_target,
    delete_in_source,
    id: str = None,
    description: str = None,
) -> "Linkage":
    """
    Make a new `Linkage` instance that describes a "recipe" to connect two molecules.
    A recipe is a linkage that can be applied numerically and requires numeric optimization as it does not
    have the internal coordinates of the atoms in the immediate vicinity of the newly formed bond.

    Parameters
    ----------
    atom1 : str or tuple of str
        The atom in the first (target) molecule to connect.
    atom2 : str or tuple of str
        The atom in the second (source) molecule to connect.
    delete_in_target : str or tuple of str
        The atom(s) in the first molecule to delete.
    delete_in_source : str or tuple of str
        The atom(s) in the second molecule to delete.
    id : str, optional
        The id of the linkage.
    description : str, optional
        A description of the linkage.

    Returns
    -------
    Linkage
        The new linkage.
    """
    return linkage(
        atom1,
        atom2,
        delete_in_target=delete_in_target,
        delete_in_source=delete_in_source,
        id=id,
        description=description,
    )


def linkage(
    atom1,
    atom2,
    delete_in_target,
    delete_in_source,
    internal_coordinates: dict = None,
    id: str = None,
    description: str = None,
) -> "Linkage":
    """
    Make a new `Linkage` instance to connect two molecules together.

    Parameters
    ----------
    atom1 : str or tuple of str
        The atom in the first (target) molecule to connect.
    atom2 : str or tuple of str
        The atom in the second (source) molecule to connect.
    delete_in_target : str or tuple of str, optional
        The atom(s) in the first molecule to delete.
    delete_in_source : str or tuple of str, optional
        The atom(s) in the second molecule to delete.
    internal_coordinates : dict, optional
        The internal coordinates of the atoms in the immediate vicinity of the newly formed bond.
        If provided, the link can be applied purely geometrically and will not require numeric optimization.
        If provided, this must be a dictionary where keys are tuples of four atoms ids and values tuples containing (in order):

            - the bond length between the first and second atom (first and third in case of an improper)
            - the bond length between the third and fourth atom
            - the bond angle between the first, second and third atom
            - the bond angle between the second, third and fourth atom
            - the dihedral angle between the first, second, third and fourth atom
            - True if the internal coordinate is improper, False otherwise

    id : str, optional
        The ID of the linkage.
    description : str, optional
        A description of the linkage.

    Returns
    -------
    Linkage
        The new linkage instance.
    """
    # make a new linkage
    new_linkage = Linkage(id=id, description=description)

    # add the bond
    new_linkage.add_bond(utils.abstract.AbstractBond(atom1, atom2))

    # add the atoms to delete
    if delete_in_target is not None:
        for i in delete_in_target:
            new_linkage.add_delete(i, "target")
    if delete_in_source is not None:
        for i in delete_in_source:
            new_linkage.add_delete(i, "source")

    # add the internal coordinates
    if internal_coordinates is not None:
        if not isinstance(internal_coordinates, dict):
            raise TypeError(
                "The internal coordinates must be provided as a dictionary."
            )
        for ic in _dict_to_ics(internal_coordinates):
            new_linkage.add_internal_coordinates(ic)

    # return the linkage
    return new_linkage


class Linkage(utils.abstract.AbstractEntity_with_IC):
    """
    Using the `Linkage` class, a template reaction instruction is stored for attaching molecules to one another.

    Parameters
    ----------
    id : str, optional
        The ID of the linkage.
    description : str, optional
        An additional description of the linkage.

    Attributes
    ----------
    id : str
        The ID of the linkage.
    bond : tuple of str
        The bond to form between the two molecules.
    internal_coordinates : list of InternalCoordinate
        The internal coordinates of the atoms in the immediate vicinity of the newly formed bond.
    deletes : tuple of list of str
        The atom IDs to delete
        in a tuple of lists where the first list
        contains the atom IDs to delete from the
        first structure (target) and the second one from the second structure (source)
    atoms : list of str
        The atom IDs of the atoms in the linkage.
    """

    def __init__(self, id=None, description: str = None) -> None:
        super().__init__(id)
        self._delete_ids = []
        self.description = description

    @property
    def bond(self) -> tuple:
        """
        The bond to form between the two molecules.
        """
        if len(self.bonds) == 0:
            return None
        return self.bonds[0]

    @property
    def _ref_atoms(self) -> tuple:
        """
        Reference atoms with 1,2 prefix for the patcher
        """
        if not self.bond:
            return None, None
        a = f"1{self.bond[0]}" if not self.bond[0].startswith("1") else self.bond[0]
        b = f"2{self.bond[1]}" if not self.bond[1].startswith("2") else self.bond[1]
        return a, b

    @classmethod
    def from_json(cls, filename: str):
        """
        Make a new `Linkage` instance from a JSON file.

        Parameters
        ----------
        filename : str
            The JSON filename.
        """
        _dict = utils.json.read(filename)
        return cls._from_dict(_dict)

    @classmethod
    def _from_dict(cls, _dict):
        """
        Make a new `Linkage` instance from a JSON dictionary.

        Parameters
        ----------
        _dict : dict
            The JSON dictionary.
        """
        new = cls(id=_dict["id"])
        new.add_bond(
            utils.abstract.AbstractBond(
                _dict["bond"]["target"], _dict["bond"]["source"]
            ),
        )
        for i in _dict["to_delete"]["target"]:
            new.add_delete(i, "target")
        for i in _dict["to_delete"]["source"]:
            new.add_delete(i, "source")
        for i in _dict["ics"]:
            new.add_internal_coordinates(
                utils.ic.InternalCoordinates._from_dict(i),
            )
        has_two_residues = False
        for i in new.internal_coordinates:
            if any([j.startswith("2") for j in i.atoms]):
                has_two_residues = True
        if not has_two_residues:
            raise ValueError(
                "The linkage contains only internal coordinates for one residue. It must contain internal coordinates spanning both residues!"
            )
        return new

    @property
    def deletes(self):
        """
        Returns the atom IDs to delete
        in a tuple of lists where the first list
        contains the atom IDs to delete from the
        first structure (target) and the second one from the second structure (source)
        """
        deletes = (
            [i[1:] for i in self._delete_ids if i[0] == "1"],
            [i[1:] for i in self._delete_ids if i[0] == "2"],
        )
        return deletes

    def add_delete(self, id, _from: str = None):
        """
        Add an atom ID to delete

        Parameters
        ----------
        id : str
            The atom ID to delete.
        _from : str, optional
            The structure from which to delete the atom.
            Can be either "source" or "target". If not provided,
            the structure is inferred from the atom ID, in which case either `1` (target) or `2` (source) must be the first character of the ID.
        """
        if _from is None:
            if not id[0] in ["1", "2"]:
                raise ValueError(
                    "The atom ID must start with either 1 or 2 to indicate from which structure it is to be deleted."
                )
        else:
            if _from == "source":
                if isinstance(id, str):
                    id = "2" + id
                else:
                    id = ("2", *id)
            elif _from == "target":
                if isinstance(id, str):
                    id = "1" + id
                else:
                    id = ("1", *id)
            else:
                raise ValueError(
                    "The _from argument must be either 'source' or 'target'."
                )
        self._delete_ids.append(id)

    add_id_to_delete = add_delete

    def add_internal_coordinates(self, ic):
        if isinstance(ic, dict):
            ic = utils.ic.InternalCoordinates._from_dict(ic)
        elif isinstance(ic, neighbors.Quartet):
            ic = utils.ic.InternalCoordinates.from_quartet(ic)
        elif isinstance(ic, (tuple, list)) and len(ic) == 11:
            ic = utils.ic.InternalCoordinates(*ic)
        elif not isinstance(ic, utils.ic.InternalCoordinates):
            raise TypeError(
                "The internal coordinate must be an instance of the InternalCoordinates class."
            )

        # now we need to vet that the internal coordinates span both molecules (by spanning different residues)
        if isinstance(ic.atom1, str):
            _residues = list(set(a[0] for a in ic.atoms))
            use_bare_strings = True
        else:
            _residues = list(set(a.get_parent() for a in ic.atoms))
            use_bare_strings = False

        if use_bare_strings:
            _residues.sort()
        else:
            _residues.sort(key=lambda x: x.id[1])

        _residues_dict = {_residues[0]: "1"}
        if len(_residues) == 2:
            _residues_dict[_residues[1]] = "2"

        if not use_bare_strings:
            prefix = (
                lambda x: _residues_dict[x.get_parent()] + x.id
                if not x.id[0] in ("1", "2")
                else x.id
            )

            ic.atom1 = prefix(ic.atom1)
            ic.atom2 = prefix(ic.atom2)
            ic.atom3 = prefix(ic.atom3)
            ic.atom4 = prefix(ic.atom4)

        return super().add_internal_coordinates(ic)

    def to_json(self, filename: str):
        """
        Write the `Linkage` instance to a JSON file.

        Parameters
        ----------
        filename : str
            The JSON filename.
        """
        utils.json.write_linkage(self, filename)

    def __getitem__(self, index):
        if len(self.bonds) == 0:
            raise IndexError("The linkage does not contain a bond.")
        return self.bonds[0][index]

    def __iter__(self):
        if len(self.bonds) == 0:
            raise StopIteration
        return iter(self.bonds[0])


def _dict_to_ics(_dict):
    """
    Convert a dictionary of internal coordinates to a list of `InternalCoordinate` instances.
    """
    ics = []
    for key, value in _dict.items():
        if len(key) == 4 and len(value) == 6:
            improper = value[-1]

            ic = utils.ic.InternalCoordinates(
                *key,
                bond_length_12=value[0] if not improper else None,
                bond_length_13=value[0] if improper else None,
                bond_length_34=value[1],
                bond_angle_123=value[2],
                bond_angle_234=value[3],
                dihedral=value[4],
                improper=improper,
            )
        elif len(key) != 4:
            raise ValueError(
                "The internal coordinate must be provided as a tuple of four atom IDs."
            )
        elif len(value) != 6:
            raise ValueError(
                "The internal coordinate must be provided as a tuple of six values."
            )
        ics.append(ic)

    return ics


if __name__ == "__main__":
    link = linkage(
        "C1",
        "O4",
        ["H1"],
        ["HO4"],
        internal_coordinates={
            ("1C1", "1C2", "2O4", "2C4"): [1.1, 1.2, 1.3, 1.4, 1.5, False]
        },
    )
    link.to_json("link.json")
    link2 = Linkage.from_json("link.json")
    pass
