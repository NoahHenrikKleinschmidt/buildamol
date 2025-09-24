import buildamol.base_classes as base_classes
import buildamol.core as core


class Reactivity:
    """
    Define a reactivity pattern for a molecule or functional group by specifying functions that identify nucleophilic and electrophilic atoms, as well as optional functions to determine which atoms to delete during the reaction.
    The class also allows for steric constraints to be applied when selecting reactive sites.

    Parameters
    ----------
    nucleophile_linker : callable, optional
        A function that takes a `Molecule` and returns an `Atom` object representing the nucleophilic site (i.e. the atom in the source molecule that will form a bond with the target molecule).
    electrophile_linker : callable, optional
        A function that takes a `Molecule` and returns a single or list of multiple `Atom` objects representing the electrophilic sites (i.e. the atoms in the target molecule that can form a bond with the source molecule).
    nucleophile_deleter : callable, optional
        A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the source molecule when the nucleophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the nucleophilic atom.
    electrophile_deleter : callable, optional
        A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the target molecule when the electrophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the electrophilic atom.
    """

    def __init__(self, *args, **kwargs):
        self._nucleophile_linker = None
        self._electrophile_linker = None
        self._nucleophile_deleter = None
        self._electrophile_deleter = None
        self._steric_constraints_func = None
        self._steric_distance = None
        self._steric_max_neighbors = None
        self._steric_n_target_sites = None
        self._serves_target = False

        if args or kwargs:
            self.set_reactivity(*args, **kwargs)
        else:
            self.set_reactivity(
                nucleophile_linker=getattr(self, "nucleophile_linker", None),
                electrophile_linker=getattr(self, "electrophile_linker", None),
                nucleophile_deleter=getattr(self, "nucleophile_deleter", None),
                electrophile_deleter=getattr(self, "electrophile_deleter", None),
            )
        self.set_steric_constraints()

    def set_steric_constraints(
        self,
        distance: float = 4.0,
        max_neighbors: int = None,
        n_target_sites: int = "all",
        func: callable = None,
    ):
        """
        Specify steric constraints to be applied when selecting reactive sites. Note that these are shared between the nucleophile and electrophile.
        If different constraints are needed, directly specify them in the linker/deleter functions.

        Parameters
        ----------
        distance : float, optional
            The distance (in angstroms) within which to count neighboring atoms for steric hindrance. Default is 4.0 Å.
        max_neighbors : int, optional
            The maximum number of neighboring atoms allowed within the specified distance for a site to be considered accessible. If None, no limit is applied. Default is None.
        n_target_sites : int or "all", optional
            The number of the most accessible sites to return after applying steric constraints. If "all", all sites that meet the steric criteria are returned. Default is "all".
        func : callable, optional
            A custom function that takes a `Molecule` and a list of `Atom` objects (in that order) and returns a filtered list of `Atom` objects based on custom steric criteria. If provided, this function is applied *before* the default steric constraints.
            But it will **not** override the default steric constraints, which will still be applied after this function.
        """
        if func is not None:
            self._steric_constraints_func = func
        self._steric_distance = distance
        self._steric_max_neighbors = max_neighbors
        self._steric_n_target_sites = n_target_sites
        return self

    def set_reactivity(
        self,
        nucleophile_linker: callable = None,
        electrophile_linker: callable = None,
        nucleophile_deleter: callable = None,
        electrophile_deleter: callable = None,
    ):
        """
        Set or update the reactivity functions for the nucleophile and electrophile. This is an in-place operation. Use `with_reactivity` to create a new instance with modified reactivity.

        Parameters
        ----------
        nucleophile_linker : callable, optional
            A function that takes a `Molecule` and returns an `Atom` object representing the nucleophilic site (i.e. the atom in the source molecule that will form a bond with the target molecule).
        electrophile_linker : callable, optional
            A function that takes a `Molecule` and returns a single or list of multiple `Atom` objects representing the electrophilic sites (i.e. the atoms in the target molecule that can form a bond with the source molecule).
        nucleophile_deleter : callable, optional
            A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the source molecule when the nucleophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the nucleophilic atom.
        electrophile_deleter : callable, optional
            A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the target molecule when the electrophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the electrophilic atom.
        """
        if nucleophile_linker is not None:
            self._nucleophile_linker = nucleophile_linker
        if electrophile_linker is not None:
            self._electrophile_linker = electrophile_linker
        if nucleophile_deleter is not None:
            self._nucleophile_deleter = nucleophile_deleter
        if electrophile_deleter is not None:
            self._electrophile_deleter = electrophile_deleter

        return self

    def with_reactivity(
        self,
        nucleophile_linker: callable = None,
        electrophile_linker: callable = None,
        nucleophile_deleter: callable = None,
        electrophile_deleter: callable = None,
    ):
        """
        Create a new instance of the class with modified reactivity functions. This does not modify the original instance.

        Parameters
        ----------
        nucleophile_linker : callable, optional
            A function that takes a `Molecule` and returns an `Atom` object representing the nucleophilic site (i.e. the atom in the source molecule that will form a bond with the target molecule).
        electrophile_linker : callable, optional
            A function that takes a `Molecule` and returns a single or list of multiple `Atom` objects representing the electrophilic sites (i.e. the atoms in the target molecule that can form a bond with the source molecule).
        nucleophile_deleter : callable, optional
            A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the source molecule when the nucleophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the nucleophilic atom.
        electrophile_deleter : callable, optional
            A function that takes an `Atom` and its parent `Molecule` (in that order) and returns an `Atom` or list thereof representing the atoms to be deleted from the target molecule when the electrophilic atom forms a bond. If not provided, the default behavior is to delete a hydrogen atom bonded to the electrophilic atom.
        """
        new = self.__class__()
        new.set_reactivity(
            nucleophile_linker=nucleophile_linker or self._nucleophile_linker,
            electrophile_linker=electrophile_linker or self._electrophile_linker,
            nucleophile_deleter=nucleophile_deleter or self._nucleophile_deleter,
            electrophile_deleter=electrophile_deleter or self._electrophile_deleter,
        )
        return new

    def find_atoms(
        self,
        mol: core.Molecule,
        role: str,
        serves_target: bool,
    ):
        """
        Find reactive atoms in the given molecule based on the specified reactivity pattern.

        Parameters
        ----------
        mol : Molecule
            The molecule in which to find reactive atoms.
        serves_target : bool
            If True, the linker function will return all identified nucleophilic sites. If False, it will return only the most accessible site based on steric constraints.
        role : str, optional
            Specify whether to find 'nucleophile' or 'electrophile' atoms.

        Returns
        -------
        linkers : list of Atom or Atom
            The identified reactive atoms based on the specified role and steric constraints.
        deleters : list of Atom or Atom
            The atoms to be deleted during the reaction based on the specified role.
        """
        if role == "nucleophile":
            linker, deleter = self.as_nucleophile(serves_target)
        elif role == "electrophile":
            linker, deleter = self.as_electrophile(serves_target)
        else:
            raise ValueError("Role must be either 'nucleophile' or 'electrophile'")
        linkers = linker(mol)
        if serves_target:
            deleters = [deleter(a, mol) for a in linkers]
        else:
            deleters = deleter(linkers, mol)
        return linkers, deleters

    def find_nucleophilic_atoms(self, mol: core.Molecule, serves_target: bool):
        """
        Find nucleophilic atoms in the given molecule based on the specified reactivity pattern.

        Parameters
        ----------
        mol : core.Molecule
            The molecule in which to find nucleophilic atoms.
        serves_target : bool
            If True, the linker function will return all identified nucleophilic sites. If False, it will return only the most accessible site based on steric constraints.

        Returns
        -------
        linkers : list of base_classes.Atom or base_classes.Atom
            The identified nucleophilic atoms based on steric constraints.
        deleters : list of base_classes.Atom or base_classes.Atom
            The atoms to be deleted during the reaction.
        """
        return self.find_atoms(mol, role="nucleophile", serves_target=serves_target)

    def find_electrophilic_atoms(self, mol: core.Molecule, serves_target: bool):
        """
        Find electrophilic atoms in the given molecule based on the specified reactivity pattern.

        Parameters
        ----------
        mol : core.Molecule
            The molecule in which to find electrophilic atoms.
        serves_target : bool
            If True, the linker function will return all identified nucleophilic sites. If False, it will return only the most accessible site based on steric constraints.

        Returns
        -------
        linkers : list of base_classes.Atom or base_classes.Atom
            The identified electrophilic atoms based on steric constraints.
        deleters : list of base_classes.Atom or base_classes.Atom
            The atoms to be deleted during the reaction.
        """
        return self.find_atoms(mol, role="electrophile", serves_target=serves_target)

    def as_nucleophile(self, serves_target: bool):
        """
        Get the nucleophilic linker and deleter functions.

        Parameters
        ----------
        serves_target : bool
            If True, the linker function will return all identified nucleophilic sites. If False, it will return only the most accessible site based on steric constraints.
        """
        self._serves_target = serves_target
        return self._nucleophile_linker_call, self._nucleophile_deleter_call

    def as_electrophile(self, serves_target: bool):
        """
        Get the electrophilic linker and deleter functions.

        Parameters
        ----------
        serves_target : bool
            If True, the linker function will return all identified nucleophilic sites. If False, it will return only the most accessible site based on steric constraints.
        """
        self._serves_target = serves_target
        return self._electrophile_linker_call, self._electrophile_deleter_call

    def _linker_call_wrapper(self, func: callable, mol: core.Molecule):
        atoms = func(mol)
        if self._steric_constraints_func is not None:
            atoms = self._steric_constraints_func(mol, atoms)
        atoms = self._default_steric_constraint_func(mol, atoms)
        if isinstance(atoms, base_classes.Atom):
            return [atoms]
        if self._serves_target:
            return list(atoms)
        else:
            return next(iter(atoms))

    def _nucleophile_linker_call(self, mol: core.Molecule):
        if self._nucleophile_linker is None:
            raise NotImplementedError("Nucleophile linker function not defined")
        linker = self._linker_call_wrapper(self._nucleophile_linker, mol)
        return linker

    def _electrophile_linker_call(self, mol: core.Molecule):
        if self._electrophile_linker is None:
            raise NotImplementedError("Electrophile linker function not defined")
        return self._linker_call_wrapper(self._electrophile_linker, mol)

    def _nucleophile_deleter_call(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        if self._nucleophile_deleter is None:
            return mol.get_hydrogen(atom)
        return self._nucleophile_deleter(atom, mol)

    def _electrophile_deleter_call(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        if self._electrophile_deleter is None:
            return mol.get_hydrogen(atom)
        return self._electrophile_deleter(atom, mol)

    def _default_steric_constraint_func(
        self, mol: core.Molecule, atoms: list[base_classes.Atom]
    ):
        if not isinstance(atoms, (list, set, tuple)):
            atoms = [atoms]
        # count close by atoms to get the most accessible site
        if self._steric_max_neighbors is not None:
            atoms = filter(
                lambda a: len(mol.get_atoms_within(a, self._steric_distance))
                <= self._steric_max_neighbors,
                atoms,
            )
        atoms = sorted(
            atoms, key=lambda n: len(mol.get_atoms_within(n, self._steric_distance))
        )
        if self._steric_n_target_sites == "all":
            return atoms
        else:
            return atoms[: self._steric_n_target_sites]


class ReactionError(Exception):
    pass


from buildamol.structural.neighbors import constraints_v2 as constraints


class Carboxyl(Reactivity):
    """
    Predefined reactivity pattern for carboxylic acids
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        C = self.electrophile_linker(mol)
        O = []
        for c in C:
            O.extend(
                mol.get_neighbors(
                    c,
                    filter=constraints.and_(
                        constraints.has_element("O"),
                        constraints.has_single_bond_with("C"),
                    ),
                )
            )
        if len(O) == 0:
            raise ReactionError("No carboxylic acid found")
        return O

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_double_bond_with("O"),
            constraints.has_single_bond_with("O"),
            constraints.not_(constraints.neighbors_any("N", "S", "P")),
        )

        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No carboxylic acid found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.and_(
            constraints.has_element("O"),
            constraints.has_single_bond_with("C"),
        )
        return mol.get_neighbors(atom, filter=filter).pop()


class Amide(Reactivity):
    """
    Predefined reactivity pattern for amide groups
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        C = self.electrophile_linker(mol)
        N = []
        for c in C:
            N.extend(
                mol.get_neighbors(
                    c,
                    filter=constraints.and_(
                        constraints.has_element("N"),
                        constraints.has_single_bond_with("C"),
                    ),
                )
            )
        if len(N) == 0:
            raise ReactionError("No amide group found")
        return N

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_double_bond_with("O"),
            constraints.has_single_bond_with("N"),
            constraints.not_(constraints.neighbors_any("S", "P")),
        )

        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No amide group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.has_element("N")
        return mol.get_neighbors(atom, filter=filter).pop()


class Ester(Carboxyl):
    """
    Predefined reactivity pattern for ester groups
    Can act as electrophile.
    """

    def nucleophile_linker(self, mol):
        raise NotImplementedError("Ester cannot act as nucleophile")


class Aldehyde(Reactivity):
    """
    Predefined reactivity pattern for aldehyde groups
    Can act as electrophile.
    """

    def nucleophile_linker(self, mol):
        raise NotImplementedError("Aldehyde cannot act as nucleophile")

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_double_bond_with("O"),
            constraints.has_single_bond_with("C"),
            constraints.not_(constraints.neighbors_any("N", "S", "P")),
        )

        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No aldehyde group found")
        return C


class Ketone(Reactivity):
    """
    Predefined reactivity pattern for ketone groups
    Can act as electrophile.
    The deleter function will remove the smaller of the two alkyl substituents or the one with more heteroatoms close-by.
    """

    def nucleophile_linker(self, mol):
        raise NotImplementedError("Ketone cannot act as nucleophile")

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_double_bond_with("O"),
            constraints.has_neighbor_hist({"C": 2, "O": 1}),
        )

        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No ketone group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        neighbors = list(mol.get_neighbors(atom, filter=constraints.has_element("C")))
        if len(neighbors) != 2:
            raise ReactionError("Ketone carbon does not have two carbon neighbors")
        a, b = neighbors[:2]

        a_neighbors_hist = sum(n.atomic_number for n in mol.get_neighbors(a, n=3))
        b_neighbors_hist = sum(n.atomic_number for n in mol.get_neighbors(b, n=3))
        if a_neighbors_hist != b_neighbors_hist:
            if a_neighbors_hist < b_neighbors_hist:
                return b
            return a

        descendants_a = mol.get_descendants(atom, a)
        descendants_b = mol.get_descendants(atom, b)
        if len(descendants_a) < len(descendants_b):
            return a
        return b


class Hydroxyl(Reactivity):
    """
    Predefined reactivity pattern for hydroxyl groups
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_single_bond_with("C"),
            constraints.neighbors_exactly("H", "C"),
        )
        O = mol.get_atoms("O", by="element", filter=filter)
        if len(O) == 0:
            raise ReactionError("No hydroxyl group found")
        return O

    def electrophile_linker(self, mol: core.Molecule):
        nucleophile_O = self.nucleophile_linker(mol)
        filter = constraints.and_(
            constraints.has_element("C"),
            constraints.not_(constraints.neighbors_any("N", "S", "P")),
        )
        C = []
        for o in nucleophile_O:
            c = mol.get_neighbors(o, filter=filter).pop()
            C.append(c)
        if len(C) == 0:
            raise ReactionError("No hydroxyl group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.has_element("O")
        return mol.get_neighbors(atom, filter=filter).pop()


class Amine(Reactivity):
    """
    Predefined reactivity pattern for amine groups
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_single_bond_with("C"),
            constraints.not_(constraints.neighbors_any("O", "S", "P")),
        )
        N = mol.get_atoms("N", by="element", filter=filter)
        if len(N) == 0:
            raise ReactionError("No amine group found")
        return N

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_single_bond_with("N"),
            constraints.not_(constraints.neighbors_any("O", "S", "P")),
        )
        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No amine group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.has_element("N")
        return mol.get_neighbors(atom, filter=filter).pop()


class Thiol(Reactivity):
    """
    Predefined reactivity pattern for thiol groups
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_single_bond_with("C"),
            constraints.not_(constraints.neighbors_any("N", "O", "P")),
        )
        S = mol.get_atoms("S", by="element", filter=filter)
        if len(S) == 0:
            raise ReactionError("No thiol group found")
        return S

    def electrophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_single_bond_with("S"),
            constraints.not_(constraints.has_double_bonds()),
            constraints.not_(constraints.neighbors_any("N", "O", "P")),
        )
        C = mol.get_atoms("C", by="element", filter=filter)
        if len(C) == 0:
            raise ReactionError("No thiol group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.has_element("S")
        return mol.get_neighbors(atom, filter=filter).pop()


class AlkylHalide(Reactivity):
    """
    Predefined reactivity pattern for alkyl halides (F, Cl, Br, I)
    Can act as electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        raise NotImplementedError("Alkyl halide cannot act as nucleophile")

    def electrophile_linker(self, mol: core.Molecule):
        halides = ("F", "CL", "BR", "I")
        filter = constraints.neighbors_any(*halides)
        C = mol.get_atoms("C", by="element", filter=filter)
        halide_neighbors = [
            mol.get_neighbors(c, filter=constraints.has_any_element_of(*halides)).pop()
            for c in C
        ]
        C = zip(C, halide_neighbors)
        C = sorted(C, key=lambda c: halides.index(c[1].element))
        C = [c[0] for c in C]
        if len(C) == 0:
            raise ReactionError("No alkyl halide group found")
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.has_any_element_of(*{"F", "CL", "BR", "I"})
        return mol.get_neighbors(atom, filter=filter).pop()


class Phosphate(Reactivity):
    """
    Predefined reactivity pattern for phosphate groups
    Can act as both nucleophile and electrophile.
    """

    def nucleophile_linker(self, mol: core.Molecule):
        filter = constraints.and_(
            constraints.has_double_bond_with("O"),
            constraints.has_single_bond_with("O"),
            constraints.not_(constraints.neighbors_any("N", "S", "C")),
        )

        P = mol.get_atoms("P", by="element", filter=filter)
        O = []
        for p in P:
            O.extend(
                mol.get_neighbors(
                    p,
                    filter=constraints.and_(
                        constraints.has_element("O"),
                        constraints.has_single_bond_with("P"),
                    ),
                )
            )
        if len(O) == 0:
            raise ReactionError("No phosphate group found")
        return O

    def electrophile_linker(self, mol: core.Molecule):
        O = self.nucleophile_linker(mol)
        C = []
        for o in O:
            C.extend(
                mol.get_neighbors(
                    o,
                    filter=constraints.has_element("C"),
                ),
            )
        if len(O) == 0:
            raise ReactionError("No phosphate group found")
        if len(C) == 0:
            return O
        return C

    def electrophile_deleter(
        self,
        atom: base_classes.Atom,
        mol: core.Molecule,
    ):
        filter = constraints.and_(
            constraints.has_element("O"),
            constraints.has_single_bond_with("P"),
        )
        return mol.get_neighbors(atom, filter=filter).pop()
