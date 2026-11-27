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

    Automatic site detection
    ------------------------
    For automatic reaction inference (see `buildamol.core.autoreact.AutoReact`), reactivities can additionally define:

    - `centers` : tuple of element symbols (upper case) of the atoms that anchor the pattern (e.g. ("C",) for a carbonyl).
      Only atoms of these elements are tested, which keeps detection fast.
    - `supersedes` : tuple of reactivity names that this reactivity takes precedence over. Atoms claimed by a
      superseding reactivity are not available to the superseded ones (e.g. a carboxyl's OH is not a hydroxyl).
    - `fallback` : if True, the reactivity is scanned last (after all non-fallback reactivities).
    - `radius` : the number of bonds around a center atom that `match` looks at (default 3). This is used to
      incrementally update cached sites after a molecule was modified (see `SiteMap`); only atoms within this
      distance of a modification are re-tested. Ring membership checks are exempt since connecting two
      molecules never creates or breaks rings.
    - `match(ctx, atom)` : test a single center atom using a `Perception` context. Return None if the atom does
      not match, otherwise a tuple `(footprint, specs)` where `footprint` is a tuple of the atoms that make up the
      matched group (these are claimed) and `specs` is a list of `SiteSpec` describing the reactive sites.

    Reactivities that do not define `centers` are detected using their linker/deleter functions instead (slower).
    """

    centers: tuple = None
    supersedes: tuple = ()
    fallback: bool = False
    radius: int = 3

    @property
    def name(self) -> str:
        """
        The name of the reactivity (used in the pairing table and for `supersedes` references).
        """
        return type(self).__name__

    def match(self, ctx: "Perception", atom: base_classes.Atom):
        """
        Test whether a center atom matches this reactivity pattern.

        Parameters
        ----------
        ctx : Perception
            The perception context of the molecule (cached graph lookups).
        atom : Atom
            The center atom to test (its element is one of `centers`).

        Returns
        -------
        tuple or None
            None if the atom does not match. Otherwise `(footprint, specs)` where footprint is a tuple of atoms
            belonging to the matched group and specs is a list of `SiteSpec` objects.
        """
        return None

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
        atoms = list(atoms)
        if len(atoms) == 0:
            return atoms
        # count close by atoms to get the most accessible site
        # (one batched distance computation for all candidate sites)
        counts = mol.count_atoms_within(atoms, self._steric_distance)
        pairs = list(zip(atoms, counts))
        if self._steric_max_neighbors is not None:
            pairs = [p for p in pairs if p[1] <= self._steric_max_neighbors]
        pairs.sort(key=lambda p: p[1])  # stable: ties keep input order
        atoms = [p[0] for p in pairs]
        if self._steric_n_target_sites == "all":
            return atoms
        else:
            return atoms[: self._steric_n_target_sites]


class ReactionError(Exception):
    pass


from collections import namedtuple
from dataclasses import dataclass

SiteSpec = namedtuple(
    "SiteSpec", ["role", "bonder", "deletes", "score", "leaving_hydrogens"]
)
SiteSpec.__doc__ = """
A reactive site as reported by `Reactivity.match`.

role : str
    Either "nucleophile" or "electrophile".
bonder : Atom
    The atom that forms the new bond.
deletes : tuple of Atom
    The atoms to delete when the bond is formed (atoms downstream of them are deleted automatically by the linkage).
score : float
    The intrinsic reactivity of the site in this role (0-1, higher is more reactive).
leaving_hydrogens : bool
    If True, `deletes` are hydrogens of which only as many as the bond order requires are removed.
"""


@dataclass(frozen=True)
class Site:
    """
    A detected reactive site in a molecule.
    """

    reactivity: str
    role: str
    bonder: base_classes.Atom
    deletes: tuple
    score: float
    footprint: tuple
    leaving_hydrogens: bool = False

    def __repr__(self):
        return (
            f"Site({self.reactivity}, {self.role}, bonder={self.bonder.id}"
            f"@{self.bonder.parent.id[1]}, deletes={[d.id for d in self.deletes]}, score={self.score:.2f})"
        )


class Perception:
    """
    A lightweight, cached view on a molecule's graph used for reactive site detection.

    Parameters
    ----------
    mol : Molecule
        The molecule to perceive.
    atoms : iterable of Atom, optional
        Restrict the candidate center atoms to these (e.g. the atoms of a single residue).
        Neighborhood lookups still see the full molecule.
    """

    __slots__ = ("mol", "adj", "atoms", "by_element", "_el")

    def __init__(self, mol: core.Molecule, atoms=None):
        self.mol = mol
        self.adj = mol._AtomGraph.adj
        atoms = list(mol.get_atoms()) if atoms is None else list(atoms)
        atoms.sort(key=lambda a: a.serial_number)
        self.atoms = atoms
        self._el = {}
        by_element = {}
        for a in atoms:
            by_element.setdefault(self.el(a), []).append(a)
        self.by_element = by_element

    def el(self, atom) -> str:
        """Upper case element symbol of an atom"""
        e = self._el.get(atom)
        if e is None:
            e = atom.element.upper()
            self._el[atom] = e
        return e

    def bonds(self, atom) -> list:
        """List of (neighbor, bond_order) tuples"""
        return [(n, d.get("bond_order", 1)) for n, d in self.adj[atom].items()]

    def hydrogens(self, atom) -> list:
        """The hydrogen neighbors of an atom (sorted by serial number)"""
        return sorted(
            (n for n in self.adj[atom] if self.el(n) == "H"),
            key=lambda a: a.serial_number,
        )

    def heavy(self, atom) -> list:
        """List of (neighbor, bond_order) tuples of non-hydrogen neighbors"""
        return [(n, o) for n, o in self.bonds(atom) if self.el(n) != "H"]

    def has_multiple_bond(self, atom, to: tuple = None) -> bool:
        """Whether the atom has a double or triple bond (optionally only to the given elements)"""
        for n, o in self.bonds(atom):
            if o > 1 and (to is None or self.el(n) in to):
                return True
        return False

    def terminal_double(self, atom, element: str):
        """Return a neighbor of the given element that is double bonded and has no other neighbors (e.g. a carbonyl O)"""
        for n, o in self.bonds(atom):
            if o == 2 and self.el(n) == element and len(self.adj[n]) == 1:
                return n
        return None

    def in_ring(self, atom1, atom2) -> bool:
        """Whether the bond atom1-atom2 is part of a ring"""
        seen = {atom2}
        stack = [atom2]
        while stack:
            current = stack.pop()
            for n in self.adj[current]:
                if current is atom2 and n is atom1:
                    continue
                if n is atom1:
                    return True
                if n not in seen:
                    seen.add(n)
                    stack.append(n)
        return False


from buildamol.structural.neighbors import constraints_v2 as constraints


class Carboxyl(Reactivity):
    """
    Predefined reactivity pattern for carboxylic acids
    Can act as both nucleophile and electrophile.
    """

    centers = ("C",)
    supersedes = ("Hydroxyl", "Aldehyde", "Ketone")

    def match(self, ctx, atom):
        O_double = ctx.terminal_double(atom, "O")
        if O_double is None:
            return None
        O_h = None
        for n, o in ctx.bonds(atom):
            if n is O_double:
                continue
            e = ctx.el(n)
            if e == "O" and o == 1 and O_h is None:
                H = ctx.hydrogens(n)
                if len(H) == 1 and len(ctx.adj[n]) == 2:
                    O_h = n
                    continue
                return None
            if e not in ("C", "H") or o != 1:
                return None
        if O_h is None:
            return None
        H = ctx.hydrogens(O_h)
        footprint = (atom, O_double, O_h, *H)
        specs = [
            SiteSpec("electrophile", atom, (O_h,), 1.0, False),
            SiteSpec("nucleophile", O_h, tuple(H), 0.4, True),
        ]
        return footprint, specs

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
    (For automatic detection only the N-H is considered reactive, as a weak nucleophile.)
    """

    centers = ("C",)
    supersedes = ("Amine", "Aldehyde", "Ketone")

    def match(self, ctx, atom):
        O_double = ctx.terminal_double(atom, "O")
        if O_double is None:
            return None
        N = [n for n, o in ctx.bonds(atom) if o == 1 and ctx.el(n) == "N"]
        if len(N) == 0:
            return None
        footprint = [atom, O_double]
        specs = []
        for n in N:
            H = ctx.hydrogens(n)
            footprint.extend((n, *H))
            if H:
                specs.append(SiteSpec("nucleophile", n, tuple(H), 0.3, True))
        return tuple(footprint), specs

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

    centers = ("C",)
    supersedes = ("Hydroxyl", "Aldehyde", "Ketone")

    def match(self, ctx, atom):
        O_double = ctx.terminal_double(atom, "O")
        if O_double is None:
            return None
        O_alk = None
        for n, o in ctx.bonds(atom):
            if n is O_double:
                continue
            e = ctx.el(n)
            if e == "O" and o == 1 and O_alk is None:
                heavy = ctx.heavy(n)
                if len(heavy) == 2 and all(
                    ctx.el(h) == "C" and ho == 1 for h, ho in heavy
                ):
                    O_alk = n
                    continue
                return None
            if e not in ("C", "H") or o != 1:
                return None
        if O_alk is None:
            return None
        footprint = (atom, O_double, O_alk)
        specs = []
        # in a lactone the alkoxy group cannot leave without deleting the ring
        if not ctx.in_ring(atom, O_alk):
            specs.append(SiteSpec("electrophile", atom, (O_alk,), 0.6, False))
        return footprint, specs

    def nucleophile_linker(self, mol):
        raise NotImplementedError("Ester cannot act as nucleophile")


class Aldehyde(Reactivity):
    """
    Predefined reactivity pattern for aldehyde groups
    Can act as electrophile.
    (For automatic detection the carbonyl O is the leaving group, i.e. condensations such as imine formation.)
    """

    centers = ("C",)
    supersedes = ("Ketone",)

    def match(self, ctx, atom):
        O_double = ctx.terminal_double(atom, "O")
        if O_double is None:
            return None
        H = ctx.hydrogens(atom)
        if len(H) == 0:
            return None
        heavy = [(n, o) for n, o in ctx.heavy(atom) if n is not O_double]
        if len(heavy) > 1 or any(ctx.el(n) != "C" or o != 1 for n, o in heavy):
            return None
        footprint = (atom, O_double, *H)
        return footprint, [SiteSpec("electrophile", atom, (O_double,), 0.8, False)]

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
    (For automatic detection the carbonyl O is the leaving group instead, i.e. condensations such as imine formation.)
    """

    centers = ("C",)

    def match(self, ctx, atom):
        O_double = ctx.terminal_double(atom, "O")
        if O_double is None:
            return None
        heavy = [(n, o) for n, o in ctx.heavy(atom) if n is not O_double]
        if len(ctx.adj[atom]) != 3 or len(heavy) != 2:
            return None
        if any(ctx.el(n) != "C" or o != 1 for n, o in heavy):
            return None
        footprint = (atom, O_double)
        return footprint, [SiteSpec("electrophile", atom, (O_double,), 0.6, False)]

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

    centers = ("O",)

    def match(self, ctx, atom):
        bonds = ctx.bonds(atom)
        if len(bonds) != 2:
            return None
        H = [n for n, _ in bonds if ctx.el(n) == "H"]
        C = [(n, o) for n, o in bonds if ctx.el(n) == "C"]
        if len(H) != 1 or len(C) != 1 or C[0][1] != 1:
            return None
        C = C[0][0]
        # acids, carbamic acids etc. are not alcohols
        if ctx.has_multiple_bond(C, to=("O", "N", "S")):
            return None
        footprint = (atom, *H)
        if ctx.has_multiple_bond(C):
            # phenols / enols: nucleophile only, no substitution at sp2 carbon
            return footprint, [SiteSpec("nucleophile", atom, tuple(H), 0.45, True)]
        specs = [
            SiteSpec("nucleophile", atom, tuple(H), 0.5, True),
            SiteSpec("electrophile", C, (atom,), 0.15, False),
        ]
        return footprint, specs

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
    (For automatic detection only primary and secondary amines are considered, as nucleophiles.)
    """

    centers = ("N",)

    def match(self, ctx, atom):
        bonds = ctx.bonds(atom)
        if any(o != 1 for _, o in bonds):
            return None
        H = ctx.hydrogens(atom)
        heavy = [n for n, _ in bonds if ctx.el(n) != "H"]
        if len(H) == 0 or any(ctx.el(n) != "C" for n in heavy):
            return None
        # N next to a carbonyl / thiocarbonyl / imine carbon is not an amine
        if any(ctx.has_multiple_bond(n, to=("O", "S", "N")) for n in heavy):
            return None
        score = 1.0 if len(heavy) <= 1 else 0.9
        if any(ctx.has_multiple_bond(n) for n in heavy):
            score -= 0.3  # aryl / vinyl amines are weaker nucleophiles
        footprint = (atom, *H)
        return footprint, [SiteSpec("nucleophile", atom, tuple(H), score, True)]

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
    (For automatic detection only the S-H is considered reactive, as a nucleophile.)
    """

    centers = ("S",)

    def match(self, ctx, atom):
        bonds = ctx.bonds(atom)
        if len(bonds) != 2 or any(o != 1 for _, o in bonds):
            return None
        H = [n for n, _ in bonds if ctx.el(n) == "H"]
        C = [n for n, _ in bonds if ctx.el(n) == "C"]
        if len(H) != 1 or len(C) != 1:
            return None
        if ctx.has_multiple_bond(C[0], to=("O", "S", "N")):
            return None  # thioacids
        footprint = (atom, *H)
        return footprint, [SiteSpec("nucleophile", atom, tuple(H), 0.9, True)]

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
    (For automatic detection aryl/vinyl halides are included as somewhat weaker electrophiles, i.e. cross-coupling partners.)
    """

    centers = ("I", "BR", "CL", "F")
    _halide_scores = {"I": 0.9, "BR": 0.85, "CL": 0.7, "F": 0.2}

    def match(self, ctx, atom):
        bonds = ctx.bonds(atom)
        if len(bonds) != 1:
            return None
        C, o = bonds[0]
        if ctx.el(C) != "C" or o != 1:
            return None
        score = self._halide_scores[ctx.el(atom)]
        if ctx.has_multiple_bond(C):
            score *= 0.8
        return (atom,), [SiteSpec("electrophile", C, (atom,), score, False)]

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
    (For automatic detection each P-OH oxygen is a nucleophile.)
    """

    centers = ("P",)
    supersedes = ("Hydroxyl",)

    def match(self, ctx, atom):
        if ctx.terminal_double(atom, "O") is None:
            return None
        footprint = [atom]
        specs = []
        for n, o in ctx.bonds(atom):
            if ctx.el(n) != "O":
                continue
            H = ctx.hydrogens(n)
            footprint.extend((n, *H))
            if o == 1 and len(H) == 1 and len(ctx.adj[n]) == 2:
                specs.append(SiteSpec("nucleophile", n, tuple(H), 0.5, True))
        if not specs:
            return None
        return tuple(footprint), specs

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


class CH(Reactivity):
    """
    Fallback reactivity for C-H bonds: the carbon is linked after removing one of its hydrogens.
    This is not a mechanistic reactivity but makes it possible to connect molecules without any
    reactive functional groups (e.g. benzene). It is always scanned last, so it only applies to carbons
    not already claimed by a functional group. sp2 carbons (e.g. aromatic C-H) are preferred over sp3 carbons.
    """

    centers = ("C",)
    fallback = True

    def match(self, ctx, atom):
        H = ctx.hydrogens(atom)
        if len(H) == 0:
            return None
        score = 0.2 if ctx.has_multiple_bond(atom) else 0.1
        H = tuple(H)
        return (atom,), [
            SiteSpec("nucleophile", atom, H, score, True),
            SiteSpec("electrophile", atom, H, score, True),
        ]

    def nucleophile_linker(self, mol: core.Molecule):
        C = mol.get_atoms(
            "C", by="element", filter=lambda a: len(mol.get_hydrogens(a)) > 0
        )
        if len(C) == 0:
            raise ReactionError("No C-H found")
        return C

    electrophile_linker = nucleophile_linker


def default_reactivities() -> list:
    """
    Get (new instances of) the default reactivities used for automatic site detection.
    """
    return [
        Carboxyl(),
        Ester(),
        Amide(),
        Phosphate(),
        Aldehyde(),
        Ketone(),
        Thiol(),
        Amine(),
        Hydroxyl(),
        AlkylHalide(),
        CH(),
    ]


def order_reactivities(reactivities: list) -> list:
    """
    Order reactivities such that each reactivity comes before all reactivities it supersedes.
    Ties keep the input order and fallback reactivities always come last.

    Parameters
    ----------
    reactivities : list of Reactivity

    Returns
    -------
    list of Reactivity
    """
    regular = [r for r in reactivities if not r.fallback]
    fallback = [r for r in reactivities if r.fallback]

    def _sort(items):
        names = {r.name for r in items}
        index = {id(r): i for i, r in enumerate(items)}
        indegree = {id(r): 0 for r in items}
        edges = {id(r): [] for r in items}
        by_name = {}
        for r in items:
            by_name.setdefault(r.name, []).append(r)
        for r in items:
            for s in r.supersedes:
                if s not in names:
                    continue
                for other in by_name[s]:
                    if other is r:
                        continue
                    edges[id(r)].append(other)
                    indegree[id(other)] += 1
        ready = sorted(
            (r for r in items if indegree[id(r)] == 0), key=lambda r: index[id(r)]
        )
        out = []
        while ready:
            r = ready.pop(0)
            out.append(r)
            for other in edges[id(r)]:
                indegree[id(other)] -= 1
                if indegree[id(other)] == 0:
                    ready.append(other)
                    ready.sort(key=lambda r: index[id(r)])
        if len(out) != len(items):
            cyclic = [r.name for r in items if r not in out]
            raise ValueError(f"Cyclic 'supersedes' relations between: {cyclic}")
        return out

    return _sort(regular) + _sort(fallback)


_DEFAULT_ORDERED = None


def _default_ordered():
    global _DEFAULT_ORDERED
    if _DEFAULT_ORDERED is None:
        _DEFAULT_ORDERED = order_reactivities(default_reactivities())
    return _DEFAULT_ORDERED


def _linker_based_sites(reactivity, mol, atoms, claimed):
    """
    Detect sites of a reactivity that only defines linker/deleter functions (no `centers`/`match`).
    """
    allowed = set(atoms)
    sites = []
    for role in ("nucleophile", "electrophile"):
        try:
            linkers, deleters = reactivity.find_atoms(mol, role, serves_target=True)
        except (ReactionError, NotImplementedError, ValueError, AttributeError):
            continue
        for bonder, deletes in zip(linkers, deleters):
            if bonder is None or bonder not in allowed or bonder in claimed:
                continue
            if deletes is None:
                deletes = ()
            elif isinstance(deletes, base_classes.Atom):
                deletes = (deletes,)
            deletes = tuple(d for d in deletes if d is not None)
            sites.append(
                Site(
                    reactivity.name,
                    role,
                    bonder,
                    deletes,
                    0.5,
                    (bonder, *deletes),
                    False,
                )
            )
    return sites


class _Match:
    """A matched group (internal): the reactivity index, center atom, claimed footprint and its sites."""

    __slots__ = ("r_idx", "center", "footprint", "sites", "seq")

    def __init__(self, r_idx, center, footprint, sites, seq=0):
        self.r_idx = r_idx
        self.center = center
        self.footprint = footprint
        self.sites = sites
        self.seq = seq


def _scan(mol, reactivities, ctx, retained=()):
    """
    Scan the center atoms of a perception context in order of precedence.

    `retained` are previously detected matches (centers outside of `ctx.atoms`) that are re-inserted
    at their place in the precedence order, so that claiming behaves exactly as in a full scan.
    """
    retained_by_r = {}
    for m in retained:
        retained_by_r.setdefault(m.r_idx, []).append(m)

    claimed = set()
    matches = []
    for r_idx, reactivity in enumerate(reactivities):
        if reactivity.centers is None:
            new = _linker_based_sites(reactivity, mol, ctx.atoms, claimed)
            for seq, s in enumerate(new):
                claimed.update(s.footprint)
                matches.append(_Match(r_idx, s.bonder, frozenset(s.footprint), (s,), seq))
            continue

        new_claims = set()
        for m in retained_by_r.get(r_idx, ()):
            if m.center in claimed or not claimed.isdisjoint(m.footprint):
                continue
            new_claims.update(m.footprint)
            matches.append(m)

        name = reactivity.name
        for element in reactivity.centers:
            for atom in ctx.by_element.get(element, ()):
                if atom in claimed:
                    continue
                result = reactivity.match(ctx, atom)
                if result is None:
                    continue
                footprint, specs = result
                if not claimed.isdisjoint(footprint):
                    continue
                new_claims.update(footprint)
                footprint = tuple(footprint)
                sites = []
                for spec in specs:
                    residue = spec.bonder.parent
                    # linkages are defined per residue, so all atoms of a site must share the bonder's residue
                    if any(d.parent is not residue for d in spec.deletes):
                        continue
                    sites.append(
                        Site(
                            name,
                            spec.role,
                            spec.bonder,
                            tuple(spec.deletes),
                            spec.score,
                            footprint,
                            spec.leaving_hydrogens,
                        )
                    )
                # matches without sites (e.g. lactones) are kept since they still claim their atoms
                matches.append(_Match(r_idx, atom, frozenset(footprint), tuple(sites)))
        claimed |= new_claims
    return matches


def _ordered_sites(matches, reactivities) -> list:
    """Flatten matches into sites in the order of a full scan (precedence, center element, serial number)."""
    element_rank = []
    for r in reactivities:
        element_rank.append({e: i for i, e in enumerate(r.centers or ())})

    def key(m):
        center = m.center
        rank = element_rank[m.r_idx].get(center.element.upper(), 0)
        return (m.r_idx, rank, center.serial_number, m.seq)

    return [s for m in sorted(matches, key=key) for s in m.sites]


def _residue_atoms(mol, residues):
    if not isinstance(residues, (list, tuple, set)):
        residues = [residues]
    atoms = []
    for r in residues:
        res = mol.get_residue(r)
        if res is None:
            raise ValueError(f"Residue {r} not found in molecule")
        atoms.extend(res.get_atoms())
    return atoms


def detect_sites(
    mol: core.Molecule,
    reactivities: list = None,
    residues=None,
    ordered: bool = False,
) -> list:
    """
    Detect all reactive sites in a molecule.

    Reactivities are scanned in order of precedence (see `Reactivity.supersedes`). Each matched group
    claims its atoms, so that less specific reactivities cannot match the same atoms again
    (e.g. the OH of a carboxylic acid is not reported as a hydroxyl). Patterns are always evaluated
    on the full molecular graph; claiming only restricts which atoms are available as group members.

    Use `SiteMap` to keep detected sites around and update them incrementally when a molecule grows.

    Parameters
    ----------
    mol : Molecule
        The molecule to scan.
    reactivities : list of Reactivity, optional
        The reactivities to scan for. By default `default_reactivities()` is used.
    residues : Residue or list of Residues, optional
        Restrict the scan to atoms of these residues. By default all atoms are scanned.
    ordered : bool
        Set to True if the reactivities are already in order of precedence (skips sorting).

    Returns
    -------
    list of Site
    """
    if reactivities is None:
        reactivities = _default_ordered()
    elif not ordered:
        reactivities = order_reactivities(reactivities)

    atoms = None if residues is None else _residue_atoms(mol, residues)
    ctx = Perception(mol, atoms)
    return _ordered_sites(_scan(mol, reactivities, ctx), reactivities)


def _graph_signature(mol) -> tuple:
    """A cheap fingerprint of a molecule's graph to detect modifications behind a cached `SiteMap`."""
    graph = mol._AtomGraph
    bonds = mol._bonds
    return (
        id(graph),
        len(graph._node),
        len(bonds),
        sum(b.order for b in bonds),
    )


class _VolatileDict(dict):
    """A dict that is not carried over when its owner is pickled or deep-copied (used for caches)."""

    def __reduce__(self):
        return (type(self), ())

    def __deepcopy__(self, memo):
        return type(self)()

    def __copy__(self):
        return type(self)()


_CACHE_KEY = "_buildamol_site_maps"


class SiteMap:
    """
    The detected reactive sites of a molecule, which can be updated incrementally.

    Site detection is a pure function of the molecular graph. When a molecule is modified locally
    (e.g. by connecting another molecule to it) only the atoms around the modification need to be
    re-tested. This makes iterative assembly (e.g. growing a polymer) linear instead of quadratic.

    Parameters
    ----------
    mol : Molecule
        The molecule to scan.
    reactivities : list of Reactivity, optional
        The reactivities to scan for. By default `default_reactivities()` is used.
    ordered : bool
        Set to True if the reactivities are already in order of precedence (skips sorting).

    Maps can be cached on the molecule itself (see `SiteMap.cached` and `SiteMap.store`). The cache lives
    exactly as long as the molecule and is not carried over by `Molecule.copy` or pickling.

    Examples
    --------
    >>> smap = SiteMap(mol)
    >>> smap.sites()                          # same as detect_sites(mol)
    >>> mol.attach(other, ...)                # modify the molecule in place
    >>> smap = smap.update(mol, seeds=[...])  # only re-scan around the seeds
    """

    __slots__ = ("reactivities", "matches", "signature", "incremental", "margin", "_sites")

    def __init__(self, mol=None, reactivities=None, ordered=False):
        if reactivities is None:
            reactivities = _default_ordered()
        elif not ordered:
            reactivities = order_reactivities(reactivities)
        self.reactivities = reactivities
        self.incremental = all(r.centers is not None for r in reactivities)
        # matches depend on atoms within `radius` of their center; new claims can block or unblock
        # neighboring groups, which adds a margin of two footprint radii
        self.margin = max((r.radius for r in reactivities), default=0) + 4
        self.matches = []
        self.signature = None
        self._sites = None
        if mol is not None:
            self.matches = _scan(mol, reactivities, Perception(mol))
            self.signature = _graph_signature(mol)

    @staticmethod
    def _store(mol, create=False):
        graph = mol._AtomGraph.graph
        store = graph.get(_CACHE_KEY)
        if store is None and create:
            store = graph[_CACHE_KEY] = _VolatileDict()
        return store

    @staticmethod
    def _key(reactivities):
        return tuple(_default_ordered() if reactivities is None else reactivities)

    @classmethod
    def lookup(cls, mol, reactivities=None):
        """
        Get the map cached on a molecule, or None if there is none or it is out of sync with the molecule.

        Parameters
        ----------
        mol : Molecule
        reactivities : list of Reactivity, optional
            The (ordered) reactivities the map was detected with. By default the default reactivities.
        """
        store = cls._store(mol)
        if not store:
            return None
        smap = store.get(cls._key(reactivities))
        if smap is None or not smap.is_valid(mol):
            return None
        return smap

    @classmethod
    def cached(cls, mol, reactivities=None, ordered=False) -> "SiteMap":
        """
        Get the map cached on a molecule, detecting (and caching) it if needed.

        Parameters
        ----------
        mol : Molecule
        reactivities : list of Reactivity, optional
            The reactivities to scan for. By default `default_reactivities()` is used.
        ordered : bool
            Set to True if the reactivities are already in order of precedence (skips sorting).
        """
        if reactivities is not None and not ordered:
            reactivities = order_reactivities(reactivities)
        smap = cls.lookup(mol, reactivities)
        if smap is None:
            smap = cls(mol, reactivities, ordered=True).store(mol)
        return smap

    def store(self, mol) -> "SiteMap":
        """Cache this map on a molecule (it must describe that molecule). Returns the map itself."""
        self._store(mol, create=True)[self._key(self.reactivities)] = self
        return self

    @classmethod
    def forget(cls, mol):
        """Remove all maps cached on a molecule"""
        store = cls._store(mol)
        if store:
            store.clear()

    def sites(self) -> list:
        """The detected sites (in the same order as `detect_sites`)"""
        if self._sites is None:
            self._sites = _ordered_sites(self.matches, self.reactivities)
        return list(self._sites)

    def is_valid(self, mol) -> bool:
        """Whether the map is (still) in sync with the molecule's graph"""
        return self.signature is not None and self.signature == _graph_signature(mol)

    def copy_to(self, original, copy) -> "SiteMap":
        """
        Transfer the sites onto a copy of the molecule they were detected in.

        Parameters
        ----------
        original : Molecule
            The molecule the map was detected in.
        copy : Molecule
            An (unmodified) copy of it, e.g. from `original.copy()`. Atoms are matched by
            their position in `get_atoms()`, which `Molecule.copy` preserves.

        Returns
        -------
        SiteMap
        """
        a = list(original.get_atoms())
        b = list(copy.get_atoms())
        if len(a) != len(b):
            raise ValueError("The copy does not have the same atoms as the original")
        mapping = dict(zip(a, b))
        mol = copy
        new = SiteMap(None, self.reactivities, ordered=True)
        get = mapping.__getitem__
        matches = []
        for m in self.matches:
            sites = tuple(
                Site(
                    s.reactivity,
                    s.role,
                    get(s.bonder),
                    tuple(map(get, s.deletes)),
                    s.score,
                    tuple(map(get, s.footprint)),
                    s.leaving_hydrogens,
                )
                for s in m.sites
            )
            matches.append(
                _Match(m.r_idx, get(m.center), frozenset(map(get, m.footprint)), sites, m.seq)
            )
        new.matches = matches
        new.signature = _graph_signature(mol)
        return new

    def update(self, mol, seeds) -> "SiteMap":
        """
        Update the map in place after a local modification of the molecule.

        Parameters
        ----------
        mol : Molecule
            The modified molecule (the one the map was detected in, modified in place).
        seeds : iterable of Atom
            The atoms whose bonds changed (e.g. the atoms of a newly formed bond) and all newly added atoms.
            Removed atoms do not need to be listed.

        Returns
        -------
        SiteMap
            The map itself (updated).
        """
        self._sites = None
        if not self.incremental:
            self.matches = _scan(mol, self.reactivities, Perception(mol))
            self.signature = _graph_signature(mol)
            return self

        adj = mol._AtomGraph.adj
        frontier = [a for a in seeds if a in adj]
        dirty = set(frontier)
        for _ in range(self.margin):
            nxt = []
            for a in frontier:
                for n in adj[a]:
                    if n not in dirty:
                        dirty.add(n)
                        nxt.append(n)
            if not nxt:
                break
            frontier = nxt

        # removed atoms are only reachable through the seeds, so groups that lost atoms
        # are either centered on a removed atom or lie within the dirty region
        retained = [m for m in self.matches if m.center not in dirty and m.center in adj]
        self.matches = _scan(mol, self.reactivities, Perception(mol, dirty), retained)
        self.signature = _graph_signature(mol)
        return self

    def __len__(self):
        return sum(len(m.sites) for m in self.matches)

    def __repr__(self):
        return f"SiteMap({len(self)} sites)"


__all__ = [
    "Reactivity",
    "ReactionError",
    "Site",
    "SiteSpec",
    "Perception",
    "Carboxyl",
    "Amide",
    "Ester",
    "Aldehyde",
    "Ketone",
    "Hydroxyl",
    "Amine",
    "Thiol",
    "AlkylHalide",
    "Phosphate",
    "CH",
    "default_reactivities",
    "order_reactivities",
    "detect_sites",
    "SiteMap",
]
