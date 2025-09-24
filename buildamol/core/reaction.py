import buildamol.core as core
import buildamol.base_classes as base_classes

from typing import *

import inspect
from itertools import product


class Reaction:
    """
    A class representing a (pseudo-) chemical reaction between two molecules.
    It serves as a factory to create Linkages from callables rather than direct atom identifiers.

    Parameters
    ----------
    atom1 : Union[Atom, callable]
        The atom in the target molecule to which the source molecule will be connected.
        This can be an Atom object or a callable that takes a Molecule and returns an Atom
    atom2 : Union[Atom, callable]
        The atom in the source molecule which will be connected to the target molecule.
        This can be an Atom object or a callable that takes a Molecule and returns an Atom
    delete_in_target : Union[List, callable], optional
        A list of atoms in the target molecule to be deleted upon connection,
        or a callable that takes the atom1 and the target molecule and returns such a list.
        Default Hydrogen-deletion is applied if None.
    delete_in_source : Union[List, callable], optional
        A list of atoms in the source molecule to be deleted upon connection,
        or a callable that takes the atom2 and the source molecule and returns such a list.
        Default Hydrogen-deletion is applied if None.
    bond_order : int, optional
        The bond order of the new bond formed between atom1 and atom2. Default is 1 (single bond).
    """

    def __init__(
        self,
        atom1: Union[base_classes.Atom, callable],
        atom2: Union[base_classes.Atom, callable],
        delete_in_target: Union[List, callable] = None,
        delete_in_source: Union[List, callable] = None,
        bond_order: int = 1,
    ):
        self._atom1 = atom1
        self._atom2 = atom2
        self._delete_in_target = self._check_delete_callable_signature(delete_in_target)
        self._delete_in_source = self._check_delete_callable_signature(delete_in_source)
        self._bond_order = bond_order

        self._memory = {
            "target": None,
            "source": None,
            "atom1": None,
            "atom2": None,
            "delete_in_target": None,
            "delete_in_source": None,
        }

    @classmethod
    def from_reactivities(
        cls,
        nucleophile: "Reactivity",
        electrophile: "Reactivity",
        bond_order: int = 1,
        target_is_electrophile: bool = True,
    ):
        """
        Set up a Reaction from two Reactivity objects, one for the nucleophile (source) and one for the electrophile (target).
        This is a convenience method to quickly create a Reaction from predefined Reactivity patterns.

        Parameters
        ----------
        nucleophile : Reactivity
            The Reactivity object defining the nucleophilic behavior of the source molecule.
        electrophile : Reactivity
            The Reactivity object defining the electrophilic behavior of the target molecule.
        bond_order : int, optional
            The bond order of the new bond formed between the nucleophile and electrophile. Default is 1 (single bond).
        target_is_electrophile : bool, optional
            Set to False to modify the roles of nucleophile and electrophile, i.e. the target molecule is the nucleophile and the source molecule is the electrophile.
        """
        atom1, delete_in_target = electrophile.as_electrophile(
            serves_target=target_is_electrophile
        )
        atom2, delete_in_source = nucleophile.as_nucleophile(
            serves_target=not target_is_electrophile
        )
        if not target_is_electrophile:
            atom1, atom2 = atom2, atom1
            delete_in_target, delete_in_source = delete_in_source, delete_in_target
        return cls(
            atom1=atom1,
            atom2=atom2,
            delete_in_target=delete_in_target,
            delete_in_source=delete_in_source,
            bond_order=bond_order,
        )

    def apply(
        self, target: "Molecule", source: "Molecule", inplace: bool = False
    ) -> "Molecule":
        """
        Apply the reaction to two molecules, creating a new molecule with the linkage applied.
        This is the same as calling the Reaction object directly.

        Parameters
        ----------
        target : Molecule
            The target molecule to which the source molecule will be connected.
        source : Molecule
            The source molecule which will be connected to the target molecule.
        inplace : bool, optional
            If True, modify the target molecule in place. If False, create a copy of the target molecule. Default is False.

        Returns
        -------
        Molecule
            A new Molecule object with the linkage applied.
        """
        if not self.can_apply(target, source):
            raise ValueError("Cannot create linkage: preconditions not met.")
        link = self.create_linkage(target, source)
        if isinstance(link, list):
            if not inplace:
                target = target.copy()
                source = source.copy()
            for l, (atom1, atom2) in zip(
                link, zip(self._memory["atom1"], self._memory["atom2"])
            ):
                if target is source:
                    l.apply(target, source, atom1.parent, atom2.parent)
                else:
                    target.set_attach_residue(atom1.parent)
                    source.set_attach_residue(atom2.parent)
                    target = self._apply_link(
                        target, source, l, inplace_a=True, inplace_b=inplace
                    )

        else:
            if target is source:
                link.apply(
                    target,
                    source,
                    self._memory["atom1"][0].parent,
                    self._memory["atom2"][0].parent,
                )
            else:
                target.attach_residue = self._memory["atom1"][0].parent
                source.attach_residue = self._memory["atom2"][0].parent
                target = self._apply_link(
                    target, source, link, inplace_a=inplace, inplace_b=inplace
                )
        return target

    def _apply_link(self, target, source, link, inplace_a, inplace_b):
        atom1 = target.attach_residue.get_atom(link.atom1)
        atom2 = source.attach_residue.get_atom(link.atom2)
        out = core.connect(
            target,
            source,
            link,
            copy_a=not inplace_a,
            copy_b=not inplace_b,
            use_patch=False,
        )
        if self._bond_order > 1:
            out.set_bond_order(
                atom1,
                atom2,
                self._bond_order,
                adjust_hydrogens=True,
            )
        return out

    def create_linkage(self, target: "Molecule", source: "Molecule") -> "Linkage":
        """
        Create one or more Linkage object(s) based on the current reaction parameters and the provided molecules.
        This does not modify the molecules, it only creates the Linkage object(s).

        This method is automatically called by the `apply` method. It requires that `can_apply` has been called beforehand to ensure that the reaction can be applied.

        Parameters
        ----------
        target : Molecule
            The target molecule to which the source molecule will be connected.
        source : Molecule
            The source molecule which will be connected to the target molecule.

        Returns
        -------
        Linkage or List[Linkage]
            A Linkage object or a list of Linkage objects representing the connection(s) to be
            made between the target and source molecules.
        """
        links = list(self._yield_linkages(target, source))
        if len(links) == 1:
            links = links[0]
        return links

    def _yield_linkages(self, target: "Molecule", source: "Molecule") -> "Linkage":
        for i in range(len(self._memory["atom1"])):
            atom1 = self._memory["atom1"][i]
            atom2 = self._memory["atom2"][i]
            delete_in_target = self._memory["delete_in_target"][i]
            delete_in_source = self._memory["delete_in_source"][i]

            link = core.linkage(
                atom1,
                atom2,
                delete_in_target=delete_in_target,
                delete_in_source=delete_in_source,
            )
            yield link

    def can_apply(self, target: "Molecule", source: "Molecule") -> bool:
        """
        Check if the reaction can be applied to the given target and source molecules.
        This checks if the specified atoms and deletions are valid in the context of the provided molecules
        and stores the resolved atoms and deletions in memory for later use.

        This method is automatically called by the `apply` method.

        Parameters
        ----------
        target : Molecule
            The target molecule to which the source molecule will be connected.
        source : Molecule
            The source molecule which will be connected to the target molecule.

        Returns
        -------
        bool
            True if the reaction can be applied, False otherwise.
        """
        atom1 = self._apply_atom_getter(self._atom1, target)
        atom2 = self._apply_atom_getter(self._atom2, source)

        for a1 in atom1:
            if a1 is None:
                return False
        for a2 in atom2:
            if a2 is None:
                return False

        if len(atom2) > 1:
            raise ValueError(
                "atom2 resolved to multiple atoms, but only a single atom is allowed for the source atom. Reacting at multiple sites is supported but only for the target molecule (atom1). "
            )

        valid_atoms1 = []
        valid_atoms2 = []
        valid_deletes1 = []
        valid_deletes2 = []
        for a1, a2 in product(atom1, atom2):
            deletes = self._find_deletes_for_anchors(target, source, a1, a2)
            if deletes is None:
                continue
            deletes1, deletes2 = deletes
            deletes1 = list(deletes1) if deletes1 is not None else None
            deletes2 = list(deletes2) if deletes2 is not None else None
            valid_atoms1.append(a1)
            valid_atoms2.append(a2)
            valid_deletes1.append(deletes1)
            valid_deletes2.append(deletes2)

        if len(valid_deletes1) == 0:
            return False

        self._memory["target"] = target
        self._memory["source"] = source
        self._memory["atom1"] = valid_atoms1
        self._memory["atom2"] = valid_atoms2
        self._memory["delete_in_target"] = valid_deletes1
        self._memory["delete_in_source"] = valid_deletes2
        return True

    def set_reactivity(
        self,
        atom1: Union[base_classes.Atom, callable] = None,
        atom2: Union[base_classes.Atom, callable] = None,
        delete_in_target: Union[List, callable] = None,
        delete_in_source: Union[List, callable] = None,
        bond_order: int = None,
    ):
        """
        Set new parameters for the reaction.

        Parameters
        ----------
        atom1 : Union[Atom, callable], optional
            The atom in the target molecule to which the source molecule will be connected.
            This can be an Atom object or a callable that takes a Molecule and returns an Atom
        atom2 : Union[Atom, callable], optional
            The atom in the source molecule which will be connected to the target molecule.
            This can be an Atom object or a callable that takes a Molecule and returns an Atom
        delete_in_target : Union[List, callable], optional
            A list of atoms in the target molecule to be deleted upon connection,
            or a callable that takes the atom1 and the target molecule and returns such a list.
            Default Hydrogen-deletion is applied if None.
        delete_in_source : Union[List, callable], optional
            A list of atoms in the source molecule to be deleted upon connection,
            or a callable that takes the atom2 and the source molecule and returns such a list.
            Default Hydrogen-deletion is applied if None.
        bond_order : int, optional
            The bond order of the new bond formed between atom1 and atom2. Default is 1 (single bond).
        """
        if atom1 is not None:
            self._atom1 = atom1
        if atom2 is not None:
            self._atom2 = atom2
        if delete_in_target is not None:
            self._delete_in_target = self._check_delete_callable_signature(
                delete_in_target
            )
        if delete_in_source is not None:
            self._delete_in_source = self._check_delete_callable_signature(
                delete_in_source
            )
        if bond_order is not None:
            self._bond_order = bond_order

    def with_reactivity(
        self,
        atom1: Union[base_classes.Atom, callable] = None,
        atom2: Union[base_classes.Atom, callable] = None,
        delete_in_target: Union[List, callable] = None,
        delete_in_source: Union[List, callable] = None,
        bond_order: int = None,
    ):
        """
        Create a new Reaction with modified parameters.

        Parameters
        ----------
        atom1 : Union[Atom, callable], optional
            The atom in the target molecule to which the source molecule will be connected.
            This can be an Atom object or a callable that takes a Molecule and returns an Atom
        atom2 : Union[Atom, callable], optional
            The atom in the source molecule which will be connected to the target molecule.
            This can be an Atom object or a callable that takes a Molecule and returns an Atom
        delete_in_target : Union[List, callable], optional
            A list of atoms in the target molecule to be deleted upon connection,
            or a callable that takes the atom1 and the target molecule and returns such a list.
            Default Hydrogen-deletion is applied if None.
        delete_in_source : Union[List, callable], optional
            A list of atoms in the source molecule to be deleted upon connection,
            or a callable that takes the atom2 and the source molecule and returns such a list.
            Default Hydrogen-deletion is applied if None.
        bond_order : int, optional
            The bond order of the new bond formed between atom1 and atom2. Default is 1 (single bond).

        Returns
        -------
        Reaction
            A new Reaction object with the modified parameters.
        """
        from copy import deepcopy

        new_reaction = deepcopy(self)
        new_reaction.set_reactivity(
            atom1=atom1,
            atom2=atom2,
            delete_in_target=delete_in_target,
            delete_in_source=delete_in_source,
            bond_order=bond_order,
        )
        return new_reaction

    def _find_deletes_for_anchors(self, target, source, atom1, atom2):

        delete_in_target = None
        delete_in_source = None

        if self._delete_in_target is not None:
            if callable(self._delete_in_target):
                delete_in_target = self._delete_in_target(atom1, target)
            else:
                delete_in_target = self._delete_in_target
            if not isinstance(delete_in_target, (list, tuple, set)):
                delete_in_target = [delete_in_target]
            for atom in delete_in_target:
                if not target.get_atom(atom):
                    return None

        if self._delete_in_source is not None:
            if callable(self._delete_in_source):
                delete_in_source = self._delete_in_source(atom2, source)
            else:
                delete_in_source = self._delete_in_source
            if not isinstance(delete_in_source, (list, tuple, set)):
                delete_in_source = [delete_in_source]
            for atom in delete_in_source:
                if not source.get_atom(atom):
                    return None

        return delete_in_target, delete_in_source

    def _apply_atom_getter(self, atom_getter, molecule):
        if callable(atom_getter):
            atom = atom_getter(molecule)
            if isinstance(atom, (list, tuple, set)) and len(atom) == 1:
                atom = set(atom).pop()
                return (molecule.get_atom(atom),)
            return molecule.get_atoms(atom)
        else:
            atom = atom_getter
            return molecule.get_atoms(atom)

    def _check_delete_callable_signature(self, func):
        if not callable(func):
            return func
        sig = inspect.signature(func)
        params = sig.parameters
        if len(params) > 2:
            raise ValueError(
                "Callable for delete_in_target or delete_in_source must accept at most two arguments: the atom and the molecule."
            )
        if len(params) == 2:
            return func
        elif len(params) == 1:

            def wrapper(atom, molecule):
                return func(atom)

            return wrapper
        elif len(params) == 0:

            def wrapper(atom, molecule):
                return func()

            return wrapper

    def __repr__(self):
        return f"Reaction(atom1={self._atom1}, atom2={self._atom2}, delete_in_target={self._delete_in_target}, delete_in_source={self._delete_in_source}, bond_order={self._bond_order})"

    def __call__(self, target, source, inplace=False):
        return self.apply(target, source, inplace)


if __name__ == "__main__":

    from buildamol.structural import constraints_v2 as constraints

    mol = Molecule.from_smiles("OCCO")
    bnz = Molecule.from_smiles("c1ccccc1C")

    R = Reaction(
        atom1=lambda mol: mol.get_atoms(
            "C", by="element", filter=constraints.has_single_bond_with("O")
        ).pop(),
        atom2=lambda mol: 1,
        delete_in_target=lambda atom: atom.get_neighbors(
            filter=constraints.has_element("O")
        ),
        bond_order=1,
    )
    new_mol = R.apply(mol, bnz)
    new_mol.show2d()

__all__ = ["Reaction"]
