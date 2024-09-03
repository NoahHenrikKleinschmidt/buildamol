"""
The `Derivator` class can be used to create stochastic and/or combinatorical derivatives of a molecule.

Usage
-----
1. Create a `Derivator` object with a molecule as an argument.
2. Use the `element_changable`, `functional_group_addable`, `bond_order_changable`, and `modifiable` methods to specify the possible changes that can be made to the molecule.
3. Use the `sample`, `all`, or `make` method to generate one or more derivative molecule(s).

Example
-------
Let's say we have a benzene ring which we want to derivatize. Let's say further that we are interested in changing the elements of two atoms and adding one functional group to an atom.
We can do this as follows:

1. Create the benzene molecule and derivator object

.. code-block:: python

    import buildamol as bam
    bam.load_small_molecules()
    benzene = bam.molecule("benzene")

    # Create the derivator object
    D = bam.Derivator(benzene)

2. Specify the possible changes that can be made to the molecule

Let's say that `C1` can be either a carbon, a nitrogen, or an oxygen.
And `C5` can be either a carbon or a nitrogen.
Also, we are interested in possibly adding a functional group to `C3`. 
Namely, we want to check with a hydroxyl group and a phosphate group.

.. code-block:: python

    # set the possible elements for C1 and C5
    D.element_changable(
        "C1",
        ("C", "N", "O"),
        )
    D.element_changable(
        "C5",
        ("C", "N"),
        )

    # set the possible functional groups for C3
    # (the groups are functions that modify the molecule in place)
    D.functional_group_addable(
        "C3",
        (D.nothing, bam.hydroxylate, bam.phosphorylate),
        )


3. Generate derivatives

Since we only have a small number of changes, we can generate all possible derivatives using the `all` method.

.. code-block:: python

    # generate all possible derivatives
    derivatives = list(D.all())


4. Visualize the derivatives

Now we can visualize the derivatives using the `draw2d` method.

.. code-block:: python

    import matplotlib.pyplot as plt

    # compute the layout for the plots (looks nicer like this)
    rows = np.sqrt(len(derivatives))
    cols = np.ceil(len(derivatives) / rows)
    rows = np.round(rows)

    fig, axs = plt.subplots(int(rows), int(cols), figsize=(12, 12))
    for i, molecule in enumerate(derivatives):
        molecule.squash()
        try:
            img = molecule.draw2d().draw()
            axs.flat[i].imshow(img)
        except Exception as e:
            print(e)
    for ax in axs.flat:
        ax.axis("off")

    plt.show()

.. image:: examples/files/derivator_example1.png

"""

from buildamol.core import Molecule, Atom, Bond
import numpy as np
import itertools


class Derivator:
    """
    A class to generate derivatives of a molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to generate derivatives of
    """

    def __init__(self, molecule: "Molecule"):
        self.molecule = molecule
        self._element_derivables = {}
        self._functional_group_derivables = {}
        self._bond_derivables = {}
        self._modifiers = (tuple(), tuple())

    @property
    def N(self) -> int:
        """
        The number of possible derivatives
        """
        n = 1
        for elements, _ in self._element_derivables.values():
            n *= len(elements)
        for orders, _ in self._bond_derivables.values():
            n *= len(orders)
        for groups, _ in self._functional_group_derivables.values():
            n *= len(groups)
        if len(self._modifiers[0]) > 0:
            n *= len(self._modifiers[0])
        return n

    def make(
        self, elements: tuple, bonds: tuple, groups: tuple, modifiers: tuple
    ) -> Molecule:
        """
        Make a specific derivative molecule.
        Provide a tuple or other ordered iterable of the same length
        as the number of atoms/bonds specified during setup for each editable category.
        Each element must be an integer denoting the index of the desired derivative application.

        Parameters
        ----------
        elements : tuple
            The elements to apply.
        bonds : tuple
            The bond orders to apply.
        groups : tuple
            The group modifiers to apply.
        modifiers : tuple
            The molecule modifiers to apply.

        Returns
        -------
        Molecule
            The generated derivative molecule
        """
        molecule = self.molecule.copy()
        element_changeble_atoms = list(self._element_derivables.keys())
        for idx, element in enumerate(elements):
            atom = element_changeble_atoms[idx]
            element = self._element_derivables[atom][0][element]
            molecule.change_element(molecule.get_atom(atom), element)

        order_changeble_bonds = list(self._bond_derivables.keys())
        for idx, order in enumerate(bonds):
            bond = order_changeble_bonds[idx]
            order = self._bond_derivables[bond][0][order]
            molecule.set_bond_order(*bond, order)

        group_addable_atoms = list(self._functional_group_derivables.keys())
        for idx, group in enumerate(groups):
            atom = group_addable_atoms[idx]
            group = self._functional_group_derivables[atom][0][group]
            atom = molecule.get_atom(atom)
            group(molecule, atom)

        for modifier in modifiers:
            m = self._modifiers[0][modifier]
            m(molecule)
        return molecule

    def sample(
        self,
        n: int = 1,
        use_probabilities: bool = True,
        apply_all_modifiers: bool = False,
    ) -> "Generator[Molecule]":
        """
        Create a derivative molecule

        Parameters
        ----------
        n : int
            The number of derivatives to generate
        use_probabilities : bool
            Whether to use the probabilities when generating derivatives
        apply_all_modifiers : bool
            If True all modifiers are applied to the sampled molecule(s). Otherwise
            one modifier is sampled and applied to each molecule.

        Yields
        ------
        Molecule
            The generated derivative molecule
        """
        for i in range(n):
            yield self._derive(use_probabilities, apply_all_modifiers)

    def all(self) -> "Generator[Molecule]":
        """
        Generate all possible derivatives

        Yields
        ------
        Molecule
            The generated derivative molecule
        """
        has_element_derivables = len(self._element_derivables) > 0
        has_bond_derivables = len(self._bond_derivables) > 0
        has_group_derivables = len(self._functional_group_derivables) > 0
        has_global_modifiers = len(self._modifiers[0]) > 0

        element_is_endpoint = has_element_derivables and not (
            has_bond_derivables or has_group_derivables or has_global_modifiers
        )
        bond_is_endpoint = has_bond_derivables and not (
            has_group_derivables or has_global_modifiers
        )
        group_is_endpoint = has_group_derivables and not has_global_modifiers

        already_made = set()
        molecules = []

        if has_element_derivables:
            incoming = []
            for molecule in self._iterate_all_elements(self.molecule):
                smiles = molecule.to_smiles()
                if smiles not in already_made:
                    already_made.add(smiles)
                    incoming.append(molecule)

            if element_is_endpoint:
                for molecule in incoming:
                    yield molecule
            else:
                molecules = incoming

        if has_bond_derivables:
            incoming = []
            for i in range(len(molecules)):
                for molecule in self._iterate_all_bonds(molecules[i]):
                    smiles = molecule.to_smiles()
                    if smiles not in already_made:
                        already_made.add(smiles)
                        incoming.append(molecule)

            if bond_is_endpoint:
                for molecule in incoming:
                    yield molecule
            else:
                molecules = incoming

        if has_group_derivables:
            incoming = []
            for i in range(len(molecules)):
                for molecule in self._iterate_all_groups(molecules[i]):
                    smiles = molecule.to_smiles()
                    if smiles not in already_made:
                        already_made.add(smiles)
                        incoming.append(molecule)

            if group_is_endpoint:
                for molecule in incoming:
                    yield molecule
            else:
                molecules = incoming

        if has_global_modifiers:
            for i in range(len(molecules)):
                for molecule in self._iterate_all_modifiers(molecules[i]):
                    smiles = molecule.to_smiles()
                    if smiles not in already_made:
                        already_made.add(smiles)
                        molecules.append(molecule)
                        yield molecule

    def element_changable(self, atom, elements: tuple, probabilities: tuple = None):
        """
        Specify an atom to be derivable to a set of possible elements

        Parameters
        ----------
        atom : Atom
            The atom
        elements : set
            The set of possible elements
        probabilities : tuple
            The probabilities of each element being used
        """
        if not isinstance(atom, Atom):
            _atom = self.molecule.get_atom(atom)
            if _atom is None:
                raise ValueError(f"Atom {atom} not found in molecule")
            atom = _atom
        if probabilities is None:
            probabilities = [1 / len(elements)] * len(elements)
        else:
            if len(elements) != len(probabilities):
                raise ValueError(
                    "The number of elements must match the number of probabilities"
                )
            probabilities = np.array(probabilities) / np.sum(probabilities)

        self._element_derivables[atom] = (
            tuple(elements),
            probabilities,
        )
        return self

    def functional_group_addable(
        self, atom, group_modifiers: tuple, probabilities: tuple = None
    ):
        """
        Specify an atom onto which a functional group can be added

        Parameters
        ----------
        atom : Atom
            The atom
        group_modifiers : tuple
            The set of possible functional groups that can be added onto the atom.
            Each entry must be a function that will modify the molecule in place, taking
            only the molecule and target atom as arguments!
        probabilities : tuple
            The probabilities of each functional group being added
        """
        if not isinstance(atom, Atom):
            _atom = self.molecule.get_atom(atom)
            if _atom is None:
                raise ValueError(f"Atom {atom} not found in molecule")
            atom = _atom
        if probabilities is None:
            probabilities = [1 / len(group_modifiers)] * len(group_modifiers)
        else:
            if len(group_modifiers) != len(probabilities):
                raise ValueError(
                    "The number of group modifiers must match the number of probabilities"
                )
            probabilities = np.array(probabilities) / np.sum(probabilities)

        self._functional_group_derivables[atom] = (
            tuple(group_modifiers),
            probabilities,
        )
        return self

    def bond_order_changable(
        self, bond, bond_orders: tuple, probabilities: tuple = None
    ):
        """
        Specify a bond to be derivable to a set of possible bond orders

        Parameters
        ----------
        bond : Bond
            The bond
        bond_orders : set
            The set of possible bond orders
        probabilities : tuple
            The probabilities of each bond order being used
        """
        if not isinstance(bond, Bond):
            _bond = self.molecule.get_bond(*bond)
            if _bond is None:
                raise ValueError(f"Bond {bond} not found in molecule")
            bond = _bond
        if probabilities is None:
            probabilities = [1 / len(bond_orders)] * len(bond_orders)
        else:
            if len(bond_orders) != len(probabilities):
                raise ValueError(
                    "The number of bond orders must match the number of probabilities"
                )
            probabilities = np.array(probabilities) / np.sum(probabilities)

        self._bond_derivables[bond] = (
            tuple(bond_orders),
            probabilities,
        )
        return self

    def global_modifiers(self, modifiers: tuple, probabilities: tuple = None):
        """
        Specify a set of possible modifications to the molecule as a whole.

        Parameters
        ----------
        modifiers : tuple
            The set of possible modifications that can be applied to the molecule.
            Each entry must be a function that will modify the molecule in place, taking
            only the molecule as an argument!
        probabilities : tuple
            The probabilities of each modification being applied
        """
        if probabilities is None:
            probabilities = [1 / len(modifiers)] * len(modifiers)
        else:
            if len(modifiers) != len(probabilities):
                raise ValueError(
                    "The number of modifiers must match the number of probabilities"
                )
            probabilities = np.array(probabilities) / np.sum(probabilities)

        self._modifiers = (tuple(modifiers), probabilities)
        return self

    def nothing(self, *args):
        """
        A function to do nothing (can be passed as a functional group modifier)
        """
        pass

    def _iterate_all_elements(self, mol):

        atoms = list(self._element_derivables.keys())
        elements_list = [elements for elements, _ in self._element_derivables.values()]

        for combination in itertools.product(*elements_list):
            molecule = mol.copy()
            for atom, element in zip(atoms, combination):
                atom_obj = molecule.get_atom(atom)
                if element != atom_obj.element:
                    molecule.change_element(atom_obj, element)
            yield molecule

    def _iterate_all_bonds(self, mol):
        bonds = list(self._bond_derivables.keys())
        orders_list = [orders for orders, _ in self._bond_derivables.values()]

        for combination in itertools.product(*orders_list):
            molecule = mol.copy()
            for bond, order in zip(bonds, combination):
                if order != bond.order:
                    molecule.set_bond_order(*bond, order)
            yield molecule

    def _iterate_all_groups(self, mol):
        atoms = list(self._functional_group_derivables.keys())
        groups_list = [
            groups for groups, _ in self._functional_group_derivables.values()
        ]

        for combination in itertools.product(*groups_list):
            molecule = mol.copy()
            for atom, group in zip(atoms, combination):
                atom = molecule.get_atom(atom)
                group(molecule, atom)
            yield molecule

    def _iterate_all_modifiers(self, mol):
        for combination in itertools.product(*self._modifiers):
            molecule = mol.copy()
            for modifier in combination:
                modifier(molecule)
            yield molecule

    def _derive(
        self, use_probabilities: bool = False, apply_all_modifiers: bool = False
    ):
        """
        Derive the molecule
        """
        molecule = self.molecule.copy()
        for atom, (
            elements,
            probabilities,
        ) in self._element_derivables.items():
            if use_probabilities:
                element = np.random.choice(elements, p=probabilities)
            else:
                element = np.random.choice(elements)
            atom = molecule.get_atom(atom)
            if element != atom.element:
                molecule.change_element(atom, element)

        for bond, (orders, probabilities) in self._bond_derivables.items():
            if use_probabilities:
                order = np.random.choice(orders, p=probabilities)
            else:
                order = np.random.choice(orders)
            molecule.set_bond_order(*bond, order)

        for atom, (
            groups,
            probabilities,
        ) in self._functional_group_derivables.items():
            if use_probabilities:
                group = np.random.choice(groups, p=probabilities)
            else:
                group = np.random.choice(groups)

            atom = molecule.get_atom(atom)
            group(molecule, atom)

        if len(self._modifiers[0]) == 0:
            return molecule

        if apply_all_modifiers:
            for modifier in self._modifiers[0]:
                modifier(molecule)
            return molecule

        if use_probabilities:
            modifier = np.random.choice(self._modifiers[0], p=self._modifiers[1])
        else:
            modifier = np.random.choice(self._modifiers[0])
        modifier(molecule)
        return molecule


if __name__ == "__main__":

    import buildamol as bam

    bam.load_small_molecules()
    # start = bam.molecule("benzene")
    # start = bam.methylate(start, "C1")
    # start = bam.benzylate(start, "C")
    # start.reindex().autolabel()

    # D = Derivator(start)
    # do_nothing = lambda m, a: None

    # D.element_changable(
    #     start.get_atom("C1", residue=2),
    #     ("C", "N", "O"),
    #     (0.6, 0.2, 0.2),
    # )
    # D.element_changable(
    #     start.get_atom("C5", residue=3),
    #     ("C", "N", "O"),
    #     (0.6, 0.3, 0.4),
    # )

    # D.functional_group_addable(
    #     start.get_atom("C3", residue=1),
    #     (D.nothing, bam.hydroxylate, bam.acetylate, bam.amidate, bam.carboxylate),
    # )
    # D.functional_group_addable(
    #     start.get_atom("C4", residue=1),
    #     (D.nothing, bam.hydroxylate, bam.acetylate, bam.amidate, bam.carboxylate),
    # )
    # # D.functional_group_addable(
    # #     start.get_atom("C2", residue=3),
    # #     (bam.methylate, bam.hydroxylate, bam.acetylate, bam.amidate, bam.carboxylate),
    # # )

    D = Derivator(bam.molecule("benzene"))
    D.element_changable(
        "C1",
        ("N", "O"),
    )
    D.element_changable(
        "C5",
        ("C", "N"),
    )
    D.functional_group_addable(
        "C3",
        (bam.hydroxylate, bam.phosphorylate),
    )

    import matplotlib.pyplot as plt

    mols = list(D.all())

    rows = np.sqrt(len(mols))
    cols = np.ceil(len(mols) / rows)

    fig, axs = plt.subplots(int(np.round(rows)), int(cols), figsize=(12, 12))
    for i, molecule in enumerate(mols):
        molecule.squash()
        try:
            img = molecule.draw2d().draw()
            axs.flat[i].imshow(img)
        except Exception as e:
            print(e)
    for ax in axs.flat:
        ax.axis("off")
    plt.show()
