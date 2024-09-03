"""
The `Assembler` class is a class that can be to assemble molecules from a library of fragments.

It requires a list of Molecules that serve as the fragments to be assembled. The class will then
generate random molecules by randomly selecting fragments from the library and attaching them to each other.

Usage
-----
1. Create a list of fragments to be used for assembly
2. Create an instance of the Assembler class with the list of fragments
3. Use the `sample` method to generate random molecules or the `make` method to create a specific fragment from an instruction matrix

Example
-------
Let's make a little toy example

.. code-block:: python

    import buildamol as bam
    from buildamol.extensions.molecular_factories import Assembler
    import matplotlib.pyplot as plt

    # get some molecules to serve as fragments
    fragments = [
        bam.Molecule.from_smiles("C1=CC=CC=C1", id="A").autolabel(),
        bam.Molecule.from_smiles("CC=O", id="B").autolabel(),
        bam.Molecule.from_smiles("COC=C", id="C").autolabel(),
        bam.Molecule.from_smiles("C1=CCC=C1", id="D").autolabel(),
        bam.Molecule.from_smiles("C(C)N", id="E").autolabel(),
    ]

    # make the assembler
    assembler = Assembler(fragments)


    # generate some molecules from 3 fragments each
    # let's make 9 molecules
    molecules = assembler.sample(n_fragments=3, n=9)

    fig, axs = plt.subplots(3, 3, figsize=(12, 12))
    for mol, ax in zip(molecules, axs.flat):
            ax.imshow(
                mol.draw2d().draw(),
            )
            ax.axis("off")    
    plt.show()

.. image:: examples/files/assembler_example1.png


Making Molecules from Arrays
----------------------------

We can also use the `make` method to create a specific molecule from an instruction matrix
This matrix is a 2D numpy array where each row corresponds to an instruction for attaching the next
fragment onto the molecule. The columns are as follows:

.. code-block::

    [
    [incoming_fragment_global_index, incoming_atom_index, target_fragment_atom],
    [incoming_fragment_global_index, incoming_atom_index, target_fragment_atom]
    ...
    ]


The `incoming_fragment_global_index` is the index of the fragment in the fragment library (i.e. in the list).
The `incoming_atom_index` is the index of the atom in the incoming fragment that will be attached to the target fragment (i.e. the attachment point).
The `target_fragment_atom` is the index of the atom in the target fragment that will be attached to the incoming fragment.

Let's make a molecule from an instruction matrix. Let's take the fourth fragment molecule as a start. Then attach the second fragment molecule to it, by attaching the its second atom to the first atom of already present molecule.
Then attach again the fourth fragment onto the molecule by attaching its first atom to the first atom of the second fragment in the molecule.

.. code-block:: python

    matrix = np.array([
    [3, 0, 0],
    [1, 1, 0],
    [3, 0, 0],
    ])

    mol = assembler.make(matrix)
    mol.draw2d().show()
    
.. image:: examples/files/assembler_example2.jpg

If including this into an automatic pipeline or an optimization loop it is recommended to wrap the whole thing into a try-except block to catch any errors that might occur due to invalid matrices.
The clue is that the atoms used for attachment should not be used more than once in the matrix. If they are used more than once, the molecule will not be able to be assembled leading to an error.
"""

import numpy as np
from buildamol.core import linkage, Molecule


class Assembler:
    """
    The Assembler class is a class that can be to assemble molecules from a library of fragments.
    Each molecule is a linear chain of fragments that are attached to each other.

    Parameters
    ----------
    fragments : list
        A list of Molecules that serve as the fragments to be assembled.
    """

    def __init__(self, fragments: list):

        # we need to maintain a per-fragment database of possible atom-sites where another fragment can be attached
        # we also need to maintain a per-fragment database of atom-ids to make linkages
        attachment_points = []
        atom_ids = []

        # let's browse through all fragments and identify the attachment points
        # also, filter out any fragments without any attachment points (good practice)
        to_drop = []
        for fdx, fragment in enumerate(fragments):

            # we define all non-Hydrogen atoms as potential attachment points
            # but only those that have a hydrogen neighbor that can be removed
            # will be considered as attachment points
            # n_atoms = sum(1 for i in fragment.get_atoms() if i.element != "H")
            a = []
            for adx, atom in enumerate(fragment.get_atoms()):
                if atom.element == "H":
                    continue
                if fragment.get_hydrogen(atom):
                    a.append(adx)
            if len(a) == 0:
                to_drop.append(fdx)
                continue
            attachment_points.append(a)
            atom_ids.append([atom.id for atom in fragment.get_atoms()])

        for fragment in to_drop:
            del fragments[fragment]

        self.fragments = fragments
        self.attachment_points = attachment_points
        self.atom_ids = atom_ids

    def specify_attachment_points(self, fragment_or_index, points: list):
        """
        Specify the attachment points for a fragment

        Parameters
        ----------
        fragment_or_index : int or Molecule
            The fragment for which to specify the attachment points
        points : list
            The attachment points to specify. These must be the indices of the atoms in the fragment as they appear in `fragment.get_atoms()` (NOT the `serial_number`!).
        """
        if isinstance(fragment_or_index, int):
            self.attachment_points[fragment_or_index] = points
        else:
            idx = self.fragments.index(fragment_or_index)
            self.attachment_points[idx] = points

    def sample(self, n_fragments: int, n: int = 1):
        """
        Generate n random molecules from the fragment library

        Parameters
        ----------
        n_fragments : int
            The number of fragments to use for each molecule
        n : int
            The number of molecules to generate

        Yields
        ------
        Molecule
            A molecule assembled from the fragments
        """
        for _ in range(n):
            matrix = self.random(n_fragments)
            yield self.make(matrix)

    def make(self, matrix: np.ndarray) -> Molecule:
        """
        Assemble a molecule based on an instruction matrix

        Parameters
        ----------
        matrix : np.ndarray
            The matrix encoding for the molecule

        Returns
        -------
        Molecule
            The assembled molecule
        """
        _used_atoms = {i: set() for i in range(len(matrix))}

        # we start by copying the first fragment
        mol = self.fragments[matrix[0, 0]].copy()

        # we then attach all other fragments
        for i in range(1, len(matrix)):
            source, source_atom, target_atom = matrix[i, :]
            target = i - 1
            # sanity checking to ensure we are not trying to attach to the same atom twice
            if target_atom in _used_atoms[target]:
                raise ValueError("Target atom already used")
            if source_atom in _used_atoms[i]:
                raise ValueError("Source atom already used")

            # make a linkage and attach the fragment
            link = linkage(
                self.atom_ids[matrix[target, 0]][target_atom],
                self.atom_ids[source][source_atom],
            )
            mol.attach(self.fragments[source], link, at_residue=int(target + 1))

            _used_atoms[target].add(target_atom)
            _used_atoms[i].add(source_atom)

        return mol

    def random(self, n_fragments: int) -> np.ndarray:
        """
        Make a random matrix encoding for a molecule assembled from fragments.

        Parameters
        ----------
        n_fragments : int
            The number of fragments to use for the molecule

        Returns
        -------
        np.ndarray
            A matrix encoding for the molecule
        """
        # we could literally just use a single line here of np.random here, but then we run the risk of
        # making invalid matrices where attachment_points are referenced more than once so we make a more intricate
        # method here to ensure our "random" matrices are valid

        matrix = np.full((n_fragments, 3), -1, dtype=int)
        matrix[0, 0] = np.random.choice(len(self.fragments))

        # we maintain a chache to keep track over which attachment points have been used already
        # on which fragments
        _used_atoms = {i: set() for i in range(n_fragments)}
        for i in range(1, n_fragments):

            # choose an incoming fragment from the database
            matrix[i, 0] = np.random.choice(len(self.fragments))

            # choose an attachment point on the incoming fragment that was not used already
            while matrix[i, 1] == -1:
                atom = np.random.choice(self.attachment_points[matrix[i, 0]])
                if atom not in _used_atoms[i]:
                    matrix[i, 1] = atom
                    _used_atoms[i].add(atom)

            # choose a target fragment in the molecule
            # and choose an attachment point in the target that was not used already
            target = i - 1
            while matrix[i, 2] == -1:
                available = [
                    i
                    for i in self.attachment_points[matrix[target, 0]]
                    if i not in _used_atoms[target]
                ]
                if len(available) > 0:
                    matrix[i, 2] = np.random.choice(available)
                    _used_atoms[target].add(matrix[i, 2])

        return matrix


if __name__ == "__main__":
    import buildamol as bam
    import matplotlib.pyplot as plt

    fragments = [
        bam.Molecule.from_smiles("C1=CC=CC=C1", id="A").autolabel(),
        bam.Molecule.from_smiles("CC=O", id="B").autolabel(),
        bam.Molecule.from_smiles("COC=C", id="C").autolabel(),
        bam.Molecule.from_smiles("C1=CCC=C1", id="D").autolabel(),
        bam.Molecule.from_smiles("C(C)N", id="E").autolabel(),
    ]

    assembler = Assembler(fragments)
    matrix = np.array(
        [
            [3, 0, 0],
            [1, 1, 0],
            [3, 0, 0],
        ]
    )

    mol = assembler.make(matrix)
    mol.draw2d().show()
    fig, axs = plt.subplots(3, 3, figsize=(12, 12))

    for mol, ax in zip(assembler.sample(3, 9), axs.flat):
        ax.imshow(
            mol.draw2d().draw(),
        )
        ax.axis("off")

    plt.show()
