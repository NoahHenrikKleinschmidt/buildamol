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


def _resolve_index(atom_or_idx, fragment):
    """
    Normalise an attachment/deletion specifier to an integer index into
    ``list(fragment.get_atoms())``.

    Accepts either a plain ``int`` (returned unchanged) or an ``Atom`` object
    whose ``molecule`` attribute is used to derive the position.
    """
    if isinstance(atom_or_idx, int):
        return atom_or_idx
    atoms = list(atom_or_idx.molecule.get_atoms())
    return atoms.index(atom_or_idx)


def _parallel_score(molecules, scoring_fn, n_workers):
    """Score a list of molecules, in parallel when n_workers > 1."""

    def _score_one(mol):
        try:
            return {"score": float(scoring_fn(mol)), "molecule": mol}
        except Exception:
            return None

    if n_workers <= 1:
        results = [_score_one(m) for m in molecules]
    else:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_workers) as exe:
            results = list(exe.map(_score_one, molecules))

    results = [r for r in results if r is not None]
    results.sort(key=lambda r: r["score"], reverse=True)
    return results


class Assembler:
    """
    The Assembler class is a class that can be to assemble molecules from a library of fragments.
    Each molecule is a linear chain of fragments that are attached to each other.

    Parameters
    ----------
    fragments : list
        A list of Molecules that serve as the fragments to be assembled.
    """

    def __init__(self, fragments: list, n_workers: int = 1, optimization_steps: int = 30, bystander_radius: float = 8):

        # we need to maintain a per-fragment database of possible atom-sites where another fragment can be attached
        # we also need to maintain a per-fragment database of atom-ids to make linkages
        attachment_points = []
        # per-fragment dict: attachment_atom_idx -> [deletable_atom_idx, ...]
        deletion_points = []
        atom_ids = []

        # let's browse through all fragments and identify the attachment points
        # also, filter out any fragments without any attachment points (good practice)
        to_drop = []
        for fdx, fragment in enumerate(fragments):

            # we define all non-Hydrogen atoms as potential attachment points
            # but only those that have at least one deletable neighbor (defaulting to
            # hydrogen neighbors) will be considered as attachment points
            atoms_list = list(fragment.get_atoms())
            a = []  # attachment point atom indices
            dp = {}  # deletion_points dict for this fragment

            for adx, atom in enumerate(atoms_list):
                if atom.element == "H":
                    continue
                h_neighbors = fragment.get_hydrogens(atom)
                h_indices = [i for i, at in enumerate(atoms_list) if at in h_neighbors]
                if h_indices:
                    a.append(adx)
                    dp[adx] = h_indices

            if len(a) == 0:
                to_drop.append(fdx)
                continue
            attachment_points.append(a)
            deletion_points.append(dp)
            atom_ids.append([atom.id for atom in atoms_list])

        for fragment in to_drop:
            del fragments[fragment]

        self.fragments = fragments
        self.attachment_points = attachment_points
        self.deletion_points = deletion_points
        self.atom_ids = atom_ids
        self.n_workers = n_workers
        self.optimization_steps = optimization_steps
        self.bystander_radius = bystander_radius

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _resolve_fragment_index(self, fragment_or_index) -> int:
        if isinstance(fragment_or_index, int):
            return fragment_or_index
        return self.fragments.index(fragment_or_index)

    # ------------------------------------------------------------------
    # Public configuration API
    # ------------------------------------------------------------------

    def specify_attachment_points(self, fragment_or_index, points: list):
        """
        Specify the attachment points for a fragment.

        Parameters
        ----------
        fragment_or_index : int or Molecule
            The fragment for which to specify the attachment points.
        points : list
            The attachment points to specify.  Each element may be either an
            ``int`` index into ``fragment.get_atoms()`` (NOT the
            ``serial_number``!) or a literal ``Atom`` object that belongs to
            the fragment.
        """
        fdx = self._resolve_fragment_index(fragment_or_index)
        frag = self.fragments[fdx]
        self.attachment_points[fdx] = [_resolve_index(p, frag) for p in points]

    def specify_deletion_points(self, fragment_or_index, deletion_points: list):
        """
        Specify which atoms to delete for each attachment point of a fragment.

        This overrides the default behaviour of automatically selecting a
        hydrogen neighbor.  Use this when you need deterministic stereo-
        chemistry control or when non-hydrogen atoms should be removed at
        specific attachment sites.

        Parameters
        ----------
        fragment_or_index : int or Molecule
            The fragment (by library index or by object reference).
        deletion_points : list of list
            A list of lists that is **index-matched** to the current
            ``attachment_points`` of the fragment.  Each inner list contains
            the atom(s) that may be deleted when the corresponding attachment
            point is used; one is chosen at random each time the point is
            activated.  Accepted element types per inner list:

            * ``int`` — index into ``list(fragment.get_atoms())``
            * ``Atom`` object belonging to the fragment

            Pass an empty inner list ``[]`` for a given attachment point to
            fall back to automatic hydrogen deletion for that site.

        Raises
        ------
        ValueError
            If ``deletion_points`` does not have the same length as the
            fragment's current ``attachment_points``.

        Examples
        --------
        .. code-block:: python

            assembler = Assembler([frag_a, frag_b])

            # Override deletion for frag_a (has 2 attachment points):
            # - attachment point 0: always delete atom index 5
            # - attachment point 1: choose randomly between atom 8 or 9
            assembler.specify_deletion_points(0, [[5], [8, 9]])

            # Same thing using Atom objects directly:
            atoms = list(frag_a.get_atoms())
            assembler.specify_deletion_points(frag_a, [[atoms[5]], [atoms[8], atoms[9]]])
        """
        fdx = self._resolve_fragment_index(fragment_or_index)
        frag = self.fragments[fdx]
        att_pts = self.attachment_points[fdx]

        if len(deletion_points) != len(att_pts):
            raise ValueError(
                f"deletion_points must be index-matched to attachment_points "
                f"(expected {len(att_pts)} inner lists, got {len(deletion_points)})."
            )

        dp = {}
        for att_idx, candidates in zip(att_pts, deletion_points):
            dp[att_idx] = [_resolve_index(c, frag) for c in candidates]
        self.deletion_points[fdx] = dp

    # ------------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------------

    def sample(self, n_fragments: int, n: int = 1):
        """
        Generate n random molecules from the fragment library.

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
        if self.n_workers <= 1:
            for _ in range(n):
                try:
                    matrix = self.random(n_fragments)
                    yield self.make(matrix)
                except Exception:
                    pass
        else:
            from concurrent.futures import ThreadPoolExecutor

            matrices = []
            for _ in range(n):
                try:
                    matrices.append(self.random(n_fragments))
                except Exception:
                    pass

            def _safe_make(m):
                try:
                    return self.make(m)
                except Exception:
                    return None

            with ThreadPoolExecutor(max_workers=self.n_workers) as exe:
                for mol in exe.map(_safe_make, matrices):
                    if mol is not None:
                        yield mol

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

            target_frag = int(matrix[target, 0])
            source_frag = int(source)
            target_atom_i = int(target_atom)
            source_atom_i = int(source_atom)

            # Determine which atoms to delete at each end of the new bond.
            # If deletion_points has candidates for this attachment atom, sample
            # one; otherwise pass None and let the linkage auto-delete a hydrogen.
            t_candidates = self.deletion_points[target_frag].get(target_atom_i, [])
            s_candidates = self.deletion_points[source_frag].get(source_atom_i, [])

            delete_in_target = (
                [self.atom_ids[target_frag][np.random.choice(t_candidates)]]
                if t_candidates
                else None
            )
            delete_in_source = (
                [self.atom_ids[source_frag][np.random.choice(s_candidates)]]
                if s_candidates
                else None
            )

            link = linkage(
                self.atom_ids[target_frag][target_atom_i],
                self.atom_ids[source_frag][source_atom_i],
                delete_in_target=delete_in_target,
                delete_in_source=delete_in_source,
            )
            mol.attach(self.fragments[source_frag], link, at_residue=int(target + 1), optimization_steps=self.optimization_steps, bystander_radius=self.bystander_radius)

            _used_atoms[target].add(target_atom)
            _used_atoms[i].add(source_atom)

        return mol

    def score_sample(self, scoring_fn, n_fragments: int, n: int = 1) -> list:
        """
        Generate ``n`` random molecules and score them, returning a list of
        ``{"score": float, "molecule": Molecule}`` dicts sorted by score.

        Molecule assembly stays serial (thread safety), while the scoring
        function is evaluated in parallel when ``n_workers > 1``.

        Parameters
        ----------
        scoring_fn : callable
            ``scoring_fn(mol) -> float``.
        n_fragments : int
            Number of fragments per assembled molecule.
        n : int
            Number of molecules to generate and score.

        Returns
        -------
        list of dict
            ``[{"score": float, "molecule": Molecule}, ...]``, sorted best-first.
        """
        molecules = list(self.sample(n_fragments, n))
        return _parallel_score(molecules, scoring_fn, self.n_workers)

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
                    a
                    for a in self.attachment_points[matrix[target, 0]]
                    if a not in _used_atoms[target]
                ]
                if available:
                    matrix[i, 2] = np.random.choice(available)
                    _used_atoms[target].add(matrix[i, 2])
                else:
                    raise ValueError(
                        f"Fragment at chain position {target} has no remaining "
                        f"attachment points for the next fragment. Consider using "
                        f"n_fragments <= number of attachment points on each fragment."
                    )

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
