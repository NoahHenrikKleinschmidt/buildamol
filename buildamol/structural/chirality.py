"""
Chiral center analysis using CIP (Cahn-Ingold-Prelog) priority rules.
"""

import numpy as np


def _cip_priority_key(atom, visited=frozenset(), depth=6):
    """
    Build a CIP priority key for *atom* as a flat tuple of atomic numbers.

    The tuple is constructed by recursively listing, at each bonding sphere,
    the sorted (descending) atomic numbers of all reachable atoms.  Comparison
    of two such tuples is lexicographic, which correctly ranks substituents
    according to CIP rules 1-2 for acyclic systems.  Ring systems are handled
    approximately: once an atom is in *visited* its branch is cut off, which
    corresponds to the CIP "phantom atom" convention (atomic number retained,
    no further expansion).

    Parameters
    ----------
    atom :
        The atom whose CIP key is being built.
    visited : frozenset
        Atoms already seen on the path from the chiral center to *atom*
        (used to avoid cycling through rings).
    depth : int
        Maximum recursion depth (shells away from the chiral center).
    """
    if depth == 0 or atom in visited:
        return (atom.atomic_number,)

    new_visited = visited | {atom}
    neighbors = [n for n in atom.get_neighbors(n=1, mode="at") if n not in visited]
    sub_keys = sorted(
        (_cip_priority_key(n, new_visited, depth - 1) for n in neighbors),
        reverse=True,
    )
    result = [atom.atomic_number]
    for sk in sub_keys:
        result.extend(sk)
    return tuple(result)


class ChiralCenter:
    """
    A tetrahedral chiral center defined by a central atom and its four substituents.

    Substituents are labelled A, B, C, D in **descending** CIP priority order
    (A = highest, D = lowest).  The R/S configuration is determined from the
    3-D coordinates stored in the molecule.

    Parameters
    ----------
    atom :
        The central atom of the chiral center.  It must be part of a
        :class:`~buildamol.core.Molecule` and must have exactly four bonded
        neighbours.

    Raises
    ------
    ValueError
        If the atom is not part of a molecule, or does not have exactly four
        bonded neighbours.

    Examples
    --------
    >>> cc = ChiralCenter(mol.get_atom("CA"))
    >>> cc.orientation
    'S'
    >>> cc.A, cc.B, cc.C, cc.D
    (<Atom N ...>, <Atom C ...>, <Atom CB ...>, <Atom H ...>)
    >>> cc.is_S()
    True
    """

    def __init__(self, atom):
        self._atom = atom
        self._A = None
        self._B = None
        self._C = None
        self._D = None
        self._assign_priorities()

    # ------------------------------------------------------------------
    # Internal priority assignment
    # ------------------------------------------------------------------

    def _assign_priorities(self):
        if self._atom.molecule is None:
            raise ValueError(
                f"Atom {self._atom} is not part of a Molecule; "
                "chirality requires connectivity information."
            )
        neighbors = list(self._atom.get_neighbors(n=1, mode="at"))
        if len(neighbors) != 4:
            raise ValueError(
                f"A chiral center requires exactly 4 bonded neighbours, "
                f"but atom {self._atom} has {len(neighbors)}."
            )

        # Exclude the central atom from all recursive CIP traversals
        visited = frozenset({self._atom})
        ranked = sorted(
            neighbors,
            key=lambda a: _cip_priority_key(a, visited),
            reverse=True,
        )
        self._A, self._B, self._C, self._D = ranked

    # ------------------------------------------------------------------
    # Substituent properties
    # ------------------------------------------------------------------

    @property
    def A(self):
        """Highest-priority substituent (CIP rank 1)."""
        return self._A

    @property
    def B(self):
        """Second-priority substituent (CIP rank 2)."""
        return self._B

    @property
    def C(self):
        """Third-priority substituent (CIP rank 3)."""
        return self._C

    @property
    def D(self):
        """Lowest-priority substituent (CIP rank 4)."""
        return self._D

    # ------------------------------------------------------------------
    # Stereodescriptor
    # ------------------------------------------------------------------

    @property
    def orientation(self) -> str:
        """
        CIP stereodescriptor: ``'R'`` or ``'S'``.

        Determined geometrically: viewing the chiral center from the side
        opposite to **D** (lowest priority), a clockwise sequence A → B → C
        corresponds to *R*; counterclockwise corresponds to *S*.

        Raises
        ------
        ValueError
            If the four substituents are coplanar (degenerate geometry).
        """
        return self._compute_orientation()

    def is_R(self) -> bool:
        """Return ``True`` if this center has *R* configuration."""
        return self.orientation == "R"

    def is_S(self) -> bool:
        """Return ``True`` if this center has *S* configuration."""
        return self.orientation == "S"

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def _compute_orientation(self) -> str:
        center = np.array(self._atom.get_coord(), dtype=float)
        a_vec = np.array(self._A.get_coord(), dtype=float) - center
        b_vec = np.array(self._B.get_coord(), dtype=float) - center
        c_vec = np.array(self._C.get_coord(), dtype=float) - center
        d_vec = np.array(self._D.get_coord(), dtype=float) - center

        # Normal to the A→B→C plane via the right-hand rule.
        # When this normal points in the same half-space as D, A→B→C appears
        # counter-clockwise to an observer sitting opposite D  → S.
        # When it points away from D, A→B→C is clockwise → R.
        normal = np.cross(b_vec - a_vec, c_vec - a_vec)
        dot = np.dot(normal, d_vec)

        if dot > 0:
            return "S"
        elif dot < 0:
            return "R"
        else:
            raise ValueError(
                "Cannot determine chirality: substituents are coplanar or two "
                "atoms share the same coordinates."
            )

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        try:
            ori = self.orientation
        except Exception:
            ori = "?"
        return (
            f"ChiralCenter(atom={self._atom!r}, orientation={ori!r}, "
            f"A={self._A!r}, B={self._B!r}, C={self._C!r}, D={self._D!r})"
        )
