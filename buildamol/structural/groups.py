"""
Geometric descriptors for chemical Functional Groups
"""

from typing import Union
import numpy as np
from scipy.spatial.distance import cdist

import buildamol.structural.base as base
import buildamol.structural.infer as infer
import buildamol.structural.geometry as geometry
import buildamol.structural.neighbors as neighbors

constraints = neighbors.constraints

from copy import deepcopy

__H_neighbor_funcs = {}


def _H_neighbor_of_assigned_atom(n: int):
    if n in __H_neighbor_funcs:
        return __H_neighbor_funcs[n]
    __H_neighbor_funcs[n] = lambda mol, res, bonder, assignment: (
        assignment[n],
        mol.get_neighbors(assignment[n], filter=lambda x: x.element == "H").pop(),
    )
    return __H_neighbor_funcs[n]


class FunctionalGroup:
    """
    A functional group that can be described by a single geometry around one central atom.

    Parameters
    ----------
    id : str
        The identifier of the functional group.
    geometry : geometry.Geometry
        The geometry of the functional group.
    atoms : str
        The elements of the atoms. The first element is the center atom.
        Then come all other atoms.
    connectivity : list[int]
        The connectivity of the atoms. Where each entry describes the bond order of the central atom to the respective other atom.
    constraints : list
        A list of functions that describe the constraints of the functional group.
        Each function describes the constraints of the respective atom in the same order as they were provided in the atoms parameter.
    invertable : bool
        If the functional group is invertable.

    Examples
    --------
    >>> from buildamol.structural.neighbors import constraints
    >>> carboxyl = FunctionalGroup(
    ...     id="carboxyl",
    ...     rank=2,
    ...     geometry=geometry.trigonal_planar,
    ...     # the first atom is the center atom (carboxyl carbon)
    ...     atoms=("C", "O", "O"),
    ...     # the connectivity of the carbonxyl carbon to the first O is a double
    ...     # bond and to the second O is a single bond
    ...     connectivity=(2, 1),
    ...     constraints=[
    ...     None, # the carboxyl carbon has no constraints
    ...     constraints.neighbors_exactly("C"), # the first oxygen must be connected to the carboxyl carbon only
    ...     constraints.neighbors_exactly("C", "H") # the second oxygen must be connected to the carboxyl carbon and a hydrogen
    ...     ],
    ...     invertable=False
    ... )
    """

    def __init__(
        self,
        id: str,
        rank: int,
        geometry: geometry.Geometry,
        atoms: tuple[str],
        connectivity: list[int],
        constraints: list[tuple],
        invertable: bool = False,
    ):
        self.id = id
        self.geometry = geometry
        self.atoms = atoms
        self.connectivity = connectivity
        self._connectivity = [(0, i + 1) for i in range(len(connectivity))]
        self._assignment = None
        self._hist = self._element_hist(atoms)
        self._set_atoms = set(atoms)
        self.n = len(atoms)
        self.rank = rank
        self.constraints = constraints
        self.invertable = invertable
        self._nucleophile_bonder = None
        self._electrophile_bonder = None
        self._nucleophile_deletes = None
        self._electrophile_deletes = None

    def with_reactivity(
        self,
        nucleophile_bonder: int = None,
        electrophile_bonder: int = None,
        nucleophile_deletes: Union[list, callable] = None,
        electrophile_deletes: Union[list, callable] = None,
    ):
        """
        Create a copy of the FunctionalGroup with a new reactivity.
        Any reactivity that is not provided will be copied from the original FunctionalGroup.

        Parameters
        ----------
        nucleophile_bonder : int
            The index of the nucleophile bonder.
        electrophile_bonder : int
            The index of the electrophile bonder.
        nucleophile_deletes : list or callable
            The indices of atoms to delete if the functional group is used as a nucleophile.
            Or a function that will return a list of atoms to delete. This function will receive the following arguments: `molecule` (Molecule), `residue` (Residue), `bonder` (Atom), `assignment` (list[Atom] that were identified as belonging to the functional group).
        electrophile_deletes : list or callable
            The indices of atoms to delete if the functional group is used as an electrophile.
            Or a function that will return a list of atoms to delete. This function will receive the following arguments: `molecule` (Molecule), `residue` (Residue), `bonder` (Atom), `assignment` (list[Atom] that were identified as belonging to the functional group).

        Returns
        -------
        FunctionalGroup
            A new FunctionalGroup with the new reactivity.
        """
        new = deepcopy(self)
        new.set_reactivity(
            nucleophile_bonder,
            electrophile_bonder,
            nucleophile_deletes,
            electrophile_deletes,
        )
        return new

    def set_reactivity(
        self,
        nucleophile_bonder: int = None,
        electrophile_bonder: int = None,
        nucleophile_deletes: Union[list, callable] = None,
        electrophile_deletes: Union[list, callable] = None,
    ):
        """
        Set the reactivity of the FunctionalGroup.
        Any reactivity that is not provided will not be changed from the current state.

        Parameters
        ----------
        nucleophile_bonder : int
            The index of the nucleophile bonder.
        electrophile_bonder : int
            The index of the electrophile bonder.
        nucleophile_deletes : list or callable
            The indices of atoms to delete if the functional group is used as a nucleophile.
            Or a function that will return a list of atoms to delete. This function will receive the following arguments: `molecule` (Molecule), `residue` (Residue), `bonder` (Atom), `assignment` (list[Atom] that were identified as belonging to the functional group).
        electrophile_deletes : list or callable
            The indices of atoms to delete if the functional group is used as an electrophile.
            Or a function that will return a list of atoms to delete. This function will receive the following arguments: `molecule` (Molecule), `residue` (Residue), `bonder` (Atom), `assignment` (list[Atom] that were identified as belonging to the functional group).

        Returns
        -------
        FunctionalGroup
            A new FunctionalGroup with the new reactivity.
        """
        if nucleophile_bonder is not None:
            self._nucleophile_bonder = nucleophile_bonder
        if electrophile_bonder is not None:
            self._electrophile_bonder = electrophile_bonder
        if nucleophile_deletes is not None:
            self._nucleophile_deletes = nucleophile_deletes
        if electrophile_deletes is not None:
            self._electrophile_deletes = electrophile_deletes
        return self

    def set_bonders(self, nucleophile: int = None, electrophile: int = None):
        """
        Set the bonder indices for the atoms that will connect to another atom when using the functional group to make a Linkage.

        Parameters
        ----------
        nucleophile : int
            The index of the nucleophile bonder.
        electrophile : int
            The index of the electrophile bonder.
        """
        self._nucleophile_bonder = nucleophile
        self._electrophile_bonder = electrophile

    def set_deletes(
        self,
        nucleophile: Union[int, list[int]] = None,
        electrophile: Union[int, list[int]] = None,
    ):
        """
        Set the delete indices for the atoms that should be deleted when the functional group is used to infer a Linkage.
        Alternatively, a callable can be provided that will receive the following arguments: `molecule` (Molecule), `residue` (Residue), `bonder` (Atom), `assignment` (list[Atom]).
        It must return a list of Atoms to delete.

        Parameters
        ----------
        nucleophile : list or callable
            The indices of atoms to delete if the functional group is used as a nucleophile.
            Or a function that will return a list of atoms to delete.
        electrophile : list or callable
            The indices of atoms to delete if the functional group is used as an electrophile.
            Or a function that will return a list of atoms to delete.
        """
        self._nucleophile_deletes = nucleophile
        self._electrophile_deletes = electrophile

    def get_atom_assignment(self):
        """
        Get the atom assignment of the functional group.
        """
        return self._assignment

    def matches(self, molecule, atoms: list, connectivity: list) -> bool:
        """
        Check if the atoms match the functional group.

        Parameters
        ----------
        atoms : list
            The atoms to check. The first atom must be the center atom.
        connectivity : list
            The connectivity of the atoms. This list contains tuples with the indices
            of the atoms that are connected.

        Returns
        -------
        bool
            True if the atoms match the functional group, False otherwise.
        """
        if len(atoms) < self.n:
            return False
        elif atoms[0].element != self.atoms[0]:
            return False
        elif not self._match_connectivity(connectivity):
            return False
        elif not self._match_elements(atoms):
            return False

        ref_coords = np.array([atom.coord for atom in atoms])
        out_coords = self.geometry.make_coords(*ref_coords[: self.geometry.max_points])

        d = cdist(out_coords, ref_coords)
        d = d < 0.95
        if not d.sum() == len(atoms):
            return False

        assignment = self._assign_atoms(molecule, atoms)
        if len(assignment) != self.n:
            return False
        self._assignment = assignment
        return True

    def apply_connectivity(self, molecule, atoms: list):
        """
        Apply the connectivity of the functional group to a set of atoms in a molecule.

        Parameters
        ----------
        molecule : Molecule
            The molecule to apply the bonds to.
        atoms : list
            The atoms to apply the bonds to.
        """
        if self._assignment is None:
            self._assignment = self._assign_atoms(molecule, atoms)

        a = self._assignment[0]
        bonds_a = molecule._get_bonds((a,), None)
        _bonds = {}
        assigned = set()
        for cdx, c in enumerate(self.connectivity):
            b = self._assignment[cdx + 1]
            bonds_b = molecule._get_bonds((b,), None)
            _bonds[b] = bonds_b
            if infer.has_free_valence(
                a, bonds_a, needed=c - 1
            ) and infer.has_free_valence(b, bonds_b, needed=c - 1):
                molecule.get_bond(a, b).order = c
                assigned.add((a, b))

        if self.invertable and len(assigned) < len(self.connectivity):
            for cdx, c in enumerate(reversed(self.connectivity)):
                if (a, b) in assigned:
                    continue
                b = self._assignment[cdx + 1]
                bonds_b = _bonds[b]
                if infer.has_free_valence(
                    a, bonds_a, needed=c - 1
                ) and infer.has_free_valence(b, bonds_b, needed=c - 1):
                    molecule.get_bond(a, b).order = c

        self._assignment = None

    def infer_electrophile_atoms(self, molecule, residue=None):
        """
        Infer the atoms to use for a Linkage when the functional group is used as an electrophile.

        Parameters
        ----------
        molecule : Molecule
            The molecule that contains the functional group.
        residue : Residue, optional
            The residue that contains the functional group.
            By default the attached residue is used.

        Returns
        -------
        Atom
            The bonder atom.
        list
            The atoms to delete.
        """
        if self._electrophile_bonder is None:
            raise ValueError("No electrophile bonder was set.")

        if residue is not None:
            residue = molecule.get_residue(residue)
        else:
            residue = molecule.attach_residue or molecule.get_residue(-1)

        assignment = self._assign_atoms(molecule, residue.get_atoms())

        bonder = assignment[self._electrophile_bonder]
        deletes = self._infer_deletes(
            molecule, residue, assignment, bonder, self._electrophile_deletes
        )

        return bonder, deletes

    def infer_nucleophile_atoms(self, molecule, residue=None):
        """
        Infer the atoms to use for a Linkage when the functional group is used as a nucleophile.

        Parameters
        ----------
        molecule : Molecule
            The molecule that contains the functional group.
        residue : Residue, optional
            The residue that contains the functional group.
            By default the attached residue is used.

        Returns
        -------
        Atom
            The bonder atom.
        list
            The atoms to delete.
        """
        if self._nucleophile_bonder is None:
            raise ValueError("No nucleophile bonder was set.")

        if residue is not None:
            residue = molecule.get_residue(residue)
        else:
            residue = molecule.attach_residue or molecule.get_residue(-1)

        assignment = self._assign_atoms(molecule, residue.get_atoms())

        bonder = assignment[self._nucleophile_bonder]
        deletes = self._infer_deletes(
            molecule, residue, assignment, bonder, self._nucleophile_deletes
        )

        return bonder, deletes

    def _infer_deletes(self, molecule, residue, assignment, bonder, deletes):
        if deletes is None:
            return None

        _deletes = []
        if callable(deletes):
            _deletes = deletes(molecule, residue, bonder, assignment)
        else:
            _deletes = [assignment[i] for i in deletes]
        return _deletes

    def _assign_atoms(self, molecule, atoms: list):
        """
        Assign atoms from a molecule to the functional group to infer the connectivity.

        Parameters
        ----------
        molecule : Molecule
            The molecule to assign the atoms to.
        atoms : list
            The atoms to assign to the functional group.
        """
        matches = {}
        G = molecule._AtomGraph
        atoms = list(atoms)
        for cdx, c in enumerate(self.constraints):
            for a in atoms:
                if not a.element == self.atoms[cdx]:
                    continue
                if c is None:
                    matches[cdx] = a
                    atoms.remove(a)
                    break
                elif c(G, a):
                    matches[cdx] = a
                    atoms.remove(a)
                    break
        return matches

    @staticmethod
    def _element_hist(elements):
        hist = {}
        for e in elements:
            if e in hist:
                hist[e] += 1
            else:
                hist[e] = 1
        return hist

    def _match_elements(self, atoms):
        hist = self._element_hist(i.element for i in atoms)
        if not self._set_atoms.issubset(hist.keys()):
            return False
        for k, v in hist.items():
            if v < self._hist.get(k, 0):
                return False
        return True

    def _match_connectivity(self, connectivity):
        for c in self._connectivity:
            if c not in connectivity:
                return False
        return True


carbonyl = FunctionalGroup(
    "carbonyl",
    1,
    geometry.trigonal_planar,
    ("C", "O"),
    (2,),
    (
        None,
        constraints.has_n_neighbors(1),
    ),
)
carbonyl.set_bonders(None, 0)
aldehyde = carbonyl

carboxyl = FunctionalGroup(
    "carboxyl",
    2,
    geometry.trigonal_planar,
    ("C", "O", "O"),
    (2, 1),
    (
        constraints.multi_constraint(
            constraints.neighbors_all("O"),
            constraints.has_n_neighbors(3),
        ),
        constraints.has_n_neighbors(1),
        constraints.multi_constraint(
            constraints.neighbors_all("C", "H"),
            constraints.extended_neighbors_all(2, "C", "O", "H"),
        ),
    ),
)
carboxyl.set_bonders(2, 0)
carboxyl.set_deletes(
    electrophile=_H_neighbor_of_assigned_atom(2),
)

amide = FunctionalGroup(
    "amide",
    2,
    geometry.trigonal_planar,
    ("C", "O", "N"),
    (2, 1),
    (
        constraints.neighbors_all("O", "N"),
        constraints.has_n_neighbors(1),
        constraints.extended_neighbors_all(2, "C", "O"),
    ),
)
amide.set_bonders(2, 0)
amide.set_deletes(
    electrophile=_H_neighbor_of_assigned_atom(2),
)


alkene = FunctionalGroup(
    "alkene",
    1,
    geometry.trigonal_planar,
    ("C", "C"),
    (2,),
    (
        constraints.neighbors_not("O", "N", "S"),
        constraints.neighbors_not("O", "N", "S"),
    ),
)

alkyne = FunctionalGroup(
    "alkyne",
    1,
    geometry.linear,
    ("C", "C"),
    (3,),
    (
        constraints.neighbors_not("O", "N", "S"),
        constraints.neighbors_not("O", "N", "S"),
    ),
)

aromatic = FunctionalGroup(
    "aromatic",
    3,
    geometry.trigonal_planar,
    ("C", "C", "C"),
    (2, 1),
    (
        constraints.neighbors_all("C", "C"),
        constraints.neighbors_all("C", "C"),
        constraints.neighbors_all("C", "C"),
    ),
    invertable=True,
)

hydroxyl = FunctionalGroup(
    "hydroxyl",
    1,
    geometry.tetrahedral,
    ("C", "O"),
    (1,),
    (
        constraints.neighbors_all("O"),
        constraints.neighbors_all("C", "H"),
    ),
)
hydroxyl.set_bonders(1, 0)
hydroxyl.set_deletes(
    electrophile=_H_neighbor_of_assigned_atom(1),
)

amine = FunctionalGroup(
    "amine",
    1,
    geometry.tetrahedral,
    ("C", "N"),
    (1,),
    (
        constraints.neighbors_all("N"),
        constraints.neighbors_all("C", "H"),
    ),
)
amine.set_bonders(1, 0)
amine.set_deletes(
    electrophile=_H_neighbor_of_assigned_atom(1),
)

thiol = FunctionalGroup(
    "thiol",
    1,
    geometry.tetrahedral,
    ("C", "S"),
    (1,),
    (
        constraints.neighbors_all("S"),
        constraints.neighbors_all("C", "H"),
    ),
)
thiol.set_bonders(1, 0)
thiol.set_deletes(
    electrophile=_H_neighbor_of_assigned_atom(1),
)


higher_order_groups = [
    carbonyl,
    carboxyl,
    amide,
    alkene,
    alkyne,
    aromatic,
]
"""
Functional groups that contain bonds with a higher order than only single bonds.
"""

if __name__ == "__main__":

    import buildamol as bam

    bam.load_amino_acids()
    tyr = bam.get_compound("TYR")
    tyr1, tyr2 = tyr.copy(2)

    b1, d1 = carboxyl.infer_electrophile_atoms(tyr1)
    b2, d2 = amine.infer_nucleophile_atoms(tyr2)
    link = bam.linkage(
        b1.id,
        b2.id,
        [i.id for i in d1] if d1 else None,
        [i.id for i in d2] if d2 else None,
    )
    mol = tyr1 % link + tyr2
    mol.show()

    # b, d = carboxyl.infer_linkage(mol)
    # v = mol.draw()
    # v.draw_atoms(b, *d, colors="orange")
    # v.show()

    # mol = mol % "LINK" * 2
    # for b in mol.get_bonds():
    #     b.single()

    # matches = {}
    # for a in mol.get_atoms():
    #     for i in (
    #         carbonyl,
    #         carboxyl,
    #         # hydroxyl,
    #         # amine,
    #         amide,
    #         # thiol,
    #         # alkene,
    #         # alkyne,
    #         aromatic,
    #     ):
    #         neighs = mol.get_neighbors(a, 1)
    #         atoms = [a, *neighs]
    #         connetivity = [(0, i + 1) for i in range(len(neighs))]
    #         m = i.matches(mol, atoms, connetivity)
    #         # if m:
    #         #     print(
    #         #         atoms,
    #         #         i.id,
    #         #     )
    #         if m:
    #             i.apply_connectivity(mol, atoms)
    #             atoms = tuple(atoms)
    #             if atoms in matches:
    #                 if i.rank > matches[atoms].rank:
    #                     matches[atoms] = i
    #             else:
    #                 matches[atoms] = i

    # for k, v in matches.items():
    #     print(k, v.id)

    # mol.show()

    # # i = aromatic
    # # A = mol.get_atom("CG")
    # # neighs = mol.get_neighbors(A, 6)
    # # atoms = [A, *neighs]
    # # connectivity = []
    # # for a in atoms:
    # #     bonds = mol.get_bonds(a)
    # #     for b in bonds:
    # #         if b[0] in atoms and b[1] in atoms:
    # #             connectivity.append((atoms.index(b[0]), atoms.index(b[1])))
    # # m = i.matches(atoms, connectivity)
    # # if m:
    # #     print(
    # #         atoms,
    # #         i.id,
    # #     )

    # # mol.show()
