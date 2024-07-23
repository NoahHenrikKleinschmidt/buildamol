from typing import Union
import networkx as nx
import Bio.PDB as bio
import numpy as np

import buildamol.structural as struct
from buildamol.graphs.base_graph import BaseGraph
import buildamol.utils.visual as vis


class ResidueGraph(BaseGraph):
    """
    A graph representation of residues bonded together as an abstraction of a large contiguous molecule.
    """

    __idx_method__ = lambda x: (
        x.serial_number if hasattr(x, "serial_number") else x.id[1]
    )

    def __init__(self, id, bonds: list):
        super().__init__(id, bonds)
        self._AtomGraph = None
        self._molecule = None
        self._atomic_bonds = {}
        self._atomic_bonds_list = []

        self._residues = {i.id: i for i in self.nodes}

    @classmethod
    def from_molecule(cls, mol, detailed: bool = False, locked: bool = True):
        """
        Create a ResidueGraph from a molecule object.

        Parameters
        ----------
        mol : Molecule
            The molecule object
        detailed : bool
            Whether to make a "detailed" residue graph representation
            including the atomic-scale bonds between residues. If True,
            locked bonds can be directly migrated from the molecule.
        locked : bool
            Whether to migrate locked bonds from the molecule. This
            is only possible if detailed is True.

        Returns
        -------
        ResidueGraph
            The ResidueGraph representation of the molecule
        """
        if len(mol.residues) < 2:
            new = cls(mol.id, [])
            new._AtomGraph = mol._AtomGraph
            new._structure = mol.structure
            new._molecule = mol
            res = mol.residues[0]
            new.add_node(res)

            if detailed:
                for node in res.child_list:
                    new.add_edge(node, res)

            return new

        try:
            connections = mol.get_residue_connections(triplet=True)
        except:
            mol.infer_residue_connections()
            connections = mol.get_residue_connections(triplet=True)
        main_connections = [
            (i.get_parent(), j.get_parent())
            for i, j in connections
            if i.get_parent() != j.get_parent()
        ]
        new = cls(mol.id, main_connections)
        new._AtomGraph = mol._AtomGraph
        new._structure = mol.structure
        new._molecule = mol

        new._atomic_bonds_list = list(connections)

        for bond in connections:
            parent_1 = bond[0].get_parent()
            parent_2 = bond[1].get_parent()
            new._atomic_bonds.setdefault((parent_1, parent_2), []).append(bond)

        if detailed:
            new.make_detailed()
            if locked:
                new._locked_edges.update(
                    (i for i in mol._AtomGraph._locked_edges if tuple(i) in new.edges)
                )
        return new

    @classmethod
    def from_AtomGraph(cls, atom_graph, infer_connections: bool = None):
        """
        Create a ResidueGraph from an AtomGraph.

        Parameters
        ----------
        atom_graph : AtomGraph
            The AtomGraph representation of the molecule
        infer_connections : bool
            Whether to infer the bonds between residues from the atom-level bonds.
            If the AtomGraph already contains atom-level bonds that connect different residues,
            this is not necessary. If this is set to None, connections will be inferred automatically
            if no atom-level bonds are present in the AtomGraph.

        Returns
        -------
        ResidueGraph
            The ResidueGraph representation of the molecule
        """
        id = atom_graph.id

        connections = [
            i for i in atom_graph.edges if i[0].get_parent() != i[1].get_parent()
        ]
        main_connections = [
            (p1.get_parent(), p2.get_parent()) for p1, p2 in connections
        ]
        if infer_connections is None:
            infer_connections = len(connections) == 0

        if infer_connections:
            main_connections = struct.infer_residue_connections(
                atom_graph.structure, triplet=False
            )
            main_connections = [
                (p1.get_parent(), p2.get_parent()) for p1, p2 in main_connections
            ]
            connections = struct.infer_residue_connections(
                atom_graph.structure, triplet=True
            )

        if len(main_connections) < 1:
            raise ValueError("No connections between residues could be inferred!")

        new = cls(id, main_connections)
        new._AtomGraph = atom_graph

        for bond in connections:
            parent_1 = bond[0].get_parent()
            parent_2 = bond[1].get_parent()
            new._atomic_bonds.setdefault((parent_1, parent_2), []).append(bond)

        return new

    def get_residue(self, r):
        """
        Get a residue in the molecule.

        Parameters
        ----------
        r : str or Residue
            The residue or it's id

        Returns
        -------
        Residue
            The residue
        """
        return self._molecule.get_residue(r)

    def add_atomic_bonds(self, *edges):
        """
        Add atom-level bonds to the graph.

        Parameters
        ----------
        *edges
            The edges to add
        """
        for edge in edges:
            p1, p2 = edge[0].get_parent(), edge[1].get_parent()
            self.remove_edges_from(((p1, p2),))
            if p1 is not p2:
                new_edges = ((p1, edge[0]), (p2, edge[1]), edge)
            else:
                new_edges = ((p1, edge[0]), edge)

            self.add_edges_from(new_edges)
            # self._atomic_bonds.setdefault(
            #     (edge[0].get_parent(), edge[1].get_parent()), []
            # ).append(edge)

    def make_detailed(
        self,
        include_samples: bool = True,
        include_far_away: bool = False,
        include_heteroatoms: bool = False,
        include_clashes: bool = True,
        n_samples: Union[int, float] = 0.5,
        f: float = 1.0,
        no_hydrogens: bool = False,
    ) -> "ResidueGraph":
        """
        Use a detailed representation of the residues in the molecule by adding the specific atoms
        that connect the residues together. This is useful for visualization and analysis.

        Note
        ----
        This function is not reversible. It is applied in-place.

        Parameters
        ----------
        include_samples : bool
            If True, a number of atoms are sampled from each residue and included in the detailed
            representation.

        include_far_away : bool
            If True, atoms that are not involved in residue connections are also included if their
            distance to the residue's center of mass is greater than f * the 75th percentile of
            atom distances to the residue's center of mass.

        include_heteroatoms : bool
            If True, all hetero-atoms are included in the detailed representation, regardless of
            their distance to the residue center of mass.

        include_clashes : bool
            If True, all atoms that are involved in a clash are included in the detailed representation.

        n_samples : int or float
            The number or fraction of atoms to sample from each residue if include_samples is True.
            If a fraction in range (0,1) is given instead of an integer, the number of atoms to
            sample is adjusted according to the residue size.

        f : float
            The factor by which the 75th percentile of atom distances to the residue's center of mass
            is multiplied to determine the cutoff distance for outlier atoms. This is only used if
            include_outliers is True.

        no_hydrogens : bool
            If True, hydrogens are not included in the detailed representation.
        """

        # self.clear_edges()

        _added_nodes = set()
        for edge in self._atomic_bonds_list:
            if not edge[0] in self.nodes:
                self.add_edge(edge[0], edge[0].get_parent())
            if not edge[1] in self.nodes:
                self.add_edge(edge[1], edge[1].get_parent())
            if edge[1].parent is edge[0].parent and self.has_edge(
                edge[1], edge[1].parent
            ):
                self.remove_edge(edge[1], edge[1].parent)
            if self.has_edge(edge[0].parent, edge[1].parent):
                self.remove_edge(edge[0].parent, edge[1].parent)
            # if edge[0].parent == edge[1].parent and self.has_edge(
            #     edge[1], edge[1].parent
            # ):
            #     self.remove_edge(edge[1], edge[1].parent)
            self.add_edge(*edge)
            _added_nodes.update(edge)

        triplets = struct.generate_triplets(self._atomic_bonds_list)
        for triplet in triplets:
            if triplet[0].get_parent() == triplet[2].get_parent():
                continue
            e1 = (triplet[0], triplet[0].get_parent())
            e3 = (triplet[2], triplet[2].get_parent())
            self.add_edge(*e1)
            self.add_edge(*e3)
            _added_nodes.update(e1)

        if include_samples:
            for residue in self.residues:
                if no_hydrogens:
                    atoms = np.array(
                        [i for i in residue.child_list if i.element != "H"]
                    )
                else:
                    atoms = np.array(residue.child_list)

                if n_samples < 1:
                    n = max(1, int(np.floor(len(atoms) * n_samples)))
                else:
                    n = n_samples

                samples = struct.sample_atoms_around_reference(
                    residue.center_of_mass(), atoms, num_samples=n
                )
                for i in samples:
                    if i not in _added_nodes:
                        self.add_edge(i, residue)
                        _added_nodes.add(i)

                # WORKS! BUT THE OTHER ONE IS SOO MUCH NICER!!!
                # if len(atoms) > n_samples:
                #     coords = np.array([i.coord for i in atoms])
                #     dists = coords - residue.center_of_mass()
                #     dists /= np.linalg.norm(dists, axis=1)[:, None]
                #     # evenly sample atoms that are all around the residue
                #     # by taking the dot product of the distance vectors
                #     # and the center of mass vector
                #     # and taking the n_samples atoms with the highest dot product
                #     # this is a proxy for the atoms that are "furthest" from the residue
                #     # center of mass
                #     dot = np.dot(dists, residue.center_of_mass())
                #     idx = np.argsort(dot)[-n_samples:]
                #     atoms = atoms[idx]
                # for atom in atoms:
                #     if atom not in _added_nodes:
                #         self.add_edge(atom, residue)
                #         _added_nodes.add(atom)

        if include_far_away:
            for residue in self.residues:
                outliers = struct.compute_outlier_atoms(residue, f=f)
                if no_hydrogens:
                    outliers = (i for i in outliers if i.element != "H")
                for outlier in outliers:
                    if outlier not in _added_nodes:
                        self.add_edge(outlier, residue)
                        _added_nodes.add(outlier)

        if include_heteroatoms:
            for residue in self.residues:
                for atom in residue.get_atoms():
                    if atom.element not in ("C", "H"):
                        if atom not in _added_nodes:
                            self.add_edge(atom, residue)
                            _added_nodes.add(atom)

        if include_clashes:
            for atoms in self._molecule.find_clashes():
                for atom in atoms:
                    if atom not in _added_nodes:
                        self.add_edge(atom, atom.get_parent())
                        _added_nodes.add(atom)

        # prune edge triplets of atoms that are part of the same residues
        self.prune_triplets()

        return self

    def prune_triplets(self):
        """
        Prune bond triangles where two nodes from the
        same residue are connected to each other and the residue...
        """
        for triplet in nx.cycle_basis(self):
            if len(triplet) == 3:
                length_12 = np.linalg.norm(triplet[0].coord - triplet[1].coord)
                length_23 = np.linalg.norm(triplet[1].coord - triplet[2].coord)
                length_13 = np.linalg.norm(triplet[0].coord - triplet[2].coord)
                lengths = [length_12, length_23, length_13]
                # remove the longest edge
                if lengths.index(max(lengths)) == 0:
                    self.remove_edge(triplet[0], triplet[1])
                elif lengths.index(max(lengths)) == 1:
                    self.remove_edge(triplet[1], triplet[2])
                else:
                    self.remove_edge(triplet[0], triplet[2])

    def draw(self):
        v = vis.ResidueGraphViewer3D()
        v.link(self)
        return v

    def lock_centers(self):
        """
        Lock any edges that connect residue centers of mass to their constituent atoms.
        This only applies to detailed graphs.
        """
        for edge in self.edges:
            if isinstance(edge[0], bio.Residue.Residue) and isinstance(
                edge[1], bio.Atom.Atom
            ):
                self._locked_edges.add(edge)
            elif isinstance(edge[0], bio.Atom.Atom) and isinstance(
                edge[1], bio.Residue.Residue
            ):
                self._locked_edges.add(edge)

    @property
    def residues(self):
        """
        Get the residues in the molecule.

        Returns
        -------
        list
            The residues in the molecule
        """
        return list(sorted(self._residues.values()))

    @property
    def atomic_bonds(self):
        """
        Get the atomic-level bonds in the molecule.

        Returns
        -------
        dict
            The atomic-level bonds in the molecule
        """
        return self._atomic_bonds

    def to_AtomGraph(self):
        """
        Convert the ResidueGraph to an AtomGraph.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        return self._AtomGraph

    def get_atomic_bond(self, residue1, residue2) -> tuple:
        """
        Get the atomic-level bond between two residues.

        Parameters
        ----------
        residue1 : Residue or str
            The first residue or it's id
        residue2 : Residue or str
            The second residue or it's id

        Returns
        -------
        tuple
            The atomic bond between the two residues
        """
        if isinstance(residue1, str):
            residue1 = self._residues.get(residue1)
        if isinstance(residue2, str):
            residue2 = self._residues.get(residue2)

        bond = self._atomic_bonds.get((residue1, residue2))
        if bond is None:
            bond = self._atomic_bonds.get((residue2, residue1))
        return bond

    def get_neighbors(self, residue: bio.Residue.Residue, n: int = 1, mode="upto"):
        """
        Get the neighbors of a residue

        Parameters
        ----------
        residue : bio.Residue.Residue
            The target residue
        n : int, optional
            The number of connections to separate the residue from its neighbors.
        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` bonds
            - "exact": get all neighbors exactly `n` bonds away

        Returns
        -------
        set
            The neighbors of the residue
        """
        if not self._neighborhood:
            self._neighborhood = struct.ResidueNeighborhood(self)
        return self._neighborhood.get_neighbors(residue, n, mode)

    def search_by_constraints(self, constraints: list) -> list:
        if not self._neighborhood:
            self._neighborhood = struct.ResidueNeighborhood(self)
        return self._neighborhood.search_by_constraints(constraints)

    def centers_of_mass(self):
        """
        Get the centers of mass of the residues in the molecule.

        Returns
        -------
        dict
            The centers of mass of the residues in the molecule
        """
        return {residue.id: residue.center_of_mass() for residue in self.residues}

    def find_rotatable_edges(
        self,
        root_node=None,
        min_descendants: int = 1,
        min_ancestors: int = 1,
        max_descendants: int = None,
        max_ancestors: int = None,
    ):
        edges = super().find_rotatable_edges(
            root_node, min_descendants, min_ancestors, max_descendants, max_ancestors
        )
        edges = [
            i for i in edges if i[0] not in self.residues and i[1] not in self.residues
        ]
        return edges


if __name__ == "__main__":
    import buildamol as bam

    f = "support/examples/man9.pdb"
    mol = bam.Molecule.from_pdb(f)
    mol.infer_bonds(restrict_residues=False)
    b = mol.get_bonds("C1", "O4")

    man = ResidueGraph.from_molecule(mol)
    man.add_atomic_bonds(*b)
    man.make_detailed(include_far_away=True, n_samples=0.2)
    x = man.find_rotatable_edges()
    v = man.draw()
    v.draw_edges(*mol.bonds, color="blue", opacity=0.1)
    v.show()

    # _man = "support/examples/MAN9.pdb"
    # _man = bam.Molecule.from_pdb(_man)
    # _man.infer_bonds(restrict_residues=False)
    # man = ResidueGraph.from_molecule(_man)

    import matplotlib.pyplot as plt
    import networkx as nx

    man.make_detailed(True)

    import buildamol.utils.visual as vis

    v = vis.MoleculeViewer3D(man)
    v.show()

    print(len(list(man.bonds)))
    nx.draw(man, with_labels=True, font_weight="bold", pos=nx.spectral_layout(man))
    plt.show()
