"""
The basic Class for Molecular Graphs
"""

from abc import abstractmethod
import warnings

import Bio.PDB as bio
import networkx as nx
import numpy as np
from scipy.spatial.transform import Rotation


class BaseGraph(nx.Graph):
    """
    The basic class for molecular graphs
    """

    def __init__(self, id, bonds: list):
        super().__init__(bonds)
        self.id = id
        self._structure = None
        self._neighborhood = None
        self._locked_edges = set()
        self._structure_was_searched = False

    @property
    def structure(self):
        """
        Returns the underlying `bio.PDB.Structure` object
        """
        if not self._structure_was_searched:
            self._structure = self._get_structure()
            self._structure_was_searched = True
        return self._structure

    @property
    def chains(self):
        """
        Returns the chains in the molecule
        """
        if not self.structure:
            return
        return list(self.structure.get_chains())

    @property
    def residues(self):
        """
        Returns the residues in the molecule
        """
        if not self.structure:
            return
        return list(self.structure.get_residues())

    @property
    def atoms(self):
        """
        Returns the atoms in the molecule
        """
        if not self.structure:
            return
        return list(self.structure.get_atoms())

    @property
    def nodes_in_cycles(self) -> set:
        """
        Returns the nodes in cycles
        """
        return nx.cycle_basis(self)

    @property
    def bonds(self):
        """
        Returns the bonds in the molecule
        """
        return list(self.edges)

    def show(self):
        """
        Show the graph
        """
        self.draw().show()

    def draw(self):
        """
        Prepare a 3D view of the graph but do not show it yet

        Returns
        -------
        PlotlyViewer3D
            A 3D viewer
        """
        raise NotImplementedError

    @abstractmethod
    def get_neighbors(self, node, n: int = 1, mode="upto"):
        """
        Get the neighbors of a node

        Parameters
        ----------
        node
            The target node

        n : int, optional
            The number of edges to separate the node from its neighbors.

        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` edges
            - "exact": get all neighbors exactly `n` edges away

        Returns
        -------
        set
            The neighbors of the node
        """
        raise NotImplementedError

    def get_descendants(self, node_1, node_2):
        """
        Get all descendant nodes that come after a specific edge
        defined in the direction from node1 to node2 (i.e. get all
        nodes that come after node2). This method is directed
        in contrast to the `get_neighbors()` method, which will get all neighboring
        nodes of an anchor node irrespective of direction.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge

        Returns
        -------
        set
            The descendant nodes

        Examples
        --------
        In case of this graph:
        
        .. code-block::

            A---B---C---D---E
                \\
                F---H
                |
                G
        
        ```
        A---B---C---D---E
             \\
              F---H
              |
              G
        ```
        
        >>> graph.get_descendants("B", "C")
        {"D", "E"}
        >>> graph.get_descendants("B", "F")
        {"H", "G"}
        >>> graph.get_descendants("B", "A")
        set() # because in this direction there are no other nodes
        """
        neighs = set(self.adj[node_2])
        neighs.remove(node_1)
        _seen = set((node_1, node_2))
        _new_neighs = set(neighs)
        while len(_new_neighs) > 0:
            neigh = _new_neighs.pop()
            descendants = set(self.adj[neigh])
            descendants -= _seen
            _seen.add(neigh)
            neighs.add(neigh)
            if len(descendants) == 0:
                continue
            _new_neighs.update(descendants)
        return neighs

    def get_ancestors(self, node_1, node_2):
        """
        Get all ancestor nodes that come before a specific edge
        defined in the direction from node1 to node2 (i.e. get all
        nodes that comebefore node1). This method is directed
        in contrast to the `get_neighbors()` method, which will get all neighboring
        nodes of an anchor node irrespective of direction.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge

        Returns
        -------
        set
            The ancestor nodes

        Examples
        --------
        In case of this graph:
        
        .. code-block::

            A---B---C---D---E
                \\
                F---H
                |
                G
        
        ```
        A---B---C---D---E
             \\
              F---H
              |
              G
        ```
        
        >>> graph.get_ancestors("B", "C")
        {"A", "F", "G", "H"}
        >>> graph.get_ancestors("F", "B")
        {"H", "G"}
        >>> graph.get_ancestors("A", "B")
        set() # because in this direction there are no other nodes
        """
        return self.get_descendants(node_2, node_1)

    def direct_edges(self):
        """
        Sort all edges such that the first node is always earlier
        in the sequence than the second node.
        """
        raise NotImplementedError

    def lock_edge(self, node_1, node_2):
        """
        Lock an edge, preventing it from being rotated.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge
        """
        self._locked_edges.add((node_1, node_2))

    def unlock_edge(self, node_1, node_2):
        """
        Unlock an edge, allowing it to be rotated.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge
        """
        if (node_1, node_2) in self._locked_edges:
            self._locked_edges.remove((node_1, node_2))

    def is_locked(self, node_1, node_2):
        """
        Check if an edge is locked

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge

        Returns
        -------
        bool
            Whether the edge is locked
        """
        return (node_1, node_2) in self._locked_edges

    def get_locked_edges(self):
        """
        Get all locked edges

        Returns
        -------
        set
            The locked edges
        """
        return self._locked_edges

    def get_unlocked_edges(self):
        """
        Get all unlocked edges

        Returns
        -------
        set
            The unlocked edges
        """
        return set(self.edges) - self._locked_edges

    def lock_all(self):
        """
        Lock all edges
        """
        self._locked_edges = set(self.edges)

    def unlock_all(self):
        """
        Unlock all edges
        """
        self._locked_edges = set()

    def rotate_around_edge(
        self,
        node_1,
        node_2,
        angle: float,
        descendants_only: bool = False,
        update_coords: bool = True,
    ):
        """
        Rotate descending nodes around a specific edge by a given angle.

        Parameters
        ----------
        node_1, node_2
            The nodes that define the edge around which to rotate.
        angle: float
            The angle to rotate by, in radians.
        descendants_only: bool, optional
            Whether to only rotate the descending nodes, by default False, in which case the entire graph
            will be rotated.
        update_coords: bool, optional
            Whether to update the coordinates of the nodes after rotation, by default True.

        Returns
        -------
        new_coords: dict
            The new coordinates of the nodes after rotation.
        """
        # ---------- sanity checks ----------
        # We can skip these here for a little performance boost
        # since we should assume that these methods are only ever
        # called from their wrappers in the entity classes...
        # ---------- sanity checks ----------
        # if node_1 not in self.nodes or node_2 not in self.nodes:
        #     raise ValueError("One or more nodes not in graph!")
        if node_1 is node_2:
            raise ValueError("Cannot rotate around an edge with only one node!")
        elif self.is_locked(node_1, node_2):
            raise ValueError("Cannot rotate around a locked edge!")

        # we need to get a reference node index to normalise the rotated
        # coordinates to the original coordinate system
        # indices = list(self.nodes)
        idx_1 = next(idx for idx, i in enumerate(self.nodes) if i is node_1)

        # define the axis of rotation as the cross product of the edge's vectors
        edge_vector = node_2.coord - node_1.coord
        edge_vector /= np.linalg.norm(edge_vector)

        # create the rotation matrix
        r = Rotation.from_rotvec(angle * edge_vector)

        # create a numpy array of the node coordinates
        if descendants_only:
            nodes = {i: i.coord for i in self.get_descendants(node_1, node_2)}
            nodes[node_2] = node_2.coord
        else:
            nodes = {i: i.coord for i in self.nodes}

        node_coords = np.array(tuple(nodes.values()))

        # indices = list(nodes.keys())
        # idx_2 = indices.index(node_2)
        idx_2 = next(idx for idx, i in enumerate(nodes.keys()) if i is node_2)

        # apply the rotation matrix to the node coordinates
        node_coords_rotated = r.apply(node_coords)

        # now adjust for the translatisonal shift around the axis
        _diff = node_coords_rotated[idx_2] - node_coords[idx_2]
        node_coords_rotated -= _diff

        # update the node coordinates in the graph
        new_coords = {i: node_coords_rotated[idx] for idx, i in enumerate(nodes.keys())}

        if update_coords:
            for node, coord in new_coords.items():
                node.coord = coord

        _new_coords = {i: i.coord for i in self.nodes}
        _new_coords.update(new_coords)
        return _new_coords

    def _get_structure(self):
        """
        Get the underlying `bio.PDB.Structure` object
        """
        if not hasattr(list(self.nodes)[0], "get_parent"):
            warnings.warn("Nodes are not Biopython entities with linked parents!")
            return None
        structure = list(self.nodes)[0].get_parent()
        if structure is None:
            warnings.warn("Nodes do not seem to have linked parents!")
            return None
        while not isinstance(structure, bio.Structure.Structure):
            structure = structure.get_parent()
        return structure

    def __str__(self):
        lines = "\n".join(nx.generate_network_text(self))
        return lines


if __name__ == "__main__":
    import biobuild as bb
    from timeit import timeit

    import seaborn as sns
    import matplotlib.pyplot as plt

    mol = bb.molecule(
        "/Users/noahhk/GIT/biobuild/biobuild/optimizers/_testing/files/EX7.json"
    )
    v = mol.draw()
    g = BaseGraph(None, mol.bonds)
    ref_atom_graph = mol.make_atom_graph()

    a, b = mol.get_atoms(68, 65)

    x = g.get_descendants(a, b)
    for i in x:
        v.draw_atom(i, color="purple")

    measure_performance = True
    repeats = 500
    number = 800
    if measure_performance:
        test_old = lambda: BaseGraph.get_descendants_old(ref_atom_graph, a, b)
        test_new = lambda: g.get_descendants(a, b)
        # test_2 = lambda: g.get_descendants_2(a, b)

        times_old = [timeit(test_old, number=number) for _ in range(repeats)]
        times_new = [timeit(test_new, number=number) for _ in range(repeats)]
        # times_2 = [timeit(test_2, number=number) for _ in range(repeats)]

        sns.distplot(times_old, label="old", kde=True, bins=20)
        sns.distplot(times_new, label="new", kde=True, bins=20)
        # sns.distplot(times_2, label="2", kde=True, bins=20)
        plt.legend()
        plt.show()

    pass

# v.draw_edges((a, b), color="limegreen", linewidth=4, elongate=1.1)
# ref_decendants = ref_atom_graph.get_descendants(a, b)

# descendants = g.get_descendants_2(a, b)

# for i in descendants:
#     v.draw_atom(i, color="purple")
# for i in ref_decendants:
#     v.draw_atom(i, color="orange")
# v.show()
# pass
