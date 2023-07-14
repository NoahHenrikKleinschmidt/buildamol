"""
The basic Class for Molecular Graphs
"""

from abc import abstractmethod
import warnings

import Bio.PDB as bio
import networkx as nx
import numpy as np
from scipy.spatial.transform import Rotation
import biobuild.utils.visual as vis


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
        MoleculeViewer3D
            The 3D viewer
        """
        return vis.MoleculeViewer3D(self)

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
            The descendants of the node

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
        # if node_1 == node_2:
        if node_1 is node_2:
            raise ValueError(
                "Cannot get descendants if only one node is given (no direction)!"
            )

        neighbors = self.get_neighbors(node_2)
        neighbors.remove(node_1)
        if len(neighbors) == 0:
            return neighbors

        _new_neighbors = neighbors.copy()
        _seen = set((node_1, node_2))

        while len(_new_neighbors) > 0:
            neighbor = _new_neighbors.pop()
            descendants = self.get_neighbors(neighbor)

            descendants -= _seen
            _seen.add(neighbor)
            neighbors.add(neighbor)

            if len(descendants) == 0:
                continue
            _new_neighbors.update(descendants)

        return neighbors

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

        # get the node coordinates as a dictionary
        # # we do this in order to update the node attributes later...
        # node_dict = nx.get_node_attributes(self, 'coord')

        if node_1 not in self.nodes or node_2 not in self.nodes:
            raise ValueError("One or more nodes not in graph!")
        elif node_1 == node_2:
            raise ValueError("Cannot rotate around an edge with only one node!")
        elif self.is_locked(node_1, node_2):
            raise ValueError("Cannot rotate around a locked edge!")

        # we need to get a reference node index to normalise the rotated
        # coordinates to the original coordinate system
        indices = list(self.nodes)
        idx_1 = indices.index(node_1)

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

        indices = list(nodes.keys())
        idx_2 = indices.index(node_2)

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

            # set the node attributes in the graph
            nx.set_node_attributes(self, new_coords, "coord")

        return new_coords

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
