"""
The basic Class for Molecular Graphs
"""

from abc import abstractmethod
import warnings

import Bio.PDB as bio
import networkx as nx
import numpy as np

# from scipy.spatial.transform import Rotation
import buildamol.structural.base as base


class BaseGraph(nx.Graph):
    """
    The basic class for molecular graphs
    """

    def __init__(self, id, bonds: list):
        super().__init__(bonds)
        self.id = id
        self._structure = None
        self._molecule = None
        self._neighborhood = None
        self._locked_edges = set()
        self._structure_was_searched = False
        self.__descendent_cache = {}
        self.__last_cache_size = len(self.nodes)

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
    def central_node(self):
        """
        Returns the central most node of the graph.
        This is computed based on the mean of all node coordinates.
        """
        # get the central node
        center = np.mean([i.coord for i in self.nodes])
        # get the node closest to the center
        root_node = min(self.nodes, key=lambda x: np.linalg.norm(x.coord - center))
        return root_node

    @property
    def nodes_in_cycles(self) -> set:
        """
        Returns the nodes in cycles
        """
        cycles = nx.cycle_basis(self)
        if len(cycles) == 0:
            return set()
        return set.union(*[set(i) for i in cycles])

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

    def get_descendants(self, node_1, node_2, use_cache: bool = True):
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
        use_cache: bool, optional
            Whether to use the cache for the descendants, by default True.
            If True and the graph has not received new nodes since the last time 
            the cache was updated, a simple lookup is performed. Otherwise the descendant
            nodes are recursively calculated again. 

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
        if node_1 is node_2:
            raise KeyError("Cannot get descendants of a node with itself!")

        if use_cache:
            _seen = self.__descendent_cache.get((node_1, node_2))
            if _seen:
                size, _seen = _seen
                if size == len(self.nodes):
                    return _seen

        __all_nodes = set(self.nodes)

        # if use_cache:
        #     _seen = self.__descendent_cache.get((node_2, node_1))
        #     if _seen:
        #         size, _seen = _seen
        #         if size == len(self.nodes):
        #             return __all_nodes - _seen - {node_1, node_2}

        _seen = set((node_1, node_2))
        _new_neighs = {node_2}
        descendants = set()
        while _new_neighs:
            neigh = _new_neighs.pop()
            _seen.add(neigh)

            descendants.clear()
            for d in self.adj[neigh]:
                if d in _seen:
                    continue
                _desc_from_cache = self.__descendent_cache.get((neigh, d))
                if _desc_from_cache:
                    _desc_from_cache = _desc_from_cache[1]
                    _seen.add(d)
                    _seen.update(_desc_from_cache)
                else:
                    descendants.add(d)

            _new_neighs.update(descendants)
            _new_neighs.difference_update(_seen)

        _seen.difference_update((node_1, node_2))

        self.__descendent_cache[(node_1, node_2)] = len(self.nodes), _seen
        self.__descendent_cache[(node_2, node_1)] = len(self.nodes), (
            __all_nodes - _seen - {node_1, node_2}
        )

        # v = self._molecule.draw()
        # v.draw_vector("edge", node_1.coord, node_2.coord, elongate=1.3, linewidth=4, color="limegreen")

        # for i in _seen:
        #     v.draw_point("n", i.coord, color="grey", showlegend=False)
        # v.show()

        return _seen

    def get_ancestors(self, node_1, node_2, use_cache: bool = True):
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
        use_cache: bool, optional
            Whether to use the cache for the ancestors, by default True.
            If True and the graph has not received new nodes since the last time
            the cache was updated, a simple lookup is performed. Otherwise the ancestor
            nodes are recursively calculated again.

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
        return self.get_descendants(node_2, node_1, use_cache)

    def search_by_constraints(self, constraints: list) -> list:
        """
        Search for neighboring nodes that match a set of constraints.

        Parameters
        ----------
        constraints : list
            A list of constraint functions, where each entry represents the constraints for a specific node.
            All constraints must be satisfied for all nodes in the neighborhood to be considered a match.

        Returns
        -------
        list
            A list of dictionaries where each dictionary contains nodes that match the constraints. The keys represent
            the constraint index which the nodes satisfy and the values are the nodes themselves.
        """
        raise NotImplementedError

    def find_cycles(self) -> list:
        """
        Find all cycles in the graph

        Returns
        -------
        list
            A list of cycles in the graph, where each cycle is a list of nodes
        """
        return nx.cycle_basis(self)

    def find_nodes_in_cycles(self) -> set:
        """
        Find all nodes that are in cycles

        Returns
        -------
        set
            The nodes in cycles
        """
        cycles = [set(i) for i in nx.cycle_basis(self)]
        if len(cycles) == 0:
            return set()
        return set.union(*cycles)

    def find_edges_in_cycles(self) -> set:
        """
        Find all edges that connect nodes in cycles, where both nodes are in the same cycle

        Returns
        -------
        set
            The edges in cycles
        """
        nodes_in_cycles = self.find_cycles()
        return set(
            (i, j)
            for i, j in self.edges
            if self.in_same_cycle(i, j, cycles=nodes_in_cycles)
        )

    def find_rotatable_edges(
        self,
        root_node=None,
        min_descendants: int = 1,
        min_ancestors: int = 1,
        max_descendants: int = None,
        max_ancestors: int = None,
    ) -> list:
        """
        Find all edges in the graph that are rotatable (i.e. not locked, single, and not in a circular constellation).
        You can also filter and direct the edges.

        Parameters
        ----------
        root_node
            A root node by which to direct the edges (closer to further).
        min_descendants: int, optional
            The minimum number of descendants that an edge must have to be considered rotatable.
        min_ancestors: int, optional
            The minimum number of ancestors that an edge must have to be considered rotatable.
        max_descendants: int, optional
            The maximum number of descendants that an edge must have to be considered rotatable.
        max_ancestors: int, optional
            The maximum number of ancestors that an edge must have to be considered rotatable.

        Returns
        -------
        list
            A list of rotatable edges
        """
        if not max_descendants:
            max_descendants = np.inf
        if not max_ancestors:
            max_ancestors = np.inf
        circulars = [set(i) for i in nx.cycle_basis(self)]
        # we changed stuff to generators to gain some performance
        # revert if it causes issues. We know that the root_node
        # step needs a list so we unpack if needed...
        rotatable_edges = (
            i
            for i in self.edges
            if not self.is_locked(*i)
            # and (hasattr(i[0], "element") and hasattr(i[1], "element"))
            and self[i[0]][i[1]].get("bond_order", 1) == 1
            and not self.in_same_cycle(*i, circulars)
        )
        if root_node is not None:
            rotatable_edges = list(rotatable_edges)
            _directed = nx.dfs_tree(self, root_node)
            rotatable_edges = [
                i
                for i in _directed.edges
                if i in rotatable_edges or i[::-1] in rotatable_edges
            ]

        rotatable_edges = [
            i
            for i in rotatable_edges
            if min_descendants < len(self.get_descendants(*i)) < max_descendants
            and min_ancestors < len(self.get_ancestors(*i)) < max_ancestors
        ]

        return rotatable_edges

    def find_edges(
        self,
        root_node=None,
        min_descendants: int = 1,
        min_ancestors: int = 1,
        max_descendants: int = None,
        max_ancestors: int = None,
        bond_order: int = None,
        exclude_cycles: bool = False,
        only_cycles: bool = False,
        exclude_locked: bool = False,
        only_locked: bool = False,
    ) -> list:
        """
        Find edges in the graph according to the given criteria.
        This does not restrict for edges that are rotatable.

        Parameters
        ----------
        root_node
            A root node by which to direct the edges (closer to further).
        min_descendants: int, optional
            The minimum number of descendants that an edge must have to be considered rotatable.
        min_ancestors: int, optional
            The minimum number of ancestors that an edge must have to be considered rotatable.
        max_descendants: int, optional
            The maximum number of descendants that an edge must have to be considered rotatable.
        max_ancestors: int, optional
            The maximum number of ancestors that an edge must have to be considered rotatable.
        bond_order: int or tuple, optional
            The bond order to filter by. If a tuple is given, the bond order must be one of the values in the tuple.
        exclude_cycles: bool, optional
            Whether to exclude edges that are in cycles, by default False
        only_cycles: bool, optional
            Whether to only include edges that are in cycles, by default False
        exclude_locked: bool, optional
            Whether to exclude locked edges, by default False
        only_locked: bool, optional
            Whether to only include locked edges, by default False

        Returns
        -------
        list
            A list of rotatable edges
        """
        if not max_descendants:
            max_descendants = np.inf
        if not max_ancestors:
            max_ancestors = np.inf

        matching_edges = iter(self.edges)

        if only_locked and exclude_locked:
            raise ValueError(
                "Cannot exclude and include locked edges at the same time!"
            )
        if exclude_locked:
            matching_edges = (i for i in matching_edges if i not in self._locked_edges)
        elif only_locked:
            matching_edges = (i for i in matching_edges if i in self._locked_edges)

        if bond_order is not None:
            if isinstance(bond_order, int):
                matching_edges = (
                    i
                    for i in matching_edges
                    if self[i[0]][i[1]].get("bond_order", 1) == bond_order
                )
            elif isinstance(bond_order, tuple):
                matching_edges = (
                    i
                    for i in matching_edges
                    if self[i[0]][i[1]].get("bond_order", 1) in bond_order
                )
            else:
                raise ValueError(f"Invalid datatype {type(bond_order)} for bond_order!")

        if exclude_cycles and only_cycles:
            raise ValueError("Cannot exclude and include cycles at the same time!")
        elif exclude_cycles:
            circulars = self.find_edges_in_cycles()
            if len(circulars) > 0:
                matching_edges = (
                    i for i in matching_edges if not self.in_same_cycle(*i, circulars)
                )
        elif only_cycles:
            circulars = self.find_edges_in_cycles()
            if len(circulars) > 0:
                matching_edges = (
                    i for i in matching_edges if self.in_same_cycle(*i, circulars)
                )

        if root_node is not None:
            matching_edges = list(matching_edges)
            _directed = nx.dfs_tree(self, root_node)
            matching_edges = [
                i
                for i in _directed.edges
                if i in matching_edges or i[::-1] in matching_edges
            ]

        matching_edges = [
            i
            for i in matching_edges
            if min_descendants < len(self.get_descendants(*i)) < max_descendants
            and min_ancestors < len(self.get_ancestors(*i)) < max_ancestors
        ]

        return matching_edges

    def sample_edges(
        self,
        edges: list = None,
        n: int = 3,
        m: int = 3,
    ) -> list:
        """
        Sample a number of rotatable edges from the graph. This is done
        by clustering the nodes together to sample "representive" edges
        from each cluster. This is useful for subsampling the rotatable
        edges for an optimization to reduce the search space.

        Parameters
        ----------
        edges : list, optional
            The edges to sample from, by default None, in which case all rotatable edges are sampled.
        n: int
            The number of clusters to sample from.
        m : int
            The number of edges to sample from each cluster
        root_node
            A root node to direct the edges (optional)

        Returns
        -------
        list
            A list of sampled edges
        """
        center = np.mean([i.coord for i in self.nodes])
        if edges is None:
            edges = self.find_rotatable_edges()
        rotatable_edges = np.array(edges)

        if len(rotatable_edges) == 0:
            raise ValueError("No rotatable edges found!")
        elif len(rotatable_edges) < n * m:
            return rotatable_edges.tolist()

        # x, y, z, n_neighbors_a_3, n_neighbors_b_3, n_descendants, dist_to_center
        data = np.zeros((len(rotatable_edges), 7))
        for idx, edge in enumerate(rotatable_edges):
            node_a, node_b = edge
            data[idx, 0:3] = (node_a.coord + node_b.coord) / 2
            data[idx, 3] = len(self.get_neighbors(node_a, 3))
            data[idx, 4] = len(self.get_neighbors(node_b, 3))
            data[idx, 5] = len(self.get_descendants(*edge))
            data[idx, 6] = np.linalg.norm(data[idx, 0:3] - center)

        from sklearn.cluster import KMeans

        kmeans = KMeans(n_clusters=min(n, len(rotatable_edges)), n_init="auto")
        kmeans.fit(data)
        labels = kmeans.predict(data)
        _rotatable_edges = []
        for i in range(kmeans.n_clusters):
            mask = np.where(labels == i)
            cluster = rotatable_edges[mask]
            if len(cluster) > m:
                prob = (
                    0.3 * (data[mask, 3] + data[mask, 4])
                    + data[mask, 5]
                    - data[mask, 6]
                ).squeeze()
                prob += np.abs(prob.min())
                prob /= prob.sum()

                cluster = np.random.choice(
                    np.arange(len(cluster)), m, replace=False, p=prob
                )
                cluster = rotatable_edges[mask][cluster]

            _rotatable_edges.extend(cluster.tolist())

        _rotatable_edges = [tuple(i) for i in _rotatable_edges]
        return _rotatable_edges

    def in_same_cycle(self, node_1, node_2, cycles=None) -> bool:
        """
        Check if two nodes are in the same cycle

        Parameters
        ----------
        node_1, node_2
            The nodes to check

        Returns
        -------
        bool
            True if the nodes are in the same cycle, False otherwise
        """
        if not cycles:
            cycles = nx.cycle_basis(self)
        for cycle in cycles:
            if node_1 in cycle and node_2 in cycle:
                return True
        return False

    def in_cycle(self, node, cycles=None) -> bool:
        """
        Check if a node is in a cycle

        Parameters
        ----------
        node
            The node to check

        Returns
        -------
        bool
            True if the node is in a cycle, False otherwise
        """
        if not cycles:
            cycles = nx.cycle_basis(self)
        for cycle in cycles:
            if node in cycle:
                return True
        return False

    def get_cycle(self, node, cycles=None) -> set:
        """
        Get the cycle that a node is in

        Parameters
        ----------
        node
            The node to check

        Returns
        -------
        set
            The nodes in the cycle that the node is in.
            If the node is not in a cycle, None is returned.
        """
        if not cycles:
            cycles = nx.cycle_basis(self)
        for cycle in cycles:
            if node in cycle:
                return set(cycle)
        return None

    def direct_edges(self, root_node=None, edges: list = None) -> list:
        """
        Sort the edges such that the first node in each edge
        is the one closer to the root node. If no root node is provided,
        the central node is used.

        Parameters
        ----------
        root_node
            The root node to use for sorting the edges. If not provided, the central node is used.
        edges : list, optional
            The edges to sort, by default None, in which case
            all edges are sorted.

        Returns
        -------
        list
            The sorted edges
        """
        if not root_node:
            root_node = self.central_node

        if edges is None:
            edges = list(self.edges)

        if root_node not in self.nodes:
            raise ValueError(f"Root node {root_node} not in graph")

        _directed = nx.dfs_tree(self, source=root_node).edges

        _tupled_edges = [(i[0], i[1]) for i in edges]
        out = [edge if edge in _directed else edge[::-1] for edge in _tupled_edges]
        return out

    def clear_cache(self):
        """
        Clear the descendant cache
        """
        self.__descendent_cache.clear()

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
        # idx_1 = next(idx for idx, i in enumerate(self.nodes) if i is node_1)

        # define the axis of rotation as the cross product of the edge's vectors
        edge_vector = node_2.coord - node_1.coord
        edge_vector /= np.linalg.norm(edge_vector)

        # create the rotation matrix
        # r = Rotation.from_rotvec(angle * edge_vector)

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
        # node_coords_rotated = r.apply(node_coords)
        node_coords_rotated = base.rotate_coords(
            node_coords - node_coords[idx_2], angle, edge_vector
        )
        node_coords_rotated += node_coords[idx_2]

        # # now adjust for the translatisonal shift around the axis
        # _diff = node_coords_rotated[idx_2] - node_coords[idx_2]
        # node_coords_rotated -= _diff

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
    import buildamol as bam
    from timeit import timeit

    import seaborn as sns
    import matplotlib.pyplot as plt

    from functools import partial

    mol = bam.molecule(
        "/Users/noahhk/GIT/biobuild/buildamol/optimizers/_testing/files/EX8.json"
    )
    v = mol.draw()
    g = BaseGraph(None, mol.bonds)
    nx.set_edge_attributes(g, 1, "bond_order")
    _g = mol.get_atom_graph()

    _g.sample_rotatable_edges = partial(BaseGraph.sample_rotatable_edges, _g)
    edges = _g.sample_edges(_g.find_rotatable_edges(min_descendants=10), n=4, m=3)

    v.draw_edges(*edges, color="limegreen", linewidth=8, elongate=1.1)
    v.show()
    pass
    #  ref_atom_graph = mol.get_atom_graph()

    # a, b = mol.get_atoms(68, 65)

    # x = g.get_descendants(a, b)
    # for i in x:
    #     v.draw_atom(i, color="purple")

    # measure_performance = True
    # repeats = 500
    # number = 800
    # if measure_performance:
    #     test_old = lambda: BaseGraph.get_descendants_old(ref_atom_graph, a, b)
    #     test_new = lambda: g.get_descendants(a, b)
    #     # test_2 = lambda: g.get_descendants_2(a, b)

    #     times_old = [timeit(test_old, number=number) for _ in range(repeats)]
    #     times_new = [timeit(test_new, number=number) for _ in range(repeats)]
    #     # times_2 = [timeit(test_2, number=number) for _ in range(repeats)]

    #     sns.distplot(times_old, label="old", kde=True, bins=20)
    #     sns.distplot(times_new, label="new", kde=True, bins=20)
    #     # sns.distplot(times_2, label="2", kde=True, bins=20)
    #     plt.legend()
    #     plt.show()

    # pass

# v.draw_edges((a, b), color="limegreen", linewidth=4, elongate=1.1)
# ref_decendants = ref_atom_graph.get_descendants(a, b)

# descendants = g.get_descendants_2(a, b)

# for i in descendants:
#     v.draw_atom(i, color="purple")
# for i in ref_decendants:
#     v.draw_atom(i, color="orange")
# v.show()
# pass
