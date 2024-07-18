"""
Classes to handle the atom and residue neighborhoods of a structure
"""

import numpy as np

import buildamol.structural.base as base


class Neighborhood:
    """
    This class handles graph connectivity of connected nodes and is the basis
    for the AtomNeighborhood and ResidueNeighborhood classes.

    Parameters
    ----------
    graph
        A networkx graph
    """

    # the method of node objects to use for mapping to an index-based dictionary
    # this assumes that the index is a **unique** identifier for the nodes
    __index_method__ = AttributeError

    # the method of node objects to use for mapping to an id-based dictionary
    # this allows for multiple nodes with the same id (e.g. multiple 'C1' atoms...)
    __id_method__ = AttributeError

    def __init__(self, graph):
        self._src = graph
        self.adj = graph.adj

        # each implementation may call this method before the super init
        # self._validate_getters()

        # self._make_node_dicts()

    def get_neighbors(self, node, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of an atom

        Parameters
        ----------
        node : object
            The target node object.
        n : int
            The (maximal) number of edges that should
            separate the target from neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n edges away from the target or closer.
            If "at", this will get all neighbors that are exactly n edges
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the target node
        """
        if isinstance(node, list):
            return [self.get_neighbors(a, n, mode) for a in node]

        self._seen = set()
        if mode == "upto":
            return self._get_neighbors_upto(node, n) - {node}
        elif mode == "at":
            return self._get_neighbors_at(node, n) - {node}
        else:
            raise ValueError(f"Invalid mode: {mode}")

    def _get_neighbors_upto(self, node, n: int):
        """
        Get all neighbors of a node that are n edges away from the target or closer
        """
        if n == 0:
            return {node}
        else:
            neighbors = set()
            for neighbor in self._src.adj[node]:
                if neighbor in self._seen:
                    continue
                if n >= 1:
                    neighbors.update(self._get_neighbors_upto(neighbor, n - 1))
                neighbors.add(neighbor)
                self._seen.add(neighbor)
            return neighbors

    def _get_neighbors_at(self, node, n: int):
        """
        Get all neighbors of a node that are exactly n edges away from the target
        """
        if n == 0:
            return {node}
        else:
            neighbors = set()
            for neighbor in self._src.adj[node]:
                if neighbor in self._seen:
                    continue
                if n >= 1:
                    neighbors.update(self._get_neighbors_upto(neighbor, n - 1))
                self._seen.add(neighbor)
            return neighbors

    def search_by_constraints(self, constraints: list):
        """
        Search for nodes that satisfy a set of constraints

        Parameters
        ----------
        constraints : list
            A list of constraint functions. Each function should take
            a graph and a node as arguments and return a boolean.

        Returns
        -------
        nodes : set
            The nodes that satisfy the constraints
        """
        matches = []
        if len(constraints) == 1:
            node_factory = lambda node: {node}
        else:
            node_factory = lambda node: [
                node,
                *self.get_neighbors(node, n=len(constraints) - 1),
            ]

        for node in self._src.nodes:
            _nodes = node_factory(node)
            m = {}
            for i, constraint in enumerate(constraints):
                for n in _nodes:
                    if constraint is None or constraint(self, n):
                        m[i] = n
                        _nodes.remove(n)
                        break
            if len(m) == len(constraints):
                if m in matches:
                    continue
                matches.append(m)
        return matches


class AtomNeighborhood(Neighborhood):
    """
    This class handles the bond connectivity neighborhood of atoms in a structure
    and can be used to obtain bonded atom triplets.

    Parameters
    ----------
    graph
        An AtomGraph
    """

    __index_method__ = __index_method__ = lambda _, node: getattr(node, "serial_number")
    __id_method__ = lambda _, node: getattr(node, "id")

    @property
    def atoms(self):
        """
        Returns the atoms in the structure
        """
        return self._src.atoms

    @property
    def bonds(self):
        """
        Returns the bonds in the structure
        """
        return self._src.bonds

    def get_neighbors(self, atom, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of an atom

        Parameters
        ----------
        atom : Bio.PDB.Atom
            The atom whose neighbors should be returned.
        n : int
            The (maximal) number of bonds that should
            separate the target from the neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n bonds away from the target or closer.
            If "at", this will get all neighbors that are exactly n bonds
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the atom
        """
        if atom is None:
            return set()
        return super().get_neighbors(atom, n, mode)


class ResidueNeighborhood(Neighborhood):
    """
    This class handles the residue connectivity neighborhood of residues in a structure
    and can be used to obtain residue triplets.

    Parameters
    ----------
    graph
        A ResidueGraph
    """

    # a maybe little hacky way to get the residue name or atom id, depending on the node type
    # since we have mixed types in detailed residue graphs
    __index_method__ = __index_method__ = lambda _, node: (
        node.id[1] if isinstance(node.id, tuple) else node.serial_number
    )
    __id_method__ = lambda _, node: (
        getattr(node, "resname")
        if hasattr(node, "resname")
        else f"{node.id}@{node.serial_number}"
    )

    @property
    def residues(self):
        """
        Returns the residues in the structure
        """
        return self._src.residues

    @property
    def bonds(self):
        """
        Returns the bonds in the structure
        """
        return self._src.bonds

    def get_neighbors(self, residue, n: int = 1, mode: str = "upto"):
        """
        Get the neighbors of a residue

        Parameters
        ----------
        residue : Bio.PDB.Residue
            The residue whose neighbors should be returned.
        n : int
            The (maximal) number of bonds that should
            separate the target from the neighbors.
        mode : str
            The mode of the neighbor search. Can be either
            "upto" or "at". If "upto", this will get all
            neighbors that are n bonds away from the target or closer.
            If "at", this will get all neighbors that are exactly n bonds
            away from the target.

        Returns
        -------
        neighbors : set
            The neighbors of the residue
        """
        if residue is None:
            return set()
        return super().get_neighbors(residue, n, mode)


class Quartet:
    """
    An atom quartet that can be used to compute internal coordinates
    """

    def __init__(self, atom1, atom2, atom3, atom4, improper: bool = False) -> None:
        self._atoms = (atom1, atom2, atom3, atom4)
        self._improper = improper

    @property
    def atoms(self):
        return self._atoms

    @property
    def atom1(self):
        return self._atoms[0]

    @property
    def atom2(self):
        return self._atoms[1]

    @property
    def atom3(self):
        return self._atoms[2]

    @property
    def atom4(self):
        return self._atoms[3]

    @property
    def improper(self):
        return self._improper

    @property
    def center_atom(self):
        if self._improper:
            return self._atoms[2]
        else:
            raise TypeError("This quartet is not an improper")

    @property
    def dist_12(self):
        return np.linalg.norm(self._atoms[1].coord - self._atoms[0].coord)

    @property
    def dist_23(self):
        return np.linalg.norm(self._atoms[2].coord - self._atoms[1].coord)

    @property
    def dist_34(self):
        return np.linalg.norm(self._atoms[3].coord - self._atoms[2].coord)

    @property
    def dist_13(self):
        return np.linalg.norm(self._atoms[2].coord - self._atoms[0].coord)

    @property
    def angle_123(self):
        return base.compute_angle(self._atoms[0], self._atoms[1], self._atoms[2])

    @property
    def angle_234(self):
        return base.compute_angle(self._atoms[1], self._atoms[2], self._atoms[3])

    @property
    def dihedral(self):
        return base.compute_dihedral(
            self._atoms[0], self._atoms[1], self._atoms[2], self._atoms[3]
        )

    def __hash__(self) -> int:
        return hash(tuple(sorted(self._atoms))) + hash(self._improper)

    def __eq__(self, other) -> bool:
        if isinstance(other, Quartet):
            return (
                set(self._atoms) == set(other._atoms)
                and self._improper == other._improper
            )
        elif isinstance(other, (tuple, list)):
            if len(other) == 4:
                return set(self._atoms) == set(other)
            elif len(other) == 5:
                return set(self._atoms) == set(other[:4]) and self._improper == other[4]
        return False

    def __repr__(self) -> str:
        if hasattr(self._atoms[0], "id"):
            return f"Quartet({self._atoms[0].id}, {self._atoms[1].id}, {self._atoms[2].id}, {self._atoms[3].id}, improper={self._improper})"
        else:
            return f"Quartet({self._atoms[0]}, {self._atoms[1]}, {self._atoms[2]}, {self._atoms[3]}, improper={self._improper})"

    def __iter__(self):
        return iter(self._atoms)

    def __getitem__(self, item):
        return self._atoms[item]


def compute_triplets(bonds: list, unique: bool = True):
    """
    Compute all possible triplets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds
    unique : bool
        Whether to return only unique triplets. If False, triplets (1,2,3) and (3,2,1) will be returned.
        Otherwise only one of them will be returned.

    Returns
    -------
    triplets : list
        A list of triplets

    Examples
    --------
    ```
    ( 1 )---( 2 )---( 4 )
       \\
       ( 3 )
         |
       ( 5 )
    ```
    >>> bonds = [(1, 2), (1, 3), (2, 4), (3, 5)]
    >>> compute_triplets(bonds)
    [(2, 1, 3), (1, 2, 4), (1, 3, 5)]
    """
    triplets = list(generate_triplets(bonds))
    if unique:
        half_length = len(triplets) // 2
        while len(triplets) != half_length:
            triplet = triplets.pop()
            if triplet[::-1] in triplets:
                continue
            else:
                triplets.insert(0, triplet)
    return triplets

    # triplets = []
    # for i, bond1 in enumerate(bonds):
    #     atom_11, atom_12 = bond1
    #     for j, bond2 in enumerate(bonds[i + 1 :]):
    #         atom_21, atom_22 = bond2
    #         if atom_11 == atom_22:
    #             continue
    #         if atom_11 == atom_21:
    #             triplets.append((atom_12, atom_11, atom_22))
    #         elif atom_11 == atom_22:
    #             triplets.append((atom_12, atom_11, atom_21))
    #         elif atom_12 == atom_21:
    #             triplets.append((atom_11, atom_12, atom_22))
    #         elif atom_12 == atom_22:
    #             triplets.append((atom_11, atom_12, atom_21))
    # return triplets


def compute_quartets(bonds: list):
    """
    Compute all possible quartets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds

    Returns
    -------
    quartets : list
        A list of quartets

    Examples
    --------
    ```
    ( 1 )---( 2 )---( 4 )
               \\
               ( 3 )
                 |
               ( 5 )
    ```
    >>> bonds = [(1, 2), (2, 3), (2, 4), (3, 5)]
    >>> compute_quartets(bonds)
    {Quartet(1, 2, 3, 5, improper=False), Quartet(5, 3, 2, 4, improper=False), Quartet(1, 3, 2, 4, improper=True)}
    """

    triplets = compute_triplets(bonds)
    quartets = set()
    for idx, triplet1 in enumerate(triplets):
        atom_1, atom_2, atom_3 = triplet1

        for triplet2 in triplets:  # [idx + 1 :]:
            atom_4, atom_5, atom_6 = triplet2

            # decision tree to map atoms into quartets
            if (
                triplet1 is triplet2
            ):  # all((atom_1 is atom_4, atom_2 is atom_5, atom_3 is atom_6)):
                continue

            quartet = None
            # ------------------- #
            # NEW IMPLEMENTATION  #
            # ------------------- #
            if atom_1 is atom_4:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            elif atom_2 is atom_4:
                if atom_1 is atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)
            elif atom_3 is atom_4:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            ###
            elif atom_1 is atom_6:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            elif atom_2 is atom_6:
                if atom_1 is atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)
            elif atom_3 is atom_6:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_4, atom_2, atom_3, True)
            ###
            # ------------------- #
            # this needs overhauling, there seems some logic error in there somewhere...
            # if atom_2 is atom_4:
            #     if atom_1 is atom_5:
            #         quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)

            #     elif atom_3 is atom_5:
            #         quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)

            # elif atom_2 is atom_5:
            #     if atom_1 is atom_4 or atom_3 is atom_4:
            #         quartet = Quartet(atom_1, atom_3, atom_2, atom_6, True)

            # elif atom_2 is atom_6:
            #     if atom_1 is atom_5:
            #         quartet = Quartet(atom_6, atom_1, atom_4, atom_3, False)

            #     elif atom_3 is atom_5:
            #         quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)

            if quartet:
                if len(set(quartet.atoms)) == 4:
                    quartets.add(quartet)

    return quartets


def generate_triplets(bonds: list):
    """
    Compute all possible triplets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds

    Returns
    -------
    triplets : generator
        A generator of triplets

    Examples
    --------
    ```
    ( 1 )---( 2 )---( 4 )
       \\
       ( 3 )
         |
       ( 5 )
    ```
    >>> bonds = [(1, 2), (1, 3), (2, 4), (3, 5)]
    >>> list(generate_triplets(bonds))
    [(2, 1, 3), (3, 1, 2), (1, 2, 4), (4, 2, 1), (1, 2, 4)]
    """
    for bond1 in bonds:
        atom_11, atom_12 = bond1
        for bond2 in bonds:
            atom_21, atom_22 = bond2

            # we used to compare with == (Which works perfectly fine)
            # but since we use this function to generate atom bond triplets
            # where the atoms are not only supposed to be equal but should literally be the same object
            # we can use the is operator to speed up the process (hopefully.)
            # UPDATE: We use == again because we can now have multiple copies
            # of the same bond (to represent double bonds etc.)
            if bond1 == bond2:
                continue
            if atom_11 is atom_21:
                yield (atom_12, atom_11, atom_22)
            elif atom_11 is atom_22:
                yield (atom_12, atom_11, atom_21)
            elif atom_12 is atom_21:
                yield (atom_11, atom_12, atom_22)
            elif atom_12 is atom_22:
                yield (atom_11, atom_12, atom_21)


def generate_quartets(bonds: list):
    """
    Generate all possible quartets of atoms from a list of bonds.

    Parameters
    ----------
    bonds : list
        A list of bonds

    Yields
    ------
    quartet : Quartet
        A quartet of atoms
    """
    triplets = compute_triplets(bonds)
    for triplet1 in triplets:
        atom_1, atom_2, atom_3 = triplet1

        for triplet2 in triplets:
            atom_4, atom_5, atom_6 = triplet2

            # decision tree to map atoms into quartets
            if triplet1 is triplet2:
                continue

            quartet = None
            # ------------------- #
            # NEW IMPLEMENTATION  #
            # ------------------- #
            if atom_1 is atom_4:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            elif atom_2 is atom_4:
                if atom_1 is atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)
            elif atom_3 is atom_4:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            ###
            elif atom_1 is atom_6:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_6, atom_2, atom_3, True)
            elif atom_2 is atom_6:
                if atom_1 is atom_5:
                    quartet = Quartet(atom_6, atom_1, atom_2, atom_3, False)
                elif atom_3 is atom_5:
                    quartet = Quartet(atom_1, atom_2, atom_3, atom_6, False)
            elif atom_3 is atom_6:
                if atom_2 is atom_5:
                    quartet = Quartet(atom_1, atom_4, atom_2, atom_3, True)

            if quartet:
                if len(set(quartet.atoms)) == 4:
                    yield quartet


class constraints:
    """
    Neighborhood structural constraints
    """

    none = lambda graph, node: True
    """
    No constraints
    """

    has_element = lambda element: (
        lambda graph, node: node.element.lower() == element.lower()
    )
    """
    The node has the specified element
    """

    has_any_element = lambda *args: (
        lambda graph, node: node.element.lower() in (i.lower() for i in args)
    )
    """
    The node has any of the specified elements
    """

    has_not_element = lambda element: (
        lambda graph, node: node.element.lower() != element.lower()
    )
    """
    The node does not have the specified element
    """

    has_id = lambda id: lambda graph, node: node.id == id
    """
    The node has the specified id
    """

    has_any_id = lambda *args: lambda graph, node: node.id in args
    """
    The node has any of the specified ids
    """

    has_not_id = lambda id: lambda graph, node: node.id != id
    """
    The node does not have the specified id
    """

    is_residue = lambda graph, node: hasattr(node, "resname")
    """
    The node is a residue
    """

    is_atom = lambda graph, node: hasattr(node, "element")
    """
    The node is an atom
    """

    alone = lambda graph, node: len(graph.adj[node]) == 0
    """
    The node has no neighbors
    """

    neighbors_any = lambda *args: (
        lambda graph, node: any(i in (j.element for j in graph.adj[node]) for i in args)
    )
    """
    The node has at least one neighbor with any of the specified elements
    """

    neighbors_all = lambda *args: (
        lambda graph, node: all(i in (j.element for j in graph.adj[node]) for i in args)
    )
    """
    The node neighbors all the specified elements but may have more
    """

    neighbors_exactly = lambda *args: (
        lambda graph, node: set(j.element for j in graph.adj[node]) == set(args)
    )
    """
    The node neighbors all and only the specified elements
    """

    neighbors_not = lambda *args: (
        lambda graph, node: all(
            i not in (j.element for j in graph.adj[node]) for i in args
        )
    )
    """
    The node does not have any neighbors with the specified elements
    """

    has_neighbor_hist = lambda hist: (
        lambda graph, node: all(
            hist[i] == sum(1 for j in graph.adj[node] if j.element == i) for i in hist
        )
    )
    """
    The node has the specified number of neighbors for each element
    Hist is a dictionary with element symbols as keys and the number of neighbors as values
    """

    has_n_neighbors = lambda n: (lambda graph, node: len(graph.adj[node]) == n)
    has_not_n_neighbors = lambda n: (lambda graph, node: len(graph.adj[node]) != n)
    has_at_least_n_neighbors = lambda n: (lambda graph, node: len(graph.adj[node]) >= n)
    has_at_most_n_neighbors = lambda n: (lambda graph, node: len(graph.adj[node]) <= n)

    extended_neighbors_any = lambda n, *args: (
        lambda graph, node: any(
            i in (j.element for j in graph.get_neighbors(node, n)) for i in args
        )
    )
    """
    The node has at least one neighbor with any of the specified elements within n bonds
    """

    extended_neighbors_all = lambda n, *args: (
        lambda graph, node: all(
            i in (j.element for j in graph.get_neighbors(node, n)) for i in args
        )
    )
    """
    The node neighbors all the specified elements within n bonds but may have more
    """

    extended_neighbors_exactly = lambda n, *args: (
        lambda graph, node: set(j.element for j in graph.get_neighbors(node, n))
        == set(args)
    )
    """
    The node neighbors all and only the specified elements within n bonds
    """

    extended_neighbors_not = lambda n, *args: (
        lambda graph, node: all(
            i not in (j.element for j in graph.get_neighbors(node, n)) for i in args
        )
    )
    """
    The node does not have any neighbors with the specified elements within n bonds
    """

    extended_has_neighbor_hist = lambda n, hist: (
        lambda graph, node: all(
            hist[i] == sum(1 for j in graph.get_neighbors(node, n) if j.element == i)
            for i in hist
        )
    )
    """
    The node has the specified number of neighbors for each element within n bonds
    Hist is a dictionary with element symbols as keys and the number of neighbors as values
    """

    extended_has_n_neighbors = lambda n, m: (
        lambda graph, node: len(graph.get_neighbors(node, n)) == m
    )
    extended_has_not_n_neighbors = lambda n, m: (
        lambda graph, node: len(graph.get_neighbors(node, n)) != m
    )
    extended_has_at_least_n_neighbors = lambda n, m: (
        lambda graph, node: len(graph.get_neighbors(node, n)) >= m
    )

    extended_has_at_most_n_neighbors = lambda n, m: (
        lambda graph, node: len(graph.get_neighbors(node, n)) <= m
    )

    def multi_constraint(*funcs):
        """
        Combine multiple constraints into one

        Parameters
        ----------
        funcs : list
            A list of constraint functions to combine.
            Each of these must take a graph and a node as arguments and return a boolean.
        """
        return lambda graph, node: all(f(graph, node) for f in funcs)
