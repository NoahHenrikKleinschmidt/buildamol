import networkx as nx
import Bio.PDB as bio

import biobuild.structural as struct
from biobuild.graphs.BaseGraph import BaseGraph


class AtomGraph(BaseGraph):
    """
    A graph representation of atoms and bonds in a contiguous molecule.
    """

    def __init__(self, id, bonds: list):
        super().__init__(id, bonds)
        nx.set_node_attributes(self, {i: i.coord for i in self.nodes}, "coord")

    @classmethod
    def from_biopython(
        cls,
        structure,
        apply_standard_bonds: bool = True,
        infer_residue_connections: bool = True,
        infer_bonds: bool = False,
        max_bond_length: float = None,
        restrict_residues: bool = True,
        _topology=None,
    ):
        """
        Create an AtomGraph from a biopython structure

        Parameters
        ----------
        structure
            The biopython structure. This can be any biopython object that houses atoms.
        infer_residue_connections: bool
            Whether to infer residue connecting bonds based on atom distances.
        infer_bonds : bool
            Whether to infer bonds from the distance between atoms. If this is set
            to True, standard bonds cannot be also applied!
        max_bond_length : float
            The maximum distance between atoms to infer a bond.
            If none is given, a default bond length is assumed.
        restrict_residues : bool
            Whether to restrict to atoms of the same residue when inferring bonds.
            If set to False, this will also infer residue connecting bonds.
        _topology
            A specific reference topology to use when re-constructing any missing parts.
            By default the default CHARMM topology is used.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        while structure.level != "S":
            _parent = structure.get_parent()
            if _parent is None:
                break
            structure = _parent

        bonds = cls._make_bonds(
            structure,
            apply_standard_bonds,
            infer_residue_connections,
            infer_bonds,
            max_bond_length,
            restrict_residues,
            _topology,
        )
        return cls(structure.id, bonds)

    @classmethod
    def from_molecule(cls, mol, locked: bool = False):
        """
        Create an AtomGraph from a molecule

        Parameters
        ----------
        mol : biobuild.molecule.Molecule
            The molecule to convert
        locked : bool, optional
            If True, any information about locked bonds will also be
            transferred to the AtomGraph, by default False.

        Returns
        -------
        AtomGraph
            The AtomGraph representation of the molecule
        """
        new = cls(mol.id, mol.bonds)
        if locked:
            new._locked_edges.update(mol.locked_bonds)
        return new

    def migrate_bonds(self, other):
        """
        Migrate bonds from another graph

        Parameters
        ----------
        other : AtomGraph
            The other graph to migrate bonds from
        """
        self._locked_edges.update(other._locked_edges)
        bond_orders = nx.get_edge_attributes(other, "bond_order")
        self.add_edges_from(other.edges)
        nx.set_edge_attributes(self, bond_orders, "bond_order")

    def get_neighbors(self, atom: bio.Atom.Atom, n: int = 1, mode="upto"):
        """
        Get the neighbors of a node

        Parameters
        ----------
        atom : bio.Atom.Atom
            The atom
        n : int, optional
            The number of bonds to separate the atom from its neighbors.

        mode : str, optional
            The mode to use for getting the neighbors, by default "upto"
            - "upto": get all neighbors up to a distance of `n` bonds
            - "exact": get all neighbors exactly `n` bonds away

        Returns
        -------
        set
            The neighbors of the atom
        """
        if not self._neighborhood:
            self._neighborhood = struct.AtomNeighborhood(self)
        return self._neighborhood.get_neighbors(atom, n, mode)

    def direct_edges(self, root_node, edges=None):
        """
        Sort the edges such that the first atom in each edge
        is the one with the lower serial number.

        Parameters
        ----------
        root_node
            The root node to use for sorting the edges
        edges : list, optional
            The edges to sort, by default None, in which case
            all edges are sorted.

        Returns
        -------
        list
            The sorted edges
        """
        if edges is None:
            edges = list(self.edges)

        if root_node not in self.nodes:
            raise ValueError(f"Root node {root_node} not in graph")

        _directed = []
        for node1, node2 in edges:
            d1 = nx.shortest_path_length(self, source=root_node, target=node1)
            d2 = nx.shortest_path_length(self, source=root_node, target=node2)
            if d1 > d2:
                edge = (node2, node1)
            else:
                edge = (node1, node2)
            _directed.append(edge)
        return _directed

    @staticmethod
    def _make_bonds(
        structure,
        apply_standard_bonds: bool,
        infer_residue_connections: bool,
        infer_bonds: bool,
        max_bond_length: float,
        restrict_residues: bool,
        _topology=None,
    ) -> list:
        """
        Make bond tuples from a structure
        """
        bonds = []
        if infer_bonds:
            apply_standard_bonds = False
            infer_residue_connections = (
                infer_residue_connections and not restrict_residues
            )
            bonds.extend(
                struct.infer_bonds(structure, max_bond_length, restrict_residues)
            )

        if apply_standard_bonds:
            bonds.extend(struct.apply_standard_bonds(structure, _topology))

        if infer_residue_connections:
            bonds.extend(struct.infer_residue_connections(structure, max_bond_length))

        return bonds


if __name__ == "__main__":
    _man = "support/examples/man9.pdb"
    man = AtomGraph.from_pdb(_man)

    import sys

    # print(man.nodes)
    # print(man.edges)
    # # x = utils.infer_residue_connections(man.structure)
    # import matplotlib.pyplot as plt
    # import networkx as nx

    # # nx.draw(man)
    # colors = {
    #     "C": "gray",
    #     "O": "red",
    #     "N": "blue",
    #     "S": "yellow",
    #     "P": "orange",
    #     "H": "lightgray",
    #     "F": "green",
    # }
    # _colors = [colors[i.element] for i in man.nodes]
    # g = man
    # # y = [(i.get_parent(), j.get_parent()) for i, j in x]

    # # g = nx.Graph(y)
    # # colors = {"MAN": "green", "BMA": "cyan", "MNA": "orange", "GAL": "pink", "NAG": "gray"}
    # # _colors = [colors[i.resname] for i in g.nodes]
    # nx.draw(g, node_color=_colors)
    # plt.show()
