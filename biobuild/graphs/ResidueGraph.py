import networkx as nx
import Bio.PDB as bio

import biobuild.structural as struct
from biobuild.graphs.BaseGraph import BaseGraph
import biobuild.utils.visual as vis


class ResidueGraph(BaseGraph):
    """
    A graph representation of residues bonded together as an abstraction of a large contiguous molecule.
    """

    __idx_method__ = (
        lambda x: x.serial_number if hasattr(x, "serial_number") else x.id[1]
    )

    def __init__(self, id, bonds: list):
        super().__init__(id, bonds)
        self._AtomGraph = None
        self._atomic_bonds = {}
        self._atomic_bonds_list = []

        self._residues = {i.id: i for i in self.nodes}
        for r in self.nodes:
            r.coord = r.center_of_mass()

        nx.set_node_attributes(
            self, {i: i.center_of_mass() for i in self.nodes}, "coord"
        )

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
            raise ValueError(
                "Molecule must have at least 2 residues to make a ResidueGraph"
            )

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

        new._atomic_bonds_list = list(connections)

        for bond in connections:
            parent_1 = bond[0].get_parent()
            parent_2 = bond[1].get_parent()
            new._atomic_bonds.setdefault((parent_1, parent_2), []).append(bond)

        if detailed:
            new.make_detailed()
            if locked:
                new._locked_edges.update(
                    (i for i in mol._AtomGraph._locked_edges if i in new.edges)
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

    def make_detailed(
        self,
        include_outliers: bool = False,
        include_heteroatoms: bool = False,
        f: float = 1.5,
    ):
        """
        Use a detailed representation of the residues in the molecule by adding the specific atoms
        that connect the residues together. This is useful for visualization and analysis.

        Note
        ----
        This function is not reversible.

        Parameters
        ----------
        include_outliers : bool
            If True, atoms that are not involved in residue connections are also included if their
            distance to the residue's center of mass is greater than f * the 75th percentile of
            atom distances to the residue's center of mass.

        include_heteroatoms : bool
            If True, all hetero-atoms are included in the detailed representation, regardless of
            their distance to the residue center of mass.

        f : float
            The factor by which the 75th percentile of atom distances to the residue's center of mass
            is multiplied to determine the cutoff distance for outlier atoms. This is only used if
            include_outliers is True.
        """

        self.clear_edges()

        _added_nodes = set()
        for edge in self._atomic_bonds_list:
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

        if include_outliers:
            for residue in self.residues:
                outliers = struct.compute_outlier_atoms(residue, f=f)
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

    def draw(self):
        v = vis.ResidueGraphViewer3D()
        v.link(self)
        return v

    def direct_edges(self):
        """
        Sort the edges such that the first atom in each edge
        is the one with the lower serial number.
        """
        for edge in self.edges:
            if self.__idx_method__(edge[0]) > self.__idx_method__(edge[1]):
                self.remove_edge(*edge)
                self.add_edge(edge[1], edge[0])

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
        return list(self._residues.values())

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

    def centers_of_mass(self):
        """
        Get the centers of mass of the residues in the molecule.

        Returns
        -------
        dict
            The centers of mass of the residues in the molecule
        """
        return {residue.id: residue.center_of_mass() for residue in self.residues}


if __name__ == "__main__":
    import biobuild as bb

    f = "support/examples/man9.pdb"
    mol = bb.Molecule.from_pdb(f)
    mol.infer_bonds(restrict_residues=False)
    man = ResidueGraph.from_molecule(mol)

    # _man = "support/examples/MAN9.pdb"
    # _man = bb.Molecule.from_pdb(_man)
    # _man.infer_bonds(restrict_residues=False)
    # man = ResidueGraph.from_molecule(_man)

    import matplotlib.pyplot as plt
    import networkx as nx

    man.make_detailed(True, f=1)

    import biobuild.utils.visual as vis

    v = vis.MoleculeViewer3D(man)
    v.show()

    print(len(list(man.bonds)))
    nx.draw(man, with_labels=True, font_weight="bold", pos=nx.spectral_layout(man))
    plt.show()
