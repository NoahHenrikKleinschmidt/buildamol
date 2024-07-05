"""
The Circulatron class to circularize a molecule
"""

import numpy as np
import buildamol.optimizers.base_rotatron as Rotatron
from buildamol.optimizers.distance_rotatron import DistanceRotatron

__all__ = ["Circulatron"]


class Circulatron(Rotatron.Rotatron):
    """
    A special rotatron to circularize a molecule.
    This rotatron works by minimizing the distance between two target nodes in the graph in order to superimpose them.
    This rotatron does not itself optimize the conformation of the resulting structure but instead uses one of the other
    rotatrons to do so.

    Note
    ----
    In order for this environment to work, it is important that the graph is NOT already circularized!

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to circulize
    target_nodes : tuple
        Two nodes in the graph that should be superimposed.
    base_rotatron : Rotatron
        The rotatron class to use for basic conformer evaluation. By default, this is the DistanceRotatron.
    rotatable_edges : list
        A list of edges that are rotatable in the graph. If not given, all rotatable edges will be considered.
    **kwargs
        Additional keyword arguments to pass to the base rotatron
    """

    def __init__(
        self,
        graph,
        target_nodes: tuple,
        base_rotatron: Rotatron = DistanceRotatron,
        hinge_node=None,
        rotatable_edges: list = None,
        **kwargs
    ):
        hinge_node = hinge_node or target_nodes[0]
        if rotatable_edges is None:
            rotatable_edges = graph.find_rotatable_edges(hinge_node)
        else:
            rotatable_edges = graph.direct_edges(hinge_node, rotatable_edges)

        super().__init__(graph, rotatable_edges=rotatable_edges, **kwargs)
        self.target_nodes = target_nodes
        self.tdx1 = list(graph.nodes).index(target_nodes[0])
        self.tdx2 = list(graph.nodes).index(target_nodes[1])
        self.base_rotatron = base_rotatron(
            graph, rotatable_edges, setup=False, **kwargs
        )
        base_rotatron.edge_masks = self.edge_masks
        base_rotatron.edge_lengths = self.edge_lengths

        self.action_space = self.base_rotatron.action_space
        self._bounds_tuple = self.base_rotatron._bounds_tuple

    def eval(self, state):
        """
        Evaluate the state using the base rotatron

        Parameters
        ----------
        state : dict
            The state to evaluate

        Returns
        -------
        float
            The energy of the state
        """
        e = self.base_rotatron.eval(state)
        # now evaluate the distance between the target nodes
        self._target_dist = np.linalg.norm(state[self.tdx1] - state[self.tdx2])
        e += self._target_dist**2
        return e

    def done(self, state):
        """
        Check if the state is done

        Parameters
        ----------
        state : dict
            The state to check

        Returns
        -------
        bool
            Whether the state is done
        """
        return self._target_dist < 0.1


if __name__ == "__main__":
    import buildamol as bam

    bam.load_amino_acids()
    mol = bam.molecule("SER") % "LINK" * 15

    targets = mol.get_atom("N", residue=1), mol.get_atom("HXT", residue=-1)
    hinge = mol.get_atom("CA", residue=5)

    res_graph = mol.get_residue_graph(True)
    res_graph.add_edge(targets[0], targets[0].parent)
    res_graph.add_edge(targets[1], targets[1].parent)
    res_graph.add_edge(hinge, hinge.parent)

    atom_graph = mol.get_atom_graph()

    edges = mol.get_residue_connections()

    graph = atom_graph

    rotatron = Circulatron(
        graph,
        targets,
        hinge_node=hinge,
        base_rotatron=bam.DistanceRotatron,
        rotatable_edges=edges,
        unfold=4.0,
    )

    v = mol.draw(atoms=False)
    v.draw_points(mol.get_coords(*targets), colors="limegreen")
    v.show()

    out = bam.optimize(mol, rotatron, algorithm="scipy")

    out.remove_atoms(targets[1])
    out.add_bond(targets[0], mol.get_atom("OXT", residue=-1))
    v += out.draw(line_color="red", atoms=False)

    out.to_pdb("circular.pdb")
    v.show()
