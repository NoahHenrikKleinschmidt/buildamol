"""
This is the basic Rotatron environment. It provides the basic functionality for preprocessing a graph into numpy arrays, masking rotatable edges, and evaluating a possible solution.
All other Rotatron environments inherit from this class.
"""

import gym
import numpy as np

from scipy.spatial.distance import cdist

import buildamol.utils.auxiliary as aux
import buildamol.graphs.base_graph as base_graph
import buildamol.structural.base as structural
from copy import deepcopy

from multiprocessing import Pool

__all__ = ["Rotatron"]


class Rotatron(gym.Env):
    """
    The base class for rotational optimization environments.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    n_processes : int
        The number of processes to use to speed up the computation of edge masks and lengths
    setup : bool
        Whether to set up the edge masks and lengths during initialization
    numba : bool
        Whether to use numba to speed up the rotation function.
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        n_processes: int = 1,
        setup: bool = True,
        numba: bool = False,
        **kwargs
    ):
        self.graph = graph
        self.rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)

        self.node_dict = {n: i for i, n in enumerate(self.graph.nodes)}
        self.n_nodes = len(self.node_dict)
        self.n_edges = len(self.rotatable_edges)

        self.state = self._make_state_from_graph(self.graph).astype(np.float64)
        self._backup_state = self.state.copy()

        self.rotation_unit_masks = np.ones(
            (len(graph.nodes), len(graph.nodes)), dtype=bool
        )

        self.edge_lengths = np.zeros(self.n_edges)
        self.edge_masks = np.zeros((self.n_edges, self.n_nodes), dtype=bool)

        self.n_processes = n_processes

        if setup:
            self._generate_edge_masks(n_processes=n_processes)
            self._generate_edge_lengths()

        self._edge_node_coords = np.array(
            [[self.node_dict[e[0]], self.node_dict[e[1]]] for e in self.rotatable_edges]
        )

        if (
            numba
            or aux.USE_ALL_NUMBA
            or (self.n_edges * self.n_nodes > 10000 and aux.USE_NUMBA)
        ):
            self._rotate = self._numba_rotate
        else:
            self._rotate = self._normal_rotate

    def eval(self, state):
        """
        Calculate the evaluation score for a given state

        Parameters
        ----------
        state : np.ndarray
            The state of the environment

        Returns
        -------
        float
            The evaluation for the state
        """
        return np.inf

    def step(self, action):
        """
        Take a step in the environment

        Parameters
        ----------
        action : np.ndarray
            The action to take

        Returns
        -------
        np.ndarray
            The new state of the environment
        float
            The evaluation for the new state
        bool
            Whether the environment is done
        dict
            Additional information
        """
        new_state = self.state
        for edge in range(self.n_edges):
            new_state = self._rotate(
                new_state,
                edge,
                action[edge],
            )

        e = self.eval(new_state)
        done = self.is_done(new_state)
        return new_state, e, done, {}

    def is_done(self, state):
        """
        Check whether the environment is done

        Parameters
        ----------
        state : np.ndarray
            The state of the environment

        Returns
        -------
        bool
            Whether the environment is done
        """
        return False

    def reset(self, *args, **kwargs):
        """
        Reset the environment
        """
        self.state[:, :] = self._backup_state

    def blank(self):
        """
        A blank action
        """
        return np.zeros(len(self.rotatable_edges))

    def copy(self):
        """
        Make a deep copy of the environment
        """
        return deepcopy(self)

    def _make_state_from_graph(self, graph):
        """
        Set up the state of the environment
        """
        state = np.array([i.coord for i in graph.nodes])
        return state

    def _make_state_from_dict(self, dict):
        """
        Set up the state of the environment
        """
        state = np.array([v for v in dict.values()])
        return state

    def _get_rotatable_edges(self, graph, rotatable_edges):
        """
        Get the rotatable edges

        Parameters
        ----------
        graph : AtomGraph or ResidueGraph
            The graph to optimize
        rotatable_edges : list
            A list of edges that can be rotated during optimization.
            If None, all non-locked edges are used.

        Returns
        -------
        list
            The rotatable edges
        """
        if rotatable_edges is None:
            _circulars = graph.nodes_in_cycles
            rotatable_edges = [
                e
                for e in graph.edges
                if e not in graph._locked_edges
                and graph.edges[e].get("bond_order", 1) == 1
                and not "Residue" in type(e[0]).__name__
                and not "Residue" in type(e[1]).__name__
                and not (e[0] in _circulars and e[1] in _circulars)
                and len(graph.get_descendants(*e)) > 1
            ]
        return rotatable_edges

    def _generate_rotation_unit_masks(self):
        """
        Generate a boolean mask (n_nodes, n_nodes) where all nodes
        that are part of the same rotation unit are set to False
        """
        dists1 = cdist(self.state, self.state)
        for i, angle in enumerate(np.random.random(self.n_edges)):
            state2 = self._rotate(self.state, i, angle)
        dists2 = cdist(state2, state2)
        for i, angle in enumerate(np.random.random(self.n_edges)):
            state3 = self._rotate(state2, i, angle)
        dists3 = cdist(state3, state3)

        d12 = np.abs(dists1 - dists2) < 1e-4
        d13 = np.abs(dists1 - dists3) < 1e-4
        d23 = np.abs(dists2 - dists3) < 1e-4

        dists = np.sum([d12, d13, d23], axis=0) == 3
        self.rotation_unit_masks = ~dists
        self.reset()

    def _find_rotation_units(self):
        self.rotation_units = {}
        patterns = []
        rdx = 0
        for edx, mask in enumerate(self.rotation_unit_masks):
            pattern = ~mask
            if any(np.all(i == pattern) for i in patterns):
                pdx = next(
                    idx for idx, i in enumerate(patterns) if np.all(i == pattern)
                )
                self.rotation_units[pdx].add(edx)
                continue
            patterns.append(pattern)
            self.rotation_units[rdx] = {edx}
            rdx += 1
        self.rotation_units = {
            r: np.array(list(v)) for r, v in self.rotation_units.items()
        }

    def _generate_edge_lengths(self):
        """
        Compute the lengths of the edges
        """
        self.edge_lengths = np.array(
            [
                np.linalg.norm(
                    self.state[self.node_dict[e[0]]] - self.state[self.node_dict[e[1]]]
                )
                for e in self.rotatable_edges
            ]
        )

    def _generate_edge_masks(self, n_processes):
        """
        Compute the edge masks of downstream nodes
        """
        if n_processes > 1:
            p = Pool(n_processes)
            p.map(self._generate_edge_mask, [e for e in self.rotatable_edges])
            p.close()
            p.join()
        else:
            self.edge_masks = np.array(
                [self._generate_edge_mask(e) for e in self.rotatable_edges],
                dtype=bool,
            )

    def _generate_edge_mask(self, edge):
        return np.array(
            [
                1 if i in self.graph.get_descendants(*edge) else 0
                for i in self.graph.nodes
            ]
        )

    def _normal_rotate(self, state, edx, angle):
        if -1e-3 < angle < 1e-3:
            return self.state

        mask = self.edge_masks[edx]

        # vec = self._get_edge_vector(edx)
        adx, bdx = self._edge_node_coords[edx]
        vec = state[bdx] - state[adx]
        vec /= self.edge_lengths[edx]

        # ref_coord = self._get_edge_ref_coord(edx)
        ref_coord = state[adx]

        state[mask] = (
            structural.rotate_coords(state[mask] - ref_coord, angle, vec) + ref_coord
        )
        return state

    def _numba_rotate(self, state, edx, angle):
        if -1e-3 < angle < 1e-3:
            return self.state

        return _numba_wrapper_rotate(
            state,
            edx,
            angle,
            self._edge_node_coords,
            self.edge_lengths,
            self.edge_masks,
        )

    # ============================================================
    # The setup helper functions can be used by other environments
    # that inherit from this base class
    # ============================================================

    def _setup_helpers_crop_faraway_nodes(self, radius, graph=None, edges=None):
        """
        This is a helper function to remove nodes that are too far away from the rotatable edges.
        """
        if graph and edges:
            rotatable_edges = self._get_rotatable_edges(graph, edges)
        else:
            rotatable_edges = self.rotatable_edges

        edge_coords = np.array([(a.coord + b.coord) / 2 for a, b in rotatable_edges])

        nodes = list(graph.nodes)
        node_coords = np.array([node.coord for node in nodes])

        dists = cdist(edge_coords, node_coords)
        dists = dists > radius

        dists = np.apply_along_axis(np.all, 0, dists)

        if np.max(dists) != 0:
            nodes_to_drop = [nodes[i] for i, d in enumerate(dists) if d]
            graph.remove_nodes_from(nodes_to_drop)

        return graph, rotatable_edges


# ============================================================
# The numba functions are used to speed up things.
# For each function there must be a _normal_ and a _numba_
# version. The _normal_ version is used in the setup if numba is not installed
# The _numba_ version is used in the step function if numba is installed
# ============================================================


@aux.njit
def _numba_wrapper_rotate(
    state, edx, angle, edge_node_coords, edge_lengths, edge_masks
):
    """
    Rotate the graph around an edge. This is the version that is used in the step function.

    Parameters
    ----------
    edx : int
        The edge index to rotate around
    angle : float
        The angle to rotate by

    Returns
    -------
    np.ndarray
        The new state of the environment
    """
    mask = edge_masks[edx]
    adx, bdx = edge_node_coords[edx]
    vec = state[bdx] - state[adx]
    vec /= edge_lengths[edx]

    ref_coord = state[adx]

    rot = structural._numba_wrapper_rotation_matrix(vec, angle)
    rot = np.transpose(rot).astype(np.float64)

    state[mask] -= ref_coord
    _c = state[mask]
    _c = np.dot(_c, rot, out=_c)
    state[mask] = _c
    state[mask] += ref_coord

    return state


# if __name__ == "__main__":
#     import buildamol as bam

#     bam.load_sugars()
#     mol = bam.molecule("GLC") % "14bb"
#     mol *= 4

#     rot = Rotatron(mol.get_atom_graph(), n_processes=4)
#     rot._generate_rotation_unit_masks()
#     rot._find_rotation_units()


# if __name__ == "__main__":
#     import buildamol as bam

#     bam.load_sugars()
#     mol = bam.molecule("GLC") % "14bb"
#     mol *= 4

#     rot = Rotatron(mol.get_atom_graph(), n_processes=4)
#     rot._generate_rotation_unit_masks()
#     rot._find_rotation_units()
