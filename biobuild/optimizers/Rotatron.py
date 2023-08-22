"""
This is the basic Rotatron environment. It provides the basic functionality for preprocessing a graph into numpy arrays, masking rotatable edges, and evaluating a possible solution.
All other Rotatron environments inherit from this class.
"""
import gym
import numpy as np

from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

import biobuild.graphs.BaseGraph as BaseGraph


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
    """

    def __init__(
        self,
        graph: "BaseGraph.BaseGraph",
        rotatable_edges: list = None,
    ):
        self.graph = graph
        self.rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)

        self.node_dict = {n: i for i, n in enumerate(self.graph.nodes)}
        self.n_nodes = len(self.node_dict)
        self.n_edges = len(self.rotatable_edges)

        self.state = self._make_state_from_graph(self.graph)
        self._backup_state = self.state.copy()

        self.rotation_unit_masks = np.ones(
            (len(graph.nodes), len(graph.nodes)), dtype=bool
        )
        self._generate_edge_masks()
        self._generate_edge_lengths()

        self._edge_node_coords = np.array(
            [
                [self.get_node_idx(e[0]), self.get_node_idx(e[1])]
                for e in self.rotatable_edges
            ]
        )

        self._best_state = self.state.copy()
        self._best_action = self.blank()
        self._action_history = self.blank()
        self._last_eval = self._init_eval(self.state)
        self._best_eval = self._last_eval
        self._backup_eval = self._last_eval

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

    def _generate_edge_masks(self):
        """
        Compute the edge masks of downstream nodes
        """
        self.edge_masks = np.array(
            [
                [
                    1 if i in self.graph.get_descendants(*e) else 0
                    for i in self.graph.nodes
                ]
                for e in self.rotatable_edges
            ],
            dtype=bool,
        )

    @property
    def best(self):
        """
        The best state, the action that lead there, and evaluation that the environment has seen
        """
        return self._best_state, self._best_action, self._best_eval

    def get_edge_idx(self, _edge):
        return self.rotatable_edges.index(_edge)

    def get_node_idx(self, _node):
        return self.node_dict[_node]

    def get_edge_vector(self, _edge):
        adx, bdx = self.get_node_idx(_edge[0]), self.get_node_idx(_edge[1])
        vec = self.state[bdx] - self.state[adx]
        return vec

    def _get_edge_vector(self, edx):
        adx, bdx = self._edge_node_coords[edx]
        vec = self.state[bdx] - self.state[adx]
        return vec

    def _get_edge_ref_coord(self, edx):
        adx, bdx = self._edge_node_coords[edx]
        return self.state[adx]

    def get_node_coords(self, _node):
        return self.state[self.get_node_idx(_node)]

    def rotate(self, edge, angle, edx=None):
        """
        Rotate the graph around an edge

        Parameters
        ----------
        edge : tuple
            The edge to rotate around
        angle : float
            The angle to rotate by

        Returns
        -------
        np.ndarray
            The new state of the environment
        """
        if -1e-3 < angle < 1e-3:
            return self.state
        edx = edx or self.get_edge_idx(edge)
        mask = self.edge_masks[edx]
        vec = self.get_edge_vector(edge)
        # new version where the lengths are pre-computed
        # since we are only rotating the lengths should not change...
        length = self.edge_lengths[edx]
        vec /= length
        ref_coord = self.get_node_coords(edge[0])

        rot = Rotation.from_rotvec(vec * angle)
        self.state[mask] = rot.apply(self.state[mask] - ref_coord) + ref_coord
        return self.state

    def _rotate(self, edx, angle):
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
        if -1e-3 < angle < 1e-3:
            return self.state

        mask = self.edge_masks[edx]
        vec = self._get_edge_vector(edx)
        # new version where the lengths are pre-computed
        # since we are only rotating the lengths should not change...
        length = self.edge_lengths[edx]
        vec /= length
        ref_coord = self._get_edge_ref_coord(edx)

        rot = Rotation.from_rotvec(vec * angle)
        self.state[mask] = rot.apply(self.state[mask] - ref_coord) + ref_coord
        return self.state

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

    def _init_eval(self, state):
        """
        The evaluation score that is computed before the first step
        """
        return self.eval(state)

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

        for i, edge in enumerate(self.rotatable_edges):
            new_state = self.rotate(
                edge,
                action[i],
            )
        e = self.eval(new_state)
        done = self.is_done(new_state)
        self._action_history += action

        if e < self._best_eval:
            self._best_eval = e
            self._best_action = self._action_history.copy()
            self._best_state = new_state
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

    def reset(self, state: bool = True, best: bool = False):
        """
        Reset the environment
        """
        if state:
            self.state[:, :] = self._backup_state
        if best:
            self._best_state[:, :] = self._backup_state
            self._best_action[:] = 0
            self._action_history[:] = 0
            self._best_eval = self._backup_eval

    def blank(self):
        """
        A blank action
        """
        return np.zeros(len(self.rotatable_edges))

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
            state2 = self._rotate(i, angle)
        dists2 = cdist(state2, state2)
        for i, angle in enumerate(np.random.random(self.n_edges)):
            state3 = self._rotate(i, angle)
        dists3 = cdist(state3, state3)

        d12 = np.abs(dists1 - dists2) < 1e-4
        d13 = np.abs(dists1 - dists3) < 1e-4
        d23 = np.abs(dists2 - dists3) < 1e-4

        dists = np.sum([d12, d13, d23], axis=0) == 3
        self.rotation_unit_masks = ~dists
        self.reset()

    def _node_rotation_unit_mask(self, edx, ndx):
        """
        Get the rotation unit mask for a node

        Parameters
        ----------
        edx : int
            The edge index
        ndx : int
            The node index

        Returns
        -------
        np.ndarray
            The rotation unit mask
        """
        array = self.edge_masks[edx]
        if array[ndx]:
            return ~array
        return array

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
        pass


if __name__ == "__main__":
    import biobuild as bb

    bb.load_sugars()
    mol = bb.molecule("GLC") % "14bb"
    mol *= 4

    rot = Rotatron(mol.get_atom_graph())
    rot._generate_rotation_unit_masks()
    rot._find_rotation_units()
