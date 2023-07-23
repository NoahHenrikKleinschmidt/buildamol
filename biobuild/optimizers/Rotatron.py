import gym
import numpy as np

from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

import biobuild.structural as structural
import biobuild.utils.auxiliary as aux
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

        self._best_state = self.state.copy()
        self._best_action = self.blank()
        self._action_history = self.blank()
        self._best_eval = np.inf

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
        return vec / np.linalg.norm(vec)

    def get_node_coords(self, _node):
        return self.state[self.get_node_idx(_node)]

    def rotate(self, edge, angle):
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
        edx = self.get_edge_idx(edge)
        mask = self.edge_masks[edx]
        vec = self.get_edge_vector(edge)
        ref_coord = self.get_node_coords(edge[0])

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
        per_edge_evals = []
        for i, edge in enumerate(self.rotatable_edges):
            new_state = self.rotate(
                edge,
                action[i],
            )
            per_edge_evals.append(self.eval(new_state))
            pass

        e = per_edge_evals[-1]
        done = self.is_done(new_state)
        self._action_history += action

        if e < self._best_eval:
            self._best_eval = e
            self._best_action = self._action_history.copy()
            self._best_state = new_state
        return new_state, e, done, {"per_edge_evals": per_edge_evals}

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
            self.state[:, :] = 0
            self.state += self._backup_state
        if best:
            self._best_state[:, :] = 0
            self._best_state += self._backup_state
            self._best_action[:] = 0
            self._action_history[:] = 0
            self._best_eval = self.eval(self.state)

        # for idx, node in enumerate(self.graph.nodes):
        #     node.coord = self.state[idx]

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
                and not (e[0] in _circulars and e[1] in _circulars)
                and len(graph.get_descendants(*e)) > 1
            ]
        return rotatable_edges
