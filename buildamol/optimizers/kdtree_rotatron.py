"""
The KDTreeRotatron environment evaluates conformations based on steric clashes between atoms.

Only atom pairs that are closer than a clash distance are considered. These are found using KD-trees,
which makes the evaluation scale roughly linearly with the number of atoms (as opposed to the quadratic scaling of
environments that consider all pairwise distances). The atoms that do not move during the optimization are stored
in a KD-tree only once, so large static environments (e.g. a protein scaffold) are cheap.

The evaluation is a soft clash penalty:

.. math::

    e = \\sum_{i<j, d_{ij} < d_{clash}} (d_{clash} - d_{ij})^2

where only pairs of atoms are considered that can actually change their distance during the optimization
(i.e. atoms that belong to different rotation units and that are not separated by less than three bonds).
Consequently, a conformation without clashes has an evaluation of exactly zero.
"""

import gymnasium as gym
import numpy as np
from scipy.spatial import cKDTree
import scipy.sparse as sparse

import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.graphs.base_graph as base_graph
import buildamol.base_classes as base_classes

__all__ = ["KDTreeRotatron"]


class KDTreeRotatron(Rotatron.Rotatron):
    """
    A clash-based Rotatron environment using KD-trees.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize. Clashes are always evaluated between atoms, so if a ResidueGraph is given,
        the AtomGraph of the underlying molecule is used instead (the rotatable edges are kept).
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    clash_distance : float
        Atoms that are closer than this distance (in Angstrom) are penalized.
    include_hydrogens : bool
        If True, hydrogen atoms are also considered for clashes. By default only heavy atoms are considered.
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    n_processes : int
        The number of processes to use for parallel computation during edge mask generation.
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        clash_distance: float = 2.8,
        include_hydrogens: bool = False,
        bounds: tuple = (-np.pi, np.pi),
        n_processes: int = 1,
        **kwargs,
    ):
        self.hyperparameters = {
            "clash_distance": clash_distance,
            "include_hydrogens": include_hydrogens,
            "bounds": bounds,
            "n_processes": n_processes,
            **kwargs,
        }
        self.clash_distance = clash_distance
        self.include_hydrogens = include_hydrogens
        self._bounds_tuple = bounds

        # clashes are only meaningful between atoms
        if not all(isinstance(n, base_classes.Atom) for n in graph.nodes):
            if getattr(graph, "_molecule", None) is None:
                raise ValueError(
                    "The KDTreeRotatron requires an AtomGraph or a graph that is linked to a molecule."
                )
            rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)
            graph = graph._molecule.get_atom_graph()

        Rotatron.Rotatron.__init__(
            self, graph, rotatable_edges, n_processes=n_processes, **kwargs
        )
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )
        self._last_eval = 0.0
        self._setup_clash_evaluation()

    def _setup_clash_evaluation(self):
        """
        Precompute which atoms are considered, which of them move, the KD-tree of the static atoms,
        and the pairs of atoms whose distance cannot change.
        """
        nodes = list(self.graph.nodes)
        self._rows = np.array(
            [
                i
                for i, n in enumerate(nodes)
                if self.include_hydrogens or getattr(n, "element", None) != "H"
            ],
            dtype=int,
        )
        masks = self.edge_masks[:, self._rows]

        # atoms with the same movement pattern belong to the same rotation unit (fixed distances)
        if len(self._rows):
            _, units = np.unique(masks.T, axis=0, return_inverse=True)
            self._units = np.asarray(units).ravel()
        else:
            self._units = np.zeros(0, dtype=int)

        moving = masks.any(axis=0)
        self._moving = np.where(moving)[0]
        self._static = np.where(~moving)[0]
        coords = self.state[self._rows]
        self._static_tree = cKDTree(coords[self._static]) if len(self._static) else None

        # atoms that are one or two bonds apart have a fixed distance
        n = len(nodes)
        index = {node: i for i, node in enumerate(nodes)}
        edges = np.array(
            [(index[a], index[b]) for a, b in self.graph.edges], dtype=int
        ).reshape(-1, 2)
        adjacency = sparse.coo_matrix(
            (np.ones(2 * len(edges)), (edges.ravel(), edges[:, ::-1].ravel())),
            shape=(n, n),
        ).tocsr()
        near = (adjacency + adjacency @ adjacency)[self._rows][:, self._rows].tocoo()
        self._near_keys = np.unique(
            np.minimum(near.row, near.col) * len(self._rows)
            + np.maximum(near.row, near.col)
        )

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
            The evaluation for the state (0 if there are no clashes)
        """
        self._last_eval = self._clash(state[self._rows])
        return self._last_eval

    def eval_batch(self, actions):
        """
        Evaluate a batch of actions at once (without changing the state of the environment).

        Parameters
        ----------
        actions : np.ndarray
            (n_candidates, n_edges) rotation angles

        Returns
        -------
        np.ndarray
            The evaluations of the actions
        np.ndarray
            Whether each action leads to a conformation without clashes
        """
        states = self._rotate_batch(actions)
        evals = np.array([self._clash(state[self._rows]) for state in states])
        return evals, evals <= 0

    def is_done(self, state):
        return self._last_eval <= 0

    def _clash(self, coords):
        """
        The soft clash penalty for the coordinates of the considered atoms
        """
        if len(self._moving) == 0:
            return 0.0
        moving = coords[self._moving]
        tree = cKDTree(moving)

        pairs = tree.query_pairs(self.clash_distance, output_type="ndarray")
        i = self._moving[pairs[:, 0]]
        j = self._moving[pairs[:, 1]]
        d = np.linalg.norm(coords[i] - coords[j], axis=1)

        if self._static_tree is not None:
            static_pairs = tree.sparse_distance_matrix(
                self._static_tree, self.clash_distance, output_type="ndarray"
            )
            i = np.concatenate([i, self._moving[static_pairs["i"]]])
            j = np.concatenate([j, self._static[static_pairs["j"]]])
            d = np.concatenate([d, static_pairs["v"]])

        keys = np.minimum(i, j) * len(self._rows) + np.maximum(i, j)
        variable = (self._units[i] != self._units[j]) & ~np.isin(
            keys, self._near_keys
        )
        return float(np.sum(np.square(self.clash_distance - d[variable])))
