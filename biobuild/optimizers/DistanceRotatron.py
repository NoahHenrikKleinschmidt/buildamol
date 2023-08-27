"""
The DistanceRotatron environment evaulates conformations based on the pairwise distances between nodes in the optimized graph.

It uses two forces, a global "unfolding" force to maximize spacial separation between nodes, and a local "pushback" force to maximize distances between the closest nodes.

The evaluation is computed as:

.. math::

    e_i = \\sum_{j \\neq i} d_{ij}^{unfold} + pushback \\cdot \\sum_{k=1}^N \\text{sorted}(d)_{ik}

There are multiple variations of this basic formulation available (see the functions below). 
"""
import gym

import numpy as np
from scipy.spatial.distance import cdist

import biobuild.optimizers.Rotatron as Rotatron
import biobuild.graphs.BaseGraph as BaseGraph

# Rotatron = Rotatron.Rotatron


def simple_concatenation_function(self, x):
    """
    A simple concatentation function that computes the evaluation as:

    Mean distance ** unfold + (mean of n smallest distances) ** pushback
    """
    smallest = np.sort(x)[: self.n_smallest]
    e = np.power(np.mean(x), self.unfold) + np.power(np.mean(smallest), self.pushback)
    return e


def concatenation_function_with_penalty(self, x):
    """
    A concatentation function that computes the evaluation as:

    (Mean distance ** unfold + (mean of n smallest distances) ** pushback) / clash penalty
    """
    smallest = np.sort(x)[: self.n_smallest]
    penalty = np.sum(x < 1.5 * self.clash_distance)
    e = np.power(np.mean(x), self.unfold) + np.power(np.mean(smallest), self.pushback)
    e /= (1 + penalty) ** 2
    return e


def concatenation_function_no_pushback(self, x):
    """
    A concatentation function that computes the evaluation as:

    Mean distance ** unfold
    """
    e = np.power(np.mean(x), self.unfold)
    return e


def concatenation_function_no_unfold(self, x):
    """
    A concatentation function that computes the evaluation as:

    Mean distance + pushback * mean of n smallest distances
    """
    smallest = np.sort(x)[: self.n_smallest]
    e = np.power(np.mean(smallest), self.pushback)
    return e


def concatenation_function_linear(self, x):
    """
    A concatentation function that computes the evaluation as:

    Mean distance * unfold + (mean of n smallest distances) * pushback
    """
    smallest = np.sort(x)[: self.n_smallest]
    e = np.multiply(np.mean(x), self.unfold) + np.multiply(
        np.mean(smallest), self.pushback
    )
    return e


class DistanceRotatron(Rotatron):
    """
    A distance-based Rotatron environment.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    radius : float
        The radius around rotatable edges to include in the distance calculation.
        Set to -1 to disable.
    pushback : float
        Short distances between atoms are given higher weight in the evaluation using this factor.
    unfold : float
        The exponent to use when computing the mean distance to others for
        each node. Higher values give higher values to global unfolding of the graph.
    clash_distance : float
        The distance at which atoms are considered to be clashing.
    crop_nodes_further_than : float
        Nodes that are further away than this factor times the radius from any rotatable edge at the beginning
        of the optimization are removed from the graph and not considered during optimization. This speeds up
        computation. Set to -1 to disable.
    n_smallest : int
        The number of smallest distances to use when computing the evaluation for each node.
    concatenation_function : callable
        A custom function to use when computing the evaluation for each node.
        This function should take the environment (self) as first argument and a 1D array of pairwise-distances from one node to all others as second argument and return a scalar.
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    """

    def __init__(
        self,
        graph: "BaseGraph.BaseGraph",
        rotatable_edges: list = None,
        radius: float = 20,
        pushback: float = 3,
        unfold: float = 2,
        clash_distance: float = 0.9,
        crop_nodes_further_than: float = -1,
        n_smallest: int = 10,
        concatenation_function: callable = None,
        bounds: tuple = (-np.pi, np.pi),
        **kwargs,
    ):
        self.kwargs = kwargs
        self.radius = radius
        self.crop_radius = crop_nodes_further_than * radius if radius > 0 else -1
        self.clash_distance = clash_distance
        self.pushback = pushback
        self.unfold = unfold
        self.n_smallest = n_smallest

        if concatenation_function is None:
            concatenation_function = concatenation_function_with_penalty

        self._concatenation_function = concatenation_function
        self._bounds_tuple = bounds

        # =====================================
        if self.crop_radius > 0:
            rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)
            edge_coords = np.array(
                [(a.coord + b.coord) / 2 for a, b in rotatable_edges]
            )
            nodes = list(graph.nodes)
            node_coords = np.array([node.coord for node in nodes])

            dists = cdist(edge_coords, node_coords)
            dists = dists > self.crop_radius
            dists = np.apply_along_axis(np.all, 0, dists)
            if np.max(dists) != 0:
                nodes_to_drop = [nodes[i] for i, d in enumerate(dists) if d]
                graph.remove_nodes_from(nodes_to_drop)
        # =====================================

        if radius > 0:
            self._radius = radius
        else:
            self._radius = np.inf

        # =====================================

        def concatenation_wrapper(x):
            mask = x < self._radius
            mask = np.logical_and(mask, self.rotation_unit_masks[self.ndx])
            self.ndx += 1
            if not np.logical_or.reduce(mask):
                return -1
            return self._concatenation_function(self, x[mask])

        self.concatenation_function = concatenation_wrapper
        self._state_dists = np.zeros((len(graph.nodes), len(graph.nodes)))

        # =====================================

        self.edx = 0

        Rotatron.__init__(self, graph, rotatable_edges)
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )
        # =====================================

        self._best_clashes = self.count_clashes()

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
        pairwise_dists = cdist(state, state)
        np.fill_diagonal(pairwise_dists, self._radius)

        self.ndx = 0
        dist_eval = np.apply_along_axis(self.concatenation_function, 1, pairwise_dists)
        mask = dist_eval > -1

        if not np.logical_or.reduce(mask):
            return self._last_eval

        mean_dist_eval = np.divide(1.0, np.mean(dist_eval[mask]))

        final = np.log(mean_dist_eval)
        self._state_dists[:, :] = pairwise_dists
        self._last_eval = final
        return final

    def step(self, action):
        for i, edge in enumerate(self.rotatable_edges):
            self.edx = i
            new_state = self._rotate(i, action[i])

        self._last_eval = self.eval(new_state)
        clashes = self.count_clashes()
        done = clashes == 0
        self._action_history += action

        if self._last_eval < self._best_eval and clashes <= self._best_clashes:
            self._best_eval = self._last_eval
            self._best_state *= 0
            self._best_state += new_state
            self._best_action *= 0
            self._best_action += action
            self._best_clashes = clashes

        return new_state, self._last_eval, done, {}

    def is_done(self, state):
        return self.count_clashes() == 0

    def count_clashes(self):
        return np.sum(self._state_dists < self.clash_distance)


__all__ = [
    "DistanceRotatron",
    "simple_concatenation_function",
    "concatenation_function_with_penalty",
    "concatenation_function_no_pushback",
    "concatenation_function_no_unfold",
    "concatenation_function_linear",
]


if __name__ == "__main__":
    import biobuild as bb
    import matplotlib.pyplot as plt
    import seaborn as sns

    mol = bb.molecule("/Users/noahhk/GIT/biobuild/_tutorials copy/ext8_opt.pdb")

    graph = mol.get_residue_graph(True)
    edges = graph.find_rotatable_edges(min_descendants=3)

    d = DistanceRotatron(
        graph,
        edges,
        pushback=3,
        concatenation_function=concatenation_function_with_penalty,
    )
    opt = bb.optimizers.optimize(mol, d, "swarm")

    print(opt.count_clashes())
    opt.to_pdb(f"opt9_pow_pushback_{d.pushback}.pdb")

    # import stable_baselines3 as sb3

    # model = sb3.PPO("MlpPolicy", d, verbose=1)
    # model.learn(total_timesteps=10000)
    # model.save("ppo_distance_rotatron")

    x_ = d._best_eval
    x0 = d.step(d.blank())
    a = d.blank()
    a[22] = 0.1
    x = d.step(a)
    print(x[1])
    pass

    # -------------------------------------

    # bb.load_sugars()
    # glc = bb.molecule("GLC")
    # glc.repeat(2, "14bb")

    # bonds = [glc.get_bonds("O4", "C4")[0]]
    # env = DistanceRotatron(glc.make_atom_graph(), bonds, radius=20)
    # actions = np.arange(
    #     -np.pi,
    #     np.pi,
    #     np.pi / 80,
    # )
    # cmap = sns.color_palette("Blues", len(actions))

    # evals = []
    # v = glc.draw()
    # for i in actions:
    #     new_state, e, done, _ = env.step(np.array([i]))
    #     evals.append([i, e])

    #     _glc = glc.copy()
    #     _glc.rotate_around_bond(*bonds[0], np.degrees(i), descendants_only=True)

    #     color = "lightgray"
    #     if e > 0.2157:
    #         color = "red"
    #     # elif e > 0.1970:
    #     #     color = "orange"
    #     # elif e < 0.1965:
    #     #     color = "green"
    #     if color == "lightgray":
    #         opacity = 0.4
    #     else:
    #         opacity = 1.0
    #     v.draw_edges(
    #         _glc.get_bonds(_glc.residues[1]), color=color, linewidth=3, opacity=opacity
    #     )

    # evals = np.array(evals)

    # plt.plot(evals[:, 0], evals[:, 1])
    # v.show()
    # plt.show()
    # pass
