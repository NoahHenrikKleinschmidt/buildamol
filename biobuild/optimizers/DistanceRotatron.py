import gym

import numpy as np
from scipy.spatial.distance import cdist

import biobuild.optimizers.Rotatron as Rotatron
import biobuild.graphs.BaseGraph as BaseGraph


class DistanceRotatron(Rotatron.Rotatron):
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
        pushback: float = 2,
        clash_distance: float = 0.9,
        crop_nodes_further_than: float = -1,
        n_smallest: int = 5,
        concatenation_function: callable = None,
        bounds: tuple = (-np.pi, np.pi),
    ):
        self.radius = radius
        self.crop_radius = crop_nodes_further_than * radius if radius > 0 else -1
        self.clash_distance = clash_distance
        self.pushback = pushback
        self.n_smallest = n_smallest

        if concatenation_function is None:

            def concatenation_function(self, x):
                smallest = np.sum(np.sort(x)[: self.n_smallest])
                return np.mean(x) + self.pushback * smallest

        self._concatenation_function = concatenation_function

        self._bounds_tuple = bounds

        # =====================================

        rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)

        # =====================================
        if self.crop_radius > 0:
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

        Rotatron.Rotatron.__init__(self, graph, rotatable_edges)
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )

        # =====================================

        if radius > 0:
            self._radius = radius
        else:
            self._radius = np.inf

        # =====================================

        def concatenation_wrapper(x):
            mask = x < self._radius
            if not np.logical_or.reduce(mask):
                return -1
            return self._concatenation_function(self, x[mask])

        self.concatenation_function = concatenation_wrapper

        # =====================================

        self._state_dists = np.zeros((self.state.shape[0], self.state.shape[0]))
        self._last_eval = np.inf
        self._last_eval = self.eval(self.state)
        self._best_eval = self._last_eval
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

        dist_eval = np.apply_along_axis(self.concatenation_function, 1, pairwise_dists)
        mask = dist_eval > -1

        if not np.logical_or.reduce(mask):
            return self._last_eval

        mean_dist_eval = 1.0 / np.mean(dist_eval[mask])

        final = np.log(mean_dist_eval)
        self._state_dists[:, :] = pairwise_dists
        self._last_eval = final
        return final

    def step(self, action):
        for i, edge in enumerate(self.rotatable_edges):
            new_state = self.rotate(
                edge,
                action[i],
            )

        self._last_eval = self.eval(new_state)
        clashes = self.count_clashes()
        done = clashes == 0
        self._action_history += action

        if self._last_eval < self._best_eval and clashes <= self._best_clashes:
            self._best_eval = self._last_eval
            self._best_state *= 0
            self._best_state += new_state
            self._best_action += self._action_history
            self._best_clashes = clashes

        return new_state, self._last_eval, done, {}

    def is_done(self):
        return self.count_clashes() == 0

    def count_clashes(self):
        return np.sum(self._state_dists < self.clash_distance)


if __name__ == "__main__":
    import biobuild as bb
    import matplotlib.pyplot as plt
    import seaborn as sns

    mol = bb.molecule(
        "/Users/noahhk/GIT/biobuild/biobuild/optimizers/_testing/files/EX6.json"
    )

    graph = mol.get_atom_graph()
    edges = graph.find_rotatable_edges(min_descendants=10)

    d = DistanceRotatron(graph, edges, radius=-1, bounds=(0, 0.5))
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
