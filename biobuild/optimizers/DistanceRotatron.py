import gym

import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

import biobuild.optimizers.Rotatron as Rotatron
import biobuild.structural as structural
import biobuild.utils.auxiliary as aux
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
    clash_distance : float
        The distance at which atoms are considered to be clashing.
    pushback : float
        The pushback factor for the distance calculation.
        The higher this value, the shorter the retraction effect.
        Tweak this parameter only if you know what you're doing and adjusting the radius alone is not enough.
    mask_rotation_units : bool
        If True, atoms that are part of the same rotational unit (i.e. between two rotational edges) are
        masked off from each other when computing the evaluation. This prevents atoms from interferring with
        the clash detection. This has less of an impact for larger radii but is curcial for small radii!
    crop_nodes_further_than: float
        Nodes that are further away than this factor times the radius from any rotatable edge at the beginning
        of the optimization are removed from the graph and not considered during optimization. This speeds up
        computation. Set to -1 to disable.
    concatenation_function : callable
        The function to use when computing the evaluation for each node.
        This function should take a 1D array and return a scalar.
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    """

    # clash_penalty : float
    #     The penalty for a clash between atoms.

    def __init__(
        self,
        graph: "BaseGraph.BaseGraph",
        rotatable_edges: list = None,
        radius: float = -1,
        clash_distance: float = 0.9,
        pushback: float = 3,
        # clash_penalty: float = 3.0,
        mask_rotation_units: bool = True,
        crop_nodes_further_than: float = -1,
        concatenation_function: callable = np.mean,
        bounds: tuple = (-np.pi, np.pi),
    ):
        self.radius = radius
        self.crop_radius = crop_nodes_further_than * radius if radius > 0 else -1
        # self.clash_penalty = clash_penalty
        self.clash_distance = clash_distance
        self.mask_rotation_units = mask_rotation_units
        self._concatenation_function = concatenation_function

        self._bounds_tuple = bounds

        self._current_edge_idx = 0
        self._current_node_index = 0

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

        self._pushback = pushback

        # =====================================

        if radius > 0:
            self._radius = radius
        else:
            self._radius = np.inf

        # =====================================

        if not self.mask_rotation_units:

            def masker(x):
                return x < self._radius

        else:

            def masker(x):
                mask = (x < self._radius).astype(np.int8)
                if self.edge_masks[self._current_edge_idx, self._current_node_index]:
                    mask *= (1 - self.edge_masks[self._current_edge_idx]).astype(
                        np.int8
                    )
                return mask.astype(bool)

        def concatenation_wrapper(x):
            mask = masker(x)
            self._current_node_index = (self._current_edge_idx + 1) % self.n_edges
            if not np.logical_or.reduce(mask):
                return -1
            return self._concatenation_function(x[mask])

        self.concatenation_function = concatenation_wrapper

        # =====================================

        # if self.clash_penalty > 0:

        #     def clash_handler():
        #         clashes = self._state_dists < self.clash_distance
        #         self._state_dists[clashes] = (
        #             self._state_dists[clashes] ** self.clash_penalty
        #         )

        # else:

        #     def clash_handler():
        #         pass

        # self.clash_handler = clash_handler

        # =====================================

        self._state_dists = np.zeros((self.state.shape[0], self.state.shape[0]))
        self._last_eval = np.inf
        self._last_eval = self.eval(self.state)
        self._best_eval = self._last_eval
        self._best_clashes = self.count_clashes()

    def eval(self, state):  # , diff=True):
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

        # if diff:
        #     mask = np.abs(pairwise_dists - self._state_dists) < 1e-4
        #     mask *= pairwise_dists > self.clash_distance
        #     pairwise_dists[mask] = self._radius + 1

        # self.clash_handler()

        rowwise_dist_eval = np.apply_along_axis(
            self.concatenation_function, 1, pairwise_dists
        )
        mask = rowwise_dist_eval > -1
        if not np.logical_or.reduce(mask):
            return self._last_eval

        rowwise_dist_eval[mask] **= self._pushback
        rowwise_dist_eval[mask] += 1e-6
        rowwise_dist_eval[mask] /= self.n_nodes
        rowwise_dist_eval[mask] **= -1
        rowwise_dist_eval[mask] = np.log(rowwise_dist_eval[mask])

        final = np.sum(rowwise_dist_eval[mask])
        self._state_dists[:, :] = pairwise_dists
        self._last_eval = final
        return final

    def step(self, action):
        for i, edge in enumerate(self.rotatable_edges):
            self._current_edge_idx = i
            new_state = self.rotate(
                edge,
                action[i],
            )

        done = self.is_done()
        self._last_eval = self.eval(new_state)
        self._action_history += action
        clashes = self.count_clashes()

        if self._last_eval < self._best_eval and clashes <= self._best_clashes:
            self._best_eval = self._last_eval
            self._best_state *= 0
            self._best_state += new_state
            # self._best_action *= 0
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
