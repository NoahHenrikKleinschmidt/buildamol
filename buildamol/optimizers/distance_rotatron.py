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

import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.graphs.base_graph as base_graph
import buildamol.utils.auxiliary as aux
import buildamol.structural.base as structural

# Rotatron = Rotatron.Rotatron


def concatenation_wrapper(x):
    pass


def simple_concatenation_function(x, unfold, pushback, n_smallest, clash_distance):
    """
    A simple concatentation function that computes the evaluation as:

    Mean distance ** unfold + (mean of n smallest distances) ** pushback
    """
    k = min(n_smallest, len(x) - 1)
    smallest = np.partition(x, k)[:k]  # np.sort(x)[:n_smallest]
    e = np.power(np.mean(x), unfold) + np.power(np.mean(smallest), pushback)
    return e


_numba_wrapper_simple_concatenation_function = aux.njit(simple_concatenation_function)


def concatenation_function_with_penalty(
    x, unfold, pushback, n_smallest, clash_distance
):
    """
    A concatentation function that computes the evaluation as:

    [(Mean distance ** unfold + (mean of n smallest distances) ** pushback)] / clash penalty
    """
    k = min(n_smallest, len(x) - 1)
    smallest = np.partition(x, k)[:k]  # np.sort(x)[:n_smallest]
    penalty = np.sum(x < 1.5 * clash_distance)
    e = np.power(np.mean(x), unfold) + np.power(np.mean(smallest), pushback)
    e /= (1 + penalty) ** 2
    return e


_numba_wrapper_concatenation_function_with_penalty = aux.njit(
    concatenation_function_with_penalty
)


def concatenation_function_no_pushback(x, unfold, pushback, n_smallest, clash_distance):
    """
    A concatentation function that computes the evaluation as:

    Mean distance ** unfold
    """
    e = np.power(np.mean(x), unfold)
    return e


_numba_wrapper_concatenation_function_no_pushback = aux.njit(
    concatenation_function_no_pushback
)


def concatenation_function_no_unfold(x, unfold, pushback, n_smallest, clash_distance):
    """
    A concatentation function that computes the evaluation as:

    (Mean of n smallest distances) ** pushback
    """
    k = min(n_smallest, len(x) - 1)
    smallest = np.partition(x, k)[:k]  # np.sort(x)[:n_smallest]
    e = np.power(np.mean(smallest), pushback)
    return e


_numba_wrapper_concatenation_function_no_unfold = aux.njit(
    concatenation_function_no_unfold
)


def concatenation_function_linear(x, unfold, pushback, n_smallest, clash_distance):
    """
    A concatentation function that computes the evaluation as:

    Mean distance * unfold + (mean of n smallest distances) * pushback
    """
    k = min(n_smallest, len(x) - 1)
    smallest = np.partition(x, k)[:k]  # np.sort(x)[:n_smallest]
    e = np.multiply(np.mean(x), unfold) + np.multiply(np.mean(smallest), pushback)
    return e


_numba_wrapper_concatenation_function_linear = aux.njit(concatenation_function_linear)


__numba_wrappers__ = {
    simple_concatenation_function: _numba_wrapper_simple_concatenation_function,
    concatenation_function_with_penalty: _numba_wrapper_concatenation_function_with_penalty,
    concatenation_function_no_pushback: _numba_wrapper_concatenation_function_no_pushback,
    concatenation_function_no_unfold: _numba_wrapper_concatenation_function_no_unfold,
    concatenation_function_linear: _numba_wrapper_concatenation_function_linear,
}


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
        This function should take the state array as first argument and may take any additional arguments.
        These additional arguments must be passed as keyword arguments to the environment during setup.
        The function must return a float.
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    n_processes : int
        The number of processes to use for parallel computation during edge mask generation.
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        radius: float = 20,
        pushback: float = 3,
        unfold: float = 2,
        clash_distance: float = 1.2,
        crop_nodes_further_than: float = -1,
        n_smallest: int = 10,
        concatenation_function: callable = None,
        bounds: tuple = (-np.pi, np.pi),
        n_processes: int = 1,
        **kwargs,
    ):
        self.hyperparameters = {
            "pushback": pushback,
            "unfold": unfold,
            "clash_distance": clash_distance,
            "crop_nodes_further_than": crop_nodes_further_than,
            "n_smallest": n_smallest,
            "concatenation_function": concatenation_function,
            "bounds": bounds,
            "radius": radius,
            "n_processes": n_processes,
            **kwargs,
        }
        self.kwargs = kwargs
        self.radius = radius
        self.crop_radius = crop_nodes_further_than * radius if radius > 0 else -1
        self.clash_distance = clash_distance
        self.pushback = pushback
        self.unfold = unfold
        self.n_smallest = n_smallest

        # =====================================

        if self.crop_radius > 0:
            graph, rotatable_edges = self._setup_helpers_crop_faraway_nodes(
                self.crop_radius, graph, rotatable_edges
            )
        # =====================================

        if radius > 0:
            self._radius = radius
        else:
            self._radius = np.inf

        # =====================================

        # self.concatenation_function = concatenation_wrapper
        self._state_dists = np.zeros((len(graph.nodes), len(graph.nodes)))

        # =====================================

        self.edx = 0

        Rotatron.Rotatron.__init__(
            self, graph, rotatable_edges, n_processes=n_processes, **kwargs
        )
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )
        # =====================================

        self._last_eval = 0.0
        self.n_smallest = min(n_smallest, self.n_nodes)
        self.hyperparameters["n_smallest"] = self.n_smallest

        # =====================================

        self._numba_func_signature = (
            "unfold",
            "pushback",
            "n_smallest",
            "clash_distance",
        )

        if concatenation_function is None:
            concatenation_function = concatenation_function_with_penalty

        if (
            kwargs.get("numba", False)
            or aux.USE_ALL_NUMBA
            or (self.n_nodes**2 > 100000 and aux.USE_NUMBA)
        ):
            concatenation_function = __numba_wrappers__.get(
                concatenation_function, concatenation_function
            )

            self.eval = self._numba_eval
        else:
            self.eval = self._normal_eval

        self._concatenation_function = concatenation_function
        self._concatenation_function_kwargs = aux.get_args(
            self._concatenation_function, self.hyperparameters
        )
        self._bounds_tuple = bounds

        # =====================================

    def _normal_eval(self, state):
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

        mask = pairwise_dists < self._radius
        mask = np.logical_and(mask, self.rotation_unit_masks)
        dist_eval = np.zeros(len(pairwise_dists))

        _changed_entries = mask.any(axis=1)
        dist_eval[~_changed_entries] = -1

        for i in np.where(_changed_entries)[0]:
            dist_eval[i] = self._concatenation_function(
                pairwise_dists[i][mask[i]], **self._concatenation_function_kwargs
            )

        mean_dist_eval = np.divide(1.0, np.mean(dist_eval[dist_eval > -1]))

        final = np.log(mean_dist_eval)  # - self._backup_eval

        min_dist = np.min(pairwise_dists)
        self._state_dists = min_dist
        self._last_eval = final
        return final

    def _numba_eval(self, state):
        min_dist, final = _numba_wrapper_eval(
            state=state,
            concatenation_function=self._concatenation_function,
            rotation_unit_masks=self.rotation_unit_masks,
            last_eval=self._last_eval,
            radius=self._radius,
            **self._concatenation_function_kwargs,
        )
        self._state_dists = min_dist
        self._last_eval = final
        return final

    def is_done(self, state):
        return np.min(self._state_dists) > self.clash_distance

    def concatenation_function(self, x):
        mask = x < self._radius
        mask = np.logical_and(mask, self.rotation_unit_masks[self.ndx])
        self.ndx += 1
        if not np.logical_or.reduce(mask):
            return -1
        return self._concatenation_function(
            x[mask], **self._concatenation_function_kwargs
        )


@aux.njit
def _numba_wrapper_eval(
    state,
    concatenation_function,
    rotation_unit_masks,
    last_eval,
    radius,
    unfold,
    pushback,
    n_smallest,
    clash_distance,
):
    pairwise_dists = structural._numba_wrapper_euclidean_distances(state, state)
    np.fill_diagonal(pairwise_dists, radius)

    dist_eval = np.zeros(len(state))

    mask = pairwise_dists < radius
    mask = np.logical_and(mask, rotation_unit_masks)
    dist_eval = np.zeros(len(pairwise_dists))

    for i in range(len(pairwise_dists)):
        if not mask[i].any():
            dist_eval[i] = -1
        else:
            dist_eval[i] = concatenation_function(
                pairwise_dists[i][mask[i]], unfold, pushback, n_smallest, clash_distance
            )

    mask = dist_eval > -1
    mean_dist_eval = np.divide(1.0, np.mean(dist_eval[mask]))

    final = np.log(mean_dist_eval)  # - self._backup_eval

    min_dist = np.min(pairwise_dists)
    return min_dist, final


__all__ = [
    "DistanceRotatron",
    "simple_concatenation_function",
    "concatenation_function_with_penalty",
    "concatenation_function_no_pushback",
    "concatenation_function_no_unfold",
    "concatenation_function_linear",
]

if __name__ == "__main__":

    import buildamol as bam

    DistanceRotatron._backup_eval = 0.0

    mol = bam.Molecule.from_json(
        "/Users/noahhk/GIT/biobuild/buildamol/optimizers/__testing__/files/EX8.json"
    )
    print("init: ", mol.count_clashes())
    graph = mol.get_residue_graph(True)
    env = DistanceRotatron(
        graph, numba=False, concatenation_function=simple_concatenation_function
    )  # , concatenatiofn_function=simple_concatenation_function)

    from time import time

    for i in range(15):

        t0 = time()
        out = bam.optimizers.optimize(mol.copy(), env)
        print(time() - t0, out.count_clashes())

#     import matplotlib.pyplot as plt
#     import seaborn as sns

#     mol = bam.molecule("/Users/noahhk/GIT/biobuild/_tutorials copy/ext8_opt.pdb")

#     graph = mol.get_residue_graph(True)
#     edges = graph.find_rotatable_edges(min_descendants=3)

#     from time import time

#     for i in range(5, 10):
#         t0 = time()
#         d = DistanceRotatron(
#             graph,
#             edges,
#             pushback=3,
#             concatenation_function=concatenation_function_with_penalty,
#             n_processes=i,
#         )
#         # opt = bam.optimizers.optimize(mol, d, "swarm")
#         print(time() - t0)
#     exit()

#     t1 = time()
#     d = DistanceRotatron(
#         graph,
#         edges,
#         pushback=3,
#         concatenation_function=concatenation_function_with_penalty,
#         n_processes=5,
#     )
#     opt = bam.optimizers.optimize(mol, d, "swarm")

#     print(opt.count_clashes())
#     opt.to_pdb(f"opt9_pow_pushback_{d.pushback}.pdb")

# import stable_baselines3 as sb3

# model = sb3.PPO("MlpPolicy", d, verbose=1)
# model.learn(total_timesteps=10000)
# model.save("ppo_distance_rotatron")

# x_ = d._best_eval
# x0 = d.step(d.blank())
# a = d.blank()
# a[22] = 0.1
# x = d.step(a)
# print(x[1])
# pass

# -------------------------------------

# bam.load_sugars()
# glc = bam.molecule("GLC")
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
