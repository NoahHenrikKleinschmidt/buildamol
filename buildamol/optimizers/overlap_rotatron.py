"""
The OverlapRotatron is a environment that approximates molecular graphs using multi-variat Gaussian distributions. The overlap between the distributions is used as the evaluation function for the environment.
Hence, this environment tries to minimize the overlap between distributions in order to find favorable conformations.

As measure for the overlap between two distributions, the Jensen-Shannon divergence is used by default. Custom overlap functions can be passed to the environment.
"""

import gym

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal

# from sklearn.mixture import GaussianMixture
# from scipy.stats import entropy


import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.graphs.base_graph as base_graph

__all__ = [
    "OverlapRotatron",
    # "likelihood_overlap",
    # "bhattacharyya_overlap",
    "jensen_shannon_overlap",
    "MVN",
]


def MVN(points, spread: float = 1.0):
    """
    Compute a multi-variate normal distribution for a given set of points.

    Parameters
    ----------
    points : np.ndarray
        The points to compute the mean and covariance matrix for.

    Returns
    -------
    mvn : scipy.stats.multivariate_normal
        The multi-variate normal distribution for the points.
    """
    return multivariate_normal(
        mean=np.mean(points, axis=0),
        cov=spread * np.cov(points, rowvar=False),
        allow_singular=True,
    )


# def likelihood_overlap(mvn1, mvn2):
#     """
#     Compute the overlap between two gaussians using likelihoods.

#     Parameters
#     ----------
#     mvn1, mvn2 : scipy.stats.multivariate_normal
#         The two gaussians to compute the overlap for.
#     Returns
#     -------
#     overlap : float
#         The overlap between the two gaussians.
#     """
#     center1 = mvn1.mean
#     center2 = mvn2.mean

#     # Compute the overlap between the two distributions
#     # using the likelihoods of the centers of the distributions
#     dist1 = mvn1.pdf(center2)
#     dist2 = mvn2.pdf(center1)

#     if dist1 == 0 or dist2 == 0:
#         return 0

#     overlap = np.log(dist1) - np.log(dist2)
#     return overlap


# def bhattacharyya_overlap(mvn1, mvn2):
#     """
#     Compute the overlap between two gaussians using the Bhattacharyya coefficient.

#     Parameters
#     ----------
#     mvn1, mvn2 : scipy.stats.multivariate_normal
#         The two gaussians to compute the overlap for.

#     Returns
#     -------
#     overlap : float
#         The overlap between the two gaussians.
#     """

#     # Create Multivariate Normal distributions for each Gaussian
#     center1 = mvn1.mean
#     center2 = mvn2.mean

#     dist1 = mvn1.pdf(center2)
#     dist2 = mvn2.pdf(center1)

#     dist = dist1 * dist2
#     if dist == 0:
#         return dist

#     # Compute the Bhattacharyya coefficient (overlap between the distributions)
#     bhattacharyya_coefficient = np.sqrt(dist)

#     # Compute the Bhattacharyya distance
#     bhattacharyya_distance = np.log(bhattacharyya_coefficient)
#     return bhattacharyya_distance


def jensen_shannon_overlap(mvn1, mvn2):
    """
    Compute the overlap between two gaussians using the Jensen-Shannon divergence.

    Parameters
    ----------
    mvn1, mvn2 : scipy.stats.multivariate_normal
        The two gaussians to compute the overlap for.

    Returns
    -------
    overlap : float
        The overlap between the two gaussians.
    """

    # Create Multivariate Normal distributions for each Gaussian
    center1 = mvn1.mean
    center2 = mvn2.mean

    pdf1_center1 = mvn1.pdf(center1)
    pdf2_center2 = mvn2.pdf(center2)
    pdf1_center2 = mvn1.pdf(center2)
    pdf2_center1 = mvn2.pdf(center1)

    mean_center1 = (pdf1_center1 + pdf2_center1) / 2
    mean_center2 = (pdf1_center2 + pdf2_center2) / 2

    dist = _kl_divergence(pdf1_center1, mean_center1) + _kl_divergence(
        pdf2_center2, mean_center2
    )
    dist *= 0.5

    return -dist


def _kl_divergence(p, q):
    """
    Compute the Kullback-Leibler divergence between two distributions.
    """
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))


# Rotatron = Rotatron.Rotatron


class OverlapRotatron(Rotatron.Rotatron):
    """
    A distribution overlap-based Rotatron environment.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    artificial_spread : float
        The spread to use for the multi-variate normal distributions. This is used to artificially increase the spread of the distributions.
        This is useful for cases where the distributions may be too tight and far apart which makes it difficult for the overlap to be computed.
    clash_distance : float
        The distance at which two atoms are considered to be clashing.
    crop_nodes_further_than : float
        If greater than 0, crop nodes that are further than this distance from the
        rotatable edges so that they are not considered in the overlap calculation.
    distance_function : callable
        A specific distance function to use for calculating the overlap. This function
        should take two arrays of shape (1, 3) (centers) and two arrays of shape (3, 3) (covariances)
        and return a scalar.
    ignore_further_than : float
        If greater than 0, centroids that are further than this distance from each other are evaluated as 0 overlap automatically.
    n_processes : int
        The number of parallel processes to use when computing edge masks.
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        artificial_spread: float = 2.0,
        clash_distance: float = 1.2,
        crop_nodes_further_than: float = -1,
        distance_function: callable = None,
        ignore_further_than: float = -1,
        n_processes: int = 1,
        bounds: tuple = (-np.pi, np.pi),
        **kwargs,
    ):
        self.hyperparameters = {
            "artificial_spread": artificial_spread,
            "clash_distance": clash_distance,
            "crop_nodes_further_than": crop_nodes_further_than,
            "distance_function": distance_function,
            "ignore_further_than": ignore_further_than,
            "n_processes": n_processes,
            "bounds": bounds,
            **kwargs,
        }
        self.crop_radius = crop_nodes_further_than
        self.clash_distance = clash_distance
        self.ignore_further_than = ignore_further_than > 0
        self._ignore_distance = ignore_further_than
        self.distance_function = distance_function or jensen_shannon_overlap

        self._bounds_tuple = bounds

        # =====================================

        rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)
        self.graph = graph
        self.rotatable_edges = rotatable_edges
        self.n_nodes = len(self.graph.nodes)
        self.n_edges = len(self.rotatable_edges)
        self.artificial_spread = artificial_spread

        # =====================================
        if self.crop_radius > 0:
            graph, rotatable_edges = self._setup_helpers_crop_faraway_nodes(
                self.crop_radius, graph, rotatable_edges
            )
        # =====================================

        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )
        Rotatron.Rotatron.__init__(
            self, graph, rotatable_edges, n_processes=n_processes, **kwargs
        )

        # =====================================

        self._generate_rotation_unit_masks()
        self._find_rotation_units()
        self.rotation_units = {
            k: v for k, v in self.rotation_units.items() if len(v) > 1
        }

        # =====================================

        # this is the mainloop for computing pairwise overlaps
        n = 0
        for i, gmm1 in enumerate(self.rotation_units):
            for j, gmm2 in enumerate(self.rotation_units):
                if i >= j:
                    continue
                n += 1
        self.overlaps = np.zeros(n + 1)
        self.centers = np.zeros((n + 1, 3))
        self.covariances = np.zeros((n + 1, 3, 3))

        # =====================================

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
        gaussians = []
        for i, mask in self.rotation_units.items():
            gaussians.append(MVN(state[mask], spread=self.artificial_spread))

        idx = 0
        for i, G1 in enumerate(gaussians):
            for j, G2 in enumerate(gaussians):
                if i >= j:
                    continue

                if (
                    self.ignore_further_than
                    and np.linalg.norm(G1.mean - G2.mean) > self._ignore_distance
                ):
                    self.overlaps[idx] = 0
                else:
                    self.overlaps[idx] = self.distance_function(G1, G2)
                idx += 1

        return np.mean(self.overlaps)


if __name__ == "__main__":
    import buildamol as bam
    from time import time

    mol = bam.molecule(
        "/Users/noahhk/GIT/biobuild/buildamol/optimizers/_testing/files/GLYCAN.json"
    )

    graph = mol.get_residue_graph()
    graph.make_detailed(n_samples=0.5, include_far_away=True)

    edges = graph.find_rotatable_edges(root_node=graph.central_node, max_descendants=20)
    env = OverlapRotatron(graph, edges, artificial_spread=2)

    print(len(edges))

    better = bam.optimizers.optimize(mol.copy(), env)
    print(mol.count_clashes(), better.count_clashes())
    better.show()

    exit()
    edges = graph.find_rotatable_edges(mol.get_atom(168), min_descendants=10)

    # env = OverlapRotatron(
    #     graph,
    #     edges,
    #     distance_function=jensen_shannon_overlap,
    # )

    # out = bam.optimizers.optimize(mol.copy(), env, "genetic", max_generations=30)
    # out.show()

    # # v = graph.draw()
    # # v.draw_edges(*edges, color="magenta", linewidth=3, opacity=1.0)
    # # v.show()
    from alive_progress import alive_bar

    n = 10
    clashes = np.zeros((3, n))
    times = np.zeros((3, n))
    with alive_bar(n * 3) as bar:
        for i, func in enumerate(
            [jensen_shannon_overlap],  # likelihood_overlap, bhattacharyya_overlap
        ):
            env = OverlapRotatron(
                graph,
                edges,
                distance_function=func,
            )
            for j in range(n):
                t1 = time()
                out = bam.optimizers.optimize(
                    mol.copy(),
                    env,
                    "swarm",
                )
                t = time() - t1
                times[i, j] = t
                clashes[i, j] = out.count_clashes()
                bar()

    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    df = pd.DataFrame(
        clashes.T,
        columns=["likelihood", "bhattacharyya", "jensen_shannon"],
    )
    df = df.melt(var_name="overlap", value_name="clashes")
    sns.violinplot(data=df, x="overlap", y="clashes", ax=axs[0])

    df2 = pd.DataFrame(
        times.T,
        columns=["likelihood", "bhattacharyya", "jensen_shannon"],
    )
    df2 = df2.melt(var_name="overlap", value_name="time")
    sns.violinplot(data=df2, x="overlap", y="time", ax=axs[1])

    axs[0].set_ylabel("Clashes")
    axs[1].set_ylabel("Time (s)")

    sns.despine()

    plt.savefig("overlap_rotatron_dist_func_comparisons_EX6.png")
    plt.show()

    # v = out.draw()
    # v.draw_edges(*mol.bonds, color="lightblue", opacity=0.5)
    # v.show()

    # if False:
    #     # FOR MAKING THE COOL FIGURES
    #     graph = mol.get_atom_graph()

    #     # edges = graph.find_rotatable_edges(min_descendants=10)
    #     edges = [mol.get_bond(150, 152)]

    #     d = OverlapRotatron(graph, edges, bounds=(0, 0.5))

    #     angles = np.arange(-np.pi, np.pi, np.pi / 10)
    #     evals = []
    #     for i in angles:
    #         new_state, e, done, _ = d.step([i])
    #         evals.append(e)
    #         d.reset()

    #     evals = np.array(evals)
    #     import matplotlib.pyplot as plt
    #     import seaborn as sns

    #     rgba_to_hex = lambda x: "#%02x%02x%02x" % tuple([int(i * 255) for i in x])

    #     cmap = sns.color_palette("coolwarm", len(evals))

    #     # evals /= np.min(evals)

    #     v = mol.draw()
    #     v.draw_edges(edges[0], color="magenta", linewidth=3, opacity=1.0)
    #     for i, e in enumerate(evals):
    #         s = mol.copy()
    #         s.rotate_around_bond(
    #             *edges[0], angles[i], descendants_only=True, angle_is_degrees=False
    #         )
    #         v.draw_edges(*s.bonds, color=rgba_to_hex(cmap[i]), linewidth=3, opacity=0.6)

    #     v.show()

    #     plt.plot(angles, evals)

    #     best_angle = angles[np.argmin(evals)]
    #     s = mol.copy()
    #     s.rotate_around_bond(
    #         *edges[0], best_angle, descendants_only=True, angle_is_degrees=False
    #     )
    #     v.draw_edges(*s.bonds, color="limegreen", linewidth=3, opacity=1.0)

    #     pass
