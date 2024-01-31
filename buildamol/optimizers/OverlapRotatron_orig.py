import gym

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal


import buildamol.optimizers.Rotatron as Rotatron
import buildamol.graphs.BaseGraph as BaseGraph

__all__ = [
    "OverlapRotatron",
    "likelihood_overlap",
    "bhattacharyya_overlap",
    "jensen_shannon_overlap",
    "gaussian",
]


def gaussian(points):
    """
    Compute the mean and covariance matrix of a given set of points.

    Parameters
    ----------
    points : np.ndarray
        The points to compute the mean and covariance matrix for.

    Returns
    -------
    center : np.ndarray
        The mean of the points.
    covariance : np.ndarray
        The covariance matrix of the points.
    """
    # Compute the center (mean) of the points cloud
    center = np.mean(points, axis=0)

    # Compute the covariance matrix of the points cloud
    covariance = np.cov(points, rowvar=False)

    return center, covariance


def likelihood_overlap(center1, center2, cov1, cov2):
    """
    Compute the overlap between two gaussians using likelihoods.

    Parameters
    ----------
    center1 : np.ndarray
        The center of the first gaussian.
    center2 : np.ndarray
        The center of the second gaussian.
    cov1 : np.ndarray
        The covariance matrix of the first gaussian.
    cov2 : np.ndarray
        The covariance matrix of the second gaussian.

    Returns
    -------
    overlap : float
        The overlap between the two gaussians.
    """

    # Create Multivariate Normal distributions for each Gaussian
    dist1 = multivariate_normal(mean=center1, cov=cov1)
    dist2 = multivariate_normal(mean=center2, cov=cov2)

    # Compute the overlap between the two distributions
    # using the likelihoods of the centers of the distributions
    dist1 = dist1.pdf(center2)
    dist2 = dist2.pdf(center1)

    if dist1 == 0 or dist2 == 0:
        return 0

    overlap = np.log(dist1) - np.log(dist2)
    return overlap


def bhattacharyya_overlap(center1, center2, cov1, cov2):
    """
    Compute the overlap between two gaussians using the Bhattacharyya coefficient.

    Parameters
    ----------
    center1 : np.ndarray
        The center of the first gaussian.
    center2 : np.ndarray
        The center of the second gaussian.
    cov1 : np.ndarray
        The covariance matrix of the first gaussian.
    cov2 : np.ndarray
        The covariance matrix of the second gaussian.

    Returns
    -------
    overlap : float
        The overlap between the two gaussians.
    """

    # Create Multivariate Normal distributions for each Gaussian
    dist1 = multivariate_normal(mean=center1, cov=cov1)
    dist2 = multivariate_normal(mean=center2, cov=cov2)

    dist1 = dist1.pdf(center2)
    dist2 = dist2.pdf(center1)

    dist = dist1 * dist2
    if dist == 0:
        return dist

    # Compute the Bhattacharyya coefficient (overlap between the distributions)
    bhattacharyya_coefficient = np.sqrt(dist)

    # Compute the Bhattacharyya distance
    bhattacharyya_distance = np.log(bhattacharyya_coefficient)
    return bhattacharyya_distance


def jensen_shannon_overlap(center1, center2, cov1, cov2):
    """
    Compute the overlap between two gaussians using the Jensen-Shannon divergence.

    Parameters
    ----------
    center1 : np.ndarray
        The center of the first gaussian.
    center2 : np.ndarray
        The center of the second gaussian.
    cov1 : np.ndarray
        The covariance matrix of the first gaussian.
    cov2 : np.ndarray
        The covariance matrix of the second gaussian.

    Returns
    -------
    overlap : float
        The overlap between the two gaussians.
    """

    # Create Multivariate Normal distributions for each Gaussian
    dist1 = multivariate_normal(mean=center1, cov=cov1)
    dist2 = multivariate_normal(mean=center2, cov=cov2)

    pdf1_center1 = dist1.pdf(center1)
    pdf2_center2 = dist2.pdf(center2)
    pdf1_center2 = dist1.pdf(center2)
    pdf2_center1 = dist2.pdf(center1)

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


Rotatron = Rotatron.Rotatron


class OverlapRotatron(Rotatron):
    """
    A distribution overlap-based Rotatron environment.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
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
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    """

    def __init__(
        self,
        graph: "BaseGraph.BaseGraph",
        rotatable_edges: list = None,
        clash_distance: float = 0.9,
        crop_nodes_further_than: float = -1,
        distance_function: callable = None,
        ignore_further_than: float = -1,
        bounds: tuple = (-np.pi, np.pi),
    ):
        self.crop_radius = crop_nodes_further_than
        self.clash_distance = clash_distance
        self.ignore_further_than = ignore_further_than > 0
        self._ignore_distance = ignore_further_than

        self._bounds_tuple = bounds
        # =====================================
        if not distance_function:
            distance_function = bhattacharyya_overlap

        self._distance_function = distance_function
        # =====================================

        rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)
        self.graph = graph
        self.rotatable_edges = rotatable_edges

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

        self._residue_graph = self.graph.__class__.__name__ == "ResidueGraph"

        residues = list(
            set(i.parent for i in self.graph.nodes if i.__class__.__name__ == "Atom")
        )

        self._residue_masks = np.array(
            [
                np.array(
                    [
                        i.parent == r if i.__class__.__name__ == "Atom" else False
                        for i in self.graph.nodes
                    ]
                )
                for r in residues
            ]
        )

        # =====================================

        # this is the mainloop for computing pairwise overlaps
        n = 0
        for i, mask in enumerate(self._residue_masks[:-1]):
            for j, mask2 in enumerate(self._residue_masks[i + 1 :], start=i + 1):
                n += 1

        self.overlaps = np.zeros(n + 1)
        self.centers = {}
        self.covariances = {}

        # =====================================

        Rotatron.__init__(self, graph, rotatable_edges)
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )

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
        # for each residue compute the gaussian

        # for i, mask in enumerate(self._residue_masks):
        #     points = state[mask]
        #     center, covariance = gaussian(points)
        #     self._gaussian_centers[i] = center
        #     self._gaussian_covariances[i] = covariance

        # for each residue in the graph, pairwise compute the overlap between the gaussians
        # and return the sum of the overlaps
        self.overlaps *= 0
        self.centers.clear()
        self.covariances.clear()
        idx = 0
        for i, mask in enumerate(self._residue_masks[:-1]):
            if i not in self.centers:
                self.centers[i], self.covariances[i] = gaussian(state[mask])
            for j, mask2 in enumerate(self._residue_masks[i + 1 :], start=i + 1):
                if j not in self.centers:
                    self.centers[j], self.covariances[j] = gaussian(state[mask2])

                if (
                    self.ignore_further_than
                    and np.linalg.norm(self.centers[i] - self.centers[j])
                    > self._ignore_distance
                ):
                    overlap = 0
                else:
                    overlap = self._distance_function(
                        self.centers[i],
                        self.centers[j],
                        self.covariances[i],
                        self.covariances[j],
                    )
                    if not np.isfinite(overlap):
                        overlap = 0
                self.overlaps[idx] = overlap
                idx += 1

        return np.mean(self.overlaps)


if __name__ == "__main__":
    import buildamol as bam
    from time import time

    mol = bam.molecule(
        "/Users/noahhk/GIT/biobuild/biobuild/optimizers/_testing/files/EX6.json"
    )

    graph = mol.get_residue_graph()
    graph.make_detailed(n_samples=0.8)

    edges = graph.find_rotatable_edges(mol.get_atom(168), min_descendants=10)

    # v = graph.draw()
    # v.draw_edges(*edges, color="magenta", linewidth=3, opacity=1.0)
    # v.show()
    from alive_progress import alive_bar

    n = 5
    clashes = np.zeros((3, n))
    times = np.zeros((3, n))
    with alive_bar(n * 3) as bar:
        for i, func in enumerate(
            [likelihood_overlap, bhattacharyya_overlap, jensen_shannon_overlap]
        ):
            env = OverlapRotatron(
                graph, edges, distance_function=func, ignore_further_than=6
            )
            for j in range(n):
                t1 = time()
                out = bam.optimizers.optimize(
                    mol.copy(), env, "genetic", max_generations=100
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

    plt.show()

    v = out.draw()
    v.draw_edges(*mol.bonds, color="lightblue", opacity=0.5)
    v.show()

    if False:
        # FOR MAKING THE COOL FIGURES
        graph = mol.get_atom_graph()

        # edges = graph.find_rotatable_edges(min_descendants=10)
        edges = [mol.get_bond(150, 152)]

        d = OverlapRotatron(graph, edges, bounds=(0, 0.5))

        angles = np.arange(-np.pi, np.pi, np.pi / 10)
        evals = []
        for i in angles:
            new_state, e, done, _ = d.step([i])
            evals.append(e)
            d.reset()

        evals = np.array(evals)
        import matplotlib.pyplot as plt
        import seaborn as sns

        rgba_to_hex = lambda x: "#%02x%02x%02x" % tuple([int(i * 255) for i in x])

        cmap = sns.color_palette("coolwarm", len(evals))

        # evals /= np.min(evals)

        v = mol.draw()
        v.draw_edges(edges[0], color="magenta", linewidth=3, opacity=1.0)
        for i, e in enumerate(evals):
            s = mol.copy()
            s.rotate_around_bond(
                *edges[0], angles[i], descendants_only=True, angle_is_degrees=False
            )
            v.draw_edges(*s.bonds, color=rgba_to_hex(cmap[i]), linewidth=3, opacity=0.6)

        v.show()

        plt.plot(angles, evals)

        best_angle = angles[np.argmin(evals)]
        s = mol.copy()
        s.rotate_around_bond(
            *edges[0], best_angle, descendants_only=True, angle_is_degrees=False
        )
        v.draw_edges(*s.bonds, color="limegreen", linewidth=3, opacity=1.0)

        pass
