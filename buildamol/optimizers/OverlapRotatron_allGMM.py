import gym

import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal

from sklearn.mixture import GaussianMixture


import buildamol.optimizers.Rotatron as Rotatron
import buildamol.graphs.BaseGraph as BaseGraph

# __all__ = [
#     "OverlapRotatron",
#     "likelihood_overlap",
#     "bhattacharyya_overlap",
#     "jensen_shannon_overlap",
#     "gaussian",
# ]


# Rotatron = Rotatron.Rotatron


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
        ignore_further_than: float = -1,
        bounds: tuple = (-np.pi, np.pi),
    ):
        self.crop_radius = crop_nodes_further_than
        self.clash_distance = clash_distance
        self.ignore_further_than = ignore_further_than > 0
        self._ignore_distance = ignore_further_than

        self._bounds_tuple = bounds
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

        Rotatron.__init__(self, graph, rotatable_edges)
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )

        # =====================================

        self.n_components = self.n_nodes // 3
        self._gmm = GaussianMixture(n_components=self.n_components).fit(self.state)

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

        self._gmm.fit(state)
        average_volume = np.mean(
            [
                np.sqrt(np.linalg.det(cov_matrix))
                for cov_matrix in self._gmm.covariances_
            ]
        )
        return -average_volume


if __name__ == "__main__":
    import buildamol as bam
    from time import time

    mol = bam.molecule("/Users/noahhk/GIT/biobuild/docs/_tutorials/ext_final_opt.pdb")
    mol.autolabel()

    graph = mol.get_residue_graph()
    graph.make_detailed(n_samples=0.6)

    edges = graph.find_rotatable_edges(graph.central_node, min_descendants=20)

    env = OverlapRotatron(graph, edges)
    t1 = time()
    out = bam.optimizers.optimize(mol.copy(), env, "scipy")
    print(out.count_clashes())
    print(time() - t1)
    out.show()
    pass
    # v = graph.draw()
    # v.draw_edges(*edges, color="magenta", linewidth=3, opacity=1.0)
    # v.show()
    # from alive_progress import alive_bar

    # n = 5
    # clashes = np.zeros((3, n))
    # times = np.zeros((3, n))
    # with alive_bar(n * 3) as bar:
    #     for i, func in enumerate(
    #         [likelihood_overlap, bhattacharyya_overlap, jensen_shannon_overlap]
    #     ):
    #         env = OverlapRotatron(
    #             graph, edges, distance_function=func, ignore_further_than=6
    #         )
    #         for j in range(n):
    #             t1 = time()
    #             out = bam.optimizers.optimize(
    #                 mol.copy(), env, "genetic", max_generations=100
    #             )
    #             t = time() - t1
    #             times[i, j] = t
    #             clashes[i, j] = out.count_clashes()
    #             bar()

    # import matplotlib.pyplot as plt
    # import pandas as pd
    # import seaborn as sns

    # fig, axs = plt.subplots(1, 2, figsize=(10, 5))

    # df = pd.DataFrame(
    #     clashes.T,
    #     columns=["likelihood", "bhattacharyya", "jensen_shannon"],
    # )
    # df = df.melt(var_name="overlap", value_name="clashes")
    # sns.violinplot(data=df, x="overlap", y="clashes", ax=axs[0])

    # df2 = pd.DataFrame(
    #     times.T,
    #     columns=["likelihood", "bhattacharyya", "jensen_shannon"],
    # )
    # df2 = df2.melt(var_name="overlap", value_name="time")
    # sns.violinplot(data=df2, x="overlap", y="time", ax=axs[1])

    # axs[0].set_ylabel("Clashes")
    # axs[1].set_ylabel("Time (s)")

    # sns.despine()

    # plt.show()

    # v = out.draw()
    # v.draw_edges(*mol.bonds, color="lightblue", opacity=0.5)
    # v.show()

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
