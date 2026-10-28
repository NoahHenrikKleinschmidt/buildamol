"""
The ForceFieldRotatron is a rotatron that uses RDKit's MMFF94 force field to evaluate a given state. Consequently, this environment can only function if RDKIt is installed.

.. note::

    Because this environment uses an actual energy function to evaluate states, it always requires all atoms of the molecule.
    If a ResidueGraph or a graph that only contains a subset of the atoms is provided, the full AtomGraph of the molecule is used instead.

"""

import warnings

import gymnasium as gym

import numpy as np


import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.graphs.base_graph as base_graph
import buildamol.utils.auxiliary as aux

# Rotatron = Rotatron.Rotatron

__all__ = ["ForceFieldRotatron"]


class ForceFieldRotatron(Rotatron.Rotatron):
    """
    A force field based rotatron. This rotatron uses RDKit's MMFF94 force field to
    evaluate the energy of a given state.

    Parameters
    ----------
    graph : AtomGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    clash_distance : float
        The distance at which two atoms are considered to be clashing.
    crop_nodes_further_than : float
        Not supported by this environment (the force field always requires all atoms). Kept for compatibility.
    mmff_variant : str
        The MMFF variant to use. Can be one of "mmff94" or "mmff94s"
    n_processes : int
        The number of processes to use for parallelization when computing edge masks
    bounds : tuple
        The bounds for the minimal and maximal rotation angles.
    kwargs
        Additional keyword arguments to pass to the Rotatron
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        clash_distance: float = 0.9,
        crop_nodes_further_than: float = -1,
        mmff_variant: str = "mmff94",
        n_processes: int = 1,
        bounds: tuple = (-np.pi, np.pi),
        **kwargs,
    ):
        self.hyperparameters = {
            "clash_distance": clash_distance,
            "crop_nodes_further_than": crop_nodes_further_than,
            "mmff_variant": mmff_variant,
            "n_processes": n_processes,
            "bounds": bounds,
            **kwargs,
        }
        self.crop_radius = crop_nodes_further_than
        self.clash_distance = clash_distance
        self._bounds_tuple = bounds
        self.mmff_variant = mmff_variant
        # =====================================

        rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)
        molecule = graph._molecule

        # the energy can only be computed if all atoms move with their rotatable edges
        # (any atom that is not part of the graph would stay in place and break its bonds)
        if set(graph.nodes) != set(molecule.get_atoms()):
            warnings.warn(
                "The ForceFieldRotatron requires a graph with all atoms of the molecule. Using the molecule's AtomGraph instead."
            )
            graph = molecule.get_atom_graph()
        if self.crop_radius > 0:
            warnings.warn(
                "The ForceFieldRotatron requires all atoms of the molecule. 'crop_nodes_further_than' is ignored."
            )

        self.graph = graph
        self.rotatable_edges = rotatable_edges
        self.mol = molecule.to_rdkit()

        # =====================================

        Rotatron.Rotatron.__init__(
            self, graph, rotatable_edges, n_processes=n_processes, **kwargs
        )
        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )

        self._ff_props = aux.MMFFGetMoleculeProperties(
            self.mol, mmffVariant=self.mmff_variant
        )
        self._ff = aux.MMFFGetMoleculeForceField(self.mol, self._ff_props)

        # map explicitily to atoms since we cannot be sure that the graph nodes correspond 1:1
        atom_index = {a: i for i, a in enumerate(molecule.get_atoms())}
        self._state_rows = np.array(
            [i for i, n in enumerate(self.graph.nodes) if n in atom_index], dtype=int
        )
        self._rdkit_rows = np.array(
            [atom_index[n] for n in self.graph.nodes if n in atom_index], dtype=int
        )
        self._positions = np.array(self.mol.GetConformer(0).GetPositions())

        # =====================================

    def _update_positions(self, state):
        """
        Update the rdkit molecule's positions
        """
        self._positions[self._rdkit_rows] = state[self._state_rows]
        conformer = self.mol.GetConformer(0)
        for idx, coord in enumerate(self._positions):
            conformer.SetAtomPosition(idx, coord)

    def energy(self, state: np.ndarray = None):
        """
        Calculate the energy of a given state

        Parameters
        ----------
        state : np.ndarray
            The state of the environment. If None, the energy of the
            current positions of the rdkit molecule is computed.

        Returns
        -------
        float
            The energy for the state
        """
        if state is None:
            self._positions[:] = self.mol.GetConformer(0).GetPositions()
        else:
            self._positions[self._rdkit_rows] = state[self._state_rows]
        # the positions are passed explicitly because the force field does not
        # reliably pick up changes of the conformer positions
        return self._ff.CalcEnergy(self._positions.ravel().tolist())

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
        return self.energy(state)


if __name__ == "__main__":
    import buildamol as bam
    from time import time

    mol = bam.molecule(
        "/Users/noahhk/GIT/biobuild/buildamol/optimizers/__testing__/files/EX8.json"
    )
    mol.autolabel()

    graph = mol.get_atom_graph()

    edges = graph.find_rotatable_edges(
        graph.central_node, min_descendants=20
    )  # mol.get_atom(168), min_descendants=10)

    env = ForceFieldRotatron(graph, edges)
    out = bam.optimizers.optimize(mol.copy(), env, "swarm")
    out.show()

    # v = graph.draw()
    # v.draw_edges(*edges, color="magenta", linewidth=3, opacity=1.0)
    # v.show()
    # from alive_progress import alive_bar

    # n = 30
    # clashes = np.zeros((3, n))
    # times = np.zeros((3, n))
    # with alive_bar(n * 3) as bar:
    #     for i, func in enumerate(
    #         [likelihood_overlap, bhattacharyya_overlap, jensen_shannon_overlap]
    #     ):
    #         env = OverlapRotatron(graph, edges, distance_function=func)
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
