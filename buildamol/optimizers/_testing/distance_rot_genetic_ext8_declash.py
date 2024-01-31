# %% [markdown]
# This notebook template is designed for testing the performance of Rotatron environments and different solving agents of different scales.

"""
NOTE
This version of the script is only for small, medium, and large whole-molecule tests, not for scaffold tests...
"""

# %%
# =============================================================================
# Work on local biobuild in GIT repo
# =============================================================================
from collections import defaultdict
import os, sys, importlib

# for inside python scripts
base = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
# base = os.path.dirname(os.path.dirname(os.getcwd()))
sys.path.insert(0, base)


def reload_optimizers():
    importlib.reload(bam.optimizers.environments)
    importlib.reload(bam.optimizers.agents)


# =============================================================================
import files
import auxiliary
import buildamol as bam
import buildamol.optimizers as optimizers

import time
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist


# %% [markdown]
# Here we can select which tests to run and on what testing structures:

# %%
# which structures to run on
structures_to_run_on = [
    # files.GLUCOSE2,
    # files.PEPTIDE,
    # files.X_WING,
    "/Users/noahhk/GIT/biobuild/docs/_tutorials/ext8_opt.pdb"
    # files.SCAFFOLD1,
    # files.SCAFFOLD3,
]

# how many times to independently run on each structure
re_runs = 20

# visualize evaluation history
visualize_eval_history = True

# visualize time history
visualize_time_history = True

# visualise clashes in final structure
visualize_clashes = True

# clash threshold
clash_cutoff = 0.9

# visualize the final structure
visualize_final_structure = False

# visualization parameters
# for draw_edges()
visualization_params = dict(color="magenta", opacity=0.3)

# export visualizations
export_visualizations = True

# export history to csv
export_history = True


# graph building function
# provide a custom callable that generates a tuple of (graph, rotatable_edges
def make_graph(x):
    g = x.make_residue_graph()
    g.make_detailed(n_samples=0.7)
    edges = x.get_residue_connections()
    edges = g.direct_edges(edges=edges)
    return g, edges


graph_factory = make_graph

# graph building parameters
graph_params = {}

# provide a custom callable to set a custom building function for the environment
rotatron_factory = None

# the rotatron class to use
rotatron_class = optimizers.DistanceRotatron


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

    dist_eval = np.apply_along_axis(self.concatenation_function, 1, pairwise_dists)
    mask = dist_eval > -1
    # mean_dist_eval = 1.0 / np.mean(dist_eval[mask])
    # min_dist_eval = 1.0 / (np.min(dist_eval[mask] ** self._pushback) + 1e-6)

    final = np.log(
        np.mean(dist_eval[mask]) / (np.min(pairwise_dists) ** self._pushback)
    )
    self._state_dists[:, :] = pairwise_dists
    self._last_eval = final
    return final

    # rowwise_dist_eval = np.apply_along_axis(
    #     self.concatenation_function, 1, pairwise_dists
    # )
    # mask = rowwise_dist_eval > -1
    # if not np.logical_or.reduce(mask):
    #     return self._last_eval

    rowwise_dist_eval[mask] **= self._pushback
    rowwise_dist_eval[mask] += 1e-6
    rowwise_dist_eval[mask] /= self.n_nodes
    rowwise_dist_eval[mask] **= -1
    rowwise_dist_eval[mask] = np.log(rowwise_dist_eval[mask])

    final = np.sum(rowwise_dist_eval[mask])
    self._state_dists[:, :] = pairwise_dists
    self._last_eval = final
    return final


optimizers.DistanceRotatron.eval = eval


p = 4


def concat1(x):
    return 1 / np.mean(x)


# rotatron parameters
rotatron_params = {
    "radius": 20,
    "pushback": p,
    "concatenation_function": concat1,
    # "concatenation_function": lambda x: np.mean(x) - np.min(x) ** 4,
}

# export name prefix
export_name_prefix = f"xwing_declash_genetic_new_concat_rad{int(rotatron_params['radius'])}_push{int(p*10)}_"

# the agent function to use
agent = optimizers.genetic_optimize

# agent parameters
agent_params = {"threshold": 1e-6, "stop_if_done": True}


# %% [markdown]
# Perform some environment setup

# %%
if agent is None:
    raise ValueError("No agent provided")
if rotatron_class is None:
    raise ValueError("No rotatron class provided")

if graph_factory is None:
    graph_factory = auxiliary.graph_factory
if rotatron_factory is None:
    rotatron_factory = auxiliary.rotatron_factory

structures_to_run_on = (bam.molecule(s) for s in structures_to_run_on)

eval_history = defaultdict(list)
time_history = defaultdict(list)
clash_history = defaultdict(list)

initial_evals = {}
initial_clashes = {}
final_visuals = {}

v = None


def make_environment(structure):
    """
    An environment generator
    """
    graph, rotatable_edges = graph_factory(structure, **graph_params)
    return rotatron_factory(rotatron_class, graph, rotatable_edges, **rotatron_params)


# %% [markdown]
# Now start the main testing code

# %%
for structure in structures_to_run_on:
    env = make_environment(structure)
    initial_evals[structure.id] = [env._best_eval] * re_runs
    initial_clashes[structure.id] = [
        auxiliary.count_clashes(env.graph, clash_cutoff)
    ] * re_runs

    if visualize_final_structure:
        if not v:
            v = structure.draw()
            v.draw_edges(*env.rotatable_edges, color="cyan", linewidth=6)

    for r in range(re_runs):
        t1 = time.time()
        # we are interested in learning the full time to make and solve the environment
        env = make_environment(structure)
        sol, eval = agent(env, **agent_params)
        t2 = time.time()
        eval_history[structure.id].append(eval)
        time_history[structure.id].append(t2 - t1)
        final = auxiliary.apply_solution(sol, env, structure.copy())
        final.to_pdb(f"{export_name_prefix}.{structure.id}.{r}.pdb")

        clash_history[structure.id].append(auxiliary.count_clashes(final, clash_cutoff))

        if visualize_final_structure:
            v.draw_edges(*final.bonds, **visualization_params)
        env.reset()
        print(
            f"Finished run {r+1} of {re_runs} on {structure.id} ({eval}, {t2-t1}, {clash_history[structure.id][-1]})"
        )

    if visualize_final_structure:
        _best = auxiliary.apply_solution(env.best[1], env, structure.copy())
        v.draw_edges(*_best.bonds, color="green", linewidth=6)
        final_visuals[structure.id] = v
        _best.show()
        v = None

# %% [markdown]
# And now do some data collecting and visualization

# %%
eval_history = auxiliary.transform_to_df(
    eval_history, initial_evals, "final", "initial"
)
clash_history = auxiliary.transform_to_df(
    clash_history, initial_clashes, "final", "initial"
)
time_history = auxiliary.transform_to_df(time_history)

if export_history:
    if not export_name_prefix:
        export_name_prefix = rotatron_class.__name__ + "." + agent.__name__
    eval_history.to_csv(f"{export_name_prefix}.eval_history.csv", index=False)
    time_history.to_csv(f"{export_name_prefix}.time_history.csv", index=False)
    clash_history.to_csv(f"{export_name_prefix}.clash_history.csv", index=False)

if visualize_eval_history or visualize_time_history or visualize_clashes:
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, axs = plt.subplots(1, 4, figsize=(15, 3))
    if visualize_eval_history:
        sns.barplot(data=eval_history, ax=axs[0], x="key", y="final")
        axs[0].set(title="Evaluation of finals", ylabel="eval-score", xlabel="")

    if visualize_time_history:
        sns.barplot(data=time_history, ax=axs[1], x="key", y=0)
        axs[1].set(title="Computation times", ylabel="seconds", xlabel="")

    if visualize_clashes:
        clash_history["diff"] = clash_history["final"] - clash_history["initial"]
        sns.barplot(data=clash_history, ax=axs[2], x="key", y="final")
        axs[2].set(title="Clashes in finals", ylabel="clashes in final", xlabel="")
        sns.barplot(data=clash_history, ax=axs[3], x="key", y="diff")
        axs[3].set(
            title="Clashes in finals",
            ylabel="diff in clashes (final - initial)",
            xlabel="",
        )

    fig.tight_layout()

    fig.suptitle(f"{rotatron_class.__name__} + {agent.__name__}")
    plt.savefig(f"{export_name_prefix}.plots.png")

# %% [markdown]
# Here can the 3d visualizations be viewed then
# ---

# %%
final_visuals

if visualize_final_structure:
    for k, v in final_visuals.items():
        v.figure.write_html(f"{export_name_prefix}.{k}.html")
