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
import buildamol.optimizers.environments as envs
import buildamol.optimizers.agents as agents

import time
import numpy as np
import pandas as pd

# %% [markdown]
# Here we can select which tests to run and on what testing structures:

# %%
# which structures to run on
structures_to_run_on = [
    # files.GLUCOSE2,
    # files.PEPTIDE,
    # files.X_WING,
    files.X_WING2,
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
clash_cutoff = 0.5

# visualize the final structure
visualize_final_structure = True

# visualization parameters
# for draw_edges()
visualization_params = dict(color="magenta", opacity=0.3)

# export visualizations
export_visualizations = True

# export history to csv
export_history = True

# export name prefix
export_name_prefix = "xwing_declash_genetic_no_early_stop_"


# graph building function
# provide a custom callable that generates a tuple of (graph, rotatable_edges
def make_graph(x):
    g = x._AtomGraph
    edges = [
        (x.get_atom(68), x.get_atom(65)),
        (x.get_atom(79), x.get_atom(68)),
        (x.get_atom(79), x.get_atom(151)),
        (x.get_atom(151), x.get_atom(148)),
        (x.get_atom(251), x.get_atom(240)),
        (x.get_atom(251), x.get_atom(323)),
        (x.get_atom(323), x.get_atom(320)),
        (x.get_atom(240), x.get_atom(237)),
    ]  # the bonds at the clashing sites...
    return g, edges


graph_factory = make_graph

# graph building parameters
graph_params = {}

# provide a custom callable to set a custom building function for the environment
rotatron_factory = None

# the rotatron class to use
rotatron_class = envs.DistanceRotatron

# rotatron parameters
rotatron_params = {"radius": -1}

# the agent function to use
agent = agents.genetic_optimize

# agent parameters
agent_params = {"threshold": 1e-6, "stop_if_done": False}


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
        print(f"Finished run {r+1} of {re_runs} on {structure.id} ({eval}, {t2-t1})")

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
