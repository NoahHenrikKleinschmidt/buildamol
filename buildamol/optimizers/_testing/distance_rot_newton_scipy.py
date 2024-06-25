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
    files.X_WING,
    files.X_WING2,
    # files.SCAFFOLD1,
    # files.SCAFFOLD3,
]

# how many times to independently run on each structure
re_runs = 10

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
export_name_prefix = None

# graph building function
# provide a custom callable that generates a tuple of (graph, rotatable_edges)
graph_factory = lambda x: (auxiliary.graph_factory(x, **graph_params)[0], None)

# graph building parameters
graph_params = {}

# provide a custom callable to set a custom building function for the environment
rotatron_factory = None

# the rotatron class to use
rotatron_class = envs.DistanceRotatron

# rotatron parameters
rotatron_params = {}

# the agent function to use
agent = agents.scipy_optimize

# agent parameters
agent_params = {}


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

structures_to_run_on = [bam.molecule(s) for s in structures_to_run_on]

eval_history = {structure.id: [] for structure in structures_to_run_on}
time_history = {structure.id: [] for structure in structures_to_run_on}
clash_history = {structure.id: [] for structure in structures_to_run_on}
final_visuals = {structure.id: None for structure in structures_to_run_on}
initial_evals = {structure.id: None for structure in structures_to_run_on}
initial_clashes = {structure.id: None for structure in structures_to_run_on}

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
    initial_evals[structure.id] = [env.step(env.blank())[1]] * re_runs
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
        if visualize_final_structure:
            final = auxiliary.apply_solution(sol, env, structure.copy())
            clash_history[structure.id].append(
                auxiliary.count_clashes(final, clash_cutoff)
            )
            v.draw_edges(*final.bonds, **visualization_params)
        else:
            clash_history[structure.id].append(
                auxiliary.count_clashes(env.graph, clash_cutoff)
            )
        env.reset()

    if visualize_final_structure:
        _best = auxiliary.apply_solution(env.best[0], env, structure.copy())
        v.draw_edges(*_best.bonds, color="green", linewidth=6)
        final_visuals[structure.id] = v
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
