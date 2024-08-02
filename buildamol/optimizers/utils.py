"""
This module contains utility functions for the optimizers.
"""

from typing import Union, List
import numpy as np
import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.optimizers.translatron as Translatron


import buildamol.core.Molecule as Molecule
import buildamol.optimizers.algorithms as agents

from multiprocessing import Pool, cpu_count
from functools import partial


__all__ = [
    "apply_rotatron_solution",
    "apply_translatron_solution",
    "optimize",
    "auto_algorithm",
    "split_environment",
    "parallel_optimize",
]


def apply_rotatron_solution(
    sol: np.ndarray, env: "Rotatron.Rotatron", mol: "Molecule.Molecule"
) -> "Molecule.Molecule":
    """
    Apply the solution of a Rotatron environment to a Molecule object.

    Parameters
    ----------
    sol : np.ndarray
        The solution of rotational angles in radians to apply
    env : Rotatron
        The environment used to find the solution
    mol : Molecule
        The molecule to apply the solution to

    Returns
    -------
    obj
        The object with the solution applied
    """
    bonds = env.rotatable_edges

    if not len(sol) == len(bonds):
        raise ValueError(
            f"Solution and environment do not match (size mismatch): {len(sol)} != {len(bonds)}"
        )

    mol._AtomGraph.clear_cache()
    for i, bond in enumerate(bonds):
        angle = sol[i]
        # used to be full_id to account for the fact that the bond might
        # come from another molecule. But there is no reason to assume someone
        # would apply the solutions of one molecule to another.
        # so we could simply use the serial_number instead of full_id which makes
        # things faster, assuming that the serial number was not altered in some way
        # outside of the molecule object.
        a, b = bond
        a = mol.get_atom(a)
        b = mol.get_atom(b)
        if a is None or b is None:
            raise ValueError(
                f"Object and environment do not match (bond mismatch): {bond}"
            )
        mol._AtomGraph.rotate_around_edge(a, b, angle, descendants_only=True)
    return mol


def apply_translatron_solution(
    sol: np.ndarray,
    env: "Translatron.Translatron",
    mol: "Molecule.Molecule",
):
    """
    Apply the solution of a Translatron environment to a Molecule object.

    Parameters
    ----------
    sol : np.ndarray
        The solution of translational vectors to apply
    env : Translatron
        The environment used to find the solution
    mol : Molecule
        The molecule to apply the solution to
    """
    mol.rotate(sol[3], "x", angle_is_degrees=False)
    mol.rotate(sol[4], "y", angle_is_degrees=False)
    mol.rotate(sol[5], "z", angle_is_degrees=False)
    mol.move(sol[:3])
    return mol


def optimize(
    mol: "Molecule.Molecule",
    env: Union["Rotatron.Rotatron", "Translatron.Translatron"] = None,
    algorithm: Union[str, callable] = None,
    **kwargs,
) -> "Molecule.Molecule":
    """
    Quickly optimize a molecule using a specific algorithm.

    Note
    ----
    This is a convenience function that will automatically create an environment and determine edges.
    However, that means that the environment will be created from scratch every time this function is called.
    Also, the environment will likely not taylor to any specifc requirements of the situation. For better performance
    and control, it is recommended to create an environment manually and supply it to the function using the `env` argument.

    Parameters
    ----------
    mol : Molecule
        The molecule to optimize. This molecule will be modified in-place.
    env : Rotatron or Translatron, optional
        The environment to use. This needs to be a Rotatron instance or Translatron instance that is fully set up and ready to use.
    algorithm : str or callable, optional
        The algorithm to use. If not provided, an algorithm is automatically determined, depending on the molecule size.
        If provided, this can be:
        - "genetic": A genetic algorithm
        - "swarm": A particle swarm optimization algorithm
        - "anneal": A simulated annealing algorithm
        - "scipy": A gradient descent algorithm (default scipy implementation, can be changed using a 'method' keyword argument)
        - "rdkit": A force field based optimization using RDKit (if installed)
        - or some other callable that takes an environment as its first argument
    **kwargs
        Additional keyword arguments to pass to the algorithm

    Returns
    -------
    Molecule
        The optimized molecule
    """
    algorithm = algorithm or auto_algorithm(mol, env)
    if algorithm == "genetic":
        agent = agents.genetic_optimize
    elif algorithm == "swarm":
        agent = agents.swarm_optimize
    elif algorithm == "scipy":
        agent = agents.scipy_optimize
    elif algorithm == "anneal":
        agent = agents.anneal_optimize
    elif algorithm == "rdkit":
        return agents.rdkit_optimize(mol, **kwargs)
    elif not callable(algorithm):
        raise ValueError(f"Unknown algorithm: {algorithm}")

    if env is None:
        from buildamol.optimizers.distance_rotatron import DistanceRotatron

        if mol.count_atoms() > 500:
            graph = mol.make_residue_graph()
            graph.make_detailed(n_samples=0.6)
            edges = mol.get_residue_connections()
            edges = graph.direct_edges(None, edges)
        else:
            graph = mol.make_atom_graph()
            edges = graph.find_rotatable_edges(min_descendants=5, max_descendants=10)

        env = DistanceRotatron(graph, edges, radius=25)

    applier = auto_applier(env)

    sol, eval = agent(env, **kwargs)

    # in case of the genetic algorithm, the solution is a list of solutions
    # so we need to take the first one
    n_edges = getattr(env, "n_edges", sol.shape[0])
    if sol.shape[0] != n_edges or "n_best" in kwargs:
        if sol[0].shape[0] != n_edges:
            raise ValueError(
                f"Solution and environment do not match (size mismatch): {sol.shape[0]} != {env.n_edges}"
            )
        final = [applier(s, env, mol.copy()) for s in sol]
        return final

    final = applier(sol, env, mol)
    return final


def auto_algorithm(mol, env=None):
    """
    Decide which algorithm to use for a quick-optimize based on the molecule size.
    """
    if env is not None:
        if isinstance(env, Rotatron.Rotatron):
            return "swarm"
        if isinstance(env, Translatron.Translatron):
            return "scipy"
    # if mol.count_atoms() < 500:
    #     if aux.HAS_RDKIT:
    #         return "rdkit"
    return "swarm"


def auto_applier(env):
    if isinstance(env, Rotatron.Rotatron):
        return apply_rotatron_solution
    if isinstance(env, Translatron.Translatron):
        return apply_translatron_solution
    else:
        raise ValueError(f"Unknown environment type: {env}")


def split_environment(
    env: "Rotatron.Rotatron", n: int = None
) -> List["Rotatron.Rotatron"]:
    """
    Split an environment into n sub-environments which are smaller and thus easier to optimize.

    Parameters
    ----------
    env : Rotatron
        The environment to split
    n : int
        The number of sub-environments to create

    Returns
    -------
    list
        A list of sub-environments
    """
    _all_edges = np.array(env.rotatable_edges)
    all_edges = np.array([(*i.coord, *j.coord) for i, j in env.rotatable_edges])

    # cluster the edges into n clusters
    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=n, n_init="auto")
    kmeans.fit(all_edges)
    labels = kmeans.predict(all_edges)

    # create the sub-environments
    sub_envs = []
    hyperparams = dict(env.hyperparameters)
    for i in range(n):
        mask = labels == i
        edges = _all_edges[mask].tolist()

        sub_env = env.__class__(env.graph, edges, setup=False, **hyperparams)

        sub_env.edge_masks = env.edge_masks[mask]
        sub_env.edge_lengths = env.edge_lengths[mask]
        sub_envs.append(sub_env)

    return sub_envs


def parallel_optimize(
    mol: "Molecule",
    envs: List["Rotatron"],
    algorithm: Union[str, callable] = None,
    n_processes: int = None,
    unify_final: bool = True,
    **kwargs,
) -> "Molecule.Molecule":
    """
    Optimize a molecule using multiple sub-environments in parallel.

    Parameters
    ----------
    mol : Molecule
        The molecule to optimize. This molecule will be modified in-place.
    envs : list
        The sub-environments to optimize
    algorithm : str or callable, optional
        The algorithm to use. If not provided, an algorithm is automatically determined, depending on the molecule size.
        If provided, this can be:
        - "genetic": A genetic algorithm
        - "swarm": A particle swarm optimization algorithm
        - "anneal": A simulated annealing algorithm
        - "scipy": A gradient descent algorithm (default scipy implementation, can be changed using a 'method' keyword argument)
        - or some other callable that takes an environment as its first argument
    n_processes : int, optional
        The number of processes to use. If not provided, the number of processes is automatically determined.
    unify_final : bool, optional
        If True, the solutions to all sub-environments are applied onto the same final molecule. If False,
        a list of molecules is returned each with the solution of one sub-environment applied.
    **kwargs
        Additional keyword arguments to pass to the algorithm

    Returns
    -------
    Molecule
        The optimized molecule
    """
    algorithm = algorithm or auto_algorithm(envs[0])
    if algorithm == "genetic":
        agent = partial(agents.genetic_optimize, **kwargs)
    elif algorithm == "swarm":
        agent = partial(agents.swarm_optimize, **kwargs)
    elif algorithm == "scipy":
        agent = partial(agents.scipy_optimize, **kwargs)
    elif algorithm == "anneal":
        agent = partial(agents.anneal_optimize, **kwargs)
    elif algorithm == "rdkit":
        raise ValueError("RDKit optimization is not supported in parallel mode")
    elif not callable(algorithm):
        raise ValueError(f"Unknown algorithm: {algorithm}")

    if n_processes is None:
        n_processes = cpu_count()

    n_processes = min(n_processes, len(envs))

    p = Pool(n_processes)
    results = p.map(agent, envs)
    p.close()
    p.join()

    if not unify_final:

        final = []
        for res, env in zip(results, envs):
            _mol = mol.copy()
            sol, eval = res
            if sol.shape[0] != env.n_edges:
                raise ValueError(
                    f"Solution and environment do not match (size mismatch): {sol.shape[0]} != {env.n_edges}"
                )
            _mol = apply_rotatron_solution(sol, env, _mol)
            final.append(_mol)
        return final

    for res, env in zip(results, envs):
        sol, eval = res
        if sol.shape[0] != env.n_edges:
            raise ValueError(
                f"Solution and environment do not match (size mismatch): {sol.shape[0]} != {env.n_edges}"
            )
        mol = apply_rotatron_solution(sol, env, mol)

    return mol


if __name__ == "__main__":

    import buildamol as bam
    import time

    mol = bam.Molecule.from_pdb(
        "/Users/noahhk/GIT/biobuild/__figure_makery/EX8_new.pdb"
    )

    print("making first environment")

    from pathlib import Path

    rotatron = bam.optimizers.OverlapRotatron
    graph = mol.get_atom_graph()
    edges = mol.get_residue_connections(triplet=False)
    edges = graph.direct_edges(edges=edges)
    edges = graph.sample_edges(edges, m=5, n=3)
    env = rotatron(graph, edges)
    envs = [env] * 2
    better = optimize(mol.copy(), env, n_best=10)

    rotatron = bam.optimizers.ForceFieldRotatron

    f = Path(__file__).parent / "{rotatron.__class__.__name__}-rotatron.pkl"
    f.unlink()
    if f.exists():
        rotatron = bam.utils.load_pickle(f)
    else:
        graph = mol.make_residue_graph(True)
        edges = mol.get_residue_connections(triplet=False)
        edges = graph.direct_edges(graph.central_node, edges)
        rotatron = rotatron(graph, edges, n_processes=6)

        bam.utils.save_pickle(rotatron, f)

    print("splitting environment")
    N = 2
    split = split_environment(rotatron, N)
    final = parallel_optimize(mol.copy(), split)
    final.to_pdb("_overlap_better.pdb")
