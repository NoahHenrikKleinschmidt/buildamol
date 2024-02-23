"""
This module contains utility functions for the optimizers.
"""

from typing import Union
import numpy as np
import buildamol.optimizers.Rotatron as Rotatron
import buildamol.optimizers.DistanceRotatron as DistanceRotatron

import buildamol.core.Molecule as Molecule
import buildamol.optimizers.algorithms as agents
import buildamol.utils.auxiliary as aux

from multiprocessing import Pool, cpu_count
from functools import partial


def apply_solution(
    sol: np.ndarray, env: "Rotatron.Rotatron", mol: "Molecule.Molecule"
) -> "Molecule.Molecule":
    """
    Apply a solution to a Molecule object.

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

    for i, bond in enumerate(bonds):
        angle = sol[i]
        a, b = mol.get_atom(bond[0].full_id), mol.get_atom(bond[1].full_id)
        if a is None or b is None:
            raise ValueError(
                f"Object and environment do not match (bond mismatch): {bond}"
            )

        mol.rotate_around_bond(
            a, b, angle, descendants_only=True, angle_is_degrees=False
        )

    return mol


def optimize(
    mol: "Molecule.Molecule",
    env: "Rotatron.Rotatron" = None,
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
    env : Rotatron, optional
        The environment to use. This needs to be a Rotatron instance that is fully set up and ready to use.
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
    algorithm = algorithm or auto_algorithm(mol)
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
        if mol.count_atoms() > 500:
            graph = mol.make_residue_graph()
            graph.make_detailed(n_samples=0.6)
            edges = mol.get_residue_connections()
            edges = graph.direct_edges(None, edges)
        else:
            graph = mol.make_atom_graph()
            edges = graph.find_rotatable_edges(min_descendants=5, max_descendants=10)

        env = DistanceRotatron(graph, edges, radius=25)

    sol, eval = agent(env, **kwargs)

    # in case of the genetic algorithm, the solution is a list of solutions
    # so we need to take the first one
    if sol.shape[0] != env.n_edges:
        if sol[0].shape[0] != env.n_edges:
            raise ValueError(
                f"Solution and environment do not match (size mismatch): {sol.shape[0]} != {env.n_edges}"
            )
        final = [apply_solution(s, env, mol.copy()) for s in sol]
        return final

    final = apply_solution(sol, env, mol)
    return final


def auto_algorithm(mol):
    """
    Decide which algorithm to use for a quick-optimize based on the molecule size.
    """
    # if mol.count_atoms() < 500:
    #     if aux.HAS_RDKIT:
    #         return "rdkit"
    return "swarm"


def split_environment(
    env: "Rotatron.Rotatron", n: int = None, radius: float = None
) -> list["Rotatron.Rotatron"]:
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

    kmeans = KMeans(n_clusters=n)
    kmeans.fit(all_edges)
    labels = kmeans.predict(all_edges)

    # create the sub-environments
    sub_envs = []
    for i in range(n):
        mask = labels == i
        edges = _all_edges[mask]
        sub_env = DistanceRotatron(
            env.graph,
            edges,
            radius=radius or env.radius,
            n_processes=env.n_processes,
            setup=False,
        )
        sub_env.edge_masks = env.edge_masks[mask]
        sub_env.edge_lengths = env.edge_lengths[mask]
        sub_envs.append(sub_env)

    return sub_envs


def parlallel_optimize(
    mol: "Molecule.Molecule",
    envs: list["Rotatron.Rotatron"],
    algorithm: Union[str, callable] = None,
    n_processes: int = None,
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

    for res, env in zip(results, envs):
        sol, eval = res
        if sol.shape[0] != env.n_edges:
            raise ValueError(
                f"Solution and environment do not match (size mismatch): {sol.shape[0]} != {env.n_edges}"
            )
        mol = apply_solution(sol, env, mol)

    return mol


__all__ = [
    "apply_solution",
    "optimize",
    "auto_algorithm",
    "split_environment",
    "parlallel_optimize",
]

if __name__ == "__main__":

    import buildamol as bam
    import time

    mol = bam.Molecule.from_pdb(
        "/Users/noahhk/GIT/biobuild/__figure_makery/EX8_new.pdb"
    )

    print("making first environment")

    rotatron = bam.optimizers.DistanceRotatron(
        mol.make_residue_graph(True), n_processes=6
    )

    print("splitting environment")
    N = 2
    split = split_environment(rotatron, N)
    multi_final = parlallel_optimize(mol.copy(), split, n_processes=6)
