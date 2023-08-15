"""
This module contains utility functions for the optimizers.
"""

from typing import Union
import numpy as np
import biobuild.optimizers.Rotatron as Rotatron
import biobuild.optimizers.DistanceRotatron as DistanceRotatron

import biobuild.core.Molecule as Molecule
import biobuild.optimizers.algorithms as agents
import biobuild.utils.auxiliary as aux


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
        bond = mol.get_bonds(bond[0].full_id, bond[1].full_id)
        if len(bond) == 0:
            raise ValueError(
                f"Object and environment do not match (bond mismatch): {bond}"
            )
        bond = bond[0]
        mol.rotate_around_bond(
            *bond, angle, descendants_only=True, angle_is_degrees=False
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
        sol = sol[0]

    final = apply_solution(sol, env, mol)
    return final


def auto_algorithm(mol):
    """
    Decide which algorithm to use for a quick-optimize based on the molecule size.
    """
    if mol.count_atoms() < 500:
        if aux.HAS_RDKIT:
            return "rdkit"
    return "swarm"


__all__ = ["apply_solution", "optimize", "auto_algorithm"]
