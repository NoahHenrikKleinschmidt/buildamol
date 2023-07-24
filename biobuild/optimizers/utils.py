"""
This module contains utility functions for the optimizers.
"""

import numpy as np
import biobuild.optimizers.environments.Rotatron as Rotatron
import biobuild.optimizers.environments.DistanceRotatron as DistanceRotatron

import biobuild.optimizers.algorithms as agents


def apply_solution(sol: np.ndarray, env: "Rotatron.Rotatron", mol: "Molecule"):
    """
    Apply a solution to a Molecule object.

    Parameters
    ----------
    sol : np.ndarray
        The solution of rotational angles in radians to apply
    env : environments.Rotatron
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


def quick_optimize(
    mol: "Molecule",
    env: "Rotatron.Rotatron" = None,
    algorithm: str = "genetic",
) -> "Molecule":
    """
    Quickly optimize a molecule using a specific algorithm.

    Note
    ----
    This is a convenience function that will automatically create an environment and determine edges.
    However, that means that the environment will be created from scratch every time this function is called.
    Also, the environment will likely not taylor to any specifc requiremehts of the molecule. For better performance
    and control, it is recommended to create an environment manually and supply it to the function using the `env` argument.

    Parameters
    ----------
    mol : Molecule
        The molecule to optimize. This molecule will be modified in-place.
    env : Rotatron.Rotatron, optional
        The environment to use, by default None
    algorithm : str, optional
        The algorithm to use, by default "genetic". This can be:
        - "genetic": A genetic algorithm
        - "swarm": A particle swarm optimization algorithm
        - "gradient": A gradient descent algorithm (default scipy implementation)

    Returns
    -------
    Molecule
        The optimized molecule
    """
    if env is None:
        if sum(1 for i in mol.get_atoms()) > 50:
            graph = mol.make_residue_graph()
            graph.make_detailed()
            edges = mol.get_residue_connections()
            edges = graph.direct_edges(None, edges)
        else:
            graph = mol.make_atom_graph()
            edges = mol.find_rotatable_edges(min_descendants=5, max_descendants=10)

        env = DistanceRotatron(graph, edges, radius=25)

    if algorithm == "genetic":
        agent = agents.genetic_optimize
    elif algorithm == "swarm":
        agent = agents.swarm_optimize
    elif algorithm == "gradient":
        agent = agents.scipy_optimize
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")

    sol, eval = agent(env)

    final = apply_solution(sol, env, mol)
    return final


__all__ = ["apply_solution", "quick_optimize"]
