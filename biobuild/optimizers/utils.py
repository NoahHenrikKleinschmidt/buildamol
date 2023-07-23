"""
This module contains utility functions for the optimizers.
"""

import numpy as np
import biobuild.optimizers.environments.Rotatron as Rotatron
import biobuild.optimizers.agents as agents


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

    if not len(angles) == len(bonds):
        raise ValueError(
            f"Solution and environment do not match (size mismatch): {len(angles)} != {len(bonds)}"
        )

    for i, bond in enumerate(bonds):
        angle = angles[i]
        bond = mol.get_bonds(bond[0].full_id, bond[1].full_id)
        if len(bond) == 0:
            raise ValueError(
                f"Object and environment do not match (bond mismatch): {bond}"
            )
        bond = bond[0]
        mol.rotate_around_bond(
            *bond, angle, descendants_only=True, angle_is_degrees=False
        )

    return obj


__all__ = ["apply_solution"]
