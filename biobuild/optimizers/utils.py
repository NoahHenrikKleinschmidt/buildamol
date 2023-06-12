"""
This module contains utility functions for the optimizers.
"""

import numpy as np
import biobuild.optimizers.environments as environments
import biobuild.optimizers.agents as agents


def apply_solution(sol: np.ndarray, env: "environments.Rotatron", obj):
    """
    Apply a solution to an object

    Parameters
    ----------
    sol : np.ndarray
        The solution of rotational angles in radians to apply
    env : environments.Rotatron
        The environment used to find the solution
    obj
        The object to apply the solution to

    Returns
    -------
    obj
        The object with the solution applied
    """
    angles = np.degrees(sol)
    bonds = env.rotatable_edges

    if not len(angles) == len(bonds):
        raise ValueError(
            f"Solution and environment do not match (size mismatch): {len(angles)} != {len(bonds)}"
        )

    for i, bond in enumerate(bonds):
        angle = angles[i]
        bond = obj.get_bonds(bond[0].full_id, bond[1].full_id)
        if len(bond) == 0:
            raise ValueError(
                f"Object and environment do not match (bond mismatch): {bond}"
            )
        bond = bond[0]
        obj.rotate_around_bond(*bond, angle, descendants_only=True)

    return obj
