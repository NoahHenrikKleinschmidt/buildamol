"""
The main function to make simple nanotubes.
"""

# This code directly uses the stuff from the nanotube tutorial notebook
# (see documentation)

import numpy as np
import buildamol.core as core

from . import polycarbons

__all__ = ["nanotube"]


def nanotube(n: int, k: int) -> core.Molecule:
    """
    Make a graphene nanotube with n carbon atoms per ring layer

    Parameters
    ----------
    n : int
        The number of carbon atoms per ring layer
    k : int
        The number of ring layers

    Returns
    -------
    Molecule
        The nanotube
    """
    angle = np.radians(32)

    # we can approximate the ring radius by the length and number of carbon atoms
    # given the formula: u = 2 * r * pi , where the circumference u is roughly
    # length * ring_size * cos(angle) (otherwise the bond lengths would be off)
    radius = (np.cos(angle) * polycarbons._length_C_C) * (n + 1) / (2 * np.pi)

    # now we can calculate the coordinates for the carbon and hydrogen atoms
    # in the ring
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)

    carbon_x = radius * np.cos(theta)
    carbon_y = radius * np.sin(theta)

    z = np.zeros(n)
    z[1::2] = np.sin(angle) * polycarbons._length_C_H

    # assembly half the unit ring
    ring = core.Molecule.new(id="RING", resname="RNG")
    Cs, cbonds = polycarbons._make_carbons(carbon_x, carbon_y, z)
    cbonds.append((f"C1", f"C{n}"))

    ring.add_atoms(*Cs)
    ring.add_bonds(*cbonds)

    # add the rest of the ring by copying and flipping the first half
    ring2 = ring.copy().move([0, 0, polycarbons._length_C_C]).flip("xy")
    ring.add_residues(ring2.get_residue(1))
    ring.add_bonds(*ring2.get_bonds())

    for i in range(0, n, 2):
        c = f"C{i+1}"
        a = ring.get_atom(c, residue=1)
        b = ring.get_atom(c, residue=2)
        ring.add_bond(a, b, 2)

    _min = min(i.coord[2] for i in ring.atoms)
    ring.move([0, 0, -_min])
    height = max(i.coord[2] for i in ring.atoms)
    dz = polycarbons._length_C_C + height  # the distance between the rings

    # =================================

    # now we can stack the rings to form the nanotube
    nanotube = ring.copy()

    for i in range(k - 1):

        # get a new ring and move it to the correct position
        # (we move downward by the length of the ring...)
        other = ring.copy().move([0, 0, -dz * (i + 1)])

        # add it to the nanotube
        nanotube.add_residues(*other.get_residues())
        nanotube.add_bonds(*other.get_bonds())

        # and add bonds between the incoming and the previous ring
        # (careful, here we use a shifted index, compared to the previous cell)
        for j in range(1, n + 1, 2):
            c = f"C{j+1}"
            a = nanotube.get_atom(c, residue=-2)
            b = nanotube.get_atom(c, residue=-3)
            nanotube.add_bond(a, b, 2)

    return nanotube


if __name__ == "__main__":
    nanotube(10)
