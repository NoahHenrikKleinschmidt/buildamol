"""
These are constructors for simple polycarbon chains.
"""

# These functions were designed for speed and not for
# beauty. There would have been much easier ways of assembling
# the molecules at the cost of slower (probably just a little)
# computation times.

import numpy as np

import buildamol.core as core
import buildamol.structural as structural

__all__ = ["linear_alkane", "cyclic_alkane", "linear_alkene"]

_length_C_C = structural.single_bond_lengths["C"]["C"]
_length_C_H = structural.single_bond_lengths["C"]["H"]


def linear_alkane(n: int, include_hydrogens: bool = True) -> core.Molecule:
    """
    Construct a linear alkane change with n carbon atoms.

    Parameters
    ----------
    n : int
        The number of carbon atoms.
    include_hydrogens : bool, optional
        If True, infer hydrogens, by default True

    Returns
    -------
    bam.Molecule
        The linear alkane chain.
    """
    if n < 1:
        raise ValueError("The number of atoms must be at least 1")

    molecule = core.Molecule.new(id="ALK", resname="ALK")

    carbon_xs = np.arange(n) * _length_C_C
    carbon_zs = np.zeros(n, dtype=np.float64)
    angle = np.radians(35)
    carbon_zs[1::2] = np.sin(angle) * _length_C_C

    Cs, bonds = _make_carbons(carbon_xs, np.zeros(n, dtype=np.float64), carbon_zs)
    molecule.add_atoms(*Cs)
    molecule.add_bonds(*bonds)

    if include_hydrogens:
        angle = structural.geometry.Tetrahedral.dihedral / 2 + 0.2 * np.pi

        hydrogen_ys1 = np.zeros(n, dtype=np.float64)
        hydrogen_ys1 -= np.sin(angle) * _length_C_H
        hydrogen_ys2 = np.zeros(n, dtype=np.float64)
        hydrogen_ys2 += np.sin(angle) * _length_C_H

        hydrogen_zs = carbon_zs.copy()
        mask = hydrogen_zs == 0
        hydrogen_zs[mask] -= _length_C_H
        hydrogen_zs[~mask] += _length_C_H

        hydrogen_xs = np.tile(carbon_xs, 2)
        hydrogen_zs = np.tile(hydrogen_zs, 2)
        hydrogen_ys = np.stack([hydrogen_ys1, hydrogen_ys2], axis=0).flatten()

        Hs, bonds = _make_hydrogens(hydrogen_xs, hydrogen_ys, hydrogen_zs, 2)
        molecule.add_atoms(*Hs)
        molecule.add_bonds(*bonds)

        # add terminal hydrogens
        v = structural.norm_vector(*molecule.get_atoms("C1", "C2"))

        if n % 2 == 0:
            v[-1] *= -1
        H = core.Atom(
            f"H{n}3",
            coord=(molecule.get_atom(f"C{n}").coord + v * _length_C_H),
            element="H",
        )
        molecule.add_atoms(H)
        molecule.add_bond(f"H{n}3", f"C{n}")

        if n % 2 == 0:
            v[-1] *= -1
        v[0] *= -1
        H = core.Atom(
            "H13",
            coord=(molecule.get_atom("C1").coord + v * _length_C_H),
            element="H",
        )
        molecule.add_atoms(H)
        molecule.add_bond("H13", "C1")

    return molecule


def cyclic_alkane(n: int, include_hydrogens: bool = True) -> core.Molecule:
    """
    Construct a cyclic alkane change with n carbon atoms.

    Parameters
    ----------
    n : int
        The number of carbon atoms.
    include_hydrogens : bool, optional
        If True, infer hydrogens, by default True

    Returns
    -------
    bam.Molecule
        The cyclic alkane chain.
    """
    if n < 3:
        raise ValueError("The number of atoms must be at least 3")

    molecule = core.Molecule.new(id="ALK", resname="ALK")

    # a radius estimate
    angle = np.radians(35)
    radius = (np.cos(angle) * _length_C_C) * (n + 1) / (2 * np.pi)

    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xs = radius * np.cos(theta)
    ys = radius * np.sin(theta)
    zs = np.zeros(n, dtype=np.float64)
    zs[1::2] = np.sin(angle) * _length_C_C

    Cs, bonds = _make_carbons(xs, ys, zs)
    bonds.append((f"C1", f"C{n}"))

    molecule.add_atoms(*Cs)
    molecule.add_bonds(*bonds)

    if include_hydrogens:
        _hydrogenator = structural.Hydrogenator()
        _hydrogenator.infer_hydrogens(molecule, bond_length=_length_C_H)

    return molecule


def linear_alkene(n: int, include_hydrogens: bool = True) -> core.Molecule:
    """
    Construct a linear alkene change with n carbon atoms.

    Parameters
    ----------
    n : int
        The number of carbon atoms.
    include_hydrogens : bool, optional
        If True, infer hydrogens, by default True

    Returns
    -------
    bam.Molecule
        The linear alkene chain.
    """
    if n < 2:
        raise ValueError("The number of atoms must be at least 2")

    molecule = core.Molecule.new(id="ALK", resname="ALK")

    carbon_xs = np.arange(n) * _length_C_C
    carbon_zs = np.zeros(n, dtype=np.float64)
    angle = np.radians(35)
    carbon_zs[1::2] = np.sin(angle) * _length_C_C

    Cs, bonds = _make_carbons(carbon_xs, np.zeros(n, dtype=np.float64), carbon_zs)
    molecule.add_atoms(*Cs)
    molecule.add_bonds(*bonds)

    for i, b in enumerate(molecule.get_bonds()):
        if i % 2 == 0:
            b.double()

    if include_hydrogens:
        hydrogen_zs = carbon_zs.copy()
        mask = hydrogen_zs == 0
        hydrogen_zs[mask] -= _length_C_H
        hydrogen_zs[~mask] += _length_C_H

        Hs, bonds = _make_hydrogens(
            carbon_xs, np.zeros(n, dtype=np.float64), hydrogen_zs, 1
        )
        molecule.add_atoms(*Hs)
        molecule.add_bonds(*bonds)

        # add terminal hydrogens
        v = structural.norm_vector(*molecule.get_atoms("C1", "C2"))

        if n % 2 == 0:
            v[-1] *= -1
        H = core.Atom(
            f"H{n}B",
            coord=(molecule.get_atom(f"C{n}").coord + v * _length_C_H),
            element="H",
        )
        molecule.add_atoms(H)
        molecule.add_bond(f"H{n}B", f"C{n}")

        if n % 2 == 0:
            v[-1] *= -1
        v[0] *= -1
        H = core.Atom(
            "H1B",
            coord=(molecule.get_atom("C1").coord + v * _length_C_H),
            element="H",
        )
        molecule.add_atoms(H)
        molecule.add_bond("H1B", "C1")

        molecule.rename_atom("H1", "H1A")
        molecule.rename_atom(f"H{n}", f"H{n}A")

    return molecule


def _make_hydrogens(coords_x, coords_y, coords_z, n) -> list:
    """
    Make a list of hydrogen atoms and create bonds for each hydrogen atom to a carbon with the same id given the number of carbons per hydrogen (n)
    """
    k = 1
    j = 1
    Hs = []
    bonds = []
    if n > 1:
        addon = lambda x: x
    else:
        addon = lambda x: ""
    for i in range(0, len(coords_x)):

        H = core.Atom(
            f"H{j}{addon(k)}",
            coord=(coords_x[i], coords_y[i], coords_z[i]),
            element="H",
        )
        Hs.append(H)
        bonds.append((f"H{j}{addon(k)}", f"C{j}"))

        j += 1

        if n > 1:
            if j == len(coords_x) / n + 1:
                j = 1
                k += 1

    return Hs, bonds


def _make_carbons(coords_x, coords_y, coords_z) -> list:
    """
    Make a list of carbon atoms and bonds from one atom to the next.
    Returns the list of atoms and the list of bonds.
    """
    bonds = []
    i = 0
    C = core.Atom("C1", coord=(coords_x[i], coords_y[i], coords_z[i]), element="C")
    Cs = [C]
    for i in range(1, len(coords_x)):
        C = core.Atom(
            f"C{i+1}", coord=(coords_x[i], coords_y[i], coords_z[i]), element="C"
        )
        Cs.append(C)
        bonds.append((f"C{i+1}", f"C{i}"))
    return Cs, bonds


if __name__ == "__main__":
    alkane = linear_alkane(25)
    alkane.to_pdb("alkane.pdb")
    # alkane.to_pdb("alkane.pdb")
    # v = alkane.draw()
    # v.viewbox(None, 20, 20)
    # v.show()
