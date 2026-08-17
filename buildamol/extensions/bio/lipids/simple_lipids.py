import buildamol.core as core
import buildamol.resources as resources
from buildamol.structural import Hydrogenator
from buildamol.extensions.polymers import linear_alkane
from typing import Union

import numpy as np

__all__ = [
    "fatty_acid",
    "triacylglycerol",
    "phospholipid",
    "sphingolipid",
    "orient_for_membrane",
]


def fatty_acid(
    length: int,
    double_bonds: Union[int, tuple],
    cis: Union[float, tuple],
    id: str = "UNK",
) -> core.Molecule:
    """
    Create a fatty acid molecule.

    Parameters
    ----------
    length : int
        The length of the fatty acid chain.
    double_bonds : int or tuple
        Double bonds to make. This can be either a single integer
        (the number of double bonds to make at random positions)
        or a tuple of integers (the positions of the double bonds).
    cis : float or tuple
        The cis configuration of the double bonds. This can be either
        a single float (the probability of a double bond being cis)
        or a tuple of floats (the probability of each double bond being cis).
        Boolean values are also accepted instead of floats.
    id : str
        The id of the fatty acid molecule.

    Returns
    -------
    Molecule
        The fatty acid molecule.

    Examples
    --------
    \# Create a fatty acid with 20 carbons, 4 double bonds, and 50% in cis configuration:
    >>> mol = fatty_acid(20, 4, 0.5)
    # Create a fatty acid with 16 carbons, and double bonds at positions (6 and 12) in trans and cis configuration:
    >>> mol = fatty_acid(16, (6, 12), (False, True))
    # Create a fatty acid with 12 carbons and two double bonds both in cis configuration:
    >>> mol = fatty_acid(12, 2, True)
    """
    # make a linear alkane
    mol = linear_alkane(length)

    # select some random positions for double bonds if not provided
    if isinstance(double_bonds, int):
        possible = list(range(2, length - 1))
        selected = []
        for _ in range(double_bonds):
            db = np.random.choice(possible)
            selected.append(db)
            possible.remove(db)
            # don't allow two double bonds in a row
            if db + 1 in possible:
                possible.remove(db + 1)
            if db - 1 in possible:
                possible.remove(db - 1)
        double_bonds = selected

    # select some random cis configurations if not provided
    if isinstance(cis, (float, int, bool)):
        cis = np.random.choice(
            [True, False], len(double_bonds), p=[float(cis), float(1 - cis)]
        )

    # now add double bonds into the molecule
    carbons = mol.get_atoms("C", by="element")
    for i, db in enumerate(double_bonds):
        mol.set_bond_order(carbons[db], carbons[db + 1], order=2, adjust_hydrogens=True)
        prob = cis[i]
        if prob:
            mol.cis(carbons[db], carbons[db + 1])

    # Add a carboxylic acid group
    core.carboxylate(mol, at_atom="C1")

    # now be sure to rename the atoms
    for idx, carbon in enumerate(carbons):
        carbon.id = f"C{idx+1}"
        for hydrogen in mol.get_hydrogens(carbon):
            hydrogen.id = f"H{idx+1}{hydrogen.id[-1]}"

    # make it all just one residue
    mol = mol.squash()
    mol.rename_residue(1, id)
    mol.id = id
    return mol


def triacylglycerol(
    chain1: core.Molecule,
    chain2: core.Molecule,
    chain3: core.Molecule,
) -> core.Molecule:
    """
    Create a triacylglycerol molecule from three fatty acid chains.

    Parameters
    ----------
    chain1 : Molecule
        The first fatty acid chain. None can be provided, to leave this position empty.
    chain2 : Molecule
        The second fatty acid chain (the one in the middle). None can be provided, to leave this position empty.
    chain3 : Molecule
        The third fatty acid chain. None can be provided, to leave this position empty.

    Returns
    -------
    Molecule
        The triacylglycerol molecule.
    """
    resources.load_lipids()
    out = core.Molecule.from_compound("GOL")

    link = core.linkage(None, "C", delete_in_source=("OXT", "HXT"))
    _chains = (chain1, chain2, chain3)
    _i = 1
    for i in range(1, 4):
        if _chains[i - 1] is None:
            continue
        link.atom1 = f"O{i}"
        out.attach(_chains[i - 1], link, at_residue=1)
        _i += 1
    return out


def phospholipid(
    chain1: core.Molecule,
    chain2: core.Molecule,
    headgroup: core.Molecule,
    headgroup_link: core.Linkage,
    id: str = "UNK",
):
    """
    Create a phospholipid molecule from two fatty acid chains and a headgroup.

    Parameters
    ----------
    chain1 : Molecule
        The first fatty acid chain. None can be provided, to leave this position empty.
    chain2 : Molecule
        The second fatty acid chain. None can be provided, to leave this position empty.
    headgroup : Molecule
        The headgroup of the phospholipid. This will be attached to the phosphate group.
        It is assumed that the headgroup does NOT have its own phosphate group. None can be provided, to leave this position empty.
    headgroup_link: Linkage
        The linkage to use to attach the headgroup to the phosphate group. The phosphate group is treated
        as "target" the headgroup is "source".
    id : str
        The id of the phospholipid molecule.

    Returns
    -------
    Molecule
        The phospholipid molecule.
    """
    out = linear_alkane(3)

    # now attach a phosphate
    out.change_element("H33", "O").rename_atom(out.get_atom("O", by="element"), "O3")
    core.phosphorylate(out, at_atom="O3")
    H = Hydrogenator()
    H.add_hydrogens(out.get_atom("O3", residue=2), out)

    # attach the two fatty acid chains
    link = core.linkage(None, "OXT")
    _chains = (chain1, chain2)
    _i = 2
    for i in (1, 2):
        if _chains[i - 1] is None:
            continue
        link.atom1 = f"C{i}"
        out.attach(_chains[i - 1], link, at_residue=1)
        _i += 1
        OXT = out.get_atom("OXT", residue=_i)
        out.rename_atom(OXT, f"O{i}")

    if headgroup is not None:

        # now attach the headgroup (making sure to rename the target atom if necessary)
        if headgroup_link.atom1 is None or headgroup_link.atom1 not in ("O3", "O4"):
            headgroup_link.atom1 = "O3"

        out.get_atom(headgroup_link.atom1, residue=2).charge = 0
        out.attach(headgroup, headgroup_link, at_residue=2)

    out.rename_residue(1, id)
    out.id = id
    return out


def sphingolipid(
    chain: core.Molecule,
    headgroup: core.Molecule,
    headgroup_link: core.Linkage,
    id: str = "UNK",
):
    """
    Create a sphingolipid molecule from a fatty acid chain and a headgroup.

    Parameters
    ----------
    chain : Molecule
        The fatty acid chain. None can be provided, to leave this position empty.
    headgroup : Molecule
        The headgroup of the sphingolipid. None can be provided, to leave this position empty.
    headgroup_link: Linkage
        The linkage to use to attach the headgroup to the sphingosine. The sphingosine is treated
        as "target" the headgroup is "source".
    id : str
        The id of the sphingolipid molecule.

    Returns
    -------
    Molecule
        The sphingolipid molecule.
    """
    resources.load_lipids()
    out = core.Molecule.from_compound("SPH")

    if chain is not None:

        # attach the fatty acid chain
        link = core.linkage("N2", "C", delete_in_source=("OXT", "HXT"))
        out.attach(chain, link, at_residue=1)

    if headgroup is not None:

        # now attach the headgroup (making sure to rename the target atom if necessary)
        if headgroup_link.atom1 is None or headgroup_link.atom1 != "O1":
            headgroup_link.atom1 = "O1"
        out.attach(headgroup, headgroup_link, at_residue=1)

    out.rename_residue(1, id)
    out.id = id
    return out


def orient_for_membrane(
    lipid: core.Molecule,
    tail_down: bool = True,
    inplace: bool = True,
) -> core.Molecule:
    """
    Orient a lipid molecule for membrane assembly.

    Rotates the molecule so that its hydrophobic tails point along the -z axis
    (or +z if *tail_down* is False) and the polar head group points in the
    opposite direction.  The resulting orientation matches a lipid sitting in
    the upper leaflet (tail_down=True) or lower leaflet (tail_down=False) of a
    bilayer and is a prerequisite for arranging lipids into a membrane grid.

    The orientation is determined entirely from the 3-D coordinates:

    1. **Head centroid** — centroid of all phosphorus atoms (P); if none are
       present, nitrogen atoms (N) are used; if still none, all oxygen atoms (O).
    2. **Tail centroid** — centroid of the terminal methyl carbons, i.e. carbon
       atoms that have no heteroatom neighbours and at most one carbon neighbour
       (the ω-carbon at the end of each acyl chain).  Carbons closer than half
       the maximum distance to the head centroid are excluded to avoid picking
       up short branches near the glycerol backbone.
    3. A single rigid-body rotation aligns the head→tail vector to ±z.

    Parameters
    ----------
    lipid : Molecule
        An assembled lipid molecule produced by :func:`phospholipid`,
        :func:`triacylglycerol`, :func:`sphingolipid`, or similar.
    tail_down : bool
        If ``True`` (default) the tails point in the -z direction (upper leaflet).
        If ``False`` the tails point in the +z direction (lower leaflet).
    inplace : bool
        If ``True`` (default) the molecule is modified in place and returned.
        If ``False`` a deep copy is rotated and returned, leaving the original
        unchanged.

    Returns
    -------
    Molecule
        The oriented lipid molecule.

    Raises
    ------
    ValueError
        If no heteroatoms are found (cannot locate head group) or if no
        terminal chain carbons are found.

    Examples
    --------
    >>> chain1 = fatty_acid(16, 0, [])
    >>> chain2 = fatty_acid(18, 1, [True])
    >>> import buildamol.resources as resources
    >>> resources.load_sugars()
    >>> import buildamol.core as core
    >>> choline = core.Molecule.from_compound("CHO")
    >>> link = core.linkage(None, "C1", delete_in_source=["O1", "HO1"])
    >>> pl = phospholipid(chain1, chain2, choline, link, id="DPPC")
    >>> oriented = orient_for_membrane(pl)
    """
    if not inplace:
        lipid = lipid.copy()

    atoms = list(lipid.get_atoms())

    # --- locate the polar head group ---
    head_atoms = [a for a in atoms if a.element == "P"]
    if not head_atoms:
        head_atoms = [a for a in atoms if a.element == "N"]
    if not head_atoms:
        head_atoms = [a for a in atoms if a.element == "O"]
    if not head_atoms:
        raise ValueError("Cannot determine head group: no heteroatoms (P, N, O) found.")
    head_centroid = np.mean([a.coord for a in head_atoms], axis=0)

    # --- locate terminal tail carbons (ω-methyls) ---
    # Qualifying atoms: element C, no heavy heteroatom neighbours, ≤1 C neighbour.
    terminal_carbons = []
    for a in atoms:
        if a.element != "C":
            continue
        neighbors = lipid.get_neighbors(a)
        if any(n.element not in ("C", "H") for n in neighbors):
            continue
        if sum(1 for n in neighbors if n.element == "C") <= 1:
            terminal_carbons.append(a)

    if not terminal_carbons:
        raise ValueError(
            "Cannot identify terminal chain carbons. "
            "Ensure the molecule contains aliphatic chain termini."
        )

    # keep only the distal termini (distance ≥ 50 % of the maximum)
    dists = np.array(
        [np.linalg.norm(a.coord - head_centroid) for a in terminal_carbons]
    )
    tail_atoms = [a for a, d in zip(terminal_carbons, dists) if d >= dists.max() * 0.5]
    tail_centroid = np.mean([a.coord for a in tail_atoms], axis=0)

    # --- compute and apply rotation ---
    tail_vec = tail_centroid - head_centroid
    norm = np.linalg.norm(tail_vec)
    if norm < 1e-10:
        raise ValueError(
            "Head and tail centroids coincide; cannot determine orientation."
        )
    tail_vec /= norm

    target = np.array([0.0, 0.0, -1.0 if tail_down else 1.0])

    axis = np.cross(tail_vec, target)
    axis_norm = np.linalg.norm(axis)
    if axis_norm > 1e-10:
        angle = np.degrees(
            np.arccos(float(np.clip(np.dot(tail_vec, target), -1.0, 1.0)))
        )
        lipid.rotate(angle, axis / axis_norm, center=lipid.center_of_geometry)

    return lipid


if __name__ == "__main__":

    mol1 = fatty_acid(20, 4, 0.5)

    mol2 = fatty_acid(20, 4, 1.0)

    mol3 = fatty_acid(20, 4, (True, False, True, False))

    tag = triacylglycerol(mol1, mol2, mol3)
    tag2 = triacylglycerol(mol1, None, mol3)
    tag2.show()

    resources.load_sugars()
    glc = core.Molecule.from_compound("GLC")
    link = core.linkage(None, "C1", delete_in_source=["O1", "HO1"])
    phos = phospholipid(mol1, mol3, glc.copy(), link)

    sphing = sphingolipid(mol2, glc, link)
    sphing.show()

    orient_for_membrane(phos.copy(), tail_down=True).show()
