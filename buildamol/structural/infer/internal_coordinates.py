"""
Internal-coordinate computation helpers.
"""

from collections import defaultdict
import re
import warnings

import numpy as np

import buildamol.structural.base as base
import buildamol.structural.neighbors as neighbors
import buildamol.utils.auxiliary as aux
import buildamol.utils.ic as _ic

from .constants import element_connectivity


def compute_internal_coordinates(bonds: list):
    """
    Compute internal coordinates for a structure.
    """
    quartets = neighbors.compute_quartets(bonds)
    ics = []
    for quartet in quartets:
        angle_123 = base.compute_angle(quartet[0], quartet[1], quartet[2])
        angle_234 = base.compute_angle(quartet[1], quartet[2], quartet[3])
        dihedral = base.compute_dihedral(quartet[0], quartet[1], quartet[2], quartet[3])
        l_12 = base.compute_distance(quartet[0], quartet[1])
        l_13 = base.compute_distance(quartet[0], quartet[2])
        l_34 = base.compute_distance(quartet[2], quartet[3])

        ic = _ic.InternalCoordinates(
            quartet[0].id,
            quartet[1].id,
            quartet[2].id,
            quartet[3].id,
            l_12,
            l_34,
            angle_123,
            angle_234,
            dihedral,
            l_13 if quartet.improper else None,
            improper=quartet.improper,
        )
        ics.append(ic)
    return ics


def compute_atom1_from_others(coords2, coords3, coords4, ic):
    """
    Compute coordinates of the first atom from internal coordinates and the other three atoms.
    """
    ic_to_xyz = (
        base._numba_wrapper_IC_to_xyz
        if (aux.USE_NUMBA or aux.USE_ALL_NUMBA)
        else base._IC_to_xyz
    )
    if ic.is_proper:
        return ic_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=coords2,
            r=-ic.bond_length_12,
            theta=-np.radians(ic.bond_angle_123),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        _vec = ic_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=np.full(3, 0),
            r=1,
            theta=np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
        if ic.bond_length_13:
            _vec *= ic.bond_length_13
        else:
            BC = np.linalg.norm(coords2 - coords3)
            AC = ic.bond_length_12
            AB = np.sqrt(
                BC**2 + AC**2 - 2 * BC * AC * np.cos(np.radians(ic.bond_angle_123))
            )
            _vec *= AB
        final = coords3 + _vec
        return final


def compute_atom4_from_others(coords1, coords2, coords3, ic):
    """
    Compute coordinates of the fourth atom from internal coordinates and the other three atoms.
    """
    ic_to_xyz = (
        base._numba_wrapper_IC_to_xyz
        if (aux.USE_NUMBA or aux.USE_ALL_NUMBA)
        else base._IC_to_xyz
    )
    if ic.is_proper:
        return ic_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        return ic_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )


def _prune_H_triplets(bonds):
    """
    Remove erroneous bonds that connect hydrogens to multiple other atoms.
    """
    bonds_with_H = [
        bond for bond in bonds if bond[0].element == "H" or bond[1].element == "H"
    ]
    triplets = neighbors.generate_triplets(bonds_with_H)
    bond_mappings = defaultdict(int)
    for a, b in bonds_with_H:
        bond_mappings[a] += 1
        bond_mappings[b] += 1

    for triplet in triplets:
        if triplet[1].element != "H":
            continue

        non_H1, H, non_H2 = triplet

        e_non_H1 = non_H1.element.title()
        e_non_H2 = non_H2.element.title()

        if bond_mappings[non_H1] > element_connectivity[e_non_H1]:
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif bond_mappings[non_H2] > element_connectivity[e_non_H2]:
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1

        elif _H_id_match(H, non_H1):
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif _H_id_match(H, non_H2):
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1

        elif _H_dist_match(H, non_H1):
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif _H_dist_match(H, non_H2):
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1
        else:
            warnings.warn(f"Could not prune H triplet! ({triplet=})", RuntimeWarning)
    return bonds


def _H_dist_match(H, non_H):
    """
    Check if the H-nonH distance suggests a potential bond.
    """
    d = np.linalg.norm(H.coord - non_H.coord)
    if non_H.element == "C":
        return 0.94 < d < 1.1
    elif non_H.element == "S":
        return 1.15 < d < 1.42
    elif non_H.element in ("O", "N"):
        return 0.94 < d < 1.07
    else:
        return 0.7 < d < 1.4


def _H_id_match(H, non_H):
    """
    Check if H atom id suggests belonging to a particular non-H atom.
    """
    if (
        non_H.element == "C"
        and re.match(r"C\d.*", non_H.id.upper()) is not None
        and re.match(r"H\d.*", H.id.upper()) is None
    ):
        return False
    elif non_H.element != "C" and re.match(r"H\d.*", H.id.upper()) is not None:
        return False
    return re.search(non_H.id[1:].upper(), H.id[1:].upper()) is not None
