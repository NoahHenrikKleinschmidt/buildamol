"""
Constants and basic lookup helpers for structural inference.
"""

import periodictable as pt


element_connectivity = {
    "C": 4,
    "H": 1,
    "O": 2,
    "N": 3,
    "S": 2,
    "P": 5,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "B": 3,
    "Si": 4,
    "Se": 2,
    "Zn": 2,
    "Ca": 2,
    "Mg": 2,
    "Fe": 2,
    "Cu": 1,
    "Mn": 2,
}


element_to_hydrogen_bond_lengths = {
    "C": 1.09,
    "N": 1.01,
    "O": 0.96,
    "S": 1.04,
    "P": 1.10,
    "F": 0.92,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
    "B": 1.19,
    "Si": 1.17,
    "Se": 1.17,
    "Zn": 1.31,
    "Ca": 1.74,
    "Mg": 1.36,
    "Fe": 1.24,
    "Cu": 1.28,
    "Mn": 1.20,
}


single_bond_lengths = {
    "C": {
        "C": 1.54,
        "N": 1.47,
        "O": 1.43,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 1.09,
    },
    "N": {
        "C": 1.47,
        "N": 1.45,
        "O": 1.43,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 1.01,
    },
    "O": {
        "C": 1.43,
        "N": 1.43,
        "O": 1.34,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 0.96,
    },
    "S": {
        "C": 1.81,
        "N": 1.81,
        "O": 1.81,
        "S": 2.05,
        "P": 2.05,
        "F": 1.81,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.04,
    },
    "P": {
        "C": 1.87,
        "N": 1.87,
        "O": 1.87,
        "S": 2.05,
        "P": 2.05,
        "F": 1.87,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.10,
    },
    "F": {
        "C": 1.35,
        "N": 1.35,
        "O": 1.35,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 0.92,
    },
    "Cl": {
        "C": 1.77,
        "N": 1.77,
        "O": 1.77,
        "S": 2.05,
        "P": 2.05,
        "F": 1.77,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 0.99,
    },
    "Br": {
        "C": 1.94,
        "N": 1.94,
        "O": 1.94,
        "S": 2.05,
        "P": 2.05,
        "F": 1.94,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.14,
    },
    "I": {
        "C": 2.14,
        "N": 2.14,
        "O": 2.14,
        "S": 2.05,
        "P": 2.05,
        "F": 2.14,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.33,
    },
    "H": {
        "C": 1.09,
        "N": 1.01,
        "O": 0.96,
        "S": 1.04,
        "P": 1.10,
        "F": 0.92,
        "Cl": 0.99,
        "Br": 1.14,
        "I": 1.33,
    },
}


double_bond_lengths = {
    "C": {
        "C": 1.34,
        "N": 1.30,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.34,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "N": {
        "C": 1.30,
        "N": 1.28,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.30,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "O": {
        "C": 1.21,
        "N": 1.21,
        "O": 1.20,
        "S": 2.05,
        "P": 2.05,
        "F": 1.21,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "S": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "P": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "F": {
        "C": 1.34,
        "N": 1.30,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.34,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "Cl": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "Br": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "I": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
}


triple_bond_lengths = {
    "C": {
        "C": 1.20,
        "N": 1.16,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.20,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "N": {
        "C": 1.16,
        "N": 1.14,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.16,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "O": {
        "C": 1.13,
        "N": 1.13,
        "O": 1.12,
        "S": 2.05,
        "P": 2.05,
        "F": 1.13,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "S": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.51,
        "P": 2.51,
        "F": 2.05,
        "Cl": 2.51,
        "Br": 2.51,
        "I": 2.51,
    },
    "P": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.51,
        "P": 2.51,
        "F": 2.05,
        "Cl": 2.51,
        "Br": 2.51,
        "I": 2.51,
    },
    "F": {
        "C": 1.20,
        "N": 1.16,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.20,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
}


bond_length_by_order = {
    1: single_bond_lengths,
    2: double_bond_lengths,
    3: triple_bond_lengths,
}


element_vdw_radii = {
    "H": 1.2,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "P": 1.8,
    "F": 1.47,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    "B": 1.92,
    "Si": 2.1,
    "Se": 1.9,
    "Zn": 1.39,
    "Ca": 2.31,
    "Mg": 1.73,
    "Fe": 1.94,
    "Cu": 1.4,
    "Mn": 2.0,
}

_VDW_DEFAULT_RADIUS = 1.8
_VDW_SCALE = 0.45
_VDW_PADDING = 0.1
_MIN_BOND_LENGTH = 0.4


acceptable_surplus_charge = 5
"""
The maximum allowed positive charge that can be left on an atom
if no hydrogens can be further removed to balance the charge.

This is used in `change_element` to prevent the creation of
unrealistically charged atoms and thereby unlikely structures.
"""


def _get_vdw_radius(element: str) -> float:
    return element_vdw_radii.get(element, _VDW_DEFAULT_RADIUS)


def _bond_cutoff_vdw(
    atom1, atom2, scale: float = _VDW_SCALE, padding: float = _VDW_PADDING
) -> float:
    return (
        scale * (_get_vdw_radius(atom1.element) + _get_vdw_radius(atom2.element))
        + padding
    )


def _max_search_radius_vdw(
    atoms, scale: float = _VDW_SCALE, padding: float = _VDW_PADDING
) -> float:
    if not atoms:
        return 0.0
    max_r = max(_get_vdw_radius(atom.element) for atom in atoms)
    return scale * (2 * max_r) + padding


def atomic_number(element: str):
    """
    Get the atomic number of an element.

    Parameters
    ----------
    element : str
        The element symbol.

    Returns
    -------
    int
        The atomic number of the element.
    """
    return pt.elements.symbol(element.title()).number
