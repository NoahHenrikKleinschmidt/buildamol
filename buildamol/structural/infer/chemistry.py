"""
Element and valence adjustment helpers.
"""

from .constants import acceptable_surplus_charge, element_connectivity
from .hydrogens import Hydrogenator


def has_free_valence(atom, bonds, needed: int = 1):
    """
    Check if an atom has free valence.
    """
    degree = sum(bond.order for bond in bonds)
    return degree <= element_connectivity.get(atom.element, 0) - needed


def change_element(atom, new_element: str, _molecule):
    """
    Change the element of an atom and update its connectivity accordingly.
    """
    atom.set_element(new_element, adjust_id=False)
    atom.set_charge(0)
    new_connectivity = element_connectivity.get(new_element, 0)

    if _molecule.get_degree(atom) > new_connectivity:
        neighbors = _molecule.get_neighbors(atom)
        to_remove = []
        loops = 0
        while _molecule.get_degree(atom) > new_connectivity:
            for neighbor in neighbors:
                if neighbor.element == "H":
                    to_remove.append(neighbor)
                    neighbors.remove(neighbor)
                    new_connectivity += 1
                    break
            loops += 1
            if loops > 100:
                leftover_charge = _molecule.get_degree(atom) - new_connectivity
                if leftover_charge > acceptable_surplus_charge:
                    raise ValueError(
                        f"Could not change element of atom {atom} to {new_element}. Not enough hydrogen atoms available to remove, and without hydrogens the atom would have a charge of {leftover_charge}! Maybe check your input or manually use `atom.set_element` to override the element."
                    )
                else:
                    atom.set_charge(leftover_charge)
                    break

        _molecule.remove_atoms(*to_remove)
    elif _molecule.get_degree(atom) < new_connectivity:
        hydrogenator = Hydrogenator()
        hydrogenator.add_hydrogens(atom, _molecule)

        if _molecule.get_degree(atom) < new_connectivity:
            leftover_charge = new_connectivity - _molecule.get_degree(atom)
            if abs(leftover_charge) > acceptable_surplus_charge:
                raise ValueError(
                    f"Could not change element of atom {atom} to {new_element}. Not enough hydrogen atoms available to add, and without hydrogens the atom would have a charge of {leftover_charge}! Maybe check your input or manually use `atom.set_element` to override the element."
                )
            else:
                atom.set_charge(leftover_charge)
