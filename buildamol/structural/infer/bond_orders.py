"""
Bond order inference methods.
"""

import warnings

import buildamol.utils.auxiliary as aux


def _infer_bond_orders_rdkit(molecule):
    """
    Infer bond orders using RDKit.
    """
    from rdkit import Chem

    rdmol = molecule.to_rdkit()
    Chem.Kekulize(rdmol, clearAromaticFlags=True)

    for bond in rdmol.GetBonds():
        a1 = molecule.get_atom(bond.GetBeginAtomIdx() + 1)
        a2 = molecule.get_atom(bond.GetEndAtomIdx() + 1)
        order = int(bond.GetBondTypeAsDouble())
        molecule.set_bond(a1, a2, order=order)


def _infer_bond_orders_native(molecule):
    """
    Infer bond orders using registered higher-order functional groups.
    """
    import buildamol.structural.groups as groups

    group_matches = {}
    _assigned_atoms = set()

    groups_to_apply = sorted(
        groups.higher_order_groups, key=lambda x: x.rank, reverse=True
    )
    for atom in molecule.get_atoms():
        neighbors = molecule.get_neighbors(atom)
        atoms = (atom, *neighbors)

        for group in groups_to_apply:
            if atoms in group_matches:
                continue
            if all(i in _assigned_atoms for i in atoms):
                continue
            if group.matches(molecule, atoms):
                if any(
                    i in _assigned_atoms
                    for i in group._assignment_cache[molecule][-1].values()
                ):
                    continue
                group_matches[atoms] = (group, group._assignment_cache[molecule][-1])
                _assigned_atoms.update(group._assignment_cache[molecule][-1].values())

    for atoms, (group, assignment) in group_matches.items():
        group._assignment = assignment
        group.apply_connectivity(molecule, atoms)


def infer_bond_orders(molecule, method: str = "native"):
    """
    Infer bond orders of a molecule.

    Parameters
    ----------
    molecule : Molecule
    method : str, optional
        The method to use for inferring bond orders. Default is "native", which uses registered higher-order functional groups. Another option is "rdkit", which uses RDKit's Kekulization.
        If "rdkit" is selected but RDKit is not installed, a RuntimeError will be raised.
    """
    if method == "rdkit":
        try:
            _infer_bond_orders_rdkit(molecule)
        except Exception as e:
            warnings.warn(
                f"RDKit failed to infer bond orders due to: {e}. Falling back to native method."
            )
            _infer_bond_orders_native(molecule)

    elif method == "native":
        _infer_bond_orders_native(molecule)

    else:
        msg = f"Method '{method}' for inferring bond orders is not recognized."
        err = ValueError
        if method == "rdkit" and not has_rdkit:
            msg += " (RDKit is not installed.)"
            err = RuntimeError
        raise err(msg)
