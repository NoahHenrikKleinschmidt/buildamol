"""
Bond and connectivity inference helpers.
"""

import warnings

import numpy as np
from Bio.PDB import NeighborSearch

import buildamol.resources as resources
import buildamol.utils.defaults as defaults

from .constants import _MIN_BOND_LENGTH, _bond_cutoff_vdw, _max_search_radius_vdw


def infer_residue_connections(
    structure, bond_length: float = None, triplet: bool = False
):
    """
    Infer connectivity between residues from inter-atom distances.
    """
    if not triplet:
        if bond_length is None:
            min_length = _MIN_BOND_LENGTH
            max_length = None
            use_vdw = True
        elif isinstance(bond_length, (int, float)):
            min_length, max_length = defaults.DEFAULT_BOND_LENGTH / 2, bond_length
            use_vdw = False
        else:
            min_length, max_length = bond_length
            use_vdw = False

        bonds = []
        _seen_residues = set()
        for residue1 in structure.get_residues():
            for residue2 in structure.get_residues():
                if residue1 == residue2:
                    continue
                elif residue2 in _seen_residues:
                    continue

                atoms = list(i for i in residue1.get_atoms() if i.element != "H")
                atoms.extend(i for i in residue2.get_atoms() if i.element != "H")

                search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
                if search_radius == 0:
                    continue

                _neighbors = NeighborSearch(atoms).search_all(radius=search_radius)

                for atom1, atom2 in _neighbors:
                    if atom1.get_parent() == atom2.get_parent():
                        continue

                    dist = np.linalg.norm(atom1.coord - atom2.coord)
                    cutoff = _bond_cutoff_vdw(atom1, atom2) if use_vdw else max_length

                    if min_length < dist <= cutoff:
                        bonds.append((atom1, atom2))

            _seen_residues.add(residue1)
    else:
        connections = infer_residue_connections(
            structure, bond_length=bond_length, triplet=False
        )
        base_bonds = infer_bonds(
            structure, bond_length=bond_length, restrict_residues=True
        )
        _additional_bonds = []
        for atom1, atom2 in connections:
            _new = (bond for bond in base_bonds if atom1 in bond)
            _additional_bonds.extend(_new)
        bonds = connections + _additional_bonds

    return bonds


def infer_bonds(structure, bond_length: float = None, restrict_residues: bool = True):
    """
    Generate connectivity by inferring bonds from inter-atom distances.
    """
    if bond_length is None:
        min_length = _MIN_BOND_LENGTH
        max_length = None
        use_vdw = True
    elif isinstance(bond_length, (int, float)):
        min_length, max_length = (defaults.DEFAULT_BOND_LENGTH / 2, bond_length)
        use_vdw = False
    else:
        min_length, max_length = (
            defaults.DEFAULT_BOND_LENGTH / 2,
            defaults.DEFAULT_BOND_LENGTH,
        )
        use_vdw = False

    bonds = []
    if restrict_residues:
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
            if search_radius == 0:
                continue

            _neighbors = NeighborSearch(atoms).search_all(radius=search_radius)
            for atom1, atom2 in _neighbors:
                if atom1.element == "H" and atom2.element == "H":
                    continue

                dist = np.linalg.norm(atom1.coord - atom2.coord)
                cutoff = _bond_cutoff_vdw(atom1, atom2) if use_vdw else max_length

                if min_length < dist <= cutoff:
                    bonds.append((atom1, atom2))
    else:
        atoms = list(structure.get_atoms())
        search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
        if search_radius == 0:
            return []

        _neighbors = NeighborSearch(atoms).search_all(radius=search_radius)

        bonds = []
        for atom1, atom2 in _neighbors:
            if atom1.element == "H" and atom2.element == "H":
                continue

            dist = np.linalg.norm(atom1.coord - atom2.coord)
            cutoff = _bond_cutoff_vdw(atom1, atom2) if use_vdw else max_length

            if min_length < dist <= cutoff:
                bonds.append((atom1, atom2))

    return bonds


def _atom_from_residue(id, residue):
    return next((atom for atom in residue.get_atoms() if atom.id == id), None)


def apply_reference_bonds(structure, _compounds=None):
    """
    Apply residue-internal bonds according to loaded reference compounds.
    """
    if _compounds is None:
        _compounds = resources.get_default_compounds()

    if structure.level == "R":
        residue = structure

        if not _compounds.has_residue(residue.resname):
            warnings.warn(
                f"[ignoring] No reference residue found in Compounds for {residue.resname}!"
            )
            return []

        ref = _compounds.get(residue.resname)
        bonds = (
            (
                _atom_from_residue(bond.atom1.id, residue),
                _atom_from_residue(bond.atom2.id, residue),
                bond.order,
            )
            for bond in ref.get_bonds()
        )
        bonds = [bond for bond in bonds if bond[1] is not None and bond[0] is not None]
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_reference_bonds(residue, _compounds)
        bonds = [
            bond
            for bond in bonds
            if bond[0].id == structure.id or bond[1].id == structure.id
        ]
        return bonds

    else:
        bonds = []
        for residue in structure.get_residues():
            bonds.extend(apply_reference_bonds(residue, _compounds))
        return bonds


def atoms_in_area(structure, center, radius):
    """
    Yield all atoms in a given spherical area around a center point.
    """
    for atom in structure.get_atoms():
        if np.linalg.norm(atom.coord - center) < radius:
            yield atom
