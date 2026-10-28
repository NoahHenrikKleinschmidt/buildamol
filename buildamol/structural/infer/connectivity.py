"""
Bond and connectivity inference helpers.
"""

import warnings

import numpy as np

import buildamol.resources as resources
import buildamol.utils.defaults as defaults

from .constants import _MIN_BOND_LENGTH, _bond_cutoff_vdw, _max_search_radius_vdw


def _search_pairs(atoms, radius):
    """
    Return (atom_i, atom_j, dist_ij) triples for all pairs within *radius*.
    Distances are included so callers avoid recomputing them.
    """
    if len(atoms) < 2:
        return []
    coords = np.array([a.coord for a in atoms])
    diff = coords[:, np.newaxis] - coords[np.newaxis, :]
    dists = np.linalg.norm(diff, axis=-1)
    i, j = np.where((dists < radius) & (dists > 0))
    mask = i < j
    ii, jj = i[mask], j[mask]
    return [
        (atoms[int(a)], atoms[int(b)], float(dists[a, b]))
        for a, b in zip(ii, jj)
    ]


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

        # Collect all non-H atoms once and build a global KDTree, which is
        # O(N log N + M) vs the previous O(N_residues² × N_per_residue²) loop.
        atoms = [a for a in structure.get_atoms() if a.element != "H"]
        if len(atoms) < 2:
            return []

        search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
        if not search_radius:
            return []

        try:
            from scipy.spatial import KDTree
            coords = np.array([a.coord for a in atoms])
            tree = KDTree(coords)
            pair_indices = tree.query_pairs(r=search_radius, output_type="ndarray")
        except Exception:
            # Scipy unavailable or error — fall back to pairwise numpy broadcast.
            # Only feasible for small structures.
            pair_indices = None

        bonds = []
        _cutoff_cache = {}  # (elem1, elem2) → float

        if pair_indices is not None and len(pair_indices):
            for a_idx, b_idx in pair_indices:
                atom1 = atoms[int(a_idx)]
                atom2 = atoms[int(b_idx)]
                if atom1.get_parent() is atom2.get_parent():
                    continue

                dist = float(np.linalg.norm(coords[a_idx] - coords[b_idx]))

                if use_vdw:
                    key = (atom1.element, atom2.element)
                    if key not in _cutoff_cache:
                        _cutoff_cache[key] = _bond_cutoff_vdw(atom1, atom2)
                        _cutoff_cache[(key[1], key[0])] = _cutoff_cache[key]
                    cutoff = _cutoff_cache[key]
                else:
                    cutoff = max_length

                if min_length < dist <= cutoff:
                    bonds.append((atom1, atom2))

        elif pair_indices is None:
            # Fallback: pairwise broadcast per residue pair (original logic).
            _seen_residues = set()
            for residue1 in structure.get_residues():
                for residue2 in structure.get_residues():
                    if residue1 == residue2:
                        continue
                    elif residue2 in _seen_residues:
                        continue

                    res_atoms = [a for a in residue1.get_atoms() if a.element != "H"]
                    res_atoms.extend(
                        a for a in residue2.get_atoms() if a.element != "H"
                    )

                    sr = _max_search_radius_vdw(res_atoms) if use_vdw else max_length
                    if not sr:
                        continue

                    for atom1, atom2, dist in _search_pairs(res_atoms, sr):
                        if atom1.get_parent() is atom2.get_parent():
                            continue
                        key = (atom1.element, atom2.element)
                        if key not in _cutoff_cache:
                            _cutoff_cache[key] = _bond_cutoff_vdw(atom1, atom2)
                            _cutoff_cache[(key[1], key[0])] = _cutoff_cache[key]
                        cutoff = _cutoff_cache[key] if use_vdw else max_length
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

    _cutoff_cache = {}  # (elem1, elem2) → float

    bonds = []
    if restrict_residues:
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
            if search_radius == 0:
                continue

            for atom1, atom2, dist in _search_pairs(atoms, search_radius):
                if atom1.element == "H" and atom2.element == "H":
                    continue

                if use_vdw:
                    key = (atom1.element, atom2.element)
                    if key not in _cutoff_cache:
                        _cutoff_cache[key] = _bond_cutoff_vdw(atom1, atom2)
                        _cutoff_cache[(key[1], key[0])] = _cutoff_cache[key]
                    cutoff = _cutoff_cache[key]
                else:
                    cutoff = max_length

                if min_length < dist <= cutoff:
                    bonds.append((atom1, atom2))
    else:
        atoms = list(structure.get_atoms())
        search_radius = _max_search_radius_vdw(atoms) if use_vdw else max_length
        if search_radius == 0:
            return []

        for atom1, atom2, dist in _search_pairs(atoms, search_radius):
            if atom1.element == "H" and atom2.element == "H":
                continue

            if use_vdw:
                key = (atom1.element, atom2.element)
                if key not in _cutoff_cache:
                    _cutoff_cache[key] = _bond_cutoff_vdw(atom1, atom2)
                    _cutoff_cache[(key[1], key[0])] = _cutoff_cache[key]
                cutoff = _cutoff_cache[key]
            else:
                cutoff = max_length

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
