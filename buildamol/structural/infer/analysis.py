"""
Structure analysis and validation helpers.
"""

import numpy as np


def cdist(a, b):
    """Euclidean pairwise distances between rows of a and b.
    Direct subtraction formula guarantees exact 0 for identical points."""
    diff = a[:, np.newaxis, :] - b[np.newaxis, :, :]
    return np.sqrt((diff ** 2).sum(axis=-1))

import buildamol.structural.base as base
import buildamol.utils.defaults as defaults


def vet_structure(
    molecule, clash_range: tuple = (0.6, 2.7), angle_range: tuple = (90, 180)
) -> bool:
    """
    Check for basic structure integrity.
    """
    for a, b in molecule.get_bonds():
        d = base.compute_distance(a, b)
        if not clash_range[0] <= d <= clash_range[1]:
            return False
    for angle in molecule.compute_angles().values():
        if not angle_range[0] <= angle <= angle_range[1]:
            return False
    for a in molecule.get_atoms():
        for b in molecule.get_atoms():
            if a is b:
                continue
            dist = base.compute_distance(a, b)
            if dist <= clash_range[0]:
                return False
    return True


def find_clashes_between(
    mol_a,
    mol_b,
    min_dist: float = 1.0,
    ignore_hydrogens: bool = False,
    coarse_precheck: bool = True,
):
    """
    Find all clashing atoms between two sets of atoms.
    """
    residues_a = list(mol_a.get_residues())
    residues_b = list(mol_b.get_residues())
    _arr_a = np.empty(len(residues_a), dtype=object)
    _arr_a[:] = residues_a
    residues_a = _arr_a
    _arr_b = np.empty(len(residues_b), dtype=object)
    _arr_b[:] = residues_b
    residues_b = _arr_b

    if ignore_hydrogens:
        residue_atoms_a = [
            np.asarray(
                [a for a in residue.get_atoms() if a.element != "H"], dtype=object
            )
            for residue in residues_a
        ]
        residue_atoms_b = [
            np.asarray(
                [a for a in residue.get_atoms() if a.element != "H"], dtype=object
            )
            for residue in residues_b
        ]
    else:
        residue_atoms_a = [
            np.asarray(list(residue.get_atoms()), dtype=object)
            for residue in residues_a
        ]
        residue_atoms_b = [
            np.asarray(list(residue.get_atoms()), dtype=object)
            for residue in residues_b
        ]

    residue_coords_a = [
        (
            np.asarray([a.get_coord() for a in atoms], dtype=float)
            if len(atoms)
            else np.empty((0, 3), dtype=float)
        )
        for atoms in residue_atoms_a
    ]
    residue_coords_b = [
        (
            np.asarray([a.get_coord() for a in atoms], dtype=float)
            if len(atoms)
            else np.empty((0, 3), dtype=float)
        )
        for atoms in residue_atoms_b
    ]

    if coarse_precheck and len(residues_a) > 1 and len(residues_b) > 1:
        residue_centers_a = np.asarray(
            [r.center_of_mass() for r in residues_a], dtype=float
        )
        residue_centers_b = np.asarray(
            [r.center_of_mass() for r in residues_b], dtype=float
        )
        residue_dists = cdist(residue_centers_a, residue_centers_b)
        residue_edge_mask = residue_dists < 12
    else:
        residue_edge_mask = np.ones((len(residues_a), len(residues_b)), dtype=bool)

    for i in range(len(residues_a)):
        atoms_a = residue_atoms_a[i]
        coords_a = residue_coords_a[i]
        if len(atoms_a) == 0:
            continue

        close_by_idx = np.where(residue_edge_mask[i])[0]
        if close_by_idx.size == 0:
            continue

        atoms_b_parts = [
            residue_atoms_b[j] for j in close_by_idx if len(residue_atoms_b[j])
        ]
        if not atoms_b_parts:
            continue
        coords_b_parts = [
            residue_coords_b[j] for j in close_by_idx if len(residue_coords_b[j])
        ]

        atoms_b = np.concatenate(atoms_b_parts)
        coords_b = np.concatenate(coords_b_parts, axis=0)

        dists = cdist(coords_a, coords_b)
        if mol_a is mol_b and dists.shape[0] == dists.shape[1]:
            np.fill_diagonal(dists, np.inf)

        xs, ys = np.where((0 < dists) & (dists < min_dist))
        for x, y in zip(xs, ys):
            yield atoms_a[x], atoms_b[y]


def sample_atoms_around_reference(
    reference_coord: np.ndarray,
    candidates: np.ndarray,
    num_samples: int,
    max_radius: float = 10.0,
):
    """
    Sample atoms around a reference coordinate.
    """
    coordinates = np.array([i.coord for i in candidates])

    phi = np.linspace(0, 2 * np.pi, num_samples)
    theta = np.linspace(0, np.pi, num_samples)

    phi_grid, theta_grid = np.meshgrid(phi, theta)

    x = max_radius * np.sin(theta_grid) * np.cos(phi_grid) + reference_coord[0]
    y = max_radius * np.sin(theta_grid) * np.sin(phi_grid) + reference_coord[1]
    z = max_radius * np.cos(theta_grid) + reference_coord[2]

    sampled_coordinates = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

    closest_indices = np.argmin(
        np.linalg.norm(sampled_coordinates[:, None] - coordinates, axis=2),
        axis=1,
    )
    samples = candidates[closest_indices]

    return samples


def compute_residue_radius(residue):
    """
    Compute the radius of a residue from center-of-mass to furthest atom.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    radius = np.max(distances)
    return radius


def compute_outlier_atoms(residue, f: float = 1.5):
    """
    Compute which atoms of a residue are especially far from center-of-mass.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    p75 = np.percentile(distances, 75)
    f = f * p75

    outlier_atoms = [atom for atom, distance in zip(atoms, distances) if distance > f]
    return outlier_atoms


def infer_surface_residues(structure, cutoff: int = 75, fraction: float = None):
    """
    .. deprecated::
        SASA-based surface detection has been removed. This function is no longer available.
    """
    raise NotImplementedError(
        "infer_surface_residues has been removed (SASA dependency dropped). "
        "Use geometry-based surface detection or an external SASA library."
    )
