"""
Structure analysis and validation helpers.
"""

import numpy as np
from scipy.spatial.distance import cdist

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
    r_a = np.empty((len(residues_a)), dtype=object)
    r_b = np.empty((len(residues_b)), dtype=object)
    r_a[:] = residues_a
    r_b[:] = residues_b
    residues_a = r_a
    residues_b = r_b

    if coarse_precheck and len(residues_a) > 1 and len(residues_b) > 1:
        residue_coords_a = np.array([r.center_of_mass() for r in residues_a])
        residue_coords_b = np.array([r.center_of_mass() for r in residues_b])

        residue_dists = cdist(residue_coords_a, residue_coords_b)
        np.fill_diagonal(residue_dists, np.inf)
        residue_edge_mask = np.zeros(residue_dists.shape, dtype=bool)
        for i in range(len(residues_a)):
            for j in range(len(residues_b)):
                if residue_dists[i, j] < 12:
                    residue_edge_mask[i, j] = True
    else:
        residue_edge_mask = np.ones((len(residues_a), len(residues_b)), dtype=bool)

    for i in range(len(residues_a)):
        residue_a = residues_a[i]
        close_by_residues = residues_b[residue_edge_mask[i]]

        if ignore_hydrogens:
            atoms_a = [a for a in residue_a.get_atoms() if a.element != "H"]
            atoms_b = []
            for residue_b in close_by_residues:
                atoms_b.extend([a for a in residue_b.get_atoms() if a.element != "H"])
        else:
            atoms_a = list(residue_a.get_atoms())
            atoms_b = []
            for residue_b in close_by_residues:
                atoms_b.extend(residue_b.get_atoms())

        atoms_a = np.array(atoms_a, dtype=object)
        atoms_b = np.array(atoms_b, dtype=object)
        coords_a = np.array([a.get_coord() for a in atoms_a])
        coords_b = np.array([a.get_coord() for a in atoms_b])
        dists = cdist(coords_a, coords_b)
        np.fill_diagonal(dists, np.inf)

        xs, ys = np.where((0 < dists) * (dists < min_dist))
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
    Infer residues likely on the surface using SASA.
    """
    sasa = defaults.get_default_instance("bioSASA")
    sasa.compute(structure, level="R")

    sasa_values = np.array([residue.sasa for residue in structure.get_residues()])
    sasa_values = sasa_values / sasa_values.max() * 100

    if fraction is not None:
        cutoff = np.percentile(sasa_values, 100 - fraction * 100)

    surface_residues = [
        residue
        for residue, sasa in zip(structure.get_residues(), sasa_values)
        if sasa > cutoff
    ]
    return surface_residues
