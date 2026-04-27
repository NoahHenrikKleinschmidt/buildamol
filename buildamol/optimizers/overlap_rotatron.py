"""
The OverlapRotatron environment optimizes molecular conformations by treating
rotation units (groups of atoms that move together) as Gaussian-distributed entities
and minimizing both their overlap and pairwise clashes.

The evaluation score is:
  score = sum of overlap penalties between all pairs of rotation units
        + sum of clash penalties for atom pairs below the clash distance threshold

Both terms are minimized, so lower score = better conformation.
"""

import gymnasium as gym
import numpy as np
from scipy.spatial.distance import pdist, squareform

import buildamol.optimizers.base_rotatron as Rotatron
import buildamol.graphs.base_graph as base_graph

__all__ = ["OverlapRotatron", "MVN", "jensen_shannon_overlap"]


class OverlapRotatron(Rotatron.Rotatron):
    """
    Gaussian-based molecular conformation optimizer using rotation units.

    Treats each rotation unit (group of atoms moving together) as a Gaussian
    distribution and minimizes both pairwise overlap and hard steric clashes.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The molecular graph to optimize
    rotatable_edges : list
        Edges that can rotate during optimization
    gaussian_spread : float
        Spread parameter for Gaussian fitting (multiplier on covariance)
    clash_distance : float
        Distance threshold below which atoms incur a clash penalty
    clash_weight : float
        Weight factor for clash penalty term (relative to overlap penalty)
    overlap_cutoff_distance : float
        Pairs of rotation units further apart than this are not scored for overlap
        (set ≤ 0 to disable cutoff)
    crop_nodes_further_than : float
        If > 0, crop nodes far from rotatable edges to reduce computation
    n_processes : int
        Number of processes for edge mask computation
    bounds : tuple
        Rotation angle bounds (min, max)
    **kwargs
        Additional backend arguments
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        gaussian_spread: float = 1.5,
        clash_distance: float = 1.2,
        clash_weight: float = 1.0,
        overlap_cutoff_distance: float = -1.0,
        crop_nodes_further_than: float = -1,
        n_processes: int = 1,
        bounds: tuple = (-np.pi, np.pi),
        **kwargs,
    ):
        # Store hyperparameters for later reference
        self.hyperparameters = {
            "gaussian_spread": gaussian_spread,
            "clash_distance": clash_distance,
            "clash_weight": clash_weight,
            "overlap_cutoff_distance": overlap_cutoff_distance,
            "crop_nodes_further_than": crop_nodes_further_than,
            "n_processes": n_processes,
            "bounds": bounds,
            **kwargs,
        }

        self.gaussian_spread = gaussian_spread
        self.clash_distance = clash_distance
        self.clash_weight = clash_weight
        self.use_overlap_cutoff = overlap_cutoff_distance > 0
        self.overlap_cutoff_distance = overlap_cutoff_distance

        self._bounds_tuple = bounds
        self.crop_radius = crop_nodes_further_than

        # Get edges and set up base
        rotatable_edges_clean = self._get_rotatable_edges(graph, rotatable_edges)
        self.graph = graph
        self.rotatable_edges = rotatable_edges_clean
        self.n_nodes = len(self.graph.nodes)
        self.n_edges = len(self.rotatable_edges)

        # Optional node cropping
        if self.crop_radius > 0:
            graph, rotatable_edges_clean = self._setup_helpers_crop_faraway_nodes(
                self.crop_radius, graph, rotatable_edges_clean
            )

        self.action_space = gym.spaces.Box(
            low=bounds[0], high=bounds[1], shape=(len(rotatable_edges_clean),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self.graph.nodes), 3)
        )

        # Initialize base Rotatron
        Rotatron.Rotatron.__init__(
            self, graph, rotatable_edges_clean, n_processes=n_processes, **kwargs
        )

        # Build rotation units from edge masks
        self._build_rotation_units()

    def _build_rotation_units(self):
        """
        Identify which atoms move together (rotation units).

        Two atoms are in the same unit if rotating any edge affects them identically.
        This method keeps ALL atoms (including singletons) in the objective.
        """
        # Create a matrix: row = edge, col = atom, value = 1 if atom moves with this edge
        rotation_matrix = self.edge_masks.astype(int)  # (n_edges, n_atoms)

        # Find unique movement patterns
        # Two atoms are in the same unit if they have identical patterns
        unique_patterns = {}
        self.rotation_units = {}

        for atom_idx in range(self.n_nodes):
            pattern = tuple(rotation_matrix[:, atom_idx])
            if pattern not in unique_patterns:
                unit_id = len(self.rotation_units)
                unique_patterns[pattern] = unit_id
                self.rotation_units[unit_id] = []
            else:
                unit_id = unique_patterns[pattern]

            self.rotation_units[unit_id].append(atom_idx)

        # Convert lists to numpy arrays
        self.rotation_units = {
            k: np.array(v, dtype=int) for k, v in self.rotation_units.items()
        }

    def eval(self, state):
        """
        Evaluate conformation quality as sum of overlap + clash penalties.

        Parameters
        ----------
        state : np.ndarray, shape (n_atoms, 3)
            Atomic coordinates

        Returns
        -------
        float
            Total penalty score (lower = better)
        """
        score = 0.0

        # --- Overlap penalty ---
        score += self._compute_overlap_penalty(state)

        # --- Clash penalty ---
        score += self.clash_weight * self._compute_clash_penalty(state)

        return float(score)

    def _compute_overlap_penalty(self, state):
        """
        Compute pairwise Gaussian overlap penalty between all rotation units.

        Uses a simple spread-based measure: penalty increases as units get closer
        relative to their spread.
        """
        unit_ids = sorted(self.rotation_units.keys())
        gaussians = []

        for unit_id in unit_ids:
            atom_indices = self.rotation_units[unit_id]
            unit_coords = state[atom_indices]

            # Fit Gaussian: mean and covariance
            mean = np.mean(unit_coords, axis=0)

            if len(atom_indices) == 1:
                # Single atom: zero covariance, use gaussian_spread as default
                cov = np.eye(3) * (self.gaussian_spread**2)
            else:
                cov = np.cov(unit_coords.T)

                # Handle degenerate covariance (linear arrangement)
                if cov.ndim == 0:
                    cov = np.eye(3) * (self.gaussian_spread**2)
                else:
                    cov = np.atleast_2d(cov)
                    if cov.shape != (3, 3):
                        cov = np.eye(3) * (self.gaussian_spread**2)
                    else:
                        # Regularize to avoid singular matrices
                        cov = cov + np.eye(3) * (self.gaussian_spread * 0.1)

            gaussians.append((mean, cov))

        # Compute pairwise overlaps
        penalty = 0.0
        for i in range(len(gaussians)):
            for j in range(i + 1, len(gaussians)):
                mean_i, cov_i = gaussians[i]
                mean_j, cov_j = gaussians[j]

                # Distance between centers
                center_dist = np.linalg.norm(mean_i - mean_j)

                # Skip far pairs if cutoff enabled
                if (
                    self.use_overlap_cutoff
                    and center_dist > self.overlap_cutoff_distance
                ):
                    continue

                # Overlap penalty based on distance and spread
                # Smaller distance → higher penalty
                # Larger spread → higher penalty tolerance
                spread_metric = (np.trace(cov_i) + np.trace(cov_j)) / 6.0

                # Penalty increases as distance decreases relative to spread
                # If center_dist << spread, high penalty
                overlap_amount = max(0.0, spread_metric - center_dist)
                penalty += overlap_amount**2

        return penalty

    def _compute_clash_penalty(self, state):
        """
        Compute direct steric clash penalty: sum of (clash_dist - actual_dist)^2
        for all atom pairs closer than clash_distance.
        """
        penalty = 0.0

        # Pairwise distances
        if self.n_nodes < 2:
            return penalty

        all_pairs = pdist(state)
        distances = squareform(all_pairs)

        # Find pairs below threshold
        clashing = distances < self.clash_distance
        np.fill_diagonal(clashing, False)  # Ignore self-distances

        if np.any(clashing):
            # Squared penalty for violations
            violations = self.clash_distance - distances[clashing]
            penalty = np.sum(violations**2)

        return penalty

    def copy(self):
        """Make a deep copy of the environment"""
        from copy import deepcopy

        return deepcopy(self)

    def reset(self, *args, **kwargs):
        """Reset to initial coordinates"""
        super().reset(*args, **kwargs)


# ============================================================================
# Utility functions for backward compatibility and advanced features
# ============================================================================

_MULTIVARIATE_NORMAL = None


def _get_multivariate_normal():
    """Lazy load scipy's multivariate_normal"""
    global _MULTIVARIATE_NORMAL
    if _MULTIVARIATE_NORMAL is None:
        from scipy.stats import multivariate_normal

        _MULTIVARIATE_NORMAL = multivariate_normal
    return _MULTIVARIATE_NORMAL


def MVN(points, spread: float = 1.0):
    """
    Compute a multi-variate normal distribution for a given set of points.

    Parameters
    ----------
    points : np.ndarray
        The points to compute the mean and covariance matrix for.
    spread : float
        Spread multiplier for covariance

    Returns
    -------
    mvn : scipy.stats.multivariate_normal
        The multi-variate normal distribution for the points.
    """
    return _get_multivariate_normal()(
        mean=np.mean(points, axis=0),
        cov=spread * np.cov(points, rowvar=False),
        allow_singular=True,
    )


def _kl_divergence(p, q):
    """
    Compute the Kullback-Leibler divergence between two distributions.

    Parameters
    ----------
    p : float or array
        First probability value(s)
    q : float or array
        Second probability value(s)

    Returns
    -------
    float
        KL divergence
    """
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))


def jensen_shannon_overlap(mvn1, mvn2):
    """
    Compute the overlap between two gaussians using the Jensen-Shannon divergence.

    Parameters
    ----------
    mvn1, mvn2 : scipy.stats.multivariate_normal
        The two gaussians to compute the overlap for.

    Returns
    -------
    overlap : float
        The overlap between the two gaussians (lower = more separated).
    """
    # Create Multivariate Normal distributions for each Gaussian
    center1 = mvn1.mean
    center2 = mvn2.mean

    pdf1_center1 = mvn1.pdf(center1)
    pdf2_center2 = mvn2.pdf(center2)
    pdf1_center2 = mvn1.pdf(center2)
    pdf2_center1 = mvn2.pdf(center1)

    mean_center1 = (pdf1_center1 + pdf2_center1) / 2
    mean_center2 = (pdf1_center2 + pdf2_center2) / 2

    dist = _kl_divergence(pdf1_center1, mean_center1) + _kl_divergence(
        pdf2_center2, mean_center2
    )
    dist *= 0.5

    return -dist
