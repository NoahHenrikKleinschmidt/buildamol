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
from scipy.spatial.distance import pdist

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

        # flat arrays to compute the unit gaussians with vectorized reductions
        unit_ids = sorted(self.rotation_units.keys())
        self._unit_order = np.concatenate([self.rotation_units[k] for k in unit_ids])
        self._unit_sizes = np.array([len(self.rotation_units[k]) for k in unit_ids])
        self._unit_starts = np.concatenate([[0], np.cumsum(self._unit_sizes)[:-1]])

        # only pairs of atoms in different rotation units can change their distance
        unit_of = np.empty(self.n_nodes, dtype=int)
        for k in unit_ids:
            unit_of[self.rotation_units[k]] = k
        i, j = np.triu_indices(self.n_nodes, k=1)
        self._clash_pairs = unit_of[i] != unit_of[j]

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

        Each unit is described by its center and its spread (its radius of gyration). The penalty increases as the distance between two units
        becomes smaller than the sum of their spreads.
        """
        sizes = self._unit_sizes
        coords = state[self._unit_order]
        means = np.add.reduceat(coords, self._unit_starts, axis=0) / sizes[:, None]

        # trace of the (sample) covariance of each unit
        squares = np.add.reduceat(np.sum(coords**2, axis=1), self._unit_starts, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            traces = (squares - sizes * np.sum(means**2, axis=1)) / (sizes - 1)
        # single atoms get the default spread, all others are regularized to avoid zero spreads
        traces = np.where(
            sizes > 1,
            np.clip(traces, 0, None) + 3 * self.gaussian_spread * 0.1,
            3 * self.gaussian_spread**2,
        )
        spreads = np.sqrt(traces)

        i, j = np.triu_indices(len(sizes), k=1)
        center_dists = np.linalg.norm(means[i] - means[j], axis=1)
        overlap = np.clip(spreads[i] + spreads[j] - center_dists, 0, None)

        # Skip far pairs if cutoff enabled
        if self.use_overlap_cutoff:
            overlap[center_dists > self.overlap_cutoff_distance] = 0

        return float(np.sum(overlap**2))

    def _compute_clash_penalty(self, state):
        """
        Compute direct steric clash penalty: sum of (clash_dist - actual_dist)^2
        for all atom pairs closer than clash_distance that can move relative to each other.
        """
        if self.n_nodes < 2:
            return 0.0

        distances = pdist(state)[self._clash_pairs]
        violations = self.clash_distance - distances[distances < self.clash_distance]
        return float(np.sum(violations**2))

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
