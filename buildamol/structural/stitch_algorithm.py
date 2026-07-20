"""
Thread-local stitching algorithm configuration.

The stitcher optimizes bond rotation angles to minimize steric clashes after
two molecules are joined. By default it uses particle swarm optimization (PSO).
This module lets you switch algorithms globally or per-block:

    # Change the global default:
    import buildamol as bam
    bam.set_stitching_algorithm("grid", n_points=72)

    # Override for one operation only (thread-safe):
    with bam.stitching_algorithm("scipy"):
        mol.attach(fragment, link)

Available algorithms
--------------------
"pso"
    Particle swarm optimization (default).  ``steps`` controls ``max_steps``,
    ``n_particles`` defaults to 5.
"scipy"
    Gradient-based minimisation via ``scipy.optimize.minimize`` (L-BFGS-B).
    ``steps`` becomes ``maxiter``.  Fast per evaluation but may settle in a
    local minimum; best combined with a good initial guess.
"genetic"
    Genetic algorithm.  ``steps`` is ``max_generations``.  Slower but robust
    on higher-dimensional angle spaces.
"anneal"
    Simulated annealing.  ``steps`` is the step budget.
"grid"
    Exhaustive grid search over all angle combinations.  Deterministic and
    globally optimal on the grid.  ``n_points`` (default 36) controls angular
    resolution (36 → 10°, 72 → 5°); cost scales as n_points^n_bonds.
    The ``steps`` argument is ignored.
"""

import contextlib
import threading

import numpy as np

_thread_local = threading.local()
_default_algorithm = "pso"
_default_kwargs: dict = {}

_ALGORITHMS = ("pso", "scipy", "genetic", "anneal", "grid", "none")


def get_stitching_algorithm():
    """Return *(algorithm_name, kwargs)* for the current thread."""
    algo = getattr(_thread_local, "algorithm", None)
    if algo is not None:
        return algo, dict(getattr(_thread_local, "algorithm_kwargs", {}))
    return _default_algorithm, dict(_default_kwargs)


def set_stitching_algorithm(algorithm: str, **kwargs):
    """Set the global default stitching optimization algorithm.

    Parameters
    ----------
    algorithm : str
        One of ``'pso'``, ``'scipy'``, ``'genetic'``, ``'anneal'``, ``'grid'``.
    **kwargs
        Algorithm-specific defaults, e.g. ``n_particles=10`` for pso or
        ``n_points=72`` for grid.  Per-call kwargs passed to
        :func:`run_stitching_optimization` override these.
    """
    global _default_algorithm, _default_kwargs
    _validate_algorithm(algorithm)
    _default_algorithm = algorithm
    _default_kwargs = dict(kwargs)


@contextlib.contextmanager
def stitching_algorithm(algorithm: str, **kwargs):
    """Context manager that overrides the stitching algorithm for this thread.

    Parameters
    ----------
    algorithm : str
        One of ``'pso'``, ``'scipy'``, ``'genetic'``, ``'anneal'``, ``'grid'``.
    **kwargs
        Algorithm-specific keyword arguments; override global defaults inside
        the block.

    Examples
    --------
    >>> with bam.stitching_algorithm("grid", n_points=72):
    ...     product = mol.attach(fragment, link)
    """
    _validate_algorithm(algorithm)
    old_algo = getattr(_thread_local, "algorithm", None)
    old_kwargs = getattr(_thread_local, "algorithm_kwargs", {})
    _thread_local.algorithm = algorithm
    _thread_local.algorithm_kwargs = dict(kwargs)
    try:
        yield
    finally:
        _thread_local.algorithm = old_algo
        _thread_local.algorithm_kwargs = old_kwargs


def _validate_algorithm(algorithm: str):
    if algorithm == "swarm":
        algorithm = "pso"  # backwards compatibility
    if algorithm not in _ALGORITHMS:
        raise ValueError(
            f"Unknown stitching algorithm {algorithm!r}. "
            f"Choose one of {_ALGORITHMS}."
        )


def grid_optimize(env, n_points: int = 36, **_ignored):
    """Exhaustive grid search over bond rotation angles.

    Samples *n_points* evenly-spaced angles in ``[-π, π)`` for each rotatable
    bond and evaluates every combination, returning the global grid minimum.

    Parameters
    ----------
    env : Rotatron
        The optimisation environment produced by ``Stitcher._optimize``.
    n_points : int
        Grid points per angle dimension.  36 → 10° resolution, 72 → 5°.
        For two rotatable bonds the cost is ``n_points²`` environment steps.
    """
    from itertools import product as _product

    n_actions = env.action_space.shape[0]
    angles = np.linspace(-np.pi, np.pi, n_points, endpoint=False)

    best_score = np.inf
    best_angles = np.zeros(n_actions)

    for combo in _product(angles, repeat=n_actions):
        action = np.array(combo)
        score = env.step(action)[1]
        env.reset()
        if score < best_score:
            best_score = score
            best_angles = action.copy()

    return best_angles, best_score


def run_stitching_optimization(env, steps: int, **caller_kwargs):
    """Dispatch to the configured stitching algorithm.

    Reads the thread-local (or global) algorithm setting, merges it with
    *caller_kwargs* (caller wins on conflicts), and runs the chosen solver.

    Parameters
    ----------
    env : Rotatron
        The DistanceRotatron produced by ``Stitcher._optimize``.
    steps : int
        Step budget passed to the algorithm (ignored by ``"grid"``).
    **caller_kwargs
        Extra kwargs forwarded to the solver; override global defaults.
    """
    import buildamol.optimizers as optimizers

    algorithm, global_kwargs = get_stitching_algorithm()
    kwargs = {**global_kwargs, **caller_kwargs}

    if algorithm == "pso":
        return optimizers.swarm_optimize(
            env,
            n_particles=kwargs.pop("n_particles", 5),
            max_steps=int(steps),
            **kwargs,
        )
    if algorithm == "scipy":
        return optimizers.scipy_optimize(env, int(steps), **kwargs)
    if algorithm == "genetic":
        return optimizers.genetic_optimize(env, int(steps), **kwargs)
    if algorithm == "anneal":
        return optimizers.anneal_optimize(env, int(steps), **kwargs)
    if algorithm == "grid":
        return grid_optimize(env, **kwargs)

    raise ValueError(f"Unknown stitching algorithm {algorithm!r}.")
