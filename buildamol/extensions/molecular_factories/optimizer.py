"""
Optimizer — wraps a BuildAMol pipeline with numeric parameter injection and
runs a search strategy (random, PSO, or scipy) to maximise a user-supplied
scoring function.
"""

import numpy as np

from .base import ChainableBlock
from buildamol.utils.auxiliary import progress_bar


class Optimizer:
    """
    Optimise a stochastic pipeline by searching its discrete/continuous
    parameter space.

    Parameters
    ----------
    pipeline : ChainableBlock
        A fully assembled pipeline (chain of blocks) whose stochastic
        blocks expose ``n_params`` / ``param_bounds``.
    scoring_fn : callable
        ``scoring_fn(molecule) -> float`` — higher is better.

    Examples
    --------
    >>> opt = Optimizer(pipeline, scoring_fn=lambda mol: mol.num_heavy_atoms)
    >>> opt.run(strategy="pso", steps=50, population=20)
    >>> best = opt.top(5)
    """

    def __init__(self, pipeline: ChainableBlock, scoring_fn, n_workers: int = 1):
        self.pipeline = pipeline
        self.scoring_fn = scoring_fn
        self.n_workers = n_workers
        self._results: list[tuple[float, object]] = []  # (score, molecule)
        self._executor = None

    # ── public interface ──────────────────────────────────────────────────────

    def run(
        self,
        strategy: str = "random",
        steps: int = 100,
        population: int = 50,
        verbose: bool = True,
        n_workers: int = None,
        **kwargs,
    ):
        """
        Run the optimisation.

        Parameters
        ----------
        strategy : {"random", "pso"/"swarm", "scipy"}
        steps : int
            Number of iterations / evaluations.
        population : int
            Swarm size (PSO only).
        verbose : bool
            Show a progress bar (default ``True``).
        n_workers : int, optional
            Number of parallel threads for pipeline execution and scoring.
            Defaults to the value set at construction time.  Threading is used
            (not multiprocessing) so the scoring function must be thread-safe
            and ideally releases the GIL (RDKit, PyTorch, admet-ai all do).
            The param-injection mechanism is already thread-safe via
            ``threading.local()``, so no special setup is needed.
        **kwargs
            Extra keyword arguments forwarded to the strategy implementation.
        """
        dispatch = {
            "random": self._run_random,
            "pso": self._run_pso,
            "swarm": self._run_pso,
            "scipy": self._run_scipy,
        }
        if strategy not in dispatch:
            raise ValueError(
                f"Unknown strategy {strategy!r}. Choose from {list(dispatch)}."
            )
        n_workers = n_workers if n_workers is not None else self.n_workers
        if n_workers > 1:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=n_workers) as exe:
                self._executor = exe
                try:
                    dispatch[strategy](steps=steps, population=population, verbose=verbose, **kwargs)
                finally:
                    self._executor = None
        else:
            dispatch[strategy](steps=steps, population=population, verbose=verbose, **kwargs)
        self._results.sort(key=lambda x: x[0], reverse=True)

    def top(self, n: int = 10) -> list:
        """Return the top *n* molecules (best scoring first)."""
        return [mol for _, mol in self._results[:n]]

    def scores(self, n: int = None) -> list[tuple[float, object]]:
        """Return ``[(score, molecule), ...]`` sorted best-first."""
        return self._results[:n] if n is not None else list(self._results)

    # ── core evaluation ───────────────────────────────────────────────────────

    def _evaluate(self, params) -> tuple[float, object] | None:
        """Inject *params*, run the pipeline, score the result."""
        ChainableBlock._inject_params(params)
        try:
            ctx = self.pipeline()
            mol = ctx.molecule
            score = float(self.scoring_fn(mol))
            return score, mol
        except Exception:
            return None
        finally:
            ChainableBlock._clear_params()

    def _evaluate_batch(self, param_list: list) -> list:
        """Evaluate multiple param vectors, in parallel when ``n_workers > 1``.

        Each thread gets its own ``threading.local()`` param iterator so there
        is no cross-thread interference.  Results that failed are excluded from
        the returned list (no ``None`` entries).
        """
        if self._executor is not None:
            raw = list(self._executor.map(self._evaluate, param_list))
        else:
            raw = [self._evaluate(p) for p in param_list]
        return [r for r in raw if r is not None]

    # ── strategies ────────────────────────────────────────────────────────────

    def _run_random(self, steps: int, population: int = 50, verbose: bool = True, **_):
        """Uniform random search over the parameter bounds."""
        bounds = self.pipeline.param_bounds

        # Generate all param vectors up front so the batch can be submitted to
        # the thread pool in one shot (no dependency between evaluations).
        if not bounds:
            param_list = [[] for _ in range(steps)]
        else:
            lo = np.array([b[0] for b in bounds], dtype=float)
            hi = np.array([b[1] for b in bounds], dtype=float)
            param_list = [lo + np.random.random(len(bounds)) * (hi - lo) for _ in range(steps)]

        bar = progress_bar(total=steps, desc="random", disable=not verbose)
        best = float("-inf")
        try:
            for result in self._evaluate_batch(param_list):
                self._results.append(result)
                bar.update(1)
                if result[0] > best:
                    best = result[0]
                    bar.set_postfix(best=f"{best:.4f}")
        finally:
            bar.close()

    def _run_pso(
        self,
        steps: int,
        population: int = 50,
        w: float = 0.72,
        c1: float = 1.49,
        c2: float = 1.49,
        verbose: bool = True,
        **_,
    ):
        """
        Particle Swarm Optimisation (Clerc-Kennedy coefficients).

        Because the pipeline is discrete (integer indices), position values are
        rounded to int when consumed by the blocks.  Velocities and positions are
        maintained as floats for smooth gradient behaviour.
        """
        bounds = self.pipeline.param_bounds
        if not bounds:
            self._run_random(steps=steps * population, verbose=verbose)
            return

        n = len(bounds)
        lo = np.array([b[0] for b in bounds], dtype=float)
        hi = np.array([b[1] for b in bounds], dtype=float)

        pos = lo + np.random.random((population, n)) * (hi - lo)
        vel = np.zeros((population, n))
        pbest_pos = pos.copy()
        pbest_score = np.full(population, -np.inf)

        gbest_pos = pos[0].copy()
        gbest_score = -np.inf

        bar = progress_bar(range(steps), desc="pso", disable=not verbose)
        for _ in bar:
            # All population members are independent within a step — evaluate in parallel.
            # executor.map preserves submission order so particle index i stays aligned.
            params = [pos[i] for i in range(population)]
            if self._executor is not None:
                raw = list(self._executor.map(self._evaluate, params))
            else:
                raw = [self._evaluate(p) for p in params]

            for i, result in enumerate(raw):
                if result is None:
                    continue
                score, mol = result
                self._results.append((score, mol))
                if score > pbest_score[i]:
                    pbest_score[i] = score
                    pbest_pos[i] = pos[i].copy()
                if score > gbest_score:
                    gbest_score = score
                    gbest_pos = pos[i].copy()

            bar.set_postfix(best=f"{gbest_score:.4f}")
            r1 = np.random.random((population, n))
            r2 = np.random.random((population, n))
            vel = w * vel + c1 * r1 * (pbest_pos - pos) + c2 * r2 * (gbest_pos - pos)
            pos = np.clip(pos + vel, lo, hi)

    def _run_scipy(
        self,
        steps: int,
        population: int = 50,
        scipy_method: str = "Nelder-Mead",
        verbose: bool = True,
        **_,
    ):
        """Use a scipy minimiser (negated score = minimise)."""
        try:
            from scipy.optimize import minimize
        except ImportError as exc:
            raise ImportError(
                "scipy is required for strategy='scipy'. Install with: pip install scipy"
            ) from exc

        bounds = self.pipeline.param_bounds
        if not bounds:
            self._run_random(steps=steps, verbose=verbose)
            return

        lo = np.array([b[0] for b in bounds], dtype=float)
        hi = np.array([b[1] for b in bounds], dtype=float)

        bar = progress_bar(total=steps, desc="scipy", disable=not verbose)
        best = float("-inf")

        def objective(params):
            nonlocal best
            result = self._evaluate(np.clip(params, lo, hi))
            if result is None:
                return 0.0
            score, mol = result
            self._results.append((score, mol))
            bar.update(1)
            if score > best:
                best = score
                bar.set_postfix(best=f"{best:.4f}")
            return -score  # minimise negated score

        x0 = lo + np.random.random(len(bounds)) * (hi - lo)
        try:
            minimize(objective, x0, method=scipy_method, options={"maxiter": steps})
        finally:
            bar.close()
