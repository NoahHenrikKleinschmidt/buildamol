"""
SOOptimizer — single-objective variant of the BuildAMol Optimizer.

Accepts a scalar scoring function ``scoring_fn(mol) -> float`` and searches
the pipeline's parameter space using random search, PSO, or a scipy minimiser.
"""

import numpy as np

from .base import ChainableBlock
from buildamol.utils.auxiliary import progress_bar


class SOOptimizer:
    """
    Optimise a stochastic pipeline by searching its discrete/continuous
    parameter space toward a single scalar objective.

    Parameters
    ----------
    pipeline : ChainableBlock
        A fully assembled pipeline whose stochastic blocks expose
        ``n_params`` / ``param_bounds``.
    scoring_fn : callable
        ``scoring_fn(molecule) -> float`` — higher is better. Used when no
        batch scorer is supplied and as the fallback for scalar evaluation.
    scoring_batch_fn : callable, optional
        ``scoring_batch_fn(molecules) -> sequence[float]``. Receives an
        ordered batch and must return one score per molecule.
    n_workers : int
        Default number of parallel threads for pipeline execution and scoring.

    Examples
    --------
    >>> opt = SOOptimizer(pipeline, scoring_fn=lambda mol: mol.count_atoms())
    >>> opt.run("pso", steps=50, population=20)
    >>> best = opt.top(5)
    """

    def __init__(
        self,
        pipeline: ChainableBlock,
        scoring_fn=None,
        n_workers: int = 1,
        scoring_batch_fn=None,
    ):
        if scoring_fn is None and scoring_batch_fn is None:
            raise ValueError("provide scoring_fn or scoring_batch_fn")
        self.pipeline = pipeline
        self.scoring_fn = scoring_fn
        self.scoring_batch_fn = scoring_batch_fn
        self.n_workers = n_workers
        self._results: list[tuple[float, object, object]] = []
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
            Show a progress bar.
        n_workers : int, optional
            Number of parallel threads for pipeline execution and scoring.
            Defaults to the value set at construction time.
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
                    dispatch[strategy](
                        steps=steps, population=population, verbose=verbose, **kwargs
                    )
                finally:
                    self._executor = None
        else:
            dispatch[strategy](
                steps=steps, population=population, verbose=verbose, **kwargs
            )
        self._results.sort(key=lambda x: x[0], reverse=True)

    def top(self, n: int = 10) -> list:
        """Return the top *n* unique molecules (best scoring first)."""
        return [mol for _, mol, _ in self._unique_results()[:n]]

    def scores(self, n: int = None) -> list[tuple[float, object]]:
        """Return ``[(score, molecule), ...]`` sorted best-first, deduplicated."""
        results = self._unique_results()
        pairs = [(score, mol) for score, mol, _ in results]
        return pairs[:n] if n is not None else pairs

    def _unique_results(self) -> list:
        seen: set = set()
        unique = []
        for score, mol, params in self._results:
            try:
                smi = mol.to_smiles() if hasattr(mol, "to_smiles") else None
            except Exception:
                smi = None
            if smi in seen:
                continue
            seen.add(smi)
            unique.append((score, mol, params))
        return unique

    def _elite_params(self, n: int, required_len: int) -> list:
        elite = []
        for _, _, p in self._results:
            if p is not None and len(p) == required_len:
                elite.append(p)
                if len(elite) == n:
                    break
        return elite

    # ── core evaluation ───────────────────────────────────────────────────────

    def _evaluate(self, params) -> tuple[float, object, object] | None:
        """Generate and score one parameter vector."""
        generated = self._generate(params)
        if generated is None:
            return None
        mol, stored_params = generated
        try:
            if self.scoring_fn is not None:
                score = self.scoring_fn(mol)
            elif self.scoring_batch_fn is not None:
                score = list(self.scoring_batch_fn([mol]))[0]
            else:
                raise ValueError("no scoring function configured")
            return float(score), mol, stored_params
        except Exception:
            return None

    def _generate(self, params):
        ChainableBlock._inject_params(params)
        try:
            ctx = self.pipeline()
            mol = ctx.molecule
            stored_params = (
                np.array(params) if params is not None and len(params) else None
            )
            return mol, stored_params
        except Exception:
            return None
        finally:
            ChainableBlock._clear_params()

    def _evaluate_batch(self, params_list):
        """Generate an ordered parameter batch and score it in one call."""
        generated = []
        if self._executor is not None:
            generated = list(self._executor.map(self._generate, params_list))
        else:
            generated = [self._generate(params) for params in params_list]

        valid = [
            (index, item) for index, item in enumerate(generated) if item is not None
        ]
        if not valid:
            return [None] * len(params_list)

        molecules = [item[1][0] for item in valid]
        try:
            if self.scoring_batch_fn is not None:
                scores = list(self.scoring_batch_fn(molecules))
            else:
                scores = [self.scoring_fn(molecule) for molecule in molecules]
            if len(scores) != len(molecules):
                raise ValueError("batch scorer must return one score per molecule")
        except Exception:
            return [None] * len(params_list)

        results = [None] * len(params_list)
        for (index, (molecule, stored_params)), score in zip(valid, scores):
            results[index] = (float(score), molecule, stored_params)
        return results

    # ── strategies ────────────────────────────────────────────────────────────

    def _run_random(self, steps: int, population: int = 50, verbose: bool = True, **_):
        bounds = self.pipeline.param_bounds
        if not bounds:
            param_list = [[] for _ in range(steps)]
        else:
            lo = np.array([b[0] for b in bounds], dtype=float)
            hi = np.array([b[1] for b in bounds], dtype=float)
            param_list = [
                lo + np.random.random(len(bounds)) * (hi - lo) for _ in range(steps)
            ]

        bar = progress_bar(total=steps, desc="random", disable=not verbose)
        best = float("-inf")
        try:
            if self._executor is not None:
                if self.scoring_batch_fn is not None:
                    batch_results = self._evaluate_batch(param_list)
                    result_iter = batch_results
                else:
                    from concurrent.futures import as_completed

                    futures = {
                        self._executor.submit(self._evaluate, p) for p in param_list
                    }
                    result_iter = (future.result() for future in as_completed(futures))
                for result in result_iter:
                    bar.update(1)
                    if result is not None:
                        self._results.append(result)
                        if result[0] > best:
                            best = result[0]
                            bar.set_postfix(best=f"{best:.4f}")
            else:
                result_iter = (
                    self._evaluate_batch(param_list)
                    if self.scoring_batch_fn is not None
                    else (self._evaluate(p) for p in param_list)
                )
                for result in result_iter:
                    bar.update(1)
                    if result is not None:
                        self._results.append(result)
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

        elite = self._elite_params(population // 2, n)
        for i, ep in enumerate(elite):
            pos[i] = np.clip(ep, lo, hi)
            pbest_pos[i] = pos[i].copy()
            pbest_score[i] = self._results[i][0]
        if elite:
            gbest_pos = pos[0].copy()
            gbest_score = self._results[0][0]

        bar = progress_bar(range(steps), desc="pso", disable=not verbose)
        for _ in bar:
            params = [pos[i] for i in range(population)]
            raw = self._evaluate_batch(params)

            for i, result in enumerate(raw):
                if result is None:
                    continue
                score, mol, stored_params = result
                self._results.append((score, mol, stored_params))
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

        n = len(bounds)
        lo = np.array([b[0] for b in bounds], dtype=float)
        hi = np.array([b[1] for b in bounds], dtype=float)

        bar = progress_bar(total=steps, desc="scipy", disable=not verbose)
        best = float("-inf")

        def objective(params):
            nonlocal best
            result = self._evaluate(np.clip(params, lo, hi))
            if result is None:
                return 0.0
            score, mol, stored_params = result
            self._results.append((score, mol, stored_params))
            bar.update(1)
            if score > best:
                best = score
                bar.set_postfix(best=f"{best:.4f}")
            return -score

        elite = self._elite_params(1, n)
        x0 = (
            np.clip(elite[0], lo, hi) if elite else lo + np.random.random(n) * (hi - lo)
        )
        try:
            minimize(objective, x0, method=scipy_method, options={"maxiter": steps})
        finally:
            bar.close()
