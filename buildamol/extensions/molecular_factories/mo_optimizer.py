"""
MOOptimizer — multi-objective variant of the BuildAMol Optimizer.

Objectives are defined as a vector-valued function
``objectives_fn(mol) -> np.ndarray`` where **lower is better for every
element** (minimisation convention).  To maximise a property negate it in the
objectives function:

    return np.array([-pred["HIA_Hou"], pred["hERG"], ...])

Strategies
----------
random   Uniform random search; maintains a Pareto archive.
nsga2    NSGA-II with uniform crossover and random-reset mutation.
"""

import numpy as np

from .base import ChainableBlock
from buildamol.utils.auxiliary import progress_bar

# ── Pareto utilities ──────────────────────────────────────────────────────────


def _dominates(a: np.ndarray, b: np.ndarray) -> bool:
    """True if *a* dominates *b* (a ≤ b on all objectives, a < b on ≥ 1)."""
    return bool(np.all(a <= b) and np.any(a < b))


def _fast_non_dominated_sort(objectives: list) -> list[list[int]]:
    """NSGA-II fast non-dominated sort.

    Returns a list of fronts (each front is a list of indices into *objectives*).
    Front 0 is the Pareto-optimal set.
    """
    n = len(objectives)
    dominated_by = [[] for _ in range(n)]  # dominated_by[i] = indices i dominates
    n_dom = [0] * n  # how many individuals dominate i
    fronts: list[list[int]] = [[]]

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if _dominates(objectives[i], objectives[j]):
                dominated_by[i].append(j)
            elif _dominates(objectives[j], objectives[i]):
                n_dom[i] += 1
        if n_dom[i] == 0:
            fronts[0].append(i)

    k = 0
    while fronts[k]:
        next_front: list[int] = []
        for i in fronts[k]:
            for j in dominated_by[i]:
                n_dom[j] -= 1
                if n_dom[j] == 0:
                    next_front.append(j)
        fronts.append(next_front)
        k += 1

    return [f for f in fronts if f]


def _crowding_distance(objectives: list) -> np.ndarray:
    """Crowding distance for the members of a single front."""
    n = len(objectives)
    if n == 0:
        return np.array([])
    distances = np.zeros(n)
    obj_array = np.array(objectives)  # (n, m)
    m = obj_array.shape[1]
    for k in range(m):
        col = obj_array[:, k]
        order = np.argsort(col)
        distances[order[0]] = distances[order[-1]] = np.inf
        obj_range = col[order[-1]] - col[order[0]]
        if obj_range == 0:
            continue
        for i in range(1, n - 1):
            distances[order[i]] += (col[order[i + 1]] - col[order[i - 1]]) / obj_range
    return distances


def _tournament_select(
    objectives: list, ranks: list[int], crowd: list[float], n: int
) -> list[int]:
    """Binary tournament selection based on Pareto rank and crowding distance."""
    pool = list(range(len(objectives)))
    selected = []
    for _ in range(n):
        a, b = np.random.choice(pool, 2, replace=False)
        if ranks[a] < ranks[b]:
            selected.append(a)
        elif ranks[a] > ranks[b]:
            selected.append(b)
        else:
            selected.append(a if crowd[a] >= crowd[b] else b)
    return selected


def _rank_and_crowd(objectives: list) -> tuple[list[int], list[float]]:
    """Compute per-individual Pareto rank and crowding distance."""
    fronts = _fast_non_dominated_sort(objectives)
    rank = [0] * len(objectives)
    crowd = [0.0] * len(objectives)
    for r, front in enumerate(fronts):
        for i in front:
            rank[i] = r
        dist = _crowding_distance([objectives[i] for i in front])
        for pos, i in enumerate(front):
            crowd[i] = float(dist[pos])
    return rank, crowd


def _survivor_selection(objectives: list, n: int) -> list[int]:
    """Select *n* survivors by Pareto rank then crowding distance."""
    fronts = _fast_non_dominated_sort(objectives)
    selected: list[int] = []
    for front in fronts:
        if len(selected) + len(front) <= n:
            selected.extend(front)
        else:
            remaining = n - len(selected)
            dist = _crowding_distance([objectives[i] for i in front])
            order = np.argsort(-dist)[:remaining]
            selected.extend(front[int(i)] for i in order)
            break
    return selected


# ── Main class ────────────────────────────────────────────────────────────────


class MOOptimizer:
    """
    Multi-objective optimisation of a stochastic BuildAMol pipeline.

    Parameters
    ----------
    pipeline : ChainableBlock
        A fully assembled pipeline whose stochastic blocks expose
        ``n_params`` / ``param_bounds``.
    objectives_fn : callable
        ``objectives_fn(molecule) -> np.ndarray`` — a 1-D array where
        **each element is minimised**.  Negate to maximise a property.
    n_workers : int
        Default number of parallel threads.

    Examples
    --------
    >>> def objectives(mol):
    ...     pred = admet_model.predict(smiles=mol.to_smiles())
    ...     return np.array([
    ...         -pred["HIA_Hou"],          # maximise → negate
    ...         -pred["Bioavailability_Ma"],
    ...         pred["hERG"],              # minimise (bad = high)
    ...         pred["AMES"],
    ...         pred["ClinTox"],
    ...         sa_score(mol),
    ...     ])
    >>> opt = MOOptimizer(pipeline, objectives_fn=objectives)
    >>> opt.run("nsga2", steps=50, population=100)
    >>> front = opt.pareto_front()
    """

    def __init__(self, pipeline: ChainableBlock, objectives_fn, n_workers: int = 1):
        self.pipeline = pipeline
        self.objectives_fn = objectives_fn
        self.n_workers = n_workers
        self._results: list[tuple[np.ndarray, object, object]] = []
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
        Run the multi-objective optimisation.

        Parameters
        ----------
        strategy : {"random", "nsga2"}
        steps : int
            Evaluations (random) or generations (nsga2).
        population : int
            Population size for nsga2.
        verbose : bool
            Show a progress bar.
        n_workers : int, optional
            Parallel threads for pipeline execution and scoring.
        **kwargs
            Forwarded to the strategy (e.g. ``crossover_prob``, ``mutation_prob``
            for nsga2).
        """
        dispatch = {
            "random": self._run_random,
            "nsga2": self._run_nsga2,
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

    def pareto_front(self) -> list:
        """Return all molecules on the current Pareto front (non-dominated set)."""
        unique = self._unique_results()
        if not unique:
            return []
        objectives = [obj for obj, _, _ in unique]
        fronts = _fast_non_dominated_sort(objectives)
        if not fronts:
            return []
        return [unique[i][1] for i in fronts[0]]

    def top(self, n: int = None) -> list:
        """Return up to *n* Pareto-front molecules sorted by crowding distance.

        Molecules with high crowding distance are better spread across the front
        and represent more diverse trade-off points.  If *n* is ``None``, the
        entire Pareto front is returned.
        """
        unique = self._unique_results()
        if not unique:
            return []
        objectives = [obj for obj, _, _ in unique]
        fronts = _fast_non_dominated_sort(objectives)
        if not fronts:
            return []
        front_idx = fronts[0]
        dist = _crowding_distance([objectives[i] for i in front_idx])
        order = np.argsort(-dist)
        ranked = [unique[front_idx[int(i)]][1] for i in order]
        return ranked[:n] if n is not None else ranked

    def scores(self, n: int = None) -> list[tuple[np.ndarray, object]]:
        """Return ``[(objectives_array, molecule), ...]`` for the Pareto front.

        Sorted by crowding distance (most diverse first).
        """
        unique = self._unique_results()
        if not unique:
            return []
        objectives = [obj for obj, _, _ in unique]
        fronts = _fast_non_dominated_sort(objectives)
        if not fronts:
            return []
        front_idx = fronts[0]
        dist = _crowding_distance([objectives[i] for i in front_idx])
        order = np.argsort(-dist)
        pairs = [
            (objectives[front_idx[int(i)]], unique[front_idx[int(i)]][1]) for i in order
        ]
        return pairs[:n] if n is not None else pairs

    def to_dataframe(
        self,
        n: int = None,
        pareto_only: bool = False,
        molecules: bool = False,
        objective_names: list = None,
    ):
        """
        Return the results as a ``pandas.DataFrame`` with one column per
        objective plus ``'smiles'`` and ``'is_pareto'``.

        Parameters
        ----------
        n : int, optional
            Only include the first *n* rows.
        pareto_only : bool
            Only include molecules on the Pareto front (sorted by crowding
            distance). Otherwise all unique results are included.
        molecules : bool
            Also include a ``'molecule'`` column with the Molecule objects.
        objective_names : list, optional
            Names to use for the objective columns. Defaults to
            ``objective_0``, ``objective_1``, ...
        """
        import pandas as pd

        unique = self._unique_results()
        if not unique:
            return pd.DataFrame()

        front = set()
        objectives = [obj for obj, _, _ in unique]
        fronts = _fast_non_dominated_sort(objectives)
        if fronts:
            front = set(fronts[0])

        if pareto_only:
            front_idx = list(fronts[0]) if fronts else []
            dist = _crowding_distance([objectives[i] for i in front_idx])
            indices = [front_idx[int(i)] for i in np.argsort(-dist)]
        else:
            indices = range(len(unique))

        rows = []
        for i in indices:
            obj, mol, _ = unique[i]
            try:
                smiles = mol.to_smiles() if hasattr(mol, "to_smiles") else None
            except Exception:
                smiles = None
            row = {}
            for j, value in enumerate(np.atleast_1d(obj)):
                name = (
                    objective_names[j]
                    if objective_names is not None and j < len(objective_names)
                    else f"objective_{j}"
                )
                row[name] = float(value)
            row["smiles"] = smiles
            row["is_pareto"] = i in front
            if molecules:
                row["molecule"] = mol
            rows.append(row)
        return pd.DataFrame(rows[:n] if n is not None else rows)

    def _unique_results(self) -> list:
        seen: set = set()
        unique = []
        for obj, mol, params in self._results:
            try:
                smi = mol.to_smiles() if hasattr(mol, "to_smiles") else None
            except Exception:
                smi = None
            if smi in seen:
                continue
            seen.add(smi)
            unique.append((obj, mol, params))
        return unique

    # ── core evaluation ───────────────────────────────────────────────────────

    def _evaluate(self, params) -> tuple[np.ndarray, object, object] | None:
        ChainableBlock._inject_params(params)
        try:
            ctx = self.pipeline()
            mol = ctx.molecule
            obj = np.asarray(self.objectives_fn(mol), dtype=float)
            stored_params = (
                np.array(params) if params is not None and len(params) else None
            )
            return obj, mol, stored_params
        except Exception:
            return None
        finally:
            ChainableBlock._clear_params()

    def _evaluate_many(self, param_list: list) -> list:
        if self._executor is not None:
            raw = list(self._executor.map(self._evaluate, param_list))
        else:
            raw = [self._evaluate(p) for p in param_list]
        return [r for r in raw if r is not None]

    # ── strategies ────────────────────────────────────────────────────────────

    def _run_random(self, steps: int, population: int = 50, verbose: bool = True, **_):
        """Uniform random search; all results are kept in a Pareto archive."""
        bounds = self.pipeline.param_bounds
        if not bounds:
            param_list = [[] for _ in range(steps)]
        else:
            lo = np.array([b[0] for b in bounds], dtype=float)
            hi = np.array([b[1] for b in bounds], dtype=float)
            param_list = [
                lo + np.random.random(len(bounds)) * (hi - lo) for _ in range(steps)
            ]

        bar = progress_bar(total=steps, desc="random (MO)", disable=not verbose)
        try:
            if self._executor is not None:
                from concurrent.futures import as_completed

                futures = {self._executor.submit(self._evaluate, p) for p in param_list}
                for fut in as_completed(futures):
                    result = fut.result()
                    bar.update(1)
                    if result is not None:
                        self._results.append(result)
            else:
                for p in param_list:
                    result = self._evaluate(p)
                    bar.update(1)
                    if result is not None:
                        self._results.append(result)
        finally:
            bar.close()
        bar.set_postfix(front=len(self.pareto_front()))

    def _run_nsga2(
        self,
        steps: int,
        population: int = 50,
        crossover_prob: float = 0.9,
        mutation_prob: float = 0.1,
        verbose: bool = True,
        **_,
    ):
        """NSGA-II: non-dominated sorting + crowding-distance selection.

        Uses uniform crossover and random-reset mutation on the pipeline's
        discrete parameter vectors.

        Parameters
        ----------
        steps : int
            Number of generations.
        population : int
            Population size (must be even).
        crossover_prob : float
            Probability of applying uniform crossover to a parent pair.
        mutation_prob : float
            Per-gene probability of random reset within bounds.
        """
        bounds = self.pipeline.param_bounds
        if not bounds:
            self._run_random(steps=steps * population, verbose=verbose)
            return

        if population % 2 != 0:
            population += 1

        n_params = len(bounds)
        lo = np.array([b[0] for b in bounds], dtype=float)
        hi = np.array([b[1] for b in bounds], dtype=float)

        # ── initialise ────────────────────────────────────────────────────────
        pop_params = [
            lo + np.random.random(n_params) * (hi - lo) for _ in range(population)
        ]
        init_results = self._evaluate_many(pop_params)
        self._results.extend(init_results)

        # align pop_params with evaluated results (some may have failed)
        pop_params = [r[2] for r in init_results if r[2] is not None]
        pop_objs = [r[0] for r in init_results]
        pop_mols = [r[1] for r in init_results]

        bar = progress_bar(range(steps), desc="nsga2", disable=not verbose)
        for _ in bar:
            if len(pop_objs) < 2:
                break

            rank, crowd = _rank_and_crowd(pop_objs)
            parent_idx = _tournament_select(pop_objs, rank, crowd, population)

            # crossover + mutation → offspring
            offspring_params = []
            for i in range(0, len(parent_idx) - 1, 2):
                p1 = pop_params[parent_idx[i] % len(pop_params)]
                p2 = pop_params[parent_idx[i + 1] % len(pop_params)]
                if np.random.random() < crossover_prob:
                    mask = np.random.randint(0, 2, n_params).astype(bool)
                    c1 = np.where(mask, p1, p2).astype(float)
                    c2 = np.where(mask, p2, p1).astype(float)
                else:
                    c1, c2 = p1.copy(), p2.copy()
                for c in (c1, c2):
                    mut = np.random.random(n_params) < mutation_prob
                    c[mut] = lo[mut] + np.random.random(int(mut.sum())) * (
                        hi[mut] - lo[mut]
                    )
                offspring_params.extend([c1, c2])

            off_results = self._evaluate_many(offspring_params)
            self._results.extend(off_results)

            # combine parent + offspring, select next generation
            comb_objs = pop_objs + [r[0] for r in off_results]
            comb_mols = pop_mols + [r[1] for r in off_results]
            comb_params = pop_params + [r[2] for r in off_results if r[2] is not None]

            survivors = _survivor_selection(comb_objs, population)
            pop_objs = [comb_objs[i] for i in survivors]
            pop_mols = [comb_mols[i] for i in survivors]
            pop_params = [comb_params[min(i, len(comb_params) - 1)] for i in survivors]

            front_size = sum(
                1 for i in survivors if _rank_and_crowd(comb_objs)[0][i] == 0
            )
            bar.set_postfix(front=front_size, pop=len(pop_objs))
