"""
The `Generator` class combines a fragment `Assembler` with a user-supplied scoring
function and a pluggable search strategy to find high-scoring molecules automatically.

The scoring function must accept a single `Molecule` argument and return a float.
Any conversion to SMILES, fingerprints, coordinates, etc. is the responsibility of
the scoring function.

Note
----
The Generator is a general purpose lightweight factory using the ``Assembler`` class.
Therefore the generation process cannot be customized in detail. If you need a more controlled
setup, check out the stand-alone building blocks of the newer factories API.

Usage
-----
1. Prepare a list of fragment molecules (same as for ``Assembler``).
2. Define a ``scoring_fn(mol) -> float``.
3. Create a ``Generator`` and call ``.run()``.

Example
-------

.. code-block:: python

    import buildamol as bam
    from buildamol.extensions.molecular_factories import Generator
    from rdkit.Chem import QED

    bam.load_small_molecules()
    fragments = [bam.molecule(name) for name in ("benzene", "isopropanol", "dimethylamine")]
    fragments = [i[0] if isinstance(i, list) else i for i in fragments]

    def qed_score(mol):
        return QED.qed(mol.to_rdkit())

    gen = Generator(fragments, qed_score, n_fragments=3, maximize=True)
    gen.run(n_steps=50, method="random")
    gen.run(n_steps=10,  method="genetic", population_size=20)
    gen.run(n_steps=10,  method="swarm",   n_particles=15)

    print(gen.to_dataframe().head())


========== ===========================
    score   smiles
========== ===========================
 0.828336   CC(O)(Cc1ccccc1)c1ccccc1
 0.827453   CC(O)Cc1ccc(-c2ccccc2)cc1
 0.827453   CC(O)Cc1cccc(-c2ccccc2)c1
 0.827453   CC(O)Cc1ccccc1-c1ccccc1
 0.825895   OC(Cc1ccccc1)Cc1ccccc1
========== ===========================


.. image:: examples/files/generator_output_example1.png


Search methods
--------------
``"random"``
    Pure random sampling — fastest baseline, good for large fragment libraries.

``"genetic"``
    Genetic algorithm over instruction matrices.
    Key kwargs: ``population_size`` (default 20), ``mutation_rate`` (default 0.15).

``"swarm"``
    Particle swarm optimisation over the continuous relaxation of the integer
    instruction matrix; matrices are rounded at each evaluation.
    Key kwargs: ``n_particles`` (default 20), ``inertia`` (default 0.7),
    ``cognitive`` (default 1.5), ``social`` (default 1.5).

``"scipy"``
    Delegates to ``scipy.optimize.minimize`` with a rounded-matrix wrapper.
    Key kwargs: ``scipy_method`` (default ``"Nelder-Mead"``).
"""

import numpy as np
from buildamol.extensions.molecular_factories.assembler import Assembler
from buildamol.utils.auxiliary import progress_bar
import contextlib
import io


class Generator:
    """
    Automated molecule generator using fragment assembly and a scoring function.

    Parameters
    ----------
    fragments : list
        Fragment molecules to assemble from (forwarded to ``Assembler``).
    scoring_fn : callable
        ``scoring_fn(mol) -> float``.  Receives a ``Molecule``; must return a
        scalar score.  Higher is better when ``maximize=True``.
    n_fragments : int
        Number of fragments per assembled molecule.
    maximize : bool
        ``True`` → search for the highest score.
        ``False`` → search for the lowest score.
    preprocess : callable, optional
        ``preprocess(ndarray) -> ndarray``.  If provided, called on each instruction matrix
        before assembly.  Can be used to enforce constraints on the search space. Requires knowledge of the assembler's fragment library and attachment points.
    postprocess : callable, optional
        ``postprocess(mol) -> mol``.  If provided, called on each assembled molecule
        before scoring.  Can be used to perform sanitisation, cyclisation, or other
        transformations that the ``Assembler`` does not handle automatically.
    n_workers : int
        Number of parallel workers for batch scoring.  Defaults to ``1``
        (serial).  When ``> 1``, ``ThreadPoolExecutor`` is used so the
        scoring function must be thread-safe.  Most ML-based scoring
        functions (admet-ai, RDKit QED, etc.) release the GIL and
        benefit from threading without any extra setup.
    """

    def __init__(
        self,
        fragments,
        scoring_fn,
        n_fragments=3,
        maximize=True,
        postprocess=None,
        preprocess=None,
        n_workers=1,
    ):
        if isinstance(fragments, Assembler):
            self.assembler = fragments
        else:
            self.assembler = Assembler(list(fragments))
        self.scoring_fn = scoring_fn
        self.n_fragments = n_fragments
        self.maximize = maximize
        self.postprocess = postprocess
        self.preprocess = preprocess
        self.n_workers = n_workers
        self._results = []
        self._executor = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, n_steps=100, method="random", verbose=False, **kwargs) -> "Generator":
        """
        Run the molecule generation / optimisation loop.

        Parameters
        ----------
        n_steps : int
            Number of steps (random samples, generations, PSO iterations, or
            scipy ``maxiter``, depending on ``method``).
        method : str
            ``"random"``, ``"genetic"``, ``"swarm"``, or ``"scipy"``.
        verbose : bool
            Show a tqdm progress bar (requires ``tqdm`` to be installed).
        **kwargs
            Method-specific keyword arguments (see module docstring).

        Returns
        -------
        Generator
            ``self``, so ``run()`` calls can be chained.
        """
        _dispatch = {
            "random": self._run_random,
            "genetic": self._run_genetic,
            "swarm": self._run_swarm,
            "scipy": self._run_scipy,
        }
        if method not in _dispatch:
            raise ValueError(f"Unknown method '{method}'. Available: {list(_dispatch)}")
        n_workers = kwargs.pop("n_workers", self.n_workers)
        if n_workers > 1:
            from concurrent.futures import ThreadPoolExecutor

            with ThreadPoolExecutor(max_workers=n_workers) as exe:
                self._executor = exe
                try:
                    _dispatch[method](n_steps=n_steps, verbose=verbose, **kwargs)
                finally:
                    self._executor = None
        else:
            _dispatch[method](n_steps=n_steps, verbose=verbose, **kwargs)
        return self

    @property
    def best(self):
        """The single best molecule found so far, or ``None``."""
        ranked = self._ranked()
        return ranked[0]["molecule"] if ranked else None

    def top(self, k=10):
        """
        Return the top-*k* results as a list of dicts with keys
        ``'score'``, ``'molecule'``, and ``'matrix'``.
        """
        return self._ranked()[:k]

    def to_dataframe(self):
        """
        Return all valid results as a ``pandas.DataFrame`` with columns
        ``'score'`` and ``'smiles'``, sorted by score.
        """
        import pandas as pd

        rows = []
        for r in self._ranked():
            try:
                smiles = r["molecule"].to_smiles()
            except Exception:
                smiles = None
            rows.append({"score": r["score"], "smiles": smiles})
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _worst(self):
        return float("-inf") if self.maximize else float("inf")

    def _is_better(self, a, b):
        return a > b if self.maximize else a < b

    def _evaluate_batch(self, matrices):
        """Evaluate a list of matrices, returning ``[(score, mol), ...]``.

        Assembly is always serial (pure-Python, GIL-bound).  Only the
        user-supplied scoring function is parallelised — it typically releases
        the GIL (e.g. PyTorch / RDKit C extensions) and is where the real
        compute lives.  The executor is created once per ``run()`` call and
        reused across all batches, so genetic/swarm iterations don't pay
        pool-startup cost on every generation.
        """
        # --- serial assembly ---
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
            io.StringIO()
        ):
            assembled = []
            for m in matrices:
                try:
                    if self.preprocess is not None:
                        m = self.preprocess(m)
                    mol = self.assembler.make(m)
                    if self.postprocess is not None:
                        mol = self.postprocess(mol)
                    assembled.append(mol)
                except Exception:
                    assembled.append(None)

        # --- parallel scoring (reuse the run()-level executor if available) ---
        worst = self._worst()

        def _score(mol):
            if mol is None:
                return worst, None
            try:
                return float(self.scoring_fn(mol)), mol
            except Exception:
                return worst, None

        if self._executor is not None:
            return list(self._executor.map(_score, assembled))
        return [_score(mol) for mol in assembled]

    def _evaluate(self, matrix):
        """Return ``(score, mol)`` or ``(worst, None)`` on any error."""
        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(
                io.StringIO()
            ):
                if self.preprocess is not None:
                    matrix = self.preprocess(matrix)
                mol = self.assembler.make(matrix)
                if self.postprocess is not None:
                    mol = self.postprocess(mol)
            return float(self.scoring_fn(mol)), mol
        except Exception:
            return self._worst(), None

    def _record(self, score, matrix, mol):
        if mol is not None and score != self._worst():
            self._results.append(
                {"score": score, "molecule": mol, "matrix": matrix.copy()}
            )

    def _ranked(self):
        seen = set()
        unique = []
        for r in sorted(self._results, key=lambda r: r["score"], reverse=self.maximize):
            try:
                smi = r["molecule"].to_smiles()
            except Exception:
                smi = None
            if smi in seen:
                continue
            seen.add(smi)
            unique.append(r)
        return unique

    def _pos_to_matrix(self, pos, n_total_frags, max_ap):
        """Round a float PSO/scipy position array to a valid integer matrix."""
        matrix = np.round(pos).astype(int)
        matrix[:, 0] = np.clip(matrix[:, 0], 0, n_total_frags - 1)
        matrix[:, 1:] = np.clip(matrix[:, 1:], 0, max_ap - 1)
        matrix[0, 1:] = -1  # first row's attachment columns are unused
        return matrix

    # ------------------------------------------------------------------
    # Search strategies
    # ------------------------------------------------------------------

    def _run_random(self, n_steps, verbose=False, **_):
        """Sample ``n_steps`` random molecules and score each one."""
        matrices = [self.assembler.random(self.n_fragments) for _ in range(n_steps)]
        bar = progress_bar(
            zip(matrices, self._evaluate_batch(matrices)),
            total=n_steps,
            desc="random",
            disable=not verbose,
        )
        for matrix, (score, mol) in bar:
            self._record(score, matrix, mol)

    def _run_genetic(
        self, n_steps, population_size=20, mutation_rate=0.15, verbose=False, **_
    ):
        """
        Genetic algorithm over instruction matrices.

        *Selection*: binary tournament.
        *Crossover*: single-point row swap between two parents.
        *Mutation*: replace one row with a row sampled from a fresh random matrix.

        ``n_steps`` = number of generations.
        """
        # Initialise population
        pop = [self.assembler.random(self.n_fragments) for _ in range(population_size)]
        eval_cache = self._evaluate_batch(pop)
        for (score, mol), matrix in zip(eval_cache, pop):
            self._record(score, matrix, mol)

        best_score = self._worst()
        bar = progress_bar(range(n_steps), desc="genetic", disable=not verbose)
        for _ in bar:
            scores = [s for s, _ in eval_cache]
            gen_best = max(scores) if self.maximize else min(scores)
            if self._is_better(gen_best, best_score):
                best_score = gen_best
            bar.set_postfix(best=f"{best_score:.4f}")

            # Tournament selection
            parents = []
            for _ in range(population_size):
                a, b = np.random.choice(len(pop), 2, replace=False)
                parents.append(
                    pop[a] if self._is_better(scores[a], scores[b]) else pop[b]
                )

            # Crossover — single-point row swap
            new_pop = []
            for i in range(0, population_size - 1, 2):
                p1, p2 = parents[i], parents[i + 1]
                pt = np.random.randint(1, self.n_fragments)
                new_pop.append(np.vstack([p1[:pt], p2[pt:]]))
                new_pop.append(np.vstack([p2[:pt], p1[pt:]]))
            if len(new_pop) < population_size:
                new_pop.append(parents[-1].copy())

            # Mutation — replace one row with a row from a fresh random matrix
            for matrix in new_pop:
                if np.random.random() < mutation_rate:
                    row = np.random.randint(1, self.n_fragments)
                    matrix[row] = self.assembler.random(self.n_fragments)[row]

            pop = new_pop
            eval_cache = self._evaluate_batch(pop)
            for (score, mol), matrix in zip(eval_cache, pop):
                self._record(score, matrix, mol)

    def _run_swarm(
        self,
        n_steps,
        n_particles=20,
        inertia=0.7,
        cognitive=1.5,
        social=1.5,
        verbose=False,
        **_,
    ):
        """
        Particle swarm optimisation.

        Particles move through a continuous relaxation of the integer instruction
        matrix space; positions are rounded to integers before each evaluation.

        ``n_steps`` = number of PSO iterations.
        """
        n_total_frags = len(self.assembler.fragments)
        max_ap = max(len(ap) for ap in self.assembler.attachment_points)

        # Initialise positions uniformly in valid integer ranges
        pos = np.random.uniform(0, 1, (n_particles, self.n_fragments, 3))
        pos[:, :, 0] *= n_total_frags
        pos[:, :, 1:] *= max_ap

        vel = np.zeros_like(pos)
        pbest_pos = pos.copy()
        pbest_score = np.full(n_particles, self._worst())
        gbest_pos = None
        gbest_score = self._worst()

        bar = progress_bar(range(n_steps), desc="swarm", disable=not verbose)
        for _ in bar:
            matrices = [
                self._pos_to_matrix(pos[i], n_total_frags, max_ap)
                for i in range(n_particles)
            ]
            results = self._evaluate_batch(matrices)
            for i, (matrix, (score, mol)) in enumerate(zip(matrices, results)):
                self._record(score, matrix, mol)
                if self._is_better(score, pbest_score[i]):
                    pbest_score[i] = score
                    pbest_pos[i] = pos[i].copy()
                    if self._is_better(score, gbest_score):
                        gbest_score = score
                        gbest_pos = pos[i].copy()
            bar.set_postfix(best=f"{gbest_score:.4f}")

            if gbest_pos is None:
                continue

            r1 = np.random.random(pos.shape)
            r2 = np.random.random(pos.shape)
            vel = (
                inertia * vel
                + cognitive * r1 * (pbest_pos - pos)
                + social * r2 * (gbest_pos - pos)
            )
            pos = pos + vel
            pos[:, :, 0] = np.clip(pos[:, :, 0], 0, n_total_frags - 1)
            pos[:, :, 1:] = np.clip(pos[:, :, 1:], 0, max_ap - 1)

    def _run_scipy(self, n_steps, scipy_method="Nelder-Mead", **_):
        """
        Delegate to ``scipy.optimize.minimize``.

        The instruction matrix is flattened to a 1-D float vector for scipy and
        rounded back to integers at each function evaluation.
        ``n_steps`` is passed as ``maxiter`` in scipy's options dict.
        """
        from scipy import optimize

        n_total_frags = len(self.assembler.fragments)
        max_ap = max(len(ap) for ap in self.assembler.attachment_points)
        starting = self.assembler.random(self.n_fragments).astype(float)

        def objective(flat):
            matrix = self._pos_to_matrix(
                flat.reshape(self.n_fragments, 3), n_total_frags, max_ap
            )
            score, mol = self._evaluate(matrix)
            self._record(score, matrix, mol)
            return -score if self.maximize else score

        optimize.minimize(
            objective,
            starting.flatten(),
            method=scipy_method,
            options={"maxiter": n_steps},
        )


# ------------------------------------------------------------------

if __name__ == "__main__":
    import buildamol as bam
    from rdkit.Chem import QED

    bam.load_small_molecules()
    fragments = [
        bam.Molecule.from_smiles("c1ccccc1", id="BNZ").autolabel(),
        bam.Molecule.from_smiles("CCO", id="EOH").autolabel(),
        bam.Molecule.from_smiles("CNC", id="DMA").autolabel(),
        bam.Molecule.from_smiles("CC(=O)O", id="AcA").autolabel(),
        bam.Molecule.from_smiles("c1ccncc1", id="Pyr").autolabel(),
    ]

    def qed_score(mol):
        return QED.qed(mol.to_rdkit())

    gen = Generator(fragments, qed_score, n_fragments=3, maximize=True)

    print("Random search...")
    gen.run(50, method="random")
    print(f"  best so far: {gen.top(1)[0]['score']:.3f}")

    print("Genetic search...")
    gen.run(10, method="genetic", population_size=10)
    print(f"  best so far: {gen.top(1)[0]['score']:.3f}")

    print("Swarm search...")
    gen.run(10, method="swarm", n_particles=10)
    print(f"  best so far: {gen.top(1)[0]['score']:.3f}")

    print("Scipy search...")
    gen.run(20, method="scipy", scipy_method="Nelder-Mead")
    print(f"  best so far: {gen.top(1)[0]['score']:.3f}")

    print("\nTop 3 results:")
    for r in gen.top(3):
        print(f"  score={r['score']:.3f}  smiles={r['molecule'].to_smiles()}")

    df = gen.to_dataframe()
    print(f"\nTotal valid results: {len(df)}")
    print(df.head())
