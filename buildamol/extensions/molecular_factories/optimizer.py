"""
Optimizer — unified entry point for single- and multi-objective optimisation.

Delegates to :class:`SOOptimizer` (``mode="single"``) or
:class:`MOOptimizer` (``mode="multi"``) depending on the *mode* argument.

Single-objective
----------------
``scoring_fn(mol) -> float``; higher is better.

    >>> opt = Optimizer(pipeline, admet_score)
    >>> opt.run("pso", steps=100, population=50)
    >>> best = opt.top(10)

Multi-objective
---------------
``objectives_fn(mol) -> np.ndarray``; every element is **minimised**
(negate to maximise).

    >>> def objectives(mol):
    ...     pred = model.predict(smiles=mol.to_smiles())
    ...     return np.array([-pred["HIA_Hou"], pred["hERG"], sa_score(mol)])
    >>> opt = Optimizer(pipeline, objectives, mode="multi")
    >>> opt.run("nsga2", steps=50, population=100)
    >>> front = opt.pareto_front()
"""

from .so_optimizer import SOOptimizer
from .mo_optimizer import MOOptimizer


class Optimizer:
    """
    Unified optimiser wrapper.

    Parameters
    ----------
    pipeline : ChainableBlock
    scoring_fn : callable
        ``(mol) -> float`` for single-objective, or
        ``(mol) -> np.ndarray`` (minimise each element) for multi-objective.
    mode : {"single", "multi"}
        Select the underlying optimiser.  Defaults to ``"single"``.
    n_workers : int
        Default number of parallel evaluation threads.

    All methods and attributes are forwarded to the underlying
    :class:`SOOptimizer` or :class:`MOOptimizer` instance, accessible via
    ``opt.optimizer``.
    """

    def __init__(
        self,
        pipeline,
        scoring_fn,
        mode: str = "single",
        n_workers: int = 1,
    ):
        if mode == "single":
            self.optimizer = SOOptimizer(pipeline, scoring_fn, n_workers=n_workers)
        elif mode == "multi":
            self.optimizer = MOOptimizer(pipeline, scoring_fn, n_workers=n_workers)
        else:
            raise ValueError(f"mode must be 'single' or 'multi', got {mode!r}")

    def __getattr__(self, name: str):
        # forward everything not found on Optimizer itself to the inner optimizer
        return getattr(self.optimizer, name)

    def __repr__(self) -> str:
        return f"Optimizer(mode={'single' if isinstance(self.optimizer, SOOptimizer) else 'multi'}, optimizer={self.optimizer!r})"
