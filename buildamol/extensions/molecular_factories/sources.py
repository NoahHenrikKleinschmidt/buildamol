import numpy as np
from .base import ChainableBlock, Context


class Choice(ChainableBlock):
    """
    Returns a molecule from a list of candidates.

    In stochastic mode a molecule is picked at random (respecting ``p``).
    When the ``Optimizer`` injects a parameter vector, the single parameter
    is interpreted as an integer index ``[0, len(molecules)-1]``.
    """

    n_params = 1

    def __init__(self, molecules, p=None, seed=None):
        self.molecules = molecules
        self.p = p
        self.seed = seed

    @property
    def param_bounds(self) -> list:
        return [(0, len(self.molecules) - 1)]

    def __call__(self, *args, **kwargs) -> Context:
        param = self._next_param()
        if param is not None:
            idx = int(param) % len(self.molecules)
        else:
            if self.seed is not None:
                np.random.seed(self.seed)
            idx = int(np.random.choice(len(self.molecules), p=self.p))
        context = Context()
        context.molecule = self.molecules[idx].copy()
        return context


class Compound(ChainableBlock):
    """
    Always returns the same molecule (optionally copying it).
    ``n_params = 0`` — no stochastic element, nothing to optimise.
    """

    n_params = 0

    def __init__(self, molecule, copy=True):
        self.molecule = molecule
        self.copy = copy

    def __call__(self, *args, **kwargs) -> Context:
        context = Context()
        context.molecule = self.molecule.copy() if self.copy else self.molecule
        return context
