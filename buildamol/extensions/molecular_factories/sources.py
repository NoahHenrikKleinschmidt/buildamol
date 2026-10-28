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

    def _call_buildamol_native(self, *args, **kwargs) -> "Context":
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

    def _call_rdkit_accelerated(self, *args, **kwargs) -> "Context":
        from .base import _bam_to_rdkit_2d
        param = self._next_param()
        if param is not None:
            idx = int(param) % len(self.molecules)
        else:
            if self.seed is not None:
                np.random.seed(self.seed)
            idx = int(np.random.choice(len(self.molecules), p=self.p))
        bam_mol = self.molecules[idx]
        context = Context()
        context.bam_molecule = bam_mol
        context.molecule = _bam_to_rdkit_2d(bam_mol)
        return context


class Compound(ChainableBlock):
    """
    Always returns the same molecule (optionally copying it).
    ``n_params = 0`` — no stochastic element, nothing to optimise.

    Parameters
    ----------
    molecule : Molecule or None
        The molecule to return.  May be ``None`` when *param* is given, in
        which case the molecule must be supplied via a keyword argument each
        time the pipeline is called.
    copy : bool
        Whether to copy the molecule before returning it (default ``True``).
    param : str, optional
        Named-input key.  When set, the pipeline call ``p(name=mol)`` binds
        *mol* to this block.  If an input is provided it takes precedence over
        *molecule*; if no input is provided the block falls back to *molecule*.
        Raises ``ValueError`` when both *molecule* is ``None`` and no input
        is bound at call time.

    Examples
    --------
    >>> scaffold = mf.Compound(core_mol, param="core")
    >>> p = scaffold | mf.FindLinkerAtoms() | mf.Connect(side_chain) | mf.Forge()
    >>> ctx = p(core=my_other_mol)   # overrides core_mol for this call
    >>> ctx = p()                    # falls back to core_mol
    """

    n_params = 0

    def __init__(self, molecule, copy: bool = True, param: str = None):
        self.molecule = molecule
        self.copy = copy
        self.param = param

    def _resolve_molecule(self):
        """Return the molecule to use, checking for a named-input override."""
        mol = self.molecule
        if self.param is not None:
            injected = self._get_input(self.param)
            if injected is not None:
                mol = injected
        if mol is None:
            if self.param:
                raise ValueError(
                    f"Compound param={self.param!r}: no molecule set and no input "
                    f"was bound. Call the pipeline as pipe({self.param}=molecule)."
                )
            raise ValueError("Compound has no molecule.")
        return mol

    def _call_buildamol_native(self, *args, **kwargs) -> "Context":
        mol = self._resolve_molecule()
        context = Context()
        context.molecule = mol.copy() if self.copy else mol
        return context

    def _call_rdkit_accelerated(self, *args, **kwargs) -> "Context":
        from .base import _bam_to_rdkit_2d
        mol = self._resolve_molecule()
        context = Context()
        context.bam_molecule = mol
        context.molecule = _bam_to_rdkit_2d(mol)
        return context


class Input(ChainableBlock):
    """
    A named input slot whose molecule is provided at pipeline-call time.

    Use this when a position in the pipeline is *always* supplied by the
    caller rather than having a built-in default::

        p = mf.Input("R1") | mf.FindLinkerAtoms() | mf.Connect(side) | mf.Forge()
        ctx = p(R1=my_molecule)

    Raises ``ValueError`` if the pipeline is called without the required
    keyword argument.

    Parameters
    ----------
    name : str
        The keyword argument name the caller must supply.
    copy : bool
        Whether to copy the bound molecule before returning it (default ``True``).

    See Also
    --------
    Compound : Use ``Compound(default_mol, param="name")`` when a fallback
               molecule should be used if no input is injected.
    """

    n_params = 0

    def __init__(self, name: str, copy: bool = True):
        self.name = name
        self.copy = copy

    def _resolve_molecule(self):
        mol = self._get_input(self.name)
        if mol is None:
            raise ValueError(
                f"Input {self.name!r}: no molecule was provided. "
                f"Call the pipeline as pipe({self.name}=molecule)."
            )
        return mol

    def _call_buildamol_native(self, *args, **kwargs) -> "Context":
        mol = self._resolve_molecule()
        context = Context()
        context.molecule = mol.copy() if self.copy else mol
        return context

    def _call_rdkit_accelerated(self, *args, **kwargs) -> "Context":
        from .base import _bam_to_rdkit_2d
        mol = self._resolve_molecule()
        context = Context()
        context.bam_molecule = mol
        context.molecule = _bam_to_rdkit_2d(mol)
        return context
