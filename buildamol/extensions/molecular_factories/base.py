import threading


def _bam_to_rdkit_2d(bam_mol):
    """Convert a BuildAMol Molecule to a lightweight RDKit Mol (no conformer, implicit Hs)."""
    from rdkit import Chem
    rdmol = bam_mol.to_rdkit()
    # molecule_to_rdkit sets noImplicit=True on every atom (all Hs are explicit
    # bonds in the bam structure). RemoveHs on such a mol stores removed H atoms
    # as numExplicitHs on the parent atom instead of as implicit Hs. That means
    # atoms carry explicit-H counts that make Connect overestimate the valence
    # (e.g. C ends up with bond-valence 4 + numExplicitHs 1 = 5). Clearing the
    # flag first lets RemoveHs recalculate a proper implicit-H count instead.
    rw = Chem.RWMol(rdmol)
    for atom in rw.GetAtoms():
        atom.SetNoImplicit(False)
    rdmol = Chem.RemoveHs(rw.GetMol())
    rdmol.RemoveAllConformers()
    return rdmol


def _rdkit_to_bam_with_embed(rdmol, optimize: bool = True):
    """Embed a conformer-less RDKit Mol and convert to a BuildAMol Molecule."""
    from rdkit.Chem import AllChem
    from buildamol.core import Molecule
    rdmol_h = AllChem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol_h, AllChem.ETKDGv3())
    if optimize:
        AllChem.MMFFOptimizeMolecule(rdmol_h)
    return Molecule.from_rdkit(rdmol_h)


def _rdkit_to_bam_lightweight(rdmol, needs_conformer: bool = False):
    """Convert an RDKit Mol to a BuildAMol Molecule without a full 3-D embed.

    Adds explicit Hs and computes a cheap 2-D layout (all z = 0), then
    delegates to ``Molecule.from_rdkit``.  The resulting Molecule supports
    all graph-based operations (``get_atom``, element/neighbour searches,
    constraint matching) but coordinates are 2-D projections, so
    distance/angle queries are not physically meaningful.

    When *needs_conformer* is ``True`` a proper ETKDGv3 conformer is generated
    instead, enabling geometry-dependent callables (distance lookups, angular
    constraints …) at the cost of a full embed.
    """
    if needs_conformer:
        return _rdkit_to_bam_with_embed(rdmol)
    from rdkit.Chem import AllChem
    from buildamol.core import Molecule
    rdmol_h = AllChem.AddHs(rdmol)
    AllChem.Compute2DCoords(rdmol_h)
    return Molecule.from_rdkit(rdmol_h)


class ChainableBlock:
    """
    A chainable block is a block that can be chained to other blocks.

    Optimisation interface
    ----------------------
    Stochastic blocks declare ``n_params`` (int) and ``param_bounds``
    (list of ``(lo, hi)`` pairs, one per param).  The ``Optimizer`` discovers
    the full parameter space by walking the pipeline tree via ``total_params()``
    and ``param_bounds``, then injects a numeric vector for each candidate via
    ``_inject_params()``.  Blocks pull values from the thread-local iterator
    with ``_next_param()``; when no vector is injected they fall back to their
    default (usually random) behaviour.
    """

    # ── thread-local param iterator ───────────────────────────────────────────
    _param_local = threading.local()

    @classmethod
    def _inject_params(cls, values):
        """Inject a numeric parameter vector for the current thread."""
        cls._param_local.iter = iter(values)

    @classmethod
    def _clear_params(cls):
        """Remove the injected vector so blocks revert to stochastic mode."""
        cls._param_local.iter = None

    @classmethod
    def _next_param(cls):
        """Return the next value from the injected vector, or None."""
        it = getattr(cls._param_local, "iter", None)
        if it is None:
            return None
        try:
            return next(it)
        except StopIteration:
            return None

    # ── parameter-space discovery ─────────────────────────────────────────────
    n_params: int = 0

    def total_params(self) -> int:
        """Total number of optimisable scalar parameters in this block (and its sub-pipelines)."""
        return self.n_params

    @property
    def param_bounds(self) -> list:
        """List of ``(lo, hi)`` bounds, one per optimisable parameter."""
        return []

    # ── backend dispatch ──────────────────────────────────────────────────────
    def __call__(self, *args, **kwargs):
        from .backend import get_backend
        key = get_backend().replace("-", "_")
        method = getattr(self, f"_call_{key}", None) or getattr(self, "_call_default", None)
        if method is None:
            raise NotImplementedError(
                f"{type(self).__name__} does not implement backend "
                f"{get_backend()!r} and has no _call_default."
            )
        return method(*args, **kwargs)

    # ── chaining ──────────────────────────────────────────────────────────────
    def __or__(self, other):
        if isinstance(other, ChainableBlock):
            return ChainedBlock(self, other)
        raise TypeError(f"Cannot chain {type(self)} with {type(other)}")


class ChainedBlock(ChainableBlock):
    """
    The result of chaining two blocks together with ``|``.
    Calls the first block, then passes its output to the second.
    """

    def __init__(self, first, second):
        self.first = first
        self.second = second

    def __call__(self, *args, **kwargs):
        return self.second(self.first(*args, **kwargs))

    def total_params(self) -> int:
        return self.first.total_params() + self.second.total_params()

    @property
    def param_bounds(self) -> list:
        return self.first.param_bounds + self.second.param_bounds


class Context:
    """
    The context passed from one block to another.
    """

    def __init__(self):
        self.molecule = None
        self.attach_residue = None
        self.linker_atom = None
        self.deleter_atom = None
        self.backend = None
        self.bam_molecule = None  # original bam mol; set by sources in rdkit-accelerated mode


class MultiContext:
    def __init__(self):
        self.contexts = []

    def add_context(self, context):
        self.contexts.append(context)

    def __getattr__(self, name):
        return [getattr(context, name) for context in self.contexts]
