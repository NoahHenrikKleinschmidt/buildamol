import threading


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


class MultiContext:
    def __init__(self):
        self.contexts = []

    def add_context(self, context):
        self.contexts.append(context)

    def __getattr__(self, name):
        return [getattr(context, name) for context in self.contexts]
