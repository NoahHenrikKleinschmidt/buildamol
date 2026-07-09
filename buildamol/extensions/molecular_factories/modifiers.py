import random as _random

import numpy as np

from .base import ChainableBlock, Context
import inspect


class Modify(ChainableBlock):
    """
    A modifier is a block that modifies a molecule. It has a `__call__` method that takes a molecule and returns a modified molecule.
    """

    def __init__(self, func):
        self.func = func

        # check if func needs an additional "atom" or "at_atom" argument aside from the mol
        sig = inspect.signature(func)
        self._needs_atom = any(
            p.name in ("atom", "at_atom") and p.default is inspect.Parameter.empty
            for p in sig.parameters.values()
        )

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        if self._needs_atom:
            if context.linker_atom is None:
                raise ValueError(
                    "Linker atom is not set on the current context. "
                    "Use SetLinkerAtoms or FindLinkerAtoms before Modify."
                )
            mol = self.func(context.molecule, context.linker_atom, *args, **kwargs)
        else:
            mol = self.func(context.molecule, *args, **kwargs)
        context.molecule = mol
        return context


class Random(ChainableBlock):
    """
    Randomly applies an inner block with probability *p*.

    Optimisation interface
    ----------------------
    ``n_params = 1`` — the single parameter is a float in ``[0, 1]``.
    The inner block is applied when ``param > 0.5``.

    If the inner block itself has optimisable parameters, ``Random`` includes
    them in its total so the vector length stays consistent.  When the block is
    *not* applied those parameters are still consumed from the iterator so that
    subsequent blocks remain correctly aligned.
    """

    n_params = 1  # the apply/skip trigger

    def __init__(self, block: ChainableBlock, p: float = 0.5):
        self.block = block
        self.p = p

    def total_params(self) -> int:
        return 1 + self.block.total_params()

    @property
    def param_bounds(self) -> list:
        return [(0.0, 1.0)] + self.block.param_bounds

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        param = self._next_param()
        apply = (param > 0.5) if param is not None else (_random.random() < self.p)

        if apply:
            return self.block(context, *args, **kwargs)

        # Consume the inner block's params so downstream blocks stay aligned.
        for _ in range(self.block.total_params()):
            self._next_param()
        return context


class Optimize(ChainableBlock):
    def __init__(self, algorithm=None, **kwargs):
        self.algorithm = algorithm
        self.kwargs = kwargs

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        mol = context.molecule
        mol.optimize(self.algorithm, **self.kwargs, inplace=True)
        return context


class SetAttachResidue(ChainableBlock):
    """
    Sets the attach residue of a molecule. The attach residue is the residue that will be used to attach other molecules.
    """

    def __init__(self, residue):
        self.residue = residue

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        if callable(self.residue):
            attach_residue = self.residue(context.molecule, *args, **kwargs)
        else:
            attach_residue = context.molecule.get_residue(self.residue)
        context.molecule.set_attach_residue(attach_residue)
        context.attach_residue = context.molecule.get_attach_residue()
        return context


class SetLinkerAtoms(ChainableBlock):
    def __init__(self, link=None, delete=None):
        self.link = link
        self.delete = delete

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        # Keep the current linker as default so that `delete`-only calls can
        # still reference it (e.g. `_CH = lambda mol, atom: atom.get_left_hydrogen()`).
        linker_atom = context.linker_atom

        if self.link is not None:
            if callable(self.link):
                linker_atom = self.link(context.molecule, *args, **kwargs)
            else:
                linker_atom = context.molecule.get_atom(
                    self.link, residue=context.molecule.get_attach_residue()
                )
            context.linker_atom = linker_atom

        if self.delete is not None:
            if callable(self.delete):
                # Callable deleters receive (molecule, linker_atom) so they can
                # navigate relative to the chosen linker atom, e.g.:
                #   _CH = lambda mol, atom: atom.get_left_hydrogen()
                deleter_atom = self.delete(
                    context.molecule, linker_atom, *args, **kwargs
                )
            else:
                deleter_atom = context.molecule.get_atom(
                    self.delete, residue=context.molecule.get_attach_residue()
                )
            context.deleter_atom = deleter_atom
        return context


class FindLinkerAtoms(ChainableBlock):
    """
    Finds a linker atom and a hydrogen to delete.

    ``how`` parameter
    -----------------
    ``None`` (default)
        Stochastic in normal mode; **optimisable** — the ``Optimizer`` injects
        a float in ``[0, 1]`` that is mapped to an index into the list of
        H-bearing candidate atoms in the attach residue.  ``n_params = 1``.
    ``"random"``
        Always random; ``n_params = 0`` (fixed stochastic, not optimisable).
    Any named strategy (``"highest_cip"``, …)
        Always uses that strategy; ``n_params = 0``.
    callable
        Called as ``how(molecule) -> (linker_atom, deleter_atom)``; ``n_params = 0``.
    """

    _allowed_how = [
        "random",
        "highest_cip",
        "lowest_cip",
        "highest_degree",
        "lowest_degree",
        "furthest_from_center",
        "closest_to_center",
        "highest_cip_and_degree",
        "lowest_cip_and_degree",
        "amine_N",
        "hydroxyl_O",
        "carbonyl_C",
        "carboxyl_O",
        "carboxyl_C",
    ]

    def __init__(self, how=None):
        self._raw_how = how  # preserve the user's original intent

        if how is None:
            self.how = None  # handled inline in __call__
        elif callable(how):
            self.how = how
        elif how in self._allowed_how:
            self.how = getattr(self, f"_how_{how}")
        else:
            raise ValueError(
                f"Invalid how: {how!r}. Must be one of {self._allowed_how} or a callable."
            )

    @property
    def n_params(self) -> int:
        return 1 if self._raw_how is None else 0

    @property
    def param_bounds(self) -> list:
        # [0, 1] float maps onto the candidate atom list at call time.
        return [(0.0, 1.0)] if self._raw_how is None else []

    # ── shared helpers ────────────────────────────────────────────────────────

    def _get_residue(self, molecule):
        return molecule.get_attach_residue() or molecule.residues[0]

    def _get_candidates(self, molecule):
        """H-bearing atoms in the attach residue; raises if none found."""
        residue = self._get_residue(molecule)
        atoms = [a for a in residue.get_atoms() if a.get_hydrogens()]
        if not atoms:
            raise ValueError(
                f"No atoms with hydrogen neighbours found in residue {residue}."
            )
        return atoms

    @staticmethod
    def _pick_h(linker):
        return min(linker.get_hydrogens(), key=lambda h: h.id)

    # ── named strategies ──────────────────────────────────────────────────────

    def _how_random(self, molecule, *args, **kwargs):
        atoms = self._get_candidates(molecule)
        linker = _random.choice(atoms)
        deleter = _random.choice(list(linker.get_hydrogens()))
        return linker, deleter

    def _how_highest_cip(self, molecule, *args, **kwargs):
        from buildamol.structural.chirality import _cip_priority_key
        atoms = self._get_candidates(molecule)
        atoms.sort(key=_cip_priority_key, reverse=True)
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_lowest_cip(self, molecule, *args, **kwargs):
        from buildamol.structural.chirality import _cip_priority_key
        atoms = self._get_candidates(molecule)
        atoms.sort(key=_cip_priority_key)
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_highest_degree(self, molecule, *args, **kwargs):
        atoms = self._get_candidates(molecule)
        atoms.sort(
            key=lambda a: sum(1 for n in a.get_neighbors() if n.element != "H"),
            reverse=True,
        )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_lowest_degree(self, molecule, *args, **kwargs):
        atoms = self._get_candidates(molecule)
        atoms.sort(key=lambda a: sum(1 for n in a.get_neighbors() if n.element != "H"))
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_furthest_from_center(self, molecule, *args, **kwargs):
        atoms = self._get_candidates(molecule)
        center = molecule.center_of_mass
        atoms.sort(key=lambda a: np.linalg.norm(a.coord - center), reverse=True)
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_closest_to_center(self, molecule, *args, **kwargs):
        atoms = self._get_candidates(molecule)
        center = molecule.center_of_mass
        atoms.sort(key=lambda a: np.linalg.norm(a.coord - center))
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_highest_cip_and_degree(self, molecule, *args, **kwargs):
        from buildamol.structural.chirality import _cip_priority_key
        atoms = self._get_candidates(molecule)
        atoms.sort(
            key=lambda a: (
                _cip_priority_key(a),
                sum(1 for n in a.get_neighbors() if n.element != "H"),
            ),
            reverse=True,
        )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_lowest_cip_and_degree(self, molecule, *args, **kwargs):
        from buildamol.structural.chirality import _cip_priority_key
        atoms = self._get_candidates(molecule)
        atoms.sort(
            key=lambda a: (
                _cip_priority_key(a),
                sum(1 for n in a.get_neighbors() if n.element != "H"),
            )
        )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_amine_N(self, molecule, *args, **kwargs):
        residue = self._get_residue(molecule)
        atoms = [a for a in residue.get_atoms() if a.element == "N" and a.get_hydrogens()]
        if not atoms:
            raise ValueError(
                f"No amine nitrogen atoms with hydrogen neighbours found in residue {residue}."
            )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_hydroxyl_O(self, molecule, *args, **kwargs):
        from buildamol.structural import constraints_v2 as c
        residue = self._get_residue(molecule)
        cc = c.and_(c.has_element("O"), c.has_bond_of_order_with(1, "H"))
        atoms = [a for a in residue.get_atoms() if cc(a)]
        if not atoms:
            raise ValueError(
                f"No hydroxyl oxygen atoms with hydrogen neighbours found in residue {residue}."
            )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_carbonyl_C(self, molecule, *args, **kwargs):
        from buildamol.structural import constraints_v2 as c
        residue = self._get_residue(molecule)
        cc = c.and_(
            c.has_element("C"),
            c.has_bond_of_order_with(2, "O"),
            c.has_bond_of_order_with(1, "H"),
        )
        atoms = [a for a in residue.get_atoms() if cc(a)]
        if not atoms:
            raise ValueError(
                f"No carbonyl carbon atoms with hydrogen neighbours found in residue {residue}."
            )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_carboxyl_O(self, molecule, *args, **kwargs):
        from buildamol.structural import constraints_v2 as c
        residue = self._get_residue(molecule)
        c_cc = c.and_(
            c.has_element("C"),
            c.has_bond_of_order_with(1, "O"),
            c.has_bond_of_order_with(2, "O"),
        )
        c_atoms = [a for a in residue.get_atoms() if c_cc(a)]
        if not c_atoms:
            raise ValueError(
                f"No carboxyl group found in residue {residue}."
            )
        oc = c.and_(
            c.has_element("O"),
            c.has_bond_of_order_with(1, "C"),
            c.has_bond_of_order_with(1, "H"),
        )
        atoms = [a for a in c_atoms[0].get_neighbors() if oc(a)]
        if not atoms:
            raise ValueError(
                f"No carboxyl OH oxygen found in residue {residue}."
            )
        linker = atoms[0]
        return linker, self._pick_h(linker)

    def _how_carboxyl_C(self, molecule, *args, **kwargs):
        from buildamol.structural import constraints_v2 as c
        residue = self._get_residue(molecule)
        cc = c.and_(
            c.has_element("C"),
            c.has_bond_of_order_with(1, "O"),
            c.has_bond_of_order_with(2, "O"),
        )
        atoms = [a for a in residue.get_atoms() if cc(a)]
        if not atoms:
            raise ValueError(
                f"No carboxyl carbon found in residue {residue}."
            )
        linker = atoms[0]
        oc = c.and_(
            c.has_element("O"),
            c.has_bond_of_order_with(1, "C"),
            c.has_bond_of_order_with(1, "H"),
        )
        deleter = next((a for a in linker.get_neighbors() if oc(a)), None)
        return linker, deleter

    # ── call ─────────────────────────────────────────────────────────────────

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        if self._raw_how is None:
            # Optimisable: map a [0,1] float onto the candidate list.
            atoms = self._get_candidates(context.molecule)
            param = self._next_param()
            if param is not None:
                idx = int(param * len(atoms)) % len(atoms)
                linker = atoms[idx]
                deleter = self._pick_h(linker)
            else:
                linker = _random.choice(atoms)
                deleter = _random.choice(list(linker.get_hydrogens()))
        else:
            linker, deleter = self.how(context.molecule, *args, **kwargs)

        context.linker_atom = linker
        context.deleter_atom = deleter
        return context


class Connect(ChainableBlock):
    n_params = 0

    def __init__(self, other_block):
        self.other_block = other_block

    def total_params(self) -> int:
        return self.other_block.total_params()

    @property
    def param_bounds(self) -> list:
        return self.other_block.param_bounds

    def __call__(self, context: Context, *args, **kwargs) -> Context:
        if context.linker_atom is None:
            raise ValueError(
                "Linker atom is not set on the current context. "
                "Use SetLinkerAtoms or FindLinkerAtoms before Connect."
            )

        # Run the other (source) pipeline independently to obtain the fragment
        # to attach. Source pipelines typically start with Choice or Compound
        # and do not need the accumulating context as input.
        other_context = self.other_block()

        if other_context.linker_atom is None:
            raise ValueError(
                "The connected pipeline did not set a linker atom. "
                "Ensure FindLinkerAtoms or SetLinkerAtoms is included in that pipeline."
            )

        from buildamol import core

        # Pass string IDs, not Atom objects, into the delete lists.
        # The linkage validation and the Stitcher work with copies of the source
        # molecule (copy_b=True), so original Atom objects would not match
        # the copy's atoms — causing StopIteration in the Stitcher and
        # "Linkage cannot be applied" in can_be_source.
        # When deleter_atom is None we pass None to let the linkage auto-detect
        # a hydrogen on the anchor atom (valid behaviour).
        delete_in_target = (
            [getattr(context.deleter_atom, "id", context.deleter_atom)]
            if context.deleter_atom is not None
            else None
        )
        delete_in_source = (
            [getattr(other_context.deleter_atom, "id", other_context.deleter_atom)]
            if other_context.deleter_atom is not None
            else None
        )

        link = core.linkage(
            atom1=context.linker_atom,
            atom2=other_context.linker_atom,
            delete_in_target=delete_in_target,
            delete_in_source=delete_in_source,
        )

        # Derive residue identifiers from the linker atoms so that
        # can_be_target / can_be_source search the correct residue.
        # Without this, both fall back to `attach_residue or -1` (last residue),
        # which fails when the linker is in the first residue of a multi-residue
        # source molecule (e.g. aryl+R2 where aryl is residue 0).
        # at_residue_a: Residue object is safe because copy_a=False (same model).
        # at_residue_b: seqid integer, copy-safe (copies preserve seqids).
        target_residue = context.linker_atom.parent
        source_seqid = other_context.linker_atom.parent.id[1]

        out = Context()
        # copy_a=False: keep the accumulating molecule in-place so that any
        # Atom object references held by the caller remain valid after this step.
        # copy_b=True (default): always work on a fresh copy of the source.
        out.molecule = core.connect(
            context.molecule,
            other_context.molecule,
            link,
            at_residue_a=target_residue,
            at_residue_b=source_seqid,
            copy_a=False,
            copy_b=True,
        )
        return out
