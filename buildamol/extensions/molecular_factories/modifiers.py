import random as _random

import numpy as np

from .base import ChainableBlock, Context
import inspect


class Modify(ChainableBlock):
    """
    A modifier is a block that modifies a molecule. It has a `__call__` method that takes a molecule and returns a modified molecule.
    """

    def __init__(self, func, mol=None):
        self.func = func
        self.mol_type = mol

        # check if func needs an additional "atom" or "at_atom" argument aside from the mol
        sig = inspect.signature(func)
        self._needs_atom = any(
            p.name in ("atom", "at_atom") and p.default is inspect.Parameter.empty
            for p in sig.parameters.values()
        )

    def _call_default(self, context: Context, *args, **kwargs) -> Context:
        from rdkit.Chem import rdchem
        from buildamol.core import Molecule
        from .base import _bam_to_rdkit_2d, _rdkit_to_bam_with_embed

        mol = context.molecule
        is_bam = isinstance(mol, Molecule)
        is_rdkit = isinstance(mol, rdchem.Mol)

        if self.mol_type == "rdkit" and is_bam:
            converted = mol.to_rdkit()
            result = self._apply_func(converted, context, *args, **kwargs)
            context.molecule = Molecule.from_rdkit(result) if isinstance(result, rdchem.Mol) else result
        elif self.mol_type == "buildamol" and is_rdkit:
            converted = _rdkit_to_bam_with_embed(mol)
            result = self._apply_func(converted, context, *args, **kwargs)
            context.molecule = _bam_to_rdkit_2d(result) if isinstance(result, Molecule) else result
        else:
            context.molecule = self._apply_func(mol, context, *args, **kwargs)

        return context

    def _apply_func(self, mol, context, *args, **kwargs):
        if self._needs_atom:
            return self.func(mol, context.linker_atom, *args, **kwargs)
        return self.func(mol, *args, **kwargs)


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

    def _call_default(self, context: Context, *args, **kwargs) -> Context:
        param = self._next_param()
        apply = (param > 0.5) if param is not None else (_random.random() < self.p)

        if apply:
            return self.block(context, *args, **kwargs)

        # Consume the inner block's params so downstream blocks stay aligned.
        for _ in range(self.block.total_params()):
            self._next_param()
        return context


class SetAttachResidue(ChainableBlock):
    """
    Sets the attach residue of a molecule. The attach residue is the residue that will be used to attach other molecules.
    """

    def __init__(self, residue):
        self.residue = residue

    def _call_buildamol_native(self, context: Context, *args, **kwargs) -> Context:
        if callable(self.residue):
            attach_residue = self.residue(context.molecule, *args, **kwargs)
        else:
            attach_residue = context.molecule.get_residue(self.residue)
        context.molecule.set_attach_residue(attach_residue)
        context.attach_residue = context.molecule.get_attach_residue()
        return context

    def _call_rdkit_accelerated(self, context: Context, *args, **kwargs) -> Context:
        # Store the raw spec so SetLinkerAtoms can apply it to any bam mol it creates.
        context.attach_residue = self.residue
        bam_mol = context.bam_molecule
        if bam_mol is not None:
            if callable(self.residue):
                attach_residue = self.residue(bam_mol, *args, **kwargs)
            else:
                attach_residue = bam_mol.get_residue(self.residue)
            bam_mol.set_attach_residue(attach_residue)
        return context


class SetLinkerAtoms(ChainableBlock):
    """
    Sets the linker and deleter atoms on the context.

    ``mol`` controls what molecule type the ``link`` / ``delete`` callables
    receive:

    ``"auto"`` (default)
        Native type of the active backend — a BuildAMol ``Molecule`` in
        *buildamol-native*, an RDKit ``Mol`` in *rdkit-accelerated*.
    ``"buildamol"``
        Always pass a BuildAMol ``Molecule``.  In *rdkit-accelerated* mode a
        lightweight 2-D molecule is created on-the-fly (cheap, no 3-D embed)
        unless ``needs_conformer=True``, in which case a full ETKDGv3 embed
        is performed instead.
    ``"rdkit"``
        Always pass an RDKit ``Mol``.  In *buildamol-native* mode the bam
        molecule is stripped of Hs and conformers before being handed to the
        callable.

    ``needs_conformer`` only takes effect when ``mol="buildamol"`` in the
    *rdkit-accelerated* backend.  Set it to ``True`` when the callable needs
    physically meaningful 3-D coordinates (distance lookups, angular
    constraints …).

    In *rdkit-accelerated* mode the ``delete`` argument is always ignored
    — Hs are implicit and ``Connect`` handles valence automatically.
    """

    def __init__(self, link=None, delete=None, mol="auto", needs_conformer=False):
        self.link = link
        self.delete = delete
        self.mol_type = mol
        self.needs_conformer = needs_conformer

    # ── mol-type resolution helpers ───────────────────────────────────────────

    def _effective_mol_type(self, context):
        """Resolve ``"auto"`` to ``"buildamol"`` or ``"rdkit"`` from context."""
        if self.mol_type != "auto":
            return self.mol_type
        from rdkit.Chem import rdchem
        return "rdkit" if isinstance(context.molecule, rdchem.Mol) else "buildamol"

    def _get_working_mol(self, context):
        """Return *(working_mol, mol_kind, bam_ref)*.

        *working_mol* — what callables receive.
        *mol_kind*    — ``"buildamol"`` or ``"rdkit"``.
        *bam_ref*     — the bam mol used to back-map atom indices when
                        ``context.molecule`` is an RDKit Mol.
        """
        from rdkit.Chem import rdchem
        from buildamol.core import Molecule

        mol_type = self._effective_mol_type(context)
        ctx_mol = context.molecule

        if mol_type == "buildamol":
            if isinstance(ctx_mol, Molecule):
                return ctx_mol, "buildamol", ctx_mol
            # rdkit context → produce a bam mol
            bam_mol = context.bam_molecule  # use cached original if available
            if bam_mol is None:
                from .base import _rdkit_to_bam_lightweight
                bam_mol = _rdkit_to_bam_lightweight(ctx_mol, self.needs_conformer)
                # Fresh mol — propagate attach_residue spec from context so
                # callables that use mol.get_attach_residue() see the right residue.
                spec = context.attach_residue
                if spec is not None:
                    if callable(spec):
                        bam_mol.set_attach_residue(spec(bam_mol))
                    else:
                        bam_mol.set_attach_residue(spec)
            return bam_mol, "buildamol", bam_mol
        else:  # "rdkit"
            if isinstance(ctx_mol, rdchem.Mol):
                return ctx_mol, "rdkit", None
            # bam context → produce a lightweight rdkit mol for the callable
            from .base import _bam_to_rdkit_2d
            return _bam_to_rdkit_2d(ctx_mol), "rdkit", ctx_mol

    def _resolve_atom(self, val, working_mol, mol_kind, *extra):
        """Resolve *val* (callable / int / str) to a raw atom in *working_mol*.

        *extra* is forwarded to callables as positional arguments after the mol
        (used for the delete callable which receives the linker atom).
        """
        if callable(val):
            return val(working_mol, *extra)
        if isinstance(val, int):
            if mol_kind == "buildamol":
                residue = getattr(working_mol, "get_attach_residue", lambda: None)()
                return working_mol.get_atom(val, residue=residue)
            return val  # rdkit atom index — pass through
        # string atom name / id
        if mol_kind == "buildamol":
            residue = getattr(working_mol, "get_attach_residue", lambda: None)()
            return working_mol.get_atom(val, residue=residue)
        raise ValueError(
            f"SetLinkerAtoms: cannot resolve atom {val!r} by name in rdkit mode. "
            "Use mol='buildamol' or pass a callable."
        )

    def _to_context_linker(self, raw, mol_kind, bam_ref, context):
        """Map a raw callable result to the linker representation for *context.molecule*.

        * bam Atom  → rdkit int index  (when context.molecule is an RDKit Mol)
        * RDKit Atom → rdkit int index (when context.molecule is an RDKit Mol)
        * bam Atom  → bam Atom         (when context.molecule is a bam Molecule)
        * int        → int              (no conversion needed)
        """
        from rdkit.Chem import rdchem
        from buildamol.core import Molecule

        if isinstance(raw, int):
            if isinstance(context.molecule, Molecule):
                # int index → look up heavy atom in bam mol
                ref = bam_ref or context.molecule
                heavy = [a for a in ref.get_atoms() if a.element.upper() != "H"]
                return heavy[raw]
            return raw  # rdkit context: int is already a valid index

        if isinstance(raw, rdchem.Atom):
            idx = raw.GetIdx()
            if isinstance(context.molecule, Molecule):
                ref = bam_ref or context.molecule
                heavy = [a for a in ref.get_atoms() if a.element.upper() != "H"]
                return heavy[idx]
            return idx  # rdkit context: just the index

        # bam Atom
        if isinstance(context.molecule, rdchem.Mol):
            # map bam atom → rdkit index via heavy-atom list of bam_ref
            heavy = [a for a in bam_ref.get_atoms() if a.element.upper() != "H"]
            return heavy.index(raw)
        return raw  # bam context: return bam Atom directly

    # ── unified dispatch ──────────────────────────────────────────────────────

    def _call_default(self, context: Context, *args, **kwargs) -> Context:
        """Handles both backends — delete is suppressed in *rdkit-accelerated*."""
        from .backend import get_backend
        is_rdkit_mode = get_backend() == "rdkit-accelerated"

        working_mol, mol_kind, bam_ref = self._get_working_mol(context)
        linker_raw = None  # raw result in working_mol's space (used by delete)

        if self.link is not None:
            linker_raw = self._resolve_atom(self.link, working_mol, mol_kind)
            context.linker_atom = self._to_context_linker(linker_raw, mol_kind, bam_ref, context)

        if is_rdkit_mode:
            context.deleter_atom = None  # implicit Hs — Connect handles valence
        elif self.delete is not None:
            # Delete callables receive (mol, linker_in_working_mol) so they can
            # navigate from the linker, e.g.:
            #   _CH = lambda mol, atom: atom.get_left_hydrogen()
            del_raw = self._resolve_atom(
                self.delete, working_mol, mol_kind,
                linker_raw if linker_raw is not None else context.linker_atom,
            )
            context.deleter_atom = self._to_context_linker(del_raw, mol_kind, bam_ref, context)

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

    # ── named strategies (buildamol-native) ───────────────────────────────────

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

    # ── buildamol-native call ─────────────────────────────────────────────────

    def _call_buildamol_native(self, context: Context, *args, **kwargs) -> Context:
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

    # ── RDKit candidate helpers ───────────────────────────────────────────────

    def _rdkit_allowed_indices(self, context):
        """Return the set of rdkit atom indices in the attach_residue, or None.

        Uses the bam_molecule heavy-atom list as the index bridge: heavy atom *i*
        in the bam mol equals rdkit atom *i* in the H-stripped mol produced by
        ``_bam_to_rdkit_2d``.  After ``Connect``, source atoms still occupy
        their original indices 0..n1-1, so the mapping stays valid.
        """
        bam_mol = context.bam_molecule
        if bam_mol is None:
            return None
        attach_res = bam_mol.get_attach_residue()
        if attach_res is None:
            spec = context.attach_residue
            if spec is None:
                return None
            attach_res = spec(bam_mol) if callable(spec) else bam_mol.get_residue(spec)
        heavy = [a for a in bam_mol.get_atoms() if a.element.upper() != "H"]
        res_atoms = set(attach_res.get_atoms())
        return {i for i, a in enumerate(heavy) if a in res_atoms}

    def _rdkit_get_candidates(self, mol, allowed_indices=None):
        """H-bearing non-hydrogen atoms in an RDKit Mol, scoped to *allowed_indices*."""
        candidates = [a for a in mol.GetAtoms()
                      if a.GetAtomicNum() != 1 and a.GetTotalNumHs() > 0]
        if allowed_indices is not None:
            candidates = [a for a in candidates if a.GetIdx() in allowed_indices]
        return candidates

    def _rdkit_how_random(self, mol, allowed_indices=None):
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        if not atoms:
            raise ValueError("No H-bearing atoms found.")
        return _random.choice(atoms).GetIdx()

    def _rdkit_how_highest_degree(self, mol, allowed_indices=None):
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return max(atoms, key=lambda a: a.GetDegree()).GetIdx()

    def _rdkit_how_lowest_degree(self, mol, allowed_indices=None):
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return min(atoms, key=lambda a: a.GetDegree()).GetIdx()

    def _rdkit_how_highest_cip(self, mol, allowed_indices=None):
        from rdkit.Chem import AssignStereochemistry
        AssignStereochemistry(mol, cleanIt=True, force=True)
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return max(atoms, key=lambda a: int(a.GetPropsAsDict().get("_CIPRank", 0))).GetIdx()

    def _rdkit_how_lowest_cip(self, mol, allowed_indices=None):
        from rdkit.Chem import AssignStereochemistry
        AssignStereochemistry(mol, cleanIt=True, force=True)
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return min(atoms, key=lambda a: int(a.GetPropsAsDict().get("_CIPRank", 0))).GetIdx()

    def _rdkit_how_highest_cip_and_degree(self, mol, allowed_indices=None):
        from rdkit.Chem import AssignStereochemistry
        AssignStereochemistry(mol, cleanIt=True, force=True)
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return max(atoms, key=lambda a: (int(a.GetPropsAsDict().get("_CIPRank", 0)), a.GetDegree())).GetIdx()

    def _rdkit_how_lowest_cip_and_degree(self, mol, allowed_indices=None):
        from rdkit.Chem import AssignStereochemistry
        AssignStereochemistry(mol, cleanIt=True, force=True)
        atoms = self._rdkit_get_candidates(mol, allowed_indices)
        return min(atoms, key=lambda a: (int(a.GetPropsAsDict().get("_CIPRank", 0)), a.GetDegree())).GetIdx()

    def _rdkit_how_furthest_from_center(self, mol, allowed_indices=None):
        raise NotImplementedError(
            "'furthest_from_center' and 'closest_to_center' require 3D coordinates "
            "and are not available in the 'rdkit-accelerated' backend."
        )

    def _rdkit_how_closest_to_center(self, mol, allowed_indices=None):
        raise NotImplementedError(
            "'furthest_from_center' and 'closest_to_center' require 3D coordinates "
            "and are not available in the 'rdkit-accelerated' backend."
        )

    def _rdkit_how_amine_N(self, mol, allowed_indices=None):
        atoms = [a for a in mol.GetAtoms()
                 if a.GetAtomicNum() == 7 and a.GetTotalNumHs() > 0]
        if allowed_indices is not None:
            atoms = [a for a in atoms if a.GetIdx() in allowed_indices]
        if not atoms:
            raise ValueError("No amine nitrogen with H found.")
        return atoms[0].GetIdx()

    def _rdkit_how_hydroxyl_O(self, mol, allowed_indices=None):
        atoms = [a for a in mol.GetAtoms()
                 if a.GetAtomicNum() == 8 and a.GetTotalNumHs() > 0]
        if allowed_indices is not None:
            atoms = [a for a in atoms if a.GetIdx() in allowed_indices]
        if not atoms:
            raise ValueError("No hydroxyl oxygen with H found.")
        return atoms[0].GetIdx()

    def _rdkit_how_carbonyl_C(self, mol, allowed_indices=None):
        from rdkit.Chem import MolFromSmarts
        patt = MolFromSmarts("[CH](=O)")
        matches = mol.GetSubstructMatches(patt)
        if allowed_indices is not None:
            matches = [m for m in matches if m[0] in allowed_indices]
        if not matches:
            raise ValueError("No carbonyl C-H found.")
        return matches[0][0]

    def _rdkit_how_carboxyl_O(self, mol, allowed_indices=None):
        from rdkit.Chem import MolFromSmarts
        patt = MolFromSmarts("[OH][C](=O)")
        matches = mol.GetSubstructMatches(patt)
        if allowed_indices is not None:
            matches = [m for m in matches if m[0] in allowed_indices]
        if not matches:
            raise ValueError("No carboxyl OH oxygen found.")
        return matches[0][0]

    def _rdkit_how_carboxyl_C(self, mol, allowed_indices=None):
        from rdkit.Chem import MolFromSmarts
        patt = MolFromSmarts("[C](=O)[OH]")
        matches = mol.GetSubstructMatches(patt)
        if allowed_indices is not None:
            matches = [m for m in matches if m[0] in allowed_indices]
        if not matches:
            raise ValueError("No carboxyl C found.")
        return matches[0][0]

    # ── rdkit-accelerated call ────────────────────────────────────────────────

    def _call_rdkit_accelerated(self, context: Context, *args, **kwargs) -> Context:
        mol = context.molecule
        allowed_indices = self._rdkit_allowed_indices(context)

        if self._raw_how is None:
            atoms = self._rdkit_get_candidates(mol, allowed_indices)
            if not atoms:
                raise ValueError("No H-bearing atoms found in molecule.")
            param = self._next_param()
            if param is not None:
                idx = int(param * len(atoms)) % len(atoms)
                linker_idx = atoms[idx].GetIdx()
            else:
                linker_idx = _random.choice(atoms).GetIdx()
        elif callable(self._raw_how):
            result = self.how(mol, *args, **kwargs)
            linker_idx = result[0] if isinstance(result, tuple) else result
            if hasattr(linker_idx, "GetIdx"):
                linker_idx = linker_idx.GetIdx()
        else:
            rdkit_method_name = f"_rdkit_how_{self._raw_how}"
            rdkit_method = getattr(self, rdkit_method_name, None)
            if rdkit_method is None:
                raise NotImplementedError(
                    f"Strategy '{self._raw_how}' has no rdkit-accelerated implementation."
                )
            linker_idx = rdkit_method(mol, allowed_indices)

        context.linker_atom = linker_idx
        context.deleter_atom = None  # implicit Hs — handled by Connect
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

    def _call_buildamol_native(self, context: Context, *args, **kwargs) -> Context:
        if context.linker_atom is None:
            raise ValueError(
                "Linker atom is not set on the current context. "
                "Use FindLinkerAtoms or SetLinkerAtoms before Connect."
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

    def _call_rdkit_accelerated(self, context: Context, *args, **kwargs) -> Context:
        if context.linker_atom is None:
            raise ValueError(
                "Linker atom is not set on the current context. "
                "Use FindLinkerAtoms or SetLinkerAtoms before Connect."
            )

        other_context = self.other_block()

        if other_context.linker_atom is None:
            raise ValueError(
                "The connected pipeline did not set a linker atom."
            )

        from rdkit.Chem import RWMol, BondType, SanitizeMol
        from rdkit import Chem

        mol1 = context.molecule
        mol2 = other_context.molecule
        n1 = mol1.GetNumAtoms()

        combined = Chem.CombineMols(mol1, mol2)
        rw = RWMol(combined)

        # Consume one explicit H from each linker atom before forming the bond.
        # If RemoveHs stored removed Hs as numExplicitHs (which counts toward
        # explicit valence), not decrementing causes overvalence after AddBond.
        for _idx in (context.linker_atom, other_context.linker_atom + n1):
            _at = rw.GetAtomWithIdx(_idx)
            _nh = _at.GetNumExplicitHs()
            if _nh > 0:
                _at.SetNumExplicitHs(_nh - 1)

        rw.AddBond(context.linker_atom, other_context.linker_atom + n1, Chem.BondType.SINGLE)
        Chem.SanitizeMol(rw)

        out = Context()
        out.molecule = rw.GetMol()
        # Preserve the bam reference for the accumulating molecule so that
        # downstream SetLinkerAtoms calls can still resolve bam-style atom
        # references. Core atoms occupy indices 0..n1-1 in the combined mol —
        # the same indices as in the original bam_molecule rdkit conversion —
        # so the mapping remains valid even after the connection.
        out.bam_molecule = context.bam_molecule
        return out


class Forge(ChainableBlock):
    """
    Terminal pipeline block that materialises the molecule into a BuildAMol Molecule.

    In the ``buildamol-native`` backend the molecule is already assembled; Forge
    is a no-op unless ``optimize=True``, in which case ``mol.optimize()`` is called.

    In the ``rdkit-accelerated`` backend Forge embeds a 3-D conformer via RDKit's
    ETKDGv3, optionally minimises with MMFF, and converts to a BuildAMol Molecule.

    Parameters
    ----------
    optimize : bool
        Run a force-field minimisation after materialisation (default ``False``).
    """

    n_params = 0

    def __init__(self, optimize: bool = False):
        self.optimize = optimize

    def _call_buildamol_native(self, context: Context, *args, **kwargs) -> Context:
        if self.optimize:
            context.molecule.optimize(inplace=True)
        return context

    def _call_rdkit_accelerated(self, context: Context, *args, **kwargs) -> Context:
        from .base import _rdkit_to_bam_with_embed
        context.molecule = _rdkit_to_bam_with_embed(context.molecule, optimize=self.optimize)
        return context
