"""
Automatic, chemistry-grounded connection of two molecules.

`AutoReact` detects reactive sites in both molecules (see `buildamol.structural.reactivity.detect_sites`),
pairs nucleophilic and electrophilic sites using a pairing table, ranks all possible pairs and
connects the molecules using the highest ranking pair.

Example
-------
>>> import buildamol as bam
>>> benzene = bam.read_smiles("c1ccccc1")
>>> aminopropanol = bam.read_smiles("OCCCN")
>>> out = bam.AutoReact()(benzene, aminopropanol)   # N-arylation (amine beats hydroxyl)
>>> print(bam.AutoReact().explain(benzene, aminopropanol)) # explanation of the chosen pair
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple, Union

import buildamol.core as core
import buildamol.base_classes as base_classes
from buildamol.structural.reactivity import (
    Site,
    ReactionError,
    SiteMap,
    detect_sites,
    order_reactivities,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PairRule:
    """
    Rule to define a pair between nucleophile and electrophile reactivities.
    It defines how well a nucleophile reacts with an electrophile.

    Parameters
    ----------
    score : float
        The preference for this pairing (0-1, higher is preferred). This dominates the ranking of candidates.
    bond_order : int
        The order of the newly formed bond (e.g. 2 for imine formation).
    label : str
        A short human-readable name of the resulting linkage (e.g. "amide").
    fallback : bool
        Whether this is a non-mechanistic fallback pairing (excluded with `AutoReact(allow_fallback=False)`).
    """

    score: float
    bond_order: int = 1
    label: str = ""
    fallback: bool = False


def _rules(*entries) -> Dict[Tuple[str, str], PairRule]:
    return {(n, e): PairRule(*rest) for n, e, *rest in entries}


DEFAULT_PAIRS: Dict[Tuple[str, str], PairRule] = _rules(
    # (nucleophile, electrophile, score, bond_order, label, is_fallback[optional])
    ("Amine", "Carboxyl", 1.00, 1, "amide"),
    ("Amine", "Ester", 0.80, 1, "amide (aminolysis)"),
    ("Amine", "Aldehyde", 0.85, 2, "imine"),
    ("Amine", "Ketone", 0.75, 2, "imine"),
    ("Amine", "AlkylHalide", 0.70, 1, "N-alkylation"),
    ("Amine", "Hydroxyl", 0.30, 1, "N-alkylation (alcohol)"),
    ("Amine", "CH", 0.25, 1, "C-N", True),
    ("Thiol", "AlkylHalide", 0.90, 1, "thioether"),
    ("Thiol", "Carboxyl", 0.60, 1, "thioester"),
    ("Thiol", "CH", 0.20, 1, "C-S", True),
    ("Hydroxyl", "Carboxyl", 0.80, 1, "ester"),
    ("Hydroxyl", "Ester", 0.50, 1, "ester (transesterification)"),
    ("Hydroxyl", "AlkylHalide", 0.60, 1, "ether"),
    ("Hydroxyl", "Hydroxyl", 0.30, 1, "ether (dehydration)"),
    ("Hydroxyl", "CH", 0.20, 1, "C-O", True),
    ("Carboxyl", "AlkylHalide", 0.55, 1, "ester (carboxylate alkylation)"),
    ("Phosphate", "Hydroxyl", 0.60, 1, "phosphoester"),
    ("Phosphate", "AlkylHalide", 0.50, 1, "phosphoester"),
    ("Amide", "AlkylHalide", 0.40, 1, "N-alkyl amide"),
    ("CH", "AlkylHalide", 0.30, 1, "C-C (coupling)"),
    ("CH", "Carboxyl", 0.30, 1, "ketone (acylation)"),
    ("CH", "CH", 0.10, 1, "C-C", True),  # fallback for C-C coupling
)
"""
The default pairing table: (nucleophile name, electrophile name) -> PairRule.
"""


@dataclass(frozen=True)
class Candidate:
    """
    A possible connection between two molecules.

    Attributes
    ----------
    score : float
        The total score (pair rule score + 0.1 * (nucleophile site score + electrophile site score)).
    rule : PairRule
        The pairing rule that applies.
    nucleophile : Site
        The nucleophilic site.
    electrophile : Site
        The electrophilic site.
    target_is_electrophile : bool
        Whether the electrophile is in the target (first) molecule.
    steric : int
        The number of heavy atoms within a few bonds of both bonders (only used to break ties).
    """

    score: float
    rule: PairRule
    nucleophile: Site
    electrophile: Site
    target_is_electrophile: bool
    steric: int = 0

    @property
    def target_site(self) -> Site:
        return self.electrophile if self.target_is_electrophile else self.nucleophile

    @property
    def source_site(self) -> Site:
        return self.nucleophile if self.target_is_electrophile else self.electrophile

    def deletes(self, site: Site) -> tuple:
        """The atoms that are actually deleted at a site for this candidate's bond order"""
        if site.leaving_hydrogens:
            return site.deletes[: self.rule.bond_order]
        return site.deletes

    def __repr__(self):
        t, s = self.target_site, self.source_site
        return (
            f"Candidate({self.rule.label or 'link'}, score={self.score:.3f}, "
            f"target={t.reactivity}:{t.bonder.id}, source={s.reactivity}:{s.bonder.id})"
        )


def _heavy_atoms_within_bonds(adj, atom, n_bonds: int) -> int:
    """The number of non-hydrogen atoms within `n_bonds` bonds of an atom (excluding the atom itself)"""
    seen = {atom}
    frontier = [atom]
    count = 0
    for _ in range(n_bonds):
        nxt = []
        for a in frontier:
            for n in adj[a]:
                if n in seen:
                    continue
                seen.add(n)
                if n.element.upper() != "H":
                    count += 1
                    # hydrogens are terminal, only heavy atoms need to be expanded
                    nxt.append(n)
        if not nxt:
            break
        frontier = nxt
    return count


class AutoReact:
    """
    Automatically connect two molecules in a chemically (at least somewhat) reasonable and reproducible way.

    Reactive sites are detected in both molecules, all nucleophile/electrophile pairs between the two
    (in both directions) are noted and scored using the pairing table; the best pair is then used to connect
    the molecules. The ranking is deterministic: ties are broken by steric accessibility, which is measured on the bond connectivity graph and is thus independent of the molecule's current 3D conformer.
    In case of steric ties, the target is chosen to be the electrophile, if still tied then finally by atom serial numbers.

    Parameters
    ----------
    pairs : dict, optional
        By default the built-in pairs are used (see `DEFAULT_PAIRS`). This can be used to overwrite existing pairs or add new ones.
        Pair values can be set to None to ignore them. Tuples of `(score, bond_order, label, fallback)` are accepted instead of a `PairRule` object.
    reactivities : list of Reactivity, optional
        The reactivities to consider when detecting reactive sites. By default `buildamol.structural.reactivity.default_reactivities()`.
    allow_fallback : bool
        Whether fallback pairings (e.g. C-H to C-H) may be used in case of no suitable mechanistic pairings. True by default, only disable if you have a specific reason to do so.
    replace_pairs : bool
        If True, `pairs` replaces the default table instead of updating it.
        Use this for fully custom pairing rules or to restrict the search space to a subset of the default pairs for performance considerations.
    steric_bonds : int
        The number of bonds around a bonder within which heavy atoms are counted for steric tie-breaking (default 3).
        This is a topological measure of how substituted a site is (e.g. a primary carbon takes precedence over a tertiary one).
    cache : bool
        Whether to cache detected sites on the molecules (see `SiteMap`). Cached sites are transferred onto copies
        and updated incrementally after each connection, so iterative assembly (e.g. `poly = R(poly, monomer)` in a loop)
        only re-scans the atoms around each new bond. Cached sites are re-detected automatically if a molecule's
        bonds were changed in the meantime (detected by atom/bond counts and bond orders). Caches are maintained per molecule and are not automatically copied when using `.copy()`!
    """

    def __init__(
        self,
        pairs: dict = None,
        reactivities: list = None,
        allow_fallback: bool = True,
        replace_pairs: bool = False,
        steric_bonds: int = 3,
        cache: bool = True,
    ):
        self.pairs = {} if replace_pairs else dict(DEFAULT_PAIRS)
        if pairs:
            for key, rule in pairs.items():
                if rule is None:
                    self.pairs.pop(key, None)
                    continue
                if not isinstance(rule, PairRule):
                    rule = PairRule(*rule)
                self.pairs[key] = rule
        self.reactivities = (
            order_reactivities(reactivities) if reactivities is not None else None
        )
        self.allow_fallback = allow_fallback
        self.steric_bonds = steric_bonds
        self.cache = cache

    def site_map(self, mol: "core.Molecule") -> SiteMap:
        """
        Get the (cached) `SiteMap` of a molecule. It is re-detected if the molecule changed since it was cached.

        Parameters
        ----------
        mol : Molecule

        Returns
        -------
        SiteMap
        """
        if self.cache:
            return SiteMap.cached(mol, self.reactivities, ordered=True)
        return SiteMap(mol, self.reactivities, ordered=True)

    def _copy(self, mol):
        new = mol.copy()
        if self.cache:
            self.site_map(mol).copy_to(mol, new).store(new)
        return new

    def sites(self, mol: "core.Molecule", residues=None) -> List[Site]:
        """
        Detect the reactive sites of a molecule.

        Parameters
        ----------
        mol : Molecule
        residues : Residue or list, optional
            Restrict detection to these residues.

        Returns
        -------
        list of Site
        """
        if residues is None and self.cache:
            return self.site_map(mol).sites()
        return detect_sites(mol, self.reactivities, residues=residues, ordered=True)

    def candidates(
        self,
        target: "core.Molecule",
        source: "core.Molecule",
        residues_target=None,
        residues_source=None,
        n: int = None,
    ) -> List[Candidate]:
        """
        Rank all possible connections between two molecules (best first).

        Parameters
        ----------
        target : Molecule
            The target (first) molecule.
        source : Molecule
            The source (second) molecule.
        residues_target, residues_source : Residue or list, optional
            Restrict site detection to these residues.
        n : int, optional
            Only return the `n` best candidates. This is faster for large molecules since
            steric tie-breaking is only computed where it can affect the top `n`.

        Returns
        -------
        list of Candidate
        """
        limit = n
        t_sites = self.sites(target, residues_target)
        s_sites = self.sites(source, residues_source)

        # group the sites by reactivity so that only pairs with a rule are visited
        groups = []
        for t_is_e, e_sites, n_sites in (
            (True, t_sites, s_sites),
            (False, s_sites, t_sites),
        ):
            electrophiles = {}
            nucleophiles = {}
            for s in e_sites:
                if s.role == "electrophile":
                    electrophiles.setdefault(s.reactivity, []).append(s)
            for s in n_sites:
                if s.role == "nucleophile":
                    nucleophiles.setdefault(s.reactivity, []).append(s)
            for (n_name, e_name), rule in self.pairs.items():
                if rule.fallback and not self.allow_fallback:
                    continue
                ns = nucleophiles.get(n_name)
                es = electrophiles.get(e_name)
                if not ns or not es:
                    continue
                bound = rule.score + 0.1 * (
                    max(s.score for s in ns) + max(s.score for s in es)
                )
                groups.append((bound, rule, ns, es, t_is_e))

        # visit the most promising pairs first; if only the best `limit` candidates are needed,
        # pairs that cannot reach the current limit-th best score are skipped
        groups.sort(key=lambda g: -g[0])
        raw = []
        for bound, rule, ns, es, t_is_e in groups:
            if limit is not None and len(raw) >= limit:
                threshold = sorted((r[0] for r in raw), reverse=True)[limit - 1]
                if bound < threshold - 1e-9:
                    break
            bo = rule.bond_order
            for n in ns:
                if n.leaving_hydrogens and len(n.deletes) < bo:
                    continue
                for e in es:
                    if e.leaving_hydrogens and len(e.deletes) < bo:
                        continue
                    score = rule.score + 0.1 * (n.score + e.score)
                    raw.append((score, rule, n, e, t_is_e))

        if not raw:
            return []

        # the score dominates the ranking, sterics only break ties. So if only the best n are
        # requested, the sterics are only needed up to the score group that contains the n-th candidate
        raw.sort(key=lambda r: -round(r[0], 9))
        if limit is not None and limit < len(raw):
            last = round(raw[limit - 1][0], 9)
            k = limit
            while k < len(raw) and round(raw[k][0], 9) == last:
                k += 1
            raw = raw[:k]

        # steric counts only break ties
        t_bonders = {}
        s_bonders = {}
        for _, _, n, e, t_is_e in raw:
            t_site, s_site = (e, n) if t_is_e else (n, e)
            t_bonders[t_site.bonder] = None
            s_bonders[s_site.bonder] = None
        for mol, bonders in ((target, t_bonders), (source, s_bonders)):
            adj = mol._AtomGraph.adj
            for a in bonders:
                bonders[a] = _heavy_atoms_within_bonds(adj, a, self.steric_bonds)

        candidates = []
        for score, rule, n, e, t_is_e in raw:
            t_site, s_site = (e, n) if t_is_e else (n, e)
            steric = t_bonders[t_site.bonder] + s_bonders[s_site.bonder]
            candidates.append(Candidate(score, rule, n, e, t_is_e, steric))

        candidates.sort(
            key=lambda c: (
                -round(c.score, 9),
                c.steric,
                not c.target_is_electrophile,
                c.target_site.bonder.serial_number,
                c.source_site.bonder.serial_number,
            )
        )
        return candidates if limit is None else candidates[:limit]

    def best(self, target, source, residues_target=None, residues_source=None):
        """
        Get the best connection candidate between two molecules.

        Raises
        ------
        ReactionError
            If no candidate exists.
        """
        candidates = self.candidates(
            target, source, residues_target, residues_source, n=1
        )
        if not candidates:
            raise ReactionError(
                "No compatible reactive sites were found between the two molecules"
            )
        return candidates[0]

    def explain(self, target, source, n: int = 10, **kwargs) -> str:
        """
        Get a human-readable ranking of the best connection candidates.

        Parameters
        ----------
        target, source : Molecule
        n : int
            The number of candidates to show.

        Returns
        -------
        str
        """
        candidates = self.candidates(target, source, **kwargs)
        if not candidates:
            return "No compatible reactive sites found."
        lines = [
            f"{'#':>3}  {'linkage':<30} {'score':>6}  {'target site':<28} {'source site':<28} {'steric':>6}"
        ]
        for i, c in enumerate(candidates[:n]):
            t, s = c.target_site, c.source_site
            label = c.rule.label + (" [fallback]" if c.rule.fallback else "")
            lines.append(
                f"{i:>3}  {label:<30} {c.score:>6.3f}  "
                f"{(t.reactivity + ' ' + t.role[:4] + ' ' + t.bonder.id):<28} "
                f"{(s.reactivity + ' ' + s.role[:4] + ' ' + s.bonder.id):<28} {c.steric:>6}"
            )
        if len(candidates) > n:
            lines.append(f"... and {len(candidates) - n} more")
        return "\n".join(lines)

    def apply(
        self,
        target: "core.Molecule",
        source: "core.Molecule",
        candidate: Union[Candidate, int] = 0,
        copy_target: bool = True,
        copy_source: bool = True,
        residues_target=None,
        residues_source=None,
    ) -> "core.Molecule":
        """
        Connect two molecules using the best (or a chosen) candidate.

        Parameters
        ----------
        target : Molecule
            The target (first) molecule.
        source : Molecule
            The source (second) molecule.
        candidate : int
            The rank of the candidate to use (0 = best).
        copy_target, copy_source : bool
            Whether to copy the molecules before connecting.
        residues_target, residues_source : Residue or list, optional
            Restrict site detection to these residues.

        Returns
        -------
        Molecule
            The connected molecule.
        """
        if target is source and not (copy_target or copy_source):
            raise ValueError(
                "AutoReact connects two separate molecules, got the same molecule twice without copying (use .copy())"
            )
        if copy_target:
            target = self._copy(target)
        if copy_source:
            source = self._copy(source)

        if isinstance(candidate, Candidate):
            raise TypeError(
                "Pass the rank of the candidate (int); candidates refer to atoms of specific molecule copies."
            )
        candidates = self.candidates(
            target,
            source,
            residues_target,
            residues_source,
            n=candidate + 1 if candidate >= 0 else None,
        )
        if not candidates:
            raise ReactionError(
                "No compatible reactive sites were found between the two molecules"
            )
        chosen = candidates[candidate]
        t_site, s_site = chosen.target_site, chosen.source_site
        logger.info(
            "AutoReact: %s (target %s %s %s -> source %s %s %s, score %.3f)",
            chosen.rule.label or "link",
            t_site.reactivity,
            t_site.role,
            t_site.bonder.id,
            s_site.reactivity,
            s_site.role,
            s_site.bonder.id,
            chosen.score,
        )

        smap = SiteMap.lookup(target, self.reactivities) if self.cache else None
        if smap is not None:
            # keep the residues referenced so that their ids cannot be reused by new residues
            residues_before = list(target.get_residues())
            before = {id(r) for r in residues_before}

        out = self._connect(target, source, chosen)

        if smap is not None:
            if out is target:
                # only the target's bonder changed its bonds; everything else that changed is new
                # (the source's atoms are merged as new residues)
                seeds = [chosen.target_site.bonder]
                for r in out.get_residues():
                    if id(r) not in before:
                        seeds.extend(r.get_atoms())
                smap.update(out, seeds).store(out)
        return out

    def _connect(self, target, source, candidate: Candidate):
        t_site, s_site = candidate.target_site, candidate.source_site
        atom1, atom2 = t_site.bonder, s_site.bonder
        link = core.linkage(
            atom1.id,
            atom2.id,
            delete_in_target=[a.id for a in candidate.deletes(t_site)],
            delete_in_source=[a.id for a in candidate.deletes(s_site)],
            description=candidate.rule.label,
        )
        out = core.connect(
            target,
            source,
            link,
            at_residue_a=atom1.parent,
            at_residue_b=atom2.parent,
            copy_a=False,
            copy_b=False,
            use_patch=False,
        )
        if candidate.rule.bond_order > 1:
            # the source atoms may have been re-parented during the merge,
            # so we resolve the new partner through the target atom's neighbors
            partner = next(
                (
                    n
                    for n in out.get_neighbors(atom1)
                    if n.id == atom2.id and n.element == atom2.element
                ),
                None,
            )
            if partner is None:
                raise ReactionError("Could not find the newly formed bond")
            # the leaving atoms were already chosen to match the bond order
            out.set_bond_order(atom1, partner, candidate.rule.bond_order)
        return out

    def __call__(self, target, source, **kwargs):
        return self.apply(target, source, **kwargs)

    def __repr__(self):
        return (
            f"AutoReact(pairs={len(self.pairs)}, allow_fallback={self.allow_fallback})"
        )


AutoReact.default = AutoReact()
"""
The `AutoReact` instance used by `autoreact()` and by the molecule operators (`+`, `+=`, `*`) when no linkage
is defined. Replace it to customize the operators, e.g. `bam.AutoReact.default = bam.AutoReact(allow_fallback=False)`.
"""


def autoreact(target, source, **kwargs) -> "core.Molecule":
    """
    Automatically connect two molecules using the best chemically reasonable linkage.
    Shortcut for `AutoReact.default(target, source, **kwargs)`.

    Parameters
    ----------
    target : Molecule
        The target (first) molecule.
    source : Molecule
        The source (second) molecule.
    kwargs
        Passed to `AutoReact.apply`.

    Returns
    -------
    Molecule
    """
    return AutoReact.default(target, source, **kwargs)


__all__ = ["AutoReact", "PairRule", "Candidate", "DEFAULT_PAIRS", "autoreact"]
