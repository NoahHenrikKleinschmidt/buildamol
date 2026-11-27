"""
Tests for automatic reactive site detection and AutoReact
"""

import pytest
import buildamol as bam
from buildamol.structural.reactivity import (
    detect_sites,
    order_reactivities,
    default_reactivities,
    Reactivity,
    ReactionError,
    SiteMap,
)


def _sites(smiles, include_ch=False):
    mol = bam.read_smiles(smiles)
    sites = detect_sites(mol)
    if not include_ch:
        sites = [s for s in sites if s.reactivity != "CH"]
    return mol, {(s.reactivity, s.role, s.bonder.element.upper()) for s in sites}


def _formula(mol):
    counts = {}
    for a in mol.get_atoms():
        e = a.element.upper()
        counts[e] = counts.get(e, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# ordering
# ---------------------------------------------------------------------------


def test_order_reactivities_respects_supersedes():
    ordered = [r.name for r in order_reactivities(default_reactivities())]
    for r in default_reactivities():
        for s in r.supersedes:
            assert ordered.index(r.name) < ordered.index(s), (r.name, s)
    assert ordered[-1] == "CH"


def test_order_reactivities_detects_cycles():
    class A(Reactivity):
        centers = ("C",)
        supersedes = ("B",)

    class B(Reactivity):
        centers = ("C",)
        supersedes = ("A",)

    with pytest.raises(ValueError):
        order_reactivities([A(), B()])


# ---------------------------------------------------------------------------
# detection and claiming
# ---------------------------------------------------------------------------


def test_detect_carboxyl_claims_hydroxyl():
    _, sites = _sites("CC(=O)O")
    assert sites == {
        ("Carboxyl", "electrophile", "C"),
        ("Carboxyl", "nucleophile", "O"),
    }


def test_detect_amide_is_not_amine():
    _, sites = _sites("CC(=O)N")
    assert sites == {("Amide", "nucleophile", "N")}


def test_detect_ketone_is_not_aldehyde():
    _, sites = _sites("CC(=O)C")
    assert sites == {("Ketone", "electrophile", "C")}
    _, sites = _sites("CC=O")
    assert sites == {("Aldehyde", "electrophile", "C")}


def test_detect_ester_is_not_carboxyl():
    _, sites = _sites("CC(=O)OC")
    assert sites == {("Ester", "electrophile", "C")}


def test_detect_lactone_has_no_leaving_group():
    _, sites = _sites("O=C1CCCO1")
    assert sites == set()


def test_detect_phenol_is_nucleophile_only():
    _, sites = _sites("Oc1ccccc1")
    assert sites == {("Hydroxyl", "nucleophile", "O")}


def test_detect_multifunctional():
    _, sites = _sites("NC(CS)C(=O)O")  # cysteine
    assert sites == {
        ("Amine", "nucleophile", "N"),
        ("Thiol", "nucleophile", "S"),
        ("Carboxyl", "electrophile", "C"),
        ("Carboxyl", "nucleophile", "O"),
    }


def test_detect_ch_fallback_only_on_unclaimed_carbons():
    mol = bam.read_smiles("CC=O")
    ch = [s for s in detect_sites(mol) if s.reactivity == "CH"]
    # the aldehyde carbon is claimed, only the methyl remains
    assert {s.bonder.id for s in ch} == {"C1"}

    mol = bam.read_smiles("c1ccccc1")
    ch = [s for s in detect_sites(mol) if s.reactivity == "CH"]
    assert len({s.bonder for s in ch}) == 6


def test_detect_restricted_to_residue():
    mol = bam.read_smiles("CC(=O)O")
    res = mol.get_residue(1)
    assert detect_sites(mol, residues=res) == detect_sites(mol)


def test_detect_custom_linker_based_reactivity():

    class Nitrogen(Reactivity):
        def nucleophile_linker(self, mol):
            return mol.get_atoms("N", by="element")

    mol = bam.read_smiles("CCN")
    sites = detect_sites(mol, reactivities=[Nitrogen()])
    assert [(s.reactivity, s.role, s.bonder.element.upper()) for s in sites] == [
        ("Nitrogen", "nucleophile", "N")
    ]


# ---------------------------------------------------------------------------
# AutoReact
# ---------------------------------------------------------------------------


def test_autoreact_prefers_amine_over_hydroxyl():
    acid = bam.read_smiles("CC(=O)O")
    aminoalcohol = bam.read_smiles("OCCCN")
    best = bam.AutoReact().best(acid, aminoalcohol)
    assert best.rule.label == "amide"
    assert best.nucleophile.reactivity == "Amine"

    out = bam.AutoReact()(acid, aminoalcohol)
    # C2H4O2 + C3H9NO - H2O
    assert _formula(out) == {"C": 5, "H": 11, "N": 1, "O": 2}


def test_autoreact_benzene_aminopropanol():
    benzene = bam.read_smiles("c1ccccc1")
    aminopropanol = bam.read_smiles("OCCCN")
    R = bam.AutoReact()
    best = R.best(benzene, aminopropanol)
    assert best.rule.fallback
    assert best.nucleophile.reactivity == "Amine"
    out = R(benzene, aminopropanol)
    # C6H6 + C3H9NO - H2
    assert _formula(out) == {"C": 9, "H": 13, "N": 1, "O": 1}


def test_autoreact_imine_sets_double_bond():
    aldehyde = bam.read_smiles("CC=O")
    amine = bam.read_smiles("NCC")
    out = bam.AutoReact()(aldehyde, amine)
    assert _formula(out) == {"C": 4, "H": 9, "N": 1}
    N = out.get_atoms("N", by="element")[0]
    orders = [out.get_bond(N, n).order for n in out.get_neighbors(N)]
    assert sorted(orders) == [1, 2]


def test_autoreact_works_in_both_directions():
    acid = bam.read_smiles("CC(=O)O")
    amine = bam.read_smiles("NCC")
    R = bam.AutoReact()
    forward = R.best(acid, amine)
    backward = R.best(amine, acid)
    assert forward.target_is_electrophile
    assert not backward.target_is_electrophile
    assert forward.rule == backward.rule


def test_autoreact_is_deterministic():
    A = bam.read_smiles("NCC(=O)O")
    B = bam.read_smiles("NCC(=O)O")
    R = bam.AutoReact()
    runs = [
        [(c.target_site.bonder.id, c.source_site.bonder.id) for c in R.candidates(A, B)]
        for _ in range(3)
    ]
    assert runs[0] == runs[1] == runs[2]


def test_autoreact_does_not_modify_inputs():
    A = bam.read_smiles("CC(=O)O")
    B = bam.read_smiles("NCC")
    nA, nB = A.count_atoms(), B.count_atoms()
    bam.AutoReact()(A, B)
    assert A.count_atoms() == nA and B.count_atoms() == nB


def test_autoreact_without_fallback_raises():
    A = bam.read_smiles("c1ccccc1")
    B = bam.read_smiles("c1ccccc1")
    R = bam.AutoReact(allow_fallback=False)
    assert R.candidates(A, B) == []
    with pytest.raises(ReactionError):
        R(A, B)


def test_autoreact_custom_pairs():
    acid = bam.read_smiles("CC(=O)O")
    aminoalcohol = bam.read_smiles("OCCCN")
    R = bam.AutoReact(pairs={("Hydroxyl", "Carboxyl"): (2.0, 1, "ester")})
    assert R.best(acid, aminoalcohol).rule.label == "ester"

    R = bam.AutoReact(pairs={("Amine", "Carboxyl"): None})
    assert R.best(acid, aminoalcohol).rule.label != "amide"


def test_autoreact_explain():
    text = bam.AutoReact().explain(bam.read_smiles("CC(=O)O"), bam.read_smiles("NCC"))
    assert "amide" in text.splitlines()[1]


def test_autoreact_same_molecule():
    A = bam.read_smiles("NCC(=O)O")
    # copies are distinct molecules, so a molecule can react with itself
    out = bam.AutoReact()(A, A)
    assert out.count_atoms() == 2 * A.count_atoms() - 3
    # intramolecular reactions are not supported
    with pytest.raises(ValueError):
        bam.AutoReact()(A, A, copy_target=False, copy_source=False)


# ---------------------------------------------------------------------------
# caching
# ---------------------------------------------------------------------------


def _site_keys(sites):
    return [
        (s.reactivity, s.role, s.bonder, s.deletes, frozenset(s.footprint), s.score)
        for s in sites
    ]


_MONOMERS = [
    "NCC(=O)O",
    "NC(CS)C(=O)O",
    "OCCCN",
    "c1ccccc1",
    "CC(=O)OC",
    "BrCCO",
    "CC=O",
    "Oc1ccc(N)cc1",
    "OP(=O)(O)OCC",
]


@pytest.mark.parametrize("copy_target", [False, True])
def test_cached_sites_match_fresh_detection(copy_target):
    monomers = [bam.read_smiles(s) for s in _MONOMERS]
    R = bam.AutoReact()
    poly = monomers[0].copy()
    for i in range(12):
        poly = R(poly, monomers[(i * 5) % len(monomers)], copy_target=copy_target)
        assert SiteMap.lookup(poly) is not None
        assert _site_keys(R.sites(poly)) == _site_keys(detect_sites(poly))


def test_site_map_is_updated_not_rebuilt(monkeypatch):
    mono = bam.read_smiles("NCC(=O)O")
    R = bam.AutoReact()
    poly = mono.copy()
    R(poly, mono, copy_target=False)  # caches poly and mono

    calls = []
    original = SiteMap.__init__

    def counting_init(self, mol=None, *args, **kwargs):
        if mol is not None:
            calls.append(mol)
        original(self, mol, *args, **kwargs)

    monkeypatch.setattr(SiteMap, "__init__", counting_init)
    for _ in range(5):
        poly = R(poly, mono, copy_target=False)
    assert calls == []


def test_site_map_copy_to():
    mol = bam.read_smiles("NC(CS)C(=O)O")
    smap = SiteMap(mol)
    copy = mol.copy()
    moved = smap.copy_to(mol, copy)
    assert moved.is_valid(copy)
    assert _site_keys(moved.sites()) == _site_keys(detect_sites(copy))


def _hydroxyl_O(mol):
    return next(
        o
        for o in mol.get_atoms("O", by="element")
        if any(n.element.upper() == "H" for n in mol.get_neighbors(o))
    )


def test_site_map_update_matches_fresh_detection():
    mol = bam.read_smiles("NCC(=O)O")
    other = bam.read_smiles("NCC(=O)O")
    smap = SiteMap(mol)

    OH = _hydroxyl_O(mol)
    C = next(iter(mol.get_neighbors(OH) - {a for a in mol.get_neighbors(OH) if a.element.upper() == "H"}))
    HO = next(a for a in mol.get_neighbors(OH) if a.element.upper() == "H")
    N = other.get_atoms("N", by="element")[0]
    HN = next(a for a in other.get_neighbors(N) if a.element.upper() == "H")
    link = bam.linkage(C.id, N.id, [OH.id, HO.id], [HN.id])

    before = {id(r) for r in mol.get_residues()}
    out = bam.connect(mol, other, link, copy_a=False, copy_b=False, use_patch=False)
    assert out is mol
    seeds = [C] + [a for r in mol.get_residues() if id(r) not in before for a in r.get_atoms()]
    smap.update(mol, seeds)
    assert smap.is_valid(mol)
    assert _site_keys(smap.sites()) == _site_keys(detect_sites(mol))
    assert "Amide" in {s.reactivity for s in smap.sites()}


def test_cache_is_refreshed_after_external_modification():
    mol = bam.read_smiles("CC(=O)O")
    R = bam.AutoReact()
    assert "Carboxyl" in {s.reactivity for s in R.sites(mol)}
    # turn the acid into an aldehyde behind the cache's back
    OH = _hydroxyl_O(mol)
    HO = next(a for a in mol.get_neighbors(OH) if a.element.upper() == "H")
    mol.remove_atoms(OH, HO)
    assert _site_keys(R.sites(mol)) == _site_keys(detect_sites(mol))
    assert "Carboxyl" not in {s.reactivity for s in R.sites(mol)}


def test_cache_does_not_keep_molecules_alive():
    import gc
    import weakref

    R = bam.AutoReact()
    mol = bam.read_smiles("NCC(=O)O")
    R.sites(mol)
    ref = weakref.ref(mol)
    del mol
    gc.collect()
    assert ref() is None


def test_cache_is_not_copied_or_pickled():
    import pickle

    mol = bam.read_smiles("NCC(=O)O")
    SiteMap.cached(mol)
    assert SiteMap.lookup(mol) is not None
    assert SiteMap.lookup(mol.copy()) is None
    assert SiteMap.lookup(pickle.loads(pickle.dumps(mol))) is None


def test_candidates_limit_is_prefix_of_full_ranking():
    A = bam.read_smiles("NC(CS)C(=O)O")
    B = bam.read_smiles("OCCCN")
    R = bam.AutoReact()
    full = R.candidates(A, B)
    for n in (1, 2, 5, 20):
        assert R.candidates(A, B, n=n) == full[:n]


# ---------------------------------------------------------------------------
# steric tie-breaking
# ---------------------------------------------------------------------------


def _candidate_keys(candidates):
    return [
        (
            c.rule.label,
            c.target_site.bonder.serial_number,
            c.source_site.bonder.serial_number,
            c.target_is_electrophile,
            c.steric,
        )
        for c in candidates
    ]


def test_steric_ranking_is_independent_of_conformation():
    import numpy as np

    A = bam.read_smiles("NCC(CO)C(=O)O")
    B = bam.read_smiles("OCCCN")
    R = bam.AutoReact()
    ranking = _candidate_keys(R.candidates(A, B))

    rng = np.random.default_rng(0)
    for atom in A.get_atoms():
        atom.coord = rng.normal(scale=5.0, size=3)
    assert _candidate_keys(R.candidates(A, B)) == ranking


def test_steric_prefers_less_substituted_site():
    # a primary and a tertiary bromide have the same score; the primary one is less hindered
    target = bam.read_smiles("BrCCC(C)(C)Br")
    best = bam.AutoReact().best(target, bam.read_smiles("SCC"))
    assert best.rule.label == "thioether"
    bonder = best.target_site.bonder
    heavy = [n for n in target.get_neighbors(bonder) if n.element.upper() != "H"]
    assert len(heavy) == 2  # the CH2 carrying a bromine


def test_assembly_is_reproducible():
    smiles = ["NC(CS)C(=O)O", "OCCCN", "c1ccccc1", "BrCCO", "NCC(=O)O"]

    def build():
        monomers = [bam.read_smiles(s) for s in smiles]
        R = bam.AutoReact()
        poly = monomers[0].copy()
        for i in range(8):
            poly = R(poly, monomers[(i * 3 + 1) % len(monomers)])
        return sorted(
            (b.atom1.serial_number, b.atom2.serial_number, b.order)
            for b in poly.get_bonds()
        )

    assert build() == build()


# ---------------------------------------------------------------------------
# operators
# ---------------------------------------------------------------------------


def _bond_signature(mol):
    return sorted(
        (b.atom1.serial_number, b.atom2.serial_number, b.order) for b in mol.get_bonds()
    )


def test_add_without_linkage_uses_autoreact():
    acid = bam.read_smiles("CC(=O)O")
    amine = bam.read_smiles("OCCCN")
    nA, nB = acid.count_atoms(), amine.count_atoms()

    out = acid + amine
    assert out is not acid and out is not amine
    assert acid.count_atoms() == nA and amine.count_atoms() == nB
    assert _formula(out) == {"C": 5, "H": 11, "N": 1, "O": 2}
    assert _bond_signature(out) == _bond_signature(bam.autoreact(acid, amine))


def test_iadd_without_linkage_is_inplace():
    acid = bam.read_smiles("CC(=O)O")
    amine = bam.read_smiles("NCC")
    ref = acid + amine
    mol = acid.copy()
    alias = mol
    mol += amine
    assert mol is alias
    assert _bond_signature(mol) == _bond_signature(ref)


def test_mul_and_repeat_without_linkage():
    gly = bam.read_smiles("NCC(=O)O")
    n = gly.count_atoms()
    # each amide bond releases one water
    tri = gly * 3
    assert tri.count_atoms() == 3 * n - 2 * 3
    assert _formula(tri) == {"C": 6, "H": 11, "N": 3, "O": 4}

    rep = gly.repeat(3, inplace=False)
    assert _bond_signature(rep) == _bond_signature(tri)
    assert gly.count_atoms() == n


def test_explicit_linkage_wins_over_autoreact():
    # with a linkage the old behaviour is used: an ether between the two hydroxyls
    # instead of the ester AutoReact would pick
    acid = bam.read_smiles("OCC(=O)O")
    alcohol = bam.read_smiles("OCC")
    assert bam.AutoReact().best(acid, alcohol).rule.label == "ester"

    O_acid = next(
        o
        for o in acid.get_atoms("O", by="element")
        if any(n.element.upper() == "H" for n in acid.get_neighbors(o))
        and not any(
            n.element.upper() == "C"
            and any(m.element.upper() == "O" and m is not o for m in acid.get_neighbors(n))
            for n in acid.get_neighbors(o)
        )
    )
    C_acid = next(n for n in acid.get_neighbors(O_acid) if n.element.upper() == "C")
    H_acid = next(n for n in acid.get_neighbors(O_acid) if n.element.upper() == "H")
    O_alc = alcohol.get_atoms("O", by="element")[0]
    H_alc = next(n for n in alcohol.get_neighbors(O_alc) if n.element.upper() == "H")
    link = bam.linkage(C_acid.id, O_alc.id, [O_acid.id, H_acid.id], [H_alc.id])

    out = (acid % link) + alcohol
    C = out.get_atom(C_acid.id, residue=1)
    assert any(
        n.element.upper() == "O" and n.parent is not C.parent for n in out.get_neighbors(C)
    )


def test_attach_residue_restricts_autoreact():
    # a dipeptide: the free amine sits in residue 1, residue 2 only has the carboxyl and amide N
    di = bam.read_smiles("NCC(=O)O") + bam.read_smiles("NCC(=O)O")
    assert len(di.residues) == 2
    acid = bam.read_smiles("CC(=O)O")

    R = bam.AutoReact()
    free = R.best(acid, di)
    assert free.rule.label == "amide"
    assert free.source_site.bonder.parent.id[1] == 1

    restricted = R.best(acid, di, residues_source=2)
    assert restricted.source_site.bonder.parent.id[1] == 2
    assert restricted.rule.label != "amide"

    # the operator honours the attach residue
    out_free = acid + di
    out_restricted = acid + (di.copy() @ 2)
    assert _formula(out_free) == {"C": 6, "H": 10, "N": 2, "O": 4}
    assert _bond_signature(out_free) != _bond_signature(out_restricted)


def test_operator_uses_default_instance():
    original = bam.AutoReact.default
    try:
        bam.AutoReact.default = bam.AutoReact(allow_fallback=False)
        with pytest.raises(ReactionError):
            bam.read_smiles("c1ccccc1") + bam.read_smiles("c1ccccc1")
    finally:
        bam.AutoReact.default = original
    out = bam.read_smiles("c1ccccc1") + bam.read_smiles("c1ccccc1")
    assert _formula(out) == {"C": 12, "H": 10}
