"""
Functions for loading standard fragment libraries for use with the
``Assembler`` and ``Generator`` classes.

Sources
-------
``"chembl"``
    Molecules fetched live from the ChEMBL REST API.
    Requires internet access and ``requests``.

    API-side filters (applied server-side, efficient):

    - ``n`` (int, default 200) — maximum fragments to return
    - ``max_mw`` / ``min_mw`` (float) — molecular weight bounds (Da)
    - ``max_hbd`` / ``min_hbd`` (int) — H-bond donor bounds
    - ``max_hba`` / ``min_hba`` (int) — H-bond acceptor bounds
    - ``max_logp`` / ``min_logp`` (float) — logP bounds
    - ``max_psa`` / ``min_psa`` (float) — polar surface area bounds (Å²)
    - ``max_rotatable_bonds`` (int) — rotatable bond count ceiling
    - ``max_heavy_atoms`` / ``min_heavy_atoms`` (int) — total heavy-atom bounds
    - ``max_aromatic_rings`` (int) — aromatic ring count ceiling
    - ``ro3_pass`` (bool) — ``True`` restricts to Rule-of-Three compliant
      molecules; ``False`` excludes them; ``None`` (default) ignores the flag
    - ``min_qed`` (float 0–1) — minimum QED drug-likeness score
    - ``max_ro5_violations`` (int) — maximum Lipinski Ro5 violations

    Post-filters (applied locally via RDKit after fetching):

    - ``max_carbons`` / ``min_carbons`` (int) — carbon atom count bounds
    - ``allowed_elements`` (set of str) — element symbols that may appear;
      any molecule containing an unlisted element is dropped
      (e.g. ``{"C", "N", "O", "F"}``)
    - ``forbidden_elements`` (set of str) — element symbols that must not
      appear (e.g. ``{"S", "P"}`` to exclude sulfur and phosphorus)
    - ``require_ring`` (bool) — ``True`` keeps only cyclic molecules,
      ``False`` keeps only acyclic ones, ``None`` (default) ignores

``"brics"``
    BRICS retrosynthetic decomposition of a user-supplied molecule set, using
    RDKit's built-in ``rdkit.Chem.BRICS`` module.  No internet needed.
    Key kwargs: ``molecules`` (list of Molecule or SMILES strings),
    ``max_mw`` (200), ``min_heavy_atoms`` (3).

Examples
--------

.. code-block:: python

    from buildamol.extensions.molecular_factories import load_fragment_library

    # Ro3-compliant fragments from ChEMBL
    fragments = load_fragment_library("chembl", n=200, ro3_pass=True)

    # Small aliphatic C/N/O fragments — good for basic building blocks
    fragments = load_fragment_library(
        "chembl",
        n=100,
        max_heavy_atoms=12,
        max_carbons=8,
        allowed_elements={"C", "N", "O"},
        require_ring=False,
    )

    # Fragment a set of known drugs into BRICS building blocks
    import buildamol as bam
    drugs = [bam.Molecule.from_smiles(s) for s in [
        "CC(=O)Oc1ccccc1C(=O)O",     # aspirin
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O", # ibuprofen
    ]]
    fragments = load_fragment_library("brics", molecules=drugs)

    # Both lists are ready for the Assembler/Generator
    from buildamol.extensions.molecular_factories import Generator
    gen = Generator(fragments, my_scoring_fn)
"""

from buildamol.core import Molecule


def make_brics_fragments(molecules, max_mw=200, min_heavy_atoms=3):
    """
    BRICS retrosynthetic decomposition of a molecule set.

    Fragments are obtained by cutting bonds that are likely retrosynthetically
    accessible (using RDKit's BRICS rules), then stripping the dummy atoms
    left by the cuts and keeping only chemically sensible pieces.

    Parameters
    ----------
    molecules : list
        BuildAMol ``Molecule`` objects or SMILES strings to decompose.
    max_mw : float
        Maximum fragment molecular weight (Da). Default 200 keeps fragments
        small enough to serve as true building blocks.
    min_heavy_atoms : int
        Minimum number of heavy (non-H) atoms; filters out trivial
        single-atom or two-atom fragments.

    Returns
    -------
    list of Molecule
    """
    return _from_brics(molecules, max_mw=max_mw, min_heavy_atoms=min_heavy_atoms)


def get_chembl_fragments(
    n=200,
    max_mw=300,
    min_mw=None,
    max_hbd=3,
    min_hbd=None,
    max_hba=3,
    min_hba=None,
    max_logp=3,
    min_logp=None,
    max_psa=None,
    min_psa=None,
    max_rotatable_bonds=None,
    max_heavy_atoms=None,
    min_heavy_atoms=None,
    max_aromatic_rings=None,
    ro3_pass=None,
    min_qed=None,
    max_ro5_violations=None,
    max_carbons=None,
    min_carbons=None,
    allowed_elements=None,
    forbidden_elements=None,
    require_ring=None,
):
    """
    Fetch fragments from the ChEMBL REST API with flexible physicochemical
    and structural filters.

    Parameters
    ----------
    n : int
        Maximum number of fragments to return.
    max_mw, min_mw : float, optional
        Molecular weight bounds (Da).
    max_hbd, min_hbd : int, optional
        H-bond donor bounds.
    max_hba, min_hba : int, optional
        H-bond acceptor bounds.
    max_logp, min_logp : float, optional
        logP bounds.
    max_psa, min_psa : float, optional
        Polar surface area bounds (Å²).
    max_rotatable_bonds : int, optional
        Rotatable bond count ceiling.
    max_heavy_atoms, min_heavy_atoms : int, optional
        Total heavy-atom count bounds.
    max_aromatic_rings : int, optional
        Aromatic ring count ceiling.
    ro3_pass : bool, optional
        ``True`` restricts to Rule-of-Three compliant molecules (MW ≤ 300,
        HBD ≤ 3, HBA ≤ 3, logP ≤ 3); ``False`` excludes them; ``None``
        (default) ignores the flag entirely.
    min_qed : float, optional
        Minimum QED drug-likeness score (0–1).
    max_ro5_violations : int, optional
        Maximum number of Lipinski Ro5 violations.
    max_carbons, min_carbons : int, optional
        Carbon atom count bounds (applied locally via RDKit after fetching).
    allowed_elements : set of str, optional
        Element symbols that may appear.  Any molecule containing an element
        not in this set is dropped (e.g. ``{"C", "N", "O", "F"}``).
    forbidden_elements : set of str, optional
        Element symbols that must not appear (e.g. ``{"S", "P"}``).
    require_ring : bool, optional
        ``True`` keeps only cyclic molecules; ``False`` keeps only acyclic
        ones; ``None`` (default) ignores ring presence.

    Returns
    -------
    list of Molecule
    """
    return _from_chembl(
        n=n,
        max_mw=max_mw,
        min_mw=min_mw,
        max_hbd=max_hbd,
        min_hbd=min_hbd,
        max_hba=max_hba,
        min_hba=min_hba,
        max_logp=max_logp,
        min_logp=min_logp,
        max_psa=max_psa,
        min_psa=min_psa,
        max_rotatable_bonds=max_rotatable_bonds,
        max_heavy_atoms=max_heavy_atoms,
        min_heavy_atoms=min_heavy_atoms,
        max_aromatic_rings=max_aromatic_rings,
        ro3_pass=ro3_pass,
        min_qed=min_qed,
        max_ro5_violations=max_ro5_violations,
        max_carbons=max_carbons,
        min_carbons=min_carbons,
        allowed_elements=allowed_elements,
        forbidden_elements=forbidden_elements,
        require_ring=require_ring,
    )


def load_fragment_library(source="chembl", **kwargs):
    """
    Load a standard fragment library from the given source.

    Parameters
    ----------
    source : str
        ``"chembl"`` or ``"brics"``.
    **kwargs
        Source-specific keyword arguments — see the module docstring or
        :func:`get_chembl_fragments` / :func:`make_brics_fragments` for
        the full parameter reference.

    Returns
    -------
    list of Molecule
        Fragment molecules, autolabeled and ready for use with
        ``Assembler`` or ``Generator``.
    """
    _dispatch = {
        "chembl": _from_chembl,
        "brics": _from_brics,
    }
    if source not in _dispatch:
        raise ValueError(f"Unknown source '{source}'. Available: {list(_dispatch)}")
    return _dispatch[source](**kwargs)


# ------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------


def _from_chembl(
    n=200,
    max_mw=300,
    min_mw=None,
    max_hbd=3,
    min_hbd=None,
    max_hba=3,
    min_hba=None,
    max_logp=3,
    min_logp=None,
    max_psa=None,
    min_psa=None,
    max_rotatable_bonds=None,
    max_heavy_atoms=None,
    min_heavy_atoms=None,
    max_aromatic_rings=None,
    ro3_pass=None,
    min_qed=None,
    max_ro5_violations=None,
    max_carbons=None,
    min_carbons=None,
    allowed_elements=None,
    forbidden_elements=None,
    require_ring=None,
):
    try:
        import requests
    except ImportError:
        raise ImportError("'requests' is required: pip install requests")

    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    base = "https://www.ebi.ac.uk/chembl/api/data/molecule"

    params = {"format": "json"}

    def _set(key, val):
        if val is not None:
            params[key] = val

    _set("molecule_properties__mw_freebase__lte", max_mw)
    _set("molecule_properties__mw_freebase__gte", min_mw)
    _set("molecule_properties__hbd__lte", max_hbd)
    _set("molecule_properties__hbd__gte", min_hbd)
    _set("molecule_properties__hba__lte", max_hba)
    _set("molecule_properties__hba__gte", min_hba)
    _set("molecule_properties__alogp__lte", max_logp)
    _set("molecule_properties__alogp__gte", min_logp)
    _set("molecule_properties__psa__lte", max_psa)
    _set("molecule_properties__psa__gte", min_psa)
    _set("molecule_properties__rtb__lte", max_rotatable_bonds)
    _set("molecule_properties__heavy_atoms__lte", max_heavy_atoms)
    _set("molecule_properties__heavy_atoms__gte", min_heavy_atoms)
    _set("molecule_properties__aromatic_rings__lte", max_aromatic_rings)
    _set("molecule_properties__qed_weighted__gte", min_qed)
    _set("molecule_properties__num_ro5_violations__lte", max_ro5_violations)
    if ro3_pass is not None:
        params["molecule_properties__ro3_pass"] = "Y" if ro3_pass else "N"

    has_local_filters = any(
        v is not None
        for v in (
            max_carbons,
            min_carbons,
            allowed_elements,
            forbidden_elements,
            require_ring,
        )
    )

    def _passes_local(rdmol):
        if max_carbons is not None or min_carbons is not None:
            n_c = sum(1 for a in rdmol.GetAtoms() if a.GetAtomicNum() == 6)
            if max_carbons is not None and n_c > max_carbons:
                return False
            if min_carbons is not None and n_c < min_carbons:
                return False
        if allowed_elements is not None:
            for atom in rdmol.GetAtoms():
                if atom.GetSymbol() not in allowed_elements:
                    return False
        if forbidden_elements is not None:
            for atom in rdmol.GetAtoms():
                if atom.GetSymbol() in forbidden_elements:
                    return False
        if require_ring is not None:
            n_rings = rdMolDescriptors.CalcNumRings(rdmol)
            if require_ring and n_rings == 0:
                return False
            if not require_ring and n_rings > 0:
                return False
        return True

    collected = []
    offset = 0
    # Fetch in pages of 200; when local filters are active we may need to
    # paginate past n API results to collect n passing ones.
    page_size = 200

    while len(collected) < n:
        params["limit"] = page_size
        params["offset"] = offset
        resp = requests.get(base, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()

        for entry in data.get("molecules", []):
            structs = entry.get("molecule_structures")
            if not structs:
                continue
            smi = structs.get("canonical_smiles")
            if not smi:
                continue
            smi = max(smi.split("."), key=len)
            rdmol = Chem.MolFromSmiles(smi)
            if rdmol is None:
                continue
            if has_local_filters and not _passes_local(rdmol):
                continue
            collected.append(Chem.MolToSmiles(rdmol))
            if len(collected) >= n:
                break

        next_page = data.get("page_meta", {}).get("next")
        if not next_page:
            break
        offset += page_size

    return _to_molecules(collected)


def _from_brics(molecules, max_mw=200, min_heavy_atoms=3):
    """
    BRICS retrosynthetic decomposition of a molecule set.

    Fragments are obtained by cutting bonds that are likely retrosynthetically
    accessible (using RDKit's BRICS rules), then stripping the dummy atoms
    left by the cuts and keeping only chemically sensible pieces.

    Parameters
    ----------
    molecules : list
        BuildAMol ``Molecule`` objects or SMILES strings to decompose.
    max_mw : float
        Maximum fragment molecular weight (Da). Default 200 keeps fragments
        small enough to serve as true building blocks.
    min_heavy_atoms : int
        Minimum number of heavy (non-H) atoms; filters out trivial
        single-atom or two-atom fragments.

    Returns
    -------
    list of Molecule
    """
    from rdkit.Chem import BRICS, MolFromSmiles, MolToSmiles, Descriptors, GetMolFrags
    from rdkit.Chem import RWMol

    raw = []
    for mol in molecules:
        smi = mol.to_smiles() if hasattr(mol, "to_smiles") else str(mol)
        rdmol = MolFromSmiles(smi)
        if rdmol is None:
            continue
        raw.extend(BRICS.BRICSDecompose(rdmol))

    # Strip BRICS dummy atoms (atomic number 0) at the graph level.
    # Regex removal of [N*] tokens leaves empty parentheses in SMILES
    # (e.g. Fc1ccc([3*])cc1 → Fc1ccc()cc1) which RDKit rejects — causing
    # fragments like fluorobenzene to be silently dropped.
    seen = set()
    cleaned = []
    for smi in raw:
        rdmol = MolFromSmiles(smi)
        if rdmol is None:
            continue
        rw = RWMol(rdmol)
        dummy_indices = sorted(
            [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0],
            reverse=True,
        )
        for idx in dummy_indices:
            rw.RemoveAtom(idx)
        try:
            from rdkit.Chem import SanitizeMol

            SanitizeMol(rw)
            frags = GetMolFrags(rw.GetMol(), asMols=True)
            if not frags:
                continue
            rdmol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        except Exception:
            continue
        if rdmol.GetNumHeavyAtoms() < min_heavy_atoms:
            continue
        if Descriptors.MolWt(rdmol) > max_mw:
            continue
        canonical = MolToSmiles(rdmol)
        if canonical in seen:
            continue
        seen.add(canonical)
        cleaned.append(canonical)

    return _to_molecules(cleaned)


# ------------------------------------------------------------------
# Shared helper
# ------------------------------------------------------------------


def _to_molecules(smiles_list):
    """Convert a list of SMILES strings to autolabeled BuildAMol Molecules."""
    out = []
    for i, smi in enumerate(smiles_list):
        try:
            mol = Molecule.from_smiles(smi, id=f"F{i:04d}").autolabel()
            out.append(mol)
        except Exception:
            continue
    return out


if __name__ == "__main__":
    import buildamol as bam

    print("Testing ChEMBL source — basic Ro3 (fetches ~20 fragments)...")
    frags = load_fragment_library("chembl", n=20, ro3_pass=True)
    print(f"  Got {len(frags)} fragments")
    for f in frags[:5]:
        print(f"    {f.id}  {f.to_smiles()}")

    print("\nTesting ChEMBL source — small aliphatic C/N/O only...")
    frags2 = load_fragment_library(
        "chembl",
        n=10,
        max_heavy_atoms=10,
        max_carbons=6,
        allowed_elements={"C", "N", "O"},
        require_ring=False,
    )
    print(f"  Got {len(frags2)} fragments")
    for f in frags2:
        print(f"    {f.id}  {f.to_smiles()}")

    print("\nTesting BRICS source (decomposes 3 known drugs)...")
    drugs = [
        bam.Molecule.from_smiles("CC(=O)Oc1ccccc1C(=O)O"),  # aspirin
        bam.Molecule.from_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O"),  # ibuprofen
        bam.Molecule.from_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C"),  # caffeine
    ]
    frags_brics = load_fragment_library("brics", molecules=drugs)
    print(f"  Got {len(frags_brics)} fragments")
    for f in frags_brics:
        print(f"    {f.id}  {f.to_smiles()}")

__all__ = ["load_fragment_library", "make_brics_fragments", "get_chembl_fragments"]
