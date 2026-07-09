"""
Functions for loading standard fragment libraries for use with the
``Assembler`` and ``Generator`` classes.

Sources
-------
``"chembl"``
    Rule-of-Three (Ro3) compliant molecules fetched live from the ChEMBL REST
    API.  Requires internet access and ``requests`` (standard in most envs).
    Key kwargs: ``n`` (default 200), ``max_mw`` (300), ``max_hbd`` (3),
    ``max_hba`` (3), ``max_logp`` (3).

``"brics"``
    BRICS retrosynthetic decomposition of a user-supplied molecule set, using
    RDKit's built-in ``rdkit.Chem.BRICS`` module.  No internet needed.
    Key kwargs: ``molecules`` (list of Molecule or SMILES strings),
    ``max_mw`` (200), ``min_heavy_atoms`` (3).

Examples
--------

.. code-block:: python

    from buildamol.extensions.molecular_factories import load_fragment_library

    # ~200 Ro3-compliant fragments from ChEMBL (requires internet)
    fragments = load_fragment_library("chembl", n=200)

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


def get_chembl_fragments(n=200, max_mw=300, max_hbd=3, max_hba=3, max_logp=3):
    """
    Fetch Rule-of-Three compliant fragments from the ChEMBL REST API.

    The Rule of Three (Ro3) is the standard physicochemical filter for
    fragment-based drug discovery: MW ≤ 300 Da, H-bond donors ≤ 3,
    H-bond acceptors ≤ 3, logP ≤ 3.

    Parameters
    ----------
    n : int
        Maximum number of fragments to return.
    max_mw : float
        Maximum molecular weight (Da).
    max_hbd : int
        Maximum H-bond donor count.
    max_hba : int
        Maximum H-bond acceptor count.
    max_logp : float
        Maximum logP.
    """
    return _from_chembl(
        n=n, max_mw=max_mw, max_hbd=max_hbd, max_hba=max_hba, max_logp=max_logp
    )


def load_fragment_library(source="chembl", **kwargs):
    """
    Load a standard fragment library from the given source.
    This function allows for easy switching between different fragment sources.

    Parameters
    ----------
    source : str
        ``"chembl"`` or ``"brics"``.
    **kwargs
        Source-specific keyword arguments (see module docstring).

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


def _from_chembl(n=200, max_mw=300, max_hbd=3, max_hba=3, max_logp=3):
    """
    Fetch Rule-of-Three compliant fragments from the ChEMBL REST API.

    The Rule of Three (Ro3) is the standard physicochemical filter for
    fragment-based drug discovery: MW ≤ 300 Da, H-bond donors ≤ 3,
    H-bond acceptors ≤ 3, logP ≤ 3.

    Parameters
    ----------
    n : int
        Maximum number of fragments to return.
    max_mw : float
        Maximum molecular weight (Da).
    max_hbd : int
        Maximum H-bond donor count.
    max_hba : int
        Maximum H-bond acceptor count.
    max_logp : float
        Maximum logP.

    Returns
    -------
    list of Molecule
    """
    try:
        import requests
    except ImportError:
        raise ImportError("'requests' is required: pip install requests")

    from rdkit import Chem

    base = "https://www.ebi.ac.uk/chembl/api/data/molecule"
    params = {
        "molecule_properties__mw_freebase__lte": max_mw,
        "molecule_properties__hbd__lte": max_hbd,
        "molecule_properties__hba__lte": max_hba,
        "molecule_properties__alogp__lte": max_logp,
        "format": "json",
    }

    collected = []
    offset = 0
    page_size = min(n, 200)

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
            # Remove salt components — keep the largest fragment
            smi = max(smi.split("."), key=len)
            rdmol = Chem.MolFromSmiles(smi)
            if rdmol is None:
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

    print("Testing ChEMBL source (fetches ~20 fragments)...")
    frags = load_fragment_library("chembl", n=20)
    print(f"  Got {len(frags)} fragments")
    for f in frags[:5]:
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
