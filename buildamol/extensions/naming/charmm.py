"""
Mapping from PDB / wwPDB residue names to their CHARMM36(m) equivalents.

Background
----------
CHARMM residue names diverge from PDB names in several systematic ways:

1. Protonation-state disambiguation
   CHARMM uses distinct residue names for different protonation states
   (HSD/HSE/HSP for histidine, ASPP/GLUP for protonated Asp/Glu, etc.)
   whereas PDB stores only HIS/ASP/GLU.  These CANNOT be resolved from
   the name alone — you must inspect the structure.  The mapping here
   returns the most common default form; override as needed.

2. Carbohydrate stereochemistry
   CHARMM encodes anomeric configuration and sugar identity in the name
   (e.g. BGLCNA = beta-D-GlcNAc).  PDB uses a 3-character ligand code
   (NAG).  Additionally CHARMM carbohydrate names can be up to 6 chars,
   which exceeds the PDB 3-char HETATM limit.

3. Lipid names
   Largely match, but a few common lipids differ.

4. Water and ions
   Minor but consistent differences (TIP3 → TIP3, HOH → TIP3, etc.)

5. Nucleic acids (RNA/DNA)
   PDB uses A/G/C/U for RNA; CHARMM uses ADE/GUA/CYT/URA (RNA) and
   DA/DG/DC/DT (DNA in some topology versions).

6. Terminal and capping residues
   PDB deposits ACE/NMA as separate residues; CHARMM treats them as
   patches on the adjacent amino acid, but some psfgen-based workflows
   keep them as standalone residues.

Usage
-----
    from pdb_to_charmm_resnames import pdb_to_charmm, translate_resname

    charmm_name = translate_resname("NAG")          # → "BGLCNA"
    charmm_name = translate_resname("HIS")          # → "HSD"  (default neutral δ)
    charmm_name = translate_resname("HOH")          # → "TIP3"

    # If you want the full entry (with notes):
    entry = pdb_to_charmm.get("NAG")
    print(entry)
    # {"charmm": "BGLCNA", "note": "β-D-GlcNAc; alpha form is AGLCNA"}

Notes on ambiguous cases
------------------------
Several PDB names map to multiple CHARMM names depending on context:

- HIS  → HSD (Nδ protonated, most common in proteins),
         HSE (Nε protonated),
         HSP (doubly protonated / positively charged)
- CYS  → CYS (free thiol) or CYSD (deprotonated) — patch DISU for disulfide
- ASP  → ASP (charged) or ASPP (protonated)
- GLU  → GLU (charged) or GLUP (protonated)
- LYS  → LYS (protonated, +1) or LSN (neutral)
- NAG  → BGLCNA (β, most common in N-linked glycans) or AGLCNA (α)
- MAN  → BMAN (β) or AMAN (α)
- FUC  → AFUC (α, natural form) or BFUC (β)

For these residues translate_resname() returns the most biologically
common default.  Pass explicit_form=True to get a list of all options.
"""

from __future__ import annotations
from typing import Optional, Union, List, Dict, Iterable
import re

from .atom_name_engine import AtomNameGraphEngine
from .residue_name_engine import ResidueNameGraphEngine, ResidueNameLookupEngine

# ---------------------------------------------------------------------------
# Core mapping
# ---------------------------------------------------------------------------
# Each entry: "PDB_CODE": {"charmm": "NAME", "note": "..."}
# "charmm" may be a string (single mapping) or list (ambiguous — first is default).

pdb_to_charmm_mapping: Dict[str, Dict] = {
    # -----------------------------------------------------------------------
    # Standard amino acids — same name in both (included for completeness)
    # -----------------------------------------------------------------------
    "ALA": {"charmm": "ALA", "note": "identical"},
    "ARG": {"charmm": "ARG", "note": "identical"},
    "ASN": {"charmm": "ASN", "note": "identical"},
    "ASP": {
        "charmm": ["ASP", "ASPP"],
        "note": "ASP=charged(default); ASPP=protonated Asp (pKa-shifted)",
    },
    "GLN": {"charmm": "GLN", "note": "identical"},
    "GLU": {
        "charmm": ["GLU", "GLUP"],
        "note": "GLU=charged(default); GLUP=protonated Glu (pKa-shifted)",
    },
    "GLY": {"charmm": "GLY", "note": "identical"},
    "ILE": {"charmm": "ILE", "note": "identical"},
    "LEU": {"charmm": "LEU", "note": "identical"},
    "MET": {"charmm": "MET", "note": "identical"},
    "PHE": {"charmm": "PHE", "note": "identical"},
    "PRO": {"charmm": "PRO", "note": "identical"},
    "SER": {"charmm": "SER", "note": "identical"},
    "THR": {"charmm": "THR", "note": "identical"},
    "TRP": {"charmm": "TRP", "note": "identical"},
    "TYR": {"charmm": "TYR", "note": "identical"},
    "VAL": {"charmm": "VAL", "note": "identical"},
    # Ambiguous protonation state
    "HIS": {
        "charmm": ["HSD", "HSE", "HSP"],
        "note": "HSD=Nδ-H (default); HSE=Nε-H; HSP=doubly protonated (+1)",
    },
    "CYS": {
        "charmm": ["CYS", "CYSD"],
        "note": "CYS=free thiol (default); disulfide handled by DISU patch",
    },
    "LYS": {
        "charmm": ["LYS", "LSN"],
        "note": "LYS=protonated +1 (default); LSN=neutral Lys",
    },
    # N/C-terminal capping (PDB stores as separate residues; CHARMM uses patches)
    "ACE": {
        "charmm": "ACE",
        "note": "Acetyl N-cap; standalone residue in CHARMM (NTER patch also common)",
    },
    "NMA": {
        "charmm": "CT3",
        "note": "N-methylamide C-cap; CHARMM uses CT3 or CTER patch",
    },
    "NH2": {"charmm": "NH2", "note": "C-terminal amide cap"},
    "FOR": {"charmm": "FOR", "note": "Formyl N-cap"},
    # -----------------------------------------------------------------------
    # Non-standard / post-translationally modified amino acids
    # -----------------------------------------------------------------------
    "MSE": {"charmm": "MSE", "note": "Selenomethionine (Se replaces S)"},
    "SEP": {"charmm": "SP1", "note": "Phosphoserine (CHARMM: SP1 or SER+PSER patch)"},
    "TPO": {"charmm": "TP1", "note": "Phosphothreonine"},
    "PTR": {"charmm": "Y1P", "note": "Phosphotyrosine"},
    "HYP": {"charmm": "HYP", "note": "4-hydroxyproline"},
    "MLY": {"charmm": "MLY", "note": "N-methyl-lysine (mono-methylated)"},
    "M3L": {"charmm": "TML", "note": "Trimethyl-lysine"},
    "KCX": {"charmm": "KCXS", "note": "Carboxylysine"},
    "CSO": {"charmm": "CSOH", "note": "S-hydroxycysteine"},
    "CME": {"charmm": "CME", "note": "S-carboxymethyl-cysteine"},
    "OCS": {"charmm": "CYSO", "note": "Cysteine sulfinic acid"},
    "FME": {"charmm": "FME", "note": "N-formylmethionine (N-terminal)"},
    "PCA": {"charmm": "PGLU", "note": "Pyroglutamate (N-terminal cyclized Glu)"},
    # -----------------------------------------------------------------------
    # Water
    # -----------------------------------------------------------------------
    "HOH": {"charmm": "TIP3", "note": "Standard water → TIP3P model"},
    "WAT": {"charmm": "TIP3", "note": "Alternative water name"},
    "H2O": {"charmm": "TIP3", "note": "Alternative water name"},
    "DOD": {"charmm": "TIP3", "note": "Deuterium oxide — usually treated as TIP3"},
    "TIP": {"charmm": "TIP3", "note": "Already CHARMM-like"},
    # -----------------------------------------------------------------------
    # Monatomic ions
    # (PDB uses the element symbol; CHARMM uses the same in most cases)
    # -----------------------------------------------------------------------
    "NA": {"charmm": "SOD", "note": "Sodium; CHARMM uses SOD"},
    "NA+": {"charmm": "SOD", "note": "Sodium (alternative PDB notation)"},
    "CL": {"charmm": "CLA", "note": "Chloride; CHARMM uses CLA"},
    "CL-": {"charmm": "CLA", "note": "Chloride (alternative)"},
    "K": {"charmm": "POT", "note": "Potassium; CHARMM uses POT"},
    "MG": {"charmm": "MG", "note": "Magnesium; same in CHARMM"},
    "CA": {"charmm": "CAL", "note": "Calcium; CHARMM uses CAL"},
    "ZN": {"charmm": "ZN2", "note": "Zinc; CHARMM uses ZN2"},
    "FE": {"charmm": ["FE2", "FE3"], "note": "Iron; FE2=Fe²⁺ (default), FE3=Fe³⁺"},
    "MN": {"charmm": "MN2", "note": "Manganese Mn²⁺"},
    "CO": {"charmm": "CO2", "note": "Cobalt Co²⁺"},
    "CU": {"charmm": ["CU1", "CU2"], "note": "Copper; CU1=Cu⁺, CU2=Cu²⁺"},
    "NI": {"charmm": "NI2", "note": "Nickel Ni²⁺"},
    "CD": {"charmm": "CAD", "note": "Cadmium"},
    "CS": {"charmm": "CES", "note": "Caesium"},
    "RB": {"charmm": "RUB", "note": "Rubidium"},
    "LI": {"charmm": "LIT", "note": "Lithium"},
    "SR": {"charmm": "SR", "note": "Strontium"},
    "BA": {"charmm": "BAR", "note": "Barium"},
    # -----------------------------------------------------------------------
    # Nucleic acids — RNA
    # -----------------------------------------------------------------------
    "A": {"charmm": "ADE", "note": "Adenine ribonucleotide"},
    "G": {"charmm": "GUA", "note": "Guanine ribonucleotide"},
    "C": {"charmm": "CYT", "note": "Cytosine ribonucleotide"},
    "U": {"charmm": "URA", "note": "Uracil ribonucleotide"},
    # PDB uses same single-letter codes with model in mmCIF to distinguish RNA/DNA;
    # for legacy PDB files these are common explicit names:
    "RA": {"charmm": "ADE", "note": "RNA adenosine"},
    "RG": {"charmm": "GUA", "note": "RNA guanosine"},
    "RC": {"charmm": "CYT", "note": "RNA cytidine"},
    "RU": {"charmm": "URA", "note": "RNA uridine"},
    # -----------------------------------------------------------------------
    # Nucleic acids — DNA
    # -----------------------------------------------------------------------
    "DA": {"charmm": "DA", "note": "DNA deoxyadenosine (same in CHARMM)"},
    "DG": {"charmm": "DG", "note": "DNA deoxyguanosine"},
    "DC": {"charmm": "DC", "note": "DNA deoxycytidine"},
    "DT": {"charmm": "DT", "note": "DNA thymidine"},
    "DU": {"charmm": "DU", "note": "DNA deoxyuridine"},
    # -----------------------------------------------------------------------
    # Carbohydrates (CHARMM36 carb FF — top_all36_carb.rtf)
    # PDB 3-char code → CHARMM up-to-6-char stereo-specific name
    # Convention: A=alpha, B=beta; GLC=glucose, MAN=mannose, etc.
    # -----------------------------------------------------------------------
    # GlcNAc / N-acetylglucosamine
    "NAG": {
        "charmm": ["BGLCNA", "AGLCNA"],
        "note": "β-D-GlcNAc (default, N-linked glycans); AGLCNA=α form. "
        "CHARMM-GUI truncates to BGLC in 4-char PSF files.",
    },
    "NDG": {
        "charmm": "BGLCNA",
        "note": "2-acetamido-2-deoxy-β-D-glucopyranose (same as NAG β)",
    },
    # Mannose
    "MAN": {
        "charmm": ["BMAN", "AMAN"],
        "note": "β-D-mannose (default); AMAN=α-D-mannose",
    },
    "BMA": {"charmm": "BMAN", "note": "β-D-mannose (explicit β PDB code)"},
    # Glucose
    "GLC": {
        "charmm": ["BGLC", "AGLC"],
        "note": "β-D-glucose (default); AGLC=α-D-glucose",
    },
    "BGC": {"charmm": "BGLC", "note": "β-D-glucose (explicit)"},
    # Galactose
    "GAL": {
        "charmm": ["BGAL", "AGAL"],
        "note": "β-D-galactose (default); AGAL=α-D-galactose",
    },
    "GLA": {"charmm": "AGAL", "note": "α-D-galactose"},
    # GalNAc / N-acetylgalactosamine
    "GalNAc": {"charmm": ["BGALNA", "AGALNA"], "note": "β-D-GalNAc (default)"},
    "A2G": {"charmm": "BGALNA", "note": "β-D-GalNAc"},
    "NGA": {"charmm": "AGALNA", "note": "α-D-GalNAc"},
    # Fucose
    "FUC": {
        "charmm": ["AFUC", "BFUC"],
        "note": "α-L-fucose (default, biological form); BFUC=β",
    },
    # Sialic acid / Neu5Ac
    "SIA": {
        "charmm": ["ANE5AC", "BNE5AC"],
        "note": "α-2,3- or α-2,6-linked Neu5Ac (default α); BNE5AC=β",
    },
    "SLB": {"charmm": "BNE5AC", "note": "β-Neu5Ac"},
    # Glucuronic acid
    "GCU": {"charmm": ["BGLCA", "AGLCA"], "note": "β-D-glucuronic acid (default)"},
    # Iduronic acid
    "IDR": {"charmm": ["AIDA", "BIDA"], "note": "α-L-iduronic acid (default)"},
    # Xylose
    "XYS": {"charmm": ["BXYL", "AXYL"], "note": "β-D-xylose (default); AXYL=α"},
    # Ribose
    "RIB": {"charmm": ["BRIB", "ARIB"], "note": "β-D-ribose"},
    # Lactose / disaccharides — usually handled as separate monosaccharide residues + patch
    # -----------------------------------------------------------------------
    # Common lipids (CHARMM36 lipid FF — top_all36_lipid.rtf)
    # Most lipid names already match; only exceptions listed.
    # -----------------------------------------------------------------------
    "POPC": {"charmm": "POPC", "note": "identical — palmitoyl-oleoyl-PC"},
    "POPE": {"charmm": "POPE", "note": "identical — palmitoyl-oleoyl-PE"},
    "DPPC": {"charmm": "DPPC", "note": "identical"},
    "DMPC": {"charmm": "DMPC", "note": "identical"},
    "DLPC": {"charmm": "DLPC", "note": "identical"},
    "DSPC": {"charmm": "DSPC", "note": "identical"},
    "DOPC": {"charmm": "DOPC", "note": "identical"},
    "DOPE": {"charmm": "DOPE", "note": "identical"},
    "POPG": {"charmm": "POPG", "note": "identical — palmitoyl-oleoyl-PG"},
    "DPPG": {"charmm": "DPPG", "note": "identical"},
    "DOPS": {"charmm": "DOPS", "note": "identical — dioleoyl-PS"},
    "POPS": {"charmm": "POPS", "note": "identical"},
    "POPE": {"charmm": "POPE", "note": "identical"},
    "CHL1": {
        "charmm": "CHOLE",
        "note": "Cholesterol; PDB uses CHL1, CHARMM uses CHOLE",
    },
    "CLR": {"charmm": "CHOLE", "note": "Cholesterol alternative PDB code"},
    "CHOL": {"charmm": "CHOLE", "note": "Cholesterol (another common abbreviation)"},
    "PALM": {"charmm": "PALM", "note": "Palmitic acid (free fatty acid)"},
    "OLEO": {"charmm": "OLEO", "note": "Oleic acid"},
    "PSM": {"charmm": "PSM", "note": "Palmitoyl sphingomyelin (CHARMM name matches)"},
    "SSM": {"charmm": "SSM", "note": "Stearoyl sphingomyelin"},
    "PNME": {"charmm": "PNME", "note": "POPE-NME (methylated)"},
    "LPC": {"charmm": "LPPC", "note": "Lyso-PC; CHARMM uses LPPC"},
    "PI": {
        "charmm": "POPI",
        "note": "Phosphatidylinositol (approximation — check acyl chains)",
    },
    # -----------------------------------------------------------------------
    # Common cofactors and small molecules with known CHARMM names
    # (CGenFF / CHARMM small molecule library)
    # -----------------------------------------------------------------------
    "HEM": {"charmm": "HEME", "note": "Iron protoporphyrin IX; CHARMM uses HEME"},
    "HEC": {
        "charmm": "HEME",
        "note": "Chloro-heme — treat as HEME, check oxidation state",
    },
    "FAD": {"charmm": "FAD", "note": "Flavin adenine dinucleotide"},
    "FMN": {"charmm": "FMN", "note": "Flavin mononucleotide"},
    "ATP": {"charmm": "ATP", "note": "Adenosine triphosphate"},
    "ADP": {"charmm": "ADP", "note": "Adenosine diphosphate"},
    "AMP": {"charmm": "AMP", "note": "Adenosine monophosphate"},
    "GTP": {"charmm": "GTP", "note": "Guanosine triphosphate"},
    "GDP": {"charmm": "GDP", "note": "Guanosine diphosphate"},
    "NAD": {"charmm": "NAD", "note": "NAD⁺ (oxidised form)"},
    "NAI": {"charmm": "NADN", "note": "NADH (reduced form); PDB uses NAI or NAH"},
    "NDP": {"charmm": "NADP", "note": "NADP⁺"},
    "NDH": {"charmm": "NADPH", "note": "NADPH"},
    "COA": {"charmm": "COA", "note": "Coenzyme A"},
    "HEA": {"charmm": "HEMEA", "note": "Heme A (cytochrome oxidase)"},
    "PLP": {"charmm": "PLP", "note": "Pyridoxal-5-phosphate"},
    "TPP": {"charmm": "TPP", "note": "Thiamine pyrophosphate"},
    "MG": {"charmm": "MG", "note": "Magnesium ion (also listed under ions above)"},
    "POP": {"charmm": "PPI", "note": "Pyrophosphate"},
    "PO4": {
        "charmm": "H2PO4",
        "note": "Inorganic phosphate; protonation state matters",
    },
    "SO4": {"charmm": "SO4", "note": "Sulfate"},
    "GOL": {"charmm": "GLYC", "note": "Glycerol; PDB uses GOL"},
    "EDO": {"charmm": "EOH", "note": "Ethylene glycol (cryo-protectant)"},
    "MPD": {"charmm": "MPD", "note": "2-methyl-2,4-pentanediol"},
    "PEG": {"charmm": "PEG", "note": "Polyethylene glycol fragment"},
    "IMD": {"charmm": "IMID", "note": "Imidazole"},
    "DMS": {"charmm": "DMSO", "note": "Dimethylsulfoxide"},
    "ACT": {"charmm": "ACET", "note": "Acetate"},
    "ACN": {"charmm": "ACNE", "note": "Acetonitrile"},
    "EOH": {"charmm": "ETHA", "note": "Ethanol"},
    "MTH": {"charmm": "METH", "note": "Methanol"},
    "BEN": {"charmm": "BENZ", "note": "Benzene"},
    "PHN": {"charmm": "PHEN", "note": "Phenol"},
    "TRS": {"charmm": "TRIS", "note": "Tris buffer"},
    "PEE": {"charmm": "PETE", "note": "PEG ethyl ester fragment"},
}


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------


def translate_resname(
    pdb_name: str,
    explicit_form: bool = False,
) -> Union[str, List[str], None]:
    """
    Return the CHARMM residue name(s) for a given PDB residue name.

    Parameters
    ----------
    pdb_name : str
        Residue name as it appears in a PDB/mmCIF file (case-insensitive).
    explicit_form : bool
        If True and the mapping is ambiguous, return the full list of
        alternatives.  If False (default), return only the first/default
        CHARMM name.

    Returns
    -------
    str or list of str or None
        The CHARMM name (or list if explicit_form=True), or None if not found.
    """
    entry = pdb_to_charmm_mapping.get(pdb_name.upper())
    if entry is None:
        return None
    charmm = entry["charmm"]
    if explicit_form:
        return charmm if isinstance(charmm, list) else [charmm]
    return charmm[0] if isinstance(charmm, list) else charmm


def get_note(pdb_name: str) -> Optional[str]:
    """Return the explanatory note for a PDB residue name mapping."""
    entry = pdb_to_charmm_mapping.get(pdb_name.upper())
    return entry["note"] if entry else None


def is_ambiguous(pdb_name: str) -> bool:
    """Return True if the PDB name maps to multiple possible CHARMM names."""
    entry = pdb_to_charmm_mapping.get(pdb_name.upper())
    if entry is None:
        return False
    return isinstance(entry["charmm"], list)


def pdb_to_charmm(
    resname_list: List[str],
    warn_ambiguous: bool = True,
    warn_unknown: bool = True,
) -> List[str]:
    """
    Translate a list of PDB residue names to CHARMM names in bulk.

    Residues that are already CHARMM-conforming and not in the mapping
    are returned unchanged (pass-through).

    Parameters
    ----------
    resname_list : list of str
    warn_ambiguous : bool
        Print a warning for residues with multiple CHARMM possibilities.
    warn_unknown : bool
        Print a warning for residue names not found in the mapping
        (returned as-is, assuming they are already CHARMM-conforming or
        are CGenFF ligands).

    Returns
    -------
    list of str
        Translated residue names, same length as input.
    """
    out = []
    for name in resname_list:
        result = translate_resname(name)
        if result is None:
            if warn_unknown:
                print(
                    f"[pdb_to_charmm] Unknown residue: '{name}' — passed through unchanged."
                )
            out.append(name)
        else:
            if warn_ambiguous and is_ambiguous(name):
                alternatives = translate_resname(name, explicit_form=True)
                print(
                    f"[pdb_to_charmm] Ambiguous: '{name}' → defaulting to '{result}'. "
                    f"Alternatives: {alternatives}. Note: {get_note(name)}"
                )
            out.append(result)
    return out


# ---------------------------------------------------------------------------
# Reverse lookup
# ---------------------------------------------------------------------------


def charmm_to_pdb(charmm_name: str) -> List[str]:
    """
    Return all PDB names that map to a given CHARMM residue name.
    Useful for going the other direction (e.g. checking what you have).
    """
    results = []
    for pdb, entry in pdb_to_charmm_mapping.items():
        c = entry["charmm"]
        names = c if isinstance(c, list) else [c]
        if charmm_name.upper() in [n.upper() for n in names]:
            results.append(pdb)
    return results


# ---------------------------------------------------------------------------
# Quick summary / diagnostics
# ---------------------------------------------------------------------------


def print_summary(category: Optional[str] = None) -> None:
    """
    Print a human-readable table of the mapping.

    Parameters
    ----------
    category : str or None
        Filter by keyword in the note, e.g. "carb", "lipid", "ion", "water".
    """
    print(f"{'PDB':<10} {'CHARMM':<20} NOTE")
    print("-" * 72)
    for pdb, entry in sorted(pdb_to_charmm_mapping.items()):
        note = entry["note"]
        if category and category.lower() not in note.lower():
            continue
        charmm = entry["charmm"]
        display = "/".join(charmm) if isinstance(charmm, list) else charmm
        print(f"{pdb:<10} {display:<20} {note}")


def _strip_comment(line: str) -> str:
    return line.split("!", 1)[0].strip()


def _guess_element(
    atom_name: str, atom_type: Optional[str], type_to_element: Dict[str, str]
) -> str:
    if atom_type and atom_type in type_to_element:
        return type_to_element[atom_type]
    token = re.sub(r"[^A-Za-z0-9]", "", atom_name.upper())
    return token[:1] if token else ""


def _parse_rtf_residue_templates(
    filename: str,
    residue_whitelist: Optional[Iterable[str]] = None,
) -> List[dict]:
    wanted = None
    if residue_whitelist is not None:
        wanted = {i.upper() for i in residue_whitelist}

    with open(filename, "r") as handle:
        lines = [_strip_comment(line) for line in handle.readlines()]

    type_to_element: Dict[str, str] = {}
    current_residue = None
    current_atoms = []
    current_bonds = []
    templates = []

    def finalize_residue():
        if current_residue is None or not current_atoms:
            return
        if wanted is not None and current_residue not in wanted:
            return
        templates.append(
            {
                "residue_name": current_residue,
                "atoms": list(current_atoms),
                "bonds": list(current_bonds),
            }
        )

    for raw in lines:
        if not raw:
            continue
        line = raw.split()
        tag = line[0]

        if tag == "MASS" and len(line) >= 5:
            atom_type = line[2]
            element = re.sub(r"[^A-Za-z]", "", line[4]).upper()
            if element:
                type_to_element[atom_type] = element
            continue

        if tag == "RESI":
            finalize_residue()
            current_residue = line[1].upper()
            current_atoms = []
            current_bonds = []
            continue

        if current_residue is None:
            continue

        if tag in ("RESI", "PRES", "END"):
            finalize_residue()
            current_residue = None
            current_atoms = []
            current_bonds = []
            continue

        if tag == "ATOM" and len(line) >= 4:
            atom_name = line[1]
            atom_type = line[2]
            charge = float(line[3])
            current_atoms.append(
                {
                    "name": atom_name,
                    "type": atom_type,
                    "charge": charge,
                    "element": _guess_element(atom_name, atom_type, type_to_element),
                }
            )
            continue

        if tag in ("BOND", "DOUBLE", "TRIPLE") and len(line) > 2:
            pairs = line[1:]
            for i in range(0, len(pairs) - 1, 2):
                current_bonds.append((pairs[i], pairs[i + 1]))
            continue

    finalize_residue()
    return templates


class CHARMMResidueNameLookupEngine(ResidueNameLookupEngine):
    """Dictionary-based residue naming engine backed by `pdb_to_charmm_mapping`."""

    def __init__(self):
        super().__init__(pdb_to_charmm_mapping, value_key="charmm")


class CHARMMResidueNameGraphEngine(ResidueNameGraphEngine):
    """Graph-isomorphism residue naming engine built from CHARMM RTF `RESI` blocks."""

    @classmethod
    def from_file(
        cls,
        filename,
        residue_whitelist: Optional[Iterable[str]] = None,
    ) -> "CHARMMResidueNameGraphEngine":
        """
        Build an engine from one or more CHARMM topology (`.rtf`) files.

        Parameters
        ----------
        filename : str or list of str
            The filename(s) of the CHARMM topology file(s)
        residue_whitelist : iterable of str, optional
            Filter to only these residue names

        Returns
        -------
        CHARMMResidueNameGraphEngine
            The engine with templates loaded from all files
        """
        engine = cls()

        # Handle both single file and list of files
        filenames = filename if isinstance(filename, (list, tuple)) else [filename]

        for fname in filenames:
            engine.load_file(fname, residue_whitelist=residue_whitelist)

        return engine

    def load_file(
        self,
        filename: str,
        residue_whitelist: Optional[Iterable[str]] = None,
    ) -> "CHARMMResidueNameGraphEngine":
        for template in _parse_rtf_residue_templates(
            filename,
            residue_whitelist=residue_whitelist,
        ):
            self.register_template(
                residue_name=template["residue_name"],
                atoms=template["atoms"],
                bonds=template["bonds"],
                source=filename,
            )
        return self


class CHARMMAtomNameGraphEngine(AtomNameGraphEngine):
    """
    Graph-based atom name inference using CHARMM RTF residue templates.

    This parses `RESI` blocks from one or more CHARMM topology files and builds
    residue-internal connectivity templates that can be matched against molecules.
    """

    __tags__ = {
        "MASS",
        "RESI",
        "PRES",
        "GROUP",
        "ATOM",
        "BOND",
        "DOUBLE",
        "TRIPLE",
        "IMPROPER",
        "IMPHI",
        "DIHE",
        "PATCH",
        "END",
    }

    @classmethod
    def from_file(
        cls,
        filename,
        residue_whitelist: Optional[Iterable[str]] = None,
    ) -> "CHARMMAtomNameGraphEngine":
        """
        Build an engine from one or more CHARMM topology (`.rtf`) files.

        Parameters
        ----------
        filename : str or list of str
            The filename(s) of the CHARMM topology file(s)
        residue_whitelist : iterable of str, optional
            Filter to only these residue names

        Returns
        -------
        CHARMMAtomNameGraphEngine
            The engine with templates loaded from all files
        """
        engine = cls()

        # Handle both single file and list of files
        filenames = filename if isinstance(filename, (list, tuple)) else [filename]

        for fname in filenames:
            engine.load_file(fname, residue_whitelist=residue_whitelist)

        return engine

    def load_file(
        self,
        filename: str,
        residue_whitelist: Optional[Iterable[str]] = None,
    ) -> "CHARMMAtomNameGraphEngine":
        for template in _parse_rtf_residue_templates(
            filename,
            residue_whitelist=residue_whitelist,
        ):
            self.register_template(
                residue_name=template["residue_name"],
                atoms=template["atoms"],
                bonds=template["bonds"],
                source=filename,
            )
        return self


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    tests = [
        ("NAG", "BGLCNA"),
        ("HOH", "TIP3"),
        ("HIS", "HSD"),
        ("NA", "SOD"),
        ("CHL1", "CHOLE"),
        ("BMA", "BMAN"),
        ("FUC", "AFUC"),
        ("ALA", "ALA"),
        ("XYZ", None),  # unknown → pass-through
    ]

    print("=== translate_resname tests ===")
    all_pass = True
    for pdb, expected in tests:
        result = translate_resname(pdb)
        status = "✓" if result == expected else "✗"
        if result != expected:
            all_pass = False
        print(f"  {status}  {pdb:6s} → {str(result):<12s}  (expected {expected})")

    print()
    print("=== Ambiguous residues ===")
    for pdb in ["HIS", "NAG", "MAN", "FE", "ASP"]:
        alts = translate_resname(pdb, explicit_form=True)
        print(f"  {pdb:<6s}: {alts}  —  {get_note(pdb)}")

    print()
    print("=== Bulk translation ===")
    chain = ["ALA", "HIS", "NAG", "HOH", "CHL1", "WEIRDNAME"]
    translated = pdb_to_charmm(chain)
    for old, new in zip(chain, translated):
        print(f"  {old:12s} → {new}")

    print()
    print("=== Reverse lookup: what PDB names map to BGLCNA? ===")
    print(" ", charmm_to_pdb("BGLCNA"))

    print()
    print(f"All tests passed: {all_pass}")
