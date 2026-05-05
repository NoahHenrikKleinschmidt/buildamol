from pathlib import Path

import buildamol as bam
import pytest
import tests.base as base
from buildamol.extensions.naming.atom_name_engine import AtomNameGraphEngine
from buildamol.extensions.naming.residue_name_engine import (
    ResidueNameGraphEngine,
    ResidueNameLookupEngine,
)
from buildamol.resources.pdbe_compounds import PDBECompounds
from buildamol.extensions.naming.charmm import (
    CHARMMAtomNameGraphEngine,
    CHARMMResidueNameLookupEngine,
    CHARMMResidueNameGraphEngine,
    pdb_to_charmm_mapping,
)
from buildamol.extensions.molecular_dynamics.atom_typing.charmm_typer import (
    CHARMMTyper,
    type_with_charmm,
)


class DummyAtom:
    def __init__(self, atom_id, element):
        self.id = atom_id
        self.name = atom_id
        self.element = element
        self.parent = None

    def __repr__(self):
        return f"DummyAtom({self.id})"


class DummyResidue:
    def __init__(self, resname, atoms, bonds):
        self.resname = resname
        self._atoms = atoms
        self._bonds = bonds
        for atom in atoms:
            atom.parent = self

    def get_atoms(self):
        return list(self._atoms)

    def get_bonds(self, residue_internal=True):
        return list(self._bonds)


def _write_simple_rtf(path: Path):
    path.write_text(
        "\n".join(
            [
                "MASS 1 CTYPE 12.011 C",
                "MASS 2 CTYPE2 12.011 C",
                "MASS 3 OTYPE 15.999 O",
                "RESI BGLCNA 0.0",
                "GROUP",
                "ATOM C   CTYPE  0.0",
                "ATOM CT  CTYPE2 0.0",
                "ATOM O1  OTYPE -0.5",
                "BOND C CT  CT O1",
                "END",
            ]
        )
        + "\n"
    )


def _first_residue(mol):
    return next(iter(mol.get_residues()))


def _scramble_residue_atom_names(residue):
    """Deterministically rename atoms to non-template names while preserving element."""
    counts = {}
    for atom in residue.get_atoms():
        element = (getattr(atom, "element", None) or atom.id[0]).upper()
        counts[element] = counts.get(element, 0) + 1
        new_name = f"{element}{90 + counts[element]}"
        atom.id = new_name
        atom.name = new_name


def _heavy_atom_ids(residue):
    return sorted(
        atom.id for atom in residue.get_atoms() if getattr(atom, "element", None) != "H"
    )


def _write_rtf_from_residue(path: Path, residue_name: str, residue):
    """Write a minimal RTF containing one RESI built from a residue graph."""
    atoms = list(residue.get_atoms())
    type_by_atom_name = {}
    mass_lines = []
    seen_type = {}
    idx = 1

    for atom in atoms:
        atom_name = str(atom.id)
        element = (getattr(atom, "element", None) or atom_name[0]).upper()
        atom_type = f"{element}T"
        type_by_atom_name[atom_name] = atom_type
        if atom_type not in seen_type:
            seen_type[atom_type] = True
            mass_lines.append(f"MASS {idx} {atom_type} 12.000 {element}")
            idx += 1

    atom_lines = [
        f"ATOM {str(atom.id):<4} {type_by_atom_name[str(atom.id)]:<4} 0.00"
        for atom in atoms
    ]

    bond_tokens = []
    for bond in residue.get_bonds(residue_internal=True):
        atom1, atom2 = bond
        bond_tokens.extend([str(atom1.id), str(atom2.id)])
    bond_line = "BOND " + " ".join(bond_tokens)

    lines = [
        *mass_lines,
        f"RESI {residue_name} 0.00",
        "GROUP",
        *atom_lines,
        bond_line,
        "END",
    ]
    path.write_text("\n".join(lines) + "\n")


def _make_test_compounds():
    mol = bam.Molecule.from_pdb(base.MANPDB)
    mol.id = "TESTMAN"
    compounds = PDBECompounds({}, id="test-compounds")
    compounds.add(mol, type="SACCHARIDE", names=["mannose"])
    compounds._compounds[mol.id]["name"] = "mannose"
    return mol, compounds


def _assert_atom_engine_matches_reference(engine, reference_mol):
    target = bam.Molecule.from_pdb(base.MANPDB)
    target_residue = _first_residue(target)
    _scramble_residue_atom_names(target_residue)

    engine.rename(target_residue, on_missing="error")

    assert _heavy_atom_ids(target_residue) == _heavy_atom_ids(
        _first_residue(reference_mol)
    )


def _assert_residue_engine_matches_reference(engine, reference_mol):
    target = bam.Molecule.from_pdb(base.MANPDB)
    target_residue = _first_residue(target)
    target_residue.resname = "UNK"

    inferred = engine.infer_residue_name(target_residue)
    assert inferred == _first_residue(reference_mol).resname


def test_base_graph_engines_construct_from_molecule_sources():
    template = bam.Molecule.from_pdb(base.MANPDB)

    atom_engine = AtomNameGraphEngine.from_molecule(template)
    residue_engine = ResidueNameGraphEngine.from_molecule(template)
    multi_atom_engine = AtomNameGraphEngine.from_molecules([template])
    multi_residue_engine = ResidueNameGraphEngine.from_molecules([template])

    assert atom_engine.has_templates(_first_residue(template).resname)
    assert _first_residue(template).resname in residue_engine.residue_names
    _assert_atom_engine_matches_reference(atom_engine, template)
    _assert_residue_engine_matches_reference(residue_engine, template)
    _assert_atom_engine_matches_reference(multi_atom_engine, template)
    _assert_residue_engine_matches_reference(multi_residue_engine, template)


def test_base_graph_engines_construct_from_pdbe_compounds():
    template, compounds = _make_test_compounds()

    atom_engine = AtomNameGraphEngine.from_compounds(compounds)
    residue_engine = ResidueNameGraphEngine.from_compounds(compounds)

    _assert_atom_engine_matches_reference(atom_engine, template)
    _assert_residue_engine_matches_reference(residue_engine, template)


@pytest.mark.parametrize(
    ("factory_name", "suffix", "writer_name"),
    [
        ("from_json", ".json", "to_json"),
        ("from_xml", ".xml", "to_xml"),
        ("from_pickle", ".pkl", "save"),
    ],
)
def test_base_graph_engines_construct_from_pdbe_compound_files(
    tmp_path, factory_name, suffix, writer_name
):
    template, compounds = _make_test_compounds()
    filename = tmp_path / f"compounds{suffix}"
    getattr(compounds, writer_name)(str(filename))

    atom_engine = getattr(AtomNameGraphEngine, factory_name)(str(filename))
    residue_engine = getattr(ResidueNameGraphEngine, factory_name)(str(filename))

    _assert_atom_engine_matches_reference(atom_engine, template)
    _assert_residue_engine_matches_reference(residue_engine, template)


def test_base_graph_engines_construct_from_engine_pickles(tmp_path):
    template = bam.Molecule.from_pdb(base.MANPDB)

    atom_engine = AtomNameGraphEngine.from_molecule(template)
    atom_pickle = tmp_path / "atom_engine.pkl"
    atom_engine.save(str(atom_pickle))

    residue_engine = ResidueNameGraphEngine.from_molecule(template)
    residue_pickle = tmp_path / "residue_engine.pkl"
    residue_engine.save(str(residue_pickle))

    loaded_atom_engine = AtomNameGraphEngine.from_pickle(str(atom_pickle))
    loaded_residue_engine = ResidueNameGraphEngine.from_pickle(str(residue_pickle))

    _assert_atom_engine_matches_reference(loaded_atom_engine, template)
    _assert_residue_engine_matches_reference(loaded_residue_engine, template)


def test_charmm_atom_name_engine_infers_names_from_connectivity(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    engine = CHARMMAtomNameGraphEngine.from_file(str(rtf_file))

    a = DummyAtom("C7", "C")
    b = DummyAtom("C8", "C")
    c = DummyAtom("O9", "O")
    residue = DummyResidue("BGLCNA", [a, b, c], [(a, b), (b, c)])

    inferred = engine.infer_atom_names(residue)

    assert inferred[a] == "C"
    assert inferred[b] == "CT"
    assert inferred[c] == "O1"


def test_charmm_typing_wrapper_renames_before_lookup(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))

    a = DummyAtom("C7", "C")
    b = DummyAtom("C8", "C")
    c = DummyAtom("O9", "O")
    residue = DummyResidue("BGLCNA", [a, b, c], [(a, b), (b, c)])

    with pytest.raises(KeyError):
        typer.get_data(a)

    typed = type_with_charmm(residue, filename=str(rtf_file), typer=typer)

    assert typed[a] == "CTYPE"
    assert typed[b] == "CTYPE2"
    assert typed[c] == "OTYPE"
    assert a.id == "C"
    assert b.id == "CT"
    assert c.id == "O1"


def test_charmm_typer_prefers_existing_type_and_charge_independently(tmp_path, capsys):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))

    atom = DummyAtom("C", "C")
    DummyResidue("BGLCNA", [atom], [])
    atom.type = "EXISTING_TYPE"
    atom.pqr_charge = 0.25

    data = typer.get_data(
        atom,
        keep_existing_type=True,
        keep_existing_charge=False,
    )
    assert data["type"] == "EXISTING_TYPE"
    assert data["charge"] == 0.0

    data = typer.get_data(
        atom,
        keep_existing_type=False,
        keep_existing_charge=True,
    )
    assert data["type"] == "CTYPE"
    assert data["charge"] == 0.25

    out = capsys.readouterr().out
    assert "keeping existing type" in out
    assert "keeping existing charge" in out


def test_charmm_typer_can_keep_existing_when_lookup_is_missing(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))

    atom = DummyAtom("UNKNOWN", "C")
    DummyResidue("BGLCNA", [atom], [])
    atom.type = "EXISTING_TYPE"
    atom.pqr_charge = -0.75

    assert typer.get_type(atom, keep_existing_type=True) == "EXISTING_TYPE"
    assert typer.get_charge(atom, keep_existing_charge=True) == -0.75

    with pytest.raises(KeyError):
        typer.get_type(atom)

    with pytest.raises(KeyError):
        typer.get_charge(atom)


def test_charmm_typer_keep_existing_messages_are_reported_once(tmp_path, capsys):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))

    atom = DummyAtom("C", "C")
    DummyResidue("BGLCNA", [atom], [])
    atom.type = "EXISTING_TYPE"

    typer.get_type(atom, keep_existing_type=True)
    typer.get_type(atom, keep_existing_type=True)

    out = capsys.readouterr().out
    assert out.count("keeping existing type") == 1


def test_charmm_typer_assign_types_does_not_keep_charge_when_only_type_requested(
    tmp_path, capsys
):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))

    atom = DummyAtom("C", "C")
    residue = DummyResidue("BGLCNA", [atom], [])
    atom.type = "EXISTING_TYPE"
    atom.pqr_charge = 1.23

    typer.assign_types(residue, keep_existing=True)
    typer.assign_charges(residue, keep_existing=False)

    out = capsys.readouterr().out
    assert "keeping existing type" in out
    assert "keeping existing charge" not in out
    assert atom.type == "EXISTING_TYPE"
    assert atom.pqr_charge == 0.0


def test_charmm_typer_keep_existing_type_skips_lookup_and_fallback_message(
    tmp_path, capsys
):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))
    typer._pres_dict["HT1"] = {"type": "HC", "charge": 0.1}

    atom = DummyAtom("HT1", "H")
    DummyResidue("ALA", [atom], [])
    atom.type = "HC"

    assert typer.get_type(atom, keep_existing=True) == "HC"

    out = capsys.readouterr().out
    assert "keeping existing type" in out
    assert "using PRES fallback" not in out


def test_charmm_typer_fallback_message_is_contextual_for_charge_pass(tmp_path, capsys):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    typer = CHARMMTyper.from_file(str(rtf_file))
    typer._pres_dict["HT1"] = {"type": "HC", "charge": 0.1}

    atom = DummyAtom("HT1", "H")
    residue = DummyResidue("ALA", [atom], [])
    atom.type = "HC"

    typer.assign_types(residue, keep_existing=True)
    typer.assign_charges(residue, keep_existing=False)

    out = capsys.readouterr().out
    assert "keeping existing type" in out
    assert "using PRES fallback for charge" in out


def test_charmm_residue_lookup_engine_uses_mapping_dict():
    lookup = CHARMMResidueNameLookupEngine()

    assert lookup.translate("NAG") == "BGLCNA"
    assert lookup.translate("HIS") == "HSD"
    assert lookup.translate("HIS", explicit_form=True) == ["HSD", "HSE", "HSP"]


def test_charmm_residue_graph_engine_infers_residue_name_from_connectivity(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)

    engine = CHARMMResidueNameGraphEngine.from_file(str(rtf_file))

    a = DummyAtom("C7", "C")
    b = DummyAtom("C8", "C")
    c = DummyAtom("O9", "O")
    residue = DummyResidue("UNK", [a, b, c], [(a, b), (b, c)])

    inferred_name = engine.infer_residue_name(residue)
    assert inferred_name == "BGLCNA"


def test_atom_graph_engine_partial_match_respects_thresholds(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)
    engine = CHARMMAtomNameGraphEngine.from_file(str(rtf_file))

    # Partial target: only C-CT fragment (2/3 atoms, 1/2 bonds from template).
    a = DummyAtom("C7", "C")
    b = DummyAtom("C8", "C")
    residue = DummyResidue("BGLCNA", [a, b], [(a, b)])

    inferred = engine.infer_atom_names(
        residue,
        min_node_coverage=0.60,
        min_edge_coverage=0.50,
    )
    assert inferred[a] == "C"
    assert inferred[b] == "CT"

    strict = engine.infer_atom_names(
        residue,
        min_node_coverage=0.90,
        min_edge_coverage=0.90,
    )
    assert strict == {}


def test_residue_graph_engine_partial_match_respects_thresholds(tmp_path):
    rtf_file = tmp_path / "tiny.rtf"
    _write_simple_rtf(rtf_file)
    engine = CHARMMResidueNameGraphEngine.from_file(str(rtf_file))

    a = DummyAtom("C7", "C")
    b = DummyAtom("C8", "C")
    residue = DummyResidue("UNK", [a, b], [(a, b)])

    inferred = engine.infer_residue_name(
        residue,
        min_node_coverage=0.60,
        min_edge_coverage=0.50,
    )
    assert inferred == "BGLCNA"

    strict = engine.infer_residue_name(
        residue,
        min_node_coverage=0.90,
        min_edge_coverage=0.90,
    )
    assert strict is None


def test_lookup_engine_rename_molecule_residue():
    mol = bam.Molecule.from_pdb(base.GLCPDB)
    residue = _first_residue(mol)
    residue.resname = "NAG"

    lookup = ResidueNameLookupEngine(pdb_to_charmm_mapping, value_key="charmm")
    lookup.rename(mol)

    assert residue.resname == "BGLCNA"


def test_atom_graph_engine_renames_real_glycan_atoms():
    template = bam.Molecule.from_pdb(base.MANPDB)
    target = bam.Molecule.from_pdb(base.MANPDB)

    template_residue = _first_residue(template)
    target_residue = _first_residue(target)

    original_atom_names = sorted(atom.id for atom in target_residue.get_atoms())
    _scramble_residue_atom_names(target_residue)
    scrambled_names = sorted(atom.id for atom in target_residue.get_atoms())
    assert scrambled_names != original_atom_names

    engine = AtomNameGraphEngine()
    engine.register_molecule(template)
    engine.rename(target_residue)

    assert _heavy_atom_ids(target_residue) == _heavy_atom_ids(template_residue)


def test_residue_graph_engine_renames_real_glycan_residue_name():
    template = bam.Molecule.from_pdb(base.GLCPDB)
    target = bam.Molecule.from_pdb(base.GLCPDB)
    target_residue = _first_residue(target)
    target_residue.resname = "UNK"

    engine = ResidueNameGraphEngine()
    engine.register_molecule(template)
    engine.rename(target)

    assert target_residue.resname == _first_residue(template).resname


def test_charmm_graph_engines_full_rename_on_real_glycan(tmp_path):
    template = bam.Molecule.from_pdb(base.MANPDB)
    template_residue = _first_residue(template)

    rtf_file = tmp_path / "man_template.rtf"
    _write_rtf_from_residue(rtf_file, "BMAN", template_residue)

    target = bam.Molecule.from_pdb(base.MANPDB)
    target_residue = _first_residue(target)
    target_residue.resname = "MAN"
    _scramble_residue_atom_names(target_residue)

    res_engine = CHARMMResidueNameGraphEngine.from_file(str(rtf_file))
    atom_engine = CHARMMAtomNameGraphEngine.from_file(str(rtf_file))

    inferred_resname = res_engine.infer_residue_name(target_residue)
    assert inferred_resname == "BMAN"

    target_residue.resname = inferred_resname
    atom_engine.rename(target_residue)

    assert _heavy_atom_ids(target_residue) == _heavy_atom_ids(template_residue)


def test_atom_graph_engine_can_rename_whole_real_glycan_molecule():
    template = bam.Molecule.from_pdb(base.MANPDB)
    target = bam.Molecule.from_pdb(base.MANPDB)

    template_residue = _first_residue(template)
    target_residue = _first_residue(target)
    _scramble_residue_atom_names(target_residue)

    engine = AtomNameGraphEngine()
    engine.register_molecule(template)
    engine.rename(target)

    assert _heavy_atom_ids(target_residue) == _heavy_atom_ids(template_residue)


def test_atom_graph_engine_on_missing_warn_and_error():
    template = bam.Molecule.from_pdb(base.MANPDB)
    target = bam.Molecule.from_pdb(base.MANPDB)
    target_residue = _first_residue(target)

    _scramble_residue_atom_names(target_residue)

    engine = AtomNameGraphEngine()
    engine.register_molecule(template)

    # Make residue non-matchable by residue-name gating.
    target_residue.resname = "UNK"

    with pytest.warns(UserWarning, match="No atom naming template could be matched"):
        engine.rename(target_residue, on_missing="warn")

    with pytest.raises(KeyError, match="No atom naming template could be matched"):
        engine.rename(target_residue, on_missing="error")


def test_residue_graph_engine_on_missing_warn_and_error():
    template = bam.Molecule.from_pdb(base.GLCPDB)
    target = bam.Molecule.from_pdb(base.MANPDB)

    engine = ResidueNameGraphEngine()
    engine.register_molecule(template)

    target_residue = _first_residue(target)
    with pytest.warns(UserWarning, match="No residue naming template could be matched"):
        engine.rename_residue(target_residue, on_missing="warn")

    with pytest.raises(KeyError, match="No residue naming template could be matched"):
        engine.rename_residue(target_residue, on_missing="error")


def test_residue_lookup_engine_on_missing_warn_and_error():
    mol = bam.Molecule.from_pdb(base.GLCPDB)
    residue = _first_residue(mol)
    residue.resname = "QQQ"

    lookup = ResidueNameLookupEngine(pdb_to_charmm_mapping, value_key="charmm")

    with pytest.warns(UserWarning, match="No residue-name mapping could be found"):
        lookup.rename(mol, on_missing="warn")

    residue.resname = "QQQ"
    with pytest.raises(KeyError, match="No residue-name mapping could be found"):
        lookup.rename(mol, on_missing="error")


def test_charmm_glycan():
    residue_engine = CHARMMResidueNameLookupEngine()
    atom_engine = CHARMMAtomNameGraphEngine.from_file(
        "/Users/noahhk/Downloads/toppar36/top_all36_carb.rtf"
    )
    bam.load_sugars()
    mol = bam.get_compound("GlcNAc")[0]
    old_res_names = [residue.resname for residue in mol.get_residues()]
    old_atom_names = [atom.id for atom in mol.get_atoms()]
    residue_engine.rename(mol, on_missing="error")
    new_res_names = [residue.resname for residue in mol.get_residues()]
    assert new_res_names == ["BGLCNA"] * len(new_res_names)
    atom_engine.rename(mol, on_missing="error")
    new_atom_names = [atom.id for atom in mol.get_atoms()]
    assert "C" in new_atom_names
    assert "CT" in new_atom_names
    assert "C7" not in new_atom_names


def test_charmm_glycan_multiple_res():
    residue_engine = CHARMMResidueNameLookupEngine()
    atom_engine = CHARMMAtomNameGraphEngine.from_file(
        "/Users/noahhk/Downloads/toppar36/top_all36_carb.rtf"
    )
    bam.load_sugars()
    mol = bam.get_compound("GlcNAc")[0] % "14bb" * 3
    old_res_names = [residue.resname for residue in mol.get_residues()]
    old_atom_names = [atom.id for atom in mol.get_atoms()]
    residue_engine.rename(mol, on_missing="error")
    new_res_names = [residue.resname for residue in mol.get_residues()]
    assert new_res_names == ["BGLCNA"] * len(new_res_names)
    atom_engine.rename(mol, on_missing="error")
    new_atom_names = [atom.id for atom in mol.get_atoms()]
    assert new_atom_names.count("C") == 3
    assert new_atom_names.count("CT") == 3
    assert "C7" not in new_atom_names
