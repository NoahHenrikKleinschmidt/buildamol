"""
Write and read pdbqt files.
"""

import os
import subprocess
import tempfile

import buildamol.utils.auxiliary as aux
import buildamol.utils.pdb as _pdb

has_meeko = aux.has_package("meeko")


def _encode_pdbqt_receptor(molecule: "Molecule", box_center=None, box_size=None) -> str:
    """
    Write a rigid receptor (multi-chain / multi-fragment) to PDBQT via
    mk_prepare_receptor (installed alongside meeko). This is the correct tool
    for proteins and nucleic acid receptors — MoleculePreparation is
    ligand-only and rejects multi-fragment input.
    """
    pdb_tmp = tempfile.mktemp(suffix=".pdb")
    pdbqt_tmp = pdb_tmp.replace(".pdb", ".pdbqt")
    cmd = "mk_prepare_receptor"
    try:
        bare_cmd_is_available = subprocess.run(
            [cmd, "--help"],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        cmd += ".py"
        try:
            with_py_cmd_is_available = subprocess.run(
                [cmd, "--help"],
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            raise RuntimeError(
                "mk_prepare_receptor is not available. Please ensure that meeko is installed and that mk_prepare_receptor is in your PATH."
            )
    try:
        pdb_content = _pdb.encode_pdb(molecule, reindex=True)
        with open(pdb_tmp, "w") as _f:
            _f.write(pdb_content)
        cmd_args = [cmd, "--read_pdb", pdb_tmp, "--write_pdbqt", pdbqt_tmp]
        if box_center is not None:
            cmd_args.extend(
                [
                    "--box_center",
                    str(box_center[0]),
                    str(box_center[1]),
                    str(box_center[2]),
                ]
            )
        if box_size is not None:
            cmd_args.extend(
                [
                    "--box_size",
                    str(box_size[0]),
                    str(box_size[1]),
                    str(box_size[2]),
                ]
            )
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"mk_prepare_receptor failed (exit {result.returncode}):\n"
                f"{result.stderr.strip()}"
            )
        with open(pdbqt_tmp) as f:
            return f.read()
    finally:
        for p in (pdb_tmp, pdbqt_tmp):
            if os.path.exists(p):
                os.unlink(p)


def encode_pdbqt(molecule: "Molecule"):
    """
    Encode a molecule to a pdbqt string.

    Single-fragment molecules (ligands) are prepared with meeko's
    MoleculePreparation. Multi-fragment molecules (protein receptors,
    nucleic acids with multiple chains) are prepared with
    mk_prepare_receptor, which meeko installs as a CLI command.

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    molecule : Molecule
        The molecule to encode

    Returns
    -------
    str
        The pdbqt string
    """
    if not has_meeko:
        raise ImportError("PDBQT encoding requires Meeko")

    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy

    # Check connectivity via the graph — O(1) vs the O(N) to_rdkit() path.
    # A multi-fragment molecule (receptor) has more connected components than 1.
    import networkx as nx
    if nx.number_connected_components(molecule._AtomGraph) > 1:
        return _encode_pdbqt_receptor(molecule)

    rdmol = molecule.to_rdkit()
    molprep = MoleculePreparation()
    pdbqt = molprep.prepare(rdmol)
    if len(pdbqt) == 0:
        raise ValueError("PDBQT encoding failed! Check the input molecule...")
    pdbqt = PDBQTWriterLegacy.write_string(pdbqt[0])
    return pdbqt[0]


def write_pdbqt(molecule: "Molecule", filename: str):
    """
    Write a molecule to a pdbqt file

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    molecule : Molecule
        The molecule to write
    filename : str
        The filename of the pdbqt file
    """
    pdbqt = encode_pdbqt(molecule)
    with open(filename, "w") as f:
        f.write(pdbqt)


def read_pdbqt(filename: str):
    """
    Read a pdbqt file into an array of atoms

    Note
    ----
    This function requires `meeko` to be installed.

    Parameters
    ----------
    filename : str
        The filename of the pdbqt file

    Returns
    -------
    list
        The list of atoms. Each entry is a tuple of form:
        ('idx', 'serial', 'name/element', 'resid', 'resname', 'chain', 'xyz/coord', 'partial_charges', 'atom_type')
    """
    if not has_meeko:
        raise ImportError("PDBQT reading requires Meeko")

    from meeko import PDBQTMolecule

    pdbqt_string = open(filename, "r").read()
    atoms = PDBQTMolecule(pdbqt_string=pdbqt_string).atoms()
    return atoms
