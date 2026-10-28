"""
OpenDock docking backend for BuildAMol.

OpenDock is a PyTorch-based protein-ligand docking framework with multiple
scoring functions (Vina, OnionNet-SFCT, RTMScore) and sampling strategies
(Monte Carlo, gradient-based optimisation, etc.).

The central conversion functions here bridge BuildAMol Molecule objects to
OpenDock's ``LigandConformation`` / ``ReceptorConformation`` classes, both
of which require PDBQT files as input.

Reference:
    https://github.com/guyuehuo/opendock-open
    https://opendock-readthedocs.readthedocs.io
"""

import os
import tempfile
from pathlib import Path

import torch

import buildamol.core as core
import buildamol.utils.auxiliary as aux
from buildamol.utils.pdbqt import encode_pdbqt, _encode_pdbqt_receptor


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _require_opendock():
    if not aux.has_package("opendock"):
        raise ImportError(
            "OpenDock is required. Install it with:\n"
            "    pip install opendock\n"
            "or from source: https://github.com/guyuehuo/opendock-open"
        )


def _require_meeko():
    if not aux.has_package("meeko"):
        raise ImportError(
            "meeko is required for PDBQT preparation. Install with:\n"
            "    pip install meeko"
        )


def _write_pdbqt(content: str, directory: str, filename: str) -> str:
    """Write a PDBQT string to *directory/filename* and return the full path."""
    path = os.path.join(directory, filename)
    with open(path, "w") as f:
        f.write(content)
    return path


# ---------------------------------------------------------------------------
# Core conversion: ligand
# ---------------------------------------------------------------------------


def to_ligand_conformation(ligand, outdir=None):
    """
    Convert a ligand to an OpenDock ``LigandConformation`` object.

    This is the primary entry point for making a BuildAMol molecule usable
    with OpenDock.  The ligand is written to a PDBQT file (with the
    ROOT/ENDROOT/BRANCH torsion tree that OpenDock's parser expects) and then
    parsed into a ``LigandConformation``.

    Parameters
    ----------
    ligand : Molecule, str, or rdkit.Chem.Mol
        The ligand to prepare.  Accepts a BuildAMol Molecule, a SMILES string,
        or an RDKit Mol.  For SMILES, a 3-D conformer is generated via
        RDKit's ETKDG before PDBQT preparation.
    outdir : str or Path, optional
        Directory where the temporary ``ligand.pdbqt`` file is written.
        If *None*, a temporary directory is created automatically.  The file
        is **not** deleted automatically — use the returned
        :class:`OpenDockSystem` (or call :func:`to_opendock`) to manage
        lifetime.  When calling this function in isolation, pass an explicit
        *outdir* and manage the directory yourself.

    Returns
    -------
    opendock.core.conformation.LigandConformation
        Ready to pass to a scoring function or sampler.
    str
        Absolute path of the PDBQT file that was written (useful for
        inspecting or reusing it).

    Examples
    --------
    >>> from buildamol.extensions.docking.backend_opendock import to_ligand_conformation
    >>> lig_cnfr, pdbqt_path = to_ligand_conformation("CC(=O)Oc1ccccc1C(=O)O")
    """
    _require_opendock()
    _require_meeko()

    from opendock.core.conformation import LigandConformation

    # --- resolve ligand to an RDKit mol so meeko can prepare it -------------
    rdmol = _ligand_to_rdmol(ligand)

    # --- meeko: prepare + write PDBQT ----------------------------------------
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    molprep = MoleculePreparation()
    setup_list = molprep.prepare(rdmol)
    if not setup_list:
        raise ValueError(
            "meeko MoleculePreparation returned nothing for this ligand. "
            "Check that the molecule has 3-D coordinates and is properly sanitised."
        )
    pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(setup_list[0])

    if outdir is None:
        outdir = tempfile.mkdtemp()
    outdir = str(outdir)
    pdbqt_path = _write_pdbqt(pdbqt_str, outdir, "ligand.pdbqt")

    cnfr = LigandConformation(pdbqt_path)
    return cnfr, pdbqt_path


# ---------------------------------------------------------------------------
# Core conversion: receptor
# ---------------------------------------------------------------------------


def to_receptor_conformation(protein, box_center, ligand_conformation=None, outdir=None):
    """
    Convert a protein to an OpenDock ``ReceptorConformation`` object.

    Parameters
    ----------
    protein : Molecule or str or Path
        The receptor to prepare.  Accepts a BuildAMol Molecule or a path to a
        PDB or PDBQT file.  PDB files are converted to PDBQT via
        ``mk_prepare_receptor`` (from meeko); PDBQT files are used directly.
    box_center : array-like of float, length 3
        Centre of the docking box in (x, y, z) Å.
    ligand_conformation : LigandConformation, optional
        A pre-built ``LigandConformation``; its heavy-atom coordinates are
        passed to the receptor so it can pre-select the binding pocket.  Not
        strictly required but improves performance.
    outdir : str or Path, optional
        Directory for the receptor PDBQT file.  Defaults to a new temp dir.

    Returns
    -------
    opendock.core.conformation.ReceptorConformation
    str
        Absolute path of the receptor PDBQT file.
    """
    _require_opendock()

    from opendock.core.conformation import ReceptorConformation

    xyz_center = torch.tensor(list(box_center), dtype=torch.float32).reshape(1, 3)

    # --- resolve protein to a PDBQT file -------------------------------------
    if isinstance(protein, (str, Path)):
        fpath = str(protein)
        if fpath.endswith(".pdbqt"):
            rec_pdbqt_path = fpath
            _managed_dir = None
        else:
            # Assume PDB — convert via meeko
            if outdir is None:
                outdir = tempfile.mkdtemp()
            _managed_dir = str(outdir)
            mol = core.Molecule.from_pdb(fpath)
            pdbqt_str = _encode_pdbqt_receptor(mol)
            rec_pdbqt_path = _write_pdbqt(pdbqt_str, _managed_dir, "receptor.pdbqt")
    else:
        # BuildAMol Molecule
        _require_meeko()
        if outdir is None:
            outdir = tempfile.mkdtemp()
        _managed_dir = str(outdir)
        pdbqt_str = _encode_pdbqt_receptor(protein)
        rec_pdbqt_path = _write_pdbqt(pdbqt_str, _managed_dir, "receptor.pdbqt")

    init_lig_xyz = (
        ligand_conformation.init_lig_heavy_atoms_xyz
        if ligand_conformation is not None
        else None
    )

    rec_cnfr = ReceptorConformation(
        rec_pdbqt_path,
        xyz_center,
        init_lig_heavy_atoms_xyz=init_lig_xyz,
    )
    return rec_cnfr, rec_pdbqt_path


# ---------------------------------------------------------------------------
# Combined conversion
# ---------------------------------------------------------------------------


def to_opendock(protein, ligand, box_center, box_size=None):
    """
    Convert a protein–ligand pair to an :class:`OpenDockSystem`.

    This is the primary "convert to OpenDock" function.  It handles PDBQT
    preparation for both the ligand and the receptor, returns a small container
    object that holds the two ``*Conformation`` objects together with the
    temporary directory, and supports use as a context manager for automatic
    cleanup.

    Parameters
    ----------
    protein : Molecule or str or Path
        The receptor.
    ligand : Molecule, str, or rdkit.Chem.Mol
        The ligand.
    box_center : array-like of float, length 3
        Docking box centre in (x, y, z) Å.
    box_size : array-like of float, length 3, optional
        Docking box dimensions (x, y, z) in Å.  Not required for building
        the conformation objects themselves, but stored on the returned system
        for convenience when constructing samplers.

    Returns
    -------
    OpenDockSystem
        Container with ``.ligand``, ``.receptor``, ``.box_center``,
        ``.box_size``, and ``.pdbqt_dir`` attributes, plus context-manager
        support for cleanup.

    Examples
    --------
    >>> import buildamol as bam
    >>> from buildamol.extensions.docking.backend_opendock import to_opendock
    >>> prot = bam.read_pdb("receptor.pdb")
    >>> with to_opendock(prot, "CC(=O)Oc1ccccc1C(=O)O", box_center=(10, 5, 3)) as sys:
    ...     from opendock.scorer.vina import VinaSF
    ...     sf = VinaSF(sys.receptor, sys.ligand)
    ...     print(sf.scoring())
    """
    tmpdir = tempfile.mkdtemp()

    lig_cnfr, lig_path = to_ligand_conformation(ligand, outdir=tmpdir)
    rec_cnfr, rec_path = to_receptor_conformation(
        protein, box_center, ligand_conformation=lig_cnfr, outdir=tmpdir
    )

    return OpenDockSystem(
        ligand=lig_cnfr,
        receptor=rec_cnfr,
        box_center=box_center,
        box_size=box_size,
        pdbqt_dir=tmpdir,
        ligand_pdbqt=lig_path,
        receptor_pdbqt=rec_path,
    )


# ---------------------------------------------------------------------------
# System container
# ---------------------------------------------------------------------------


class OpenDockSystem:
    """
    Container for an OpenDock docking setup.

    Holds the ``LigandConformation``, ``ReceptorConformation``, box
    parameters, and the directory containing the associated PDBQT files.
    Supports use as a context manager so the temporary PDBQT directory is
    cleaned up automatically on exit.

    Attributes
    ----------
    ligand : LigandConformation
    receptor : ReceptorConformation
    box_center : tuple of float
    box_size : tuple of float or None
    pdbqt_dir : str
        Directory containing ``ligand.pdbqt`` and ``receptor.pdbqt``.
    ligand_pdbqt : str
        Full path to the ligand PDBQT file.
    receptor_pdbqt : str
        Full path to the receptor PDBQT file.
    """

    def __init__(
        self,
        ligand,
        receptor,
        box_center,
        box_size,
        pdbqt_dir,
        ligand_pdbqt,
        receptor_pdbqt,
    ):
        self.ligand = ligand
        self.receptor = receptor
        self.box_center = box_center
        self.box_size = box_size
        self.pdbqt_dir = pdbqt_dir
        self.ligand_pdbqt = ligand_pdbqt
        self.receptor_pdbqt = receptor_pdbqt

    def cleanup(self):
        """Remove the temporary PDBQT directory and its contents."""
        import shutil
        if self.pdbqt_dir and os.path.isdir(self.pdbqt_dir):
            shutil.rmtree(self.pdbqt_dir, ignore_errors=True)
            self.pdbqt_dir = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.cleanup()

    def __repr__(self):
        return (
            f"OpenDockSystem("
            f"box_center={self.box_center}, "
            f"box_size={self.box_size}, "
            f"pdbqt_dir={self.pdbqt_dir!r})"
        )


# ---------------------------------------------------------------------------
# Top-level dock
# ---------------------------------------------------------------------------


def dock(
    protein,
    ligand,
    box_center,
    box_size,
    steps=500,
    n_out=10,
    random_start=True,
    rmsd_cutoff=2.0,
    out=None,
):
    """
    Dock a ligand to a protein using OpenDock's Monte Carlo / Vina pipeline.

    Parameters
    ----------
    protein : Molecule or str or Path
        The receptor.
    ligand : Molecule, str, or rdkit.Chem.Mol
        The ligand.
    box_center : array-like of float, length 3
        Docking box centre (x, y, z) in Å.
    box_size : array-like of float, length 3
        Docking box dimensions (x, y, z) in Å.
    steps : int, optional
        Number of Monte Carlo steps (default 500).
    n_out : int, optional
        Maximum number of output poses to return (default 10).
    random_start : bool, optional
        Randomise the ligand starting position before sampling (default True).
    rmsd_cutoff : float, optional
        RMSD threshold for pose clustering (default 2.0 Å).
    out : str or Path, optional
        If given, write the docked trajectory to this PDB file.

    Returns
    -------
    list of tuple
        Each element is ``(score, conformation_tensor)`` for the top poses,
        sorted best score first.
    """
    _require_opendock()

    from opendock.scorer.vina import VinaSF
    from opendock.sampler.monte_carlo import MonteCarloSampler
    from opendock.core.clustering import BaseCluster

    with to_opendock(protein, ligand, box_center, box_size) as system:
        sf = VinaSF(system.receptor, system.ligand)

        sampler = MonteCarloSampler(
            system.ligand,
            system.receptor,
            sf,
            box_center=list(box_center),
            box_size=list(box_size),
            random_start=random_start,
            minimizer=None,
        )
        sampler.sampling(steps)

        if out is not None:
            from opendock.core import io
            io.write_ligand_traj(
                system.ligand.cnfrs_history,
                system.receptor,
                system.ligand,
                str(out),
                scores=system.ligand.scores_history,
            )

        cluster = BaseCluster(
            system.ligand.cnfrs_history,
            None,
            system.ligand.scores_history,
            system.ligand,
            rmsd_cutoff,
        )
        scores, cnfrs, _ = cluster.clustering()

        return list(zip(scores[:n_out], cnfrs[:n_out]))


# ---------------------------------------------------------------------------
# Internal: ligand format resolution
# ---------------------------------------------------------------------------


def _ligand_to_rdmol(ligand):
    """
    Convert various ligand representations to an RDKit Mol with 3-D coords.
    meeko requires an embedded RDKit Mol for PDBQT preparation.
    """
    from rdkit.Chem import MolFromSmiles, AddHs, AllChem

    if isinstance(ligand, str):
        # SMILES
        mol = MolFromSmiles(ligand)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES: {ligand!r}")
        mol = AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result != 0:
            raise ValueError(
                f"3-D conformer generation failed for SMILES: {ligand!r}"
            )
        AllChem.MMFFOptimizeMolecule(mol)
        return mol

    if hasattr(ligand, "to_smiles"):
        # BuildAMol Molecule — go through SMILES for a clean embedding
        return _ligand_to_rdmol(ligand.to_smiles())

    if hasattr(ligand, "GetAtoms"):
        # Already an RDKit Mol — ensure it has 3-D coords
        import copy
        from rdkit.Chem import AddHs
        mol = AddHs(copy.deepcopy(ligand))
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol)
        return mol

    raise TypeError(
        f"Unsupported ligand type {type(ligand)}. "
        "Pass a SMILES string, BuildAMol Molecule, or RDKit Mol."
    )


__dock__ = dock
__all__ = [
    "to_ligand_conformation",
    "to_receptor_conformation",
    "to_opendock",
    "OpenDockSystem",
    "dock",
]
