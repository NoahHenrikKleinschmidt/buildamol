"""
DiffDock docking backend for BuildAMol.

Uses the DiffDock SE(3)-equivariant diffusion model for blind protein-ligand
docking (no binding-site specification required).

Reference:
    Corso et al. (2023). DiffDock: Diffusion Steps, Twists, and Turns for
    Molecular Docking. ICLR 2023.
"""

import copy
import logging
import os
import sys
import tempfile
from argparse import Namespace
from functools import partial
from pathlib import Path

import numpy as np
import torch
import yaml

import buildamol.core as core

logger = logging.getLogger(__name__)


# Optimised temperature-annealing parameters from default_inference_args.yaml
# Order: [translation, rotation, torsion]
_DEFAULT_TEMP_SAMPLING = [1.170050527854316, 2.06391612594481, 7.044261621607846]
_DEFAULT_TEMP_PSI = [0.727287304570729, 0.9022615585677628, 0.5946212391366862]
_DEFAULT_TEMP_SIGMA_DATA = [0.9299802531572672, 0.7464326999906034, 0.6943254174849822]
_DEFAULT_INITIAL_NOISE_STD = 1.4601642460337794


def dock(
    protein,
    ligand,
    weights,
    return_type="coords",
    samples=10,
    **kwargs,
):
    """
    Dock a ligand to a protein using DiffDock.

    Parameters
    ----------
    protein : Molecule or str
        The protein to dock to. Can be a BuildAMol Molecule or a path to a PDB file.
    ligand : Molecule, str, or rdkit.Chem.Mol
        The ligand to dock. Can be a BuildAMol Molecule, a SMILES string, or an
        RDKit Mol object.
    weights : str or Path
        Path to the DiffDock weights directory. Expected to contain
        ``score_model/`` and optionally ``confidence_model/`` subdirectories,
        each with a ``model_parameters.yml`` and a checkpoint ``.pt`` file.
    return_type : str, optional
        ``'coords'`` (default) – return a numpy coordinate array of shape
        ``(samples, n_atoms, 3)``.
        ``'molecule'`` – return a list of BuildAMol Molecule objects.
    samples : int, optional
        Number of docked poses to generate (default 10).
    **kwargs
        Extra keyword arguments forwarded to :class:`Docker`.

    Returns
    -------
    coords_or_molecules : np.ndarray or list[Molecule]
        Docked poses sorted by confidence score (highest first).
    scores : np.ndarray
        Confidence scores, shape ``(samples,)``.
    """
    docker = Docker(
        protein, weights, return_type=return_type, samples=samples, **kwargs
    )
    return docker(ligand)


class Docker:
    """
    DiffDock protein-ligand docker.

    The protein is preprocessed once during ``__init__`` (ESM-2 language-model
    embeddings + receptor graph construction) and the result is cached so that
    successive calls to :meth:`dock` with different ligands are fast.

    Parameters
    ----------
    protein : Molecule or str or Path
        The protein to dock to.  Accepts a BuildAMol Molecule or a path to a
        PDB file.
    weights : str or Path
        Path to the DiffDock weights directory.  It should contain
        ``score_model/`` and optionally ``confidence_model/``
        subdirectories with ``model_parameters.yml`` and checkpoint files.
        Alternatively, pass the score-model directory directly.
    return_type : str, optional
        ``'coords'`` (default) or ``'molecule'``.
    samples : int, optional
        Number of sampled poses per call (default 10).
    device : str or torch.device, optional
        Compute device (default: CUDA if available, otherwise CPU).
    inference_steps : int, optional
        Total diffusion schedule length (default 20).
    actual_steps : int, optional
        Steps actually executed during reverse diffusion (default 19).
    batch_size : int, optional
        Batch size for the diffusion loop (default 10).
    diffdock_root : str or Path, optional
        Path to the cloned DiffDock repository.
        Defaults to ``/Users/noahhk/GIT/DiffDock``.

    Examples
    --------
    >>> import buildamol as bam
    >>> from buildamol.extensions.docking.backend_diffdock import Docker
    >>> prot = bam.read_pdb("ribosome_fragment.pdb")
    >>> docker = Docker(prot, weights="workdir/v1.1")
    >>> coords, scores = docker("CC(=O)Oc1ccccc1C(=O)O")  # aspirin SMILES
    >>> best_pose_coords = coords[0]          # shape (n_atoms, 3)
    >>> print(f"Top confidence score: {scores[0]:.3f}")
    """

    def __init__(
        self,
        protein,
        weights,
        return_type="coords",
        samples=10,
        device=None,
        inference_steps=20,
        actual_steps=19,
        batch_size=10,
        diffdock_root=None,
    ):
        self.return_type = return_type
        self.samples = samples
        self.inference_steps = inference_steps
        self.actual_steps = actual_steps
        self.batch_size = batch_size
        self._protein_tmpdir = None

        root = Path(diffdock_root) if diffdock_root else Path(weights).parent
        _ensure_diffdock_importable(root)

        if device is None:
            self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self._device = torch.device(device)

        logger.info(f"DiffDock will run on {self._device}")

        self._load_models(Path(weights))
        self._preprocess_protein(protein)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def __call__(self, ligand):
        """Alias for :meth:`dock`."""
        return self.dock(ligand)

    def dock(self, ligand):
        """
        Run DiffDock on one or more ligands against the pre-processed protein.

        Parameters
        ----------
        ligand : Molecule, str, rdkit.Chem.Mol, or list thereof
            The ligand(s) to dock.

        Returns
        -------
        coords_or_molecules : np.ndarray or list[Molecule]
            Shape ``(samples, n_atoms, 3)`` when ``return_type='coords'``.
            List of BuildAMol Molecules when ``return_type='molecule'``.
            Sorted by confidence (highest first).
        scores : np.ndarray
            Confidence scores, shape ``(samples,)``.
            Returns zeros if no confidence model was loaded.
        """
        if isinstance(ligand, (list, tuple)):
            return [self.dock(lig) for lig in ligand]

        # These imports succeed because DiffDock root is already on sys.path
        from utils.sampling import randomize_position, sampling
        from utils.diffusion_utils import get_t_schedule

        mol = self._prepare_mol(ligand)
        complex_graph = self._build_complex_graph(mol)

        # N independent copies for parallel diffusion sampling
        data_list = [copy.deepcopy(complex_graph) for _ in range(self.samples)]

        # Confidence model needs its own copy list (positions get overwritten
        # by sampling() with the score-model output before computing scores)
        confidence_data_list = (
            [copy.deepcopy(complex_graph) for _ in range(self.samples)]
            if self._confidence_model is not None
            else None
        )

        randomize_position(
            data_list,
            self._score_model_args.no_torsion,
            False,
            self._score_model_args.tr_sigma_max,
            initial_noise_std_proportion=_DEFAULT_INITIAL_NOISE_STD,
        )

        tr_schedule = get_t_schedule(
            inference_steps=self.inference_steps, sigma_schedule="expbeta"
        )

        data_list, confidence = sampling(
            data_list=data_list,
            model=self._score_model,
            inference_steps=self.actual_steps,
            tr_schedule=tr_schedule,
            rot_schedule=tr_schedule,
            tor_schedule=tr_schedule,
            device=self._device,
            t_to_sigma=self._t_to_sigma,
            model_args=self._score_model_args,
            confidence_model=self._confidence_model,
            confidence_data_list=confidence_data_list,
            confidence_model_args=self._confidence_model_args,
            batch_size=self.batch_size,
            no_final_step_noise=True,
            temp_sampling=_DEFAULT_TEMP_SAMPLING,
            temp_psi=_DEFAULT_TEMP_PSI,
            temp_sigma_data=_DEFAULT_TEMP_SIGMA_DATA,
        )

        # Re-add protein centroid offset → absolute Cartesian coordinates
        orig_center = complex_graph.original_center.cpu().numpy()
        ligand_pos = np.stack(
            [g["ligand"].pos.cpu().numpy() + orig_center for g in data_list]
        )

        # Sort poses by confidence score (descending)
        if confidence is not None:
            if self._confidence_model_args is not None and isinstance(
                getattr(
                    self._confidence_model_args, "rmsd_classification_cutoff", None
                ),
                list,
            ):
                confidence = confidence[:, 0]
            confidence = confidence.cpu().numpy()
            order = np.argsort(confidence)[::-1]
            confidence = confidence[order]
            ligand_pos = ligand_pos[order]
        else:
            confidence = np.zeros(self.samples)

        if self.return_type == "molecule":
            return self._to_molecules(ligand_pos, mol), confidence
        return ligand_pos, confidence

    # ------------------------------------------------------------------
    # Internal: model loading
    # ------------------------------------------------------------------

    def _load_models(self, weights: Path):
        from utils.utils import get_model
        from utils.diffusion_utils import t_to_sigma as _t_to_sigma_fn

        # Resolve score / confidence directories
        if (weights / "score_model").is_dir():
            score_dir = weights / "score_model"
            conf_dir = (
                weights / "confidence_model"
                if (weights / "confidence_model").is_dir()
                else None
            )
        else:
            score_dir = weights
            conf_dir = None

        if not score_dir.is_dir():
            raise FileNotFoundError(
                f"DiffDock score model directory not found: {score_dir}\n"
                "Pass the path to the weights directory (containing score_model/ "
                "and confidence_model/) or directly to the score model folder."
            )

        # Load score model hyper-parameters and instantiate model
        with open(score_dir / "model_parameters.yml") as f:
            self._score_model_args = Namespace(**yaml.full_load(f))

        self._t_to_sigma = partial(_t_to_sigma_fn, args=self._score_model_args)
        # knn_only_graph matches inference.py: False unless not_knn_only_graph=False
        self._knn_only_graph = not getattr(
            self._score_model_args, "not_knn_only_graph", True
        )

        logger.info("Loading DiffDock score model …")
        score_model = get_model(
            self._score_model_args,
            self._device,
            t_to_sigma=self._t_to_sigma,
            no_parallel=True,
            old=False,
        )
        ckpt_path = score_dir / "best_ema_inference_epoch_model.pt"
        if not ckpt_path.exists():
            raise FileNotFoundError(f"Score model checkpoint not found: {ckpt_path}")
        state_dict = torch.load(ckpt_path, map_location="cpu")
        score_model.load_state_dict(state_dict, strict=True)
        score_model.to(self._device).eval()
        self._score_model = score_model
        logger.info("Score model loaded.")

        # Load confidence model if present
        if conf_dir is not None:
            with open(conf_dir / "model_parameters.yml") as f:
                self._confidence_model_args = Namespace(**yaml.full_load(f))

            logger.info("Loading DiffDock confidence model …")
            conf_model = get_model(
                self._confidence_model_args,
                self._device,
                t_to_sigma=self._t_to_sigma,
                no_parallel=True,
                confidence_mode=True,
                old=True,
            )
            # Try common checkpoint names
            for ckpt_name in ("best_model_epoch75.pt", "best_model.pt"):
                conf_ckpt = conf_dir / ckpt_name
                if conf_ckpt.exists():
                    break
            else:
                raise FileNotFoundError(
                    f"No confidence model checkpoint found in {conf_dir}"
                )
            state_dict = torch.load(conf_ckpt, map_location="cpu")
            conf_model.load_state_dict(state_dict, strict=True)
            conf_model.to(self._device).eval()
            self._confidence_model = conf_model
            logger.info("Confidence model loaded.")
        else:
            self._confidence_model = None
            self._confidence_model_args = None

    # ------------------------------------------------------------------
    # Internal: protein preprocessing (runs once)
    # ------------------------------------------------------------------

    def _preprocess_protein(self, protein):
        """
        Build and cache the receptor graph.

        Computes ESM-2 embeddings and constructs the HeteroData receptor
        subgraph so that ligand docking calls only need to add the ligand graph.
        """
        from torch_geometric.data import HeteroData
        from utils.inference_utils import (
            get_sequences_from_pdbfile,
            compute_ESM_embeddings,
        )
        from datasets.process_mols import moad_extract_receptor_structure
        from esm import pretrained

        # Materialise protein as a PDB file if it's a Molecule object
        if isinstance(protein, (str, Path)):
            protein_file = str(protein)
        else:
            self._protein_tmpdir = tempfile.TemporaryDirectory()
            protein_file = os.path.join(self._protein_tmpdir.name, "receptor.pdb")
            protein.to_pdb(protein_file)

        # ESM-2 language model embeddings (expensive — done once)
        logger.info("Computing ESM-2 receptor embeddings (this runs once) …")
        sequence = get_sequences_from_pdbfile(protein_file)
        esm_model, alphabet = pretrained.load_model_and_alphabet("esm2_t33_650M_UR50D")
        esm_model.eval()
        if torch.cuda.is_available():
            esm_model = esm_model.cuda()

        chains = sequence.split(":")
        labels = [f"receptor_chain_{j}" for j in range(len(chains))]
        all_embeddings = compute_ESM_embeddings(esm_model, alphabet, labels, chains)
        lm_embeddings = [all_embeddings[lbl] for lbl in labels]
        del esm_model  # release GPU memory

        # Build receptor HeteroData graph (expensive — done once)
        logger.info("Building receptor graph …")
        receptor_graph = HeteroData()
        moad_extract_receptor_structure(
            path=protein_file,
            complex_graph=receptor_graph,
            neighbor_cutoff=self._score_model_args.receptor_radius,
            max_neighbors=getattr(
                self._score_model_args, "c_alpha_max_neighbors", None
            ),
            lm_embeddings=lm_embeddings,
            knn_only_graph=self._knn_only_graph,
            all_atoms=self._score_model_args.all_atoms,
            atom_cutoff=getattr(self._score_model_args, "atom_radius", 5),
            atom_max_neighbors=getattr(
                self._score_model_args, "atom_max_neighbors", None
            ),
        )

        # Shift receptor to its own centroid (ligand will be centred separately)
        protein_center = torch.mean(receptor_graph["receptor"].pos, dim=0, keepdim=True)
        receptor_graph["receptor"].pos -= protein_center
        if self._score_model_args.all_atoms:
            receptor_graph["atom"].pos -= protein_center
        receptor_graph.original_center = protein_center

        self._receptor_template = receptor_graph
        logger.info("Protein preprocessing complete.")

    # ------------------------------------------------------------------
    # Internal: per-call graph construction
    # ------------------------------------------------------------------

    def _build_complex_graph(self, mol, name="ligand"):
        """
        Assemble a full HeteroData complex graph from the cached receptor
        template and a freshly featurised ligand.
        """
        from datasets.process_mols import get_lig_graph_with_matching

        complex_graph = copy.deepcopy(self._receptor_template)
        complex_graph["name"] = name

        get_lig_graph_with_matching(
            mol,
            complex_graph,
            popsize=None,
            maxiter=None,
            matching=False,
            keep_original=False,
            num_conformers=1,
            remove_hs=getattr(self._score_model_args, "remove_hs", False),
        )

        # Centre ligand at its own centroid (receptor already centred)
        lig_center = torch.mean(complex_graph["ligand"].pos, dim=0, keepdim=True)
        complex_graph["ligand"].pos -= lig_center

        complex_graph.mol = mol
        complex_graph["success"] = True
        return complex_graph

    @staticmethod
    def _prepare_mol(ligand):
        """
        Convert the ligand input to an RDKit Mol with a 3D conformer suitable
        for DiffDock featurisation.
        """
        from rdkit.Chem import MolFromSmiles, AddHs
        from datasets.process_mols import generate_conformer, read_molecule

        if isinstance(ligand, str):
            # Try as SMILES first, then as a file path
            mol = MolFromSmiles(ligand)
            if mol is None:
                mol = read_molecule(ligand, remove_hs=False, sanitize=True)
                if mol is None:
                    raise ValueError(
                        f"Could not parse ligand — not a valid SMILES or "
                        f"readable molecule file: {ligand!r}"
                    )
                mol.RemoveAllConformers()
            mol = AddHs(mol)
            generate_conformer(mol)

        elif hasattr(ligand, "to_smiles"):
            # BuildAMol Molecule — go via SMILES for a clean 2D→3D pipeline
            smiles = ligand.to_smiles()
            mol = MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(
                    f"Could not parse SMILES obtained from ligand: {smiles!r}"
                )
            mol = AddHs(mol)
            generate_conformer(mol)

        elif hasattr(ligand, "GetAtoms"):
            # RDKit Mol — regenerate conformer from scratch (DiffDock uses it
            # as a starting topology, not a starting pose)
            mol = AddHs(copy.deepcopy(ligand))
            mol.RemoveAllConformers()
            generate_conformer(mol)

        else:
            raise TypeError(
                f"Unsupported ligand type {type(ligand)}. "
                "Pass a SMILES string, BuildAMol Molecule, or RDKit Mol."
            )

        return mol

    def _to_molecules(self, ligand_pos, mol):
        """Convert docked coordinate arrays to BuildAMol Molecule objects."""
        from rdkit.Chem import RemoveAllHs

        remove_hs = getattr(self._score_model_args, "remove_hs", False)
        out = []
        for pos in ligand_pos:
            rdmol = copy.deepcopy(mol)
            if remove_hs:
                rdmol = RemoveAllHs(rdmol)
            conf = rdmol.GetConformer()
            for i, (x, y, z) in enumerate(pos):
                conf.SetAtomPosition(i, (float(x), float(y), float(z)))
            out.append(core.Molecule.from_rdkit(rdmol))
        return out

    def __del__(self):
        if self._protein_tmpdir is not None:
            try:
                self._protein_tmpdir.cleanup()
            except Exception:
                pass


# ------------------------------------------------------------------
# Module helpers
# ------------------------------------------------------------------


def _ensure_diffdock_importable(root: Path):
    """Add the DiffDock repository root to sys.path if not already present."""
    root_str = str(root.resolve())
    if not root.is_dir():
        raise FileNotFoundError(
            f"DiffDock repository not found at {root}. "
            "Clone the repo or pass diffdock_root= to Docker."
        )
    if root_str not in sys.path:
        sys.path.insert(0, root_str)


__dock__ = dock
__all__ = ["dock", "Docker"]
