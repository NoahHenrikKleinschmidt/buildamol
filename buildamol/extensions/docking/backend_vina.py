"""
Docking module for BuildAMol using the AutoDock Vina Python bindings directly (in-process).

This backend is designed for high-throughput use, e.g. in generative pipelines with
optimization loops where many candidates must be evaluated against the same pocket:

- The receptor is prepared once per `Docker`. Affinity maps for all common Vina (XS)
  atom types are computed at setup and re-used for all subsequent ligands; they are only
  re-computed (~1 s) when a ligand introduces a new XS atom type (e.g. macrocycle pseudo atoms).
- Besides a full global search (``mode="dock"``), ligands that were already built
  inside the pocket can be locally refined (``mode="optimize"``) or just scored
  (``mode="score"``), which takes milliseconds instead of seconds to minutes.
- Batches of ligands are docked in parallel, one ligand per worker process, which
  scales much better than giving Vina more threads for a single ligand.

Requires the `vina` and `meeko` packages.
"""

import os
import re
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np

import buildamol.core as core
import buildamol.utils.auxiliary as aux

MODES = ("dock", "optimize", "score")

# Small molecules that together contain all common Vina (XS) atom types.
# The affinity maps are computed for all of these at setup so that typical
# ligands can be evaluated without re-computing maps. Maps are only re-computed
# when a ligand introduces an XS type that is not covered yet.
_REFERENCE_SMILES = (
    "CCO",
    "CC(=O)N(C)C",
    "CC(=O)NC",
    "CN",
    "C[NH3+]",
    "c1ccncc1",
    "CSC",
    "CS",
    "CP(=O)(O)O",
    "CF",
    "CCl",
    "CBr",
    "CI",
)

# AutoDock atom types for which Vina can compute affinity maps.
# Vina terminates the entire Python process (instead of raising an error) when it
# encounters a ligand atom type without a map (e.g. Si or metals), so ligands with
# other atom types are rejected before they are passed to Vina.
_SUPPORTED_ATOM_TYPES = {
    "C", "A", "N", "NA", "O", "OA", "S", "SA", "P", "F", "Cl", "Br", "I", "H", "HD"
}  # fmt: skip


def dock(
    protein: core.Molecule,
    ligand: core.Molecule,
    center: tuple,
    box_size: tuple,
    mode: str = "dock",
    exhaustiveness: int = 8,
    n_poses: int = 5,
    max_evals: int = 0,
    ncpu: int = None,
    seed: int = 0,
    scoring_function: str = "vina",
    **kwargs,
):
    """
    Dock a ligand to a protein using AutoDock Vina (in-process via the `vina` Python bindings).

    Note
    ----
    When docking many ligands against the same pocket (e.g. in an optimization loop) use
    a `Docker` directly. It keeps the prepared receptor and the affinity maps around so they
    are not re-computed for every call.

    Parameters
    ----------
    protein : Molecule or str
        The protein to dock the ligand to. This may also be the path to an already prepared receptor PDBQT file.
    ligand : Molecule or list of Molecules
        The ligand(s) to dock. RDKit molecules and SMILES strings are also accepted.
        Ligands should be protonated (Vina needs at least the polar hydrogens).
    center : tuple
        The center of the docking box in (x, y, z) coordinates
    box_size : tuple
        The size of the docking box in (x, y, z)
    mode : str, optional
        - "dock": global search for the best pose (default)
        - "optimize": local refinement of the ligand's current pose (fast, the ligand must be in the box)
        - "score": score the ligand's current pose without moving it (fastest, the ligand must be in the box)
    exhaustiveness : int, optional
        The exhaustiveness of the global search (default is 8). Lower values are faster.
    n_poses : int, optional
        The (maximal) number of poses to return when docking (default is 5)
    max_evals : int, optional
        The maximal number of evaluations per Monte Carlo run (default is 0 which uses Vina's heuristic).
        Setting this caps the runtime for large and flexible ligands.
    ncpu : int, optional
        The number of CPUs to use (default is all available). A single ligand is docked with `ncpu` threads,
        multiple ligands are docked in parallel (one ligand per CPU).
    seed : int, optional
        The random seed (default is 0, which lets Vina pick a random seed)
    scoring_function : str, optional
        The scoring function to use, either 'vina' or 'vinardo' (default is 'vina')
    **kwargs
        Additional keyword arguments are passed to `Docker.dock`, e.g. `min_rmsd` or `energy_range`.

    Returns
    -------
    Molecule or list of Molecules
        The docked ligand. The returned molecule has one `model` per pose, sorted by score (best first).
        The best score is available as `docking_score` on the molecule, the score of each pose as `docking_score` on each model.
        If multiple ligands were docked a list of Molecules is returned (with `None` for ligands that failed).

    Examples
    --------
    >>> import buildamol as bam
    >>> from buildamol.extensions.docking import dock
    >>> protein = bam.Molecule.from_pdb("protein.pdb")
    >>> ligand = bam.Molecule.from_molfile("ligand.sdf")
    >>> docked = dock(protein, ligand, center=(10, 12, 5), box_size=(20, 20, 20))
    >>> docked.docking_score
    -7.3
    """
    with Docker(scoring_function=scoring_function, seed=seed) as docker:
        docker.setup(protein, center, box_size)
        return docker.dock(
            ligand,
            mode=mode,
            exhaustiveness=exhaustiveness,
            n_poses=n_poses,
            max_evals=max_evals,
            ncpu=ncpu,
            **kwargs,
        )


class Docker:
    """
    Dock (or score) many ligands against the same protein pocket with AutoDock Vina.

    The receptor is prepared once and the affinity maps are re-used for all ligands
    (they are only re-computed if a ligand introduces a new Vina atom type).
    Use this class directly (rather than the `dock` function) inside optimization loops.

    Parameters
    ----------
    scoring_function : str, optional
        The scoring function to use, either 'vina' or 'vinardo' (default is 'vina')
    seed : int, optional
        The random seed (default is 0, which lets Vina pick a random seed)

    Examples
    --------
    >>> from buildamol.extensions.docking.backend_vina import Docker
    >>> with Docker() as docker:
    ...     docker.setup(protein, center=(10, 12, 5), box_size=(20, 20, 20))
    ...     for candidate in candidates:
    ...         # quick local refinement of a candidate that was built inside the pocket
    ...         refined = docker.optimize(candidate)
    ...         print(refined.docking_score)
    ...     # full docking of a batch of candidates in parallel
    ...     docked = docker.dock(candidates, exhaustiveness=4)
    """

    def __init__(self, scoring_function: str = "vina", seed: int = 0):
        if not aux.has_package("vina") or not aux.has_package("meeko"):
            raise ImportError("Docking requires the `vina` and `meeko` packages")
        if scoring_function not in ("vina", "vinardo"):
            raise ValueError(
                f"Unsupported scoring function '{scoring_function}': use 'vina' or 'vinardo'"
            )
        self.scoring_function = scoring_function
        self.seed = seed
        self.working_directory = None
        self.receptor_file = None
        self.center = None
        self.box_size = None
        self._session = None
        self._pool = None
        self._pool_key = None

    def setup(self, protein: core.Molecule, center: tuple, box_size: tuple, dir=None):
        """
        Set up the docking environment for a specific protein pocket.

        Parameters
        ----------
        protein : Molecule or str
            The protein to dock ligands to. This may also be the path to an already prepared receptor PDBQT file.
        center : tuple
            The center of the docking box in (x, y, z) coordinates
        box_size : tuple
            The size of the docking box in (x, y, z)
        dir : str, optional
            The directory to store the docking files (default is None, which creates a temporary directory that is removed on `close`)
        """
        self.close()
        self.center = [float(i) for i in center]
        self.box_size = [float(i) for i in box_size]
        if dir is None:
            self.working_directory = Path(tempfile.mkdtemp(prefix="bam_vina_"))
            self._owns_directory = True
        else:
            self.working_directory = Path(dir).resolve()
            self.working_directory.mkdir(parents=True, exist_ok=True)
            self._owns_directory = False

        if isinstance(protein, (str, Path)):
            self.receptor_file = str(Path(protein).resolve())
        else:
            from buildamol.utils.pdbqt import _encode_pdbqt_receptor

            pdbqt = _encode_pdbqt_receptor(
                protein, box_center=self.center, box_size=self.box_size
            )
            self.receptor_file = str(self.working_directory / "receptor.pdbqt")
            with open(self.receptor_file, "w") as f:
                f.write(pdbqt)
        return self

    def dock(
        self,
        ligand,
        mode: str = "dock",
        exhaustiveness: int = 8,
        n_poses: int = 5,
        max_evals: int = 0,
        ncpu: int = None,
        min_rmsd: float = 1.0,
        energy_range: float = 3.0,
    ):
        """
        Dock, locally optimize, or score one or more ligands.

        Parameters
        ----------
        ligand : Molecule or list of Molecules
            The ligand(s). RDKit molecules and SMILES strings are also accepted.
        mode : str, optional
            - "dock": global search for the best pose (default)
            - "optimize": local refinement of the ligand's current pose (the ligand must be in the box)
            - "score": score the ligand's current pose without moving it (the ligand must be in the box)
        exhaustiveness : int, optional
            The exhaustiveness of the global search (default is 8). Lower values are faster.
        n_poses : int, optional
            The (maximal) number of poses to return when docking (default is 5)
        max_evals : int, optional
            The maximal number of evaluations per Monte Carlo run (default is 0 which uses Vina's heuristic)
        ncpu : int, optional
            The number of CPUs to use (default is all available). A single ligand is docked with `ncpu` threads,
            multiple ligands are docked in parallel (one ligand per CPU).
        min_rmsd : float, optional
            The minimal RMSD between returned poses (default is 1.0)
        energy_range : float, optional
            The maximal energy difference (kcal/mol) of returned poses to the best pose (default is 3.0)

        Returns
        -------
        Molecule or list of Molecules
            The docked ligand(s) with one `model` per pose (best first) and a `docking_score` attribute.
            If a list of ligands was given, a list is returned that contains `None` for ligands that failed.
        """
        if self.receptor_file is None:
            raise ValueError("Docker must be set up before docking")
        if mode not in MODES:
            raise ValueError(f"Unknown mode '{mode}': use one of {MODES}")

        single = not isinstance(ligand, (list, tuple, np.ndarray))
        ligands = [ligand] if single else list(ligand)
        ncpu = ncpu or os.cpu_count() or 1
        params = dict(
            exhaustiveness=exhaustiveness,
            n_poses=n_poses,
            max_evals=max_evals,
            min_rmsd=min_rmsd,
            energy_range=energy_range,
        )

        prepared = []
        for lig in ligands:
            try:
                p = _prepare_ligand(lig)
                self._check_box(p[2], mode)
                prepared.append(p)
            except Exception as e:
                if single:
                    raise
                warnings.warn(f"Could not prepare ligand {_ligand_name(lig)}: {e}")
                prepared.append(None)

        jobs = [i for i, p in enumerate(prepared) if p is not None]
        results = [None] * len(ligands)
        errors = {}

        # global searches of multiple ligands run in parallel (one ligand per process),
        # everything else is fast enough to just run in-process
        if mode == "dock" and len(jobs) > 1 and ncpu > 1:
            n_workers = min(ncpu, len(jobs))
            pool = self._get_pool(n_workers, max(1, ncpu // n_workers))
            futures = {
                i: pool.submit(_worker_run, prepared[i][2], mode, params) for i in jobs
            }
            from concurrent.futures.process import BrokenProcessPool

            for i, future in futures.items():
                try:
                    results[i] = future.result()
                except BrokenProcessPool as e:
                    errors[i] = e
                    self._pool = None
                    self._pool_key = None
                except Exception as e:
                    errors[i] = e
        else:
            cpu = min(ncpu, exhaustiveness) if mode == "dock" else None
            session = self._get_session(cpu)
            for i in jobs:
                try:
                    results[i] = session.run(prepared[i][2], mode, **params)
                except Exception as e:
                    errors[i] = e

        out = []
        for i, lig in enumerate(ligands):
            if i in errors:
                if single:
                    raise errors[i]
                warnings.warn(f"Docking failed for ligand {_ligand_name(lig)}: {errors[i]}")
                out.append(None)
                continue
            if results[i] is None:
                out.append(None)
                continue
            template, rdmol, pdbqt = prepared[i]
            energies, poses = results[i]
            out.append(_poses_to_molecule(template, rdmol, pdbqt, poses, energies))

        return out[0] if single else out

    def optimize(self, ligand, **kwargs):
        """
        Locally refine the current pose of one or more ligands (which must be placed inside the box).
        This is a shortcut for `dock(ligand, mode="optimize")`.
        """
        return self.dock(ligand, mode="optimize", **kwargs)

    def score(self, ligand, **kwargs):
        """
        Score the current pose of one or more ligands (which must be placed inside the box) without moving them.
        This is a shortcut for `dock(ligand, mode="score")`.
        """
        return self.dock(ligand, mode="score", **kwargs)

    def close(self):
        """
        Shut down any worker processes and remove the temporary working directory.
        """
        if self._pool is not None:
            self._pool.shutdown(cancel_futures=True)
            self._pool = None
            self._pool_key = None
        if self._session is not None:
            self._session.close()
            self._session = None
        if self.working_directory is not None and getattr(
            self, "_owns_directory", False
        ):
            shutil.rmtree(self.working_directory, ignore_errors=True)
        self.working_directory = None
        self.receptor_file = None

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _check_box(self, pdbqt, mode):
        xyz = _pdbqt_pose_coords(pdbqt)[0]
        box = np.asarray(self.box_size)
        suggestion = (
            "Use `box_around(ligand)` to get a box that fits the ligand, "
            "e.g. `docker.setup(protein, *box_around(ligand))`."
        )
        if mode == "dock":
            # compare the extents along the principal axes since the ligand may be rotated
            centered = xyz - xyz.mean(axis=0)
            axes = np.linalg.svd(centered, full_matrices=False)[2]
            extent = np.sort(np.ptp(centered @ axes.T, axis=0))
            if np.any(extent > np.sort(box)):
                warnings.warn(
                    f"The ligand (extent {np.round(extent, 1).tolist()} Å) does not fit into the docking box "
                    f"{box.tolist()} in its current conformation. Docking will be slow and the poses unreliable. "
                    + suggestion
                )
            return
        lower = np.asarray(self.center) - box / 2
        upper = lower + box
        outside = np.any((xyz < lower) | (xyz > upper), axis=1)
        if outside.any():
            raise ValueError(
                f"{outside.sum()} of {len(xyz)} ligand atoms are outside the docking box "
                f"(ligand spans {np.round(xyz.min(0), 1).tolist()} to {np.round(xyz.max(0), 1).tolist()}, "
                f"box spans {np.round(lower, 1).tolist()} to {np.round(upper, 1).tolist()}). "
                f"Mode '{mode}' requires the ligand to lie inside the box. " + suggestion
            )

    def _session_args(self, cpu):
        return (
            self.receptor_file,
            self.center,
            self.box_size,
            self.scoring_function,
            self.seed,
            cpu,
        )

    def _get_session(self, cpu=None):
        # the number of threads is fixed when Vina is created; scoring and local optimization do not care about it
        if self._session is None or (cpu is not None and self._session.cpu != cpu):
            if self._session is not None:
                self._session.close()
            self._session = _VinaSession(*self._session_args(cpu or 1))
        return self._session

    def _get_pool(self, n_workers, cpu_per_worker):
        key = (n_workers, cpu_per_worker)
        if self._pool is None or self._pool_key != key:
            if self._pool is not None:
                self._pool.shutdown(cancel_futures=True)
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor

            self._pool = ProcessPoolExecutor(
                max_workers=n_workers,
                mp_context=multiprocessing.get_context("spawn"),
                initializer=_init_worker,
                initargs=self._session_args(cpu_per_worker),
            )
            self._pool_key = key
        return self._pool


class _VinaSession:
    """
    A Vina instance with a loaded receptor and pre-computed affinity maps.
    """

    def __init__(self, receptor_file, center, box_size, scoring_function, seed, cpu):
        self.cpu = cpu
        self.center = center
        self.box_size = box_size
        self._vina_args = dict(sf_name=scoring_function, cpu=cpu, seed=seed, verbosity=0)
        self._receptor_file = receptor_file
        self._tmpdir = tempfile.mkdtemp(prefix="bam_vina_session_")
        self._map_ligands = [_smiles_to_pdbqt(s, center) for s in _REFERENCE_SMILES]
        self._map_types = set()
        self._compute_maps()

    def _compute_maps(self):
        from vina import Vina

        # Once maps exist, Vina terminates the process when a ligand with uncovered
        # types is set, so extending the maps requires a fresh Vina instance
        v = Vina(**self._vina_args)
        v.set_receptor(rigid_pdbqt_filename=self._receptor_file)
        v.set_ligand_from_string(self._map_ligands)
        v.compute_vina_maps(center=self.center, box_size=self.box_size)
        self._vina = v
        self._map_types = set().union(*(_pdbqt_xs_types(p) for p in self._map_ligands))

    def run(self, pdbqt, mode, exhaustiveness=8, n_poses=5, max_evals=0, min_rmsd=1.0, energy_range=3.0):
        # Vina terminates the entire process (rather than raising) if a ligand has an
        # XS atom type without an affinity map, so maps must be extended beforehand.
        # XS types depend on the bonding context, not only on the AutoDock type
        # (e.g. a carbonyl and a hydroxyl oxygen are both OA but O_A vs O_DA).
        _check_atom_types(pdbqt)
        if not _pdbqt_xs_types(pdbqt) <= self._map_types:
            self._map_ligands.append(pdbqt)
            self._compute_maps()

        v = self._vina
        v.set_ligand_from_string(pdbqt)
        try:
            if mode == "dock":
                v.dock(
                    exhaustiveness=exhaustiveness,
                    n_poses=n_poses,
                    min_rmsd=min_rmsd,
                    max_evals=max_evals,
                )
                energies = v.energies(n_poses=n_poses, energy_range=energy_range)[:, 0]
                poses = v.poses(n_poses=n_poses, energy_range=energy_range)
            elif mode == "optimize":
                energies = np.array([v.optimize()[0]])
                posefile = os.path.join(self._tmpdir, "pose.pdbqt")
                v.write_pose(posefile, overwrite=True)
                with open(posefile) as f:
                    poses = f.read()
            else:
                energies = np.array([v.score()[0]])
                poses = None
        except RuntimeError as e:
            if "outside the grid box" in str(e):
                raise ValueError(
                    f"The ligand is outside the docking box. Use mode='dock' or adjust the box (mode was '{mode}')."
                ) from None
            raise
        return energies, poses

    def close(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)


_WORKER_SESSION = None


def _init_worker(*session_args):
    global _WORKER_SESSION
    _WORKER_SESSION = _VinaSession(*session_args)


def _worker_run(pdbqt, mode, params):
    return _WORKER_SESSION.run(pdbqt, mode, **params)


def _ligand_name(ligand):
    return getattr(ligand, "id", None) or str(ligand)[:40]


def _prepare_ligand(ligand):
    """
    Convert a ligand to (template Molecule, RDKit Mol with hydrogens, PDBQT string).
    The atom order of the template matches the first atoms of the RDKit Mol.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if isinstance(ligand, core.Molecule):
        template = ligand
        rdmol = ligand.to_rdkit()
    elif isinstance(ligand, str):
        rdmol = Chem.MolFromSmiles(ligand)
        if rdmol is None:
            raise ValueError(f"Could not parse SMILES '{ligand}'")
        rdmol = Chem.AddHs(rdmol)
        if AllChem.EmbedMolecule(rdmol, randomSeed=42) != 0:
            raise ValueError(f"Could not embed SMILES '{ligand}'")
        template = None
    elif hasattr(ligand, "GetAtoms"):
        rdmol = Chem.Mol(ligand)
        if rdmol.GetNumConformers() == 0:
            rdmol = Chem.AddHs(rdmol)
            if AllChem.EmbedMolecule(rdmol, randomSeed=42) != 0:
                raise ValueError("Could not embed the RDKit molecule")
        else:
            rdmol = Chem.AddHs(rdmol, addCoords=True)
        template = None
    else:
        raise ValueError(
            f"Unsupported ligand type {type(ligand)}: use a Molecule, RDKit Mol, or SMILES string"
        )

    if not any(a.GetAtomicNum() == 1 for a in rdmol.GetAtoms()):
        warnings.warn(
            "The ligand has no hydrogens. Vina requires (at least) polar hydrogens, consider calling `add_hydrogens()` first."
        )
    if template is None:
        template = core.Molecule.from_rdkit(rdmol)

    pdbqt = _rdkit_to_pdbqt(rdmol)
    _check_atom_types(pdbqt)
    return template, rdmol, pdbqt


def _rdkit_to_pdbqt(rdmol):
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation().prepare(rdmol)
    if len(setups) != 1:
        raise ValueError(
            f"Expected a single connected ligand but got {len(setups)} fragments"
        )
    # the index map links the PDBQT atoms back to the atoms of the input molecule
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0], add_index_map=True)
    if not ok:
        raise ValueError(f"PDBQT preparation failed: {err}")
    return pdbqt


def _smiles_to_pdbqt(smiles, center):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    conf = mol.GetConformer()
    xyz = conf.GetPositions()
    xyz += np.asarray(center) - xyz.mean(axis=0)
    for i, pos in enumerate(xyz):
        conf.SetAtomPosition(i, pos.tolist())
    return _rdkit_to_pdbqt(mol)


def _pdbqt_atom_types(pdbqt):
    return {
        line[77:].strip()
        for line in pdbqt.splitlines()
        if line.startswith(("ATOM", "HETATM"))
    }


# covalent radii used by Vina to perceive bonds (atom_constants.h)
_COVALENT_RADII = {
    "C": 0.77, "A": 0.77, "N": 0.75, "NA": 0.75, "O": 0.73, "OA": 0.73,
    "P": 1.06, "S": 1.02, "SA": 1.02, "H": 0.37, "HD": 0.37, "F": 0.71,
    "I": 1.33, "Cl": 0.99, "Br": 1.14,
}  # fmt: skip
_XS_FIXED = {
    "S": "S_P", "SA": "S_P", "P": "P_P", "F": "F_H",
    "Cl": "Cl_H", "Br": "Br_H", "I": "I_H",
}  # fmt: skip


def _pdbqt_xs_types(pdbqt):
    """
    Get the Vina (XS) atom types of a ligand PDBQT, mirroring `model::assign_types` of Vina.
    """
    types, coords = [], []
    for line in pdbqt.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            types.append(line[77:].strip())
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not types:
        return set()
    coords = np.array(coords)
    radii = np.array([_COVALENT_RADII.get(t, 0.77) for t in types])
    dists = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)
    bonded = dists < 1.1 * (radii[:, None] + radii[None, :])
    np.fill_diagonal(bonded, False)

    xs = set()
    for i, t in enumerate(types):
        if t in ("H", "HD"):
            continue
        neighbors = [types[j] for j in np.flatnonzero(bonded[i])]
        if t in ("C", "A") or re.fullmatch(r"CG\d", t):
            hetero = any(n not in ("C", "A", "H", "HD") for n in neighbors)
            x = "C_P" if hetero else "C_H"
            if t.startswith("CG"):
                x += "_" + t
            xs.add(x)
        elif t in ("N", "NA", "O", "OA"):
            acceptor = t in ("NA", "OA")
            donor = "HD" in neighbors
            suffix = "DA" if acceptor and donor else "A" if acceptor else "D" if donor else "P"
            xs.add(f"{t[0]}_{suffix}")
        else:
            xs.add(_XS_FIXED.get(t, t))
    return xs


def _check_atom_types(pdbqt):
    types = _pdbqt_atom_types(pdbqt)
    # CG0/G0 etc. are pseudo atoms that meeko introduces to keep macrocycles flexible
    unsupported = {
        t for t in types - _SUPPORTED_ATOM_TYPES if not re.fullmatch(r"C?G\d", t)
    }
    if unsupported:
        raise ValueError(
            f"The ligand contains atom types that are not supported by Vina: {sorted(unsupported)}"
        )
    return types


def _pdbqt_index_map(pdbqt):
    """
    Get the mapping of PDBQT atom numbers (1-based, in file order) to input atom indices (0-based)
    from the `REMARK INDEX MAP` lines written by meeko.
    """
    mapping = {}
    for line in pdbqt.splitlines():
        if line.startswith("REMARK INDEX MAP"):
            values = [int(i) for i in line[len("REMARK INDEX MAP") :].split()]
            for input_idx, pdbqt_number in zip(values[::2], values[1::2]):
                mapping[pdbqt_number] = input_idx - 1
    return mapping


def _pdbqt_pose_coords(poses):
    """
    Get the coordinates of each pose (in PDBQT atom order) from a (multi-model) PDBQT string.
    """
    blocks, current = [], []
    for line in poses.splitlines():
        if line.startswith("MODEL"):
            current = []
        elif line.startswith(("ATOM", "HETATM")):
            current.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        elif line.startswith("ENDMDL"):
            blocks.append(np.array(current))
            current = []
    if current:
        blocks.append(np.array(current))
    return blocks


def _kabsch(source, target):
    """
    Get the rotation and translation that superimpose `source` onto `target`.
    """
    source_center, target_center = source.mean(axis=0), target.mean(axis=0)
    H = (source - source_center).T @ (target - target_center)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, target_center - R @ source_center


def _map_pose(rdmol, ref_xyz, pose_xyz, index_map):
    """
    Map PDBQT pose coordinates back onto the atoms of the input RDKit molecule.
    Hydrogens that meeko merged into their parent atom are placed by moving them rigidly
    with their parent's local environment (bond lengths and angles are not changed by Vina).
    """
    xyz = np.full((rdmol.GetNumAtoms(), 3), np.nan)
    for number, input_idx in index_map.items():
        xyz[input_idx] = pose_xyz[number - 1]

    for idx in np.flatnonzero(np.isnan(xyz[:, 0])):
        atom = rdmol.GetAtomWithIdx(int(idx))
        if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
            raise ValueError("Could not map the docked poses back onto the ligand")
        parent = atom.GetNeighbors()[0]
        frame = [parent.GetIdx()]
        frame += [n.GetIdx() for n in parent.GetNeighbors() if not np.isnan(xyz[n.GetIdx(), 0])]
        if len(frame) < 3:
            # e.g. methyl groups: include the next shell to fix the rotation around the bond
            for i in list(frame[1:]):
                frame += [
                    n.GetIdx()
                    for n in rdmol.GetAtomWithIdx(i).GetNeighbors()
                    if n.GetIdx() not in frame and not np.isnan(xyz[n.GetIdx(), 0])
                ]
        frame = np.array(frame)
        if np.isnan(xyz[frame[0], 0]):
            raise ValueError("Could not map the docked poses back onto the ligand")
        R, t = _kabsch(ref_xyz[frame], xyz[frame])
        xyz[idx] = R @ ref_xyz[idx] + t
    return xyz


def _poses_to_molecule(template, rdmol, pdbqt, poses, energies):
    """
    Make a copy of the template molecule with one model per pose.
    """
    out = template.copy()
    # keep only the active conformation of the template, poses are added as new models
    for model in list(out.models)[1:]:
        out.remove_model(model)
    out._active_conformer_id = out._model.id

    atoms = list(out.get_atoms())
    if poses is None:
        coords = [np.array([a.coord for a in atoms])]
    else:
        index_map = _pdbqt_index_map(pdbqt)
        ref_xyz = rdmol.GetConformer().GetPositions()
        coords = [
            _map_pose(rdmol, ref_xyz, pose_xyz, index_map)[: len(atoms)]
            for pose_xyz in _pdbqt_pose_coords(poses)
        ]

    for i, xyz in enumerate(coords):
        for atom, pos in zip(atoms, xyz):
            atom.coord = np.asarray(pos, dtype=float)
        if i > 0:
            out.add_model()
            out.models[-1].docking_score = float(energies[i])
    # the canonical model holds the best pose
    for atom, pos in zip(atoms, coords[0]):
        atom.coord = np.asarray(pos, dtype=float)
    out.models[0].docking_score = float(energies[0])
    out.docking_score = float(energies[0])
    return out


def box_around(ligand, padding: float = 4.0):
    """
    Get a docking box that encloses a ligand.

    Parameters
    ----------
    ligand : Molecule or array-like
        The ligand (or its coordinates) to enclose
    padding : float, optional
        The space (in Å) to add on each side of the ligand (default is 4.0)

    Returns
    -------
    tuple
        The center and the size of the box

    Examples
    --------
    >>> docker = Docker().setup(protein, *box_around(ligand))
    """
    xyz = ligand.get_coords() if hasattr(ligand, "get_coords") else np.asarray(ligand)
    lower, upper = xyz.min(axis=0), xyz.max(axis=0)
    center = (lower + upper) / 2
    size = upper - lower + 2 * padding
    return tuple(center.tolist()), tuple(size.tolist())


__dock__ = dock
__all__ = ["dock", "Docker", "MODES", "box_around"]
