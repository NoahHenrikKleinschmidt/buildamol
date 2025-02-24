"""
Docking module for BuildAMol using the `easydock` package.
"""

import buildamol.core as core
import buildamol.utils.auxiliary as aux

from tempfile import TemporaryDirectory
from pathlib import Path
import yaml
import numpy as np

EASYDOCK_VINA_BACKEND = "vina"

def set_vina_backend(backend: str):
    """
    Set the default vina backend to use

    Parameters
    ----------
    backend : str
        The backend to use, either 'vina' or 'gnina'
    """
    global EASYDOCK_VINA_BACKEND
    EASYDOCK_VINA_BACKEND = backend

def get_vina_backend():
    """
    Get the default vina backend

    Returns
    -------
    str
        The default vina backend
    """
    return EASYDOCK_VINA_BACKEND

def dock(protein: core.Molecule, ligand: core.Molecule, center: tuple, box_size: tuple, backend: str = None, n_poses: int = 5, exhaustiveness: int = 8, ncpu: int = 1, **kwargs):
    """
    Dock a ligand to a protein using AutoDock Vina via the `easydock` package.

    Parameters
    ----------
    protein : Molecule
        The protein to dock the ligand to
    ligand : Molecule or list of Molecules
        The ligand to dock. This can be a single molecule or a list of molecules
    center : tuple
        The center of the docking box in (x, y, z) coordinates
    box_size : tuple
        The size of the docking box in (x, y, z)
    backend : str, optional
        The docking backend to use, either 'vina' or 'gnina' (default is 'vina')
    n_poses : int, optional
        The number of poses to generate (default is 5)
    exhaustiveness : int, optional
        The exhaustiveness of the docking search (default is 8)
    ncpu : int, optional
        The number of CPUs to use (default is 1)
    **kwargs
        Additional keyword arguments to pass to the docking backend

    Returns
    -------
    Molecule or list of Molecules
        The docked ligand. The molecule that is returned has multiple `models`, one for each pose. If multiple ligands were docked, a list of Molecules is returned.

    Examples
    --------
    Assuming we want to dock a single ligand to a protein to a binding site that is roughly in the center of the protein:
    >>> import buildamol as bam
    >>> from buildamol.extensions.docking import dock
    >>> protein = bam.Molecule.from_pdb("protein.pdb")
    >>> ligand = bam.Molecule.from_molfile("ligand.sdf")
    >>> where = protein.center_of_geometry # or some more intelligent choice...
    >>> docked = dock(protein, ligand, where, box_size=(10, 10, 10))
    >>> print(docked)
    Molecule(id='ligand.sdf')
    >>> print(docked.models)
    [Model(1), Model(2), Model(3), Model(4), Model(5)]
    """
    if not aux.has_package("easydock"):
        raise ImportError("Docking requires the `easydock` package")
    if backend is None:
        backend = EASYDOCK_VINA_BACKEND
    docker = Docker(backend)
    docker.setup(protein, exhaustiveness=exhaustiveness, n_poses=n_poses, ncpu=ncpu, **kwargs)
    return docker.dock(ligand, center, box_size)

class Docker:
    """
    A class to dock ligands to a protein using AutoDock Vina via the `easydock` package.

    Parameters
    ----------
    backend : str, optional
        The docking backend to use, either 'vina' or 'gnina' (default is 'vina')
    """

    def __init__(self, backend: str = "vina"):
        self.backend = backend
        self.working_directory = None
        if backend == "vina":
            from easydock.vina_dock import mol_dock
        elif backend == "gnina":
            from easydock.gnina_dock import mol_dock
        else:
            raise ValueError(
                f"Backend {backend} not supported by `easydock`: use 'vina' or 'gnina'"
            )
        self._core_docking_fn = mol_dock
        from easydock.run_dock import docking

        self._docking_fn = docking

    def setup(
        self,
        protein: core.Molecule,
        exhaustiveness: int = 8,
        seed: int = 0,
        n_poses: int = 5,
        ncpu: int = 1,
        dir=None,
        **kwargs,
    ):
        """
        Set up the docking environment for a specific protein

        Parameters
        ----------
        protein : Molecule
            The protein to dock ligands to
        exhaustiveness : int, optional
            The exhaustiveness of the docking search (default is 8)
        seed : int, optional
            The seed for the random number generator (default is 0)
        n_poses : int, optional
            The number of poses to generate (default is 5)
        ncpu : int, optional
            The number of CPUs to use (default is 1)
        dir : str, optional
            The directory to store the docking files (default is None, which creates a temporary directory)
        **kwargs
            Additional keyword arguments to pass to the docking backend
        """
        outdir = _prepare_outdir(dir)
        proteinfile = outdir / "protein.pdbqt"
        protein.to_pdbqt(proteinfile)
        _write_easydock_config(exhaustiveness, seed, n_poses, ncpu, outdir, **kwargs)
        self.working_directory = outdir

    def dock(self, ligand: core.Molecule, center: tuple, box_size: tuple):
        """
        Dock a ligand to the protein

        Parameters
        ----------
        ligand : Molecule or list of Molecules
            The ligand(s) to dock
        center : tuple
            The center of the docking box in (x, y, z) coordinates
        box_size : tuple
            The size of the docking box in (x, y, z).

        Returns
        -------
        Molecule or list of Molecules
            The docked ligand(s). If one ligand was docked only one Molecule is returned, otherwise a list of Molecules is returned.
            Each molecule that is returned as multiple `models`, one for each pose.
        """
        if self.working_directory is None:
            raise ValueError("Docker must be set up before docking")

        _write_gridfile(center, box_size, self.working_directory)
        config = self.working_directory / "config.yaml"
        ncpu = yaml.safe_load(open(config))["ncpu"]

        if isinstance(ligand, (list, tuple, np.ndarray)) and hasattr(
            ligand[0], "to_rdkit"
        ):
            ligand = [l.to_rdkit() for l in ligand]
        elif hasattr(ligand, "to_rdkit"):
            ligand = ligand.to_rdkit()
        elif "rdkit" in str(type(ligand)):
            pass
        else:
            raise ValueError(
                "Ligand must be a BuildAMol Molecule or a list of Molecules"
            )

        if not isinstance(ligand, (list, tuple, np.ndarray)):
            ligand = [ligand]

        docked = list(
            self._docking_fn(
                ligand,
                dock_func=self._core_docking_fn,
                dock_config=str(self.working_directory / "config.yaml"),
                ncpu=ncpu,
            )
        )

        # postprocess docked poses
        out = []
        for id, data in docked:
            pdb_string = data["pdb_block"]
            with open(self.working_directory / f"{id}.pdb", "w") as f:
                f.write(pdb_string)
            new = core.Molecule.from_pdb(
                    self.working_directory / f"{id}.pdb", id=id, model="all"
                )
            new.docking_score = data["docking_score"]
            out.append(new)

        # if only one ligand was docked, reduce the dimensionality of the output
        if len(out) == 1:
            out = out[0]

        # remove the temporary directory
        self.working_directory.rmdir()
        return out


def _write_easydock_config(exhaustiveness, seed, n_poses, ncpu, outdir, **kwargs):
    configfile = outdir / "config.yaml"
    with open(configfile, "w") as f:
        f.write(f"protein: {outdir}/protein.pdbqt\n")
        f.write(f"protein_setup: {outdir}/grid.txt\n")
        f.write(f"exhaustiveness: {exhaustiveness}\n")
        f.write(f"seed: {seed}\n")
        f.write(f"n_poses: {n_poses}\n")
        f.write(f"ncpu: {ncpu}\n")
        for key, value in kwargs.items():
            f.write(f"{key}: {value}\n")


def _write_gridfile(center, box_size, outdir):
    with open(outdir / "grid.txt", "w") as f:
        f.write(f"center_x = {center[0]}\n")
        f.write(f"center_y = {center[1]}\n")
        f.write(f"center_z = {center[2]}\n")
        f.write(f"size_x = {box_size[0]}\n")
        f.write(f"size_y = {box_size[1]}\n")
        f.write(f"size_z = {box_size[2]}\n")


def _prepare_outdir(dir):
    if dir is None:
        dir = Path(TemporaryDirectory(dir=".").name)
    else:
        dir = Path(dir)
    dir.mkdir(exist_ok=True)
    return dir


__dock__ = dock
__all__ = ["dock", "set_vina_backend", "get_vina_backend", "Docker"]