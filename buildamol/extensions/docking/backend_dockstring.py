"""
Docking module for BuildAMol using the `dockstring` package.
"""

import buildamol.core as core
import buildamol.utils.auxiliary as aux

from tempfile import TemporaryDirectory
from pathlib import Path
from copy import deepcopy

def dock(protein: core.Molecule, ligand: core.Molecule, center: tuple, box_size: tuple, ph: float = 7.4, num_cpus: int = 1):
    """
    Dock a ligand to a protein using the `dockstring` package.

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
    ph : float, optional
        The pH to use for docking (default is 7.4)
    num_cpus : int, optional
        The number of CPUs to use (default is 1)
    
    Returns
    -------
    Molecule or list of Molecules
        The docked ligand(s). For each ligand, a Molecule with the docked poses (models) is returned. Each model comes with a `docking_score` attribute.
        If only one ligand is docked, a single Molecule is returned. If a list of ligands is docked, a list of Molecules is returned.
    """
    docker = Docker()
    docker.setup(protein)
    return docker.dock(ligand, center, box_size, ph, num_cpus)

class Docker:
    """
    A class for docking molecules using the `dockstring` package.
    """
    def __init__(self):
        if not aux.has_package("dockstring") or not aux.has_package("openbabel"):
            raise ImportError("Docking requires `dockstring` and `openbabel` to be installed")
        self.working_direrctory = None
        self.protein_id = None

        import dockstring as ds
        self._target = ds.target.Target
        
    def setup(self, protein: core.Molecule, dir=None):
        """
        Setup a protein for docking
        
        Parameters
        ----------
        protein : Molecule
            The protein to setup
        """
        outdir = _prepare_outdir(dir)
        protein.to_pdbqt(outdir / (protein.id + "_target.pdbqt"))
        self.working_direrctory = outdir
        self.protein_id = protein.id
        self.target = self._target(name=protein.id, working_dir=outdir)

    def dock(self, ligand: core.Molecule, center: tuple, box_size: tuple, ph: float = 7.4, num_cpus: int = 1):
        """
        Run docking of one or more ligand(s) to a protein

        Parameters
        ----------
        ligand : Molecule or list of Molecules
            The ligand(s) to dock. This can be a single molecule or a list of molecules
        center : tuple
            The center of the docking box in (x, y, z) coordinates
        box_size : tuple
            The size of the docking box in (x, y, z)
        ph : float, optional    
            The pH to use for docking (default is 7.4)
        num_cpus : int, optional
            The number of CPUs to use (default is 1)
        
        Returns
        -------
        Molecule or list of Molecules
            The docked ligand(s). For each ligand, a Molecule with the docked poses (models) is returned. Each model comes with a `docking_score` attribute.
            If only one ligand is docked, a single Molecule is returned. If a list of ligands is docked, a list of Molecules is returned.
        """
        if isinstance(ligand, str):
            pass
        elif isinstance(ligand, core.Molecule):
            ligand = ligand.to_smiles()
        elif "Chem.Mol" in str(type(ligand)):
            from rdkit import Chem
            ligand = Chem.MolToSmiles(ligand)
        elif isinstance(ligand, (list, tuple, set)):
            return [self.dock(lig, center, box_size, ph, num_cpus) for lig in ligand]

        if self.working_direrctory is None:
            raise ValueError("No protein setup")
        _write_conf(self.working_direrctory, self.protein_id, center, box_size)
        out = self.target.dock(ligand, pH=ph, num_cpus=num_cpus)

        # postprocess the results
        best_docking_score, data = out
        docked = _get_ligands_from_rdkit(data["ligand"])
        affinities = data["affinities"]
        for i, pose in enumerate(docked.models):
            pose.docking_score = affinities[i]
        return docked

def _get_ligands_from_rdkit(rdkit_ligand):
    mols = []
    for i in range(rdkit_ligand.GetNumConformers()):
        lig = deepcopy(rdkit_ligand)
        conf_ids = list(range(lig.GetNumConformers()))
        conf_ids.remove(i)
        for conf_id in conf_ids:
            lig.RemoveConformer(conf_id)

        new = core.Molecule.from_rdkit(lig)
        mols.append(new)
    
    final = mols[0]
    for mol in mols[1:]:
        final.add_model(mol.model)
    return final

def _write_conf(dir, protein_id, center, box_size):
    with open(dir / f"{protein_id}_conf.txt", "w") as f:
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
__all__ = ["dock", "Docker"]