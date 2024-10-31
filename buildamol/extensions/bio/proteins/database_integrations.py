"""
Functions to query databases for protein structures.
"""

import buildamol.core as core
import requests
from pathlib import Path

__all__ = ["query_rcsb", "query_pdbe", "query_alphafold"]


def query_rcsb(id: str) -> core.Molecule:
    """
    Query the RCSB database for a protein structure.

    Parameters
    ----------
    id : str
        The ID of the protein to query.

    Returns
    -------
    core.Molecule
        The protein structure.
    """
    url = f"https://files.rcsb.org/download/{id}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    return core.Molecule._from_pdb_string(response.text, id=id)


def query_pdbe(id: str) -> core.Molecule:
    """
    Query the PDBe database for a protein structure.

    Parameters
    ----------
    id : str
        The ID of the protein to query.

    Returns
    -------
    core.Molecule
        The protein structure.
    """
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/download/{id}.cif"
    response = requests.get(url)
    response.raise_for_status()
    ciffile = Path(f"{id}.cif")
    with open(ciffile, "w") as f:
        f.write(response.text)
    mol = core.Molecule.from_cif(ciffile)
    ciffile.unlink()
    return mol


def query_alphafold(id: str, model: str = "F1", version: int = 4) -> core.Molecule:
    """
    Query the AlphaFold database for a protein structure.

    Parameters
    ----------
    id : str
        The ID of the protein to query.
    model : str, optional
        The model to query, by default "F1".
    version : int, optional
        The version of the model to query, by default 4.

    Returns
    -------
    core.Molecule
        The protein structure.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-{model}-model_v{version}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    return core.Molecule._from_pdb_string(response.text, id=id)
