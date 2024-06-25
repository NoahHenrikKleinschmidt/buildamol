import buildamol.resources.charmm as charmm
import buildamol.resources.pdbe_compounds as pdbe_compounds
from pathlib import Path


def export_custom_resources(directory_name: str):
    """
    Export any custom default settings for CHARMM topology and PDBE compounds to XML files in the specified directory.

    Parameters
    ----------
    directory_name : str
        The directory to which the XML files will be written
    """
    pdbe_compounds.export_compounds(directory_name + "/pdbe_compounds.xml")
    charmm.export_topology(directory_name + "/charmm_topology.xml")


def import_custom_resources(directory_name: str, overwrite: bool = False):
    """
    Import any custom default settings for CHARMM topology and PDBE compounds from XML or JSON backup files in the specified directory.

    Parameters
    ----------
    directory_name : str
        The directory containing the XML or JSON backup files.
        This directory must contain a "pdbe_compounds.xml" or
        "pdbe_compounds.json" file and a "charmm_topology.xml"
        or "charmm_topology.json" file.
    overwrite : bool
        Whether to overwrite the default settings permanently with the custom settings (all future sessions).
    """
    directory = Path(directory_name)
    compounds = directory / "pdbe_compounds.xml"
    topology = directory / "charmm_topology.xml"

    if not compounds.exists():
        compounds = directory / "pdbe_compounds.json"
        if not compounds.exists():
            raise FileNotFoundError(
                "No PDBe compounds file in JSON or XML found in the specified directory"
            )

    if not topology.exists():
        topology = directory / "charmm_topology.json"
        if not topology.exists():
            raise FileNotFoundError(
                "No CHARMM topology file in JSON or XML found in the specified directory"
            )

    compounds = pdbe_compounds.read_compounds(compounds)
    topology = charmm.read_topology(topology)
    if overwrite:
        pdbe_compounds.set_default_compounds(compounds, overwrite=True)
        charmm.set_default_topology(topology, overwrite=True)


__all__ = ["export_custom_resources", "import_custom_resources"]
