import biobuild.resources.charmm as charmm
import biobuild.resources.pdbe_compounds as pdbe_compounds


def export_custom_resources(directory_name: str):
    """
    Export any custom default settings for CHARMM topology and PDBE compounds to JSON files in the specified directory.

    Parameters
    ----------
    directory_name : str
        The directory to which the JSON files will be written
    """
    pdbe_compounds.export_compounds(directory_name + "/pdbe_compounds.json")
    charmm.export_topology(directory_name + "/charmm_topology.json")


def import_custom_resources(directory_name: str, overwrite: bool = False):
    """
    Import any custom default settings for CHARMM topology and PDBE compounds from JSON files in the specified directory.

    Parameters
    ----------
    directory_name : str
        The directory containing the JSON files
    overwrite : bool
        Whether to overwrite the default settings permanently with the custom settings (all future sessions).
    """
    compounds = pdbe_compounds.read_compounds(directory_name + "/pdbe_compounds.json")
    topology = charmm.read_topology(directory_name + "/charmm_topology.json")
    if overwrite:
        pdbe_compounds.set_default_compounds(compounds, overwrite=True)
        charmm.set_default_topology(topology, overwrite=True)
