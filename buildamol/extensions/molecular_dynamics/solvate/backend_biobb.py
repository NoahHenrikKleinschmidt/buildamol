import buildamol.utils.auxiliary as aux
from buildamol.core import Molecule
from pathlib import Path
from uuid import uuid4


if not aux.has_package("biobb_amber"):
    raise ImportError(
        "biobb_amber is not available. Please install biobb_amber to use this feature."
    )

from biobb_amber.leap.leap_solvate import leap_solvate


def solvate(
    mol: "Molecule",
    properties: dict,
    outfile_prefix: str = None,
    return_molecule: bool = True,
    **kwargs,
) -> "Molecule":
    """
    Solvate a molecule using the leap_solvate function from biobb_amber.
    This function uses the leap_solvate function from biobb_amber to solvate a molecule with water and ions.
    The leap_solvate function is a wrapper around the tleap program from AmberTools.

    Parameters
    ----------
    mol : Molecule
        The molecule to solvate
    properties : dict
        The properties to use for solvation. See https://biobb-amber.readthedocs.io/en/latest/leap.html#leap.leap_solvate.LeapSolvate
        for more information on the available properties, possible values, and default settings.
    outfile_prefix : str
        The prefix to use for the output files. The output files will be named
        {outfile_prefix}.pdb, {outfile_prefix}.top, and {outfile_prefix}.crd.
        If None, no output files are written.
    **kwargs
        Any additional keyword arguments are interpreted as "properties" and will be merged with the dictionary passed in the properties parameter.
    """

    # Check if the molecule is a valid input
    if not isinstance(mol, Molecule) or not hasattr(mol, "to_pdb"):
        raise ValueError(
            "The input molecule must be of type Molecule or at least have a 'to_pdb' method that will write a pdb file of the molecule structure."
        )

    if not return_molecule and outfile_prefix is None:
        raise ValueError(
            "Either outfile_prefix must be set or return_molecule must be True, otherwise this function call is pointless."
        )

    if properties is None:
        properties = {}
    properties = {**properties, **kwargs}

    if len(properties) == 0:
        properties = None

    if outfile_prefix is None:
        should_remove_output_files = True
        # Generate a unique name for the output files
        outfile_prefix = str(uuid4())
    else:
        should_remove_output_files = False

    outfile_prefix = Path(outfile_prefix).stem
    input_pdb = outfile_prefix.with_suffix(".input.pdb")
    out_pdb = outfile_prefix.with_suffix(".pdb")
    out_top = outfile_prefix.with_suffix(".top")
    out_crd = outfile_prefix.with_suffix(".crd")

    mol.to_pdb(input_pdb)

    leap_solvate(
        input_pdb_path=str(input_pdb),
        output_pdb_path=str(out_pdb),
        output_top_path=str(out_top),
        output_crd_path=str(out_crd),
        properties=properties,
    )

    if return_molecule:
        out = Molecule.Molecule.from_pdb(str(out_pdb))
    else:
        out = None

    if should_remove_output_files:
        # Remove the output files
        for file in [input_pdb, out_pdb, out_top, out_crd]:
            if file.exists():
                file.unlink()

    return out
