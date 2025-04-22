"""
Solvate a molecular system with water using PDBFixer.
"""

from io import StringIO
import buildamol.utils.auxiliary as aux
import buildamol.utils.pdb as pdb

if not aux.has_package("pdbfixer"):
    raise ImportError(
        "PDBFixer is not available. Please install PDBFixer to use this feature."
    )

from pdbfixer import PDBFixer
from openmm.unit import molar, nanometer
from openmm.vec3 import Vec3
from openmm.app import PDBFile


class Ions:
    """
    Available ions for solvation in PDBFixer.
    """

    Cs = "Cs+"
    """
    Cesium cation (Cs+)
    """
    K = "K+"
    """
    Potassium cation (K+)
    """

    Li = "Li+"
    """
    Lithium cation (Li+)
    """

    Na = "Na+"
    """
    Sodium cation (Na+)
    """

    Rb = "Rb+"
    """
    Rubidium cation (Rb+)
    """

    Cl = "Cl-"
    """
    Chloride anion (Cl-)
    """

    Br = "Br-"
    """
    Bromide anion (Br-)
    """

    F = "F-"
    """
    Fluoride anion (F-)
    """

    I = "I-"
    """
    Iodide anion (I-)
    """

    def __repr__(self):
        return f"Ions({self.__dict__})"


def solvate(
    mol: "Molecule",
    box_shape: str = "cube",
    box_padding: float = 10.0,
    ions: list = None,
    ionic_strength: float = 0.0,
    symmetric_box: bool = True,
    outfile: str = None,
    return_molecule: bool = True,
):
    """
    Solvate a molecular system with water using PDBFixer.

    Parameters
    ----------
    mol : Molecule
        The molecule to solvate.
    box_shape : str
        The shape of the box to use for solvation. Can be any of 'cube', 'dodecahedron', and 'octahedron'.
    box_padding : float
        The padding to add to the box size in angstroms.
    ions : list
        The list of ions to use for solvation. Can contain any two of 'Na+', 'Cl-', 'K+', 'Br-', 'Cs+', 'Li+', 'Rb+', 'I-', and 'F-'.
        Or use the Ions enum class to get the available ions. By default, the list will be ['Na+', 'Cl-'].
    ionic_strength : float
        The ionic strength of the solution in molar.
    symmetric_box : bool
        Whether to use a symmetric box or not. If True, the box size will be the same in all directions. If False, the box size will be determined by the molecule's bounding box.
    outfile : str
        The file to write the solvated molecule to. If None, the molecule will not be written to a file.
    return_molecule : bool
        Whether to return the solvated molecule or not. If True, the solvated molecule will be returned. If False, the function will return None.

    Returns
    -------
    Molecule
        The solvated molecule, which will contain a new chain with the water molecules. Water molecules are named HOH and do not have any bonds.
    """
    pdb_string = StringIO(pdb.encode_pdb(mol))
    box_size = _determine_box_size(mol)
    if symmetric_box:
        box_size[:] = max(box_size)
        box_size += box_padding
    else:
        box_size += box_padding
        box_size = Vec3(*box_size)
    fixer = PDBFixer(
        pdbfile=pdb_string,
    )
    if ions is None:
        ions = [Ions.Na, Ions.Cl]

    positive_ion = [ion for ion in ions if ion.endswith("+")]
    positive_ion = positive_ion[0] if len(positive_ion) > 0 else None
    negative_ion = [ion for ion in ions if ion.endswith("-")]
    negative_ion = negative_ion[0] if len(negative_ion) > 0 else None
    fixer.addSolvent(
        boxShape=box_shape,
        boxSize=box_size / 10 * nanometer,
        ionicStrength=ionic_strength * molar,
        positiveIon=positive_ion,
        negativeIon=negative_ion,
    )

    if outfile is not None:
        with open(outfile, "w") as f:
            PDBFile.writeFile(
                fixer.topology,
                fixer.positions,
                f,
            )
    if return_molecule:
        tmp = StringIO()
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            tmp,
        )
        tmp.seek(0)
        out = mol.__class__._from_pdb_string(tmp.read())
        tmp.close()
        return out


def _determine_box_size(
    mol,
):
    """
    Determine the box size for the solvation box.
    """
    # Get the coordinates of the molecule
    coords = mol.get_coords()
    # Get the bounding box of the molecule
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    # Calculate the box size
    box_size = max_coords - min_coords
    return box_size


if __name__ == "__main__":
    import buildamol as bam

    bam.load_sugars()
    mol = bam.molecule("d-glucose")[0]

    solvated_mol = solvate(
        mol,
        box_shape="cubic",
        box_padding=5,
        ions=[Ions.Na, Ions.Cl],
        ionic_strength=0.1,
    )
    pass
