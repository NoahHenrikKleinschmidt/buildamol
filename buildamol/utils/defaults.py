import os
import Bio.PDB as bio

import buildamol.utils.constants as constants

# =================================================================
# Default settings of values
# =================================================================

DEFAULT_BOND_LENGTH = 1.6
"""
The default length of a bond in Angstrom
"""

DEFAULT_CHARMM_TOPOLOGY_FILE = os.path.join(constants.RESOURCES, "topology.pkl")
"""
The path to the default CHARMM topology file
"""


DEFAULT_PDBE_COMPONENT_FILES = {
    "base": os.path.join(constants.RESOURCES, "components_base.pkl"),
    "lipids": os.path.join(constants.RESOURCES, "components_lipids.pkl"),
    "sugars": os.path.join(constants.RESOURCES, "components_sugars.pkl"),
    "nucleotides": os.path.join(constants.RESOURCES, "components_nucleotides.pkl"),
    "amino_acids": os.path.join(constants.RESOURCES, "components_amino_acids.pkl"),
    "small_molecules": os.path.join(constants.RESOURCES, "components_small.pkl"),
}
"""
The paths to the default PDBe component files
"""

DEFAULT_SASA_PROBE_RADIUS = 1.4
"""
The default probe radius for calculating solvent accessible surface area
"""

DEFAULT_SASA_N = 300
"""
The default number of points used to calculate solvent accessible surface area
"""

# =================================================================
# Default instances of auxiliary classes
# =================================================================


__default_instances__ = dict(
    bioPDBParser=bio.PDBParser(),
    bioMMCIFParser=bio.MMCIFParser(),
    bioSASA=bio.SASA.ShrakeRupley(DEFAULT_SASA_PROBE_RADIUS, DEFAULT_SASA_N),
)
"""
Default instance dictionary
"""

__bioPDBParser__ = __default_instances__["bioPDBParser"]
"""
The default instance of Bio.PDB.PDBParser
"""

# =================================================================
# Auxiliary functions
# =================================================================


def get_default_instance(key):
    """
    Get the default instance of a class

    Parameters
    ----------
    key : str
        The key of the default instance

    Returns
    -------
    obj
        The default instance of the class
    """
    return __default_instances__.get(key, None)


def set_default_instance(key, obj):
    """
    Set the default instance of a class

    Parameters
    ----------
    key : str
        The key of the default instance
    obj
        The new default instance
    """
    __default_instances__[key] = obj


# =================================================================
# Default values
# =================================================================
