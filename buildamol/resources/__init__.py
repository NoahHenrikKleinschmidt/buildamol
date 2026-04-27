"""
External data resource integrations to BuildAMol. 

PDBe Compounds
--------------

The Protein Data Bank in Europe (PDBe) provides a database of small molecules that are found in PDB structures.
They are available via the `PDBe Component Library <https://www.ebi.ac.uk/pdbe/pdb-component-library/>`_ . BuildAMol
integrates a subset of the database directly for Molecule creation. Since this database is integrated directly in BuildAMol it can be used offline.

PubChem
-------

PubChem maintains an extensive database of small molecules. BuildAMol integrates the PubChem API for Molecule creation.
Since this database is queried online and not integrated directly in BuildAMol it requires an internet connection.
PubChem follows a different atom labelling convention than the CHARMM force field. Hence, many compounds may not be compatible
with the pre-set patches in BuildAMol and may thus not work at all or produce erroneous results. To assist in importing 
PubChem compounds anyway, BuildAMol provides a function ``pubchem_to_cif`` that can be used to convert PubChem compound
data into a CIF file which can be more easily edited and subsequently loaded into BuildAMol using the `PDBECompounds` class
or using ``Molecule.from_cif()``. Also methods such as ``Molecule.autolabel()`` are designed to assist in the conversion of PubChem
compounds into CHARMM-compatible molecules.

CHARMM
------

CHARMM is a molecular simulation software that is widely used in the field of computational chemistry.
It is developed and maintained by the `CHARMM Development Project <http://www.charmm.org/>`_.
BuildAMol integrates the CHARMM force field for pre-defined molecular linkages.
"""

import buildamol.resources.pdbe_compounds as pdbe_compounds
import buildamol.resources.pubchem as pubchem
import buildamol.resources.charmm as charmm
import buildamol.resources.defaults as defaults

from buildamol.resources.pdbe_compounds import *
from buildamol.resources.charmm import *
from buildamol.resources.defaults import *

# Keep CHARMM functionality available under `buildamol.resources.charmm`
# while avoiding default star-export into the top-level `buildamol` namespace.
__all__ = [
	"pdbe_compounds",
	"pubchem",
	*pdbe_compounds.__all__,
	*defaults.__all__,
]
