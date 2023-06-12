"""
External data resource integrations to _biobuild_. 

PDBe Compounds
==============

The Protein Data Bank in Europe (PDBe) provides a database of small molecules that are found in PDB structures.
They are available via the `PDBe Component Library <https://www.ebi.ac.uk/pdbe/pdb-component-library/>`. _biobuild_
integrates a subset of the database directly for Molecule creation. Since this database is integrated directly in biobuild it can be used offline.

PubChem
=======

PubChem maintains an extensive database of small molecules. _biobuild_ integrates the PubChem API for Molecule creation.
Since this database is queried online and not integrated directly in biobuild it requires an internet connection.
PubChem follows a different atom labelling convention than the CHARMM force field. Hence, many compounds may not be compatible
with the pre-set `Patches` in _biobuild_ and may thus not work at all or produce erroneous results. To assist in importing 
PubChem compounds anyway, _biobuild_ provides a function `pubchem_to_cif` that can be used to convert PubChem compound
data into a CIF file which can be more easily edited and subsequently loaded into _biobuild_ using the `PDBECompounds` class
or using `Molecule.from_cif()`.

CHARMM
======

CHARMM is a molecular simulation software that is widely used in the field of computational chemistry.
It is developed and maintained by the `CHARMM Development Project <http://www.charmm.org/>`_.
_biobuild_ integrates the CHARMM force field for use in conformational optimization.

Sequons
=======

Sequons are short amino acid sequences that are recognized by the cell's glycosylation machinery.
They are used to identify potential glycosylation sites in proteins. _biobuild_ integrates a small database of sequons
that can be used to identify potential glycosylation sites in proteins.
"""
import biobuild.resources.pdbe_compounds as pdbe_compounds
import biobuild.resources.pubchem as pubchem
import biobuild.resources.charmm as charmm

from biobuild.resources.pdbe_compounds import *
from biobuild.resources.charmm import *
