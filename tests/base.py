"""
Basic constants and stuff for the tests
"""

import os


HOME = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

MANNOSE9 = os.path.join(HOME, "support/examples/man9.pdb")
"""
The PDB file for the example Mannose9 glycan
"""

MANNOSE = os.path.join(HOME, "support/examples/MAN.pdb")
"""
The PDB file for the Mannose sugar (only a single Mannose)
"""

GALACTOSE = os.path.join(HOME, "support/examples/GAL.pdb")
"""
The PDB file for the Galactose sugar (only a single Galactose)
"""

PROTEIN = os.path.join(HOME, "support/examples/4tvp.prot.pdb")
"""
THe PDB file of the protein 4TVP (no glycans)
"""

CHARMM_TOPOLOGY_FILE = os.path.join(HOME, "support/toppar_charmm/carbohydrates.rtf")
"""
The path to the CHARMM topology file
"""

CHARMM_PARAMETERS_FILE = os.path.join(HOME, "support/toppar_charmm/carbohydrates.prm")
"""
The path to the CHARMM parameters file
"""

PDBE_TEST_FILE = os.path.join(HOME, "support/pdbe_compounds/test.cif")
"""
The path to the test file for the PDBe compounds
"""
