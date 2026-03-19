"""
Basic constants and stuff for the tests
"""

from pathlib import Path


HOME = Path(__file__).parent.parent

FILES = Path(__file__).parent / "files"

ALLOW_VISUAL = False


MANPDB = FILES / "man.pdb"
GLCXML = FILES / "glc.xml"
GLCPDB = FILES / "glc.pdb"
GALPDB = FILES / "gal.pdb"
MAN9PDB = FILES / "man9.pdb"
PDBE_TEST_FILE = FILES / "compounds.cif"
CHARMM_TOPOLOGY_FILE = FILES / "patches.rtf"
DMBXML = FILES / "dimethylbenzene.xml"
EX8PDB = FILES / "ex8.pdb"
EX6JSON = FILES / "ex6.json"
ROTSCAFJSON = FILES / "rotaxane_scaffold.json"
PADPDB = FILES / "paddlewheel.pdb"
