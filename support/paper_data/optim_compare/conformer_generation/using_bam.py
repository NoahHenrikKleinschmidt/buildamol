"""
The script used to generate conformers for the testing molecule using BuildAMol
"""

import buildamol as bam
from pathlib import Path

DIR = Path(__file__).parent.parent
FILE = DIR / "test_mol.pdb"
DESTINATION = DIR / "bam_ov_conformers"
if not DESTINATION.exists():
    DESTINATION.mkdir()

mol = bam.molecule(str(FILE))
mol.add_hydrogens()

N = 50

# Generate N conformers using totally default settings
for i in range(N):
    mol_copy = mol.copy()
    mol_copy.optimize(rotatron="overlap", rotatron_kws=dict(spread=5))
    mol_copy.to_pdb(DESTINATION / f"conf-{i}.pdb")

print("Done")
