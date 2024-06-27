from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import buildamol as bam
from buildamol.extensions import polymers

t_start = time()

# first make the cicle as an alkane then modify to get the nitrogens
N, M = 30, 5
cycle = polymers.cyclic_alkane(n=N)
for i in range(1, N, M):
    cycle.change_element(f"C{i+1}", "N")

# write the core to a pdb file
cycle.to_pdb(parent / "bam_cycle.pdb")

print(f"{time() - t_start:.2f}")
