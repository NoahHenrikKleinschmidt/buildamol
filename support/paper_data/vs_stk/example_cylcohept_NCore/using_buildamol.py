from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import buildamol as bam

t_start = time()

# get the core molecule
core = bam.molecule("cycloheptane")

# now add methyl arms to the core (we make sure to add them equatorially)
for carbon in core.get_atoms("C", by="element"):
    H_to_delete = core.get_equatorial_hydrogen(carbon)
    bam.methylate(core, carbon, H_to_delete)

# write the core to a pdb file
core.to_pdb(parent / "bam_core.pdb")

print(f"Time taken: {time() - t_start:.2f} s")
