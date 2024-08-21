from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import buildamol as bam

t_start = time()

# get the core molecule
core = bam.molecule("cycloheptane").autolabel()

# now add methyl arms to the core (we make sure to add them equatorially)
for carbon in core.get_atoms("C", by="element"):
    H_to_delete = core.get_equatorial_hydrogen(carbon)
    bam.methylate(core, carbon, H_to_delete)

    # now amidate some methyl groups
    # the incoming methyl group will always be the last residue...
    if carbon.id in ("C1", "C3", "C5"):
        metyhl_group = core.get_residue(-1)
        carbon = core.get_atom("C", by="element", residue=metyhl_group)
        bam.aminate(core, carbon)

# write the core to a pdb file
core.to_pdb(parent / "bam_core.pdb")

print(f"{time() - t_start:.2f}")
