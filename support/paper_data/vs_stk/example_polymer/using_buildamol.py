from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent.parent

import buildamol as bam

t_start = time()

# define the repeating unit
unit = bam.molecule("N=CCCC=NCC").autolabel()

# specify the link between the units (using atom ids)
link_units = bam.linkage("C6", "N4")

# build the polymer
polymer = unit.repeat(10, link_units)

# make the end groups
# (change the nitrogen at the start to an oxygen to make an aldehyde)
polymer.change_element(polymer.get_atom("N4", residue=1), "O")
# (amidate the other end)
bam.amidate(polymer, at_atom=polymer.get_atom("C6", residue=-1))

# write the polymer to a pdb file
polymer.to_pdb(parent / "bam_polymer.pdb")

print(f"{time() - t_start:.2f}")
