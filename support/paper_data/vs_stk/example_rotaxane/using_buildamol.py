from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

from buildamol.extensions import polymers, complexes

t1 = time()

# make the ring
ring = polymers.cyclic_alkane(10)

# make the axle basis and then modify the atoms
N = 35
axle = polymers.linear_alkane(N)
axle.change_element(axle.get_hydrogen("C1"), "Br")
axle.change_element(axle.get_hydrogen(f"C{N}"), "Br")

i = 4
while i < N:
    # here we need to use change_element as the hydrogen count changes
    axle.change_element(f"C{i}", "N")
    i += 5

# make the rotaxane
rotaxane = complexes.rotaxane(axle, cycles=ring.copy(3), optimize=False)

# Write the constructed molecule to a pdb file.
rotaxane.to_pdb(parent / "bam_rotaxane.pdb")

print(time() - t1)
