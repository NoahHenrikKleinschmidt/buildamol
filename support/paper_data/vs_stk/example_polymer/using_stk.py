"""
This is literally the stk example code
"""

from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent.parent

import stk
import stko

t_start = time()

# React the amine functional groups during construction.
bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
# React the aldehyde functional groups during construction.
bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
# Build a polymer.
polymer = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(bb1, bb2),
        repeating_unit="AB",
        num_repeating_units=10,
        optimizer=stk.Collapser(scale_steps=False),
    ),
)

stk.PdbWriter().write(polymer, parent / "stk_polymer.pdb")

print(f"{time() - t_start:.2f}")

# Optimize the polymer.
polymer = stko.UFF().optimize(polymer)

stk.PdbWriter().write(polymer, parent / "stko_polymer.pdb")

print(f"{time() - t_start}")
