from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import stk
import stko

t1 = time()

cycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="[Br]CC[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="A",
        num_repeating_units=5,
        optimizer=stk.MCHammer(),
    ),
)
axle = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()]),
        ),
        repeating_unit="AB",
        num_repeating_units=7,
        optimizer=stk.MCHammer(),
    )
)
rotaxane = stk.ConstructedMolecule(
    topology_graph=stk.rotaxane.NRotaxane(
        axle=stk.BuildingBlock.init_from_molecule(axle),
        cycles=(stk.BuildingBlock.init_from_molecule(cycle),),
        repeating_unit="A",
        num_repeating_units=3,
    ),
)
# Write the constructed molecule to a pdb file.
rotaxane.write(parent / "stk_rotaxane.pdb")

print(time() - t1)

# Optimize the rotaxane.
rotaxane = stko.UFF().optimize(rotaxane)

# Write the optimized molecule to a pdb file.
rotaxane.write(parent / "stko_rotaxane.pdb")

print(time() - t1)
