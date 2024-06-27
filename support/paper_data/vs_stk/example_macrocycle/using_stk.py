from time import time
from pathlib import Path

import stk
import stko

parent = Path(__file__).resolve().parent

t_start = time()

macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="BrCCBr",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="BrCNCBr",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="AB",
        num_repeating_units=5,
        optimizer=stk.MCHammer(),
    ),
)

stk.PdbWriter().write(macrocycle, parent / "stk_macrocycle.pdb")

print(f"{time() - t_start:.2f} s")

# optimize the macrocycle
macrocycle = stko.UFF().optimize(macrocycle)

stk.PdbWriter().write(macrocycle, parent / "stko_macrocycle.pdb")

print(f"{time() - t_start:.2f}")
