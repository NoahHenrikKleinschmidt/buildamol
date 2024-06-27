import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import stk
import stko

t_start = time.time()

core = stk.BuildingBlock(
    smiles="C(Br)1C(Br)C(Br)C(Br)C(Br)C(Br)C1Br",
    functional_groups=stk.BromoFactory(),
)
arm1 = stk.BuildingBlock(
    smiles="BrC",
    functional_groups=stk.BromoFactory(),
)
arm2 = stk.BuildingBlock(
    smiles="BrCN",
    functional_groups=stk.BromoFactory(),
)

ncore = stk.ConstructedMolecule(
    stk.small.NCore(
        core_building_block=core,
        arm_building_blocks=[arm1, arm2],
        repeating_unit="ABABABA",
        optimizer=stk.MCHammer(),
    )
)

stk.PdbWriter().write(ncore, parent / "stk_ncore.pdb")

print(f"{time.time() - t_start:.2f}")

# optimize the ncore
ncore = stko.UFF().optimize(ncore)

stk.PdbWriter().write(ncore, parent / "stko_ncore.pdb")

print(f"{time.time() - t_start:.2f}")
