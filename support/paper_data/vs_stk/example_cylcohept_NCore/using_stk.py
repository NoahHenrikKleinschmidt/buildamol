import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import stk

t_start = time.time()

arm = stk.BuildingBlock(
    smiles="BrC",
    functional_groups=stk.BromoFactory(),
)
core2 = stk.BuildingBlock(
    smiles="C(Br)1C(Br)C(Br)C(Br)C(Br)C(Br)C1Br",
    functional_groups=stk.BromoFactory(),
)

ncore = stk.ConstructedMolecule(
    stk.small.NCore(
        core_building_block=core2,
        arm_building_blocks=arm,
        repeating_unit="A",
        optimizer=stk.MCHammer(),
    )
)

stk.PdbWriter().write(ncore, parent / "stk_ncore.pdb")

print(f"Time taken: {time.time() - t_start:.2f} s")
