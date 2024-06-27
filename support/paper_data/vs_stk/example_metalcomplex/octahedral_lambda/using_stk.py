"""
This is literally the stk example code
"""

from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import stk

t_start = time()

metal = stk.BuildingBlock(
    smiles="[Fe+2]",
    functional_groups=(stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(6)),
    position_matrix=[[0, 0, 0]],
)

bidentate = stk.BuildingBlock(
    smiles="C=NC/C=N/Br",
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#35]",
            bonders=(1,),
            deleters=(),
        ),
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    ],
)

complex = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.OctahedralLambda(
        metals=metal,
        ligands=bidentate,
        optimizer=stk.MCHammer(),
    ),
)

stk.PdbWriter().write(complex, parent / "stk_complex.pdb")

print(f"Time taken: {time() - t_start:.2f} s")
