"""
This is literally the stk example code
"""

from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import stk
import stko

t_start = time()

metal = stk.BuildingBlock(
    smiles="[Fe+2]",
    functional_groups=(stk.SingleAtom(stk.Fe(0, charge=2)) for i in range(6)),
    position_matrix=[[0, 0, 0]],
)

ligand = stk.BuildingBlock(
    smiles="C1=CC=NC(=C1)C2=CC(=CC=C2)C3=NC=CC=C3",
    functional_groups=[
        stk.SmartsFunctionalGroupFactory(
            smarts="[#6]~[#7X2]~[#6]",
            bonders=(1,),
            deleters=(),
        ),
    ],
)

complex = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.OctahedralDelta(
        metals=metal,
        ligands=ligand,
        optimizer=stk.MCHammer(),
    ),
)

stk.PdbWriter().write(complex, parent / "stk_complex.pdb")

print(f"{time() - t_start:.2f}")

# optimize the complex
complex = stko.UFF().optimize(complex)

stk.PdbWriter().write(complex, parent / "stko_complex.pdb")

print(f"{time() - t_start}")
