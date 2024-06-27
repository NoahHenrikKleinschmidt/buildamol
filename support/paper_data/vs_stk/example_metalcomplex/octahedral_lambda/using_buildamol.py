from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import buildamol as bam
from buildamol.extensions import complexes

t_start = time()

iron = bam.Atom.new("FE", pqr_charge=2)
complexer = complexes.MetalComplexer(iron, bam.structural.geometry.Octahedral())
ligand = bam.Molecule.from_smiles("C=[N+](CC=[N+](Br)[H])[H]", id="LIG").autolabel()

complexer.make_core(2)
complex = complexer.add_ligands(
    ligands=[ligand] * 3,
    binders=[("N1", "N3") for i in range(3)],
    acceptors=[("H1", "H5"), ("H6", "H4"), ("H2", "H3")],
    optimize_kwargs={
        "algorithm": "scipy",
    },
)

complex.to_pdb(parent / "bam_complex.pdb")

print(f"Time taken: {time() - t_start:.2f} s")

complex.show()
