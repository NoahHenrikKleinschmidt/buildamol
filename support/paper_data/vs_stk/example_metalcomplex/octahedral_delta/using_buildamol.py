from time import time
from pathlib import Path

parent = Path(__file__).resolve().parent

import buildamol as bam
from buildamol.extensions import complexes

t_start = time()

iron = bam.Atom.new("FE", pqr_charge=2)
complexer = complexes.MetalComplexer(iron, bam.structural.geometry.Octahedral())
ligand = bam.Molecule.from_smiles(
    "C1=CC=[N+](C(=C1)C2=CC(=CC=C2)C3=[N+](C=CC=C3)[H])[H]", id="LIG"
)

complexer.make_core(4)

complex = complexer.add_ligands(
    ligands=[ligand] * 3,
    binders=[[i.id for i in ligand.get_atoms("N", by="element")]] * 3,
    acceptors=[("H1", "H5"), ("H6", "H4"), ("H2", "H3")],
)

complex.to_pdb(parent / "bam_complex.pdb")

print(f"Time taken: {time() - t_start:.2f} s")

complex.show()
