
---

![](docs/_resources/logo_large.png)

---

[![cite - BMC Cheminformatics](https://img.shields.io/badge/cite-BMC_Cheminformatics-50EDEA)](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00900-6)
[![version DOI - 10.5281/zenodo.12581092](https://img.shields.io/badge/version_DOI-10.5281%2Fzenodo.12581092-blue)](https://doi.org/10.5281/zenodo.12581092)
[![Made with Python](https://img.shields.io/badge/Python->=3.8-blue?logo=python&logoColor=white)](https://python.org "Go to Python homepage")
[![Documentation Status](https://readthedocs.org/projects/biobuild/badge/?version=latest)](https://biobuild.readthedocs.io/en/latest/?badge=latest)
[![Check out - Tutorials](https://img.shields.io/badge/check_out-Tutorials-e61882)](https://biobuild.readthedocs.io/en/latest/tutorials.html)
[![PyPI version](https://badge.fury.io/py/buildamol.svg)](https://badge.fury.io/py/buildamol)
[![Downloads](https://static.pepy.tech/badge/buildamol)](https://pepy.tech/project/buildamol)
[![code style - black](https://img.shields.io/badge/code_style-black-black)](https://black.readthedocs.io/ "Go to Black homepage")
[![CodeFactor](https://www.codefactor.io/repository/github/noahhenrikkleinschmidt/buildamol/badge/main)](https://www.codefactor.io/repository/github/noahhenrikkleinschmidt/buildamol/overview/main)


BuildAMol is a molecular building suite designed to facilitate the generation and alteration of atomic models for large and small chemical structures.

It allows for an easy modeling process inside a Jupyter Notebook  or can be integrated into automated pipelines. BuildAMol offers direct integrations to [PubChem](https://pubchem.ncbi.nlm.nih.gov), and the [PDBE component library](https://www.google.com/search?client=safari&rls=en&q=pdbe+component+library&ie=UTF-8&oe=UTF-8) as well as the [CHARMM project](http://charmm-gui.org) to provide pre-defined template structures and linkages to use out-of-the-box. Quick-conversions to popular libraries such as [RDKit](https://www.rdkit.org) allow for a smooth workflow, going from modeling to analysis.

BuildAMol allows users to:
--------------------------
- build any larger molecular structure they like with full control
- automate molecular modeling tasks (e.g. see the [Ligand Design Pipeline](https://biobuild.readthedocs.io/en/latest/examples/ligand_design.html) or [Molecular Derivatives](https://biobuild.readthedocs.io/en/latest/examples/derivator_example.html) examples)
- improve the conformation of an existing structure
- create highly customizable 2D and 3D visualizations of structures
- quickly obtain molecular structures for chemical compounds
- convert data formats

BuildAMol cannot:
-----------------
- model real-life chemical reaction mechanisms in detail
- perform molecular dynamics or quantum chemistry computations
- generate molecules for the user **out of the blue** - the user needs to to have some idea of what to build or how to build it...


Installing BuildAMol
--------------------

BuildAMol can be installed via pip using:

```bash
pip install buildamol
```

Getting Started
---------------
BuildAMol has a comprehensive [documentation](https://biobuild.readthedocs.io/en/latest/) on ReadTheDocs. There you can find also also a number of **tutorials** to get you started on the API covering both basic operations as well as more complex and applied workflows such as building materials, preparing molecules for molecular dynamics, or designing protein ligands. 

Example 1 - Building A Dendrimer From Scratch
---------------------------------------------

This code will model a polyphenylene dendrimer as it was originally described by [Bauer et al. (2002)](https://doi.org/10.1002/1521-3765(20020902)8:17<3858::AID-CHEM3858>3.0.CO;2-5) using only the most barebone functionalities of the BuildAMol library with 100% control over which atoms form bonds and how.
 
```python
import buildamol as bam

bam.load_small_molecules()
benzene = bam.molecule("benzene")

# set up the linkage instructions
# always shifting the carbon at which to attach
periphery = benzene.copy()
link = bam.linkage(None, "C1")
for carbon in range(1, 6):
    link.atom1 = f"C{carbon}"
    periphery.attach(benzene, link, at_residue=1)

# assemble the dendrimer starting with the central benzene
mol = benzene.copy()
link2 = bam.linkage(None, "C4")

# and attach the periphery to the core
for carbon in mol.get_atoms("C", by="element"):
    link2.atom1 = carbon
    mol.attach(periphery, link2, at_residue=1, other_residue=2)

# optimize the conformation
mol.optimize()
mol.to_pdb("polyphenylene.pdb")
mol.show3d()
```

![](support/graphics/polyphenylene.gif)

Example 2 - Making a Glycan-Aspirin Conjugate
---------------------------------------------

There are also a number of already available extensions to make life easier when constructing certain kinds of molecules. For example, we can build Glycans directly from commonly used IUPAC notation. We can also exploit BuildAMol's various inference-level tools to determine how to connect molecular fragments together. In the example below we create a glycan-drug conjugate, automatically searching for the right atoms to use for connecting the molecules.

```python
import buildamol as bam
from buildamol.extensions.bio import glycans
from buildamol.structural import constraints_v2 as constraints
from buildamol.structural.reactivity import Carboxyl

# construct a small glycan
glycan = glycans.glycan("Neu5Ac(a2-3)Gal(b1-4)GlcNAc")

# and now create a conjugate with a drug-like molecule
# e.g. aspirin
aspirin = bam.molecule("aspirin")

# find the right atoms to define a linkage 
# (here: connect the Nitrogen atom of the last sugar residue
# to the carbonyl Carbon of the carboxyl group of aspirin, while splitting of an acetonic acid)
N = glycan.get_atom("N", by="element", residue=-1)
C_next_to_N = N.get_neighbors(filter=constraints.has_double_bond_with("O")).pop()

# find the reactive carboxyl atoms in the aspirin with help of the Carboxyl reactivity class
carboxyl = Carboxyl()
C_of_COOH, OH_of_COOH = carboxyl.find_atoms(aspirin, role="electrophile", serves_target=False)

# define the linkage between glycan and aspirin
link = bam.linkage(
    N, C_of_COOH,
    delete_in_target=[C_next_to_N],
    delete_in_source=[OH_of_COOH],
)

# now create the conjugate
conjugate = bam.connect(glycan, aspirin, link)
conjugate.draw2d().highlight_residues(-1, color="yellow").draw()
```
![](docs/_resources/glycan_drug_conjugate.png)

Example 3 - Letting Molecules "react" via Functional Groups
-----------------------------------------------------------

BuildAMol does not limit itself to "chemically" plausible connections of Molecules. However, it provides ways to conveniently connect Molecules based on reactivity patterns of functional groups. For example we can hydroxylate a carbon polymer and then let the hydroxyl groups react with the amide groups of a second molecule like so:

```python
import buildamol as bam
from buildamol.structural.reactivity import Hydroxyl, Amide
from buildamol.extensions.polymers.polycarbons import cyclic_alkane

# create two cyclic alkanes 
# (one larger that we hydroxylate at several positions,
# the other smaller that we amidate at one position)
A = cyclic_alkane(65)
sites = A.get_atoms(list(range(1, 65, 5)))
A = bam.hydroxylate(A, sites, [A.get_equatorial_hydrogen(s) for s in sites])

B = cyclic_alkane(6)
B = bam.amidate(B, 1, B.get_equatorial_hydrogen(1))

# now we let them react via the hydroxyl and amide groups
# attaching the amides and splitting off water
reaction = bam.Reaction.from_reactivities(
    nucleophile=Amide(),
    electrophile=Hydroxyl(),
)

out = reaction(A, B)
out.draw2d().draw(width=800, height=800)
```

![](docs/_resources/cyclic_alkanes.png)

Example 4 - Multi-Molecule Systems
----------------------------------

BuildAMol allows allows the free arrangement and orientation of molecules in 3D space. We can for instance create a dual ring system using the cyclic alkane from Example 3.

```python
# out from Example 3
dual_rings = out.copy()
# copy the ring, move it and rotate it
ring2 = out.copy().move([12, 0, 2]).rotate(90, "x")

# merge the two rings as one system
dual_rings.merge(ring2)
dual_rings.show3d()
```

![](docs/_resources/cyclic_alkanes_dual.GIF)


BuildAMol Paper
---------------
To learn more about the benchmarking we did and further details on the software, please check out the [BuildAMol paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00900-6). If you were using BuildAMol for your project, cite the paper as follows:

```
@article{buildamol,
	author = {Kleinschmidt, Noah and Lemmin, Thomas},
	journal = {Journal of Cheminformatics},
	number = {1},
	pages = {104},
	title = {BuildAMol: a versatile Python toolkit for fragment-based molecular design},
	volume = {16},
	year = {2024}}
```