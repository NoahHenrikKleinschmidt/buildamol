
![](docs/_resources/logo_large.png)

---

[![Documentation Status](https://readthedocs.org/projects/biobuild/badge/?version=latest)](https://buildamol.readthedocs.io/en/latest/?badge=latest)
[![CodeFactor](https://www.codefactor.io/repository/github/noahhenrikkleinschmidt/buildamol/badge/main)](https://www.codefactor.io/repository/github/noahhenrikkleinschmidt/buildamol/overview/main)


BuildAMol is a molecular building suite designed to facilitate the generation and alteration of atomic models for large and small chemical structures.

It allows for an easy modeling process inside a Jupyter Notebook  or can be integrated into automated pipelines. BuildAMol offers direct integrations to [PubChem](https://pubchem.ncbi.nlm.nih.gov), and the [PDBE component library](https://www.google.com/search?client=safari&rls=en&q=pdbe+component+library&ie=UTF-8&oe=UTF-8) as well as the [CHARMM project](http://charmm-gui.org) to provide pre-defined template structures and linkages to use out-of-the-box. Quick-conversions to popular libraries such as [RDKit](https://www.rdkit.org) allow for a smooth workflow, going from modeling to analysis.

BuildAMol allows users to:
--------------------------
- build any larger molecular structure they like
- improve the conformation of an existing structure
- convert data formats
- visualize the structures as they build them
- quickly obtain molecular structures for chemical compounds

BuildAMol cannot:
-----------------
- imitate real-life chemical reaction mechanisms
- perform molecular dynamics or quantum chemistry computations
- generate molecules _for_ the user out of the blue - the user needs to know what they want to build...


Installing BuildAMol
--------------------

BuildAMol can be installed via pip using:

```bash
pip install buildamol
```

Getting Started
---------------
BuildAMol has a comprehensive [documentation](https://biobuild.readthedocs.io/en/latest/) on ReadTheDocs. There you can find also also a number of **tutorials** to get you started on the API covering both basic operations as well as more complex and applied workflows such as building materials, preparing molecules for molecular dynamics, or designing protein ligands. 


Quick Example - Building A Dendrimer
------------------------------------

This code will model a polyphenylene dendrimer as it was originally described by [Bauer et al. (2002)](https://doi.org/10.1002/1521-3765(20020902)8:17<3858::AID-CHEM3858>3.0.CO;2-5). 
```python
import buildamol as bam

bam.load_small_molecules()

benzene = bam.molecule("benzene")

# -----------------------------
#     make the periphery
# -----------------------------
periphery = benzene.copy()

# set up the linkage instructions
# always shifting the carbon at which to attach
link = bam.linkage("C1", "C1")
for carbon in range(1, 6):
    link.atom1 = f"C{carbon}"
    periphery.attach(benzene, link, at_residue=1)

# -----------------------------
#     assemble the molecule
# -----------------------------
mol = benzene.copy()
link2 = bam.linkage("C1", "C4")

# and attach the periphery to the core
for carbon in mol.get_atoms("C", by="element"):
    link2.atom1 = carbon
    mol.attach(periphery, link2, at_residue=1, other_residue=2)

# -----------------------------
#   optimize the conformation
# -----------------------------
mol.optimize()
mol.to_pdb("polyphenylene.pdb")
```

![](support/graphics/polyphenylene.gif)
