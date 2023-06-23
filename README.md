
![](docs/_resources/logo_large.png)

`biobuild` is a molecular building suite designed to facilitate the creation of large biomolecules such as glycans. 
It allows for an easy molecule creation process in a jupyter-notebook environment. `biobuild` offers direct integrations
to _pubchem_, and the _PDBE component library_ as well as the _CHARMM force field_ for pre-defined component structures and linkage types.

`biobuild` allows users to:
---------------------------
- build any larger molecular structure they like
- improve the conformation of an existing structure
- convert data formats
- visualize the structures as they build them
- quickly obtain molecular structures for chemical compounds

`biobuild` cannot:
------------------
- generate circular structures (users need to choose suitable templates with rings already present)
- imitate real-life chemical reaction mechanisms
- perform molecular dynamics or quantum chemistry computations
- generate molecules _for_ the user - the user needs to know what they want to build...


Example
-------

Generating a glycan structure is as simple as:

```python

import biobuild as bb

# get the monosaccharides
# (using their PDBE identifiers)
nag = bb.molecule("NAG") # N-acetylglucosamine, a.k.a. GlcNAc
bma = bb.molecule("BMA") # beta-mannose
man = bb.molecule("MAN") # alpha-mannose

# start by connecting two NAGs together
# 'beta 1->4' glycosydic linkage is pre-defined
# in the CHARMM force field and can be used by its name '14bb' directly
glycan = bb.connect(nag, nag,  "14bb")

# add a beta-mannose to the last NAG
glycan.attach(bma, "14bb")

# add an alpha-mannose to the beta-mannose
# using an 'alpha 1->3' linkage ('13ab' in CHARMM)
glycan.attach(man, "13ab")

# add another alpha-mannose
# at the second-to-last residue (BMA)
glycan.attach(man, "16ab", at_residue=-2)

# add one final alpha-mannose
glycan = bb.connect(glycan, man, "16ab")

# now visualise the structure
glycan.show()

# and save to a PDB file
glycan.to_pdb("my_glycan.pdb")
```

![](docs/_build/html/_images/glycan_example.gif)
