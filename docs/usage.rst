.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Overview
========

BuildAMol is a python library to create, modify, and visualize molecular structures in a simple and intuitive way.
Complex structures can be built from simple building blocks, and the library is designed to be user-friendly, extendible, and highly flexible.
BuildAMol is not designed for any specifc type or class of molecule. Whether it's small protein ligands, metal complexes, or large polymers, BuildAMol can handle it all.
BuildAMol is designed to be an expert-driven library, meaning that it is designed to be used by people who are already familiar with molecular structures and their properties.
However, BuildAMol can be used in automated frameworks to generate molecular structures without human intervention (check out our tutorial on automated ligand generation to see how this might be done).
Also, to faciliate the workflow of the user, BuildAMol provides easy interfaces with other popular libraries such as `Bio.PDB`, `RDKit`, and `OpenBabel` to allow for direct interconversion
between these libraries and BuildAMol.

Here are a few things that BuildAMol can do:
- assemble two benzene rings together to form a biphenyl
- polymerize a glucose molecule to form cellulose
- change the elements in a carboxylic acid to form an amide
- create a metal complex with three ligands
- optimize a molecular structure
- search for substructures to identify functional groups
- visualize a protein-ligand complex
- compute the bond angles in a molecule
- sample 10 different conformations of a molecule
- and much more!

Most of the functionality of BuildAMol is based on the `Molecule` class, which is equipped with a large number of methods to manipulate the structure.
However, BuildAMol also provides a top-level functional API as well as an operator-based API to allow for more concise and intuitive code depending on the user's preferences!
What is more, BuildAMol comes with a built-in database of small molecules, amino acids, sugars, and lipids from the PDB Component Dictionary, which offers thousands of fragments to be 
usable off-line directly in BuildAMol. However, if you need more molecules, BuildAMol will retrieve any unknown molecule automatically from PubChem so you have access to 50 million potential fragment molecules without leaving your Jupyter Notebook or Python script!

Quicklook Examples
==================

Building a simple polymer
-------------------------

Let's make a simple polymer with a repeating unit composed of a benzene ring and two amino acids, glycine and tyrosine.

.. code-block:: python

   import buildamol as bam

   # first we load some reference data
   # (not necessary but makes getting the molecules 
   # quicker, since we don't need to query PubChem)
   bam.load_small_molecules()
   bam.load_amino_acids()  

   # now get the molecules we want to use
   benzene = bam.molecule("benzene")
   glycine = bam.molecule("glycine")
   tyrosine = bam.molecule("tyrosine")


   # now we can connect the molecules together
   # we start by defining how we want to connect the benzene
   # ring to the amino acid tyrosine
   link = bam.linkage("C1", "N") # = connect C1 from benzene to N from tyrosine
   mol = benzene.attach(tyrosine, link1)

   # now we can add a glycine to the molecule
   link.atom1 = "C3" # just re-use the link object but change the target atom
   mol = mol.attach(glycine, link, at_residue=1) # at_residue=1 means we always attach to the first residue (=benzene)

   # now we can polymerize the molecule
   
   # we make a new linkage (luckily there are only two atoms with names C5 and OH so we don't have
   # to worry about specifying which residues to connect exactly)
   # let's use a for-loop to automate the process
   link2 = bam.linkage("C5", "OH")
   _mol = mol.copy()
   for i in range(10):
      mol = mol.attach(_mol, link2)
   
   # fun fact: we could have just used
   # mol = mol.repeat(10, link2)

   # now that we have the molecule we can perform some quick
   # optimization on it
   mol = mol.optimize()

   # save the molecule to a MOL file
   mol.to_molfile("polymer.mol")

   # and finally visualize it
   mol.show()

.. image:: _resources/simple_polymer.gif
   :width: 80%
   :align: center
   :alt: Example repeated structure with benzene ring, glycine and tyrosine.


Building a glycan
-----------------

Glycans are complex carbohydrates that are often found in biological systems such as in glycoproteins.
Since they are structurally highly diverse and can be quite large, they are a good example to showcase the capabilities of BuildAMol.
Due to their flexibility, glycans can be tricky to predict and model with fully automated deep-learning methods but they are easy to build manually using the monoscaccharides:

Actually, BuildAMol not only includes the monsaccharides in its database, but also the most common glycosidic linkages, so building a glycan is as simple as connecting the monosaccharides together!

.. code-block:: python

   import buildamol as bam

   bam.load_sugars()

   # get the monosaccharides
   # (using their PDBE identifiers)
   nag = bam.molecule("NAG") # N-acetylglucosamine, a.k.a. GlcNAc
   bma = bam.molecule("BMA") # beta-mannose
   man = bam.molecule("MAN") # alpha-mannose

   # start by connecting two NAGs together
   # the % operator specifies which linkage to use when
   # connecting the two molecules using the + operator
   # btw. the 'beta 1->4' glycosydic linkage is pre-defined
   # in the database as '14bb' and can be accessed using its ID
   glycan = nag % "14bb" + nag

   # add a beta-mannose to the last NAG
   # (where + makes a copy, += will modify the molecule in place)
   glycan += bma

   # add an alpha-mannose to the beta-mannose
   # using an 'alpha 1->3' linkage ('13ab' in the database)
   glycan.attach(man, "13ab")

   # add another alpha-mannose
   # at the second-to-last residue (BMA)
   glycan.attach(man, "16ab", at_residue=-2)

   # add one final alpha-mannose
   # (this time using the connect function instead of the attach method or the + operator)
   glycan = bam.connect(glycan, man, "16ab")

   # save the glycan to a PDB file
   glycan.to_pdb("glycan.pdb")

   # now visualise the structure
   glycan.show()

.. image:: _resources/glycan_example.gif
   :width: 80%
   :align: center
   :alt: Example glycan structure.

In the above visualization, `NAG` residues are colored in pink, `BMA` in orange, and `MAN` in green. Hetero-atoms are colored according to IUPAC conventions.

The above example demonstrates how we can use BuildAMol to create a glycan structure from scratch. The example also demonstrates how we can use the three different syntaxes
to achieve this. Using the toplevel function `connect`, using the method `attach`, or by simple "molecular arithmetics" through the `+` operator.

