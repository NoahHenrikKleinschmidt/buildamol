.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Why `biobuild`?
---------------

`biobuild` is a mocule building suite for cheminformatics and bioinformatics.
It is designed to be a simple, easy to use, and extensible tool for building
larger molecular structures from smaller components. `biobuild` tries to leverage the power
of geometric operations in combination with machine learning algorithms 
to facilitate the building of complex molecular structures where direct imputation from SMILES may not be feasible.

In order to facilitate the building of larger structures, `biobuild` provides a set of tools to not only easily
and quickly connect smaller molecules to larger structures, but also integrates a vast database of available molecules
that can be loaded without the need to download external data sources, or read files from disk. 

Why "bio"-build and not "molecule"-build, or something?
-------------------------------------------------------

`biobuild` is designed on top of `biopython <https://biopython.org/>`_ and was originally designed to build larger
biomolecules such as glycans. However, the underlying algorithms and implementations are not limited to these fields
and hence, `biobuild` was created as a stand-along framework for molecule building. `biopython` was chosen as backend
due to its great flexibility when it comes to adding and removing atoms to and from structures. However, `biobuild`
is not limited to `biopython` but exports import and export of structures to and from `rdkit <https://www.rdkit.org/>`_,
another widely used cheminformatics library.

What can `biobuild` do?
-----------------------

`biobuild` can be used to build larger molecules from smaller components - that's pretty much it. However, given its integrated
compounds database, it can also be used to quickly obtain molecular structures for other purposes such as docking. Its ability
to convert to and from `rdkit` also makes it a great tool to edit `rdkit` structures in a more intuitive way. Lastly, `biobuild`
offers a suite of conformational optimization algorithms to improve the geometry of structures. Hence, `biobuild` can be used
to sample conformers for molecular docking or the generation of a deep learning dataset. 

Who is `biobuild` for?
----------------------

Anyone who would like to build larger molecules from smaller components. `biobuild` is designed to be easy to use and intuitive. 
Most of its user-relevant functions are accessible through toplevel functions, or methods of the `Molecule` class. That in mind,
`biobuild` is designed to work interactively with the user, for instance in a Jupyter notebook - `biobuild` offers some customizable 3D visualization. 
However, `biobuild` is equally suited for scripting and integration into larger workflows such as `snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipelines.
How it is used is up to you!

