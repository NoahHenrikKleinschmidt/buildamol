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
and quickly connect smaller molecules to larger structures, but also integrates a extensible database of available molecules
that can be loaded without the need to download external data sources, or read files from disk. 

Why "bio"-build?
----------------

`biobuild` is designed on top of `biopython <https://biopython.org/>`_ and was originally designed to build larger
biomolecules, specifically glycans. However, the underlying algorithms and implementations are not limited to these glycobiology
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

What can `biobuild` not do?
---------------------------

`biobuild` is not a quantum chemistry package. It does not offer any quantum chemical calculations, nor does it offer any tools
to analyze the electronic structure of molecules (aside from whatever is available from `biopython`, of course). 
It is also not a molecular dynamics package. Also, currently `biobuild` is limited to creating single bond-connections between
molecules. Hence, `biobuild` cannot create (small) ring structures. However, using fragments with the desired ring structures already present (e.g. from pubchem)
large circular molecules can be built regardless. However, `biobuild` is currently not optimized for circular molecules of any kind, hence, users wishing to 
create circular molecules will likely want to use molecular dynamics to optimize the geometry of any structure they create.
Finally, `biobuild` does not generate molecules `for` you! It is **not** a `de novo` molecule generator that tries to find molecules with certain properties - e.g. molecules that bind a certain protein pocket. 
The user needs to know what they want to build and have some idea of the steps involved - `biobuild` will help them build it, that's the deal.

Who is `biobuild` for?
----------------------

Anyone who would like to build larger molecules from smaller components. `biobuild` is designed to be easy to use and intuitive. 
Most of its user-relevant functions are accessible through toplevel functions, or methods of the `Molecule` class. That in mind,
`biobuild` is designed to work interactively with the user, for instance in a Jupyter notebook - `biobuild` offers some customizable 3D visualization. 
However, `biobuild` is equally suited for scripting and integration into larger workflows such as `snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipelines.
How it is used is up to you!

