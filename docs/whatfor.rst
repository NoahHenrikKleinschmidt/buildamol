.. BuildAMol documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Why BuildAMol?
--------------

BuildAMol is a mocule building suite for cheminformatics and bioinformatics.
It is designed to be a simple, easy to use, and extensible tool for building
larger molecular structures from smaller components. BuildAMol tries to leverage the power
of geometric operations in combination with machine learning algorithms 
to facilitate the building of complex molecular structures where direct imputation from SMILES may not be feasible.

In order to facilitate the building of larger structures, BuildAMol provides a set of tools to not only easily
and quickly connect smaller molecules to larger structures, but also integrates a extensible database of available molecules
that can be loaded without the need to download external data sources, or read files from disk.

Especially for molecules with repeated elements such as polymers, BuildAMol can be a handy tool allowing molecular assembly in a `for` loop
or similar construct to build large structures incrementally and thus much faster as compared to all-at-once approaches such as SMILES imputation.



What can BuildAMol do?
----------------------

BuildAMol can be used to build larger molecules from smaller components - that's pretty much it. It offers a number of ways to do so, however, ranging from 100% user-controlled atom-level connections to automated chemistry-inspired reaction templates. 

Also, given its integrated compounds database, it can also be used to quickly obtain molecular structures for other purposes such as docking. 
Molecules not available in the database are automatically fetched from `PubChem <https://pubchem.ncbi.nlm.nih.gov/>`_ providing 50 million compounds at your fingertips.
The ability to convert to and from `rdkit` on-the-fly also makes BuildAMol a great tool to edit `rdkit` structures in a more intuitive way. BuildAMol
offers a suite of conformational optimization algorithms to improve the geometry of structures and methods to manually adjust torsion angles and bond lengths among a myriad of other features, ensuring that *your* molecule has the geometry *you* want.
Finally, BuildAMol comes with a set of extensions that provide more applied functionalities for particular use cases such as building polymers, biomolecules, or complexes.


What can BuildAMol not do?
--------------------------

BuildAMol is not a quantum chemistry package. It does not offer any quantum chemical calculations, nor does it offer any tools
to analyze the electronic structure of molecules (aside from whatever is available from `biopython`, of course). 
It is also not a molecular dynamics package either. Circular structures are supported but not the primary focus of BuildAMol.
If you are looking to build circular molecules you will definitely want to install RDKit as well to optimize circular conformations.
Most importantly, BuildAMol is not a *de novo* molecular generator. It can be used to create automated molecule design pipelines but it needs some coding effort to do so. 


Who is BuildAMol for?
---------------------

Anyone who would like to build larger molecules from smaller components. BuildAMol is designed to be easy to use and intuitive. 
Most of its user-relevant functions are accessible through toplevel functions, or methods of the `Molecule` class. That in mind,
BuildAMol is designed to work interactively with the user, for instance in a Jupyter notebook - BuildAMol offers some customizable 3D visualization. 
However, BuildAMol is equally suited for scripting and integration into larger workflows such as `snakemake <https://snakemake.readthedocs.io/en/stable/>`_ pipelines.
How it is used is up to you!

