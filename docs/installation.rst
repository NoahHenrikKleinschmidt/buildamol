.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
------------


`BuildAMol` can be directly installed via `pip` from the python package index:

.. code-block:: bash

   pip install buildamol

If you want to install BuildAMol from source, you can clone the repository and install it using the following command:

.. code-block:: bash

   git clone https://github.com/NoahHenrikKleinschmidt/biobuild
   cd biobuild
   pip install .

Optional dependencies
---------------------

BuildAMol relies only on `biopython <https://biopython.org/>`_ for its core functions. However, many additional features may require other libraries.
Especially, `RDKit <https://www.rdkit.org/>`_ provides many useful features that BuildAMol can use if it is installed. Apart from RDKit,
BuildAMol supports data exchange between `Openbabel <http://openbabel.org/wiki/Main_Page>`_ and `OpenMM <http://openmm.org/>`_. If you want to install 
these libraries you can use the following commands:

.. code-block:: bash

   pip install buildamol[rdkit]
   pip install buildamol[openbabel]
   pip install buildamol[openmm]
   # or all at once
   pip install buildamol[full]

We recommend to install at least `buildamol[rdkit]` if you want to use most features of `BuildAMol`.


Updating BuildAMol
------------------

.. warning:: 

   When updating BuildAMol to a newer version, make sure to export any custom default `PDBECompounds` or `CHARMMTopology` settings
   before updating, as the update will overwrite the pickle files where defaults are stored. You can re-set the custom settings after the update.

In order to export the custom data it is recommended not to pickle the objects (as a pickled object may not be compatbile with a newer version of BuildAMol).
Instead we recommend to export your data to an agnostic format before updating. BuildAMol offers quick-exports to JSON and XML formats for such purposes.
Both the `PDBECompounds` and `CHARMMTopology` classes have methods to export their data to JSON and XML formats. Alternatively, to export all custom resources at once, you can use the following functions:

.. code-block:: python

   import buildamol as bam

   bam.export_custom_resources('my_settings')
   
   # update ...

   # after updating buildamol
   bam.import_custom_resources('my_settings')

