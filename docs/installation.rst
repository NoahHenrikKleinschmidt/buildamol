.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
------------

.. note:: 

   Biobuild is not yet released officially on PyPI. You can currently install it from the github repository.
   A release is coming shortly, so stay tuned!


.. `Biobuild` can be directly installed via `pip` from the python package index:

.. .. code-block:: bash

..    pip install biobuild

If you want to install Biobuild from source, you can clone the repository and install it using the following command:

.. code-block:: bash

   git clone https://github.com/NoahHenrikKleinschmidt/biobuild
   cd biobuild
   pip install .

Optional dependencies
---------------------

Biobuild relies only on `biopython <https://biopython.org/>`_ for its core functions. However, many additional features may require other libraries.
Especially, `RDKit <https://www.rdkit.org/>`_ provides many useful features that Biobuild can use if it is installed. Apart from RDKit,
Biobuild supports data exchange between `Openbabel <http://openbabel.org/wiki/Main_Page>`_ and `OpenMM <http://openmm.org/>`_. If you want to install 
these libraries you can use the following commands:

.. code-block:: bash

   pip install biobuild[rdkit]
   pip install biobuild[openbabel]
   pip install biobuild[openmm]
   # or all at once
   pip install biobuild[full]

We recommend to install at least `biobuild[rdkit]` if you want to use most features of `Biobuild`.


Updating `biobuild`
-------------------

.. warning:: 

   When updating `biobuild` to a newer version, make sure to export any custom default `PDBECompounds` or `CHARMMTopology` settings
   before updating, as the update will overwrite the pickle files where defaults are stored. You can re-set the custom settings after the update.

In order to export the custom data it is recommended not to pickle the objects (as a pickled object may not be compatbile with a newer version of Biobuild).
Instead, you can export the data as a JSON file using the `to_json` method of the `PDBECompounds` and `CHARMMTopology` classes or the `export_custom_resources` function:

.. code-block:: python

   import biobuild as bb

   bb.export_custom_resources('my_settings')
   # after updating biobuild
   bb.import_custom_resources('my_settings')

