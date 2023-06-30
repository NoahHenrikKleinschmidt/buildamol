.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Installation
------------

`biobuild` can be directly installed via `pip` from the python package index:

.. code-block:: bash

   pip install biobuild

If you want to install `biobuild` from source, you can clone the repository and install it using the following command:

.. code-block:: bash

   git clone https://github.com/NoahHenrikKleinschmidt/biobuild
   cd biobuild
   pip install .


Updating `biobuild`
-------------------

.. warning:: 

   When updating `biobuild` to a newer version, make sure to export any custom default `PDBECompounds` or `CHARMMTopology` settings
   before updating, as the update will overwrite the pickle files where defaults are stored. You can re-set the custom settings after the update.

In order to export the custom data it is recommended not to pickle the objects (as a pickled object may not be compatbile with a newer version of `biobuild`).
Instead, you can export the data as a JSON file using the `to_json` method of the `PDBECompounds` and `CHARMMTopology` classes:

.. code-block:: python

   my_pdbecompounds = PDBECompounds(...)
   my_charmm_topology = CHARMMTopology(...)

   my_pdbecompounds.to_json('my_pdbecompounds.json')
   my_charmm_topology.to_json('my_charmm_topology.json')

After updating `biobuild` you can re-set the custom data using the `from_json` method of the `PDBECompounds` and `CHARMMTopology` classes:

.. code-block:: python

   my_pdbecompounds = PDBECompounds.from_json('my_pdbecompounds.json')
   my_charmm_topology = CHARMMTopology.from_json('my_charmm_topology.json')

For more convenience, biobuild also offers the `export_custom_resources` and `import_custom_resources` functions to perform both actions at once:

.. code-block:: python

   import biobuild as bb

   bb.export_custom_resources('my_settings')
   # after updating biobuild
   bb.import_custom_resources('my_settings')

