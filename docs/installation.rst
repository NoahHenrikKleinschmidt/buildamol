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
   before updating, as the update will overwrite the default settings file. You can re-set the custom settings after the update.