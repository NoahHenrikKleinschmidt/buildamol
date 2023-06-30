The structural package
======================

The `structural` package contains modules for dealing with molecular structures - it is the functional core that performs the real action behind most of the methods that are available in the `Molecule` class's methods. 
This includes operations such as computing angles, rotating around bonds, and assembling molecules together.

The `structural` package is divided into several submodules, each of which is described below.

Basic structural operations
---------------------------

The `base` module contains the basic structural operations that are used by the other modules.

.. automodule:: biobuild.structural.base
   :members:
   :undoc-members:
   :show-inheritance:

Inferring structural properties
-------------------------------

The `infer` module contains methods for inferring structural properties from a molecule's atoms and connectivity.
Among these are connections, or atom labels.

.. automodule:: biobuild.structural.infer
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: biobuild.structural.smiles
   :members:
   :undoc-members:
   :show-inheritance:


Assembling molecules
--------------------

In order to assemble two molecules together there are two modules that are used: `stitch` and `patch`.
The `patch` module defines a `patch` function and corresponding `Patcher` class that can be used to assemble two molecules together using
a `patch` (i.e. a `Linkage` with internal coordinates). The `stitch` module defines a `stitch` function and corresponding `Stitcher` class that can be used to assemble two molecules together using a `recipe` (i.e. a `Linkage` without internal coordinates).


.. tab-set::

   .. tab-item:: patch module

      .. automodule:: biobuild.structural.patch
         :members:
         :undoc-members:
         :show-inheritance:
   
   .. tab-item:: stitch module

      .. automodule:: biobuild.structural.stitch
         :members:
         :undoc-members:
         :show-inheritance:

   
   .. tab-item:: connector module

      .. automodule:: biobuild.structural.connector
         :members:
         :undoc-members:
         :show-inheritance:


Graph neighborhood in molecular structures
------------------------------------------

In order to infer the neighborhoods from connectivity information in molecules, the `neighbors` module defines a `Neighborhood` class that can be used to infer the neighborhoods of atoms in a molecule.


.. automodule:: biobuild.structural.neighbors
   :members:
   :undoc-members:
   :show-inheritance:

