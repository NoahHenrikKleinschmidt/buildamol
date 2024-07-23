The core package
================

.. automodule:: buildamol.core
   :members:
   :undoc-members:
   :show-inheritance:

The `Molecule` module
---------------------

The `Molecule` is the heart of the buildamol package. It is a class that
handles molecular structures and allows the user to easily assemble them into 
larger constructs. The `Molecule` is a child of the `BaseEntity` class that provides
most of the functionality. The `Molecule` module defines additionally a number of toplevel
functions to easily create new molecules by querying databases or reading files.

.. automodule:: buildamol.core.Molecule
   :members:
   :undoc-members:
   :show-inheritance:

The `Linkage` module
--------------------

The `Linkage` module defines the `Linkage` class that is used to connect two molecules
in a specific way. 

.. automodule:: buildamol.core.Linkage
   :members:
   :undoc-members:
   :show-inheritance:

The `base` module
-----------------

The `entity` module defines the `BaseEntity` class that is the base class for `Molecules`
(and whatever other classes a user may wish to define that are similar in concept).

.. dropdown:: The `BaseEntity` class

   .. autoclass:: buildamol.core.entity.BaseEntity
      :members:
      :undoc-members:
      :show-inheritance:

In addition the `base_classes` module defines wrappers around native `BioPython` classes
such as `Atom`, `Residue`, etc. These classes are used by buildamol in order to facilitate atom identifcation
in situations where multiple identical molecules are connected to each other. All these classes support a `from_biopython` and `to_biopython` conversion.

.. dropdown:: The `base_classes` module

   .. automodule:: buildamol.base_classes
      :members:
      :undoc-members:
      :show-inheritance:

      

.. toctree::
   :maxdepth: 2

   buildamol.core.Molecule
   buildamol.core.Linkage
   buildamol.core.base

