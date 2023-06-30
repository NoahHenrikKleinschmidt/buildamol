The core package
================

.. automodule:: biobuild.core
   :members:
   :undoc-members:
   :show-inheritance:

The `Molecule` module
---------------------

The `Molecule` is the heart of the biobuild package. It is a class that
handles molecular structures and allows the user to easily assemble them into 
larger constructs. The `Molecule` is a child of the `BaseEntity` class that provides
most of the functionality. The `Molecule` module defines additionally a number of toplevel
functions to easily create new molecules by querying databases or reading files.

.. automodule:: biobuild.core.Molecule
   :members:
   :undoc-members:
   :show-inheritance:

The `Linkage` module
--------------------

The `Linkage` module defines the `Linkage` class that is used to connect two molecules
in a specific way. 

.. automodule:: biobuild.core.Linkage
   :members:
   :undoc-members:
   :show-inheritance:

The `base` module
-----------------

The `entity` module defines the `BaseEntity` class that is the base class for `Molecules`
(and whatever other classes a user may wish to define that are similar in concept).

.. dropdown:: The `BaseEntity` class

   .. autoclass:: biobuild.core.entity.BaseEntity
      :members:
      :undoc-members:
      :show-inheritance:

In addition the `base_classes` module defines biobuil'd wrappers around native `BioPython` classes
such as `Atom`, `Residue`, etc. These classes are used by biobuild in order to facilitate atom identifcation
in situations where multiple identical molecules are connected to each other. All these classes support a `from_biopython` and `to_biopython` conversion.

.. dropdown:: The `base_classes` module

   .. automodule:: biobuild.core.base_classes
      :members:
      :undoc-members:
      :show-inheritance:

      

.. toctree::
   :maxdepth: 2

   biobuild.core.Molecule
   biobuild.core.Linkage
   biobuild.core.base

