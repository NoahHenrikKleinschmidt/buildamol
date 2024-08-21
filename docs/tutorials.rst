.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Tutorials
=========

On this page you will find some tutorials on how to use BuildAMol. Feel free to browse through them to learn more.

Basic Usage
-----------

.. nbgallery::
   examples/getting_started
   examples/building_workflow
   examples/drawing
   examples/syntax_flavours
   examples/optimization
   examples/converting_formats

Resource-oriented Tutorials
---------------------------

.. nbgallery::
   examples/built-in
   examples/default_settings


BuildAmol and Chemistry
-----------------------

.. nbgallery::
   examples/functional_groups
   examples/geometry_scratch
   examples/derivator_example


BuildAmol and Biology
---------------------

.. nbgallery::
   examples/bio_extension
   examples/glycosylation
   examples/ligand_design

BuildAmol and Materials
-----------------------

.. nbgallery::
   examples/nanotube
   examples/metal_organic_framework

BuildAmol and Molecular Dynamics
--------------------------------

.. nbgallery::
   examples/molecular_dynamics
   examples/conformation_sampling
   examples/solvation_box

   
More applied tutorials
-----------------------

.. nbgallery::

   examples/making_modifiers
   examples/rotaxane_easy
   examples/rotaxane_hard
   examples/building_circular
   examples/building_polyphenylene
   examples/irregular_dendrimer
   examples/fragment_library

.. .. _cards-clickable::

.. .. card:: Getting Started
..    :link: example_getting_started
..    :link-type: ref

..    This tutorial covers the fundamentals of buildamol. It is slightly more detailed than the "Basic Usage" page, however, and covers topics that the other tutorials elaborate on.


.. .. card:: A typical workflow
..    :link: example_building_workflow
..    :link-type: ref

..    In this tutorial we will build a molecule to get the hang of a "typical" BuildAMol workflow.


.. .. card:: Syntax Flavors
..    :link: example_syntax
..    :link-type: ref

..    This tutorial covers the different syntax flavors that BuildAMol supports.


.. .. card:: Visualizing Molecules
..    :link: example_visualizing_molecules
..    :link-type: ref

..    This tutorial covers how to visualize molecules using buildamol.


.. .. card:: Optimizing Molecules
..    :link: example_optimization
..    :link-type: ref

..    This tutorial covers how to optimize molecular conformations in buildamol.


.. .. card:: Other Libraries and Formats
..    :link: example_converting_formats
..    :link-type: ref

..    This tutorial covers how to convert between different file formats and how to use other libraries with buildamol.


.. Resource-oriented Tutorials
.. ---------------------------

.. .. card:: Built-in Resources
..    :link: example_built_in_resources
..    :link-type: ref

..    This tutorial covers how to use the built-in resources of buildamol.

.. .. card:: Adding Resources
..    :link: example_defaults
..    :link-type: ref

..    This tutorial covers how to add your own resources to BuildAMol so they are available to all your projects.


.. More applied Tutorials
.. ----------------------

.. .. card:: Using Functional Groups
..    :link: example_functional_groups
..    :link-type: ref

..    In this tutorial we will explore how to use functional groups to connect molecules together in BuildAMol.
..    We also cover how to create your own custom functional groups.

.. .. card:: Building circular Molecules
..    :link: example_building_circular
..    :link-type: ref

..    This tutorial covers how to build circular molecules in BuildAMol together with RDKit.

.. .. card:: Building Polyphenylene
..    :link: example_building_polyphenylene
..    :link-type: ref

..    In this tutorial we will build a polyphenylene dendrimer using some automization.

.. .. card:: Glycosylation
..    :link: example_glycosylation
..    :link-type: ref

..    In this tutorial we will glycosylate a protein.


.. .. card:: Molecular Dynamics
..    :link: example_md
..    :link-type: ref

..    In this tutorial we will build a molecule and run a molecular dynamics simulation on it.

.. .. card:: Conformation Sampling
..    :link: example_conformation_sampling
..    :link-type: ref

..    In this tutorial we will generate multiple conformers of a molecule using BuildAMol's optimization methods.

.. .. card:: Building a Rotaxane - The easy way
..    :link: example_rotaxane_easy
..    :link-type: ref

..    In this tutorial we will build a rotaxane using BuildAMol and learn how to combine multiple molecules into one system by aligning and merging them.

.. .. card:: Building a Rotaxane - The hard way
..    :link: example_rotaxane_hard
..    :link-type: ref

..    In this tutorial we will build the same rotaxane as above but create our own optimization setup to spacially arrange the ring around the axle molecule.


.. .. card:: Building a Solvation Box
..    :link: example_solvationbox
..    :link-type: ref

..    In this tutorial we will build a solvation box around a molecule using BuildAMol.

.. .. card:: Molecules from scratch - building PF5
..    :link: example_geometry
..    :link-type: ref

..    In this tutorial we will build a PF5 molecule using BuildAMol's molecular geometries to automatically generate coordinates.
   
.. .. card:: Molecules from scratch - building a Nanotube
..    :link: example_nanotube
..    :link-type: ref

..    In this tutorial we will build a nanotube 100% from scratch, atom by atom, using numpy to build coordinates and BuildAMol to assemble a molecule.
