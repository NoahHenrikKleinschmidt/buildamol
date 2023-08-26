.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. image:: _resources/biobuild_docs_header.png
   :width: 80%
   :align: center
   :alt: biobuild

===================================
Welcome to Biobuild's documentation
===================================

`Biobuild` is a fragment-based molecular assembly toolkit for the generation of atomic models for medium to large chemical compounds.
It is designed to help researchers working on structrual biology generate realistic structures of molecules
that may not be immediatelly imputable from SMILES alone and for which no 3D references are immediately downloadable
from databases. `Biobuild` is based on `biopython <http://biopython.org/wiki/Main_Page>`_ and accessible as a `python package`
offering a simple API to generate, manipulate, visualize, and export 3D structures of molecules.

.. admonition:: Example
      
   .. image:: _resources/large.gif
      :width: 80%
      :align: center
      :alt: dendrimer

   This dendrimer was chemically described by `Pedro-Hernandez et al. (2022) <http://benthamscience.com/article/119156>`_ and generated with `Biobuild`.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :hidden:

   whatfor
   installation
   usage
   tutorials
   documentation
