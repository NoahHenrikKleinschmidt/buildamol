.. biobuild documentation master file, created by
   sphinx-quickstart on Tue Jun 13 14:40:03 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Using `biobuild`
================

`biobuild` offers essentially three possible syntaxes:
- a toplevel functional API
- a method-based API
- an operator-based API

Functional API
--------------

The toplevel functional API is the most straightforward one, and is for many operations the easiest one.
Functions such as `connect` or `polymerize` are intuitive and easy to use. However, the functional API is 
not as extensive as the method-based API, since it simply calls the method-based API under the hood.

Method-based API
----------------

The method-based API allows for more flexibility. At the heart of `biobuild` is the `Molecule` class, which
has a great number of methods to modify the structure or add/remove parts of it. Methods such as `get_atom`,
`add_bond`, or `remove_residue` are basic operations to modify a structure in-place, while methods such as `attach`
allow to expand the structure by adding new parts to it.

Operator-based API
------------------

The operator-based API is a short-hand proxy of the method-based API (just as the functional API is a proxy).
It is essentially restricted to operations that regard connecting two molecules together. However, it is the most
condensed way to write `biobuild` code - sometimes at the expense of readability. Available operators are:

- `+` for connecting two molecules together
- `*` for polymerizing a molecule
- `%` for specifying the linkage between two molecules
- `@` for specifying the residue at which to create a connection between two molecules
- `^` for specifying the atom to use for a connection (more detailed than `@`)

.. note::

   In-place versions of the operators are also available, e.g. `+=` for connecting two molecules in-place, or `*=` for in-place polymerization.


Built-in-resources
==================

`biobuild` comes with a number of built-in resources. Namely, `biobuild` integrates the `PDBE component library <https://www.ebi.ac.uk/pdbe/pdb-component-library/#:~:text=The%20PDB%20Component%20Library%20is,and%20related%20protein%20structural%20data.>` for
components up to 40 atoms in size by default - naturally, the full library can be loaded if desired. Furthermore, `biobuild` integrates parts of the `CHARMM force field <https://www.charmm.org/>`_ for
references of molecular connections. Finally, `biobuild` integrates `pubchempy` for the direct retrieval of molecules from PubChem.

Toplevel functions exist to access these resources, e.g. `biobuild.available_linkages()` to get a list of pre-defined linkages,
or `biobuild.has_compound("alpha-mannose")` to check if a particular compound is available in the loaded PDBE component library.

Example
=======

`biobuild` was originally conceptualized with the aim of creating glycan structures - so, please, forgive if the example below produces a glycan. The following example demonstrates
how we can create a larger structure from single monosaccharides using `biobuild` (using all three syntaxes intermixed):

.. code-block:: python

   import biobuild as bb

   # get the monosaccharides
   # (using their PDBE identifiers)
   nag = bb.molecule("NAG") # N-acetylglucosamine, a.k.a. GlcNAc
   bma = bb.molecule("BMA") # beta-mannose
   man = bb.molecule("MAN") # alpha-mannose

   # start by connecting two NAGs together
   # 'beta 1->4' glycosydic linkage is pre-defined
   # in the CHARMM force field and can be used by its name '14bb' directly
   glycan = nag % "14bb" + nag

   # add a beta-mannose to the last NAG
   glycan += bma

   # add an alpha-mannose to the beta-mannose
   # using an 'alpha 1->3' linkage ('13ab' in CHARMM)
   glycan.attach(man, "13ab")

   # add another alpha-mannose
   # at the second-to-last residue (BMA)
   glycan.attach(man, "16ab", at_residue=-2)

   # add one final alpha-mannose
   glycan = bb.connect(glycan, man, "16ab")

   # now visualise the structure
   glycan.show()

