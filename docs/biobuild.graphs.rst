The graphs package
==================

This module defines classes to represent PDB structures as graph objects.
These store the molecule connectivity information. Two graphs are available, 
the `AtomGraph`, as an atomic-scale representation of a molecule, and the `ResidueGraph` which represents each residue as a single node in the graph, thereby abstracting away the atomic details.

The graphs contain many of the analytical methods used for molecule structure analysis and manipulation that does not concern adding or removing atoms.

The graphs also serve an important function in biobuild's conformational optimization process as they are the main data structure to which the optimization algorithms are applied.


The BaseGraph module
--------------------

The `BaseGraph` is at the basis of both the `AtomGraph` and `ResidueGraph`.

.. dropdown:: The `BaseGraph` class

   .. autoclass:: buildamol.graphs.BaseGraph.BaseGraph
      :members:
      :undoc-members:
      :show-inheritance:



The AtomGraph module
--------------------

The `AtomGraph` handles the atom connectivity within `Molecule` objects.
It provides the bulk of connectivity related methods such as `get_neighbors`.

.. dropdown:: The `AtomGraph` class

   .. autoclass:: buildamol.graphs.AtomGraph
      :members:
      :undoc-members:
      :show-inheritance:

The ResidueGraph module
-----------------------

The `ResidueGraph` handles the residue connectivity within `Molecule` objects.
It is an abstraction of the `AtomGraph` and provides the many of the same methods.
`ResidueGraph` objects serve as primary input for structural optimization algorithms in the `optimizers`
package of buildamol. 

.. dropdown:: The `ResidueGraph` class

   .. autoclass:: buildamol.graphs.ResidueGraph
      :members:
      :undoc-members:
      :show-inheritance:
