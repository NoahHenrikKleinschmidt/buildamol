

Reaction
--------

The `Reaction` class defines a linkage-generator engine that can be used to automatically determine how two molecules should be connected based on a functional rule set.
The `Reaction` class works hand in hand with the :py:class:`Reactivity` classes that define how to find suitable attachment points in a molecule.
However, the `Reaction` class can also be used without a :py:class:`Reactivity` class by directly providing functions that define how to find attachment points in the source and target molecules.
Be sure to check out the Reaction tutorial and Reactivities tutorial for more information.

For example to let BuildAMol figure out how to connect a carboxylic acid to an amine to form an amide bond you can use the following code:

.. code-block:: python

      import buildamol as bam
      from buildamol.structural.reactivity import Carboxyl, Amine

      # start with two molecules
      acid = bam.read_smiles("CC(=O)O")
      amine = bam.read_smiles("CCN")

      # create a reaction that uses the Carboxyl reactivity to find attachment points in the acid
      # and the Amine reactivity to find attachment points in the amine
      reaction = bam.Reaction.from_reactivities(
         electrophile=Carboxyl(),
         nucleophile=Amine()
      )

      # now we can use the reaction to connect the two molecules
      product = reaction(acid, amine)

Make sure to check out the Tutorials for more examples on how to use the `Reaction` class and the various :py:class:`Reactivity` classes.


.. autoclass:: buildamol.core.Reaction
   :members:
   :undoc-members:
   :show-inheritance:
