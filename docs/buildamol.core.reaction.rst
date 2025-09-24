

buildamol Reaction
------------------

The `Reaction` class defines a linkage-generator engine that can be used to automatically determine how two molecules should be connected based on a functional rule set.
The `Reaction` class works hand in hand with the `Reactivity` classes that define how to find suitable attachment points in a molecule.
However, the `Reaction` class can also be used without a `Reactivity` class by directly providing functions that define how to find attachment points in the source and target molecules.
Be sure to check out the Reaction tutorial and Reactivities tutorial for more information.

.. autoclass:: buildamol.core.Reaction
   :members:
   :undoc-members:
   :show-inheritance:
