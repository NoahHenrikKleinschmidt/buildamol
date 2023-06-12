"""
This module defines classes to represent PDB structures as graph objects.
These store the molecule connectivity information. Two graphs are available, 
the `AtomGraph`, as an atomic-scale representation of a molecule, and the `ResidueGraph` which represents each residue as a single node in the graph, thereby abstracting away the atomic details.

The graphs contain most of the analytical methods used for molecule structure analysis and manipulation that does not concern adding or removing atoms.

The graphs also serve an important function in _biobuild's_ conformational optimization process as they are the main data structure to which the optimization algorithms are applied.
"""

from biobuild.graphs.AtomGraph import AtomGraph
from biobuild.graphs.ResidueGraph import ResidueGraph

__all__ = ["AtomGraph", "ResidueGraph"]
