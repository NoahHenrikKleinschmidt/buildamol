"""
The base class for the Patcher and Stitcher classes.
"""


class Connector:
    def __init__(self, copy_target: bool = False, copy_source: bool = False):
        self.copy_target = copy_target
        self.copy_source = copy_source

        self.target = None
        self.source = None

        self._anchors = (None, None)
        self._target_residue = None
        self._source_residue = None

    def set_target(self, target: "Molecule.Molecule"):
        """
        Set the target molecule

        Parameters
        ----------
        target : Molecule
            The target molecule
        """
        self.target = target

    def set_source(self, source: "Molecule.Molecule"):
        """
        Set the source molecule

        Parameters
        ----------
        source : Molecule
            The source molecule
        """
        self.source = source

    def get_anchors(self, _ref_atoms, target_residue=None, source_residue=None):
        """
        Find appropriate anchor atoms in the source and target molecules

        Parameters
        ----------
        _ref_atoms : tuple
            The two atoms that form the bond between the two molecules.
            These may be specified by name or by index.
        target_residue : int or str or Residue
            The residue in the target molecule to which the source molecule will be patched.
            By default, the last residue in the target molecule will be used if it contains
            appropriate anchor atoms.

        source_residue : int or str or Residue
            The residue in the source molecule that will be patched into the target molecule.
            By default, the first residue in the source molecule will be used if it contains
            appropriate anchor atoms.

        Returns
        -------
        tuple
            The two atoms that form the bond between the two molecules
        """
        if not self.target:
            raise AttributeError("No target set")
        if not self.source:
            raise AttributeError("No source set")

        if not _ref_atoms[0]:
            ref_atom_1 = [self.target.root_atom]
        else:
            # the whole point of this stuff is to allow the stitcher (which uses this)
            # to also work with patches (which would nromally have a 1/2 prefix) which would
            # otherwise prevent anchor finding...
            if isinstance(_ref_atoms[0], str) and _ref_atoms[0].startswith("1"):
                ref_atom_1 = self.target.get_atoms(_ref_atoms[0][1:])
            else:
                ref_atom_1 = self.target.get_atoms(_ref_atoms[0])

        if not _ref_atoms[1]:
            ref_atom_2 = [self.source.root_atom]
        else:
            if isinstance(_ref_atoms[1], str) and _ref_atoms[1].startswith("2"):
                ref_atom_2 = self.source.get_atoms(_ref_atoms[1][1:])
            else:
                ref_atom_2 = self.source.get_atoms(_ref_atoms[1])

        if target_residue:
            target_residue = self.target.get_residue(target_residue)
            ref_atom_1 = [i for i in ref_atom_1 if i.parent == target_residue]

        if source_residue:
            source_residue = self.source.get_residue(source_residue)
            ref_atom_2 = [i for i in ref_atom_2 if i.parent == source_residue]

        if len(ref_atom_1) == 0:
            raise ValueError("No anchor atom found in target molecule")
        if len(ref_atom_2) == 0:
            raise ValueError("No anchor atom found in source molecule")

        ref_atom_1 = ref_atom_1[-1]
        ref_atom_2 = ref_atom_2[0]

        self._anchors = (ref_atom_1, ref_atom_2)
        self._target_residue = ref_atom_1.parent
        self._source_residue = ref_atom_2.parent
        return ref_atom_1, ref_atom_2
