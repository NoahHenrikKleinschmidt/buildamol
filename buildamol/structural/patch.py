"""
Functions to patch molecules together
"""

import numpy as np

# import buildamol.utils.abstract as abstract
import buildamol.structural.connector as base
import buildamol.structural.infer as infer
import buildamol.core.Molecule as Molecule
import buildamol.core.Linkage as Linkage


class PatchError(Exception):
    pass


class Patcher(base.Connector):
    """
    This class is responsible for patching molecules together
    based on a `Linkage` object and two `Molecule` objects,
    of which one is the target and one is the source. Residues from the source
    are integrated into the target molecule.
    The Linkage object must define a "patch", meaning
    it must contain data on internal coodrinates (ICs) of the resulting molecule in the vicinity
    of the newly formed bond.

    Parameters
    ----------
    copy_target : bool
        Whether to copy the target molecule before patching
    copy_source : bool
        Whether to copy the source molecule before patching
    """

    def __init__(self, copy_target: bool = False, copy_source: bool = False):
        super().__init__(copy_target=copy_target, copy_source=copy_source)
        self.patch = None

    def set_patch(self, patch: "Linkage.Linkage"):
        """
        Set the patch

        Parameters
        ----------
        patch : Linkage
            The patch
        """
        self.patch = patch

    def apply(
        self,
        patch: "Linkage.Linkage" = None,
        target: "Molecule.Molecule" = None,
        source: "Molecule.Molecule" = None,
        target_residue=None,
        source_residue=None,
    ):
        """
        Patch two molecules together.
        This will merge the source molecule's residues into the target molecule.
        If no patch, target, and source have been set already, they can be provided
        as arguments to this method.

        Parameters
        ----------
        patch : Linkage
            The patch to apply (this needs to have ICs)
        target : Molecule
            The target molecule
        source : Molecule
            The source molecule
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
            The patched target and source molecules (not yet merged into one molecule)
        """
        if patch:
            self.patch = patch
        if target:
            self.target = target
        if source:
            self.source = source

        if self.copy_target:
            self.target = self.target.copy()
        if self.copy_source:
            self.source = self.source.copy()

        # a dictionary to store the anchor atoms in the source
        # molecule and their original coordinates
        self._source_computed_anchors = {}

        if not target_residue:
            target_residue = self.target.attach_residue
        if not source_residue:
            source_residue = self.source.attach_residue
        self.get_anchors(target_residue, source_residue)

        # compute internal coordinates to rebuild the source molecule
        # _ics = infer.compute_internal_coordinates(self.source.bonds)
        # self._source_ICs = abstract.AbstractEntity()
        # self._source_ICs.internal_coordinates = _ics

        # self._v.draw_edges(self.source.bonds, color="red", opacity=0.5)

        # align the source molecule to the target molecule
        self._align_anchors()

        # now start imputing the source molecule into the target molecule
        # starting with the atoms that are described in the patch's internal coordinates
        self._infer_patch_IC_neighbors()

        # and now trans-rotate the source molecule such that
        # the inferred neighboring atoms overlay their imputed counterparts
        self._transpose_source()

        self._delete_atoms()

        return self.target, self.source

    def merge(self):
        """
        Merge the source molecule into the target molecule

        Returns
        -------
        Molecule
            The target molecule with the source molecule merged into it
        """
        # self.target.adjust_indexing(self.source)
        self.target.add_residues(*self.source.residues)
        self.target._bonds.extend(self.source.bonds)
        self.target._AtomGraph.migrate_bonds(self.source._AtomGraph)
        self.target.add_bond(*self._anchors)
        return self.target

    def get_anchors(self, target_residue=None, source_residue=None):
        """
        Returns the two atoms that will be used for anchoring structure alignment.
        These atoms form a bond between the two molecules.

        Parameters
        ----------
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
            the first atom is from molecule1 and the second is from molecule2.
        """
        if not self.patch:
            raise AttributeError("No patch set")

        _ref_atoms = {int(i[0]) - 1: i[1:] for i in self.patch._ref_atoms}
        ref_atom_1, ref_atom_2 = super().get_anchors(
            _ref_atoms, target_residue, source_residue
        )
        return ref_atom_1, ref_atom_2

    def _align_anchors(self):
        """
        Align the source to the target molecule such that the anchors are in the correct distance
        """
        ic_34 = self._match_IC_34()
        if len(ic_34) == 0:
            raise PatchError("No matching IC found")
        ic = ic_34[0]

        atom1 = self._get_ref_atom(ic.atom1)
        atom2 = self._get_ref_atom(ic.atom2)

        # atom1 = self.target.get_atoms(ic.atom1[1:])
        # atom2 = self.target.get_atoms(ic.atom2[1:])
        # if self._target_residue:
        #     atom1 = [i for i in atom1 if i.get_parent() == self._target_residue]
        #     atom2 = [i for i in atom2 if i.get_parent() == self._target_residue]
        # if len(atom1) == 0:
        #     raise PatchError("No anchor atom found in target molecule")
        # if len(atom2) == 0:
        #     raise PatchError("No anchor atom found in target molecule")
        # atom1 = atom1[-1]
        # atom2 = atom2[-1]

        new_anchor_source = infer.compute_atom4_from_others(
            atom1.coord,
            atom2.coord,
            self._anchors[0].coord,
            ic,
        )

        # self._v.draw_point("target ref1", atom1.coord, color="teal")
        # self._v.draw_point("target ref2", atom2.coord, color="teal")

        # store the old coordinates and set the new coordinates
        self._source_computed_anchors[self._anchors[1]] = self._anchors[1].coord
        self._anchors[1].set_coord(new_anchor_source)

        # self._v.draw_point("anchor_source (new)", self._anchors[1].coord, color="green")
        # self._v.draw_point("anchor_target", self._anchors[0].coord, color="green")

    def _infer_patch_IC_neighbors(self):
        """
        Infer the neighbors of the internal coordinates in the patch.
        This is required to rebuild the source molecule.
        """
        for ic in self.patch.internal_coordinates:
            if ic.atom4[0] == "1":
                continue

            atom4 = self._get_ref_atom(ic.atom4)
            if atom4 in self._source_computed_anchors.keys():
                continue

            atom3 = self._get_ref_atom(ic.atom3)
            atom2 = self._get_ref_atom(ic.atom2)
            atom1 = self._get_ref_atom(ic.atom1)

            _new_coords = infer.compute_atom4_from_others(
                atom1.coord, atom2.coord, atom3.coord, ic
            )

            # self._v.draw_point(
            #     atom4.id + " (old)",
            #     atom4.coord,
            #     color="brown",
            #     opacity=0.6,
            # )

            if np.isnan(_new_coords).any():
                continue

            self._source_computed_anchors[atom4] = atom4.coord
            atom4.set_coord(_new_coords)

            # self._v.draw_point(
            #     atom4.id
            #     + " (new should-be)"
            #     + f"tb-deleted={atom4.id in self.patch.deletes[1]}",
            #     atom4.coord,
            #     color="blue",
            #     opacity=0.6,
            # )

    def _transpose_source(self):
        """
        Transpose the source molecule to overlay the
        """

        if len(self._source_computed_anchors) < 3:
            raise PatchError("Not enough anchors to transpose")

        _old_coords = np.stack(list(self._source_computed_anchors.values()))
        _new_coords = np.array([i.coord for i in self._source_computed_anchors.keys()])

        # for n in _old_coords:
        #     self._v.draw_point("old", n, color="red")

        # for n in _new_coords:
        #     self._v.draw_point("new", n, color="limegreen")

        # compute translation vector
        old_centroid = _old_coords.mean(axis=0)
        new_centroid = _new_coords.mean(axis=0)
        # translation_vector = new_centroid - old_centroid

        # self._v.draw_vector(
        #     "translation vector",
        #     old_centroid,
        #     translation_vector,
        #     color="teal",
        # )

        _relative_old_coords = _old_coords - old_centroid
        _relative_new_coords = _new_coords - new_centroid

        H = (_relative_old_coords).T.dot(_relative_new_coords)
        U, S, VT = np.linalg.svd(H)
        R = VT.T @ U.T

        atoms = [
            i
            for i in self.source.get_atoms()
            if i not in self._anchors and i not in self._source_computed_anchors.keys()
        ]
        atom_coords = np.array([i.coord for i in atoms])
        atom_coords -= old_centroid
        atom_coords = atom_coords @ R.T
        atom_coords += new_centroid  # old_centroid + translation_vector

        for atom, coord in zip(atoms, atom_coords):
            atom.set_coord(coord)

        # THE OLD STUFF HERE WORKED JUST FINE BUT IS LESS EFFICIENT SINCE WE
        # KEEP PERFORMING THE MATRIX MULTIPLICATIONS MANY TIMES IN A PYTHON
        # LOOP. THE NEW STUFF ABOVE IS MORE EFFICIENT SINCE WE ONLY PERFORM
        # THE MATRIX MULTIPLICATIONS ONCE AND THEN APPLY THE TRANSFORMATION
        # TO ALL ATOMS AT ONCE.

        # for atom in self.source.get_atoms():
        # self._v.draw_point(
        #     atom.id + " (old)",
        #     atom.coord,
        #     color="brown",
        #     opacity=0.3,
        #     showlegend=False,
        # )

        # if atom in self._source_computed_anchors.keys():
        #     continue

        # vec = self._source_computed_anchors[atom] - old_centroid
        # new_coord = (R @ vec) + old_centroid + translation_vector

        # self._v.draw_point(
        #     atom.id + " (new computed from old)",
        #     new_coord,
        #     color="purple",
        #     opacity=0.6,
        # )

        # atom.set_coord(current[idx])
        # idx += 1

        # vec = atom.coord - old_centroid
        # new_coord = (R @ vec) + old_centroid + translation_vector
        # atom.set_coord(new_coord)

        # self._v.draw_point(
        #     atom.id + " (new)",
        #     atom.coord,
        #     color="lightblue",
        #     opacity=0.6,
        # )

    # self._v.draw_edges(self.source.bonds, color="blue", opacity=0.6)

    def _match_IC(self, n_target: int, n_source: int):
        """
        Get appropriate internal coordinates from the patch
        that include the anchor atoms of the target and source molecule
        at specific positions (starting at 1).

        Parameters
        ----------
        n_target : int
            The position of the target anchor atom in the internal coordinate
        n_source : int
            The position of the source anchor atom in the internal coordinate

        Returns
        -------
        list
            All available internal coordinates that include the anchor atoms at the specified positions
        """

        ids = [None, None, None, None]
        ids[n_target - 1] = "1" + self._anchors[0].id
        ids[n_source - 1] = "2" + self._anchors[1].id

        ics = self.patch.get_internal_coordinates(
            *ids,
            mode="partial",
        )

        if len(ics) == 0:
            raise PatchError(
                "No internal coordinates found for anchor atoms at positions {} and {}".format(
                    n_target, n_source
                )
            )

        return ics

    def _match_IC_34(self, ic=None):
        """
        Get appropriate internal coordinates from the patch
        that include the anchor atoms of the target and source molecule
        at positions 3 and 4 and have three atoms from the target molecule.
        Or check if a given internal coordinate matches the criteria.

        Returns
        -------
        list or bool
            All available internal coordinates that include the anchor atoms at the specified positions
            or True if the given internal coordinate matches the criteria.
        """
        if ic:
            return ic.atom1[0] == "1" and ic.atom2[0] == "1" and ic.atom3[0] == "1"

        ids = ["1" + i.id for i in self.target.get_atoms()]
        ics = self._match_IC(3, 4)
        ics = [i for i in ics if i.atom1 in ids and i.atom2 in ids and i.atom3 in ids]
        if len(ics) == 0:
            raise PatchError(
                "No internal coordinates found for anchor atoms at positions 3 and 4"
            )
        return ics

    def _delete_atoms(self):
        """
        Delete all atoms that need to be deleted according to the patch
        """
        delete_from_target, delete_from_source = self.patch.deletes
        if self._target_residue:
            delete_from_target = [
                i for i in self._target_residue.child_list if i.id in delete_from_target
            ]
        if self._source_residue:
            delete_from_source = [
                i for i in self._source_residue.child_list if i.id in delete_from_source
            ]
        self.target._remove_atoms(*delete_from_target)
        self.source._remove_atoms(*delete_from_source)

    @property
    def _objs(self):
        return {
            # target, object, specific residue, index policy if no residue
            "1": (self.target, self._target_residue, -1),
            "2": (self.source, self._source_residue, 0),
        }

    def _get_ref_atom(self, atom: str):
        """
        Get the atom that is used as a reference for the given atom
        when computing the relative coordinates.

        Parameters
        ----------
        atom : str
            The atom for which to get the reference atom

        Returns
        -------
        Atom
            The reference atom
        """
        res, id = atom[0], atom[1:]
        _obj, _res, _idx = self._objs[res]

        atoms = _obj.get_atoms(id, by="id")
        if _res:
            atoms = [i for i in atoms if i.get_parent() is _res]
        if len(atoms) == 0:
            raise PatchError("No atom found with id {}".format(atom))
        atom = atoms[_idx]
        return atom


__default_keep_keep_patcher__ = Patcher(copy_target=False, copy_source=False)
"""
Default instance of the Patcher that will modify the target molecule in place
"""

# __default_keep_copy_patcher__ = Patcher(copy_target=False, copy_source=True)
# """
# Default instance of the Patcher that will modify the target molecule in place
# but use a copy of the source molecule to leave it intact.
# """

# __default_copy_copy_patcher__ = Patcher(copy_target=True, copy_source=True)
# """
# Default instance of the Patcher that will copy both the target and source molecules
# to leave the originals intact.
# """


def patch(
    target: "Molecule.Molecule",
    source: "Molecule.Molecule",
    patch: "Linkage.Linkage",
    target_residue=None,
    source_residue=None,
    copy_target: bool = False,
    copy_source: bool = False,
) -> "Molecule.Molecule":
    """
    Patch two molecules together.
    This will merge the source molecule's residues into the target molecule.

    Parameters
    ----------
    patch : Linkage
        The patch to apply
    target : Molecule
        The target molecule
    source : Molecule
        The source molecule
    target_residue : int or str or Residue
        The residue in the target molecule to which the source molecule will be patched.
        By default, the last residue in the target molecule will be used if it contains
        appropriate anchor atoms.
    source_residue : int or str or Residue
        The residue in the source molecule that will be patched into the target molecule.
        By default, the first residue in the source molecule will be used if it contains
        appropriate anchor atoms.
    copy_target : bool
        Whether to copy the target molecule before patching it
    copy_source : bool
        Whether to copy the source molecule before patching it

    Returns
    -------
    Molecule
        The patched molecule
    """
    if copy_target:
        target = target.copy()
    if copy_source:
        source = source.copy()
    __default_keep_keep_patcher__.apply(
        patch=patch,
        target=target,
        source=source,
        target_residue=target_residue,
        source_residue=source_residue,
    )
    return __default_keep_keep_patcher__.merge()


if __name__ == "__main__":
    import buildamol as bam

    bam.load_sugars()

    glc = bam.get_compound("GLC")

    patcher = Patcher(copy_target=False, copy_source=False)
    v = glc.draw()
    patcher._v = v

    link = bam.get_linkage("14bb")
    a, b = patcher.apply(link, glc, glc.copy())
    merged = patcher.merge()
    merged.show()
    # ser = bam.molecule("ser.json")
    # his = bam.molecule("his.json")

    # link = bam.Linkage.from_json("peptide_linkage.json")

    # patcher = Patcher(copy_target=False, copy_source=False)
    # v = ser.draw()
    # patcher._v = v
    # patcher.apply(link, ser, ser.copy())
    # merged = patcher.merge()

    # patcher.apply(link, merged, ser.copy())
    # merged = patcher.merge()

    # merged.show()

    # man = "support/examples/MAN.pdb"
    # man1 = bam.Molecule.from_pdb(man)
    # man1.infer_bonds()

    # man2 = man1.copy()

    # # now make sure that man2 has some different coordinates
    # man2.rotate_around_bond(1, 2, 35)
    # r = np.random.rand(3) * 0.1
    # for i in man2.atoms:
    #     i.coord += r

    # man1.lock_all()
    # man2.lock_all()

    # top = bam.get_default_topology()

    # colors = ["red", "green", "blue", "orange", "purple", "pink", "black"]

    # patcher = Patcher(False, True)
    # for i in ("14bb", "14bb", "12ab"):
    #     # v2 = bam.utils.visual.MoleculeViewer3D(man1)
    #     # patcher._v = v2

    #     patch_or_recipe = top.get_patch(i)
    #     patcher.apply(patch_or_recipe, man1, man2)
    #     man1 = patcher.merge()

    #     new = man1
    #     seen_atoms = set()
    #     for atom in new.atoms:
    #         assert atom.serial_number not in seen_atoms
    #         seen_atoms.add(atom.serial_number)

    #     res_con = bam.structural.infer_residue_connections(new.chain, triplet=True)

    #     v2 = bam.utils.visual.MoleculeViewer3D(man1)
    #     for idx, residue in enumerate(man1.residues):
    #         bonds = [
    #             i
    #             for i in man1.bonds
    #             if i[0] in residue.child_list and i[1] in residue.child_list
    #         ]
    #         v2.draw_edges(edges=bonds, color=colors[idx], linewidth=2)

    #     # v2.draw_edges(edges=list(new.bonds), color="blue", opacity=1)
    #     # v2.draw_edges(edges=list(new._locked_bonds), color="red", linewidth=3)
    #     # v2.draw_edges(edges=res_con, color="limegreen", linewidth=4)
    #     v2.show()
