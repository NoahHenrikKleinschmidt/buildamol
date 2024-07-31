import buildamol.core as core
import buildamol.structural as structural
import buildamol.optimizers as optimizers

from typing import List

import numpy as np
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA

__all__ = ["rotaxane", "RotaxaneBuilder"]


def rotaxane(
    axle: core.Molecule,
    cycles: List[core.Molecule],
    optimize: bool = True,
    copy_axle: bool = False,
    copy_cycles: bool = False,
) -> core.Molecule:
    """
    Create a rotaxane from an axle and one or more cyclic molecules.

    Parameters
    ----------
    axle : core.Molecule
        The axle molecule.
    cycles : list[core.Molecule]
        One or more cyclic molecules. These will be threaded onto the axle.
    optimize: bool, optional
        Whether to optimize the spacial arrangement of each cycle around the axle.
        This is helpful if you have a non-linear axle molecule that requires a positional fitting
        to minimize the risk of clashes.
    copy_axle : bool, optional
        If True, the axle molecule will be copied, by default False
    copy_cycles : bool, optional
        If True, the cyclic molecules will be copied, by default False

    Returns
    -------
    core.Molecule
        The rotaxane where each cycle is its own chain.
    """
    if isinstance(cycles, core.Molecule):
        cycles = [cycles]
    elif not isinstance(cycles, (list, tuple)):
        raise ValueError(
            f"cycles must be a list or tuple of Molecules. Got {type(cycles)} instead..."
        )

    if copy_axle:
        axle = axle.copy()
    if copy_cycles:
        cycles = [cycle.copy() for cycle in cycles]

    builder = RotaxaneBuilder()
    builder.distribute_along_axle(axle, cycles)
    if optimize:
        builder.optimize()
    return builder.merge()


class RotaxaneBuilder:
    """
    A rotaxane builder class. Use this to retain more control over the rotaxane building process, compared to the toplevel `rotaxane` function.
    """

    def __init__(self):
        self.axle = None
        self.cycles = []

    def distribute_along_axle(self, axle: core.Molecule, cycles: List[core.Molecule]):
        """
        Distribute the cycles along the axle molecule.
        This will spacially position the cycles in equal intervals along the longest axis of the axle molecule and rotate them perpendicular to the longest axis.
        No molecules will be merged at this stage.

        Parameters
        ----------
        axle : core.Molecule
            The axle molecule
        cycles : List[core.Molecule]
            The cyclic molecules to distribute along the axle.
        """
        self.axle = axle
        if isinstance(cycles, core.Molecule):
            cycles = [cycles]
        elif not isinstance(cycles, (list, tuple)):
            raise ValueError(
                f"cycles must be a list or tuple of Molecules. Got {type(cycles)} instead..."
            )
        self.cycles = cycles

        self.axle_coords = self.axle.get_coords()
        self._set_initial_cycle_positions()
        self._rotate_cycles_perpendicular_to_main_axis()
        self.cycle_coords = np.array([cycle.get_coords() for cycle in self.cycles])

    def optimize(
        self,
        translation: bool = False,
        rotation: bool = True,
        translation_bounds: tuple = (-10, 10),
    ):
        """
        Optimize the spacial arrangement of each cycle around the axle.
        This is helpful if you have a non-linear axle molecule that requires a positional fitting
        to minimize the risk of clashes.

        This step requires that the cycles have been distributed along the axle using the `distribute_along_axle` method.

        Parameters
        ----------
        translation : bool, optional
            Whether to optimize the translation of the cycles, by default False
        rotation : bool, optional
            Whether to optimize the rotation of the cycles, by default True
        translation_bounds : tuple, optional
            The bounds for the minimal and maximal translation and rotation values.
            If a tuple of length two this is interpreted as the low and high bounds for translation only.
            Otherwise provide a tuple of length 6 for the bounds of translation and rotation. In this case values can be either
            singular (int/float) in which case they are interpreted as symmetric extrama (min=-value, max=+value) or as tuples with (min=value[0], max=value[1]).
            Mixed inputs are allowed.
        """
        if len(self.cycles) == 0:
            raise ValueError("No cycles are available for optimization!")

        bounds = [
            translation_bounds,
            translation_bounds,
            translation_bounds,
            np.pi,
            np.pi,
            np.pi,
        ]
        if not translation:
            bounds[0] = bounds[1] = bounds[2] = 0
        if not rotation:
            bounds[3] = bounds[4] = bounds[5] = 0

        for cdx, cycle in enumerate(self.cycles):

            anchor = self.cycle_anchors[cdx]

            if len(self.cycles) > 1:
                other_cycles = self.cycle_coords[
                    [False if c == cdx else True for c in range(len(self.cycles))]
                ]

                def constraint(env, coords):
                    dist_anchor = ((np.mean(coords, axis=0) - anchor) ** 2).sum()
                    dist_others = 0  # TBA

                    clashes = (cdist(coords, self.axle_coords) < 3).sum() ** 5

                    return dist_anchor + dist_others + clashes

            else:

                def constraint(env, coords):
                    dist_anchor = ((np.mean(coords, axis=0) - anchor) ** 2).sum()
                    clashes = (cdist(coords, self.axle_coords) < 3).sum() ** 5

                    return dist_anchor + clashes

                translatron = optimizers.Translatron(
                    cycle._AtomGraph, constraint, bounds=bounds
                )
                translatron.step(translatron.action_space.sample())

                solution, _ = optimizers.scipy_optimize(translatron)

                if rotation:
                    cycle.rotate(solution[3], "x", angle_is_degrees=False)
                    cycle.rotate(solution[4], "y", angle_is_degrees=False)
                    cycle.rotate(solution[5], "z", angle_is_degrees=False)
                if translation:
                    cycle.move_to(solution[:3])

                self.cycle_coords[cdx] = cycle.get_coords()

    def merge(self) -> core.Molecule:
        """
        Merge the cycles and axle into a single molecule.
        Specifically, the cycles are merges as new chains into the axle molecule.
        """
        for cycle in self.cycles:
            self.axle.merge(cycle)
        return self.axle

    def _set_initial_cycle_positions(self):
        axis = self._find_longest_axis()
        self.main_axis = axis

        # get some anchor atom of the axle
        anchor = next(self.axle.get_atoms()).get_coord()

        # get the maximum distance along the axis
        max_distance = np.max(cdist(self.axle_coords, axis.reshape(1, -1)))

        # Generate the broad anchor positions for the cycles
        self.cycle_anchors = np.array(
            [
                (anchor + (i + 1) * max_distance / (len(self.cycles) + 1) * axis)
                for i in range(len(self.cycles))
            ]
        )
        for cycle, pos in zip(self.cycles, self.cycle_anchors):
            cycle.move_to(pos)

    def _rotate_cycles_perpendicular_to_main_axis(self):
        for anchor, cycle in zip(self.cycle_anchors, self.cycles):

            cycle_coords = cycle.get_coords()
            cycle_plane_vector = structural.plane_of_points(cycle_coords)

            # compute the angle between the plane vector and longest axis
            angle = np.pi / 2 + np.dot(self.main_axis, cycle_plane_vector)

            # compute the axis to rotate around perpendicularly to the other two vectors
            axis = np.cross(self.main_axis, cycle_plane_vector)
            axis /= np.linalg.norm(axis)

            # now rotate the cycle coords by the angle around the axis
            cycle_coords = (
                structural.rotate_coords(cycle_coords - anchor, angle, axis) + anchor
            )

            for atom, coord in zip(cycle.get_atoms(), cycle_coords):
                atom.coord = coord

    def _find_longest_axis(self):

        # Perform Principal Component Analysis (PCA) on the coordinates
        pca = PCA(n_components=3)
        pca.fit(self.axle_coords)

        # Get the eigenvectors (axes) of the PCA
        axes = pca.components_

        # Find the longest axis
        longest_axis = axes[np.argmax(pca.explained_variance_)]

        return longest_axis


if __name__ == "__main__":

    import buildamol.extensions.polymers as p

    # make the ring
    ring = p.cyclic_alkane(10)

    # make the axle basis and then modify the atoms
    N = 35
    axle = p.linear_alkane(N)
    axle.get_hydrogen("C1").set_element("Br")
    axle.get_hydrogen(f"C{N}").set_element("Br")
    i = 4
    while i < 35:
        # here we need to use change_element as the hydrogen count changes
        axle.change_element(f"C{i}", "N")
        i += 5

    # cycle.rotate(90, "y")

    builder = RotaxaneBuilder()
    builder.distribute_along_axle(axle, ring.copy(2))
    # builder.optimize()

    builder.merge().to_pdb("rotaxane_builder.pdb")

    # out = rotaxane(axle, cycle.copy(2))
    # out.show()
