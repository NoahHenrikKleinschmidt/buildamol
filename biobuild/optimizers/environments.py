import numpy as np
from copy import deepcopy

# import gymnasium as gym
import gym

from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation

import biobuild.structural as structural
import biobuild.utils.auxiliary as aux


class Rotatron(gym.Env):
    """
    The base class for rotation environments

    Parameters
    ----------
    graph : biobuild.graphs.Graph,
        The graph to optimize
    rotatable_edges : list, optional,
        A list of edges to rotate around
    """

    def __init__(self, graph, rotatable_edges=None) -> None:
        super().__init__()
        if rotatable_edges is None:
            rotatable_edges = [
                edge for edge in graph.edges if not graph.is_locked(*edge)
            ]

        self.graph = graph
        self.rotatable_edges = [i for i in rotatable_edges if i in graph.edges]

        self.graph = deepcopy(self.graph)
        # self.rotatable_edges = deepcopy(self.rotatable_edges)

        self.n = len(self.rotatable_edges)
        self.reward = 0

        self._make_coord_array()
        self._make_descendant_masks()
        self._make_edge_ref_masks()

        self._effector_mask = np.ones(len(self._nodes), dtype=bool)
        self.__orig_coords = deepcopy(self._coords)
        self.__renderer__ = None

    @classmethod
    def load(cls, path):
        """
        Load the environment from a file

        Parameters
        ----------
        path : str
            The path to the file
        """
        return aux.load_pickle(path)

    @property
    def effector_coords(self):
        """
        Get the effector coordinates
        """
        return self._coords[self._effector_mask]

    @property
    def effector_nodes(self):
        """
        Get the effector nodes
        """
        return np.array(self._nodes)[self._effector_mask]

    @property
    def state(self):
        """
        Get the current state of the environment
        """
        return self._coords

    @state.setter
    def state(self, state):
        """
        Set the current state of the environment
        """
        self._coords = deepcopy(state)

    def save(self, path):
        """
        Save the environment to a file

        Parameters
        ----------
        path : str
            The path to the file
        """
        aux.save_pickle(self, path)

    def set_state(self, state):
        """
        Set the current state of the environment
        """
        self.state = state

    def get_state(self):
        """
        Get the current state of the environment
        """
        return self.state

    def reset(self, graph: bool = False):
        """
        Reset the environment
        """
        self._coords = deepcopy(self.__orig_coords)
        if graph:
            self.apply_to_graph()
        return self.state

    def render(self, **kwargs):
        """
        Render the environment
        """
        raise NotImplementedError

    def step(self, action):
        """
        Take a step in the environment
        """
        raise NotImplementedError

    def apply_to_graph(self):
        """
        Apply the current state to the graph
        """
        raise NotImplementedError
        for i, node in enumerate(self._nodes):
            node.coord = deepcopy(self._coords[i])

    def mask_effector_coords(self, mask):
        """
        Mask the effector coordinates for which the reward is calculated.
        By default, all coordinates are used.

        Parameters
        ----------
        mask : np.ndarray
            The mask to apply
        """
        self._effector_mask = mask

    def get_descendant_coords(self, bond: int):
        """
        Get the descendant coordinates for a given bond

        Parameters
        ----------
        bond : int
            The bond to get the descendants for

        Returns
        -------
        descendant_coords : np.ndarray
            The descendant coordinates
        """
        return self._coords[self._descendant_masks[bond]]

    def rotate(self, bond: int, angle: float):
        """
        Rotate a coordinate array around an axis by a specified angle

        Parameters
        ----------
        bond : int
            The bond to rotate around (as sampled from the action space)
        angle : float
            The angle to rotate by in radians

        Returns
        -------
        rotated_coords : np.ndarray
            The rotated coordinates
        desc_mask : np.ndarray
            The mask of the descendant nodes
        """
        desc_mask = self._descendant_masks[bond]
        coords = self._coords[desc_mask]

        # Get the axis to rotate around
        edge_masks = self._edge_ref_masks[bond]
        a = self._coords[edge_masks[0]][0]
        b = self._coords[edge_masks[1]][0]
        axis = b - a
        axis = axis / np.linalg.norm(axis)

        # get the reference coordinate for the bond
        ref_coord = b

        # translate the coordinates so that the reference coordinate is at the origin
        coords = coords - ref_coord

        axis = angle * axis

        # Get the rotation matrix
        R = Rotation.from_rotvec(axis)

        # Rotate the coordinates
        rotated_coords = R.apply(coords)

        # Translate the coordinates back to the original position
        rotated_coords = rotated_coords + ref_coord

        return rotated_coords, desc_mask

    def blank(self):
        """
        Return a blank action
        """
        raise NotImplementedError

    def compute_reward(self, coords=None):
        """
        Compute the reward of the given or current coordinate array
        """
        raise NotImplementedError

    def _make_coord_array(self):
        """
        Make the coordinate array
        """
        self._nodes = list(self.graph.nodes)
        self._coords = np.array([node.coord for node in self._nodes], dtype=float)

    def _make_descendant_masks(self):
        """
        Map the descendant nodes for each rotatable edge
        """
        # masks = {
        #     # the list comprehension is pretty slow
        #     # idx: np.array([j in self.graph.get_descendants(*edge) for j in self._nodes])
        #     idx: np.fromiter(
        #         # WE CAN IMPROVE THIS LINE BY DOING SOMETHING LIKE
        #         # FOR J IN SELF._NODES IF J NOT AN INTERNAL NODE
        #         # ELSE JUST USE THE ALL FALSE MASK...
        #         (j in self.graph.get_descendants(*edge) for j in self._nodes),
        #         dtype=bool,
        #         count=len(self._nodes),
        #     )
        #     for idx, edge in enumerate(self.rotatable_edges)
        # }
        masks = {}
        for idx, edge in enumerate(self.rotatable_edges):
            mask = np.zeros(len(self._nodes), dtype=bool)
            for j in self.graph.get_descendants(*edge):
                mask[self._nodes.index(j)] = True
            masks[idx] = mask
        self._descendant_masks = masks

    def _make_edge_ref_masks(self):
        """
        Map the coordinates of one end of each rotatable edge from the _coords array
        """
        masks = {
            idx: (
                np.array([j == node1 for j in self._nodes]),
                np.array([j == node2 for j in self._nodes]),
            )
            for idx, (node1, node2) in enumerate(self.rotatable_edges)
        }
        self._edge_ref_masks = masks


class MultiBondRotatron(Rotatron):
    """
    This environment samples an angle for each rotatable bond and evaluates
    the reward of the resulting structure

    Parameters
    ----------
    graph
        A detailed ResidueGraph object
    rotatable_edges
        A list of rotatable edges. If None, all non-locked edges from the graph are used.
    mask_same_residues: bool
        Whether to mask the coordinates of the same residue (ignore the distances between atoms of the same residues, since they never change).
    mask_max_distance: float
        Ignore the distances between atoms if they are greater than this value.
    """

    def __init__(
        self,
        graph,
        rotatable_edges=None,
        mask_same_residues: bool = True,
        mask_max_distance: float = 10.0,
    ):
        super().__init__(graph, rotatable_edges)
        self.action_space = gym.spaces.Box(
            low=-np.pi, high=np.pi, shape=(len(self.rotatable_edges),)
        )
        self.observation_space = gym.spaces.Box(
            low=-np.inf, high=np.inf, shape=(len(self._nodes), 3)
        )
        self.mask_same_residues = mask_same_residues
        self.mask_max_distance = mask_max_distance

        # Mask the effector coordinates
        if mask_same_residues:
            self._mask_residues()

    def blank(self):
        """
        Return a blank action
        """
        return np.zeros(len(self.rotatable_edges))

    def step(self, action):
        """
        Take a step in the environment
        """
        self.reset()
        angles = action
        for idx in range(len(angles)):
            angle = angles[idx]
            new_coords, mask = self.rotate(idx, angle)
            self._coords[mask] = new_coords
        reward = self.compute_reward()

        return self.state, reward, False, {}

    def _regional_mask(self, dists):
        """
        Compute a mask for the distances that are greater than the max distance
        """
        return dists > self.mask_max_distance

    def _reduce_func(self, dists):
        """
        Compute the mean of the distances that are not masked
        """
        return dists[dists > 0].mean()

    def compute_reward(self):
        """
        Compute the reward of the given or current coordinate array
        """

        coords = self.effector_coords

        # Compute the inter-residue distances
        dists = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

        # Mask the residue distances
        np.fill_diagonal(dists, -1)
        dists[self._residue_masks] = -1

        # Ignore distances greater than 10 angstroms
        regional_mask = np.apply_along_axis(self._regional_mask, 1, dists)
        dists[regional_mask] = -1

        # Reduce the distances to a single value per node
        dists = np.apply_along_axis(self._reduce_func, 1, dists)

        # energy = (1 / dists) ** 12 - (1 / dists) ** 6
        reward = -4 * np.sum((1 / dists) ** 4)

        return reward

    def _mask_residues(self):
        """
        For each node in the graph, mask which nodes are part of the same
        residue as these are not relevant for reward computation and could
        interfere with the learning process
        """
        self._residue_masks = np.empty((len(self._nodes), len(self._nodes)), dtype=bool)
        for residue in self.graph.residues:
            node_mask = np.array(
                [node is residue or node in residue.child_list for node in self._nodes]
            )
            node_mask[self._nodes.index(residue)] = True
            for node in residue.child_list:
                if node not in self._nodes:
                    continue
                self._residue_masks[self._nodes.index(node)] = node_mask
            self._residue_masks[self._nodes.index(residue)] = node_mask

        for edge in self.rotatable_edges:
            a, b = edge
            self._residue_masks[self._nodes.index(a), self._nodes.index(b)] = True
            self._residue_masks[self._nodes.index(b), self._nodes.index(a)] = True


class SphereRotatron(MultiBondRotatron):
    """
    This environment samples an angle for each rotatable bond and evaluates
    the reward of the resulting structure. It uses a spherical representation
    of residues to evaluate the reward.

    Parameters
    ----------
    graph
        A detailed ResidueGraph object
    rotatable_edges
        A list of rotatable edges. If None, all non-locked edges from the graph are used.
    """

    def __init__(self, graph, rotatable_edges=None):
        super().__init__(graph, rotatable_edges)
        self._effector_mask = np.array(
            [node in self.graph.residues for node in self._nodes]
        )
        self._residue_radii = np.array(
            [
                structural.compute_residue_radius(node)
                for node in self._nodes
                if node in self.graph.residues
            ]
        )
        self._residue_radii_sums = (
            self._residue_radii[:, np.newaxis] + self._residue_radii
        )

    def compute_reward(self):
        """
        Compute the reward of the given or current coordinate array
        """

        coords = self.effector_coords

        # Compute the inter-residue distances
        dists = cdist(coords, coords) - self._residue_radii_sums
        reward = -np.sum((1 / dists) ** 4)
        return reward


class DiscreteRotatron(MultiBondRotatron):
    """
    This environment uses a discrete sample space for rotational angles instead of a continuous one.
    This helps reduce the complexity of the environment and makes it easier to learn.

    Parameters
    ----------
    graph
        A detailed ResidueGraph object
    rotatable_edges
        A list of rotatable edges. If None, all non-locked edges from the graph are used.
    mask_same_residues: bool
        Whether to mask the nodes of the same residue from the reward computation.
    d: int
        The angle discretization to use, given in degrees.
        1 degree by default, meaning 360 possible actions between -180 and +180 degrees.
    """

    def __init__(
        self,
        graph,
        rotatable_edges=None,
        mask_same_residues: bool = True,
        d: int = 1,
    ) -> None:
        super().__init__(graph, rotatable_edges, mask_same_residues)

        d = np.radians(d)
        self._angles = np.arange(-np.pi, np.pi, d)
        self.action_space = gym.spaces.MultiDiscrete(
            [len(self._angles) for _ in self.rotatable_edges]
        )

    def step(self, action):
        angles = self._angles[action.astype(int)]
        return super().step(angles)


# if __name__ == "__main__":
#     import biobuild as bb
#     import jax.numpy as np
#     import jax

#     mol = bb.Molecule.load("/Users/noahhk/GIT/biobuild/test.mol")

#     connections = sorted(mol.get_residue_connections())
#     graph = mol.make_residue_graph()
#     graph.make_detailed(True, True, f=0.8)

#     env = MultiBondRotatron(graph, connections)

#     reward_f = lambda x: env.step(x)[1]

#     def jax_reward(x_jax):
#         reward = jax.pure_callback(
#             reward_f, jax.core.ShapedArray((1, 1), np.float32), x_jax
#         )
#         return np.array(reward)

#     grad_reward = grad(jax_reward)

#     grad_reward(np.array(env.action_space.sample()))
