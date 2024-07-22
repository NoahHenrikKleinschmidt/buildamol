"""
This is the Translatron environment that can be used to place a molecule according to some constraints.
"""

import gym
from gym import spaces
import numpy as np

import buildamol.utils.auxiliary as aux
import buildamol.structural.base as structural

__all__ = ["Translatron"]


class Translatron(gym.Env):
    """
    The Translatron environment can be used to place a molecule in space according to some constraints.
    It will produce a vector of 6 values, the first 3 are the translation in x, y, and z, and the last 3 are the
    rotations around the x, y, and z axes.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph of the molecule to be placed.
    constraint_func : callable
        A function that takes the environment and the new coordinates and returns a reward value.
    finish_func : callable
        A function that takes the environment and the new coordinates and returns a boolean value indicating
        whether the optimization is finished.
    bounds : tuple
        The bounds for the minimal and maximal translation and rotation values.
        If a tuple of length two this is interpreted as the low and high bounds for translation only.
        Otherwise provide a tuple of length 6 for the bounds of translation and rotation. In this case values can be either
        singular (int/float) in which case they are interpreted as symmetric extrama (min=-value, max=+value) or as tuples with (min=value[0], max=value[1]).
        Mixed inputs are allowed.
    """

    def __init__(
        self,
        graph,
        constraint_func: callable,
        finish_func: callable = None,
        bounds: tuple = (-10, 10),
    ):
        self.graph = graph
        self.constraint_func = constraint_func
        self.finish_func = finish_func or self.__no_finish

        if len(bounds) == 2:
            low = [bounds[0], bounds[0], bounds[0], -np.pi, -np.pi, -np.pi]
            high = [bounds[1], bounds[1], bounds[1], np.pi, np.pi, np.pi]
        elif len(bounds) == 6:
            bounds = [(-i, i) if isinstance(i, (int, float)) else i for i in bounds]
            low = [
                bounds[0][0],
                bounds[1][0],
                bounds[2][0],
                bounds[3][0],
                bounds[4][0],
                bounds[5][0],
            ]
            high = [
                bounds[0][1],
                bounds[1][1],
                bounds[2][1],
                bounds[3][1],
                bounds[4][1],
                bounds[5][1],
            ]
        else:
            raise ValueError(
                f"Expected 'bounds' to be of length 2 (only translation min/max) or length 6 (min/max for translation and rotations), but got length {len(bounds)}"
            )

        self.action_space = spaces.Box(
            low=np.array(low),
            high=np.array(high),
            shape=(6,),
            dtype=np.float32,
        )
        self._setup(graph)

        if aux.USE_NUMBA:
            self.apply = self._numba_apply
        else:
            self.apply = self._normal_apply

    def reset(self, *args, **kwargs):
        self.state[:] = self.__backup_coords

    def step(self, action):
        new_coords = self.apply(action)
        reward = self.constraint_func(self, new_coords)
        done = self.finish_func(self, new_coords)
        return new_coords, reward, done, {}

    def eval(self, state):
        return self.constraint_func(self, state)

    def _normal_apply(self, action):
        """
        Apply an action to the current state.
        """
        new = structural.rotate_coords(self.state, action[3], structural.x_axis)
        new = structural.rotate_coords(new, action[4], structural.y_axis)
        new = structural.rotate_coords(new, action[5], structural.z_axis)
        new = new + action[:3]
        return new

    def _numba_apply(self, action):
        new = structural._numba_wrapper_rotate_coords(
            self.state, action[3], structural.x_axis
        )
        new = structural._numba_wrapper_rotate_coords(new, action[4], structural.y_axis)
        new = structural._numba_wrapper_rotate_coords(new, action[5], structural.z_axis)
        new = new + action[:3]
        return new

    def blank(self):
        return np.zeros(6)

    def _setup(self, graph):
        coords = np.array([i.coord for i in graph.nodes])
        self.state = coords
        self.__backup_coords = coords.copy()

    def __no_finish(self, coords, *args, **kwargs):
        return False


if __name__ == "__main__":

    import buildamol as bam

    mol = bam.read_pdb("/Users/noahhk/GIT/biobuild/0.pdb")

    mol2 = mol.copy()
    mol2.move([10, 10, 10])
    mol2.rotate(90, [1, 0, 0])
    mol2.rotate(36, [0, 1, 0])

    ref_coords = np.array([i.coord for i in mol.get_atoms()])

    def constraint(env, coords):
        dist = (coords - ref_coords) ** 2
        dist = np.mean(dist)
        env._dist = dist
        return dist

    def finish(env, coords):
        return env._dist < 0.1

    env = Translatron(mol2.get_atom_graph(), constraint)

    sol, eval = bam.optimizers.swarm_optimize(env, n_particles=50)

    v = mol.draw(atoms=False)

    mol2.move(sol[:3])
    mol2.rotate(sol[3], "x")
    mol2.rotate(sol[4], "y")
    mol2.rotate(sol[5], "z")

    v.draw_edges(*mol2.get_bonds(), color="blue", linewidth=3, showlegend=False)
    v.show()
