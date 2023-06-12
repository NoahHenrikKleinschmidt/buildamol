"""
Intelligent agents for optimizing molecules
"""

from typing import Union
import numpy as np

import biobuild.optimizers.environments as environments

import pyswarms as ps
import scipy.optimize as opt

from alive_progress import alive_bar


def scipy_optimize(
    env, steps: int = 1e5, method: str = "L-BFGS-B", bounds=(-np.pi, np.pi), **kws
):
    """
    Optimize a Rotatron environment through
    a simple scipy optimization

    Parameters
    ----------
    env : biobuild.optimizers.environments.Rotatron
        The environment to optimize
    steps : int, optional
        The number of steps to take, by default 1000
    method : str, optional
        The optimizer to use, by default "L-BFGS-B".
        This can be any optimizer from scipy.optimize.minimize
    bounds : tuple, optional
        The bounds to use for solutions as a tuple of size 2 with a minimal and maximal angle to allow.
    kws : dict, optional
        Keyword arguments to pass as options to the optimizer

    Returns
    -------
    np.ndarray
        The optimized action
    float
        The reward for the optimized action
    """

    x0 = env.action_space.sample()

    def loss_fn(x):
        state, reward, *_ = env.step(x)
        return -reward

    if bounds:
        bounds = opt.Bounds(*bounds)

    kws["maxiter"] = int(steps)
    result = opt.minimize(loss_fn, x0, method=method, bounds=bounds, options=kws)
    return result.x, -result.fun


def swarm_optimize(env, steps: int = 30, n: int = 10, bounds=(-np.pi, np.pi), **kws):
    """
    Optimize a Rotatron environment through
    a pyswarms swarm optimization

    Parameters
    ----------
    env : biobuild.optimizers.environments.Rotatron
        The environment to optimize
    steps : int, optional
        The number of steps to take, by default 1000
    n : int, optional
        The number of particles to use, by default 10
    bounds : tuple, optional
        The bounds to use for solutions as a tuple of size 2 with a minimal and maximal angle to allow.
    kws : dict, optional
        Keyword arguments to pass as options to the optimizer

    Returns
    -------
    np.ndarray
        The optimized action
    float
        The reward for the optimized action
    """

    x0 = env.action_space.sample()

    if bounds:
        bounds = (np.full(x0.shape, bounds[0]), np.full(x0.shape, bounds[1]))

    def loss_fn(sol):
        costs = np.zeros(n)
        i = 0
        for s in sol:
            costs[i] = -env.step(s)[1]
            i += 1
        return costs

    options = kws.pop("options", {})
    options.setdefault("c1", 0.5)
    options.setdefault("c2", 0.3)
    options.setdefault("w", 0.9)

    verbose = kws.pop("verbose", False)

    optimizer = ps.single.GlobalBestPSO(
        n_particles=n, dimensions=x0.shape[0], bounds=bounds, options=options, **kws
    )
    reward, result = optimizer.optimize(loss_fn, iters=steps, verbose=verbose)
    return result, -reward


# class Rotator:
#     """
#     The base agent class for working with Rotatron problems

#     Parameters
#     ----------
#     env : biobuild.optimizers.environments.Rotatron
#         The environment to optimize
#     searcher : evotorch.algorithms.Searcher, optional
#         A search algorithm to use, by default None.
#         Each child class must define a searcher.
#     verbose : bool, optional
#         Whether to show a progress bar, by default False
#     **kwargs
#         Keyword arguments for the searcher
#     """

#     def __init__(self, env, searcher=None, verbose: bool = False, **kwargs) -> None:
#         self.problem = env
#         if searcher is not None:
#             self.searcher = searcher(self.problem, **kwargs)
#         else:
#             self.searcher = None
#         if verbose:
#             self.bar = alive_bar
#         else:
#             self.bar = aux.DummyBar

#     def render(self, *args, **kwargs) -> None:
#         """
#         Render the environment
#         """
#         self.problem.render(*args, **kwargs)

#     def reset(self) -> None:
#         """
#         Reset the environment
#         """
#         self.problem.reset()

#     def step(self, *args, **kwargs) -> None:
#         """
#         Run a single step of the environment
#         """
#         if not self.searcher:
#             raise ValueError("No searcher defined")
#         self.problem.step(*args, **kwargs)

#     def run(self, steps: int, *args, **kwargs) -> None:
#         """
#         Run the agent

#         Parameters
#         ----------
#         steps : int
#             The number of steps to run the agent
#         *args
#             Arguments for the searcher
#         **kwargs
#             Keyword arguments for the searcher

#         Returns
#         -------
#         tuple
#             The best found solution and its obtained reward
#         """
#         if not self.searcher:
#             raise ValueError("No searcher defined")

#         with self.bar(steps) as bar:
#             for _ in range(steps):
#                 self.searcher.step(*args, **kwargs)
#                 bar()
#         return self.best

#     @property
#     def best(self) -> tuple:
#         """
#         Get the best found solution and its obtained reward
#         """
#         if not self.searcher:
#             return None, None
#         else:
#             return (
#                 self.searcher.status["best"].values.numpy(),
#                 self.searcher.status["best_eval"],
#             )

#     def __repr__(self) -> str:
#         return f"{self.__class__.__name__}({self.problem})"

#     def __call__(self, steps: int, *args, **kwds):
#         self.run(steps, *args, **kwds)
#         return self.best


# class EvoRotator(Rotator):
#     """
#     This rotator uses an evolutionary algorithm to search for optimal conformations

#     Parameters
#     ----------
#     env : biobuild.optimizers.environments.Rotatron
#         The environment to optimize
#     verbose : bool, optional
#         Whether to show a progress bar, by default False
#     population_size : int, optional
#         The population of conformations to maintian, by default 20
#     mutation_rate : float, optional
#         The mutation rate, by default 0.25 (25%)
#     mutation_stdev : float, optional
#         The standard deviation of the mutation of angles, by default 0.8,
#         in radians, confinded to [-pi, pi].
#     """

#     def __init__(
#         self,
#         env,
#         verbose: bool = False,
#         population_size: int = 20,
#         mutation_rate: float = 0.25,
#         mutation_stdev: float = 0.8,
#         **kwargs,
#     ) -> None:
#         searcher = ea.Cosyne
#         kwargs.setdefault("popsize", population_size)
#         kwargs.setdefault("mutation_probability", mutation_rate)
#         kwargs.setdefault("mutation_stdev", mutation_stdev)
#         kwargs.setdefault("tournament_size", int(max(1, population_size / 5)))
#         super().__init__(env, searcher, verbose, **kwargs)


# class GradientRotator(Rotator):
#     """
#     This rotator uses gradient descent to search for optimal conformations

#     Parameters
#     ----------
#     env : biobuild.optimizers.environments.Rotatron
#         The environment to optimize
#     verbose : bool, optional
#         Whether to show a progress bar, by default False
#     lr : float, optional
#         The learning rate, by default 0.1
#     """

#     def __init__(self, env, verbose: bool = False, lr: float = 0.1, **kwargs) -> None:
#         kwargs.setdefault("lr", lr)
#         super().__init__(env, None, verbose, **kwargs)
#         self.searcher = optim.SGD()
#     def step(self, *args, **kwargs) -> None:
#         """
#         Run a single step of the environment
#         """
#         if not self.searcher:
#             raise ValueError("No searcher defined")

if __name__ == "__main__":
    import biobuild as bb
    from time import time

    glc = bb.Molecule.from_compound("GLC")
    lac = glc.repeat(2, "14bb")

    c1 = lac.get_atom("C1", residue=2)
    o4 = lac.get_atom("O4", residue=1)
    c4 = lac.get_atom("C4", residue=1)
    o5 = lac.get_atom("O5", residue=2)
    c5 = lac.get_atom("C5", residue=1)

    _dihedral_c1_o4 = (o5, c1, o4, c4)
    _dihedral_o4_c4 = (c1, o4, c4, c5)

    dihedral_c1_o4 = lac.compute_dihedral(*_dihedral_c1_o4)
    dihedral_o4_c4 = lac.compute_dihedral(*_dihedral_o4_c4)

    v = lac.draw()
    v.draw_edges(lac.bonds, color="limegreen", linewidth=5)

    lac.rotate_around_bond(c1, o4, -dihedral_c1_o4, descendants_only=True)
    lac.rotate_around_bond(o4, c4, -dihedral_o4_c4, descendants_only=True)

    edges = lac.get_residue_connections()
    v.draw_edges(lac.bonds, color="black", linewidth=4)
    v.draw_edges(edges, color="magenta", linewidth=4)

    # ====================================================================

    graph = lac.make_residue_graph()
    graph.make_detailed(True, True)
    env = bb.optimizers.MultiBondRotatron(graph, edges)

    # ====================================================================
    colors = ["red", "orange", "yellow", "green", "blue", "purple"]
    idx = 0
    for steps in (20, 30, 40, 50, 60):
        t1 = time()
        for i in range(20):
            sol, reward = swarm_optimize(env, steps=steps)

            lac.rotate_around_bond(c1, o4, np.degrees(sol[0]), descendants_only=True)
            lac.rotate_around_bond(o4, c4, np.degrees(sol[1]), descendants_only=True)

            v.draw_edges(lac.bonds, color=colors[idx], opacity=0.2, linewidth=4)

            lac.rotate_around_bond(c1, o4, -np.degrees(sol[0]), descendants_only=True)
            lac.rotate_around_bond(o4, c4, -np.degrees(sol[1]), descendants_only=True)
        idx += 1
        t2 = time()
        print(f"Steps: {steps}, Time: {t2-t1}")
    v.show()
    pass

# if __name__ == "__main__":
#     import biobuild as bb
#     from copy import deepcopy
#     from time import time
#     import stable_baselines3 as sb3
#     import halo

#     mol_ref = bb.Molecule.from_pdb(
#         "/Users/noahhk/GIT/biobuild/support/examples/man9.pdb"
#     )
#     mol_ref.infer_bonds(restrict_residues=False)

#     mol = bb.Molecule.load("/Users/noahhk/GIT/biobuild/test.mol")
#     # mol1 = deepcopy(mol)

#     # angles = np.random.randint(-180, 180, size=(len(mol.get_residue_connections())))
#     connections = sorted(mol.get_residue_connections())
#     # for i, angle in enumerate(angles):
#     #     mol.rotate_around_bond(*connections[i], angle, descendants_only=True)

#     g = mol.make_residue_graph()
#     g.make_detailed(True, True, f=1)

#     v2 = bb.utils.visual.MoleculeViewer3D(mol.make_residue_graph(True))
#     v = bb.utils.visual.MoleculeViewer3D(mol)

#     v.draw_edges(mol_ref.bonds, color="cyan")
#     v2.draw_edges(mol_ref.make_residue_graph(True).edges, color="cyan")

#     steps = (1e5,)  # , 1e4, 1e5, 1e6)
#     colors = ("green",)  # "orange", "red", "purple")
#     for step, color in zip(steps, colors):
#         for i in range(3):
#             mol1 = deepcopy(mol)

#             spinner = halo.Halo(text="Optimizing", spinner="dots")
#             spinner.start()

#             t1 = time()
#             env = environments.MultiBondRotatron(g, connections)
#             result = scipy_optimize(env, steps=step)  # , method="Nelder-Mead")
#             t2 = time()

#             spinner.succeed(f"scipy Optimization finished in {t2-t1:.2f} seconds")

#             for i, angle in enumerate(result):
#                 bond = env.rotatable_edges[i]
#                 bond = (bond[0].serial_number, bond[1].serial_number)
#                 mol1.rotate_around_bond(*bond, np.degrees(angle), descendants_only=True)

#             v.draw_edges(mol1.bonds, color=color, opacity=0.3)
#             v2.draw_edges(mol1.make_residue_graph(True).edges, color=color, opacity=0.3)

#     v.show()
#     v2.show()

#     # r = EvoRotator(env, verbose=True)
#     # best, reward = r.run(5000)

#     # for i, angle in enumerate(best):
#     #     bond = env.rotatable_edges[i]
#     #     bond = (bond[0].serial_number, bond[1].serial_number)
#     #     mol1.rotate_around_bond(*bond, np.degrees(angle), descendants_only=True)

#     # v.draw_edges(mol1.bonds, color="limegreen")
#     # v.show()
#     pass
