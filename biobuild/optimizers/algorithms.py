import numpy as np
import scipy.optimize as opt


def scipy_optimize(
    env,
    steps: int = 1e5,
    method: str = "L-BFGS-B",
    bounds=(-np.pi, np.pi),
    **kws,
):
    """
    Optimize a Rotatron environment through
    a simple scipy optimization

    Parameters
    ----------
    env : biobuild.optimizers.environments.Rotatron
        The environment to optimize
    steps : int, optional
        The number of steps to take, by default 1e5
    method : str, optional
        The optimizer to use, by default "L-BFGS-B".
        This can be any optimizer from scipy.optimize.minimize
    bounds : tuple, optional
        The bounds to use for solutions as a tuple of size 2 with a minimal and maximal angle to allow.
    kws : dict, optional
        Keyword arguments to pass as options to the optimizer

    Returns
    -------
    tuple or Molecule
        If molecule is None, a tuple of the optimized action and the evaluation for the optimized action.
        If molecule is not None, the optimized molecule.
    """

    x0 = env.action_space.sample()

    def loss_fn(x):
        state, _eval, *_ = env.step(x)
        env.reset()
        return _eval

    if bounds:
        bounds = opt.Bounds(*bounds)

    kws["maxiter"] = int(steps)
    result = opt.minimize(loss_fn, x0, method=method, bounds=bounds, options=kws)
    return (
        result.x,
        result.fun,
    )


def genetic_optimize(
    env,
    max_steps: int = 1e3,
    stop_if_done: bool = True,
    threshold: float = 1e-4,
    variation: float = 0.3,
    population_size: int = 20,
    n_parents: int = 9,
    n_children: int = 8,
    n_mutations: int = 2,
    n_newcomers: int = 1,
):
    """
    A simple genetic algorithm for optimizing a Rotatron environment.

    Parameters
    ----------
    env : biobuild.optimizers.environments.Rotatron
        The environment to optimize
    stop_if_done : bool, optional
        Stop the optimization if the environment signals it is done or the solutions have converged.
    threshold : float, optional
        A thershold to use for convergence of the best solution found.
        The algorithm will stop if the variation of the best solution evaluation history
        is less than this threshold.
    max_steps : int, optional
        The maximum number of steps to take.
    variation : float, optional
        The variation to use for the initial action.
    population_size : int, optional
        The size of the population.
    n_parents : int, optional
        The number of parents to select. The parents are selected from the best solutions.
    n_children : int, optional
        The number of children to generate from the parents.
    n_mutations : int, optional
        The number of mutations to apply. Mutations are applied to the selected parents,
        generating abarrent clones.
    n_newcomers: int, optional
        Newcomers are entirely new solution candidates.

    Returns
    -------
    solution, evaluation
        The solution and evaluation for the solution
    """
    if n_children + n_mutations + n_parents + n_newcomers != population_size:
        raise ValueError(
            "The sum of n_children, n_mutations, n_parents, and n_newcomers must equal population_size"
        )

    if any([n < 1 for n in [n_children, n_mutations, n_parents, n_newcomers]]):
        raise ValueError(
            "n_children, n_mutations, n_parents, and n_newcomers must all be at least 1"
        )

    blank = env.blank()
    population = np.stack([blank] * population_size)
    evals = np.zeros(population_size)

    bests = np.zeros(max(10, int(max_steps * 0.05)))

    pop_range = np.arange(0, population_size)
    children_range = np.arange(0, n_children)
    mutations_range = np.arange(0, n_mutations)
    newcomers_range = np.arange(0, n_newcomers)

    for i in pop_range:
        population[i] = env.action_space.sample()
        _, evals[i], *_ = env.step(population[i])
        env.reset()

    steps = 0
    while steps < max_steps:
        sorting = np.argsort(evals)
        parents = population[sorting[:n_parents]]

        children = np.stack([blank] * n_children)
        for i in children_range:
            p1, p2 = np.random.choice(np.arange(0, n_parents), size=2, replace=False)
            p1, p2 = parents[int(p1)], parents[int(p2)]
            children[i] = (p1 + p2) / 2 + np.random.normal(
                0, variation / 2, size=blank.shape
            )

        mutations = np.stack([blank] * n_mutations)
        for i in mutations_range:
            mutations[i] = parents[np.random.randint(0, n_parents)] + np.random.normal(
                0, variation, size=blank.shape
            )

        newcomers = np.stack([blank] * n_newcomers)
        for i in newcomers_range:
            newcomers[i] = env.action_space.sample()

        population = np.concatenate([parents, children, mutations, newcomers])
        evals = np.zeros(population_size)

        for i in pop_range:
            _, evals[i], done, *_ = env.step(population[i])
            env.reset(best=False)
            if done and stop_if_done:
                break

        if done and stop_if_done:
            break

        bests = np.roll(bests, -1)
        bests[-1] = env._best_eval
        if stop_if_done and np.var(bests) < threshold:
            break

        steps += 1

    best = np.argmin(evals)
    return population[best], evals[best]


def swarm_optimize(
    env,
    n_particles: int = 10,
    max_steps: int = 30,
    stop_if_done: bool = True,
    threshold: float = 1e-5,
    variation: float = 0.1,
    recycle: float = 0.3,
    recycle_every: int = 5,
):
    """
    Optimize a rotatron environment through a simple particle swarm optimization.

    Parameters
    ----------
    env : biobuild.optimizers.environments.Rotatron
        The environment to optimize
    n_particles : int, optional
        The number of particles to use.
    max_steps : int, optional
        The maximum number of steps to take.
    stop_if_done : bool, optional
        Stop the optimization if the environment signals it is done or the solutions have converged.
    threshold : float, optional
        A threshold to use for convergence of the best solution found.
        The algorithm will stop if the variation of the best solution evaluation history
        is less than this threshold.
    variation : float, optional
        The variation to use for updating particle positions.
    recycle : float, optional
        The fraction of particles to replace by variations of the best particle when updating the particle positions.
        This will remove the worst particles and replace them with the best particle + some noise.
    recycle_every : int, optional
        The number of steps to take before recycling particles.

    Returns
    -------
    solution, evaluation
        The solution and evaluation for the solution
    """
    blank = env.blank()
    particles = np.stack([blank] * n_particles)
    velocities = np.stack([blank] * n_particles)
    evals = np.zeros(n_particles)
    n_recycle = int(n_particles * recycle)

    particle_range = np.arange(0, n_particles)
    for i in particle_range:
        particles[i] = env.action_space.sample()
        _, evals[i], *_ = env.step(particles[i])
        env.reset()

    best = np.argmin(evals)
    best_eval = evals[best]
    best_solution = particles[best]

    steps = 0
    while steps < max_steps:
        for i in particle_range:
            velocities[i] += np.random.normal(0, variation, size=blank.shape)
            velocities[i] *= 0.9
            particles[i] += velocities[i]
            _, evals[i], *_ = env.step(particles[i])
            env.reset()
            if evals[i] < best_eval:
                best_eval = evals[i]
                best_solution = particles[i]

        if steps % recycle_every == 0:
            # remove the worst particles
            # and replace them with the best particle+noise
            particles[np.argsort(evals)[n_recycle:]] = best_solution + np.random.normal(
                0, variation, size=blank.shape
            )

        if stop_if_done and env.is_done():
            break

        if stop_if_done and np.var(evals) < threshold:
            break

        steps += 1

    return best_solution, best_eval


# -------------------------------------------------------------------------------------
# This is the original code for swarm_optimize using pyswarms.
# This used to work really well with the old implementations of Rotatron.
# However, it seems to be not adapting well to the new Rotatron so I decided
# to switch to the barebone implementation in numpy alone above.
# -------------------------------------------------------------------------------------
# import pyswarms as ps
# def swarm_optimize(
#     env,
#     steps: int = 30,
#     n: int = 10,
#     molecule: "Molecule" = None,
#     bounds=(-np.pi, np.pi),
#     **kws,
# ):
#     """
#     Optimize a Rotatron environment through
#     a pyswarms swarm optimization

#     Parameters
#     ----------
#     env : biobuild.optimizers.environments.Rotatron
#         The environment to optimize
#     steps : int, optional
#         The number of steps to take, by default 1000
#     n : int, optional
#         The number of particles to use, by default 10
#     molecule : Molecule, optional
#         A molecule to apply the optimized solution to directly.
#         If None, the solution is returned, by default None. The solution is
#         applied in-place to the molecule directly.
#     bounds : tuple, optional
#         The bounds to use for solutions as a tuple of size 2 with a minimal and maximal angle to allow.
#     kws : dict, optional
#         Keyword arguments to pass as options to the optimizer

#     Returns
#     -------
#     tuple or Molecule
#         If molecule is None, a tuple of the optimized action and the reward for the optimized action.
#         If molecule is not None, the optimized molecule.
#     """

#     x0 = env.action_space.sample()

#     if bounds:
#         bounds = (np.full(x0.shape, bounds[0]), np.full(x0.shape, bounds[1]))

#     def loss_fn(sol):
#         costs = np.zeros(n)
#         i = 0
#         for s in sol:
#             costs[i] = env.step(s)[1]
#             i += 1
#         return costs

#     options = kws.pop("options", {})
#     options.setdefault("c1", 0.5)
#     options.setdefault("c2", 0.3)
#     options.setdefault("w", 0.9)

#     verbose = kws.pop("verbose", False)

#     optimizer = ps.single.GlobalBestPSO(
#         n_particles=n, dimensions=x0.shape[0], bounds=bounds, options=options, **kws
#     )
#     reward, result = optimizer.optimize(loss_fn, iters=steps, verbose=verbose)
#     if molecule:
#         molecule = outils.apply_solution(result, env, molecule)
#         return molecule
#     return result, reward

__all__ = [
    "swarm_optimize",
    "scipy_optimize",
    "genetic_optimize",
]
