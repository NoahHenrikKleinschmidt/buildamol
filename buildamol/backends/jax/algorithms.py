"""
JAX backend for optimization algorithms.

This module provides JAX implementations of key optimization algorithms.
JAX enables functional transformations and preserves JAX arrays throughout computation.

Functions in this module override the NumPy defaults when backend="jax" is selected.
"""

import numpy as np
import buildamol.utils.auxiliary as aux


__all__ = [
    "swarm_optimize",
    "anneal_optimize",
    "genetic_optimize",
]


def swarm_optimize(
    env,
    n_particles: int = None,
    max_steps: int = 30,
    stop_if_done: bool = True,
    threshold: float = 1e-6,
    w: float = 0.9,
    c1: float = 0.5,
    c2: float = 0.3,
    cooldown_rate: float = 0.99,
    n_best: int = 1,
    backend: str = None,
):
    """
    JAX-native particle swarm optimization.

    Preserves JAX arrays throughout computation for compatibility with JAX transformations.

    Parameters
    ----------
    env : buildamol.optimizers.environments.Rotatron
        The environment to optimize
    n_particles : int, optional
        The number of particles to use.
    max_steps : int, optional
        The maximum number of steps to take.
    stop_if_done : bool, optional
        Stop if environment signals done or if solutions converge.
    threshold : float, optional
        Convergence threshold for fitness variance.
    w : float, optional
        Inertia parameter.
    c1 : float, optional
        Cognitive parameter.
    c2 : float, optional
        Social parameter.
    cooldown_rate : float, optional
        Rate at which inertia parameter decays.
    n_best : int, optional
        The number of best solutions to return.
    backend : str, optional
        Backend to use (ignored, always JAX).

    Returns
    -------
    solution, evaluation
        The best solution(s) and their evaluation(s).
    """
    jnp = aux.get_jax_numpy()

    if n_particles is None:
        n_particles = max(15, len(env.rotatable_edges) // 3)

    positions = jnp.array([env.action_space.sample() for i in range(n_particles)])
    velocities = jnp.array(np.random.rand(n_particles, env.action_space.shape[0]))
    best_positions = positions.copy()
    best_fitnesses = jnp.full(n_particles, jnp.inf)
    fitnesses = jnp.zeros(n_particles)

    bounds = getattr(env, "_bounds_tuple", None) or (-9999, 9999)
    bounds = jnp.array(bounds, dtype=jnp.float64)

    best_fitness = np.inf
    best_solution = np.zeros(env.action_space.shape[0])

    steps = 0
    while steps < max_steps:
        for i in range(n_particles):
            fitnesses = fitnesses.at[i].set(env.step(np.asarray(positions[i]))[1])
            if fitnesses[i] < best_fitnesses[i]:
                best_fitnesses = best_fitnesses.at[i].set(fitnesses[i])
                best_positions = best_positions.at[i].set(positions[i])
            if float(fitnesses[i]) < best_fitness:
                best_fitness = float(fitnesses[i])
                best_solution[:] = np.asarray(positions[i])
            env.reset()

        velocities = _jax_update_velocities(
            velocities, w, c1, c2, best_positions, best_solution, positions
        )
        positions = _jax_update_positions(positions, velocities, bounds)

        if stop_if_done:
            if float(jnp.var(best_fitnesses)) < threshold:
                break

        w *= cooldown_rate
        steps += 1

    if n_best == 1:
        return best_solution, best_fitness
    else:
        best = jnp.argsort(best_fitnesses)[:n_best]
        return np.array([np.asarray(best_positions[b]) for b in best]), np.array(
            [float(best_fitnesses[b]) for b in best]
        )


def anneal_optimize(
    env,
    n_particles: int = None,
    max_steps: int = 100,
    stop_if_done: bool = True,
    threshold: float = 1e-6,
    variance: float = 0.3,
    cooldown_rate: float = 0.98,
    n_best: int = 1,
    backend: str = None,
):
    """
    JAX-native simulated annealing optimization.

    Preserves JAX arrays throughout computation for compatibility with JAX transformations.

    Parameters
    ----------
    env : buildamol.optimizers.environments.Rotatron
        The environment to optimize
    n_particles : int, optional
        The number of particles to use.
    max_steps : int, optional
        The maximum number of steps to take.
    stop_if_done : bool, optional
        Stop the optimization if the environment signals it is done or convergence achieved.
    threshold : float, optional
        Convergence threshold for fitness variance.
    variance : float, optional
        The variation to use for updating particle positions.
    cooldown_rate : float, optional
        Rate at which temperature and variance decay.
    n_best : int, optional
        The number of best solutions to return.
    backend : str, optional
        Backend to use (ignored, always JAX).

    Returns
    -------
    solution, evaluation
        The solution(s) and evaluation(s).
    """
    jnp = aux.get_jax_numpy()

    if n_particles is None:
        n_particles = max(15, len(env.rotatable_edges) // 2)

    particles = jnp.array([env.action_space.sample() for i in range(n_particles)])
    fitnesses = jnp.full(n_particles, 9999.0)
    best_fitness = np.inf
    best_solution = np.zeros(env.action_space.shape[0])

    bounds = getattr(env, "_bounds_tuple", None) or (-9999, 9999)

    temperature = 1.0

    steps = 0
    while steps < max_steps:
        for i in range(n_particles):
            position = np.asarray(particles[i]) + np.random.uniform(
                -variance, variance, size=np.asarray(particles[i]).shape
            )
            position = np.clip(position, bounds[0], bounds[1])
            fitness = env.step(position)[1]

            if _jax_accept(float(fitnesses[i]), fitness, temperature):
                fitnesses = fitnesses.at[i].set(fitness)
                particles = particles.at[i].set(jnp.array(position))

                if fitness < best_fitness:
                    best_fitness = fitness
                    best_solution[:] = position

            env.reset()

        temperature *= cooldown_rate
        variance *= cooldown_rate * 0.95
        steps += 1

        if stop_if_done:
            if float(jnp.var(fitnesses)) < threshold:
                break

    if n_best == 1:
        return best_solution, best_fitness
    else:
        best = jnp.argsort(fitnesses)[:n_best]
        return np.array([np.asarray(particles[b]) for b in best]), np.array(
            [float(fitnesses[b]) for b in best]
        )


def genetic_optimize(
    env,
    max_generations: int = 5e2,
    stop_if_done: bool = True,
    threshold: float = 1e-6,
    variation: float = 0.2,
    population_size: int = 50,
    parents: int = 0.25,
    children: int = 0.3,
    mutants: int = 0.3,
    newcomers: int = 0.15,
    variation_cooldown: float = 1,
    n_best: int = 1,
    backend: str = None,
):
    """
    JAX-native genetic algorithm optimization.

    Preserves JAX arrays throughout computation for compatibility with JAX transformations.

    Parameters
    ----------
    env : buildamol.optimizers.environments.Rotatron
        The environment to optimize
    max_generations : int, optional
        The maximum number of generations.
    stop_if_done : bool, optional
        Stop if environment signals done or convergence achieved.
    threshold : float, optional
        Convergence threshold for fitness variance.
    variation : float, optional
        Initial variation for mutation.
    population_size : int, optional
        The size of the population.
    parents : int or float, optional
        The number or fraction of parents to select.
    children : int or float, optional
        The number or fraction of children to generate.
    mutants : int or float, optional
        The number or fraction of mutants to generate.
    newcomers : int or float, optional
        The number or fraction of new random solutions.
    variation_cooldown : float, optional
        Rate at which variation decreases per generation.
    n_best : int, optional
        The number of best solutions to return.
    backend : str, optional
        Backend to use (ignored, always JAX).

    Returns
    -------
    solution, evaluation
        The best solution(s) and their evaluation(s).
    """
    jnp = aux.get_jax_numpy()

    max_generations = int(max_generations)

    # Convert fractions to counts
    if isinstance(parents, float):
        parents = round(parents * population_size)
    if isinstance(children, float):
        children = round(children * population_size)
    if isinstance(mutants, float):
        mutants = round(mutants * population_size)
    if isinstance(newcomers, float):
        newcomers = round(newcomers * population_size)

    # Adjust counts to match population size
    while children + mutants + parents + newcomers > population_size:
        if newcomers > 1:
            newcomers -= 1
        elif children > 1:
            children -= 1
        elif mutants > 1:
            mutants -= 1
        elif parents > 1:
            parents -= 1
    while children + mutants + parents + newcomers < population_size:
        mutants += 1

    if any([n <= 0 for n in [children, mutants, parents, newcomers]]):
        raise ValueError(
            "n_children, n_mutations, n_parents, and n_newcomers must all be at least 1"
        )

    n_parents = parents
    n_children = children
    n_mutants = mutants
    n_newcomers = newcomers

    blank = env.blank()
    if hasattr(env, "_bounds_tuple"):
        min_angle, max_angle = env._bounds_tuple
    else:
        min_angle, max_angle = -np.pi, np.pi

    population = jnp.array(np.stack([blank] * population_size))
    evals = jnp.zeros(population_size)

    bests = jnp.zeros(max(10, int(max_generations * 0.3)))

    pop_range = np.arange(0, population_size)
    parents_range = np.arange(0, n_parents)
    children_range = np.arange(0, n_children)
    mutations_range = np.arange(0, n_mutants)
    newcomers_range = np.arange(0, n_newcomers)

    for i in pop_range:
        population = population.at[i].set(jnp.array(env.action_space.sample()))
        _, evals = evals.at[i].set(env.step(np.asarray(population[i]))[1])
        env.reset()

    steps = 0
    while steps < max_generations:
        parents_pop, children_pop, mutations_pop = _jax_make_generation(
            evals,
            population,
            n_parents,
            n_children,
            n_mutants,
            variation,
            blank,
            parents_range,
            children_range,
            mutations_range,
        )

        newcomers_pop = jnp.array(np.stack([blank] * n_newcomers))
        for i in newcomers_range:
            newcomers_pop = newcomers_pop.at[i].set(
                jnp.array(env.action_space.sample())
            )

        population = jnp.concatenate(
            [parents_pop, children_pop, mutations_pop, newcomers_pop]
        )
        population = jnp.clip(population, min_angle, max_angle)

        evals = jnp.zeros(population_size)

        done = False
        for i in pop_range:
            _, eval_val, done, *_ = env.step(np.asarray(population[i]))
            evals = evals.at[i].set(eval_val)
            env.reset()
            if done and stop_if_done:
                break

        if done and stop_if_done:
            break

        bests, should_break = _jax_converged_break(
            bests, evals, threshold, steps, max_generations, stop_if_done
        )
        if should_break:
            break

        variation *= variation_cooldown
        steps += 1

    best = jnp.argsort(evals)[:n_best]
    if n_best == 1:
        return np.asarray(population[best[0]]), float(evals[best[0]])
    return np.array([np.asarray(population[b]) for b in best]), np.array(
        [float(evals[b]) for b in best]
    )


# =====================================================================================
# HELPER FUNCTIONS
# =====================================================================================


def _jax_update_velocities(
    velocities, inertia, cognitive, social, best_particle, best_position, positions
):
    """Update particle velocities using JAX arrays."""
    jnp = aux.get_jax_numpy()
    dist_best = best_particle - positions
    dist_social = best_position - positions
    velocities = (
        inertia * velocities
        + cognitive * np.random.rand() * dist_best
        + social * np.random.rand() * dist_social
    )
    return velocities


def _jax_update_positions(positions, velocities, bounds):
    """Update particle positions using JAX arrays."""
    jnp = aux.get_jax_numpy()
    positions = jnp.clip(positions + velocities, bounds[0], bounds[1])
    return positions


def _jax_accept(old, new, temp):
    """Accept or reject new position based on simulated annealing criterion."""
    jnp = aux.get_jax_numpy()
    return float(jnp.exp((old - new) / temp)) > np.random.rand()


def _jax_make_generation(
    evals,
    population,
    n_parents,
    n_children,
    n_mutants,
    variation,
    blank,
    parents_range,
    children_range,
    mutations_range,
):
    """Create new generation using genetic algorithm operations."""
    jnp = aux.get_jax_numpy()

    sorting = jnp.argsort(evals)
    parents = population[sorting[:n_parents]]

    children = jnp.zeros((n_children, blank.shape[0]))
    for i in children_range:
        p1, p2 = np.random.choice(parents_range, size=2, replace=False)
        p1, p2 = parents[int(p1)], parents[int(p2)]
        children = children.at[i].set(
            (p1 + p2) / 2
            + jnp.array(
                np.random.uniform(-variation / 3, variation / 3, size=blank.shape)
            )
        )

    mutations = jnp.zeros((n_mutants, blank.shape[0]))
    for i in mutations_range:
        mutations = mutations.at[i].set(
            parents[np.random.randint(0, n_parents)]
            + jnp.array(np.random.uniform(-variation, variation, size=blank.shape))
        )

    return parents, children, mutations


def _jax_converged_break(bests, evals, threshold, steps, max_generations, stop_if_done):
    """Check convergence and decide whether to break."""
    jnp = aux.get_jax_numpy()
    bests = jnp.roll(bests, -1)
    bests = bests.at[-1].set(jnp.min(evals))
    if (
        stop_if_done
        and steps > max_generations / 10
        and float(jnp.var(bests[bests != 0])) < threshold
    ):
        return bests, True
    return bests, False
