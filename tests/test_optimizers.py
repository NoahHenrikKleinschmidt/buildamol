"""
Tests for the optimizers
"""

import tests.base as base
import numpy as np
import pytest
import buildamol as bam
import buildamol.optimizers as opt


def test_distance_rotatron():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_distance_rotatron_resgraph():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_overlap_rotatron():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_overlap_rotatron_resgraph():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_forcefield_rotatron():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_forcefield_rotatron_resgraph():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.swarm
def test_swarm():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.swarm_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None


@pytest.mark.slow
def test_genetic():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.genetic_optimize(env, max_generations=3)
    assert sol is not None
    assert _eval is not None


@pytest.mark.slow
def test_anneal():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.anneal_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None


# =====================================================================================
# BACKEND-SPECIFIC TESTS
# =====================================================================================


@pytest.mark.parametrize("backend", ["numpy", "numba", "jax"])
@pytest.mark.swarm
def test_swarm_explicit_backend(backend):
    """Test swarm optimization with explicit backend parameter."""
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.swarm_optimize(env, n_particles=3, max_steps=3, backend=backend)
    assert sol is not None
    assert _eval is not None
    assert isinstance(sol, np.ndarray)
    assert isinstance(_eval, (float, np.floating))


@pytest.mark.parametrize("backend", ["numpy", "numba", "jax"])
@pytest.mark.slow
def test_anneal_explicit_backend(backend):
    """Test anneal optimization with explicit backend parameter."""
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.anneal_optimize(env, n_particles=3, max_steps=3, backend=backend)
    assert sol is not None
    assert _eval is not None
    assert isinstance(sol, np.ndarray)
    assert isinstance(_eval, (float, np.floating))


@pytest.mark.parametrize("backend", ["numpy", "numba", "jax"])
@pytest.mark.slow
def test_genetic_explicit_backend(backend):
    """Test genetic optimization with explicit backend parameter."""
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, backend=backend)
    assert env is not None

    sol, _eval = opt.genetic_optimize(env, max_generations=3, backend=backend)
    assert sol is not None
    assert _eval is not None
    assert isinstance(sol, np.ndarray)
    assert isinstance(_eval, (float, np.floating))


@pytest.mark.slow
def test_scipy():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.scipy_optimize(env)
    assert sol is not None
    assert _eval is not None


def test_apply():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.swarm_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None

    out = opt.apply_rotatron_solution(sol, env, mol.copy())
    assert out is not None

    assert out.count_atoms() == mol.count_atoms()
    assert out.count_bonds() == mol.count_bonds()

    before = np.array([i.coord for i in mol.get_atoms()])
    after = np.array([i.coord for i in out.get_atoms()])

    assert not np.allclose(before, after, atol=1e-2)

    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.swarm
def test_optim_distance_swarm():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "swarm", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_distance_anneal():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "anneal", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_distance_genetic():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "genetic", max_generations=500)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_distance_scipy():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "scipy")

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.swarm
def test_optim_overlap_swarm():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.OverlapRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "swarm", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_overlap_anneal():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.OverlapRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "anneal", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_overlap_genetic():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.OverlapRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "genetic", max_generations=500)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_overlap_scipy():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.OverlapRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "scipy")

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.swarm
def test_optim_forcefield_swarm():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "swarm", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_forcefield_anneal():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "anneal", n_particles=50)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_forcefield_genetic():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "genetic", max_generations=500)

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.optim
@pytest.mark.slow
def test_optim_forcefield_scipy():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None

    out = opt.optimize(mol.copy(), env, "scipy")

    print(out.count_clashes())
    if base.ALLOW_VISUAL:
        out.show()


@pytest.mark.numba
def test_use_numba_setting():
    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA

    bam.use_numba()

    assert bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA

    bam.dont_use_numba()

    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA

    bam.use_all_numba()

    assert not bam.utils.auxiliary.USE_NUMBA
    assert bam.utils.auxiliary.USE_ALL_NUMBA

    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.swarm
def test_optim_numba_distance_swarm():
    bam.use_all_numba()
    test_optim_distance_swarm()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_anneal():
    bam.use_all_numba()
    test_optim_distance_anneal()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_genetic():
    bam.use_all_numba()
    test_optim_distance_genetic()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_scipy():
    bam.use_all_numba()
    test_optim_distance_scipy()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.swarm
def test_optim_numba_overlap_swarm():
    bam.use_all_numba()
    test_optim_overlap_swarm()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_anneal():
    bam.use_all_numba()
    test_optim_overlap_anneal()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_genetic():
    bam.use_all_numba()
    test_optim_overlap_genetic()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_scipy():
    bam.use_all_numba()
    test_optim_overlap_scipy()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.swarm
def test_optim_numba_forcefield_swarm():
    bam.use_all_numba()
    test_optim_forcefield_swarm()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_anneal():
    bam.use_all_numba()
    test_optim_forcefield_anneal()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_genetic():
    bam.use_all_numba()
    test_optim_forcefield_genetic()
    bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_scipy():
    bam.use_all_numba()
    test_optim_forcefield_scipy()
    bam.dont_use_numba()


# def test_benchmark_numba():
#     mol = bam.read_pdb("/Users/noahhk/GIT/biobuild/0.pdb")
#     g = mol.get_atom_graph()
#     edges = g.find_rotatable_edges(g.central_node)
#     edges = g.sample_edges(edges)

#     from time import time

#     f = "numba_benchmark.csv"
#     open(f, "w").close()

#     for mode in ("normal", "numba"):

#         env = opt.DistanceRotatron(g, edges, n_processes=5)

#         for i in range(3):
#             m = mol.copy()
#             start = time()
#             out = opt.optimize(m, env, algorithm="swarm", n_particles=50)
#             dt = time() - start
#             print(out.count_clashes(), dt, mode, sep=",", file=open(f, "a"))

#         bam.use_all_numba()


# def test_benchmark_numba_hybrid():
#     mol = bam.read_pdb("/Users/noahhk/GIT/biobuild/0.pdb")
#     g = mol.get_atom_graph(True)
#     edges = g.find_rotatable_edges(g.central_node)
#     edges = g.sample_edges(edges)

#     from time import time

#     f = "numba_benchmark.csv"
#     open(f, "w").close()

#     env = opt.DistanceRotatron(g, edges, n_processes=5, numba=True)

#     for i in range(10):
#         m = mol.copy()
#         start = time()
#         out = opt.optimize(m, env, algorithm="swarm", max_steps=100, n_particles=50)
#         dt = time() - start
#         print(out.count_clashes(), dt, "hybrid", sep=",", file=open(f, "a"))


# #         start = time()
# #         sol, _eval = opt.swarm_optimize(env2, n_particles=50)
# #         dt = time() - start
# #         print(dt, ", normal", file=open("rotate_benchmark_numba.csv", "a"))


def test_translatron_optimize():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)

    ref_coords = mol.get_coords()

    def constraint(env, coords):
        dist = (coords - ref_coords) ** 2
        dist = np.mean(dist)
        env._dist = dist
        return dist

    def finish(env, coords):
        return env._dist < 0.1

    mol.move([10, 10, 10])
    mol.rotate(34, [0.6, 0.12, 0.5])
    mol.rotate(90, [0.1, 0.2, 0.3])

    env = opt.Translatron(mol.get_atom_graph(), constraint, finish)

    out = opt.optimize(mol.copy(), env)

    assert out.count_clashes() == mol.count_clashes()

    if base.ALLOW_VISUAL:
        v = mol.draw(atoms=False)
        v += out.draw(atoms=False, line_color="red")
        v.show()


def test_constraint_rotatron_step():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    base_env = opt.DistanceRotatron(g, edges)

    def constraint(rotatron, state, **kwargs):
        ref = kwargs["ref"]
        coords = np.asarray(state)
        return float(np.mean((coords - ref) ** 2))

    def finisher(rotatron, state, **kwargs):
        return constraint(rotatron, state, **kwargs) < 1e-2

    wrapped = opt.ConstraintRotatron(
        base_env,
        constraint,
        finisher,
        ref=base_env._backup_state.copy(),
    )

    action = wrapped.action_space.sample()
    state, score, done, info = wrapped.step(action)
    assert state is not None
    assert np.isfinite(score)
    assert isinstance(done, bool)
    assert isinstance(info, dict)
    wrapped.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_constraint_rotatron_step_jax_wrapped_env():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    base_env = opt.DistanceRotatron(g, edges, backend="jax")

    def constraint(rotatron, state, **kwargs):
        ref = kwargs["ref"]
        coords = np.asarray(state)
        return float(np.mean((coords - ref) ** 2))

    wrapped = opt.ConstraintRotatron(
        base_env,
        constraint,
        ref=base_env._backup_state.copy(),
    )

    action = wrapped.action_space.sample()
    state, score, done, info = wrapped.step(action)
    assert state is not None
    assert np.isfinite(score)
    assert isinstance(done, bool)
    assert isinstance(info, dict)
    wrapped.reset()


def test_apply_inplace():
    mol = bam.read_smiles("C1=CC(=O)C(CCCOC)CCC1")
    mol.autolabel()
    mol.optimize()
    if base.ALLOW_VISUAL:
        mol.show()


# ==== JAX TESTS ====


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_distance_rotatron_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_distance_rotatron_resgraph_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_overlap_rotatron_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_overlap_rotatron_resgraph_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_forcefield_rotatron_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


@pytest.mark.skipif(not bam.utils.auxiliary.HAS_JAX, reason="JAX not available")
def test_forcefield_rotatron_resgraph_jax():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges, backend="jax")
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()
