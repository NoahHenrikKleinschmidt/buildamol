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
    bam.dont_use_jax()
    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA
    assert not bam.utils.auxiliary.USE_JAX

    bam.use_numba()

    assert bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA
    assert not bam.utils.auxiliary.USE_JAX

    bam.dont_use_numba()

    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA
    assert not bam.utils.auxiliary.USE_JAX

    bam.use_all_numba()

    assert not bam.utils.auxiliary.USE_NUMBA
    assert bam.utils.auxiliary.USE_ALL_NUMBA
    assert not bam.utils.auxiliary.USE_JAX

    bam.dont_use_numba()


def test_use_jax_setting():
    bam.dont_use_numba()
    bam.dont_use_jax()

    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA
    assert not bam.utils.auxiliary.USE_JAX

    bam.use_jax()
    if bam.utils.auxiliary.HAS_JAX:
        assert bam.utils.auxiliary.USE_JAX
    else:
        assert not bam.utils.auxiliary.USE_JAX
    assert not bam.utils.auxiliary.USE_NUMBA
    assert not bam.utils.auxiliary.USE_ALL_NUMBA

    bam.use_numba()
    assert not bam.utils.auxiliary.USE_JAX

    bam.dont_use_numba()
    bam.dont_use_jax()


def test_optim_swarm_jax_fallback():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, eval = opt.swarm_optimize(env, max_steps=5, n_particles=5, backend="jax")
    assert sol is not None
    assert np.isfinite(eval)


def test_optim_numba_distance_swarm():
    bam.use_all_numba()
    try:
        test_optim_distance_swarm()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_anneal():
    bam.use_all_numba()
    try:
        test_optim_distance_anneal()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_genetic():
    bam.use_all_numba()
    try:
        test_optim_distance_genetic()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_distance_scipy():
    bam.use_all_numba()
    try:
        test_optim_distance_scipy()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.swarm
def test_optim_numba_overlap_swarm():
    bam.use_all_numba()
    try:
        test_optim_overlap_swarm()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_anneal():
    bam.use_all_numba()
    try:
        test_optim_overlap_anneal()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_genetic():
    bam.use_all_numba()
    try:
        test_optim_overlap_genetic()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_overlap_scipy():
    bam.use_all_numba()
    try:
        test_optim_overlap_scipy()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.swarm
def test_optim_numba_forcefield_swarm():
    bam.use_all_numba()
    try:
        test_optim_forcefield_swarm()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_anneal():
    bam.use_all_numba()
    try:
        test_optim_forcefield_anneal()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_genetic():
    bam.use_all_numba()
    try:
        test_optim_forcefield_genetic()
    finally:
        bam.dont_use_numba()


@pytest.mark.numba
@pytest.mark.optim
@pytest.mark.slow
def test_optim_numba_forcefield_scipy():
    bam.use_all_numba()
    try:
        test_optim_forcefield_scipy()
    finally:
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


def _scrambled_man9(seed=0):
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node, min_descendants=3)
    env = opt.DistanceRotatron(g, edges)
    angles = np.random.default_rng(seed).uniform(-np.pi, np.pi, env.n_edges)
    return opt.apply_rotatron_solution(angles, env, mol)


def test_evolution_never_worse():
    mol = _scrambled_man9()
    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node, min_descendants=3)
    env = opt.DistanceRotatron(g, edges)

    start = env.step(env.blank())[1]
    env.reset()

    sol, evaluation = opt.evolution_optimize(env, max_generations=10, seed=1)
    assert sol.shape == (env.n_edges,)
    assert evaluation <= start
    assert np.isclose(env.step(sol)[1], evaluation)
    env.reset()

    # with blocks only the indices in the blocks are changed
    blocks = [list(range(0, 5)), list(range(5, 10))]
    sol, evaluation = opt.evolution_optimize(
        env, blocks=blocks, sweeps=2, max_generations=5, seed=1
    )
    assert evaluation <= start
    assert np.all(sol[10:] == 0)

    sols, evals = opt.evolution_optimize(env, max_generations=5, n_best=3, seed=1)
    assert sols.shape == (3, env.n_edges)
    assert np.all(np.diff(evals) >= 0)


def test_evolution_optimize_function():
    mol = _scrambled_man9()
    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node, min_descendants=3)
    env = opt.DistanceRotatron(g, edges)
    clashes = mol.count_clashes()
    out = opt.optimize(mol.copy(), env, "evolution", max_generations=20, seed=1)
    assert out.count_clashes() <= clashes


def test_optimize_with_callable_algorithm():
    mol = _scrambled_man9()
    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node, min_descendants=3)
    env = opt.DistanceRotatron(g, edges)
    out = opt.optimize(mol.copy(), env, lambda env: (env.blank(), 0.0))
    assert np.allclose(out.get_coords(), mol.get_coords())


def test_translatron_solution_matches_environment():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    for g in (mol.get_atom_graph(), mol.get_residue_graph(detailed=True)):
        env = opt.Translatron(g, lambda env, coords: 0.0)
        sol = np.array([1.0, 2.0, 3.0, 0.4, -0.3, 0.8])
        state = env.step(sol)[0].copy()
        env.reset()

        out = opt.apply_translatron_solution(sol, env, mol.copy())
        atoms = [
            (i, out.get_atom(n.serial_number))
            for i, n in enumerate(g.nodes)
            if isinstance(n, bam.Atom)
        ]
        coords = np.array([a.coord for _, a in atoms])
        assert np.allclose(state[[i for i, _ in atoms]], coords)


def test_forcefield_rotatron_energy_matches_rdkit():
    from rdkit.Chem import AllChem

    def energy(mol):
        rdmol = mol.to_rdkit()
        props = AllChem.MMFFGetMoleculeProperties(rdmol)
        return AllChem.MMFFGetMoleculeForceField(rdmol, props).CalcEnergy()

    bam.load_sugars()
    mol = (bam.molecule("GLC") % "14bb") * 3
    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node, min_descendants=3)[:4]
    env = opt.ForceFieldRotatron(g, edges)

    rng = np.random.default_rng(1)
    # the evaluation must not depend on previously evaluated states
    for angles in (rng.uniform(-np.pi, np.pi, env.n_edges), env.blank()):
        evaluation = env.step(angles)[1]
        env.reset()
        out = opt.apply_rotatron_solution(angles, env, mol.copy())
        assert np.isclose(evaluation, energy(out))


def test_crop_does_not_modify_graph():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    g = mol.get_atom_graph()
    n_nodes = len(g.nodes)
    edges = g.find_rotatable_edges(g.central_node, min_descendants=5)[:3]

    env = opt.DistanceRotatron(g, edges, radius=5, crop_nodes_further_than=1.0)
    assert env.n_nodes < n_nodes
    assert len(g.nodes) == n_nodes

    # no edges given
    env = opt.DistanceRotatron(g, None, radius=5, crop_nodes_further_than=1.0)
    assert len(g.nodes) == n_nodes


def test_circulatron_and_constraint_rotatron_is_done():
    mol = bam.read_pdb(base.MAN9PDB)
    mol.infer_bonds(restrict_residues=False)
    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node, min_descendants=5)[:4]

    base_env = opt.DistanceRotatron(g, edges)
    env = opt.ConstraintRotatron(
        base_env, lambda rotatron, state: 0.0, lambda rotatron, state: True
    )
    assert env.n_edges == base_env.n_edges
    assert env.graph is base_env.graph
    assert env.state is base_env.state
    assert env.is_done(env.state) == base_env.is_done(env.state)

    nodes = list(g.nodes)
    env = opt.Circulatron(g, (nodes[0], nodes[0]), rotatable_edges=edges)
    # the target nodes are identical, so the environment is done immediately
    assert env.step(env.blank())[2]


def test_kdtree_rotatron():
    mol = _scrambled_man9()
    for graph in (mol.get_atom_graph(), mol.get_residue_graph(detailed=True)):
        edges = graph.find_rotatable_edges(graph.central_node, min_descendants=3)
        env = opt.KDTreeRotatron(graph, edges)
        # clashes are always evaluated between atoms
        assert all(isinstance(n, bam.Atom) for n in env.graph.nodes)

        actions = np.random.default_rng(0).uniform(-np.pi, np.pi, (5, env.n_edges))
        actions[0] = 0
        stepwise = []
        for action in actions:
            stepwise.append(env.step(action)[1])
            env.reset()
        batch, done = env.eval_batch(actions)
        assert np.allclose(stepwise, batch)
        assert np.all(done == (batch <= 0))

        # the evaluation of an applied solution matches the environment
        out = opt.apply_rotatron_solution(actions[2], env, mol.copy())
        out_graph = out.get_atom_graph()
        out_edges = [
            tuple(out.get_atom(a.serial_number) for a in edge)
            for edge in env.rotatable_edges
        ]
        out_env = opt.KDTreeRotatron(out_graph, out_edges)
        assert np.isclose(out_env.step(out_env.blank())[1], batch[2])


def test_kdtree_rotatron_ignores_fixed_pairs():
    import networkx as nx

    mol = _scrambled_man9()
    graph = mol.get_atom_graph()
    edges = graph.find_rotatable_edges(graph.central_node, min_descendants=3)
    env = opt.KDTreeRotatron(graph, edges)
    evaluation = env.step(env.blank())[1]

    # brute force: only heavy atom pairs that can change their distance
    # (different rotation units and more than two bonds apart) are penalized
    nodes = list(graph.nodes)
    heavy = [i for i, n in enumerate(nodes) if n.element != "H"]
    units = {
        i: tuple(env.edge_masks[:, i]) for i in heavy
    }
    expected = 0
    for x, i in enumerate(heavy):
        for j in heavy[x + 1 :]:
            d = np.linalg.norm(env.state[i] - env.state[j])
            if d >= env.clash_distance or units[i] == units[j]:
                continue
            if nx.shortest_path_length(graph, nodes[i], nodes[j]) <= 2:
                continue
            expected += (env.clash_distance - d) ** 2
    assert evaluation > 0
    assert np.isclose(evaluation, expected)


def test_optimize_default_environment():
    mol = _scrambled_man9()
    clashes = mol.count_clashes()
    out = opt.optimize(mol.copy())
    assert out.count_clashes() < clashes
