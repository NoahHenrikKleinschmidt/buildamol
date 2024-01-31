"""
Tests for the optimizers
"""

import tests.base as base
import numpy as np
import buildamol as bam
import buildamol.optimizers as opt


def test_distance_rotatron():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_distance_rotatron_resgraph():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges, unfold=3, pushback=10)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_overlap_rotatron():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_overlap_rotatron_resgraph():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.OverlapRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_forcefield_rotatron():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_forcefield_rotatron_resgraph():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_residue_graph(detailed=True)
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.ForceFieldRotatron(g, edges)
    assert env is not None
    env.step(env.action_space.sample())
    env.reset()


def test_swarm():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.swarm_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None


def test_genetic():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.genetic_optimize(env, max_generations=3)
    assert sol is not None
    assert _eval is not None


def test_anneal():
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)
    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.anneal_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None


def test_scipy():
    mol = bam.read_pdb(base.MANNOSE9)
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
    mol = bam.read_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)
    assert mol is not None

    g = mol.get_atom_graph()
    edges = g.find_rotatable_edges(g.central_node)

    env = opt.DistanceRotatron(g, edges)
    assert env is not None

    sol, _eval = opt.swarm_optimize(env, n_particles=3, max_steps=3)
    assert sol is not None
    assert _eval is not None

    out = opt.apply_solution(sol, env, mol.copy())
    assert out is not None

    assert out.count_atoms() == mol.count_atoms()
    assert out.count_bonds() == mol.count_bonds()

    before = np.array([i.coord for i in mol.get_atoms()])
    after = np.array([i.coord for i in out.get_atoms()])

    assert not np.allclose(before, after, atol=1e-2)
