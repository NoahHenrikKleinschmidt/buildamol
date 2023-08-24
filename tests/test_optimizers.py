"""
Tests for the optimizers
"""

import tests.base as bas
import numpy as np
import biobuild as bb
import biobuild.optimizers as opt


def test_rotatron_constraints():
    mol = bb.Molecule.from_pdb(base.MANNOSE9)
    mol.infer_bonds(restrict_residues=False)

    connections = mol.get_residue_connections()

    env = opt.MultiBondRotatron(mol.make_residue_graph(detailed=True), connections)

    initial_coords = env._coords.copy()
    nodes = list(env.graph.nodes)

    initial_distances = []
    for node in nodes:
        neighbors = env.graph.get_neighbors(node)
        dists = np.array(
            [np.linalg.norm(node.coord - neighbor.coord) for neighbor in neighbors]
        )
        initial_distances.append(dists)

    for run in range(3):
        v = bb.utils.visual.MoleculeViewer3D(mol.make_residue_graph(detailed=True))

        for i in range(200):
            action = env.action_space.sample()
            descendant_mask = env._descendant_masks[action[0]]
            descendant_coords = env._coords[descendant_mask].copy()
            ancestor_coords = env._coords[~descendant_mask].copy()

            env.step(action)

            assert not np.allclose(descendant_coords, env._coords[descendant_mask])
            assert np.allclose(ancestor_coords, env._coords[~descendant_mask])

        env.apply_to_graph()
        new_distances = []
        for node in nodes:
            dists = []
            neighbors = env.graph.get_neighbors(node)
            for neighbor in neighbors:
                dists.append(np.linalg.norm(node.coord - neighbor.coord))
            new_distances.append(np.asarray(dists))

        v.draw_edges(env.graph.edges, color="cyan", linewidth=2)
        v.show()

        for i in range(len(nodes)):
            assert np.all(np.abs(initial_distances[i] - new_distances[i]) < 1e-3)

        env.reset()
        assert np.allclose(initial_coords, env.state)
