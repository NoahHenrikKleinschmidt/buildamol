import biobuild as bb

bb.load_small_molecules()

benzene = bb.molecule("benzene")

# -----------------------------
#     make the periphery
# -----------------------------
periphery = benzene.copy()

# set up the linkage instructions
# always shifting the carbon at which to attach
link = bb.linkage("C1", "C1")
for carbon in range(1, 6):
    link.atom1 = f"C{carbon}"
    periphery.attach(benzene, link, at_residue=1)

# -----------------------------
#     assemble the molecule
# -----------------------------
mol = benzene.copy()

link2 = bb.linkage("C1", "C4")
periphery.set_attach_residue(2)

# and attach the periphery to the core
for carbon in core.get_atoms("C", by="element"):
    link2.atom1 = carbon.id
    mol.attach(periphery, link2, at_residue=1)

# -----------------------------
#     optimize the molecule
# -----------------------------
graph = core.get_atom_graph()
edges = core.get_residue_connections()
edges = graph.direct_edges(graph.central_node, edges)

env = bb.optimizers.DistanceRotatron(graph, edges)

mol_opt = bb.optimizers.optimize(mol, env, algorithm="swarm")

mol_opt.to_pdb("polyphenylene.pdb")
