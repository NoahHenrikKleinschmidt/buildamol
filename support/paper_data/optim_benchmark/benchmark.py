"""
This is for the paper benchmark optimization
"""

import buildamol as bam
from time import time
from pathlib import Path

from argparse import ArgumentParser

if __name__ == "__main__":

    parser = ArgumentParser()
    parser.add_argument("--mol", type=str, help="molecule input")
    parser.add_argument("--repeats", type=int, help="number of repeats")
    parser.add_argument("--reset", action="store_true", help="reset the output file")
    parser.add_argument(
        "--outdir", type=str, help="output directory for results", default="."
    )
    args = parser.parse_args()

    mol = bam.molecule(args.mol)

    residue_graph = mol.get_residue_graph(True)
    atom_graph = mol.get_atom_graph(True)

    edges = residue_graph.find_rotatable_edges(
        residue_graph.central_node, min_descendants=10
    )
    edges = residue_graph.sample_edges(edges, n=5)

    REPEATS = args.repeats

    init_clashes = str(mol.count_clashes())
    dir = Path(args.outdir)
    if not dir.exists():
        dir.mkdir()

    outfile = dir / f"{mol.id}-results.tsv"

    if (args.reset and outfile.exists()) or not outfile.exists():
        with open(outfile, "w") as f:
            f.write(
                "\t".join(
                    [
                        "env",
                        "detail",
                        "time",
                        "init_clashes",
                        "final_clashes",
                    ]
                )
                + "\n"
            )

    def record(env, outs, time, detail):
        data = [
            env.__name__,
            detail,
            str(time),
            init_clashes,
            None,
        ]
        with open(outfile, "a") as f:
            for clash in (str(i.count_clashes()) for i in outs):
                data[-1] = clash
                f.write("\t".join(data) + "\n")

    for env in [
        bam.optimizers.DistanceRotatron,
        bam.optimizers.OverlapRotatron,
        bam.optimizers.ForceFieldRotatron,
    ]:

        print(env.__name__)
        atom_start = time()
        atom_env = env(
            atom_graph,
            edges,
        )
        envs = [atom_env] * REPEATS
        atom_outs = bam.optimizers.parallel_optimize(mol, envs, unify_final=False)
        _time = time() - atom_start
        record(env, atom_outs, _time, "atom")

        final = bam.Molecule.empty()
        for i in atom_outs:
            final.merge(i)

        final.to_pdb(outfile.stem + ("." + env.__name__ + ".atom.pdb"))

        residue_start = time()
        residue_env = env(residue_graph, edges)
        envs = [residue_env] * REPEATS
        residue_outs = bam.optimizers.parallel_optimize(mol, envs, unify_final=False)
        _time = time() - residue_start
        record(env, residue_outs, _time, "residue")

        final = bam.Molecule.empty()
        for i in residue_outs:
            final.merge(i)

        final.to_pdb(outfile.stem + ("." + env.__name__ + ".residue.pdb"))
        print("done")
