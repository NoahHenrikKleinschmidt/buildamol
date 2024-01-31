import buildamol as bam
from time import time
import numpy as np

import Bio.PDB as PDB

import argparse
from alive_progress import alive_bar


def setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "file",
        type=str,
        help="The file to optimize. This can be any Biobuild readable file.",
    )
    parser.add_argument("-n", type=int, default=10)
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output file prefix. Used for both the stats and the visualizations.",
    )
    parser.add_argument(
        "--atomgraph", action="store_true", help="Use an AtomGraph for optimization."
    )
    parser.add_argument(
        "--no-genetic",
        nargs="+",
        help="Do not use the genetic algorithm for optimization (because it takes forever on larger molecules). Specify which environments to skip using their names DistanceRotatron, OverlapRotatron, ForceFieldRotatron. All three or only some can be specified.",
    )
    parser.add_argument(
        "--no-scipy",
        nargs="+",
        help="Do not use the scipy algorithm for optimization (because it takes forever on larger molecules). Specify which environments to skip using their names DistanceRotatron, OverlapRotatron, ForceFieldRotatron. All three or only some can be specified.",
    )
    parser.add_argument(
        "--residue-connections",
        action="store_true",
        help="Use residue connections instead of finding rotatable edges.",
    )
    parser.add_argument(
        "--min-descendants",
        type=int,
        default=5,
        help="Minimum number of descendants for a rotatable edge to be considered.",
    )
    parser.add_argument(
        "--max-descendants",
        type=int,
        default=500,
        help="Maximum number of descendants for a rotatable edge to be considered.",
    )
    return parser.parse_args()


algorithms = [
    bam.optimizers.genetic_optimize,
    bam.optimizers.swarm_optimize,
    bam.optimizers.anneal_optimize,
    bam.optimizers.scipy_optimize,
]

environments = [
    bam.optimizers.DistanceRotatron,
    # bam.optimizers.OverlapRotatron,
    # bam.optimizers.ForceFieldRotatron,
]


def write_pdb(mol, filename):
    mol._model.id = mdx
    mdx += 1
    io.set_structure(mol.to_pdb_structure())
    io.save(filename)
    bam.utils.pdb.write_connect_lines(mol, filename)


if __name__ == "__main__":
    args = setup()

    # class Args:
    #     file = "/Users/noahhk/GIT/biobuild/biobuild/optimizers/_testing/files/GL2.json"
    #     n = 10
    #     o = None
    #     atomgraph = False

    # args = Args()
    args.output = args.output if args.output else ".".join(args.file.split(".")[:-1])

    mol = bam.molecule(args.file)
    mol.autolabel()
    mol._model.id = 1
    biostruct = mol.to_biopython()
    mdx = 2

    if args.atomgraph:
        graph = mol.get_atom_graph()
    else:
        graph = mol.get_residue_graph()
        graph.make_detailed(n_samples=0.7, include_far_away=True)
    if args.residue_connections:
        edges = mol.get_residue_connections()
        edges = graph.directed_edges(None, edges)
    else:
        edges = graph.find_rotatable_edges(
            graph.central_node,
            min_descendants=args.min_descendants,
            max_descendants=args.max_descendants,
        )

    statsfile = open(args.output + ".stats.tsv", "w")
    statsline = "{file}\t{env}\t{algo}\t{atomgraph}\t{time}\t{clashes}\t{eval}"
    statsfile.write(
        statsline.format(
            file="file",
            env="env",
            atomgraph="atomgraph",
            algo="algo",
            time="time",
            clashes="clashes",
            eval="eval",
        )
        + "\n"
    )

    def concatenation_function(self, x):
        smallest = np.sort(x)[: self.n_smallest]
        penalty = np.sum(x < 1.5 * self.clash_distance)
        e = np.mean(x) ** self.unfold + np.mean(smallest) ** self.pushback
        e /= (1 + penalty) ** 2
        return e

    # return 0.5 * np.mean(x) + self.pushback * smallest

    with alive_bar(args.n * len(environments) * len(algorithms)) as bar:
        for i, env in enumerate(environments):
            E = env(
                graph,
                edges,
                # n_smallest=10,
                pushback=3,
                unfold=2,
                # radius=12,
                concatenation_function=concatenation_function,
            )

            line = statsline.format(
                file=args.file,
                env=E.__class__.__name__,
                atomgraph=args.atomgraph,
                algo="None",
                time="0",
                clashes=mol.count_clashes(),
                eval=str(E._best_eval).replace("[", "").replace("]", ""),
            )
            statsfile.write(line + "\n")

            for j, algo in enumerate(algorithms):
                if (
                    algo.__name__ == "genetic_optimize"
                    and args.no_genetic is not None
                    and E.__class__.__name__ in args.no_genetic
                ) or (
                    algo.__name__ == "scipy_optimize"
                    and args.no_scipy is not None
                    and E.__class__.__name__ in args.no_scipy
                ):
                    for k in range(args.n):
                        line = statsline.format(
                            file=args.file,
                            env=E.__class__.__name__,
                            atomgraph=args.atomgraph,
                            algo=algo.__name__,
                            time="NaN",
                            clashes="NaN",
                            eval="NaN",
                        )
                        statsfile.write(line + "\n")
                        bar()
                    continue
                for k in range(args.n):
                    t1 = time()
                    sol, _eval = algo(
                        E,
                        # c1=0.2,
                        # mutants=0.5,
                        # parents=0.25,
                        # children=0.1,
                        # newcomers=0.15,
                        # population_size=50,
                        # max_generations=20,
                        # stop_if_done=False,
                        # threshold=1e-6,
                    )
                    delta_t = time() - t1
                    _eval = str(_eval).replace("[", "").replace("]", "")

                    if sol.shape[0] != len(E.rotatable_edges):
                        sol = sol[0]
                    out = bam.optimizers.apply_solution(sol, E, mol.copy())
                    clashes = out.count_clashes()
                    print(
                        f"{E.__class__.__name__} - {algo.__name__} - {delta_t:.2f}s, {clashes}"
                    )
                    # if clashes >= 6:
                    #     out.show()

                    line = statsline.format(
                        file=args.file,
                        env=E.__class__.__name__,
                        atomgraph=args.atomgraph,
                        algo=algo.__name__,
                        time=delta_t,
                        clashes=clashes,
                        eval=_eval,
                    )
                    statsfile.write(line + "\n")
                    out._model.id = mdx
                    mdx += 1
                    biostruct.add(out.to_biopython().child_list[0])
                    bar()
                    E.reset()

        io = PDB.PDBIO()
        io.set_structure(biostruct)
        io.save(f"{args.output}.pdb")
        bam.utils.pdb.write_connect_lines(mol, f"{args.output}.pdb")
        with open(f"{args.output}.pdb", "r") as f:
            c = f.read()
        with open(f"{args.output}.pdb", "w") as f:
            c = c.replace("MDL\n", "ENDMDL\n")
            f.write(c)

    statsfile.close()
    print("Done!")
