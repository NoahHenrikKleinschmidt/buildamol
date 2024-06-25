"""
This script takes one PDBECompounds pickle file and filters out all molecules that are larger than a certain size.
"""

import argparse
import buildamol as bam
import os


def setup():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=str, help="The input PDBECompounds pickle file")
    parser.add_argument("output", type=str, help="The output PDBECompounds pickle file")
    parser.add_argument(
        "--max_size",
        type=int,
        default=80,
        help="The maximum size of the molecule to keep. Default is 80.",
    )
    return parser.parse_args()


def main(args):
    if args.input == args.output:
        args.output += ".tmp"
    c = bam.read_compounds(args.input, False)
    ids = list(c.ids)
    for i in ids:
        pdb = c._pdb[i]
        if len(pdb["atoms"]["ids"]) > args.max_size:
            c.remove(i)
    c.save(args.output)
    if args.output.endswith(".tmp"):
        args.output = args.output[:-4]
        os.rename(args.output + ".tmp", args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main(setup())
    # class Args:
    #     input = (
    #         "/Users/noahhk/GIT/biobuild/support/pdbe_compounds/components.cif.comp.pkl"
    #     )
    #     output = (
    #         "/Users/noahhk/GIT/biobuild/support/pdbe_compounds/components.max40.pkl"
    #     )
    #     max_size = 40

    # main(Args())
