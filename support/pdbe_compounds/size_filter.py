"""
This script takes one PDBE CIF file and and removes any entries that describe molecules exceeding a certain size.
"""

import argparse
import re
import os


def setup():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=str, help="The input PDBE CIF file")
    parser.add_argument("output", type=str, help="The output PDBE CIF file")
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
    with open(args.input, "r") as fi:
        with open(args.output, "w") as fo:
            compound_fragment = []
            fragment_size = 0
            start_counting = False
            may_write = True
            for line in fi:
                if may_write:
                    compound_fragment.append(line)
                if line.startswith("_chem_comp_atom.pdbx_ordinal"):
                    start_counting = True
                    continue
                if start_counting:
                    fragment_size += 1
                if fragment_size > 0 and (
                    line.startswith("#") or line.startswith("loop_")
                ):
                    start_counting = False
                    fragment_size -= 1
                    if fragment_size <= args.max_size:
                        may_write = True
                    else:
                        may_write = False
                        compound_fragment.clear()
                        fragment_size = 0
                if line.startswith("data_") and fragment_size > 0:
                    for line in compound_fragment:
                        fo.write(line)
                    compound_fragment.clear()
                    fragment_size = 0
                    start_counting = False
                    may_write = True

    if args.output.endswith(".tmp"):
        args.output = args.output[:-4]
        os.rename(args.output + ".tmp", args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main(setup())
    # class Args:
    #     input = "/Users/noahhk/GIT/biobuild/support/pdbe_compounds/components.cif"
    #     output = (
    #         "/Users/noahhk/GIT/biobuild/support/pdbe_compounds/components.beta.cif"
    #     )

    # main(Args())
