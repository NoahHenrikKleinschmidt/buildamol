"""
This script can be used to generate a PDBECompounds database pickle object from a PDBE CIF file.
"""

import argparse
import biobuild as bb


def setup():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=str, help="The input PDBE CIF file")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="The output pickle file. By default, the same as the input filename, but with '.comp.pkl'. Set to 'None' to disable saving (in case you only want to set the default database)",
        default=None,
    )
    parser.add_argument(
        "-s",
        "--setdefault",
        action="store_true",
        help="Set the database as the default database for biobuild (overwriting the current one in the process)",
        default=False,
    )
    return parser.parse_args()


def main(args):
    db = bb.resources.PDBECompounds.from_file(args.input)
    if args.output is None:
        args.output = args.input + ".comp.pkl"
    if args.output != "None":
        db.save(args.output)
    if args.setdefault:
        bb.set_default_compounds(db, overwrite=True)


if __name__ == "__main__":
    main(setup())
