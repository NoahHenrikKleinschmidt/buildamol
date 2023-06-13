"""
This script takes one pubchem JSON file and corresponding SDF structure file and converts them to a PDBE-CIF file that can be imported into biobuild. 
The created CIF file is not a full-format CIF file, but a CIF file that contains only the information needed for biobuild!
"""

import argparse
from biobuild.resources import pubchem


def setup():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("json", help="Input JSON file")
    parser.add_argument("sdf", help="Input SDF file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output CIF file. By default the output filename is the same as the JSON file with the extension changed to .cif",
        default=None,
    )
    parser.add_argument(
        "-id",
        "--compound-id",
        help="Set a custom compound id. By default an ID is inferred from the JSON file.",
        default=None,
    )
    args = parser.parse_args()
    return args


def main(args):
    pubchem.pubchem_to_cif(args.json, args.sdf, args.compound_id, args.output)


if __name__ == "__main__":
    main(setup())
