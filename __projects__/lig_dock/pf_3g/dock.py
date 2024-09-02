import dockstring as ds
import argparse as ap

parser = ap.ArgumentParser(description="Docking simulation")
parser.add_argument("ligand", type=str, help="Ligand file")
parser.add_argument("target", type=str, help="Target file")
parser.add_argument("output", type=str, help="Output file")
args = parser.parse_args()
