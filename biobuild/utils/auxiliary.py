"""
Utility and auxiliary functions
"""

import math
import os
import re
import string

import pickle
import numpy as np

# =================================================================
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger

    HAS_RDKIT = True
except ImportError:
    AllChem = None
    Chem = None
    RDLogger = None
    HAS_RDKIT = False


try:
    from openbabel import pybel

    HAS_PYBEL = True
except ImportError:
    pybel = None
    HAS_PYBEL = False


try:
    import openmm.app

    HAS_OPENMM = True
except ImportError:
    openmm = None
    HAS_OPENMM = False


# =================================================================
def load_pickle(filename):
    """
    Load an object from a pickle file
    """
    with open(filename, "rb") as f:
        return pickle.load(f)


def save_pickle(obj, filename):
    """
    Save an object to a pickle file
    """
    with open(filename, "wb") as f:
        pickle.dump(obj, f)


def filename_to_id(filename):
    """
    Extract the id from a filename

    Parameters
    ----------
    filename : str
        The filename

    Returns
    -------
    str
        The id
    """
    base, suffix = os.path.splitext(os.path.basename(filename))
    return base


def change_suffix(filename, suffix):
    """
    Change the suffix of a filename

    Parameters
    ----------
    filename : str
        The filename
    suffix : str
        The new suffix

    Returns
    -------
    str
        The filename with the new suffix
    """
    base, _ = os.path.splitext(os.path.basename(filename))
    if suffix[0] != ".":
        suffix = "." + suffix
    return base + suffix


def chain_id_maker(cdx: int):
    """
    Make a string chain id from a counting integer

    Parameters
    ----------
    cdx : int
        The counting integer

    Returns
    -------
    str
        The chain id
    """
    if cdx < 26:
        return chr(cdx + 65)
    else:
        # recursive call
        return chain_id_maker(cdx // 26 - 1) + chain_id_maker(cdx % 26)


def make_formula(atoms):
    """
    Make a chemical formula from a list of atoms

    Parameters
    ----------
    atoms : list of Atoms
        The atoms. Each atom must have an element attribute (str).

    Returns
    -------
    str
        The chemical formula
    """
    formula = {}
    for atom in atoms:
        e = atom.element.title()
        if e in formula:
            formula[e] += 1
        else:
            formula[e] = 1
    return "".join([f"{atom}{formula[atom]}" for atom in sorted(formula.keys())])


def remove_nonprintable(text):
    """
    Remove non-printable characters from a string

    Parameters
    ----------
    text : str
        The string

    Returns
    -------
    str
        The string without non-printable characters
    """
    return "".join(filter(lambda x: x in string.printable, text))


class DummyStructure:
    """
    A dummy pdb structure

    Used for the surface residue inference

    Parameters
    ----------
    residues : list of Bio.PDB.Residue.Residue
        The residues to include in the structure
    """

    def __init__(self, residues) -> None:
        self.atoms = []
        self.atoms.extend(a for r in residues for a in r.get_atoms())
        self.residues = residues
        self.chains = []
        self.chains.extend(
            r.get_parent() for r in residues if r.get_parent() not in self.chains
        )
        self.models = []
        self.models.extend(
            c.get_parent() for c in self.chains if c.get_parent() not in self.models
        )
        self.level = "S"

    def get_atoms(self):
        return iter(self.atoms)

    def get_residues(self):
        return iter(self.residues)

    def get_chains(self):
        return iter(self.chains)

    def get_models(self):
        return iter(self.models)


class DummyBar:
    """
    A dummy progress bar
    """

    def __init__(self, *args, **kwargs) -> None:
        pass

    def __call__(self, *args, **kwargs) -> None:
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        pass


# # =================================================================
# # OLD STUFF
# # =================================================================


# def readLinesFromFile(fileName):
#     """Reads all lines in a file
#     Parameters:
#         fileName: path to file
#     Returns:
#         lines: list with all the lines in a file
#     """
#     with open(fileName, "r") as f:
#         lines = f.readlines()
#     return lines


# def topological_sort(unsorted_graph):
#     """Topological sorting of a graph
#     Parameters:
#         unsorted_graph: dictionary representation of a graph
#     Returns:
#         sorted_graph: list of nodes and corresponding edges
#     """
#     sorted_graph = []
#     # sort graph
#     while unsorted_graph:
#         acyclic = False

#         for node, edges in unsorted_graph.copy().items():
#             for edge in edges:
#                 if edge in unsorted_graph:
#                     break
#             else:
#                 acyclic = True
#                 del unsorted_graph[node]
#                 sorted_graph.append((node, edges))

#         if not acyclic:
#             print(
#                 "WARNING! Cyclique dependency occurred in ICs. Impossible to build residue"
#             )
#             print(unsorted_graph)
#             print(sorted_graph)
#             return ""
#     return sorted_graph[::-1]


# def rotation_matrix2(angle, direction, point=None):
#     """Return matrix to rotate about axis defined by point and direction."""
#     sina = math.sin(angle)
#     cosa = math.cos(angle)
#     direction = np.array(direction[:3])
#     direction = direction / math.sqrt(np.dot(direction, direction))
#     # rotation matrix around unit vector
#     R = np.diag([cosa, cosa, cosa])
#     R += np.outer(direction, direction) * (1.0 - cosa)
#     direction *= sina
#     R += np.array(
#         [
#             [0.0, -direction[2], direction[1]],
#             [direction[2], 0.0, -direction[0]],
#             [-direction[1], direction[0], 0.0],
#         ]
#     )
#     M = np.identity(4)
#     M[:3, :3] = R
#     if point is not None:
#         # rotation not around origin
#         point = np.array(point[:3], dtype=np.float64, copy=False)
#         M[:3, 3] = point - np.dot(R, point)
#     return M


# def alphanum_sort(l):
#     """Alphanumerical sort of a list from
#     https://arcpy.wordpress.com/2012/05/11/sorting-alphanumeric-strings-in-python/
#     Parameter:
#         l: list
#     Returns:
#         alphanumerically sorted list
#     """
#     convert = lambda text: int(text) if text.isdigit() else text
#     alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
#     return sorted(l, key=alphanum_key)


# aaa2a = {
#     "CYS": "C",
#     "ASP": "D",
#     "SER": "S",
#     "GLN": "Q",
#     "LYS": "K",
#     "ILE": "I",
#     "PRO": "P",
#     "THR": "T",
#     "PHE": "F",
#     "ASN": "N",
#     "GLY": "G",
#     "HIS": "H",
#     "LEU": "L",
#     "ARG": "R",
#     "TRP": "W",
#     "ALA": "A",
#     "VAL": "V",
#     "GLU": "E",
#     "TYR": "Y",
#     "MET": "M",
# }

# biobuild_PATH = os.path.dirname(os.path.realpath(__file__))
# # SELF_BIN = os.path.dirname(os.path.realpath(sys.argv[0]))
# # sys.path.insert(0, SELF_BIN + '/support')
