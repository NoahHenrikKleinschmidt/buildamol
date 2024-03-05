"""
Utility and auxiliary functions
"""

import os

# import re
import string

import pickle
import numpy as np

# =================================================================
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    from rdkit.Chem import Draw

    HAS_RDKIT = True
except ImportError:
    AllChem = None
    Chem = None
    Draw = None
    RDLogger = None
    HAS_RDKIT = False


try:
    from openbabel import pybel

    HAS_PYBEL = True
except ImportError:
    pybel = None
    HAS_PYBEL = False


try:
    import openmm.app as openmm

    HAS_OPENMM = True
except ImportError:
    openmm = None
    HAS_OPENMM = False


try:
    from numba import njit

    HAS_NUMBA = True
    USE_NUMBA = False
    USE_ALL_NUMBA = False
except ImportError:
    njit = None
    HAS_NUMBA = False
    USE_NUMBA = False
    USE_ALL_NUMBA = False


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


def get_args(func, namespace):
    """
    Filter a dictionary based on the argument-namespace of a function

    Parameters
    ----------
    func : function
        The function whose argument-namespace will be used for filtering
    namespace : dict
        The dictionary to be filtered

    Returns
    -------
    dict
        The filtered dictionary
    """
    arg_names = func.__code__.co_varnames[: func.__code__.co_argcount]
    kwargs = {k: v for k, v in namespace.items() if k in arg_names}
    return kwargs


def element_range(symbol: str, n: int, start: int = 1):
    """
    Generate a range of systemic atom ids based on element symbol and number of atoms

    Parameters
    ----------
    symbol : str
        The element symbol
    n : int
        The number of elements
    start : int, optional
        The start number, by default 1

    Returns
    -------
    list of str
        The element names

    Examples
    --------
    >>> element_range("C", 3)
    ['C1', 'C2', 'C3']
    """
    return [f"{symbol}{i}" for i in range(start, start + n)]


def coord_array(*objs) -> np.ndarray:
    """
    Creates a numpy array of coordinates from objects with a get_coord() method or a coord attribute
    """
    coords = []
    for obj in objs:
        if hasattr(obj, "get_coord"):
            coords.append(obj.get_coord())
        elif hasattr(obj, "coord"):
            coords.append(obj.coord)
        else:
            raise ValueError(f"Object {obj} has no get_coord() or coord attribute")
    return np.array(coords)


def use_numba():
    """
    Use Numba if available to speed up some functions
    """
    global USE_NUMBA
    if HAS_NUMBA:
        USE_NUMBA = True


def use_all_numba():
    """
    Use Numba if available to speed up all functions
    """
    global USE_ALL_NUMBA
    if HAS_NUMBA:
        USE_ALL_NUMBA = True


def dont_use_numba():
    """
    Don't use Numba
    """
    global USE_NUMBA
    USE_NUMBA = False


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
