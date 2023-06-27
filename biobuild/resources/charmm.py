"""
This module defines File parsers and data classes to work with the CHARMM force field

Reading CHARMM files
====================

CHARMM files can be read using the `CHARMMTopology` class.

.. code-block:: python

    from biobuild.resources import CHARMMTopology

    charmm_topology_file = "top_all36_prot.rtf"

    # load a CHARMM topology file
    top = CHARMMTopology.from_file(charmm_topology_file)


Because parsing can be a time intensive task, it is recommended to `save` the parsed objects to a pickle file for later use.
In this case a pre-parsed topology object can be loaded using the `load` method.

.. code-block:: python

    # save the parsed objects to a pickle file
    top.save("top_all36_prot.pkl")

    # load a CHARMM topology file
    top = CHARMMTopology.load("top_all36_prot.pkl")

    
Working with CHARMM objects
===========================

The `CHARMMTopology` class includes methods to work with the parsed data, 
such as `get_residue` or `get_mass`. They also support adding new data via the corresponding `add_{...}` methods.


Setting default CHARMM objects
==============================

The `biobuild.utils.defaults` module pre-loads default CHARMM topology and parameters objects for convenience. Many methods that make use of 
these objects such as the `attach` method of `Molecule` instances also accept arguments for custom topology objects.
For convenience, a custom object can be set as the default, however, using the `set_default_topology` functions.

.. code-block:: python

    from biobuild.resources import set_default_topology
    
    # set a custom topology object as the default
    set_default_topology(top)


.. warning::

    The `set_default_{...}` functions include an argument `overwrite` that is set to `False` by default. If set to `True` the default object is permanently overwritten
    and will be used for all future sessions. If set to `False` (default), the default object is only set for the current session and any future sessions will use the
    original defaults. In case the defaults were overwritten, they can be reset to the original using the `restore_default_{...}` functions.

    .. code-block:: python

        # set a custom topology object as the default
        # for the current session only
        set_default_topology(top)

        # set a custom topology object as the default
        # for all future sessions
        set_default_topology(top, overwrite=True)

        # restore the original default topology object
        # for all future sessions
        restore_default_topology()
"""

import os
import pickle

import re
import warnings

import biobuild.utils as utils
import biobuild.core.Linkage as Linkage
import biobuild.utils.defaults as _defaults

# ===================================================================
# Base Parser
# ===================================================================


def load_topology(filename: str, set_default: bool = True) -> "CHARMMTopology":
    """
    Load a CHARMM topology from a pickle file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    set_default: bool
        If `True`, the topology is set as the default topology

    Returns
    -------
    CHARMMTopology
        The topology object
    """
    top = CHARMMTopology.load(filename)
    if set_default:
        set_default_topology(top)
    return top


def read_topology(filename: str, set_default: bool = True) -> "CHARMMTopology":
    """
    Read a CHARMM topology file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    set_default: bool
        If `True`, the topology is set as the default topology

    Returns
    -------
    CHARMMTopology
        The parsed topology object
    """
    top = CHARMMTopology.from_file(filename)
    if set_default:
        set_default_topology(top)
    return top


def save_topology(filename: str, topology: "CHARMMTopology" = None):
    """
    Save a CHARMM topology to a pickle file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    topology: CHARMMTopology
        The topology object. If `None`, the default topology is used.
    """
    if topology is None:
        topology = get_default_topology()
    topology.save(filename)


def set_default_topology(obj, overwrite: bool = False):
    """
    Set a CHARMMTopology object as the new default.

    Parameters
    ----------
    obj : CHARMMTopology
        The CHARMMTopology object to set as the new default
    overwrite : bool
        If set to `True`, the new object will be permanently saved as
        the default. Otherwise, the new object will only be used for
        the current session.
    """
    if obj.__class__.__name__ != "CHARMMTopology":
        raise TypeError("The object must be a CHARMMTopology instance.")
    if overwrite:
        current = _defaults.__default_instances__.get("Topology", None)
        if current:
            current.save(_defaults.DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak")
    _defaults.__default_instances__["Topology"] = obj
    if overwrite:
        obj.save(_defaults.DEFAULT_CHARMM_TOPOLOGY_FILE)


def get_default_topology() -> "CHARMMTopology":
    """
    Get the default CHARMMTopology object

    Returns
    -------
    CHARMMTopology
        The default CHARMMTopology object
    """
    return _defaults.__default_instances__.get("Topology", None)


def has_patch(name: str) -> bool:
    """
    Check if a patch is defined in the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the patch

    Returns
    -------
    bool
        `True` if the patch is defined, `False` otherwise
    """
    return get_default_topology().has_patch(name)


has_linkage = has_patch


def available_patches():
    """
    Get a list of available patches.

    Returns
    -------
    list
        A list of available patches
    """
    return get_default_topology().patches


available_linkages = available_patches


def add_patch(patch, overwrite: bool = False):
    """
    Add a patch to the CHARMM topology.

    Parameters
    ----------
    patch: Patch
        The patch object
    overwrite: bool
        If `True`, the topology with the added patch is saved as the new default topology
    """
    get_default_topology().add_patch(patch)
    if overwrite:
        set_default_topology(get_default_topology(), overwrite=True)


add_linkage = add_patch


def has_residue(name: str) -> bool:
    """
    Check if a residue is defined in the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the residue

    Returns
    -------
    bool
        `True` if the residue is defined, `False` otherwise
    """
    return get_default_topology().has_residue(name)


def in_charmm_topology(name: str) -> bool:
    """
    Check if a residue or patch is defined in the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the residue or patch

    Returns
    -------
    bool
        `True` if the residue or patch is defined, `False` otherwise
    """
    has_patch = get_default_topology().has_patch(name)
    has_residue = get_default_topology().has_residue(name)
    return has_patch or has_residue


def available_residues():
    """
    Get a list of available residues.

    Returns
    -------
    list
        A list of available residues
    """
    return get_default_topology().residues


def restore_default_topology():
    """
    Restore the default CHARMMTopology object from the backup file
    """
    _defaults.__default_instances__["Topology"] = CHARMMTopology.load(
        DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak"
    )
    os.remove(DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak")


def get_patch(name: str):
    """
    Get a patch from the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name of the patch

    Returns
    -------
    Patch
        The patch object
    """
    return get_default_topology().get_patch(name)


get_linkage = get_patch


def add_patch(patch, overwrite: bool = False):
    """
    Add a patch to the CHARMM topology.

    Parameters
    ----------
    patch: Patch
        The patch object
    overwrite: bool
        If `True`, the topology with the added patch is saved to a pickle file and will be used as the default topology for all future sessions.
    """
    top = get_default_topology()
    top.add_patch(patch)
    set_default_topology(top, overwrite)


add_linkage = add_patch


class CHARMMParser:
    """
    The base class to parse CHARMM files
    """

    def __init__(self) -> None:
        self.id = None
        self._dict = {
            "masses": {},
        }
        self._file = None
        self._was_loaded_from_pickle = False

    @classmethod
    def load(cls, filename: str):
        """
        Load a pre-parsed object from a pickle file.

        Parameters
        ----------
        filename: str
            The path to the pickle file
        """
        new = cls()
        _data = pickle.load(open(filename, "rb"))
        new._vet_load(_data)
        new._dict = _data
        new._file = filename
        new.id = os.path.basename(filename)
        new._was_loaded_from_pickle = True
        return new

    @classmethod
    def from_file(cls, filename: str):
        """
        Parse a CHARMM file and return a new object.

        Parameters
        ----------
        filename: str
            The path to the CHARMM file
        """
        new = cls()
        new._parse(filename)
        new._file = filename
        new.id = os.path.basename(filename)
        return new

    @property
    def masses(self):
        return self._dict["masses"]

    def save(self, filename: str = None):
        """
        Save the data dictionary as binary file that can be directly loaded again.

        Parameters
        ----------
        filename: str
            The path to the pickle file to save in. By default, this will be the same filename as the one from which the data was loaded or parsed (adding the file-suffix `.pkl`)
        """
        if not filename:
            if not self._file:
                raise ValueError(
                    "No filename was given and no filename from a source file is available!"
                )
            filename = self._file
            if not self._was_loaded_from_pickle and not filename.endswith(".pkl"):
                filename += ".pkl"
        pickle.dump(self._dict, open(filename, "wb"))

    def get_mass(self, id):
        """
        Get a mass by its atom type ID

        Parameters
        ----------
        id : str
            The ID of the atom type

        Returns
        -------
        float
            The mass
        """
        if isinstance(id, (list, tuple)):
            return [self.get_mass(i) for i in id]
        return self._dict["masses"].get(id)

    def add_mass(self, key, mass):
        """
        Add a mass to the topology

        Parameters
        ----------
        key : str
            The ID of the atom type
        mass : float
            The mass of the atom type
        """
        if isinstance(key, (list, tuple)):
            for k, m in zip(key, mass):
                self.add_mass(k, m)
            return
        self._dict["masses"][key] = float(mass)

    def _parse(self, filename: str):
        """
        Parse a CHARMM file and store the data in a dictionary

        Parameters
        ----------
        filename: str
            The path to the CHARMM file
        """
        raise NotImplementedError(
            "This method needs to be implemented by the subclass!"
        )

    def _vet_load(self, _data):
        """
        An optional method to vet the loaded data for specific properties
        """
        pass

    @staticmethod
    def _read_line(line: str):
        """
        Reads a line of the CHARMM file

        Parameters
        ----------
        line : str
            The line of the file

        Returns
        -------
        list
            The line split into its components
        """
        return line.strip().split("!")[0].split()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.id})"


# ===================================================================
# Topology
# ===================================================================


class CHARMMTopology(CHARMMParser):
    """
    A data class and file parser for CHARMM topology files.

    Attributes
    ----------
    id : str
        The ID of the topology file
    residues : list of AbstractResidue
        The residues in the topology file
    patches : list of Linkage
        The Linkages (patches, in the CHARMM nomenclature, prefixed with 'PRES') in the topology file
    masses : list of floats
        The masses in the topology file
    """

    __tags__ = (
        "MASS",
        "RESI",
        "PRES",
        "GROUP",
        "ATOM",
        "BOND",
        "DOUBLE",
        "IMPROPER",
        "IMPHI",
        "DIHE",
        "PATCH",
        "END",
    )

    def __init__(self, id=None):
        super().__init__()
        self._dict["residues"] = {}
        self._dict["patches"] = {}
        self.id = id

    @property
    def residues(self):
        return list(self._dict["residues"].values())

    @property
    def patches(self):
        return list(self._dict["patches"].values())

    def has_residue(self, id):
        """
        Check if a residue is in the topology

        Parameters
        ----------
        id : str
            The ID of the residue

        Returns
        -------
        bool
            True if the residue is in the topology
        """
        if isinstance(id, (list, tuple)):
            return all(self.has_residue(i) for i in id)
        return id in self._dict["residues"]

    def get_residue(self, id):
        """
        Get a residue by its ID

        Parameters
        ----------
        id : str
            The ID of the residue

        Returns
        -------
        AbstractResidue
            The residue
        """
        if isinstance(id, (list, tuple)):
            return [self.get_residue(i) for i in id]
        return self._dict["residues"].get(id)

    def find_residue(self, *atoms):
        """
        Find a residue by its atoms. This will require the
        residue to have all atoms in the list.

        Parameters
        ----------
        atoms : tuple of AbstractAtom
            The atoms to search for

        Returns
        -------
        AbstractResidue
            The residue
        """
        for residue in self.residues:
            if all(residue.has_atom(atom) for atom in atoms):
                return residue
        return None

    def add_residue(self, residue):
        """
        Add a residue to the topology

        Parameters
        ----------
        residue : AbstractResidue
            The residue to add
        """
        if isinstance(residue, (list, tuple)):
            for r in residue:
                self.add_residue(r)
            return
        elif not isinstance(residue, utils.abstract.AbstractResidue):
            raise TypeError("The residue must be an instance of AbstractResidue")
        self._dict["residues"][residue.id] = residue

    def get_patch(self, id):
        """
        Get a patch by its ID

        Parameters
        ----------
        id : str
            The ID of the patch

        Returns
        -------
        Linkage
            The patch
        """
        if isinstance(id, (list, tuple)):
            return [self.get_patch(i) for i in id]
        return self._dict["patches"].get(id)

    def has_patch(self, id):
        """
        Check if a patch is in the topology

        Parameters
        ----------
        id : str
            The ID of the patch

        Returns
        -------
        bool
            True if the patch is in the topology
        """
        if isinstance(id, (list, tuple)):
            return all(self.has_patch(i) for i in id)
        return id in self._dict["patches"]

    def add_patch(self, patch):
        """
        Add a patch to the topology

        Parameters
        ----------
        patch : Linkage
            The patch to add
        """
        if isinstance(patch, (list, tuple)):
            for p in patch:
                self.add_patch(p)
            return
        elif not isinstance(patch, Linkage):
            raise TypeError("The patch must be an instance of Linkage")
        self._dict["patches"][patch.id] = patch

    def _parse(self, filename: str):
        """
        Reads and parses the data from a CHARMM Topology file

        Parameters
        ----------
        filename : str
            The path to the CHARMM Topology file
        """
        with open(filename, "r") as file:
            lines = file.read().split("\n")  # readlines but remove the endlines
            lines = [
                line.strip().split("!") for line in lines
            ]  # get rid of all comments
            lines = [
                line[0] for line in lines if line[0] != ""
            ]  # get rid of all empty lines

        idx = 0
        while idx < len(lines):
            line = lines[idx]
            if line.startswith("MASS"):
                self._parse_mass(line)
            elif line.startswith("RESI"):
                idx = self._parse_residue(lines, idx)
            elif line.startswith("PRES"):
                idx = self._parse_patch(lines, idx)
            idx += 1

        self._adopt_atom_masses()
        self._make_ICs_improper()
        self._file = filename

    def _parse_mass(self, line: str):
        """
        Parses the mass information from a line of the topology file

        Parameters
        ----------
        line : str
            The line of the topology file
        """
        line = self._read_line(line)
        id = line[2]
        mass = float(line[3])
        self._dict["masses"][id] = mass

    def _parse_residue(self, lines: list, idx: int):
        """
        Parses the residue information from a line of the topology file

        Parameters
        ----------
        lines: list
            The list of lines of the topology file
        idx: int
            The index of the current line
        """
        line = self._read_line(lines[idx])

        residue = utils.abstract.AbstractResidue(id=line[1])

        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            start = line[0]

            if start == "" or re.match("GROU(P| )", start):
                idx += 1
                continue

            if start == "ATOM":
                self._parse_atom(line, residue)

            elif start == "BOND":
                self._parse_bond(line, residue)

            elif start == "IC":
                self._parse_ic(line, residue)

            elif start in self.__tags__:
                idx -= 1
                break
            idx += 1

        self._dict["residues"][residue.id] = residue
        return idx

    def _parse_patch(self, lines: list, idx: int):
        """
        Parses the patch information from a line of the topology file

        Parameters
        ----------
        lines: list
            The list of lines of the topology file
        idx: int
            The index of the current line
        """
        line = self._read_line(lines[idx])

        patch = Linkage(id=line[1])

        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            start = line[0]

            if start == "" or re.match("GROU(P| )", start):
                idx += 1
                continue

            if start == "ATOM":
                self._parse_atom(line, patch)

            elif re.match("dele(te)?", start.lower()):
                self._parse_delete(line, patch)

            elif start == "BOND":
                self._parse_bond(line, patch)

            elif start == "IC":
                self._parse_ic(line, patch)

            elif start in self.__tags__:
                idx -= 1
                break
            idx += 1

        self._dict["patches"][patch.id] = patch
        return idx

    def _adopt_atom_masses(self):
        """
        Adopt the atom masses from the topology file
        """
        for residue in self.residues:
            for atom in residue.atoms:
                if atom.mass is None:
                    atom.mass = self._dict["masses"][atom.type]

    def _make_ICs_improper(self):
        """
        Ensure that improper internal coordinates are also labelled as such
        based on whether their third atom is connected to the first atom
        """
        for residue in self.residues:
            for ic in residue.internal_coordinates:
                if not ic.improper and residue.get_bond(ic.atom1, ic.atom3):
                    ic.improper = True
                    # Actually a bad idea to do this,
                    # because the file still stores 1-2 lengths instead of 1-3 lengths!
                    # if not ic.bond_length_13:
                    #     ic.bond_length_13 = ic.bond_length_12
                    #     ic.bond_length_12 = None

    def _vet_load(self, _data):
        """
        Checks that the data loaded from a file is valid
        """
        if not isinstance(_data, dict):
            raise TypeError("The file must contain a dictionary object")
        if (
            "residues" not in _data.keys()
            or "patches" not in _data.keys()
            or "masses" not in _data.keys()
        ):
            raise KeyError(
                "The dictionary must contain 'residues', 'masses', and 'patches' keys"
            )

    @staticmethod
    def _parse_ic(line: list, obj):
        """
        Parses the internal coordinate information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : Linkage
            The object to which the internal coordinate should be added
        """
        is_improper = "*" in line[3]
        if is_improper:
            line[3] = line[3][1:]

        atom1 = line[1]
        atom2 = line[2]
        atom3 = line[3]
        atom4 = line[4]

        if isinstance(obj, utils.abstract.AbstractResidue):
            atom1 = obj.get_atom(atom1)
            atom2 = obj.get_atom(atom2)
            atom3 = obj.get_atom(atom3)
            atom4 = obj.get_atom(atom4)

            if atom1 is None or atom2 is None or atom3 is None or atom4 is None:
                warnings.warn(
                    f"[ignoring line] Found an invalid internal coordinate in {line}"
                )
                return

        if is_improper:
            _bond_lengths = {
                "bond_length_12": None,
                "bond_length_13": float(line[5]),
                "bond_length_34": float(line[9]),
            }
        else:
            _bond_lengths = {
                "bond_length_12": float(line[5]),
                "bond_length_13": None,
                "bond_length_34": float(line[9]),
            }

        ic = utils.ic.InternalCoordinates(
            atom1=atom1,
            atom2=atom2,
            atom3=atom3,
            atom4=atom4,
            bond_angle_123=float(line[6]),
            dihedral=float(line[7]),
            bond_angle_234=float(line[8]),
            improper=is_improper,
            **_bond_lengths,
        )

        obj.add_internal_coordinates(ic)

    @staticmethod
    def _parse_atom(line: list, obj):
        """
        Parses the atom information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : Linkage
            The object to which the atom belongs
        """
        atom = utils.abstract.AbstractAtom(
            id=line[1], type=line[2], charge=float(line[3])
        )
        obj.add_atom(atom)

    @staticmethod
    def _parse_bond(line: list, obj):
        """
        Parses the bond information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        obj : Linkage
            The object to which the bond belongs
        """
        line = line[1:]

        # split line into tuple pairs
        line = [(line[i], line[i + 1]) for i in range(0, len(line), 2)]

        if isinstance(obj, utils.abstract.AbstractResidue):
            for a1, a2 in line:
                atom1 = obj.get_atom(a1)
                atom2 = obj.get_atom(a2)
                if atom1 is None or atom2 is None:
                    warnings.warn(f"[ignoring bond] Found an invalid bond in {line}")
                    return

                bond = utils.abstract.AbstractBond(atom1, atom2)
                obj.add_bond(bond)

        elif isinstance(obj, Linkage):
            for a1, a2 in line:
                bond = utils.abstract.AbstractBond(a1, a2)
                obj.add_bond(bond)

    @staticmethod
    def _parse_delete(line: list, patch):
        """
        Parses the delete information from a line of the topology file

        Parameters
        ----------
        line : list
            The line of the topology file, split into a list
        patch : Linkage
            The patch to which the delete belongs
        """
        id = line[2]
        patch.add_delete(id)


# Set the default topology to the CHARMM topology
set_default_topology(CHARMMTopology.load(_defaults.DEFAULT_CHARMM_TOPOLOGY_FILE))


__all__ = [
    "CHARMMTopology",
    "get_default_topology",
    "set_default_topology",
    "restore_default_topology",
    "has_patch",
    "has_residue",
    "available_patches",
    "available_linkages",
    "available_residues",
    "add_patch",
    "add_linkage",
    "has_linkage",
    "get_patch",
    "get_linkage",
    "read_topology",
    "load_topology",
    "save_topology",
]

if __name__ == "__main__":
    _carbs = "/Users/noahhk/GIT/biobuild/support/toppar_charmm/carbohydrates.rtf"
    _top = CHARMMTopology.from_file(_carbs)
    print(_top)

    from biobuild.utils.defaults import DEFAULT_CHARMM_TOPOLOGY_FILE

    _save_to = "/Users/noahhk/GIT/biobuild/biobuild/resources/"
    _top.save(_save_to + os.path.basename(DEFAULT_CHARMM_TOPOLOGY_FILE))
