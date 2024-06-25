"""
This module defines File parsers and data classes to work with the CHARMM force field

CHARMM files can be read using the `CHARMMParser` class. Which implements a child-class `CHARMMTopology` that is used to read and store CHARMM topology files.
Specifically, the `CHARMMTopology` stores the linkage (patches) data from the topology file.

Reading CHARMM files
====================

The quickest way to read a CHARMM topology file is to use the `read_topology` function.

.. code-block:: python

    import buildamol as bam 

    charmm_topology_file = "top_all36_prot.rtf"

    # read a CHARMM topology file
    top = bam.read_topology(charmm_topology_file)

There is an optional argument `set_default` that is set to `True` by default. If set to `True` the parsed object is set as the default object for the current session.

Alternatively, the `CHARMMTopology` class can be used to read a CHARMM topology file.

.. code-block:: python

    from buildamol.resources import CHARMMTopology
    
    # read a CHARMM topology file
    top = CHARMMTopology.from_file(charmm_topology_file)
    # works the same as
    # top = bam.read_topology(charmm_topology_file)

    
.. note::

    There is currently no implementation of a CHARMM parameter file parser, because BuildAMol does not use this data.
    However, the `CHARMMParser` can be used as a base class to implement a custom parser for CHARMM parameter files.


Saving and loading CHARMM objects
=================================
        
Because parsing can be a time intensive task, it is recommended to `save` the parsed objects to a pickle file for later use.
In this case a pre-parsed topology object can be loaded using the `load` method. Both methods accept a file path as argument,
and both also have a functional equivalent called `save_topology` and `read_topology` (yes, that's the same as above).

.. code-block:: python

    # save the parsed objects to a pickle file
    bam.save_topology(top, "top_all36_prot.pkl")
    # or 
    top.save("top_all36_prot.pkl")

    # load a CHARMM topology file
    top = bam.read_topology("top_all36_prot.pkl")
    # or
    top = CHARMMTopology.load("top_all36_prot.pkl")

    
Topologies also support an encoding agnosting export to JSON format using the `to_json` method or XML format using the `to_xml` method (both formats work for the function `export_topology`). This is 
useful for sharing topologies with other users who may be using a different version of BuildAMol and therefore could potentially have issues
with pickled objects, or for programmatic creation of new topologies (since especially XML is an easier format to work with than the CHARMM topology file format itself).

.. code-block:: python

    # export a CHARMM topology file to JSON
    top.to_json("top_all36_prot.json")

    # or to an XML file with the function
    bam.export_topology("top_all36_prot.xml", top)
    
    
Working with CHARMM objects
===========================

Naturally, the `CHARMMTopology` class includes methods to work with the parsed data, such as `has_linkage` or `get_linkage`. 
They also support adding new data via the the `add_linkage` method. For these as well (you guessed it) there are functional equivalents, which only work for the currently loaded default topology, however!


.. code-block:: python

    # check if a given topology has a linkage
    top.has_linkage("DISU") # use a specific topology
    # or
    bam.has_linkage("DISU") # use the default topology

    # get a linkage
    top.get_linkage("DISU") # use a specific topology
    # or
    bam.get_linkage("DISU") # use the default topology

    # add a linkage
    my_linkage = bam.linkage(...)
    top.add_linkage(my_linkage) # use a specific topology
    # or
    bam.add_linkage(my_linkage) # use the default topology



Setting default CHARMM objects
==============================

BuildAMol pre-loads a default CHARMM topology for convenience. Many methods that make use of 
this object such as the `attach` method of `Molecule` instances also accept arguments for custom topology objects.
For convenience, a custom object can be set as the default, however, using the `set_default_topology` function.

.. code-block:: python
    
    # set a custom topology object as the default
    bam.set_default_topology(top)


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

import buildamol.utils as utils
import buildamol.core.Linkage as Linkage
import buildamol.utils.defaults as _defaults

# ===================================================================
# Base Parser
# ===================================================================


def read_topology(filename: str, set_default: bool = True) -> "CHARMMTopology":
    """
    Make a CHARMMTopology from a CHARMM topology file, a JSON, an XML, or a pickle file.

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
    if filename.endswith(".pkl"):
        top = CHARMMTopology.load(filename)
    elif filename.endswith(".json"):
        top = CHARMMTopology.from_json(filename)
    elif filename.endswith(".xml"):
        top = CHARMMTopology.from_xml(filename)
    else:
        top = CHARMMTopology.from_file(filename)
    if set_default:
        set_default_topology(top)
    return top


def save_topology(filename: str, topology: "CHARMMTopology" = None):
    """
    Save a CHARMMTopology to a pickle file.

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


def export_topology(filename: str, topology: "CHARMMTopology" = None):
    """
    Export a CHARMM topology to a JSON or XML file.

    Parameters
    ----------
    filename: str
        The name of the topology file
    topology: CHARMMTopology
        The topology object. If `None`, the default topology is used.
    """
    if topology is None:
        topology = get_default_topology()
    if filename.endswith(".xml"):
        topology.to_xml(filename)
    elif filename.endswith(".json"):
        topology.to_json(filename)
    else:
        raise ValueError("Strange file extension. Please use .xml or .json")


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
            if not os.path.exists(_defaults.DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak"):
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


def restore_default_topology(overwrite: bool = True):
    """
    Restore the default CHARMMTopology object from the backup file

    Parameters
    ----------
    overwrite : bool
        If set to `True`, the backup is permanently set as the default again.
    """
    _defaults.__default_instances__["Topology"] = CHARMMTopology.load(
        _defaults.DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak"
    )
    if overwrite:
        os.rename(
            _defaults.DEFAULT_CHARMM_TOPOLOGY_FILE + ".bak",
            _defaults.DEFAULT_CHARMM_TOPOLOGY_FILE,
        )


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
    if patch.id is None:
        raise ValueError("The linkage must have an ID")
    get_default_topology().add_patch(patch)
    if overwrite:
        set_default_topology(get_default_topology(), overwrite=True)


add_linkage = add_patch


def get_patch(name: str):
    """
    Get a patch from the CHARMM topology file.

    Parameters
    ----------
    name: str
        The name/id of the patch

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
            The path to the pickle file to save in. By default, this will be
            the same filename as the one from which the data was loaded or
            parsed (adding the file-suffix `.pkl`)
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
    patches : list of Linkage
        The Linkages (patches, in the CHARMM nomenclature, prefixed with 'PRES') in the topology file
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
        self._dict["patches"] = {}
        self.id = id

    @property
    def patches(self):
        return list(self._dict["patches"].values())

    @property
    def linkages(self):
        return self.patches

    @classmethod
    def from_json(cls, filename: str) -> "CHARMMTopology":
        """
        Make a CHARMMTopology from a previously exported JSON file.

        Parameters
        ----------
        filename: str
            The path to the JSON file
        """
        _dict = utils.json.read(filename)
        new = cls()
        for link in _dict["patches"]:
            new.add_patch(Linkage._from_dict(link))
        new.id = _dict["id"]
        new._file = filename
        return new

    def to_json(self, filename: str = None):
        """
        Export the topology as JSON file.

        Parameters
        ----------
        filename: str
            The path to the JSON file to save in. By default, this will be
            the same filename as the one from which the data was loaded or
            parsed (adding the file-suffix `.json`)
        """
        if not filename:
            if not self._file:
                raise ValueError(
                    "No filename was given and no filename from a source file is available!"
                )
            filename = self._file
            if not filename.endswith(".json"):
                filename += ".json"
        utils.json.write_charmm_topology(self, filename)

    @classmethod
    def from_xml(cls, filename: str) -> "CHARMMTopology":
        """
        Make a CHARMMTopology from a previously exported XML file.

        Parameters
        ----------
        filename: str
            The path to the XML file
        """
        xml = utils.xml.read_xml(filename)
        new = cls()
        for patch in xml["patches"]:
            new.add_patch(Linkage._from_xml(patch))
        new.id = xml["id"]
        new._file = filename
        return new

    def to_xml(self, filename: str = None):
        """
        Export the topology as XML file.

        Parameters
        ----------
        filename: str
            The path to the XML file to save in. By default, this will be
            the same filename as the one from which the data was loaded or
            parsed (adding the file-suffix `.xml`)
        """
        if not filename:
            if not self._file:
                raise ValueError(
                    "No filename was given and no filename from a source file is available!"
                )
            filename = self._file
            if not filename.endswith(".xml"):
                filename += ".xml"
        xml = utils.xml.encode_topology(self)
        utils.xml.write_xml(filename, xml)

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
            lines = [line.strip() for line in lines]
            # lines = [
            #     line[0] for line in lines if line[0] != ""
            # ]  # get rid of all empty lines

        idx = 0
        while idx < len(lines):
            line = lines[idx]
            if line.startswith("PRES"):
                idx = self._parse_patch(lines, idx)
            idx += 1

        self._file = filename

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
        descr = lines[idx].split("!")
        if len(descr) == 2:
            descr = descr[1].strip()
        else:
            descr = None
        line = self._read_line(lines[idx])

        patch = Linkage(id=line[1], description=descr)

        idx += 1
        while idx < len(lines):
            line = self._read_line(lines[idx])
            if len(line) == 0:
                idx += 1
                continue
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

    def _vet_load(self, _data):
        """
        Checks that the data loaded from a file is valid
        """
        if not isinstance(_data, dict):
            raise TypeError("The file must contain a dictionary object")
        if "patches" not in _data.keys() or not isinstance(
            _data.get("patches", None), dict
        ):
            raise KeyError("Data must contain a 'patches' key with a dictionary value")

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

        # if isinstance(obj, utils.abstract.AbstractResidue):
        #     atom1 = obj.get_atom(atom1)
        #     atom2 = obj.get_atom(atom2)
        #     atom3 = obj.get_atom(atom3)
        #     atom4 = obj.get_atom(atom4)

        #     if atom1 is None or atom2 is None or atom3 is None or atom4 is None:
        #         warnings.warn(
        #             f"[ignoring line] Found an invalid internal coordinate in {line}"
        #         )
        #         return

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
set_default_topology(
    read_topology(_defaults.DEFAULT_CHARMM_TOPOLOGY_FILE, set_default=False)
)


__all__ = [
    "CHARMMTopology",
    "get_default_topology",
    "set_default_topology",
    "restore_default_topology",
    "has_patch",
    "available_patches",
    "available_linkages",
    "add_patch",
    "add_linkage",
    "has_linkage",
    "get_patch",
    "get_linkage",
    "read_topology",
    "export_topology",
    "save_topology",
]

if __name__ == "__main__":
    # _carbs = "/Users/noahhk/GIT/biobuild/support/toppar_charmm/carbohydrates.rtf"
    # _top = CHARMMTopology.from_file(_carbs)
    patches = "/Users/noahhk/GIT/biobuild/support/charmm_topology/patches.rtf"
    _top = CHARMMTopology.from_file(patches)
    print(_top)

    from buildamol.utils.defaults import DEFAULT_CHARMM_TOPOLOGY_FILE

    _save_to = "/Users/noahhk/GIT/biobuild/biobuild/resources/"
    _top.to_json(_save_to + os.path.basename(DEFAULT_CHARMM_TOPOLOGY_FILE))
