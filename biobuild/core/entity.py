"""
The base class for classes storing and manipulating biopython structures
This houses most of the essential functionality of the library for most users. 
The ``Molecule`` class adds additional features on top. 
"""

from copy import deepcopy
from typing import Union
import warnings

import Bio.PDB as bio
import numpy as np

import biobuild.core.Linkage as Linkage
import biobuild.structural as structural
import biobuild.utils as utils
import biobuild.graphs as graphs
import biobuild.resources as resources
import biobuild.core.base_classes as base_classes


class BaseEntity:
    """
    THe Base class for all classes that store and handle biopython structures, namely the Molecule class.

    Parameters
    ----------
    structure : Bio.PDB.Structure
        The biopython structure
    model : int
        The index of the model to use (default: 0)
    """

    def __init__(self, structure, model: int = 0):
        if isinstance(structure, bio.Structure.Structure):
            structure = base_classes.Structure.from_biopython(structure)
        self._base_struct = structure
        self._id = structure.id

        self._model = self._base_struct.child_list[model]
        if len(self._model.child_list) == 0:
            raise ValueError("The model is empty")

        self._bonds = []
        # self._locked_bonds = set()

        self._AtomGraph = graphs.AtomGraph(self.id, [])
        self._AtomGraph.add_nodes_from(self._model.get_atoms())

        # let the molecule store a patch to use for attaching other
        # molecules to it, or which Recipe to use by default for stitching
        # we use the same attribute for this since they are mutually exclusive and equivalent
        self._linkage = None

        self._working_chain = None
        self._root_atom = None
        self._attach_residue = None

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        id: str = None,
    ):
        """
        Read a Molecule from a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        new = cls(struct)
        bonds = utils.pdb.parse_connect_lines(filename)
        if len(bonds) != 0:
            for b in bonds:
                new.add_bond(*b)
        return new

    @classmethod
    def from_cif(cls, filename: str, id: str = None):
        """
        Load a Molecule from a CIF file

        Parameters
        ----------
        filename : str
            Path to the CIF file
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        try:
            struct = utils.defaults.get_default_instance(
                "bioMMCIFParser"
            ).get_structure(id, filename)
            new = cls(struct)
            bonds = utils.cif.parse_bond_table(filename)
            if len(bonds) != 0:
                for b in bonds:
                    new.add_bond(*b)
            return new
        except KeyError:
            try:
                c = resources.PDBECompounds.from_file(filename)
                if len(c) == 0:
                    raise ValueError(f"No compounds found in {filename}")
                return c.get(c.ids[0])
            except Exception as e:
                raise e

    @classmethod
    def from_json(cls, filename: str):
        """
        Make a Molecule from a JSON file

        Parameters
        ----------
        filename : str
            Path to the JSON file
        """
        _dict = utils.json.read(filename)
        return cls._from_dict(_dict)

    @classmethod
    def _from_dict(cls, _dict):
        """
        Make a Molecule from a JSON dictionary
        """
        _struct = base_classes.Structure(_dict["id"])
        structure_dict = _dict["structure"]

        model = base_classes.Model(structure_dict["model"]["id"])
        _struct.add(model)

        chain_dict = {}
        for chain in structure_dict["chains"]["id"]:
            chain = base_classes.Chain(chain)
            model.add(chain)
            chain_dict[chain.id] = chain

        residue_dict = {}
        for i in range(len(structure_dict["residues"]["serial"])):
            residue = base_classes.Residue(
                structure_dict["residues"]["name"][i],
                " ",
                structure_dict["residues"]["serial"][i],
            )
            parent = chain_dict[structure_dict["residues"]["parent"][i]]
            parent.add(residue)
            residue_dict[residue.serial_number] = residue

        for i in range(len(structure_dict["atoms"]["serial"])):
            atom = base_classes.Atom(
                id=structure_dict["atoms"]["id"][i],
                coord=np.array(structure_dict["atoms"]["coords_3d"][i]),
                serial_number=structure_dict["atoms"]["serial"][i],
                element=structure_dict["atoms"]["element"][i].upper(),
            )
            parent = residue_dict[structure_dict["atoms"]["parent"][i]]
            parent.add(atom)

        new = cls(_struct)
        for bond, order in zip(
            structure_dict["bonds"]["serial"], structure_dict["bonds"]["order"]
        ):
            for i in range(order):
                new.add_bond(bond[0], bond[1])
        return new

    @classmethod
    def from_pybel(cls, mol):
        """
        Load a Molecule from a Pybel molecule

        Parameters
        ----------
        mol : pybel.Molecule
            The Pybel molecule
        """
        conv = utils.convert.PybelBioPythonConverter()
        new = conv.pybel_molecule_to_biopython(mol)
        new = base_classes.Structure.from_biopython(new)
        new = cls(new)
        return new

    @classmethod
    def from_rdkit(cls, mol, id: str = None):
        """
        Load a Molecule from an RDKit molecule

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The RDKit molecule
        id : str
            The id of the Molecule. By default an id is inferred from the "_Name" property of the mol object (if present).
        """
        if id is None:
            if not mol.HasProp("_Name"):
                id = "unnamed"
            else:
                id = mol.GetProp("_Name")

        conv = utils.convert.RDKITBiopythonConverter()
        conv.rdkit_to_pdbio(mol)
        new = cls.from_pdb(conv.__fileio__, id=id)
        return new

    @classmethod
    def load(cls, filename: str):
        """
        Load a Molecule from a pickle file

        Parameters
        ----------
        filename : str
            Path to the file
        """
        obj = utils.load_pickle(filename)
        if obj.__class__.__name__ != cls.__name__:
            raise TypeError(
                f"Object loaded from {filename} is not a {cls.__name__} but a {type(obj)}"
            )
        return obj

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, new_id):
        self._id = new_id
        self._base_struct.id = new_id

    @property
    def structure(self):
        """
        The biopython structure
        """
        return self._base_struct

    @property
    def model(self):
        """
        The biopython model
        """
        return self._model

    @property
    def linkage(self):
        """
        The patch or recipe to use for attaching other molecules to this one
        """
        return self._linkage

    @property
    def patch(self):
        """
        The patch to use for attaching other molecules to this one (synonym for recipe)
        """
        return self._linkage

    @property
    def recipe(self):
        """
        The recipe to use for stitching other molecules to this one (synonym for patch)
        """
        return self._linkage

    @recipe.setter
    def recipe(self, value):
        # self._patch = value
        self.set_linkage(value)

    @patch.setter
    def patch(self, value):
        # self._patch = value
        self.set_linkage(value)

    @linkage.setter
    def linkage(self, value):
        # self._patch = value
        self.set_linkage(value)

    @property
    def bonds(self):
        """
        All bonds in the structure
        """
        return self._bonds
        # return list(self._AtomGraph.edges)

    @bonds.setter
    def bonds(self, value):
        if value is None or len(value) == 0:
            self._bonds.clear()
            self._AtomGraph.clear_edges()
        self._bonds = value
        self._AtomGraph.clear_edges()
        self._AtomGraph.add_edges_from(value)

    @property
    def locked_bonds(self):
        """
        All bonds that are locked and cannot be rotated around.
        """
        return self._AtomGraph._locked_edges

    @locked_bonds.setter
    def locked_bonds(self, value):
        if value is None or len(value) == 0:
            self._AtomGraph._locked_edges.clear()
        elif isinstance(value, set):
            self._AtomGraph._locked_edges = value
        else:
            raise TypeError("locked_bonds must be a set")

    @property
    def chains(self):
        """
        A sorted list of all chains in the structure
        """
        return sorted(self._model.get_chains(), key=lambda x: len(x.id))

    @property
    def residues(self):
        """
        A sorted list of all residues in the structure
        """
        return sorted(self._model.get_residues(), key=lambda x: x.id[1])

    @property
    def atoms(self):
        """
        A sorted list of all atoms in the structure
        """
        return sorted(self._AtomGraph.nodes)
        # return list(self._model.get_atoms())

    @property
    def atom_triplets(self):
        """
        Compute triplets of three consequtively bonded atoms
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_triplets(self.bonds)

    @property
    def atom_quartets(self):
        """
        Compute quartets of four consequtively bonded atoms
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_quartets(self.bonds)

    @property
    def angles(self):
        """
        Compute all angles of consecutively bonded atom triplets within the molecule.

        Returns
        -------
        angles : dict
            A dictionary of the form {atom_triplet: angle}
        """
        return {
            triplet: structural.compute_angle(*triplet)
            for triplet in structural.generate_triplets(self.bonds)
        }

    @property
    def dihedrals(self):
        """
        Compute all dihedrals of consecutively bonded atom quartets within the molecule.

        Returns
        -------
        dihedrals : dict
            A dictionary of the form {atom_quartet: dihedral}
        """
        return {
            quartet: structural.compute_dihedral(*quartet)
            for quartet in self.atom_quartets
        }

    @property
    def _chain(self):
        if self._working_chain is None:
            return self.chains[-1]
        else:
            return self._working_chain

    @_chain.setter
    def _chain(self, value):
        self._working_chain = value

    @property
    def root_atom(self):
        """
        The root atom of this molecule/scaffold at which it is attached to another molecule/scaffold
        """
        return self._root_atom

    @root_atom.setter
    def root_atom(self, value):
        if value is None:
            self._root_atom = None
        else:
            self._root_atom = self.get_atom(value)

    @property
    def root_residue(self):
        """
        The residue of the root atom
        """
        if self._root_atom:
            return self.root_atom.get_parent()

    @property
    def attach_residue(self):
        """
        The residue at which to attach other molecules to this one.
        """
        return self._attach_residue

    @attach_residue.setter
    def attach_residue(self, value):
        if value is None:
            self._attach_residue = None
        else:
            self._attach_residue = self.get_residue(value)

    def get_root(self) -> bio.Atom.Atom:
        """
        Get the root atom of the molecule. The root atom is the atom
        at which it is attached to another molecule.
        """
        return self.root_atom

    def set_root(self, atom):
        """
        Set the root atom of the molecule

        Parameters
        ----------
        atom : Atom or int or str or tuple
            The atom to be used as the root atom.
            This may be an Atom object, an atom serial number, an atom id (must be unique), or the full-id tuple.
        """
        self.root_atom = atom

    def get_linkage(self):
        """
        Get the linkage that is currently set as default attachment specication for this molecule
        """
        return self._linkage

    def set_linkage(
        self,
        link: Union[str, "Linkage.Linkage"] = None,
        _topology=None,
    ):
        """
        Set a linkage to be used for attaching other molecules to this one

        Parameters
        ----------
        link : str or Linkage
            The linkage to be used. Can be either a string with the name of a known Linkage in the loaded topology, or an instance of the Linkage class.
            If None is given, the currently loaded default linkage is removed.
        _topology
            The topology to use for referencing the link.
        """
        if link is None:
            self._linkage = None
            return

        if isinstance(link, str):
            if not _topology:
                _topology = resources.get_default_topology()
            self._linkage = _topology.get_patch(link)
        elif isinstance(link, Linkage.Linkage):
            self._linkage = link
        else:
            raise ValueError(f"Unknown linkage type {type(link)}")

    def save(self, filename: str):
        """
        Save the object to a pickle file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        utils.save_pickle(self, filename)

    def show(self, residue_graph: bool = False):
        """
        Open a browser window to view the molecule in 3D

        Parameters
        ----------
        residue_graph : bool
            If True, a residue graph is shown instead of the full structure.
        """
        v = self.draw(residue_graph)
        v.show()

    def nglview(self):
        """
        View the molecule in 3D through nglview
        """
        utils.visual.NglViewer(self).show()

    def draw(self, residue_graph: bool = False):
        """
        Prepare a view of the molecule in 3D
        but do not open a browser window.

        Parameters
        ----------
        residue_graph : bool
            If True, a residue graph is shown instead of the full structure.

        Returns
        -------
        viewer : MoleculeViewer3D
            The viewer object
        """
        if residue_graph:
            return utils.visual.MoleculeViewer3D(self.make_residue_graph())
        else:
            return utils.visual.MoleculeViewer3D(self)

    def vet(
        self, clash_range: tuple = (0.6, 1.7), angle_range: tuple = (90, 180)
    ) -> bool:
        """
        Vet the structural integrity of a molecule.
        This will return True if there are no clashes and all angles
        of adjacent atom triplets are within a tolerable range, False otherwise.

        Parameters
        ----------
        clash_range : tuple, optional
            The minimal and maximal allowed distances for two bonded atoms (in Angstrom).
        angle_range : tuple, optional
            The minimal and maximal allowed angles between tree adjacent bonded atoms (in degrees).

        Returns
        -------
        bool
            True if the structure is alright, False otherwise.
        """
        return structural.vet_structure(self, clash_range, angle_range)

    def copy(self):
        """
        Create a deepcopy of the molecule
        """
        new = deepcopy(self)
        new._base_struct._new_id()
        new._AtomGraph.clear()
        for model in new._base_struct.child_list:
            new._base_struct.child_dict.pop(model.get_id())
            model._new_id()
            new._base_struct.child_dict[model.get_id()] = model
            for chain in model.child_list:
                model.child_dict.pop(chain.get_id())
                chain._new_id()
                model.child_dict[chain.get_id()] = chain
                for residue in chain.child_list:
                    chain.child_dict.pop(residue.get_id())
                    residue._new_id()
                    chain.child_dict[residue.get_id()] = residue
                    for atom in residue.child_list:
                        residue.child_dict.pop(atom.get_id())
                        atom._new_id()
                        residue.child_dict[atom.get_id()] = atom
        new._AtomGraph.add_nodes_from(new.get_atoms())
        new._AtomGraph.add_edges_from(new.get_bonds())
        for b in new.get_bonds():
            new._AtomGraph.edges[b]["bond_order"] = (
                new._AtomGraph.edges[b].get("bond_order", 0) + 1
            )
        return new

    def get_attach_residue(self):
        """
        Get the residue that is used for attaching other molecules to this one.
        """
        return self._attach_residue

    def set_attach_residue(self, residue: Union[int, bio.Residue.Residue] = None):
        """
        Set the residue that is used for attaching other molecules to this one.

        Parameters
        ----------
        residue
            The residue to be used for attaching other molecules to this one
        """
        if residue is None:
            self._attach_residue = None
        else:
            residue = self.get_residue(residue)
            self._attach_residue = residue

    def rotate_descendants(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        angle: float,
    ):
        """
        Rotate all descendant atoms (atoms after atom2) of a bond.

        Parameters
        ----------
        atom1 : Union[str, int, bio.Atom.Atom]
            The first atom
        atom2 : Union[str, int, bio.Atom.Atom]
            The second atom (whose downstream neighbors are rotated)
        angle : float
            The angle to rotate by in degrees
        """
        self.rotate_around_bond(atom1, atom2, angle, descendants_only=True)

    def rotate_ancestors(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        angle: float,
    ):
        """
        Rotate all ancestor atoms (atoms before atom1) of a bond

        Parameters
        ----------
        atom1 : Union[str, int, bio.Atom.Atom]
            The first atom (whose upstream neighbors are rotated)
        atom2 : Union[str, int, bio.Atom.Atom]
            The second atom
        angle : float
            The angle to rotate by in degrees
        """
        self.rotate_around_bond(atom2, atom1, angle, descendants_only=True)

    def rotate_around_bond(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        angle: float,
        descendants_only: bool = False,
    ):
        """
        Rotate the structure around a bond

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        angle
            The angle to rotate by in degrees
        descendants_only
            Whether to only rotate the descendants of the bond, i.e. only atoms that come after atom2
            (sensible only for linear molecules, or bonds that are not part of a circular structure).
        
        Examples
        --------
        For a molecule starting as:
        ```
                     OH
                    /
        (1)CH3 --- CH 
                    \\
                    CH2 --- (2)CH3
        ```
        we can rotate around the bond `(1)CH3 --- CH` by 180° using

        >>> import numpy as np
        >>> angle = 180
        >>> mol.rotate_around_bond("(1)CH3", "CH", angle)

        and thus achieve the following:
        ```
                     CH2 --- (2)CH3
                    /
        (1)CH3 --- CH 
                    \\
                     OH
        ```

        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if (atom1, atom2) in self.locked_bonds:
            raise RuntimeError("Cannot rotate around a locked bond")

        angle = np.radians(angle)

        self._AtomGraph.rotate_around_edge(atom1, atom2, angle, descendants_only)

    def get_ancestors(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
    ):
        """
        Get the atoms upstream of a bond. This will return the set
        of all atoms that are connected before the bond atom1-atom2 in the direction of atom1,
        the selection can be reversed by reversing the order of atoms (atom2-atom1).

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom

        Returns
        -------
        set
            A set of atoms

        Examples
        --------
        For a molecule
        ```
                     OH
                    /
        (1)CH3 --- CH 
                    \\
                    CH2 --- (2)CH3
        ```
        >>> mol.get_ancestors("(1)CH3", "CH")
        set()
        >>> mol.get_ancestors("CH", "CH2")
        {"(1)CH3", "OH"}
        >>> mol.get_ancestors("CH2", "CH")
        {"(2)CH3"}
        """
        return self.get_descendants(atom2, atom1)

    def get_descendants(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
    ):
        """
        Get the atoms downstream of a bond. This will return the set
        of all atoms that are connected after the bond atom1-atom2 in the direction of atom2,
        the selection can be reversed by reversing the order of atoms (atom2-atom1).

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom

        Returns
        -------
        set
            A set of atoms

        Examples
        --------
        For a molecule
        ```
                     OH
                    /
        (1)CH3 --- CH 
                    \\
                    CH2 --- (2)CH3
        ```
        >>> mol.get_descendants("(1)CH3", "CH")
        {"OH", "CH2", "(2)CH3"}
        >>> mol.get_descendants("CH", "CH2")
        {"(2)CH3"}
        >>> mol.get_descendants("CH2", "CH")
        {"OH", "(1)CH3"}
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        return self._AtomGraph.get_descendants(atom1, atom2)

    def get_neighbors(
        self,
        atom: Union[int, str, tuple, bio.Atom.Atom],
        n: int = 1,
        mode: str = "upto",
    ):
        """
        Get the neighbors of an atom.

        Parameters
        ----------
        atom
            The atom
        n
            The number of bonds that may separate the atom from its neighbors.
        mode
            The mode to use. Can be "upto" or "at". If `upto`, all neighbors that are at most `n` bonds away
            are returned. If `at`, only neighbors that are exactly `n` bonds away are returned.

        Returns
        -------
        set
            A set of atoms


        Examples
        --------
        For a molecule
        ```
                     O --- (2)CH2
                    /         \\
        (1)CH3 --- CH          OH
                    \\
                    (1)CH2 --- (2)CH3
        ```
        >>> mol.get_neighbors("(2)CH2", n=1)
        {"O", "OH"}
        >>> mol.get_neighbors("(2)CH2", n=2, mode="upto")
        {"O", "OH", "CH"}
        >>> mol.get_neighbors("(2)CH2", n=2, mode="at")
        {"CH"}
        """
        atom = self.get_atom(atom)
        return self._AtomGraph.get_neighbors(atom, n, mode)

    def reindex(
        self, start_chainid: int = 1, start_resid: int = 1, start_atomid: int = 1
    ):
        """
        Reindex the atoms and residues in the structure.
        You can use this method if you made substantial changes
        to the molecule and want to be sure that there are no gaps in the
        atom and residue numbering.

        Parameters
        ----------
        start_chainid : int
            The starting chain id (default: 1=A, 2=B, ..., 26=Z, 27=AA, 28=AB, ...)
        start_resid : int
            The starting residue id
        start_atomid : int
            The starting atom id
        """
        cdx = start_chainid - 1
        rdx = start_resid
        adx = start_atomid

        chains = list(self.chains)
        for chain in chains:
            chain._id = utils.auxiliary.chain_id_maker(cdx)
            cdx += 1

            for residue in chain.child_list:
                # chain.child_dict.pop(residue.id)
                residue.serial_number = rdx
                # residue.id = (residue.id[0], rdx, *residue.id[2:])
                rdx += 1
                # chain.child_dict[residue._id] = residue

                for atom in residue.child_list:
                    atom.serial_number = adx
                    adx += 1

        # update the atom graph
        # let's see iff all breaks if we don't do this
        # self.update_atom_graph()

    def adjust_indexing(self, mol):
        """
        Adjust the indexing of a molecule to match the scaffold index

        Parameters
        ----------
        mol : Molecule
            The molecule to adjust the indexing of
        """
        cdx = len(self.chains)
        rdx = len(self.residues)
        adx = sum(1 for i in self._model.get_atoms())
        mol.reindex(cdx + 1, rdx + 1, adx + 1)

    def get_chains(self):
        return self._model.get_chains()

    def get_models(self):
        return self._base_struct.get_models()

    def get_residues(
        self,
        *residues: Union[int, str, tuple, bio.Residue.Residue],
        by: str = None,
        chain=None,
    ):
        """
        Get residues from the structure either based on their
        name, serial number or full_id.

        Parameters
        ----------
        residues
            The residues' id, seqid or full_id tuple. If None is passed, the iterator over all residues is returned.
        by : str
            The type of parameter to search for. Can be either 'name', 'seqid' or 'full_id'
            By default, this is inferred from the datatype of the residue parameter.
            If it is an integer, it is assumed to be the sequence identifying number,
            if it is a string, it is assumed to be the residue name and if it is a tuple, it is assumed
            to be the full_id.
        chain : str
            Further restrict to residues from a specific chain.

        Returns
        -------
        list or generator
            The residue(s)
        """
        if len(residues) == 0:
            if chain is not None:
                chain = self.get_chain(chain)
                return chain.get_residues()
            return self._model.get_residues()

        _residues = []
        for residue in residues:
            if isinstance(residue, bio.Residue.Residue):
                if residue in self.residues:
                    return residue
                else:
                    return self.get_residue(residue.id[1], by="seqid", chain=chain)

            if by is None:
                if isinstance(residue, int):
                    by = "seqid"
                elif isinstance(residue, str):
                    by = "name"
                elif isinstance(residue, tuple):
                    by = "full_id"
                else:
                    raise ValueError(
                        f"Cannot infer search parameter from residue query '{residue}', provide `by` manually: 'name', 'seqid' or 'full_id'"
                    )

            if by == "name":
                _residue = [
                    i for i in self._model.get_residues() if i.resname == residue
                ]
            elif by == "seqid":
                if residue < 0:
                    residue = len(self.residues) + residue + 1
                _residue = [i for i in self._model.get_residues() if i.id[1] == residue]
            elif by == "full_id":
                _residue = [
                    i for i in self._model.get_residues() if i.full_id == residue
                ]
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'name', 'seqid' or 'full_id'"
                )
            if chain is not None:
                chain = self.get_chain(chain)
                _residue = [i for i in _residue if i.get_parent() == chain]

            _residues.extend(_residue)

        return _residues

    def get_atoms(self, *atoms: Union[int, str, tuple], by: str = None) -> list:
        """
        Get one or more atoms from the structure either based on their
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues all of them are returned in a list.
        It is a safer option to use the full_id or serial number to retrieve a specific atom.
        If no search parameters are provided, the underlying iterator of the structure is returned.

        Parameters
        ----------
        atoms
            The atom id, serial number or full_id tuple,
            this supports multiple atoms to search for. However,
            only one type of parameter is supported per call. If left empty, the underlying generator is returned.

        by : str
            The type of parameter to search for. Can be either 'id', 'serial' or 'full_id'
            If None is given, the parameter is inferred from the datatype of the atoms argument
            'serial' in case of `int`, 'id' in case of `str`, `full_id` in case of a tuple.

        Returns
        -------
        atom : list or generator
            The atom(s)
        """
        if len(atoms) == 0:
            return self._model.get_atoms()

        if by is None:
            if isinstance(atoms[0], int):
                by = "serial"
            elif isinstance(atoms[0], str):
                by = "id"
            elif isinstance(atoms[0], tuple):
                by = "full_id"
            elif isinstance(atoms[0], bio.Atom.Atom):
                return atoms
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
                )

        if by == "id":
            atoms = [i for i in self._model.get_atoms() if i.id in atoms]
        elif by == "serial":
            atoms = [i for i in self._model.get_atoms() if i.serial_number in atoms]
        elif by == "full_id":
            atoms = [i for i in self._model.get_atoms() if i.full_id in atoms]
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
            )

        return atoms

    def get_atom(
        self,
        atom: Union[int, str, tuple],
        by: str = None,
        residue: Union[int, bio.Residue.Residue] = None,
    ):
        """
        Get an atom from the structure either based on its
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues, only the first one is returned.
        To get all atoms matching the criteria, use the `get_atoms` method.

        Parameters
        ----------
        atom
            The atom id, serial number or full_id tuple
        by : str
            The type of parameter to search for. Can be either 'id', 'serial' or 'full_id'
            Because this looks for one specific atom, this parameter can be inferred from the datatype
            of the atom parameter. If it is an integer, it is assumed to be the serial number,
            if it is a string, it is assumed to be the atom id and if it is a tuple, it is assumed
            to be the full_id.
        residue: int or Residue
            A specific residue to search in. If None, the entire structure is searched.

        Returns
        -------
        atom : bio.Atom.Atom
            The atom
        """
        if residue is not None:
            residue = self.get_residue(residue)
            atom_gen = residue.get_atoms
        else:
            atom_gen = self._model.get_atoms

        if isinstance(atom, base_classes.Atom):
            if atom in atom_gen():
                return atom
            else:
                return self.get_atom(atom.full_id, by="full_id")

        if by is None:
            if isinstance(atom, (int, np.int64)):
                by = "serial"
            elif isinstance(atom, str):
                by = "id"
            elif isinstance(atom, tuple):
                by = "full_id"
            else:
                raise ValueError(
                    "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
                )

        if by == "id":
            _atom = (i for i in atom_gen() if i.id == atom)
        elif by == "serial":
            _atom = (i for i in atom_gen() if i.serial_number == atom)
        elif by == "full_id":
            _atom = (i for i in atom_gen() if i.full_id == atom)
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'id', 'serial' or 'full_id'"
            )

        return next(_atom, None)

    def get_bonds(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom, bio.Residue.Residue] = None,
        atom2: Union[int, str, tuple, bio.Atom.Atom] = None,
        either_way: bool = True,
        residue_internal: bool = True,
    ):
        """
        Get one or multiple bonds from the molecule. If only one atom is provided, all bonds
        that are connected to that atom are returned.

        Parameters
        ----------
        atom1
            The atom id, serial number or full_id tuple of the first atom.
            This may also be a residue, in which case all bonds between atoms in that residue are returned.
        atom2
            The atom id, serial number or full_id tuple of the second atom
        either_way : bool
            If True, the order of the atoms does not matter, if False, the order of the atoms
            does matter. By setting this to false, it is possible to also search for bonds that have
            a specific atom in position 1 or 2 depending on which argument was set, while leaving the other input as none.
        residue_internal : bool
            If True, only bonds where both atoms are in the given residue (if `atom1` is a residue) are returned.
            If False, all bonds where either atom is in the given residue are returned.

        Returns
        -------
        bond : list
            The bond(s). If no input is given, all bonds are returned.
        """
        if atom1 is None and atom2 is None:
            return self.bonds

        if isinstance(atom1, bio.Residue.Residue):
            if residue_internal:
                return [
                    i
                    for i in self.bonds
                    if i[0].get_parent() is atom1 and i[1].get_parent() is atom1
                ]
            else:
                return [
                    i
                    for i in self.bonds
                    if i[0].get_parent() is atom1 or i[1].get_parent() is atom1
                ]

        if atom1:
            atom1 = self.get_atom(atom1)
        if atom2:
            atom2 = self.get_atom(atom2)
        return self._get_bonds(atom1, atom2, either_way)

    def get_residue(
        self,
        residue: Union[int, str, tuple, bio.Residue.Residue],
        by: str = None,
        chain=None,
    ):
        """
        Get a residue from the structure either based on its
        name, serial number or full_id. Note, if multiple residues match the requested criteria,
        for instance there are multiple 'MAN' from different chains, only the first one is returned.

        Parameters
        ----------
        residue
            The residue id, seqid or full_id tuple
        by : str
            The type of parameter to search for. Can be either 'name', 'seqid' or 'full_id'
            By default, this is inferred from the datatype of the residue parameter.
            If it is an integer, it is assumed to be the sequence identifying number,
            if it is a string, it is assumed to be the residue name and if it is a tuple, it is assumed
            to be the full_id.
        chain : str
            Further restrict to a residue from a specific chain.

        Returns
        -------
        residue : bio.Residue.Residue
            The residue
        """
        if isinstance(residue, bio.Residue.Residue):
            if residue in self.residues:
                return residue
            else:
                return self.get_residue(residue.id[1], by="seqid", chain=chain)

        if by is None:
            if isinstance(residue, int):
                by = "seqid"
            elif isinstance(residue, str):
                by = "name"
            elif isinstance(residue, tuple):
                by = "full_id"
            else:
                raise ValueError(
                    f"Cannot infer search parameter from residue query '{residue}', provide `by` manually: 'name', 'seqid' or 'full_id'"
                )

        if by == "name":
            _residue = (i for i in self._model.get_residues() if i.resname == residue)
        elif by == "seqid":
            if residue < 0:
                residue = len(self.residues) + residue + 1
            _residue = (i for i in self._model.get_residues() if i.id[1] == residue)
        elif by == "full_id":
            _residue = (i for i in self._model.get_residues() if i.full_id == residue)
        else:
            raise ValueError(
                "Unknown search parameter, must be either 'name', 'seqid' or 'full_id'"
            )
        if chain is not None:
            chain = self.get_chain(chain)
            _residue = (i for i in _residue if i.get_parent() == chain)
        return next(_residue, None)

    def get_chain(self, chain: str):
        """
        Get a chain from the structure either based on its
        name.

        Parameters
        ----------
        chain
            The chain id
        Returns
        -------
        chain : bio.Chain.Chain
            The chain
        """
        if isinstance(chain, bio.Chain.Chain):
            if chain in self.chains:
                return chain
            else:
                return self.get_chain(chain.id)
        return next((i for i in self.chains if i.id == chain), None)

    def add_residues(
        self,
        *residues: bio.Residue.Residue,
        adjust_seqid: bool = True,
        _copy: bool = False,
    ):
        """
        Add residues to the structure

        Parameters
        ----------
        residues : bio.Residue.Residue
            The residues to add
        adjust_seqid : bool
            If True, the seqid of the residues is adjusted to
            match the current number of residues in the structure
            (i.e. a new residue can be given seqid 1, and it will be adjusted
            to the correct value of 3 if there are already two other residues in the molecule).
        _copy : bool
            If True, the residues are copied before adding them to the molecule.
            This is useful if you want to add the same residue to multiple molecules, while leaving
            them and their original parent structures intakt.
        """
        rdx = len(self.residues)
        adx = sum(1 for i in self._model.get_atoms())
        for residue in residues:
            p = residue.get_parent()
            if p:
                residue.detach_parent()
            if _copy and p is not None:
                r = residue.copy()
                residue.set_parent(p)
                residue = r

            rdx += 1
            if adjust_seqid and residue.id[1] != rdx:
                residue.serial_number = rdx
                # the above line with serial_number works for the
                # biobuild class derivatives, the below line
                # is for the original biopython classes
                # residue.id = (residue.id[0], rdx, *residue.id[2:])

            self._chain.add(residue)

            for atom in residue.child_list:
                if adjust_seqid:
                    adx += 1
                    atom.set_serial_number(adx)
                self._AtomGraph.add_node(atom)

    def remove_residues(self, *residues: Union[int, bio.Residue.Residue]) -> list:
        """
        Remove residues from the structure

        Parameters
        ----------
        residues : int or bio.Residue.Residue
            The residues to remove, either the object itself or its seqid

        Returns
        -------
        list
            The removed residues
        """
        _residues = []
        for residue in residues:
            if isinstance(residue, int):
                residue = self._chain.child_list[residue - 1]

            for atom in residue.child_list:
                self._AtomGraph.remove_node(atom)
                self._purge_bonds(atom)

            # keep the memory of the parent in the residue that is removed...
            chain = residue.get_parent()
            chain.detach_child(residue.get_id())
            residue.set_parent(chain)

            _residues.append(residue)
        return _residues

    def rename_residue(self, residue: Union[int, bio.Residue.Residue], name: str):
        """
        Rename a residue

        Parameters
        ----------
        residue : int or bio.Residue.Residue
            The residue to rename, either the object itself or its seqid
        name : str
            The new name
        """
        if isinstance(residue, int):
            residue = self._chain.child_list[residue - 1]
        residue.resname = name

    def rename_atom(
        self,
        atom: Union[int, bio.Atom.Atom],
        name: str,
        residue: Union[int, bio.Residue.Residue] = None,
    ):
        """
        Rename an atom

        Parameters
        ----------
        atom : int or bio.Atom.Atom
            The atom to rename, either the object itself or its serial number
        name : str
            The new name (id)
        residue : int or bio.Residue.Residue
            The residue to which the atom belongs, either the object itself or its seqid.
            Useful when giving a possibly redundant id as identifier in multi-residue molecules.
        """
        atom = self.get_atom(atom, residue=residue)
        # p = atom.get_parent()
        # _old = atom.id
        atom.id = name
        atom.name = name
        # if p:
        #     p.child_dict.pop(_old)
        #     p.child_dict[name] = atom

    def add_atoms(self, *atoms: bio.Atom.Atom, residue=None, _copy: bool = False):
        """
        Add atoms to the structure. This will automatically adjust the atom's serial number to
        fit into the structure.

        Parameters
        ----------
        atoms : bio.Atom.Atom
            The atoms to add
        residue : int or str
            The residue to which the atoms should be added,
            this may be either the seqid or the residue name,
            if None the atoms are added to the last residue.
            Note, that if multiple identically named residues
            are present, the first one is chosen, so using the
            seqid is a safer option!
        _copy : bool
            If True, the atoms are copied and then added to the structure.
            This will leave the original atoms (and their parent structures) untouched.
        """
        if residue is not None:
            if isinstance(residue, int):
                residue = self._chain.child_list[residue - 1]
            elif isinstance(residue, str):
                residue = next(i for i in self._chain.child_list if i.resname == i)
            target = residue
        else:
            target = self._chain.child_list[-1]

        _max_serial = sum(1 for i in self._model.get_atoms())
        for atom in atoms:
            if _copy:
                atom = deepcopy(atom)
            _max_serial += 1
            atom.set_serial_number(_max_serial)
            target.add(atom)
            self._AtomGraph.add_node(atom)

    def remove_atoms(self, *atoms: Union[int, str, tuple, bio.Atom.Atom]) -> list:
        """
        Remove one or more atoms from the structure

        Parameters
        ----------
        atoms
            The atoms to remove, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        Returns
        -------
        list
            The removed atoms
        """
        _atoms = []
        for atom in atoms:
            _atom = self.get_atoms(atom)
            if len(_atom) == 0:
                continue
            for atom in _atom:
                self._AtomGraph.remove_node(atom)
                self._purge_bonds(atom)
                p = atom.get_parent()
                if p:
                    p.detach_child(atom.get_id())
                    atom.set_parent(
                        p
                    )  # this is necessary to avoid a bug where atoms remain longer in the memory atoms list than they should
                _atoms.append(atom)

        # reindex the atoms
        adx = 0
        for atom in self._model.get_atoms():
            adx += 1
            atom.serial_number = adx

        return _atoms

    def add_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
    ):
        """
        Add a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        if not atom1:
            raise ValueError("Atom1 not found!")
        if not atom2:
            raise ValueError("Atom2 not found!")

        # if (atom1, atom2) not in self._bonds:
        self._bonds.append((atom1, atom2))
        if not self._AtomGraph.has_edge(atom1, atom2):
            self._AtomGraph.add_edge(atom1, atom2, bond_order=1)
        else:
            self._AtomGraph.edges[atom1, atom2]["bond_order"] = (
                self._AtomGraph.edges[atom1, atom2].get("bond_order", 0) + 1
            )

    def add_bonds(self, *bonds):
        """
        Add multiple bonds at once

        Parameters
        ----------
        bonds
            The bonds to add, each bond is a tuple of two atoms.
            Each atom may be specified directly (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        for bond in bonds:
            self.add_bond(*bond)

    def remove_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        either_way: bool = True,
    ):
        """
        Remove a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        either_way : bool
            If True, the bond will be removed in both directions, i.e. if the bond is (1, 2)
            it will be removed if either (1, 2) or (2, 1) is provided.
        """
        if either_way:
            self.remove_bond(atom1, atom2, either_way=False)
            self.remove_bond(atom2, atom1, either_way=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        self._remove_bond(atom1, atom2)

    def purge_bonds(self, atom: Union[int, str, bio.Atom.Atom] = None):
        """
        Remove all bonds connected to an atom

        Parameters
        ----------
        atom
            The atom to remove the bonds from, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms. If None, all bonds
            are removed.
        """
        if atom is None:
            self.bonds = []
            return

        atom = self.get_atom(atom)
        self._purge_bonds(atom)

    def lock_all(self, both_ways: bool = True):
        """
        Lock all bonds in the structure

        Parameters
        ----------
        both_ways : bool
            If True, the bond is locked in both directions
            i.e. atom1 --- atom2 direction will be unavailable for
            rotation as well as atom2 --- atom1 direction as well.
        """
        self.locked_bonds = set(self.bonds)
        if both_ways:
            self.locked_bonds.update(b[::-1] for b in self.bonds)

    def unlock_all(self):
        """
        Unlock all bonds in the structure
        """
        self._AtomGraph.unlock_all()

    def lock_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        both_ways: bool = False,
    ):
        """
        Lock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        both_ways : bool
            If True, the bond is locked in both directions. By default the bond is only locked in the specified direction.
        """
        if both_ways:
            self.lock_bond(atom2, atom1, both_ways=False)
            self.lock_bond(atom2, atom1, both_ways=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._AtomGraph.lock_edge(atom1, atom2)

    def unlock_bond(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
        both_ways: bool = False,
    ):
        """
        Unlock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        both_ways : bool
            If True, the bond is unlocked in both directions. By default the bond is only unlocked in the specified direction.
        """
        if both_ways:
            self.unlock_bond(atom1, atom2, both_ways=False)
            self.unlock_bond(atom2, atom1, both_ways=False)
            return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if (atom1, atom2) in self._AtomGraph.edges:
            self._AtomGraph.unlock_edge(atom1, atom2)
        if both_ways and (atom2, atom1) in self._AtomGraph.edges:
            self._AtomGraph.unlock_edge(atom2, atom1)

    def is_locked(
        self,
        atom1: Union[int, str, tuple, bio.Atom.Atom],
        atom2: Union[int, str, tuple, bio.Atom.Atom],
    ):
        """
        Check if a bond is locked

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        Returns
        -------
        bool
            True if the bond is locked, False otherwise
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)

        bond = (atom1, atom2)
        return bond in self.locked_bonds

    def infer_bonds(
        self, max_bond_length: float = None, restrict_residues: bool = True
    ) -> list:
        """
        Infer bonds between atoms in the structure

        Parameters
        ----------
        max_bond_length : float
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.
        restrict_residues : bool
            Whether to restrict bonds to only those in the same residue.
            If False, bonds between atoms in different residues are also inferred.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.infer_bonds(
            self._base_struct, max_bond_length, restrict_residues
        )
        to_add = [b for b in bonds if b not in self.bonds]
        self._bonds.extend(to_add)
        self._AtomGraph.add_edges_from(to_add, bond_order=1)
        return bonds

    def get_residue_connections(
        self, residue_a=None, residue_b=None, triplet: bool = True
    ):
        """
        Get bonds between atoms that connect different residues in the structure
        This method is different from `infer_residue_connections` in that it works
        with the already present bonds in the molecule instead of computing new ones.

        Parameters
        ----------
        residue_a, residue_b : Union[int, str, tuple, bio.Residue.Residue]
            The residues to consider. If None, all residues are considered.
            Otherwise, only between the specified residues are considered.
        triplet : bool
            Whether to include bonds between atoms that are in the same residue
            but neighboring a bond that connects different residues. This is useful
            for residues that have a side chain that is connected to the main chain.
            This is mostly useful if you intend to use the returned list for some purpose,
            because the additionally returned bonds are already present in the structure
            from inference or standard-bond applying and therefore do not actually add any
            particular information to the Molecule object itself.

        Returns
        -------
        list
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = (i for i in self.bonds if i[0].get_parent() != i[1].get_parent())

        if residue_a is not None and residue_b is None:
            residue_a = self.get_residue(residue_a)
            bonds = (
                i
                for i in bonds
                if i[0].get_parent() is residue_a or i[1].get_parent() is residue_a
            )
        elif residue_b is not None and residue_a is None:
            residue_b = self.get_residue(residue_b)
            bonds = (
                i
                for i in bonds
                if i[0].get_parent() is residue_b or i[1].get_parent() is residue_b
            )
        elif residue_a is not None and residue_b is not None:
            residue_a, residue_b = self.get_residue(residue_a), self.get_residue(
                residue_b
            )
            bonds = (
                i
                for i in bonds
                if (i[0].get_parent() is residue_a and i[1].get_parent() is residue_b)
                or (i[1].get_parent() is residue_a and i[0].get_parent() is residue_b)
            )
        if triplet:
            bonds = self._make_bond_triplets(bonds)
        return list(bonds)

    def _make_bond_triplets(self, bonds) -> set:
        """
        Make triplets of bonds for bonds that are connected to either side of a bond.

        Parameters
        ----------
        bonds : iterable
            A set of tuples of atom pairs that are bonded

        Returns
        -------
        set
            A set of tuples of atom pairs that are bonded and connect to either partner of the original bonds
            These are added to the original bonds set.
        """
        bonds = set(bonds)
        _new = set()
        for atom1, atom2 in bonds:
            neighs = self.get_neighbors(atom1)
            neighs.remove(atom2)
            neighs -= set(i for i in neighs if i.element == "H")
            if len(neighs) == 1:
                neigh = neighs.pop()
                if neigh.get_parent() is atom1.get_parent():
                    _new.add((atom1, neigh))
                continue

            neighs = self.get_neighbors(atom2)
            neighs.remove(atom1)
            neighs -= set(i for i in neighs if i.element == "H")
            if len(neighs) == 1:
                neigh = neighs.pop()
                if neigh.get_parent() is atom2.get_parent():
                    _new.add((atom2, neigh))
                continue

        bonds.update(_new)
        return bonds

    def infer_residue_connections(
        self, bond_length: Union[float, tuple] = None, triplet: bool = True
    ) -> list:
        """
        Infer bonds between atoms that connect different residues in the structure

        Parameters
        ----------
        bond_length : float or tuple
            If a float is given, the maximum distance between atoms to consider them bonded.
            If a tuple, the minimal and maximal distance between atoms. If None, the default value is min 0.8 Angstrom, max 1.6 Angstroms.
        triplet : bool
            Whether to include bonds between atoms that are in the same residue
            but neighboring a bond that connects different residues. This is useful
            for residues that have a side chain that is connected to the main chain.
            This is mostly useful if you intend to use the returned list for some purpose,
            because the additionally returned bonds are already present in the structure 
            from inference or standard-bond applying and therefore do not actually add any 
            particular information to the Molecule object itself.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded and considered residue connections.

        Examples
        --------
        For a molecule with the following structure:
        ```     
             connection -->  OA     OB --- H
                            /  \\   /
              (1)CA --- (2)CA   (1)CB 
               /         \\        \\
           (6)CA         (3)CA    (2)CB --- (3)CB
               \\         /
             (5)CA --- (4)CA
        ``` 
        The circular residue A and linear residue B are connected by a bond between
        `(1)CA` and the oxygen `OA` and `(1)CB`. By default, because OA originally is associated
        with residue A, only the bond `OA --- (1)CB` is returned. However, if `triplet=True`,
        the bond `OA --- (1)CA` is also returned, because the entire connecting "bridge" between residues
        A and B spans either bond around `OA`.
        >>> mol.infer_residue_connections(triplet=False)
        [("OA", "(1)CB")]
        >>> mol.infer_residue_connections(triplet=True)
        [("OA", "(1)CB"), ("OA", "(2)CA")]
        """
        bonds = structural.infer_residue_connections(
            self._base_struct, bond_length, triplet
        )
        self._bonds.extend(b for b in bonds if b not in self._bonds)
        self._AtomGraph.add_edges_from(bonds)
        return bonds

    def apply_standard_bonds(self, _topology=None) -> list:
        """
        Get the standard bonds for the structure

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.apply_reference_bonds(self._base_struct, _topology)
        to_add = [b for b in bonds if b not in self.bonds]
        self._bonds.extend(to_add)
        self._AtomGraph.add_edges_from(to_add)
        return bonds

    def autolabel(self):
        """
        Automatically label atoms in the structure to match the CHARMM force field
        atom nomenclature. This is useful if you want to use some pre-generated
        PDB file that may have used a different labelling scheme for atoms.

        Note
        ----
        The labels are infererred and therefore may occasionally not be "correct".
        It is advisable to check the labels after using this method.
        """
        self = structural.autolabel(self)

    def relabel_hydrogens(self):
        """
        Relabel hydrogen atoms in the structure to match the standard labelling according
        to the CHARMM force field. This is useful if you want to use some pre-generated
        PDB file that may have used a different labelling scheme for atoms.
        """
        structural.relabel_hydrogens(self)

    def compute_angle(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        atom3: Union[str, int, bio.Atom.Atom],
    ):
        """
        Compute the angle between three atoms where atom2 is the middle atom.

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        atom3
            The third atom

        Returns
        -------
        float
            The angle in degrees
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        atom3 = self.get_atom(atom3)
        return structural.compute_angle(atom1, atom2, atom3)

    def compute_dihedral(
        self,
        atom1: Union[str, int, bio.Atom.Atom],
        atom2: Union[str, int, bio.Atom.Atom],
        atom3: Union[str, int, bio.Atom.Atom],
        atom4: Union[str, int, bio.Atom.Atom],
    ):
        """
        Compute the dihedral angle between four atoms

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        atom3
            The third atom
        atom4
            The fourth atom

        Returns
        -------
        float
            The dihedral angle in degrees
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        atom3 = self.get_atom(atom3)
        atom4 = self.get_atom(atom4)
        return structural.compute_dihedral(atom1, atom2, atom3, atom4)

    def make_atom_graph(self, locked: bool = True):
        """
        Generate an AtomGraph for the Molecule

        Parameters
        ----------
        locked : bool
            If True, the graph will also migrate the information on any locked bonds into the graph.

        Returns
        -------
        AtomGraph
            The generated graph
        """
        graph = deepcopy(self._AtomGraph)
        if not locked:
            graph.unlock_all()
        return graph

    def update_atom_graph(self):
        """
        Generate a new up-to-date `AtomGraph` after any manual changes were done to the Molecule's underlying biopython structure.
        """
        self._AtomGraph.clear()
        self._AtomGraph.add_nodes_from(self.get_atoms())
        self._AtomGraph.add_edges_from(self._bonds)

    def make_residue_graph(self, detailed: bool = False, locked: bool = True):
        """
        Generate a ResidueGraph for the Molecule

        Parameters
        ----------
        detailed : bool
            If True the graph will include the residues and all atoms that form bonds
            connecting different residues. If False, the graph will only include the residues
            and their connections without factual bonds between any existing atoms.
        locked : bool
            If True, the graph will also migrate the information on any locked bonds into the graph.
            This is only relevant if detailed is True.
        """
        graph = graphs.ResidueGraph.from_molecule(self, detailed, locked)
        return graph

    def to_pdb(self, filename: str, symmetric: bool = True):
        """
        Write the molecule to a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        symmetric : bool
            If True, bonds are written symmetrically - i.e. if atom A is bonded to atom B, then atom B is also bonded to atom A,
            and both atoms will get an entry in the "CONECT" section. If False, only one of the atoms will get an entry in the
            "CONECT" section.
        """
        utils.pdb.write_pdb(self, filename, symmetric=symmetric)

        # io = bio.PDBIO()
        # io.set_structure(self._base_struct)
        # io.save(filename)
        # utils.pdb.write_connect_lines(self, filename)
        # with open(filename, "r") as f:
        #     content = f.read()
        # content = utils.remove_nonprintable(content)
        # with open(filename, "w") as f:
        #     f.write(content)

    def to_cif(self, filename: str):
        """
        Write the molecule to a CIF file

        Parameters
        ----------
        filename : str
            Path to the CIF file
        """
        io = bio.MMCIFIO()
        io.set_structure(self._base_struct.to_biopython())
        io.save(filename)
        utils.cif.write_bond_table(self, filename)
        with open(filename, "r") as f:
            content = f.read()
        content = utils.remove_nonprintable(content)
        with open(filename, "w") as f:
            f.write(content)

    def to_json(
        self,
        filename: str,
        type: str = None,
        names: list = None,
        identifiers: list = None,
        one_letter_code: str = None,
        three_letter_code: str = None,
    ):
        """
        Write the molecule to a JSON file

        Parameters
        ----------
        filename : str
            Path to the JSON file
        type : str
            The type of the molecule to be written to the JSON file (e.g. "protein", "ligand", etc.).
        names : list
            A list of names of the molecules to be written to the JSON file.
        identifiers : list
            A list of identifiers of the molecules to be written to the JSON file (e.g. SMILES, InChI, etc.).
        one_letter_code : str
            A one-letter code for the molecule to be written to the JSON file.
        three_letter_code : str
            A three-letter code for the molecule to be written to the JSON file.
        """
        utils.json.write_molecule(
            self, filename, type, names, identifiers, one_letter_code, three_letter_code
        )

    def to_pybel(self):
        """
        Convert the molecule to a Pybel molecule

        Returns
        -------
        pybel.Molecule
            The Pybel molecule
        """
        conv = utils.convert.PybelBioPythonConverter()
        mol = conv.biobuild_to_pybel(self)
        mol.title = self.id
        return mol

    def to_rdkit(self):
        """
        Convert the molecule to an RDKit molecule

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The RDKit molecule
        """
        conv = utils.convert.RDKITBiopythonConverter()
        conv.molecule_to_pdbio(self)
        mol = conv.pdbio_to_rdkit()
        mol.SetProp("_Name", self.id)
        return mol

    def to_biopython(self):
        """
        Convert the molecule to a Biopython structure

        Returns
        -------
        Bio.PDB.Structure.Structure
            The Biopython structure
        """
        return self._base_struct.to_biopython()

    def infer_missing_atoms(self, _topology=None, _compounds=None):
        """
        Infer missing atoms in the structure based on a reference topology or compounds database.
        By default, if a residue is not available in the topology, the compounds database is used.

        Parameters
        ----------
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        _compounds
            A specific compounds object to use for referencing.
            If None, the default compounds object is used.
        """
        structural.fill_missing_atoms(self._base_struct, _topology, _compounds)
        for atom in self._base_struct.get_atoms():
            if atom not in self._AtomGraph:
                self._AtomGraph.add_node(atom)
        self.infer_bonds()

    def _get_bonds(
        self,
        atom1,
        atom2,
        either_way: bool = True,
    ):
        """
        The core function of `get_bonds` which expects atoms to be provided as Atom objects
        """
        if atom1 and atom2:
            if either_way:
                return [
                    i
                    for i in self.bonds
                    if (atom1 is i[0] and atom2 is i[1])
                    or (atom1 is i[1] and atom2 is i[0])
                    # if (atom1.full_id == i[0].full_id and atom2.full_id == i[1].full_id)
                    # or (atom1.full_id == i[1].full_id and atom2.full_id == i[0].full_id)
                ]
            else:
                return [
                    i
                    for i in self.bonds
                    if (atom1 is i[0] and atom2 is i[1])
                    # if (atom1.full_id == i[0].full_id and atom2.full_id == i[1].full_id)
                ]
        elif atom1:
            if either_way:
                return [i for i in self.bonds if (atom1 is i[0] or atom1 is i[1])]
            else:
                return [i for i in self.bonds if atom1 is i[0]]
        elif atom2:
            if either_way:
                return [i for i in self.bonds if (atom2 is i[0] or atom2 is i[1])]
            else:
                return [i for i in self.bonds if atom2 is i[1]]
        else:
            raise ValueError("No atom provided")

    def _purge_bonds(self, atom):
        """
        The core function of `purge_bonds` which expects atoms to be provided as Atom objects.
        """
        # the full_id thing seems to prevent some memory leaky-ness
        bonds = [
            i
            for i in self._bonds
            if atom is i[0] or atom is i[1]
            # if atom.full_id == i[0].full_id or atom.full_id == i[1].full_id
        ]
        for bond in bonds:
            self._remove_bond(*bond)

    def _remove_bond(self, atom1, atom2, either_way: bool = False):
        """
        The core function of `remove_bond` which expects atoms to be provided as Atom objects.
        """
        b = (atom1, atom2)
        if b in self._bonds:
            self._bonds.remove(b)
        if self._AtomGraph.has_edge(atom1, atom2):
            if self._AtomGraph[atom1][atom2].get("bond_order", 1) == 1:
                self._AtomGraph.remove_edge(atom1, atom2)
            else:
                self._AtomGraph[atom1][atom2]["bond_order"] = (
                    self._AtomGraph[atom1][atom2].get("bond_order", 1) - 1
                )
        self.locked_bonds.discard(b)

        if either_way:
            if b[::-1] in self._bonds:
                self._bonds.discard(b[::-1])
            self.locked_bonds.discard(b[::-1])
            if self._AtomGraph.has_edge(atom2, atom1):
                self._AtomGraph.remove_edge(atom2, atom1)

    def _remove_atoms(self, *atoms):
        """
        The core alternative of `remove_atoms` which expects atoms to be provided as Atom objects.
        """
        for atom in atoms:
            self._AtomGraph.remove_node(atom)
            self._purge_bonds(atom)
            p = atom.get_parent()
            p.detach_child(atom.get_id())
            atom.set_parent(p)

        # reindex the atoms
        adx = 0
        for atom in self._model.get_atoms():
            adx += 1
            atom.serial_number = adx

    def _direct_bonds(
        self,
        bonds: list,
        by: str = "resid",
        direct_connections: set = None,
        save: bool = True,
    ):
        """
        Sort a given list of bonds such that the first atom is the "earlier" atom
        in the bond.

        Parameters
        ----------
        bonds : list
            A list of tuples of atom pairs that are bonded and connect different residues
        by : str
            The attribute to sort by. Can be either "serial", "resid" or "root".
            In the case of "serial", the bonds are sorted by the serial number of the first atom.
            In the case of "resid", the bonds are sorted by the residue number of the first atom.
            In this case, bonds connecting atoms from the same residue are sorted by the serial number of the first atom.
            In the case of "root" the bonds are sorted based on the graph distances to the root atom,
            provided that the root atom is set (otherwise the atom with serial 1 is used).
        direct_connections : set
            A set of atom pairs that are bonded and connect different residues. Only used if by == "resid".
        save : bool
            Whether to save the sorted bonds as the new bonds in the Molecule.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded and connect different residues
        """
        if by == "serial":
            directed = [
                bond if bond[0].serial_number < bond[1].serial_number else bond[::-1]
                for bond in bonds
            ]
        elif by == "root":
            root = self._root_atom
            if root is None:
                root = self.atoms[0]
            directed = self._AtomGraph.direct_edges(root, bonds)
        elif by == "resid":
            if not direct_connections:
                raise ValueError("direct_connections must be provided if by == 'resid'")
            directed = [
                bond if not should_invert(bond, direct_connections) else bond[::-1]
                for bond in bonds
            ]
        if save:
            for old, new in zip(bonds, directed):
                self.remove_bond(*old)
                self.add_bond(*new)
        return directed

    def __mod__(self, patch):
        """
        Add a patch to the molecule using the % operator (i.e. mol % patch)
        """
        self.set_linkage(patch)
        return self

    def __xor__(self, atom):
        """
        Set the root atom for the molecule using the ^ operator (i.e. mol ^ atom)
        """
        self.set_root(atom)
        return self

    def __matmul__(self, residue):
        """
        Set the residue at which the molecule should be attached to another molecule
        using the @ operator (i.e. mol @ 1, for residue 1)
        """
        self.set_attach_residue(residue)
        return self


def should_invert(bond, direct_connecting_atoms):
    """
    Check if a given bond should be inverted during bond direction

    Parameters
    ----------
    bond : tuple
        A tuple of two atoms that are bonded
    direct_connecting_atoms : set
        A set of atoms that directly participate in bonds connecting different residues

    Returns
    -------
    bool
        Whether the bond should be inverted
    """
    atom1, atom2 = bond
    if atom1.parent.id[1] > atom2.parent.id[1]:
        return True
    elif atom1.parent.id[1] == atom2.parent.id[1]:
        if atom1 in direct_connecting_atoms:
            return True
        elif atom2 not in direct_connecting_atoms:
            return False
        elif atom1.serial_number > atom2.serial_number:
            return True
    return False


if __name__ == "__main__":
    # f = "/Users/noahhk/GIT/biobuild/support/examples/4tvp.prot.pdb"
    # e = BaseEntity.from_pdb(f)
    # e.infer_bonds(restrict_residues=False)
    # e.get_residue_connections()
    # pass

    BaseEntity.from_json("man.json").show()

    # import biobuild as bb

    # man = bb.Molecule.from_compound("MAN")
    # man.repeat(3, "14bb")
    # cs = man.get_residue_connections()

    # v = utils.visual.MoleculeViewer3D(man)
    # for c in cs:
    #     v.draw_vector(
    #         None,
    #         c[0].coord,
    #         1.1 * (c[1].coord - c[0].coord),
    #         color="limegreen",
    #     )
    # v.show()
    # pass
