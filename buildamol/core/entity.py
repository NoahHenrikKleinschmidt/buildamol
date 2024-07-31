"""
The base class for classes storing and manipulating molecular structures
This houses most of the essential functionality of the library for most users. 
The ``Molecule`` class adds additional features on top. 
"""

from copy import deepcopy
from typing import Union
import warnings

import Bio.PDB as bio
import numpy as np

import buildamol.base_classes as base_classes
import buildamol.core.Linkage as Linkage
import buildamol.structural as structural
import buildamol.utils as utils
import buildamol.graphs as graphs
import buildamol.resources as resources


class BaseEntity:
    """
    The Base class for all classes that store and handle molecular structures.
    This class is not meant to be used directly but serves as the base for the Molecule class.

    Parameters
    ----------
    structure : Structure or Bio.PDB.Structure
        A BuildAMol or Biopython structure
    model : int
        The index of the model to use (default: 0)
    """

    __slots__ = (
        "_base_struct",
        "_model",
        "_id",
        "_bonds",
        "_AtomGraph",
        "_linkage",
        "_working_chain",
        "_root_atom",
        "_attach_residue",
    )

    def __init__(self, structure, model: int = 0):
        if not isinstance(structure, base_classes.Structure):
            if isinstance(structure, bio.Structure.Structure):
                structure = base_classes.Structure.from_biopython(structure)
            else:
                raise TypeError(
                    f"structure must be a Structure object, got {type(structure)}"
                )
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
        model: int = 0,
        has_atom_ids: bool = True,
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
        model : int
            The index of the model to use (default: 0)
        has_atom_ids : bool
            If the PDB file provides no atom ids, set this to False in order to autolabel the atoms.
        """
        if id is None:
            id = utils.filename_to_id(filename)
        f = open(filename)
        content = f.read()
        f.close()
        return cls._from_pdb_string(content, id=id)

    @classmethod
    def _from_pdb_string(cls, string, id: str = None):
        """
        Read a Molecule from a PDB string

        Parameters
        ----------
        string : str
            The PDB string
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        """
        content = string.split("\n")
        atoms = utils.pdb._parse_atom_lines(content)

        if len(atoms) == 1:
            atoms[0] = atoms[-1]
        if 0 not in atoms:
            atoms[0] = next(iter(atoms.values()))
        atoms.pop(-1)

        structure = base_classes.Structure(id)
        chains = {}
        residues = {}
        for model, _atoms in atoms.items():
            model = base_classes.Model(model)
            structure.add(model)
            chains.clear()
            residues.clear()
            for atom_info in _atoms:
                if atom_info["chain"] not in chains:
                    chains[atom_info["chain"]] = base_classes.Chain(atom_info["chain"])
                    model.add(chains[atom_info["chain"]])

                res_seq = (atom_info["chain"], atom_info["res_seq"])
                if res_seq not in residues:
                    residues[res_seq] = base_classes.Residue(
                        atom_info["residue"], " ", atom_info["res_seq"]
                    )
                    chains[atom_info["chain"]].add(residues[res_seq])

                atom = base_classes.Atom.new(
                    atom_info["id"],
                    generate_id=False,
                    altloc=atom_info["alt_loc"],
                    serial_number=atom_info["serial"],
                    coord=(atom_info["x"], atom_info["y"], atom_info["z"]),
                    occupancy=atom_info["occ"],
                    pqr_charge=eval(f"{atom_info['charge']}+0"),
                    element=atom_info["element"],
                )
                residues[res_seq].add(atom)

        new = cls(structure)
        bonds = utils.pdb._parse_connect_lines(content)
        if len(bonds) != 0:
            atom_mapping = {a.serial_number: a for a in new.get_atoms()}
            for a1, a2, o in bonds:
                new._set_bond(atom_mapping[a1], atom_mapping[a2], o)
            del atom_mapping
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
    def from_xml(cls, filename: str):
        """
        Make a Molecule from an XML file

        Parameters
        ----------
        filename : str
            Path to the XML file
        """
        tree = utils.xml.read_xml(filename)
        _struct = base_classes.Structure(tree.attributes["id"])
        _model = base_classes.Model(0)
        _struct.add(_model)
        _model._atom_index_mapping = {}
        for chain in tree.get_child("structure").children:
            _chain = base_classes.Chain(chain.get_attribute("id"))
            _model.add(_chain)
            for residue in chain.children:
                _residue = base_classes.Residue(
                    residue.get_attribute("name"),
                    " ",
                    residue.get_attribute("serial", int),
                )
                _chain.add(_residue)
                for atom in residue.children:
                    _atom = base_classes.Atom(
                        atom.get_attribute("id"),
                        coord=np.zeros(3),
                        serial_number=atom.get_attribute("serial", int),
                        element=atom.get_attribute("element"),
                    )
                    for attr in atom.attributes:
                        if attr not in ["id", "serial", "element"]:
                            setattr(_atom, attr, atom.get_attribute(attr))
                    _residue.add(_atom)
                    _model._atom_index_mapping[_atom.serial_number] = _atom
        new = cls(_struct)

        for bond in tree.get_child("connectivity").children:
            atom1 = new._model._atom_index_mapping[bond.get_attribute("atom1", int)]
            atom2 = new._model._atom_index_mapping[bond.get_attribute("atom2", int)]
            new._set_bond(
                atom1,
                atom2,
                bond.get_attribute("order", int),
            )

        for model in tree.get_child("coordinates").children:
            if model.get_attribute("id", int) != 0:
                new_model = _model.copy()
                new_model.id = model.get_attribute("id", int)
                new.structure.add(new_model)
                new._model = new_model

            for atom in model.children:
                _atom = new._model._atom_index_mapping[
                    atom.get_attribute("serial", int)
                ]
                _atom.coord = np.array(
                    [
                        atom.get_attribute("x", float),
                        atom.get_attribute("y", float),
                        atom.get_attribute("z", float),
                    ]
                )

        new.set_model(0)
        for model in new.get_models():
            del model._atom_index_mapping
        return new

    @classmethod
    def from_molfile(cls, filename: str):
        """
        Make a Molecule from a molfile

        Parameters
        ----------
        filename : str
            Path to the molfile
        """
        rdmol = utils.sdmol.read_mol(filename)
        return cls.from_rdkit(rdmol)

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
        atom_mapping = {a.serial_number: a for a in new.get_atoms()}
        for bond, order in zip(
            structure_dict["bonds"]["serial"], structure_dict["bonds"]["order"]
        ):
            for i in range(order):
                new._set_bond(atom_mapping[bond[0]], atom_mapping[bond[1]], order)

        del atom_mapping
        return new

    @classmethod
    def from_openmm(cls, topology, positions):
        """
        Load a Molecule from an OpenMM topology and positions

        Parameters
        ----------
        topology : simtk.openmm.app.Topology
            The OpenMM topology
        positions : simtk.unit.Quantity
            The OpenMM positions
        """
        conv = utils.convert.OpenMMBioPythonConverter()
        conv._openmm_to_pdbio(topology, positions)
        new = cls.from_pdb(conv.__fileio__)

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
        atom_mapping = {a.serial_number: a for a in new.get_atoms()}
        for i in range(mol.OBMol.NumBonds()):
            bond = mol.OBMol.GetBondById(i)
            a = atom_mapping[bond.GetBeginAtomIdx()]
            b = atom_mapping[bond.GetEndAtomIdx()]
            new._set_bond(a, b, bond.GetBondOrder())
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
                id = "Untitled"
            else:
                id = mol.GetProp("_Name")

        conv = utils.convert.RDKITBiopythonConverter()
        conv._rdkit_to_pdbio(mol)
        new = cls.from_pdb(conv.__fileio__, id=id)
        return new

    @classmethod
    def from_stk(cls, obj):
        """
        Load a Molecule from an stk ConstructedMolecule

        Parameters
        ----------
        obj : stk.ConstructedMolecule
            The stk ConstructedMolecule
        """
        conv = utils.convert.STKBuildAMolConverter()
        conv.stk_to_pdbio(obj)
        new = cls.from_pdb(conv.__fileio__)
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
        The buildamol base-structure
        """
        return self._base_struct

    @property
    def model(self):
        """
        The working model of the structure
        """
        return self._model

    @property
    def models(self):
        """
        A list of all models in the base-structure
        """
        return list(self._base_struct.get_models())

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
        else:
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
    def _chains(self):
        """
        Get the list of chains as returned by the get_chains() method
        """
        return list(self._model.get_chains())

    @property
    def residues(self):
        """
        A sorted list of all residues in the structure
        """
        return sorted(self._model.get_residues(), key=lambda x: x.id[1])

    @property
    def _residues(self):
        """
        Get the list of residues as returned by the get_residues() method
        """
        return list(self._model.get_residues())

    @property
    def atoms(self):
        """
        A sorted list of all atoms in the structure
        """
        return sorted(self._AtomGraph.nodes)
        # return list(self._model.get_atoms())

    @property
    def _atoms(self):
        """
        Get the list of atoms as returned by the get_atoms() method
        """
        return list(self.get_atoms())

    @property
    def center_of_mass(self):
        """
        The center of mass of the molecule
        """
        return structural.center_of_gravity(
            np.array([a.mass for a in self.atoms]),
            np.array([a.coord for a in self.atoms]),
        )

    @property
    def center_of_geometry(self):
        """
        The center of geometry of the molecule
        """
        return np.array([a.coord for a in self.atoms]).mean(axis=0)

    @property
    def mass(self):
        """
        The total mass of the molecule
        """
        return sum(a.mass for a in self.atoms)

    def get_atom_triplets(self):
        """
        Compute triplets of three consequtively bonded atoms
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return structural.compute_triplets(self.bonds)

    def get_atom_quartets(self) -> list:
        """
        Compute quartets of four consequtively bonded atoms

        Returns
        -------
        atom_quartets : list
            A list of atom quartets
        """
        if len(self.bonds) == 0:
            warnings.warn("No bonds found (yet), returning empty list")
            return []
        return list(structural.compute_quartets(self.bonds))

    def compute_angles(self):
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

    def compute_dihedrals(self):
        """
        Compute all dihedrals of consecutively bonded atom quartets within the molecule.

        Returns
        -------
        dihedrals : dict
            A dictionary of the form {atom_quartet: dihedral}
        """
        return {
            quartet: structural.compute_dihedral(*quartet)
            for quartet in self.get_atom_quartets()
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

    def get_root(self) -> base_classes.Atom:
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
        return self

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
            return self

        if isinstance(link, str):
            if not _topology:
                _topology = resources.get_default_topology()
            self._linkage = _topology.get_patch(link)
        elif isinstance(link, Linkage.Linkage):
            self._linkage = link
        else:
            raise ValueError(f"Unknown linkage type {type(link)}")
        return self

    def save(self, filename: str):
        """
        Save the object to a pickle file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        """
        utils.save_pickle(self, filename)

    def nglview(self):
        """
        View the molecule in 3D through nglview
        """
        return utils.visual.NglViewer(self)

    def py3dmol(self, style: str = "stick", color: str = None, size: tuple = None):
        """
        View the molecule in 3D through py3Dmol

        Parameters
        ----------
        style : str
            The style to use for the visualization. Can be
            "line", "stick", "sphere", "cartoon", "surface", or "label"
        color : str
            A specific color to use for the visualization
        size : tuple
            The size of the view as a tuple of (width, height) in pixels.
        """
        size = size or (600, 500)
        viewer = utils.visual.Py3DmolViewer(self, *size)
        if color:
            color = {"color": color}
        else:
            color = {}
        viewer.style = {style: color}
        viewer.view.setStyle(viewer.style)
        return viewer

    def chem2dview(self, *args, **kwargs):
        """
        View the molecule in 2D through RDKit
        """
        return utils.visual.Chem2DViewer(self)

    draw2d = chem2dview

    def show2d(self, *args, **kwargs):
        """
        View the molecule in 2D through RDKit
        """
        viewer = self.draw2d(*args, **kwargs)
        viewer.show()

    def plotly(
        self,
        residue_graph: bool = False,
        atoms: bool = True,
        line_color: str = "black",
    ):
        """
        Prepare a view of the molecule in 3D using Plotly
        but do not open a browser window.

        Parameters
        ----------
        residue_graph : bool
            If True, a residue graph is shown instead of the full structure.
        atoms : bool
            Whether to draw the atoms (default: True)
        line_color : str
            The color of the lines connecting the atoms

        Returns
        -------
        viewer : MoleculeViewer3D
            The viewer object
        """
        if residue_graph:
            return self.make_residue_graph().draw()
        else:
            v = utils.visual.MoleculeViewer3D()
            v.bond_color = line_color
            v.link(self)
            v.setup(atoms)

            return v

    def draw(self, *args, **kwargs):
        backend = utils.visual.DEFAULT_BACKEND
        return getattr(self, backend)(*args, **kwargs)

    draw3d = draw

    def show(self, *args, **kwargs):
        self.draw(*args, **kwargs).show()

    show3d = show

    # def vet(
    #     self, clash_range: tuple = (0.7, 1.7), angle_range: tuple = (90, 180)
    # ) -> bool:
    #     """
    #     Vet the structural integrity of a molecule.
    #     This will return True if there are no clashes and all angles
    #     of adjacent atom triplets are within a tolerable range, False otherwise.

    #     Parameters
    #     ----------
    #     clash_range : tuple, optional
    #         The minimal and maximal allowed distances for two bonded atoms (in Angstrom).
    #         The minimal distance is also used for non-bonded atoms.

    #     angle_range : tuple, optional
    #         The minimal and maximal allowed angles between tree adjacent bonded atoms (in degrees).

    #     Returns
    #     -------
    #     bool
    #         True if the structure is alright, False otherwise.
    #     """
    #     return structural.vet_structure(self, clash_range, angle_range)

    def find_clashes(
        self,
        clash_threshold: float = 1.0,
        ignore_hydrogens: bool = True,
        coarse_precheck: bool = True,
    ) -> list:
        """
        Find all clashes in the molecule.

        Parameters
        ----------
        clash_threshold : float, optional
            The minimal allowed distance between two atoms (in Angstrom).
        ignore_hydrogens : bool, optional
            Whether to ignore clashes with hydrogen atoms (default: True)
        coarse_precheck : bool, optional
            If set to True a coarse-grained pre-screening on residue-level is done
            to speed up the computation. This may cause the sytem to overlook clashes if
            individual residues are particularly large, however (e.g. lipids with long carbon chains).

        Returns
        -------
        list
            A list of tuples of atoms that clash.
        """
        return [
            i
            for i in structural.find_clashes_between(
                self, self, clash_threshold, ignore_hydrogens, coarse_precheck
            )
        ]

    def find_clashes_with(
        self,
        other,
        clash_threshold: float = 1.0,
        ignore_hydrogens: bool = True,
        coarse_precheck: bool = True,
    ) -> list:
        """
        Find all clashes between this molecule and another one.

        Parameters
        ----------
        other : Molecule
            The other molecule to compare with
        clash_threshold : float, optional
            The minimal allowed distance between two atoms (in Angstrom).
        ignore_hydrogens : bool, optional
            Whether to ignore clashes with hydrogen atoms (default: True)
        coarse_precheck : bool, optional
            If set to True a coarse-grained pre-screening on residue-level is done
            to speed up the computation. This may cause the sytem to overlook clashes if
            individual residues are particularly large, however (e.g. lipids with long carbon chains).

        Returns
        -------
        list
            A list of tuples of atoms that clash.
        """
        return [
            i
            for i in structural.find_clashes_between(
                self, other, clash_threshold, ignore_hydrogens, coarse_precheck
            )
        ]

    def count_clashes(
        self,
        clash_threshold: float = 1.0,
        ignore_hydrogens: bool = True,
        coarse_precheck: bool = True,
    ) -> int:
        """
        Count all clashes in the molecule.

        Parameters
        ----------
        clash_threshold : float, optional
            The minimal allowed distance between two atoms (in Angstrom).
        ignore_hydrogens : bool, optional
            Whether to ignore clashes with hydrogen atoms (default: True)
        coarse_precheck : bool, optional
            If set to True a coarse-grained pre-screening on residue-level is done
            to speed up the computation. This may cause the sytem to overlook clashes if
            individual residues are particularly large, however (e.g. lipids with long carbon chains).

        Returns
        -------
        int
            The number of clashes.
        """
        return sum(
            1
            for i in structural.find_clashes_between(
                self, self, clash_threshold, ignore_hydrogens, coarse_precheck
            )
        )

    def has_clashes(
        self,
        clash_threshold: float = 1.0,
        ignore_hydrogens: bool = True,
        coarse_precheck: bool = True,
    ) -> bool:
        """
        Check if the molecule has any clashes.

        Parameters
        ----------
        clash_threshold : float, optional
            The minimal allowed distance between two atoms (in Angstrom).
        ignore_hydrogens : bool, optional
            Whether to ignore clashes with hydrogen atoms (default: True)
        coarse_precheck : bool, optional
            If set to True a coarse-grained pre-screening on residue-level is done
            to speed up the computation. This may cause the sytem to overlook clashes if
            individual residues are particularly large, however (e.g. lipids with long carbon chains).

        Returns
        -------
        bool
            True if there are clashes, False otherwise.
        """
        return (
            next(
                structural.find_clashes_between(
                    self, self, clash_threshold, ignore_hydrogens, coarse_precheck
                ),
                None,
            )
            is not None
        )

    def copy(self, n: int = 1) -> list:
        """
        Create one or multiple deepcopy of the molecule

        Parameters
        ----------
        n : int, optional
            The number of copies to make, by default 1

        Returns
        -------
        Molecule or list
            The copied molecule(s)
        """
        if n > 1:
            return [self.copy() for i in range(n)]
        else:
            new = deepcopy(self)
            new._base_struct._new_id()
            new._AtomGraph.clear()

            new._base_struct.child_dict.clear()
            for model in new._base_struct.child_list:
                model._new_id()
                new._base_struct.child_dict[model.get_id()] = model

                model.child_dict.clear()
                for chain in model.child_list:
                    chain._new_id()
                    model.child_dict[chain.get_id()] = chain

                    chain.child_dict.clear()
                    for residue in chain.child_list:
                        residue._new_id()
                        chain.child_dict[residue.get_id()] = residue

                        residue.child_dict.clear()
                        for atom in residue.child_list:
                            atom._new_id()
                            residue.child_dict[atom.get_id()] = atom

            new._AtomGraph.add_nodes_from(new.get_atoms())
            new._AtomGraph.add_edges_from(new.get_bonds())
            for b in new.get_bonds():
                # I don't think this is necessary. Add again if it causes problems...
                # b.atom1 = new.get_atom(b.atom1.serial_number)
                # b.atom2 = new.get_atom(b.atom2.serial_number)
                new._AtomGraph.edges[b]["bond_order"] = b.order
                new._AtomGraph.edges[b]["bond_obj"] = b
            return new

    def merge(self, other, adjust_indexing: bool = True):
        """
        Merge another molecule into this one. This will simply add all chains, residues, and atoms of the other molecule to this one.
        It will NOT perform any kind of geometrical alignment or anything like that.

        Parameters
        ----------
        other : Molecule
            The other molecule to merge into this one
        adjust_indexing : bool
            Whether to adjust the indexing of the atoms and residues in the merged molecule
        """
        if adjust_indexing:
            self.adjust_indexing(other)

        self.add_chains(*other.chains)
        self._set_bonds(*other.get_bonds())
        return self

    def clear(self):
        """
        Clear the molecule of all models, chains, residues, and atoms.
        """
        self._base_struct.child_dict.clear()
        self._base_struct.child_list.clear()
        self._AtomGraph.clear()
        self._AtomGraph.clear_cache()
        self._bonds.clear()

    def squash(self, chain_id: str = "A", resname: str = "UNK"):
        """
        Turn the entire molecule into a single chain with a single residue.
        """
        chain = base_classes.Chain(chain_id)
        residue = base_classes.Residue(resname, " ", 1)
        chain.add(residue)
        for atom in self.get_atoms():
            residue.add(atom)

        self._model.child_dict.clear()
        self._model.child_list.clear()
        self.add_chains(chain, adjust_seqid=False)
        self.reindex()
        return self

    def squash_chains(self, chain_id: str = "A"):
        """
        Turn all chains of the molecule into a single chain but preserve the residues.
        """
        chain = base_classes.Chain(chain_id)
        for residue in self.get_residues():
            chain.add(residue)
        self._model.child_dict.clear()
        self._model.child_list.clear()
        self.add_chains(chain, adjust_seqid=False)
        self.reindex()
        return self

    def get_attach_residue(self):
        """
        Get the residue that is used for attaching other molecules to this one.
        """
        return self._attach_residue

    def set_attach_residue(self, residue: Union[int, base_classes.Residue] = None):
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

    def move_to(self, pos: np.ndarray):
        """
        Move the molecule to a specific position in 3D space

        Parameters
        ----------
        pos : np.ndarray
            The position to move the molecule to. This will be the new center of geometry.
        """
        vector = pos - self.center_of_geometry
        self.move(vector)
        return self

    place = move_to

    def move(self, vector: np.ndarray):
        """
        Move the molecule in 3D space

        Parameters
        ----------
        vector : np.ndarray
            The vector to move the molecule by
        """
        for atom in self.get_atoms():
            atom.move(vector)
        return self

    def rotate(
        self,
        angle: float,
        axis: np.ndarray,
        center: np.ndarray = None,
        angle_is_degrees: bool = True,
    ):
        """
        Rotate the molecule around an axis

        Parameters
        ----------
        angle : float
            The angle to rotate by
        axis : np.ndarray or str
            The axis to rotate around. This must be a unit vector.
            Alternatively, it may be one of the strings "x", "y", or "z" to rotate around the respective axes.
        center : np.ndarray
            The center of the rotation. By default the center of geometry is used to achieve relative rotations (i.e. without translation). Use "absolute" if you want to rotate around the literal axes.
        angle_is_degrees : bool
            Whether the angle is given in degrees (default) or radians
        """
        if isinstance(axis, (list, tuple)):
            axis = np.array(axis)
        if isinstance(axis, np.ndarray):
            axis = axis / np.linalg.norm(axis)
        elif isinstance(axis, str):
            axis = axis.lower().strip()
            sign = None
            if len(axis) == 2:
                sign, axis = axis[0], axis[1]
            if axis == "x":
                axis = structural.x_axis
            elif axis == "y":
                axis = structural.y_axis
            elif axis == "z":
                axis = structural.z_axis
            else:
                raise ValueError(f"Unknown axis input: {axis=}")
            if sign == "-":
                axis = -axis
        else:
            raise ValueError(f"Unknown axis type {type(axis)}")

        if center is None:
            center = self.center_of_geometry
        elif isinstance(center, str) and center.lower() == "absolute":
            center = np.zeros(3)
        elif not isinstance(center, np.ndarray):
            center = np.array(center)

        if angle_is_degrees:
            angle = np.radians(angle)

        structural.rotate_molecule(self, angle, axis, center)

        return self

    def flip(self, plane_vector: np.ndarray, center: np.ndarray = None):
        """
        Flip the molecule around an axis

        Parameters
        ----------
        plane_vector : np.ndarray or str
            The vector defining the plane to flip around. This must be a unit vector.
            It may also be one of the strings "xy", "xz", or "yz" to flip around the respective planes.
        center : np.ndarray
            The center of the flip
        """
        if not isinstance(plane_vector, np.ndarray):
            if isinstance(plane_vector, str):
                if plane_vector.lower() == "xy":
                    plane_vector = np.array([0, 0, 1])
                elif plane_vector.lower() == "xz":
                    plane_vector = np.array([0, 1, 0])
                elif plane_vector.lower() == "yz":
                    plane_vector = np.array([1, 0, 0])
            else:
                plane_vector = np.array(plane_vector)

        structural.flip_molecule(self, plane_vector, center)

        return self

    def transpose(
        self,
        vector: np.ndarray,
        angle: float,
        axis: np.ndarray,
        center: np.ndarray = None,
        angle_is_degrees: bool = True,
    ):
        """
        Transpose the molecule in 3D space

        Parameters
        ----------
        vector : np.ndarray
            The vector to move the molecule by
        angle : float
            The angle to rotate by
        axis : np.ndarray
            The axis to rotate around. This must be a unit vector.
        center : np.ndarray
            The center of the rotation
        angle_is_degrees : bool
            Whether the angle is given in degrees (default) or radians
        """
        if not isinstance(vector, np.ndarray):
            vector = np.array(vector)
        if not isinstance(axis, np.ndarray):
            axis = np.array(axis)
        if center is not None and not isinstance(center, np.ndarray):
            center = np.array(center)

        coords = np.array([a.coord for a in self.get_atoms()])
        if center is None:
            center = np.zeros(3)
        coords -= center

        if angle_is_degrees:
            angle = np.radians(angle)

        new_coords = structural.rotate_coords(coords=coords, angle=angle, axis=axis)
        new_coords += center

        for atom, coord in zip(self.get_atoms(), new_coords):
            atom.coord = coord + vector

        return self

    def superimpose_to_bond(
        self,
        ref_bond: Union[tuple, base_classes.Bond],
        other_bond: Union[tuple, base_classes.Bond],
    ):
        """
        Superimpose the molecule to another molecule based on two bonds.
        This will move this molecule so that the atoms in ref_bond are superimposed to the atoms in other_bond.

        Parameters
        ----------
        ref_bond : tuple or Bond
            The bond to reference in this molecule
        other_bond : tuple or Bond
            The bond to superimpose to in the other molecule
        """
        if not isinstance(ref_bond[0], np.ndarray):
            ref_bond = (ref_bond[0].get_coord(), ref_bond[1].get_coord())
        if not isinstance(other_bond[0], np.ndarray):
            other_bond = (other_bond[0].get_coord(), other_bond[1].get_coord())

        new_coords = structural.superimpose_points(
            self.get_coords(), ref_bond, other_bond
        )
        for adx, atom in enumerate(self.get_atoms()):
            atom.coord = new_coords[adx]

        return self

    def superimpose_to_pair(self, pair1, pair2):
        """
        Superimpose the molecule to another molecule based on two atom pairs (they do not need to be bonded).
        This will move this molecule so that the atoms in pair1 are superimposed to the atoms in pair2.

        Parameters
        ----------
        pair1 : tuple
            The pair to superimpose in this molecule. These may either be Atom objects or any input which can be used to get atoms in this molecule.
        pair2 : tuple
            The pair to superimpose to. These must be either Atom objects or arbitrary coordinates (np.ndarray).
        """
        if not isinstance(pair1[0], np.ndarray):
            if isinstance(pair1[0], (int, str)):
                pair1 = (self.get_atom(a) for a in pair1)
            pair1 = tuple(a.get_coord() for a in pair1)
        if not isinstance(pair2[0], np.ndarray):
            pair2 = tuple(a.get_coord() for a in pair2)

        if len(pair1) != 2 or len(pair2) != 2:
            raise ValueError("Both pairs must have exactly 2 elements")

        new_coords = structural.superimpose_points(self.get_coords(), pair1, pair2)
        for adx, atom in enumerate(self.get_atoms()):
            atom.coord = new_coords[adx]

        return self

    def superimpose_to_atom(
        self,
        ref_atom: Union[base_classes.Atom, int, str],
        other_atom: Union[base_classes.Atom, np.ndarray],
    ):
        """
        Superimpose the molecule to another molecule based on two atoms.
        This will move this molecule so that the atom in ref_atom is superimposed to the atom in other_atom.

        Parameters
        ----------
        ref_atom : Atom or int or str
            The atom to superimpose in this molecule
        other_atom : Atom or np.ndarray
            The atom to superimpose to in the other molecule or an arbitrary coordinate
        """
        if not isinstance(other_atom, np.ndarray):
            other_atom = other_atom.get_coord()
        if not isinstance(ref_atom, np.ndarray):
            if isinstance(ref_atom, (int, str)):
                ref_atom = self.get_atom(ref_atom)
            ref_atom = ref_atom.get_coord()

        vec = other_atom - ref_atom
        self.move(vec)
        return self

    def superimpose_to_triplet(self, ref_triplet: tuple, other_triplet: tuple):
        """
        Superimpose the molecule to another molecule based on two atom triplets.
        This will move this molecule so that the atoms in ref_triplet are superimposed to the atoms in other_triplet.

        Parameters
        ----------
        ref_triplet : tuple
            The triplet to superimpose to. These may either be Atom objects or any input which can be used to get atoms in this molecule.
        other_triplet : tuple
            The triplet to superimpose from.. These must be either Atom objects or arbitrary coordinates (np.ndarray).
        """
        if not isinstance(ref_triplet[0], np.ndarray):
            if isinstance(ref_triplet[0], (int, str)):
                ref_triplet = (self.get_atom(a) for a in ref_triplet)
            ref_triplet = tuple(a.get_coord() for a in ref_triplet)
        if not isinstance(other_triplet[0], np.ndarray):
            other_triplet = tuple(a.get_coord() for a in other_triplet)

        if len(ref_triplet) != 3 or len(other_triplet) != 3:
            raise ValueError("Both triplets must have exactly 3 elements")

        new_coords = structural.superimpose_points(
            self.get_coords(), ref_triplet, other_triplet
        )
        for adx, atom in enumerate(self.get_atoms()):
            atom.coord = new_coords[adx]

        return self

    def superimpose_to_residue(self, ref_residue, other_residue):
        """
        Superimpose the molecule to another molecule based on two residues.
        This will move this molecule so that the residues are superimposed.

        Parameters
        ----------
        ref_residue : Residue or int or str
            The residue to superimpose to in this molecule
        other_residue : Residue
            The residue to superimpose to in the other molecule
        """
        if isinstance(ref_residue, (int, str)):
            ref_residue = self.get_residue(ref_residue)
        if not isinstance(other_residue, base_classes.Residue):
            raise ValueError("other_residue must be a Residue object")

        ref_atoms = tuple(ref_residue.get_atoms())
        ref_triplet = []
        other_triplet = []
        for r in ref_atoms:
            o = other_residue.get_atom(r.id)
            if o is not None:
                ref_triplet.append(r)
                other_triplet.append(o)
            if len(ref_triplet) == 3:
                break
        if len(ref_triplet) < 3:
            raise ValueError(
                f"Could not find corresponding atoms in the other residue. For this method to work 3 corresponding atoms with the same IDs are required but only {len(ref_triplet)} were found."
            )

        return self.superimpose_to_triplet(ref_triplet, other_triplet)

    def stack(self, axis: Union[str, np.ndarray], n: int, pad: float = 0):
        """
        Stack the molecule along an axis. This will create n copies of the molecule along the axis with a padding of pad between them.
        This method is a convenience wrapper for `move` and `merge` and will not perform any kind of alignment or rotation.

        Parameters
        ----------
        axis : str or np.ndarray
            The axis to stack along. This can be either a unit vector or one of the strings "x", "y", or "z" to stack along the respective axes.
        n : int
            The number of copies to stack
        pad : float
            The padding between the copies
        """
        if isinstance(axis, (list, tuple)):
            axis = np.array(axis)
        if isinstance(axis, np.ndarray):
            axis = axis / np.linalg.norm(axis)
        elif isinstance(axis, str):
            axis = axis.lower().strip()
            sign = None
            if len(axis) == 2:
                sign, axis = axis[0], axis[1]
            if axis == "x":
                axis = structural.x_axis
            elif axis == "y":
                axis = structural.y_axis
            elif axis == "z":
                axis = structural.z_axis
            else:
                raise ValueError(f"Unknown axis input: {axis=}")
            if sign == "-":
                axis = -axis
        else:
            raise ValueError(f"Unknown axis type {type(axis)}")

        _self = self.copy()
        length = self.compute_length_along_axis(axis)
        for i in range(1, n + 1):
            incoming = _self.copy()
            incoming.move(i * axis * (length + pad))
            self.merge(incoming)
        del _self
        return self

    def align_to(self, axis: Union[str, np.ndarray]):
        """
        Align the structure (via it's primary axis, i.e. the axis perpendicular to the main plane) to some other axis.
        This will rotate the molecule so that the primary axis is aligned with the given axis. This only works for (more or less) planar molecules.

        Parameters
        ----------
        axis : str or np.ndarray
            The axis to align to. This can be either a unit vector or one of the strings "x", "y", or "z" to align to the respective axes.
        """
        if isinstance(axis, (list, tuple)):
            axis = np.array(axis)
        if isinstance(axis, np.ndarray):
            axis = axis / np.linalg.norm(axis)
        elif isinstance(axis, str):
            axis = axis.lower().strip()
            sign = None
            if len(axis) == 2:
                sign, axis = axis[0], axis[1]
            if axis == "x":
                axis = structural.x_axis
            elif axis == "y":
                axis = structural.y_axis
            elif axis == "z":
                axis = structural.z_axis
            else:
                raise ValueError(f"Unknown axis input: {axis=}")
            if sign == "-":
                axis = -axis
        else:
            raise ValueError(f"Unknown axis type {type(axis)}")

        # compute the main plane of the molecule
        plane = self.compute_normal_axis()

        # compute the rotation axis
        rot_axis = np.cross(plane, axis)
        rot_angle = np.arccos(np.dot(plane, axis))

        # rotate the molecule
        self.rotate(rot_angle, rot_axis, angle_is_degrees=False)

        return self

    def compute_normal_axis(self) -> np.ndarray:
        """
        Compute the normal axis of the molecule. This is the axis that is perpendicular to the main plane of the molecule.
        This can be computed on any molecule but will only be meaningful for (more or less) planar molecules.
        """
        return structural.plane_of_points(self.get_coords())

    def compute_principal_axis(self) -> np.ndarray:
        """
        Compute the principal axis of the molecule. This is the axis that shows the most variance in the coordinates.
        This can be computed on any molecule but will only be meaningful for (more or less) linear molecules.
        """
        return structural.principal_axis(self.get_coords())

    def compute_length_along_axis(self, axis: Union[str, np.ndarray]) -> float:
        """
        Compute the length of the molecule along a specific axis. This can be computed on any molecule but may not be meaningful in all cases (e.g. circular or branched molecules).

        Parameters
        ----------
        axis : str or np.ndarray
            The axis to compute the length along. This can be either a unit vector or one of the strings "x", "y", or "z" to align to the respective axes.
        """
        if isinstance(axis, (list, tuple)):
            axis = np.array(axis)
        if isinstance(axis, np.ndarray):
            axis = axis / np.linalg.norm(axis)
        elif isinstance(axis, str):
            axis = axis.lower().strip()
            sign = None
            if len(axis) == 2:
                sign, axis = axis[0], axis[1]
            if axis == "x":
                axis = structural.x_axis
            elif axis == "y":
                axis = structural.y_axis
            elif axis == "z":
                axis = structural.z_axis
            else:
                raise ValueError(f"Unknown axis input: {axis=}")
            if sign == "-":
                axis = -axis
        else:
            raise ValueError(f"Unknown axis type {type(axis)}")
        return structural.length_along_axis(self.get_coords(), axis)

    def bend_at_bond(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        angle: float,
        neighbor: Union[str, int, base_classes.Atom] = None,
        angle_is_degrees: bool = True,
    ):
        """
        Bend the molecule at a specific bond. This will rotate the atoms downstream of the bond in direction atom1->atom2 by the given angle.
        The axis of rotation will be the plane vector specified by the two atoms and one neighboring atom. A specific neighbor can be provided
        to ensure a specific plane is used (recommended), otherwise a random neighbor of atom1 will be used
        (preference is given to non-Hydrogens but a Hydrogen will be used if no other neighbor is found).

        Parameters
        ----------
        atom1 : Union[str, int, base_classes.Atom]
            The first atom of the bond
        atom2 : Union[str, int, base_classes.Atom]
            The second atom of the bond
        angle : float
            The angle to bend by
        neighbor : Union[str, int, base_classes.Atom], optional
            The atom to use as a neighbor for the plane vector, by default None, in which case a random neighbor of atom1 will be used.
            It is recommended to specify this to ensure a specific plane is used.
        angle_is_degrees : bool, optional
            Whether the angle is given in degrees (default) or radians
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if neighbor is not None:
            neighbor = self.get_atom(neighbor)
        else:
            H = self.get_hydrogen(atom1)
            neighbor = self.get_neighbors(atom1) - {atom2, H}
            if len(neighbor) == 0:
                if H is None:
                    raise ValueError(
                        "No neighbor found and no Hydrogen available to construct a plane vector! Specify a neighbor manually."
                    )
                neighbor = H
            else:
                neighbor = neighbor.pop()

        if angle_is_degrees:
            angle = np.radians(angle)

        descendants = list(self.get_descendants(atom1, atom2)) + [atom2]
        plane = structural.plane_of_points([a.coord for a in [neighbor, atom1, atom2]])
        coords = np.array([a.coord for a in descendants])
        coords -= atom1.coord
        new_coords = structural.rotate_coords(coords, angle, plane)
        new_coords += atom1.coord
        for adx, atom in enumerate(descendants):
            atom.coord = new_coords[adx]

        return self

    def rotate_descendants(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        angle: float,
        angle_is_degrees: bool = True,
    ):
        """
        Rotate all descendant atoms (atoms after atom2) of a bond.

        Parameters
        ----------
        atom1 : Union[str, int, base_classes.Atom]
            The first atom
        atom2 : Union[str, int, base_classes.Atom]
            The second atom (whose downstream neighbors are rotated)
        angle : float
            The angle to rotate by
        angle_is_degrees : bool
            Whether the angle is given in degrees (default) or radians
        """
        self.rotate_around_bond(
            atom1,
            atom2,
            angle,
            descendants_only=True,
            angle_is_degrees=angle_is_degrees,
        )
        return self

    def rotate_ancestors(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        angle: float,
        angle_is_degrees: bool = True,
    ):
        """
        Rotate all ancestor atoms (atoms before atom1) of a bond

        Parameters
        ----------
        atom1 : Union[str, int, base_classes.Atom]
            The first atom (whose upstream neighbors are rotated)
        atom2 : Union[str, int, base_classes.Atom]
            The second atom
        angle : float
            The angle to rotate by
        angle_is_degrees : bool
            Whether the angle is given in degrees (default) or radians
        """
        self.rotate_around_bond(
            atom2,
            atom1,
            angle,
            descendants_only=True,
            angle_is_degrees=angle_is_degrees,
        )
        return self

    def rotate_around_bond(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        angle: float,
        descendants_only: bool = False,
        angle_is_degrees: bool = True,
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
        angle_is_degrees
            Whether the angle is given in degrees (default) or radians
        
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
        if angle_is_degrees:
            angle = np.radians(angle)
        self._rotate_around_bond(atom1, atom2, angle, descendants_only)
        return self

    def get_ancestors(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
    ) -> set:
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
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
    ) -> set:
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
        atom: Union[int, str, tuple, base_classes.Atom],
        n: int = 1,
        mode: str = "upto",
        filter: callable = None,
    ) -> set:
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
        filter
            A filter function that is applied to the neighbors. If the filter returns True, the atom is included in the result.

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
        out = self._AtomGraph.get_neighbors(atom, n, mode)
        if filter:
            return {a for a in out if filter(a)}
        return out

    def get_equatorial_neighbor(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the equatorial neighbor of an atom, if the atom is in a ring structure.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The equatorial neighbor, if it exists, None otherwise
        """
        atom = self.get_atom(atom)
        return structural.get_equatorial_neighbor(self, atom)

    def get_axial_neighbor(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the axial neighbor of an atom, if the atom is in a ring structure.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The axial neighbor, if it exists, None otherwise
        """
        atom = self.get_atom(atom)
        return structural.get_axial_neighbor(self, atom)

    def get_equatorial_hydrogen(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the equatorial hydrogen neighbor of an atom, if the atom is in a ring structure.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The equatorial hydrogen, if it exists, None otherwise
        """
        atom = self.get_atom(atom)
        return structural.get_equatorial_hydrogen_neighbor(self, atom)

    def get_axial_hydrogen(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the axial hydrogen neighbor of an atom, if the atom is in a ring structure.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The axial hydrogen, if it exists, None otherwise
        """
        atom = self.get_atom(atom)
        return structural.get_axial_hydrogen_neighbor(self, atom)

    def get_left_hydrogen(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the "left-protruding" hydrogen neighbor of an atom with two hydrogens and two non-hydrogen neighbors.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The left hydrogen, if it exists, None otherwise

        Example
        -------
        In a molecule:
        ```
                   H_B
                   |
            CH3 -- C -- CH2 -- OH
                   |
                   H_A
        ```
        We want to get the left and right hydrogens of the central C atom (labeled only C).
        Using part of the logic behind R/S nomenclature for chiral centers, we prioritize the non-H neighbors
        and then rotate the molecule such that the highest order non-H neighbor points toward the user and the other
        non-H neighbor points away. The left and right hydrogens are then determined based on their orientation in this view.

        In this case, the left hydrogen is H_A and the right hydrogen is H_B.
        """
        atom = self.get_atom(atom)
        return structural.get_left_hydrogen_neighbor(self, atom)

    def get_right_hydrogen(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get the "right-protruding" hydrogen neighbor of an atom with two hydrogens and two non-hydrogen neighbors.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The right hydrogen, if it exists, None otherwise

        Example
        -------
        In a molecule:
        ```
                   H_B
                   |
            CH3 -- C -- CH2 -- OH
                   |
                   H_A
        ```
        We want to get the left and right hydrogens of the central C atom (labeled only C).
        Using part of the logic behind R/S nomenclature for chiral centers, we prioritize the non-H neighbors
        and then rotate the molecule such that the highest order non-H neighbor points toward the user and the other
        non-H neighbor points away. The left and right hydrogens are then determined based on their orientation in this view.

        In this case, the left hydrogen is H_A and the right hydrogen is H_B.
        """
        atom = self.get_atom(atom)
        return structural.get_right_hydrogen_neighbor(self, atom)

    def get_hydrogen(
        self, atom: Union[int, str, tuple, base_classes.Atom]
    ) -> base_classes.Atom:
        """
        Get any hydrogen neighbor of an atom.

        Parameters
        ----------
        atom
            The atom

        Returns
        -------
        Atom
            The hydrogen, if it exists, None otherwise
        """
        Hs = self.get_neighbors(atom, n=1, filter=lambda a: a.element == "H")
        if len(Hs) == 0:
            return None
        return Hs.pop()

    def get_hydrogens(
        self, atom: Union[int, str, tuple, base_classes.Atom] = None
    ) -> set:
        """
        Get multiple hydrogen atoms

        Parameters
        ----------
        atom
            A specific atom whose hydrogen neighbors should be returned. If None, all hydrogen atoms in the molecule are returned.

        Returns
        -------
        set
            A set of hydrogen atoms
        """
        if atom is not None:
            return self.get_neighbors(atom, n=1, filter=lambda a: a.element == "H")
        else:
            return set(self.get_atoms("H", by="element"))

    def trans(self, *bond: Union[base_classes.Atom, tuple, base_classes.Bond]):
        """
        Rotate the molecule such that the atoms in the bond are in a trans configuration.

        Parameters
        ----------
        *bond : Atom or tuple or Bond
            The bond to rotate
        """
        if len(bond) == 1 and isinstance(bond[0], (tuple, base_classes.Bond)):
            bond = bond[0]
        a = self.get_atom(bond[0])
        b = self.get_atom(bond[1])
        neighbor1 = self.get_neighbors(a, filter=lambda x: x.element != "H") - {b}
        neighbor2 = self.get_neighbors(b, filter=lambda x: x.element != "H") - {a}
        if len(neighbor1) == 0:
            neighbor1 = self.get_neighbors(a) - {b}
        if len(neighbor2) == 0:
            neighbor2 = self.get_neighbors(b) - {a}
        if len(neighbor1) == 0 or len(neighbor2) == 0:
            raise ValueError("Both atoms must have at least one neighbor!")

        neighbor1 = neighbor1.pop()
        neighbor2 = neighbor2.pop()

        angle = structural.angle_between(neighbor1.coord, a.coord, neighbor2.coord)
        if angle < 100:
            self._rotate_around_bond(a, b, np.pi, descendants_only=True)
        return self

    def cis(self, *bond: Union[base_classes.Atom, tuple, base_classes.Bond]):
        """
        Rotate the molecule such that the atoms in the bond are in a cis configuration.

        Parameters
        ----------
        *bond : Atom or tuple or Bond
            The bond to rotate
        """
        if len(bond) == 1 and isinstance(bond[0], (tuple, base_classes.Bond)):
            bond = bond[0]
        a = self.get_atom(bond[0])
        b = self.get_atom(bond[1])
        neighbor1 = self.get_neighbors(a, filter=lambda x: x.element != "H") - {b}
        neighbor2 = self.get_neighbors(b, filter=lambda x: x.element != "H") - {a}
        if len(neighbor1) == 0:
            neighbor1 = self.get_neighbors(a) - {b}
        if len(neighbor2) == 0:
            neighbor2 = self.get_neighbors(b) - {a}
        if len(neighbor1) == 0 or len(neighbor2) == 0:
            raise ValueError("Both atoms must have at least one neighbor!")

        neighbor1 = neighbor1.pop()
        neighbor2 = neighbor2.pop()
        angle = structural.angle_between(neighbor1.coord, a.coord, neighbor2.coord)
        if angle > 100:
            self._rotate_around_bond(a, b, np.pi, descendants_only=True)
        return self

    def is_cis(self, *bond: Union[base_classes.Atom, tuple, base_classes.Bond]) -> bool:
        """
        Check if the atoms in the bond are in a cis configuration.

        Parameters
        ----------
        *bond : Atom or tuple or Bond
            The bond to check

        Returns
        -------
        bool
            Whether the bond is in a cis configuration
        """
        if len(bond) == 1 and isinstance(bond[0], (tuple, base_classes.Bond)):
            bond = bond[0]
        a = self.get_atom(bond[0])
        b = self.get_atom(bond[1])
        neighbor1 = self.get_neighbors(a, filter=lambda x: x.element != "H") - {b}
        neighbor2 = self.get_neighbors(b, filter=lambda x: x.element != "H") - {a}
        if len(neighbor1) == 0:
            neighbor1 = self.get_neighbors(a) - {b}
        if len(neighbor2) == 0:
            neighbor2 = self.get_neighbors(b) - {a}
        if len(neighbor1) == 0 or len(neighbor2) == 0:
            raise ValueError("Both atoms must have at least one neighbor!")

        neighbor1 = neighbor1.pop()
        neighbor2 = neighbor2.pop()
        angle = structural.angle_between(neighbor1.coord, a.coord, neighbor2.coord)
        return angle > 100

    def is_trans(
        self, *bond: Union[base_classes.Atom, tuple, base_classes.Bond]
    ) -> bool:
        """
        Check if the atoms in the bond are in trans configuration

        Parameters
        ----------
        *bond : Atom or tuple or Bond
            The bond to check

        Returns
        -------
        bool
            Whether the bond is in a trans configuration
        """
        return not self.is_cis(*bond)

    def search_by_constraints(self, constraints: list) -> list:
        """
        Search for atoms based on a list of constraints. The constraints must be constraint functions
        from `structural.neighbors.constraints`. Each entry in the constraints list represents the constraints for one specific atom.
        Constraints apply to atom neighborhoods not the atom graph as a whole! This means that constraints are applied to the neighbors of the atoms when searching!

        Parameters
        ----------
        constraints : list
            A list of constraint functions
        
        Returns
        -------
        list
            A list of matching atoms. Each entry in this list will be a dictionary mapping the atoms (values) to the constraint function index for which they match (key).

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
        we can search for the metyhl groups by using the following constraints:
        
        >>> from buildamol.core.structural.neighbors import constraints
        >>> constraints = [
        ...     # the first atom must be a carbon and have three hydrogen neighbors
        ...     # we only search for the methyl-carbons...
        ...     constraints.multi_constraint(
        ...         constraints.has_element("C"),
        ...         constraints.has_neighbor_hist({"H": 3}),
        ...     ),
        ...     ]
        >>> mol.search_by_constraints(constraints)
        [{0: (1)C}, {0: (2)C}]
        """
        return self._AtomGraph.search_by_constraints(constraints)

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
        for model in self.get_models():

            cdx = start_chainid - 1
            rdx = start_resid
            adx = start_atomid

            for chain in model.child_list:
                chain._id = utils.auxiliary.chain_id_maker(cdx)
                cdx += 1

                for residue in chain.child_list:
                    residue.serial_number = rdx
                    rdx += 1

                    for atom in residue.child_list:
                        atom.serial_number = adx
                        adx += 1
        return self

    def index_by_chain(self):
        """
        Reindex the residues in the structure by chain. This will let each chain start with a residue 1.
        This will not reindex the atoms, only the residues.
        """
        for model in self.get_models():
            for chain in model.child_list:
                rdx = 1
                for residue in chain.child_list:
                    residue.serial_number = rdx
                    rdx += 1
        return self

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
        return mol.reindex(cdx + 1, rdx + 1, adx + 1)

    def get_chains(self):
        return self._model.get_chains()

    def get_models(self):
        return self._base_struct.get_models()

    def split_models(self) -> list:
        """
        Split the molecule into multiple molecules, each containing one of the models.
        """
        models = []
        for model in self.get_models():
            new = self.__class__.empty(id=self.id)
            new._base_struct.child_list.clear()
            new._base_struct.child_dict.clear()
            model.parent = new._base_struct
            model.id = 0
            new.add_model(model)
            models.append(new)
        return models

    def set_model(self, model: int):
        """
        Set the current working model of the molecule

        Parameters
        ----------
        model : Int
            The id of the model to set as active
        """
        if isinstance(model, int):
            if model < len(self._base_struct.child_list):
                new_model = self._base_struct.child_list[model]
            else:
                raise IndexError(
                    f"Model {model} not in molecule. Available models: {self.models}. First add the model to the molecule to set it as the active model!"
                )
        elif isinstance(model, base_classes.Model):
            if model in self._base_struct.child_list:
                new_model = model
            else:
                raise ValueError(
                    f"Model {model} not in molecule. Available models: {self.models}. First add the model to the molecule to set it as the active model!"
                )
        # we have to update the bond references to the new model
        # since each model stores copies of the atoms
        atom_mapping = {i.serial_number: i for i in new_model.get_atoms()}
        for bond in self.get_bonds():
            bond.atom1 = atom_mapping[bond.atom1.serial_number]
            bond.atom2 = atom_mapping[bond.atom2.serial_number]
        del atom_mapping
        self._model = new_model
        return self

    def get_model(self, model: int = None) -> base_classes.Model:
        """
        Get a model from the molecule.

        Parameters
        ----------
        model : Int
            The id of the model to get. If not provided the current working model is returned.

        Returns
        -------
        Model
            The model
        """
        if model is None:
            return self._model

        if isinstance(model, int):
            return self._base_struct.child_list[model]
        elif isinstance(model, base_classes.Model):
            if model in self._base_struct.child_list:
                return model
            else:
                raise ValueError(
                    f"Model {model} not in molecule. Available models: {self.get_models()}. First add the model to the molecule to set it as the active model!"
                )

    def add_model(self, model: Union[int, base_classes.Model] = None):
        """
        Add a new model to the molecule's structure

        Parameters
        ----------
        model : int or Model
            If not given, a new completely blank model is created. If an integer is given, an existing model
            is copied and added to the molecule. If a Model object is given, it is added to the molecule.
        """
        if isinstance(model, int):
            new = self.get_model(model).copy()
            new.id = len(self._base_struct.child_list)
        elif isinstance(model, base_classes.Model):
            new = model
            new.id = len(self._base_struct.child_list)
        elif isinstance(model, BaseEntity):
            new = model._model
            new.id = len(self._base_struct.child_list)
        elif model is None:
            new = base_classes.Model(len(self._base_struct.child_list))
        else:
            raise ValueError(f"Unknown model type {type(model)}")

        self._base_struct.add(new)
        return self

    def remove_model(self, model: Union[int, base_classes.Model]):
        """
        Remove a model from the molecule

        Parameters
        ----------
        model : int or Model
            The model to remove
        """
        if isinstance(model, int):
            model = self.get_model(model)
        self._base_struct.child_list.remove(model)
        self._base_struct.child_dict.pop(model.get_id())
        return self

    def get_structure(self) -> base_classes.Structure:
        return self._base_struct

    def get_residues(
        self,
        *residues: Union[int, str, tuple, base_classes.Residue],
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
            The type of parameter to search for. Can be either 'name', 'seqid' (or 'serial') or 'full_id'
            By default, this is inferred from the datatype of the residue parameter.
            If it is an integer, it is assumed to be the sequence identifying number (serial number),
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
            if isinstance(residue, base_classes.Residue):
                if residue in self._model.get_residues():
                    _residues.append(residue)
                    continue
                else:
                    _residues.append(
                        self.get_residue(residue.id[1], by="seqid", chain=chain)
                    )
                    continue

            if by is None:
                by = infer_search_param(residue)
            if by == "id":
                by = "name"
            elif by == "seqid":
                by = "serial"

            if by == "name":
                _residue = [
                    i for i in self._model.get_residues() if i.resname == residue
                ]
            elif by == "serial":
                if residue < 0:
                    residue = len(self.residues) + residue + 1
                _residue = [i for i in self._model.get_residues() if i.id[1] == residue]
            elif by == "full_id":
                _residue = [
                    i for i in self._model.get_residues() if i.full_id == residue
                ]
            else:
                raise ValueError(
                    f"Unknown search parameter '{by}', must be either 'name', 'seqid' or 'full_id' -> erroneous input: {residue=}"
                )
            if chain is not None:
                chain = self.get_chain(chain)
                _residue = [i for i in _residue if i.get_parent() == chain]

            _residues.extend(_residue)

        return _residues

    def count_bonds(self) -> int:
        """
        Count the number of bonds in the structure

        Returns
        -------
        int
            The number of bonds
        """
        return len(self._bonds)

    def count_atoms(self) -> int:
        """
        Count the number of atoms in the structure

        Returns
        -------
        int
            The number of atoms
        """
        return sum(1 for i in self._model.get_atoms())

    def count_residues(self) -> int:
        """
        Count the number of residues in the structure

        Returns
        -------
        int
            The number of residues
        """
        return sum(1 for i in self._model.get_residues())

    def count_chains(self) -> int:
        """
        Count the number of chains in the structure

        Returns
        -------
        int
            The number of chains
        """
        return sum(1 for i in self._model.get_chains())

    def count_models(self) -> int:
        """
        Count the number of models in the structure

        Returns
        -------
        int
            The number of models
        """
        return sum(1 for i in self.get_models())

    def get_atoms(
        self,
        *atoms: Union[int, str, tuple],
        by: str = None,
        keep_order: bool = False,
        residue: Union[int, base_classes.Residue] = None,
        filter: callable = None,
    ) -> list:
        """
        Get one or more atoms from the structure either based on their
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues all of them are returned in a list.
        It is a safer option to use the full_id or serial number to retrieve a specific atom.
        If no search parameters are provided, the underlying atom-generator of the structure is returned.

        Note
        ----
        This does not support mixed queries. I.e. you cannot query for an atom with id 'C1' and serial number 1 at the same time.
        Each call can only query for one type of parameter.

        Parameters
        ----------
        atoms
            The atom id, serial number, full_id tuple, or element string symbol.
            This supports multiple atoms to search for. However,
            only one type of parameter is supported per call. If left empty, the underlying generator is returned.

        by : str
            The type of parameter to search for. Can be either 'id', 'serial', 'full_id', or 'element'
            If None is given, the parameter is inferred from the datatype of the atoms argument
            'serial' in case of `int`, 'id' in case of `str`, `full_id` in case of a tuple.

        keep_order : bool
            Whether to return the atoms in the order they were queried. If False, the atoms are returned in the order they appear in the structure.

        residue: int or Residue
            A specific residue to search in. If None, the entire structure is searched.

        filter : callable
            A filter function that is applied to the atoms. If the filter returns True, the atom is included in the result.
            The filter function must take an atom as its only argument and return a boolean.

        Returns
        -------
        atom : list or generator
            The atom(s)
        """
        if len(atoms) == 0:
            if residue is not None:
                _residue = self.get_residue(residue)
                if _residue is None:
                    raise ValueError(f"Residue {residue} not found")
                atom_gen = _residue.get_atoms
            else:
                atom_gen = self._model.get_atoms
            if filter is not None:
                return (a for a in atom_gen() if filter(a))
            return atom_gen()

        elif len(atoms) == 1 and isinstance(atoms[0], (list, set, tuple)):
            atoms = atoms[0]

        if isinstance(atoms[0], base_classes.Atom):
            atoms_to_keep = [i for i in atoms if i.parent.parent.parent is self._model]
            atoms_to_get = [i.full_id for i in set(atoms) - set(atoms_to_keep)]
            if len(atoms_to_get) > 0:
                atoms_to_get = self.get_atoms(
                    atoms_to_get, by="full_id", residue=residue
                )
            _atoms = atoms_to_keep + atoms_to_get
            if filter:
                _atoms = [a for a in _atoms if filter(a)]
            if keep_order:
                return sorted(_atoms, key=lambda x: atoms.index(x))
            return _atoms

        if by is None:
            by = infer_search_param(atoms[0])

        if residue is not None:
            _residue = self.get_residue(residue)
            if residue is None:
                raise ValueError(f"Residue '{residue}' not found")
            atom_gen = _residue.get_atoms
        else:
            atom_gen = self._model.get_atoms

        # these used to be list-comprehensions. Revert to them if it causes problems
        # that the generators are used here...
        if by == "id":
            _atoms = (i for i in atom_gen() if i.id in atoms)
        elif by == "serial":
            _atoms = (i for i in atom_gen() if i.serial_number in atoms)
        elif by == "full_id":
            _atoms = (self.get_atom(i, by="full_id", residue=residue) for i in atoms)
        elif by == "element":
            atoms = [i.upper() for i in atoms]
            _atoms = (i for i in atom_gen() if i.element in atoms)
        else:
            raise ValueError(
                f"Unknown search parameter '{by}', must be either 'id', 'serial', 'full_id', or 'element' -> erroneous input: {atoms=}"
            )

        if filter:
            _atoms = (a for a in _atoms if filter(a))

        # finally make sure that the order of the atoms is the same as the input
        if keep_order:
            if by == "id":
                _atoms = sorted(_atoms, key=lambda x: atoms.index(x.id))
            elif by == "serial":
                _atoms = sorted(_atoms, key=lambda x: atoms.index(x.serial_number))
            elif by == "full_id":
                _atoms = sorted(_atoms, key=lambda x: atoms.index(x.full_id))
            elif by == "element":
                _atoms = sorted(_atoms, key=lambda x: atoms.index(x.element))

        if not isinstance(_atoms, list):
            _atoms = list(_atoms)
        return _atoms

    def get_atom(
        self,
        atom: Union[int, str, tuple],
        by: str = None,
        residue: Union[int, base_classes.Residue] = None,
    ):
        """
        Get an atom from the structure either based on its
        id, serial number or full_id. Note, if multiple atoms match the requested criteria,
        for instance there are multiple 'C1' from different residues, only the first one is returned.
        To get all atoms matching the criteria, use the `get_atoms` method.

        Parameters
        ----------
        atom
            The atom id, serial number, full_id tuple, or element symbol.
        by : str
            The type of parameter to search for. Can be either 'id', 'serial', 'full_id', or 'element'.
            Because this looks for one specific atom, this parameter can be inferred from the datatype
            of the atom parameter. If it is an integer, it is assumed to be the serial number,
            if it is a string, it is assumed to be the atom id and if it is a tuple, it is assumed
            to be the full_id.

        residue: int or Residue
            A specific residue to search in. If None, the entire structure is searched.

        Returns
        -------
        atom : base_classes.Atom
            The atom
        """
        if isinstance(atom, base_classes.Atom):
            if atom.parent.parent.parent is self._model:
                return atom
            else:
                return self.get_atom(atom.full_id, by="full_id")

        if residue is not None:
            _residue = self.get_residue(residue)
            if _residue is None:
                raise ValueError(f"Residue {residue} not found")
            residue = _residue
            atom_gen = residue.get_atoms
        else:
            atom_gen = self._model.get_atoms

        if by is None:
            by = infer_search_param(atom)

        if by == "id":
            _atom = (i for i in atom_gen() if i.id == atom)
        elif by == "serial":
            if atom < 0:
                atom = self.count_atoms() + atom + 1
            _atom = (i for i in atom_gen() if i.serial_number == atom)
        elif by == "full_id":
            if residue is None:
                try:
                    _model = next(i for i in self.get_models() if i.id == atom[1])
                    _chain = next(i for i in _model.get_chains() if i.id == atom[2])
                    _residue = next(
                        i
                        for i in _chain.get_residues()
                        if i.serial_number == atom[3][1]
                    )
                    _atom = (i for i in _residue.get_atoms() if i.id == atom[4][0])
                except StopIteration:
                    _atom = (i for i in atom_gen() if i.full_id == atom)
            else:
                if residue.serial_number != atom[3][1]:
                    return None
                _atom = (i for i in atom_gen() if i.id == atom[4][0])
        elif by == "element":
            _atom = (i for i in atom_gen() if i.element == atom.upper())
        else:
            raise ValueError(
                f"Unknown search parameter '{by}', must be either 'id', 'serial', 'full_id', or 'element' -> erroneous input: {atom=}"
            )

        return next(_atom, None)

    def get_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
        add_if_not_present: bool = True,
    ) -> base_classes.Bond:
        """
        Get/make a bond between two atoms.

        Parameters
        ----------
        atom1: str or int or tuple or Atom
            The first atom
        atom2: str or int or tuple or Atom
            The second atom
        add_if_not_present : bool
            Whether to add the bond if it is not present

        Returns
        -------
        bond : Bond
            The bond object. If the bond is not present and add_if_not_present is False, None is returned.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        return self._get_bond(atom1, atom2, add_if_not_present)

    def _get_bond(self, atom1, atom2, add_if_not_present=True):
        """
        Get/make a bond between two atoms.
        This is the core method that is used to get bonds between atoms.
        It expects the atoms to be Atom objects which are present in the molecule.
        """
        has_edge = self._AtomGraph.has_edge(atom1, atom2)
        if add_if_not_present and not has_edge:
            self._set_bond(atom1, atom2)
        elif not has_edge:
            return None
        return self._AtomGraph.edges[atom1, atom2]["bond_obj"]

    def get_bonds(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom, base_classes.Residue] = None,
        atom2: Union[int, str, tuple, base_classes.Atom] = None,
        residue_internal: bool = True,
        either_way: bool = True,
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
        residue_internal : bool
            If True, only bonds where both atoms are in the given residue (if `atom1` is a residue) are returned.
            If False, all bonds where either atom is in the given residue are returned.
        either_way : bool
            If True, the order of the atoms does not matter, if False, the order of the atoms
            does matter. By setting this to false, it is possible to also search for bonds that have
            a specific atom in position 1 or 2 depending on which argument was set, while leaving the other input as none.

        Returns
        -------
        bond : list or generator
            The bond(s). If no input is given, all bonds are returned as a generator.
        """
        if atom1 is None and atom2 is None:
            return iter(self._bonds)

        if isinstance(atom1, base_classes.Residue):
            if residue_internal:
                return [
                    i
                    for i in self.get_bonds()
                    if i[0].parent is atom1 and i[1].parent is atom1
                ]
            else:
                return [
                    i
                    for i in self.get_bonds()
                    if i[0].parent is atom1 or i[1].parent is atom1
                ]

        if atom1:
            atom1 = self.get_atoms(atom1)
        if atom2:
            atom2 = self.get_atoms(atom2)

        return self._get_bonds(atom1, atom2, either_way)

    def set_bond_order(self, atom1, atom2, order: int, adjust_hydrogens: bool = False):
        """
        Set the order of a bond between two atoms

        Parameters
        ----------
        atom1
            The first atom
        atom2
            The second atom
        order : int
            The order of the bond
        adjust_hydrogens : bool
            Whether to adjust the number of hydrogens on the atoms based on the bond order
        """
        bond = self.get_bond(atom1, atom2)
        self._AtomGraph.edges[atom1, atom2]["bond_order"] = order
        bond.order = order
        if adjust_hydrogens:
            H = structural.Hydrogenator()
            self.remove_atoms(self.get_hydrogens(atom1) | self.get_hydrogens(atom2))
            H.add_hydrogens(atom1, self)
            H.add_hydrogens(atom2, self)
        return self

    def get_residue(
        self,
        residue: Union[int, str, tuple, base_classes.Residue],
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
            The type of parameter to search for. Can be either 'name', 'serial' (or 'seqid') or 'full_id'
            By default, this is inferred from the datatype of the residue parameter.
            If it is an integer, it is assumed to be the sequence identifying number (serial number),
            if it is a string, it is assumed to be the residue name and if it is a tuple, it is assumed
            to be the full_id.
        chain : str
            Further restrict to a residue from a specific chain.

        Returns
        -------
        residue : base_classes.Residue
            The residue
        """
        if isinstance(residue, base_classes.Residue):
            if residue.parent.parent is self._model:
                return residue
            else:
                return self.get_residue(residue.id[1], by="serial", chain=chain)

        if by is None:
            by = infer_search_param(residue)
        if by == "id":
            by = "name"
        elif by == "seqid":
            by = "serial"

        if by == "name":
            _residue = (i for i in self._model.get_residues() if i.resname == residue)
        elif by == "serial":
            if residue < 0:
                residue = len(self.residues) + residue + 1
            _residue = (i for i in self._model.get_residues() if i.id[1] == residue)
        elif by == "full_id":
            _model = next(i for i in self.get_models() if i.id == residue[1])
            _chain = next(i for i in _model.get_chains() if i.id == residue[2])
            _residue = (
                i for i in _chain.get_residues() if i.serial_number == residue[3][1]
            )
        else:
            raise ValueError(
                f"Unknown search parameter, must be either 'name', 'seqid'/'serial', or 'full_id' -> erroneous input: {residue=}"
            )
        if chain is not None:
            chain = self.get_chain(chain)
            _residue = (i for i in _residue if i.parent == chain)
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
        chain : base_classes.Chain
            The chain
        """
        if isinstance(chain, base_classes.Chain):
            if chain in self.chains:
                return chain
            else:
                return self.get_chain(chain.id)
        return next((i for i in self.chains if i.id == chain), None)

    def add_chains(
        self,
        *chains: base_classes.Chain,
        adjust_seqid: bool = True,
        _copy: bool = False,
    ):
        """
        Add chains to the structure

        Parameters
        ----------
        chains : base_classes.Chain
            The chains to add

        adjust_seqid : bool
            If True, the seqid of the chains is adjusted to
            match the current number of chains in the structure
            (i.e. a new chain can be given seqid A, and it will be adjusted
            to the correct value of C if there are already two other chains in the molecule).

        _copy : bool
            If True, the chains are copied before adding them to the molecule.
            This is useful if you want to add the same chain to multiple molecules, while leaving
            them and their original parent structures intakt.
        """
        for chain in chains:
            if _copy:
                chain = chain.copy()
            if adjust_seqid:
                chain._id = utils.auxiliary.chain_id_maker(len(self.chains))
            self._model.add(chain)
        return self

    def add_residues(
        self,
        *residues: base_classes.Residue,
        adjust_seqid: bool = True,
        _copy: bool = False,
    ):
        """
        Add residues to the structure

        Parameters
        ----------
        residues : base_classes.Residue
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
                # buildamol class derivatives, the below line
                # is for the original biopython classes
                # residue.id = (residue.id[0], rdx, *residue.id[2:])

            self._chain.add(residue)

            for atom in residue.child_list:
                if adjust_seqid:
                    adx += 1
                    atom.set_serial_number(adx)
                self._AtomGraph.add_node(atom)
        return self

    def remove_residues(self, *residues: Union[int, base_classes.Residue]) -> list:
        """
        Remove residues from the structure

        Parameters
        ----------
        residues : int or base_classes.Residue
            The residues to remove, either the object itself or its seqid

        Returns
        -------
        list
            The removed residues
        """
        _residues = []
        for residue in residues:
            if isinstance(residue, int):
                residue = self.get_residue(residue)

            for atom in residue.child_list:
                self._purge_bonds(atom)
                self._AtomGraph.remove_node(atom)

            # keep the memory of the parent in the residue that is removed...
            chain = residue.parent
            chain.detach_child(residue.get_id())
            residue.set_parent(chain)

            _residues.append(residue)
        return _residues

    def rename_chain(self, chain: Union[str, base_classes.Chain], name: str):
        """
        Rename a chain

        Parameters
        ----------
        chain : str or Chain
            The chain to rename, either the object itself or its id
        name : str
            The new name
        """
        chain = self.get_chain(chain)
        chain._id = name
        return self

    def rename_residue(self, residue: Union[int, base_classes.Residue], name: str):
        """
        Rename a residue

        Parameters
        ----------
        residue : int or Residue
            The residue to rename, either the object itself or its seqid
        name : str
            The new name
        """
        residue = self.get_residue(residue)
        residue.resname = name
        return self

    def rename_atom(
        self,
        atom: Union[int, base_classes.Atom],
        name: str,
        residue: Union[int, base_classes.Residue] = None,
    ):
        """
        Rename an atom

        Parameters
        ----------
        atom : int or base_classes.Atom
            The atom to rename, either the object itself or its serial number
        name : str
            The new name (id)
        residue : int or base_classes.Residue
            The residue to which the atom belongs, either the object itself or its seqid.
            Useful when giving a possibly redundant id as identifier in multi-residue molecules.
        """
        atom = self.get_atom(atom, residue=residue)
        atom.id = name
        atom.name = name
        return self

    def rename_residues(self, old_name: str, new_name: str):
        """
        Rename multiple residues to the same name

        Parameters
        ----------
        old_name : str
            The name of the residues to rename
        new_name : str

        """
        for residue in self.get_residues(old_name):
            residue.resname = new_name
        return self

    def rename_atoms(self, old_name: str, new_name: str, residue_name: str = None):
        """
        Rename multiple atoms to the same name

        Parameters
        ----------
        old_name : str
            The name of the atoms to rename
        new_name : str
            The new name
        residue_name : str
            The name of the residue of the atoms to rename (if only atoms from a specific type of residue should be renamed).
        """
        if residue_name:
            filter = lambda x: x.parent.resname == residue_name
        else:
            filter = None
        for atom in self.get_atoms(old_name, by="id", filter=filter):
            atom.name = new_name
            atom.id = new_name
        return self

    def change_element(
        self,
        atom: Union[int, base_classes.Atom],
        element: str,
        adjust_bond_length: bool = True,
    ):
        """
        Change the element of an atom. This will automatically add or remove hydrogens
        if the new element has a different valency.

        Parameters
        ----------
        atom : int or base_classes.Atom
            The atom to rename, either the object itself or its serial number
        element : str
            The new element
        adjust_bond_length : bool
            If True, adjust the bond length to match the new element. This may slow down the process if the atom is central in a very large molecule.
        """
        atom = self.get_atom(atom)

        current_element = atom.element

        structural.change_element(atom, element, self)

        if adjust_bond_length:
            neighbor = self.get_neighbors(atom).pop()
            dist = structural.single_bond_lengths.get(neighbor.element.title(), {}).get(
                atom.element.title(), None
            )
            if dist:
                self.adjust_bond_length(
                    neighbor, atom, dist, move_descendants=current_element != "H"
                )

        return self

    def set_charge(
        self,
        atom: Union[str, int, tuple, base_classes.Atom],
        charge: int,
        adjust_protonation: bool = True,
    ):
        """
        Set the charge of an atom. This will automatically adjust the number of protons on the atom
        if the charge is changed.

        Parameters
        ----------
        atom : str or int or tuple or Atom
            The atom whose charge should be changed
        charge : int
            The new charge. This is NOT the charge difference to apply but the final charge of the atom.
        adjust_protonation : bool
            If True, adjust the number of protons on the atom to match the charge.
        """
        atom = self.get_atom(atom)
        if adjust_protonation:
            structural.adjust_protonation(self, atom, charge)
        else:
            atom.charge = charge
        return self

    def add_atoms(self, *atoms: base_classes.Atom, residue=None, _copy: bool = False):
        """
        Add atoms to the structure. This will automatically adjust the atom's serial number to
        fit into the structure.

        Parameters
        ----------
        atoms : base_classes.Atom
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
        if len(atoms) == 1 and isinstance(atoms[0], (list, tuple)):
            atoms = iter(atoms[0])

        if residue is not None:
            target = self.get_residue(residue)
        else:
            target = self._chain.child_list[-1]

        _max_serial = sum(1 for i in self._model.get_atoms())
        for atom in atoms:
            if _copy:
                atom = atom.copy()
            _max_serial += 1
            atom.set_serial_number(_max_serial)
            target.add(atom)
            self._AtomGraph.add_node(atom)
        return self

    def link_atoms(self, *atoms: base_classes.Atom, residue=None):
        """
        Softlink atoms to the structure. This will add the atoms to the index of the maintained structure
        but it will not adjust the atoms' own parent references. This is useful if you want to have atoms
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        atoms : base_classes.Atom
            The atoms to link
        """
        if len(atoms) == 1 and isinstance(atoms[0], (list, tuple)):
            atoms = iter(atoms[0])
        if residue is not None:
            target = self.get_residue(residue)
            if target is None:
                raise ValueError(f"Residue '{residue}' not found")
        else:
            target = self._chain.child_list[-1]

        for atom in atoms:
            target.link(atom)
            self._AtomGraph.add_node(atom)
        return self

    def unlink_atoms(self, *atoms: base_classes.Atom):
        """
        Unlink atoms from the structure. This will remove the atoms from the index of the maintained structure
        but it will not adjust the atoms' own parent references. This is useful if you want to have atoms
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        atoms : base_classes.Atom
            The atoms to unlink
        """
        if len(atoms) == 1 and isinstance(atoms[0], (list, tuple)):
            atoms = iter(atoms[0])

        for atom in atoms:
            for residue in self.get_residues():
                if atom in residue.child_list:
                    self._purge_bonds(atom)
                    self._AtomGraph.remove_node(atom)
                    residue.unlink(atom)
                    break

        return self

    def link_residues(self, *residues: base_classes.Residue, chain=None):
        """
        Softlink residues to the structure. This will add the residues to the index of the maintained structure
        but it will not adjust the residues' own parent references. This is useful if you want to have residues
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        residues : base_classes.Residue
            The residues to link
        chain : str or base_classes.Chain
            The chain to which the residues should be linked.
            If None, the residues are linked to the current working chain.
        """
        if len(residues) == 1 and isinstance(residues[0], (list, tuple)):
            residues = iter(residues[0])

        if chain is not None:
            chain = self.get_chain(chain)
        else:
            chain = self._chain
        for residue in residues:
            chain.link(residue)
            for atom in residue.child_list:
                self._AtomGraph.add_node(atom)
        return self

    def unlink_residues(self, *residues: base_classes.Residue):
        """
        Unlink residues from the structure. This will remove the residues from the index of the maintained structure
        but it will not adjust the residues' own parent references. This is useful if you want to have residues
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        residues : base_classes.Residue
            The residues to unlink
        """
        if len(residues) == 1 and isinstance(residues[0], (list, tuple)):
            residues = iter(residues[0])

        for residue in residues:
            for chain in self.get_chains():
                if residue in chain.child_list:
                    chain.unlink(residue)
                    for atom in residue.child_list:
                        self._purge_bonds(atom)
                        self._AtomGraph.remove_node(atom)
                    break
        return self

    def link_chains(self, *chains: base_classes.Chain):
        """
        Softlink chains to the structure. This will add the chains to the index of the maintained structure
        but it will not adjust the chains' own parent references. This is useful if you want to have chains
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        chains : base_classes.Chain
            The chains to link
        """
        if len(chains) == 1 and isinstance(chains[0], (list, tuple)):
            chains = iter(chains[0])

        for chain in chains:
            self._model.link(chain)
            for atom in chain.get_atoms():
                self._AtomGraph.add_node(atom)
        return self

    def unlink_chains(self, *chains: base_classes.Chain):
        """
        Unlink chains from the structure. This will remove the chains from the index of the maintained structure
        but it will not adjust the chains' own parent references. This is useful if you want to have chains
        be accessible from multiple Molecule objects.

        Parameters
        ----------
        chains : base_classes.Chain
            The chains to unlink
        """
        if len(chains) == 1 and isinstance(chains[0], (list, tuple)):
            chains = iter(chains[0])

        for chain in chains:
            for model in self.get_models():
                if chain in model.child_list:
                    model.unlink(chain)
                    for atom in chain.get_atoms():
                        self._purge_bonds(atom)
                        self._AtomGraph.remove_node(atom)
                    break
        return self

    def remove_atoms(self, *atoms: Union[int, str, tuple, base_classes.Atom]) -> list:
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
        if len(atoms) == 1 and isinstance(atoms[0], (list, set, tuple)):
            atoms = iter(atoms[0])

        _atoms = []
        for atom in atoms:
            _atom = self.get_atoms(atom)
            if len(_atom) == 0:
                continue
            for atom in _atom:
                self._purge_bonds(atom)
                self._AtomGraph.remove_node(atom)
                p = atom.get_parent()
                if p:
                    p.detach_child(atom.get_id())
                    atom.set_parent(
                        p
                    )  # this is necessary to avoid a bug where atoms remain longer in the memory atoms list than they should
                _atoms.append(atom)

        # (same as in the _remove_atoms) I don't see why
        # this is necessary, but it is in the original code
        # reindex the atoms
        adx = 0
        for atom in self._model.get_atoms():
            adx += 1
            atom.serial_number = adx

        return _atoms

    def add_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
        order: int = 1,
    ):
        """
        Add a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        order : int
            The order of the bond, i.e. 1 for single, 2 for double, 3 for triple, etc.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._add_bond(atom1, atom2, order)
        return self

    def set_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
        order: int = 1,
    ):
        """
        Specify a bond between two atoms. The difference between this method and `add_bond` is that
        the latter can be used to incrementally add bond orders (i.e. make a double bond out of a single bond
        by calling the method twice). This method will always set the bond order to the provided value.

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        order : int
            The order of the bond, i.e. 1 for single, 2 for double, 3 for triple, etc.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._set_bond(atom1, atom2, order)

    def add_bonds(self, *bonds):
        """
        Add multiple bonds at once.

        Parameters
        ----------
        bonds
            The bonds to add, each bond is a tuple of two atoms.
            Each atom may be specified directly (BuildAMol object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        if len(bonds) == 1 and isinstance(bonds[0], (list, set, tuple)):
            bonds = iter(bonds[0])

        for bond in bonds:
            if isinstance(bond, base_classes.Bond):
                bond = bond.to_tuple()
            self.add_bond(*bond)
        return self

    def set_bonds(self, *bonds):
        """
        Specify multiple bonds at once. The difference between this method and `add_bonds` is that
        the latter can be used to incrementally add bond orders (i.e. make a double bond out of a single bond
        by calling the method twice or certain bonds are specified multiple times in the arguments).
        This method will always set the bond order to the provided value.

        Parameters
        ----------
        bonds
            The bonds to add, each bond is a tuple of two atoms.
            Each atom may be specified directly (BuildAMol object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        if len(bonds) == 1 and isinstance(bonds[0], (list, set, tuple)):
            bonds = iter(bonds[0])

        for bond in bonds:
            if isinstance(bond, base_classes.Bond):
                bond = bond.to_tuple()
            self.set_bond(*bond)
        return self

    def _add_bonds(self, *bonds):
        """
        Add multiple bonds at once. This requires that the tuple objects are indeed Atoms in the structure!
        """
        if len(bonds) == 1 and isinstance(bonds[0], (list, set, tuple)):
            bonds = iter(bonds[0])
        for bond in bonds:
            if isinstance(bond, base_classes.Bond):
                bond = bond.to_tuple()
            self._add_bond(*bond)
        return self

    def _set_bonds(self, *bonds):
        """
        Specify multiple bonds at once. This requires that the tuple objects are indeed Atoms in the structure!
        The difference between this method and `add_bonds` is that
        the latter can be used to incrementally add bond orders (i.e. make a double bond out of a single bond
        by calling the method twice or certain bonds are specified multiple times in the arguments).
        This method will always set the bond order to the provided value.
        """
        if len(bonds) == 1 and isinstance(bonds[0], (list, set, tuple)):
            bonds = iter(bonds[0])
        for bond in bonds:
            if isinstance(bond, base_classes.Bond):
                bond = bond.to_tuple()
            self._set_bond(*bond)
        return self

    def remove_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
    ):
        # either_way: bool = True,
        """
        Remove a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.
        """
        # either_way : bool
        #     If True, the bond will be removed in both directions, i.e. if the bond is (1, 2)
        #     it will be removed if either (1, 2) or (2, 1) is provided.
        # if either_way:
        #     self.remove_bond(atom1, atom2, either_way=False)
        #     self.remove_bond(atom2, atom1, either_way=False)
        #     return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._remove_bond(atom1, atom2)
        return self

    def purge_bonds(self, atom: Union[int, str, base_classes.Atom] = None):
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
        return self

    def get_degree(self, atom: Union[int, str, base_classes.Atom]):
        """
        Get the degree of an atom in the structure

        Parameters
        ----------
        atom
            The atom to get the degree of, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        Returns
        -------
        int
            The degree of the atom's connectivity as the sum of the bond orders that connect it to its neighbors
        """
        atom = self.get_atom(atom)
        return sum(b.order for b in self.bonds if atom in b)

    def adjust_bond_length(
        self, atom1, atom2, length: float, move_descendants: bool = False
    ):
        """
        Adjust the bond length between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (Atom object)
            or by providing the serial number, the full_id or the id of the atoms.
        length : float
            The new bond length
        move_descendants : bool
            If True, this method will infer all descendant atoms and move them accordingly to preserve the overall geometry of the
            molecule. It will make things slower, however!
        """
        atom1, atom2 = self.get_atom(atom1), self.get_atom(atom2)

        diff = length - (atom1 - atom2)
        structural.adjust_bond_length((atom1, atom2), length)
        if move_descendants:
            vec = structural.norm_vector(atom1, atom2)
            descendants = self.get_descendants(atom1, atom2)
            for d in descendants:
                d.coord += vec * diff

        return self

    def lock_all(self):  # , both_ways: bool = True):
        """
        Lock all bonds in the structure so they cannot be rotated around
        """

        # Parameters
        # ----------
        # both_ways : bool
        #     If True, the bond is locked in both directions
        #     i.e. atom1 --- atom2 direction will be unavailable for
        #     rotation as well as atom2 --- atom1 direction as well.
        self._AtomGraph.lock_all()
        # if both_ways:
        #     self.locked_bonds.update(b[::-1] for b in self.bonds)
        return self

    def unlock_all(self):
        """
        Unlock all bonds in the structure
        """
        self._AtomGraph.unlock_all()
        return self

    def lock_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
    ):
        # both_ways: bool = False,
        """
        Lock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        """
        # both_ways : bool
        #     If True, the bond is locked in both directions. By default the bond is only locked in the specified direction.
        # if both_ways:
        #     self.lock_bond(atom2, atom1, both_ways=False)
        #     self.lock_bond(atom2, atom1, both_ways=False)
        #     return

        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        self._AtomGraph.lock_edge(atom1, atom2)
        return self

    def unlock_bond(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
    ):
        # both_ways: bool = False,
        """
        Unlock a bond between two atoms

        Parameters
        ----------
        atom1, atom2
            The atoms to bond, which can either be directly provided (biopython object)
            or by providing the serial number, the full_id or the id of the atoms.

        """
        # both_ways : bool
        #     If True, the bond is unlocked in both directions. By default the bond is only unlocked in the specified direction.
        # if both_ways:
        #     self.unlock_bond(atom1, atom2, both_ways=False)
        #     self.unlock_bond(atom2, atom1, both_ways=False)
        #     return
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        if self._AtomGraph.has_edge(atom1, atom2):
            self._AtomGraph.unlock_edge(atom1, atom2)
        # if both_ways and (atom2, atom1) in self._AtomGraph.edges:
        #     self._AtomGraph.unlock_edge(atom2, atom1)
        return self

    def is_locked(
        self,
        atom1: Union[int, str, tuple, base_classes.Atom],
        atom2: Union[int, str, tuple, base_classes.Atom],
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
        self,
        max_bond_length: float = None,
        restrict_residues: bool = True,
        infer_bond_orders: bool = False,
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
        infer_bond_orders : bool
            Whether to infer the bond orders (double and tripple bonds) based on registered functional groups.
            This will slow the inference down, however.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.infer_bonds(
            self._base_struct, max_bond_length, restrict_residues
        )
        self._set_bonds(*bonds)
        if infer_bond_orders:
            structural.infer_bond_orders(self)

        return bonds

    def infer_bonds_for(self, *residues, max_bond_length: float = None):
        """
        Infer bonds between atoms in the structure for a specific set of residues

        Parameters
        ----------
        residues
            The residues to consider
        max_bond_length : float
            The maximum distance between atoms to consider them bonded.
            If None, the default value is 1.6 Angstroms.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = []
        for res in self.get_residues(*residues):
            bonds.extend(
                structural.infer_bonds(res, max_bond_length, restrict_residues=False)
            )
        self._set_bonds(*bonds)
        return bonds

    def get_residue_connections(
        self,
        residue_a=None,
        residue_b=None,
        triplet: bool = True,
        rotatable_only: bool = False,
    ):
        """
        Get bonds between atoms that connect different residues in the structure
        This method is different from `infer_residue_connections` in that it works
        with the already present bonds in the molecule instead of computing new ones.

        Parameters
        ----------
        residue_a, residue_b : Union[int, str, tuple, base_classes.Residue]
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
        rotatable_only : bool
            Whether to only return bonds that are rotatable. This is useful if you want to
            use the returned bonds for optimization.

        Returns
        -------
        list
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = (i for i in self._bonds if i[0].get_parent() != i[1].get_parent())

        if residue_a is not None and residue_b is None:
            residue_a = self.get_residues(residue_a)
            bonds = (
                i
                for i in bonds
                if i[0].get_parent() in residue_a or i[1].get_parent() in residue_a
            )
        elif residue_b is not None and residue_a is None:
            residue_b = self.get_residues(residue_b)
            bonds = (
                i
                for i in bonds
                if i[0].get_parent() in residue_b or i[1].get_parent() in residue_b
            )
        elif residue_a is not None and residue_b is not None:
            residue_a, residue_b = self.get_residues(residue_a), self.get_residues(
                residue_b
            )
            bonds = (
                i
                for i in bonds
                if (i[0].get_parent() in residue_a and i[1].get_parent() in residue_b)
                or (i[1].get_parent() in residue_a and i[0].get_parent() in residue_b)
            )
        if triplet:
            bonds = self._make_bond_triplets(bonds)

        if rotatable_only:
            bonds = [
                bond
                for bond in bonds
                if bond not in self._locked_bonds
                and not self._AtomGraph.in_same_cycle(*bond)
            ]
            return bonds
        return [b for b in bonds]

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
        for bond in bonds:
            atom1, atom2 = bond
            neighs = self.get_neighbors(atom1)
            neighs.discard(atom2)
            neighs -= set(i for i in neighs if i.element == "H")
            if len(neighs) == 1:
                neigh = neighs.pop()
                if neigh.parent is atom1.parent:
                    _new.add((atom1, neigh))
                continue

            neighs = self.get_neighbors(atom2)
            neighs.discard(atom1)
            neighs -= set(i for i in neighs if i.element == "H")
            if len(neighs) == 1:
                neigh = neighs.pop()
                if neigh.parent is atom2.parent:
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
            A list of bonds that link atoms from different residues.

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
        _bonds = [base_classes.Bond(*b) for b in bonds]
        _bonds = [b for b in _bonds if b not in self._bonds]
        self._bonds.extend(_bonds)
        self._AtomGraph.add_edges_from(bonds)
        for b in _bonds:
            self._AtomGraph.edges[b[0], b[1]]["bond_order"] = 1
            self._AtomGraph.edges[b[0], b[1]]["bond_obj"] = b
        return bonds

    def apply_standard_bonds(self, _compounds=None) -> list:
        """
        Use reference compounds to infer bonds in the structure. This will be exclusively based on the
        residue and atom ids and not on the actual distances between atoms.

        Parameters
        ----------
        _compounds
            The compounds to use for the standard bonds. If None, the default compounds are used.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        bonds = structural.apply_reference_bonds(self._base_struct, _compounds)
        self._add_bonds(*bonds)
        return bonds

    def apply_standard_bonds_for(self, *residues, _compounds=None) -> list:
        """
        Use reference compounds to infer bonds in the structure for specific residues. This will be exclusively based on the
        residue and atom ids and not on the actual distances between atoms.

        Parameters
        ----------
        residues
            The residues to consider
        _compounds
            The compounds to use for the standard bonds. If None, the default compounds are used.

        Returns
        -------
        list
            A list of tuples of atom pairs that are bonded
        """
        if len(residues) == 1 and isinstance(residues[0], (set, list, tuple)):
            residues = residues[0]
        bonds = []
        for res in self.get_residues(*residues):
            bonds.extend(structural.apply_reference_bonds(res, _compounds))
        self._set_bonds(*bonds)
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
        return self

    def relabel_hydrogens(self):
        """
        Relabel hydrogen atoms in the structure to match the standard labelling according
        to the CHARMM force field. This is useful if you want to use some pre-generated
        PDB file that may have used a different labelling scheme for atoms.
        """
        structural.relabel_hydrogens(self)
        return self

    def add_hydrogens(self, *atoms: Union[int, str, base_classes.Atom]):
        """
        Infer missing hydrogens in the structure.

        Parameters
        ----------
        atoms
            The atoms to infer hydrogens for. If None, all atoms are considered.
        """
        if len(atoms) == 1 and isinstance(atoms[0], (list, tuple, set)):
            atoms = atoms[0]
        H = structural.infer.Hydrogenator()
        if atoms:
            atoms = self.get_atoms(*atoms)
            for a in atoms:
                H.add_hydrogens(a, self)
        else:
            H.infer_hydrogens(self, bond_length=1.05)
        return self

    def remove_hydrogens(self, *atoms: Union[int, str, base_classes.Atom]):
        """
        Remove all hydrogens in the structure.

        Parameters
        ----------
        atoms
            The atoms to remove hydrogens from. If None, all atoms are considered.
        """
        if len(atoms) == 1 and isinstance(atoms[0], (list, tuple, set)):
            atoms = atoms[0]
        if atoms:
            atoms = self.get_atoms(*atoms)
            hydrogens = set()
            for atom in atoms:
                hydrogens.update(self.get_hydrogens(atom))
            self._remove_atoms(*hydrogens)
        else:
            self._remove_atoms(*self.get_atoms("H", by="element"))
        return self

    def get_quartets(self):
        """
        A generator for all atom quartets in the structure
        """
        yield from structural.neighbors.generate_quartets(self.bonds)

    def quartet(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        atom3: Union[str, int, base_classes.Atom],
        atom4: Union[str, int, base_classes.Atom],
    ):
        """
        Make an atom quartet from four atoms.

        Parameters
        ----------
        atom1, atom2, atom3, atom4
            The four atoms that make up the quartet.
        """
        atom1 = self.get_atom(atom1)
        atom2 = self.get_atom(atom2)
        atom3 = self.get_atom(atom3)
        atom4 = self.get_atom(atom4)
        atoms = [atom1, atom2, atom3, atom4]
        for quartet in self.get_quartets():
            if quartet == atoms:
                return quartet
        raise ValueError("The given atoms do not form a quartet.")

    def compute_angle(
        self,
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        atom3: Union[str, int, base_classes.Atom],
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
        atom1: Union[str, int, base_classes.Atom],
        atom2: Union[str, int, base_classes.Atom],
        atom3: Union[str, int, base_classes.Atom],
        atom4: Union[str, int, base_classes.Atom],
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

    def get_atom_graph(self, _copy: bool = True) -> graphs.AtomGraph:
        """
        Get an AtomGraph for the Molecule

        Parameters
        ----------
        _copy : bool
            If True, not the "original" AtomGraph object that the Molecule relies on is returned but a new one.
            However, the molecule will still be linked to the new graph. This is useful if you want to make changes
            to the graph itself (not including changes to the graph nodes, i.e. the atoms itself, such as rotations).

        Returns
        -------
        AtomGraph
            The generated graph
        """
        if not _copy:
            return self._AtomGraph
        atom_graph = graphs.AtomGraph(self.id, [])
        atom_graph.add_nodes_from(self._AtomGraph.nodes)
        atom_graph.migrate_bonds(self._AtomGraph)
        atom_graph._locked_edges.update(self._AtomGraph._locked_edges)
        atom_graph._molecule = self
        return atom_graph

    def update_atom_graph(self):
        """
        Generate a new up-to-date `AtomGraph` after any manual changes were done to the Molecule's underlying biopython structure.
        """
        self._AtomGraph.clear()
        self._AtomGraph.add_nodes_from(self.get_atoms())
        self._AtomGraph.add_edges_from(self._bonds)
        return self

    def make_residue_graph(
        self, detailed: bool = False, locked: bool = True
    ) -> graphs.ResidueGraph:
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
        graph._molecule = self
        return graph

    get_residue_graph = make_residue_graph
    make_atom_graph = get_atom_graph

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

    def to_molfile(self, filename: str):
        """
        Write the molecule to a Molfile

        Parameters
        ----------
        filename : str
            Path to the Mol file
        """
        utils.sdmol.write_mol(self, filename)

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
            mol=self,
            filename=filename,
            type=type,
            names=names,
            identifiers=identifiers,
            one_letter_code=one_letter_code,
            three_letter_code=three_letter_code,
        )

    def to_xml(self, filename: str, atom_attributes: list = None):
        """
        Write the molecule to an XML file

        Parameters
        ----------
        filename : str
            Path to the XML file
        atom_attributes : list
            A list of attributes to include in the XML file. Always included are:
                - serial_number
                - id
                - element
        """
        xml = utils.xml.encode_molecule(self, atom_attributes)
        utils.xml.write_xml(filename, xml)

    def to_openmm(self):
        """
        Convert the molecule to an OpenMM Topology

        Returns
        -------
        openmm.app.PDBFile
        """
        # since we are lazy loading the modules it is possible that they raise an exception
        # the first time they are accessed, so we try to catch that here
        try:
            return utils.convert.OpenMMBioPythonConverter().buildamol_to_openmm(self)
        except:
            try:
                return utils.convert.OpenMMBioPythonConverter().buildamol_to_openmm(
                    self
                )
            except Exception as e:
                raise ValueError(
                    "Could not convert the molecule to an OpenMM Topology."
                ) from e

    def to_pybel(self):
        """
        Convert the molecule to a Pybel molecule

        Returns
        -------
        pybel.Molecule
            The Pybel molecule
        """
        conv = utils.convert.PybelBioPythonConverter()
        mol = conv.buildamol_to_pybel(self)
        if self.id is not None:
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
        mol = conv._pdbio_to_rdkit()
        if self.id is not None:
            mol.SetProp("_Name", self.id)
        return mol

    def to_stk(self):
        """
        Convert the molecule to a STK molecule

        Returns
        -------
        stk.BuildingBlock
            The STK molecule
        """
        # we use this setup because of the lazy loading of stk
        # which may raise a ValueError the first time the module is accessed
        try:
            conv = utils.convert.STKBuildAMolConverter()
            return conv.buildamol_to_stk(self)
        except:
            try:
                conv = utils.convert.STKBuildAMolConverter()
                return conv.buildamol_to_stk(self)
            except Exception as e:
                raise ValueError(
                    "Could not convert the molecule to a STK molecule."
                ) from e

    def to_biopython(self):
        """
        Convert the molecule to a Biopython structure

        Returns
        -------
        Bio.PDB.Structure.Structure
            The Biopython structure
        """
        return self._base_struct.to_biopython()

    def to_numpy(self, export_bonds: bool = True):
        """
        Convert the molecule to numpy arrays

        Parameters
        ----------
        export_bonds : bool
            If True, the bonds are also exported.
            If False, the bond array will remain empty.

        Returns
        -------
        tuple
            The atomic numbers and atomic coordinates in one array
            and the bonds with atom serial numbers and bond order in a second array
        """
        return utils.convert.mol_to_numpy_array(self, export_bonds)

    def get_coords(self, *atom_selector, **atom_selectors) -> np.ndarray:
        """
        Get the coordinates of the atoms in the molecule

        Parameters
        ----------
        atom_selectors
            Arguments or keyword arguments to pass to get_atoms(). If None, all atoms are selected.

        Returns
        -------
        np.ndarray
            The coordinates
        """
        return np.array(
            [atom.coord for atom in self.get_atoms(*atom_selector, **atom_selectors)]
        )

    def get_bond_array(self) -> np.ndarray:
        """
        Get the bonds of the atoms in the molecule
        as an array of atom1, atom2, bond_order

        Returns
        -------
        np.ndarray
            The bonds
        """
        return np.array(
            [
                [bond.atom1.serial_number, bond.atom2.serial_number, bond.order]
                for bond in self.get_bonds()
            ]
        )

    def get_bond_mask(self) -> np.ndarray:
        """
        Get the bonds of the atoms in the molecule
        as a 2D mask where fields with 1 indicate a bond between the atoms
        of row and column.

        Returns
        -------
        np.ndarray
            The bond mask
        """
        n = len(self.get_atoms())
        mask = np.zeros((n, n), dtype=int)
        for bond in self.get_bonds():
            mask[bond.atom1.serial_number - 1, bond.atom2.serial_number - 1] = 1
            mask[bond.atom2.serial_number - 1, bond.atom1.serial_number - 1] = 1
        return mask

    # def infer_missing_atoms(self, _topology=None, _compounds=None):
    #     """
    #     Infer missing atoms in the structure based on a reference topology or compounds database.
    #     By default, if a residue is not available in the topology, the compounds database is used.

    #     Parameters
    #     ----------
    #     _topology
    #         A specific topology to use for referencing.
    #         If None, the default CHARMM topology is used.
    #     _compounds
    #         A specific compounds object to use for referencing.
    #         If None, the default compounds object is used.
    #     """
    #     structural.fill_missing_atoms(self._base_struct, _topology, _compounds)
    #     for atom in self._base_struct.get_atoms():
    #         if atom not in self._AtomGraph:
    #             self._AtomGraph.add_node(atom)
    #     self.infer_bonds()

    def _get_bonds(
        self,
        atom1,
        atom2,
        either_way: bool = True,
    ):
        """
        The core function of `get_bonds` which expects atoms to be provided as Atom objects
        """
        iterable_a = isinstance(atom1, (list, tuple))
        iterable_b = isinstance(atom2, (list, tuple))

        if iterable_a and iterable_b:
            bonds = []
            for a1 in atom1:
                for a2 in atom2:
                    bonds.extend(self._get_bonds(a1, a2, either_way=either_way))
            return bonds

        elif iterable_a:
            bonds = []
            for a1 in atom1:
                bonds.extend(self._get_bonds(a1, atom2, either_way=either_way))
            return bonds

        elif iterable_b:
            bonds = []
            for a2 in atom2:
                bonds.extend(self._get_bonds(atom1, a2, either_way=either_way))
            return bonds

        if atom1 and atom2:
            if either_way:
                return [
                    i
                    for i in self._bonds
                    # ALL OF THESE USED TO BE is COMPARISONS!
                    if (atom1 == i[0] and atom2 == i[1])
                    or (atom1 == i[1] and atom2 == i[0])
                    # if (atom1.full_id == i[0].full_id and atom2.full_id == i[1].full_id)
                    # or (atom1.full_id == i[1].full_id and atom2.full_id == i[0].full_id)
                ]
            else:
                return [
                    i
                    for i in self._bonds
                    if (atom1 == i[0] and atom2 == i[1])
                    # if (atom1.full_id == i[0].full_id and atom2.full_id == i[1].full_id)
                ]
        elif atom1:
            if either_way:
                return [i for i in self._bonds if (atom1 == i[0] or atom1 == i[1])]
            else:
                return [i for i in self._bonds if atom1 == i[0]]
        elif atom2:
            if either_way:
                return [i for i in self._bonds if (atom2 == i[0] or atom2 == i[1])]
            else:
                return [i for i in self._bonds if atom2 == i[1]]
        else:
            raise ValueError("No atom provided")

    def _purge_bonds(self, atom):
        """
        The core function of `purge_bonds` which expects atoms to be provided as Atom objects.
        """
        # the full_id thing seems to prevent some memory leaky-ness
        bonds = (
            idx
            for idx, i in enumerate(self._bonds)
            if atom is i[0] or atom is i[1]
            # if atom.full_id == i[0].full_id or atom.full_id == i[1].full_id
        )
        _bonds = {idx: i for idx, i in enumerate(self._bonds)}
        for b in bonds:
            _bonds.pop(b)
        self._bonds = list(_bonds.values())

    def _add_bond(self, atom1, atom2, order=1):
        """
        Add a bond between two atoms
        This method expects the atoms to be present in the structure!
        """
        if not atom1:
            raise ValueError("Atom1 not found!")
        if not atom2:
            raise ValueError("Atom2 not found!")

        if not self._AtomGraph.has_edge(atom1, atom2):
            bond = base_classes.Bond(atom1, atom2, order)
            self._AtomGraph.add_edge(atom1, atom2, bond_order=order, bond_obj=bond)
            self._bonds.append(bond)
        else:
            self._AtomGraph.edges[atom1, atom2]["bond_order"] = (
                self._AtomGraph.edges[atom1, atom2].get("bond_order", 0) + 1
            )
            self._AtomGraph.edges[atom1, atom2]["bond_obj"].order += 1

    def _set_bond(self, atom1, atom2, order=1):
        """
        The core function of `set_bond` which expects atoms to be provided as Atom objects.
        """
        if not self._AtomGraph.has_edge(atom1, atom2):
            bond = base_classes.Bond(atom1, atom2, order)
            self._AtomGraph.add_edge(atom1, atom2, bond_order=order, bond_obj=bond)
            self._bonds.append(bond)
        else:
            self._AtomGraph.edges[atom1, atom2]["bond_order"] = order
            self._AtomGraph.edges[atom1, atom2]["bond_obj"].order = order

        return self

    def _remove_bond(self, atom1, atom2):  # , either_way: bool = False):
        """
        The core function of `remove_bond` which expects atoms to be provided as Atom objects.
        """
        if self._AtomGraph.has_edge(atom1, atom2):
            bond_obj = self._AtomGraph[atom1][atom2]["bond_obj"]
            self._bonds.remove(bond_obj)
            if self._AtomGraph[atom1][atom2].get("bond_order", 1) == 1:
                self._AtomGraph._locked_edges.discard(bond_obj)
                self._AtomGraph.remove_edge(atom1, atom2)
            else:
                self._AtomGraph[atom1][atom2]["bond_order"] = (
                    self._AtomGraph[atom1][atom2].get("bond_order", 1) - 1
                )
                bond_obj.order -= 1
        else:
            bond_obj = base_classes.Bond(atom1, atom2)
            self._bonds.remove(bond_obj)

        # if either_way:
        #     if bond_obj[::-1] in self._bonds:
        #         self._bonds.discard(b[::-1])
        #     self.locked_bonds.discard(b[::-1])
        #     if self._AtomGraph.has_edge(atom2, atom1):
        #         self._AtomGraph.remove_edge(atom2, atom1)

    def _remove_atoms(self, *atoms):
        """
        The core alternative of `remove_atoms` which expects atoms to be provided as Atom objects.
        """
        for atom in atoms:
            self._purge_bonds(atom)
            self._AtomGraph.remove_node(atom)
            p = atom.get_parent()
            p.detach_child(atom.get_id())
            atom.set_parent(p)

        # I think we should drop this because it's not really necessary...
        # but it might cause problems down the line if we do so and I don't
        # want to work it out now...
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
        # FUTURE TODO: This currently still works with tuples of atoms instead of the new Bond objects.
        # This should be changed to use the Bond objects instead, which can be inverted and do not need to be
        # pre-processed into new lists of tuples.
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

    def _rotate_around_bond(
        self,
        atom1: base_classes.Atom,
        atom2: base_classes.Atom,
        angle: float,
        descendants_only: bool = False,
    ):
        """
        Rotate the structure around a bond. This is the same as the `rotate_around_bond` method,
        but expects the atoms to be provided as Atom objects. And it does not check if a bond is locked.

        Note
        ----
        This function expects the angle to be in RADIANS! Contrary to the `rotate_around_bond` method!
        """
        self._AtomGraph.rotate_around_edge(atom1, atom2, angle, descendants_only)

    def __mod__(self, patch):
        """
        Add a patch to the molecule using the % operator (i.e. mol % patch)
        """
        # we also allow setting functional groups as patches
        # using the operator
        if patch.__class__.__name__ == "FunctionalGroup":
            self._linkage = patch
        else:
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


def infer_search_param(input):
    """
    Infer the search parameter 'by' for the get_atoms/residues etc. methods
    """
    if isinstance(input, (int, np.int64, np.int16, np.int32, np.int8)):
        return "serial"
    elif isinstance(input, str):
        return "id"
    elif isinstance(input, (tuple, list)):
        return "full_id"
    else:
        raise ValueError(f"Could not infer search parameter for {input}")


if __name__ == "__main__":
    f = "/Users/noahhk/GIT/biobuild/__figure_makery__/e8_opt.pdb"
    e = BaseEntity.from_pdb(f)
    # e.infer_bonds(restrict_residues=False)
    # e.get_residue_connections()
    # pass

    import buildamol as bam

    bam.load_amino_acids()
    # ser = molecule("SER")
    ser = BaseEntity.from_compound("SER")

    mol = ser % "LINK" * 6

    print(mol.residues)

    mol.squash()

    print(mol.residues)

    mol.show()

    # import buildamol as bam

    # man = bam.Molecule.from_compound("MAN")
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
