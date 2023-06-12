"""
Functions to infer structural data such as missing atom coordinates or bond connectivity
"""

import pandas as pd
import networkx as nx
from collections import defaultdict
from copy import deepcopy
import warnings
import numpy as np

from Bio.PDB import NeighborSearch
from Bio import SVDSuperimposer
import periodictable as pt

import biobuild.utils.ic as _ic
import biobuild.utils.abstract as _abstract
import biobuild.utils.defaults as defaults
import biobuild.resources as resources
import biobuild.structural.base as base
import biobuild.structural.neighbors as neighbors

element_connectivity = {
    "C": 4,
    "H": 1,
    "O": 2,
    "N": 3,
    "S": 2,
    "P": 5,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "B": 3,
    "Si": 4,
    "Se": 2,
    "Zn": 2,
    "Ca": 2,
    "Mg": 2,
    "Fe": 2,
    "Cu": 1,
    "Mn": 2,
}


class AutoLabel:
    """
    A automatic atom labeller

    Parameters
    ----------
    atom_graph : nx.Graph
        The molecule's atom graph
    """

    def __init__(self, atom_graph):
        self.graph = atom_graph
        self._bond_orders = nx.get_edge_attributes(self.graph, "bond_order")
        self._df = self._make_df()

    @property
    def carbons(self):
        """
        All carbon atoms in the molecule.
        """
        return [n for n in self.graph.nodes if n.element == "C"]

    def relabel(self):
        """
        Relabel the atoms in the molecule

        Returns
        -------
        pd.DataFrame
            A dataframe with the atom objects and their new labels.
        """
        self._parse_carbon_labels()
        self._parse_hetero_labels()
        self._parse_hydrogen_labels()
        self._df["label"] = self._df["element"] + self._df["label"]
        return self._df[["atom", "label"]]

    def _neighbors(self, atoms):
        """
        Get the neighbors of a list of atoms.
        """
        neighbors = []
        neighbor_element_sum = []
        for atom in atoms:
            ndx = 0
            edx = 0
            for neighbor in self.graph.neighbors(atom):
                ndx += 1
                n_num = pt.elements.symbol(neighbor.element).number
                n_num *= self._bond_orders.get((atom, neighbor), 1)
                edx += n_num
            neighbors.append(ndx)
            neighbor_element_sum.append(edx)
        return neighbors, neighbor_element_sum

    def _make_df(self):
        """
        Make a dataframe of the molecule connectivity.
        """

        neighbors, neighbor_element_sum = self._neighbors(self.graph.nodes)
        self._df = pd.DataFrame(
            {
                "atom": list(self.graph.nodes),
                "element": [a.element for a in self.graph.nodes],
                "neighbors": neighbors,
                "neighbor_element_sum": neighbor_element_sum,
            }
        )
        self._df["total"] = self._df.neighbors * self._df.neighbor_element_sum
        self._df["label"] = "none"
        return self._df

    def _parse_c1(self, carbons=None):
        if carbons is not None:
            carbons = self._df[self._df.atom.isin(carbons)]
        else:
            carbons = self._df[self._df.element == "C"]
        carbons = carbons[carbons.label == "none"]
        carbons = carbons.sort_values("total", ascending=False)
        # carbons = carbons.sort_values(
        #     ["neighbors", "neighbor_element_sum"], ascending=False
        # )
        carbons = carbons.reset_index(drop=True)
        return carbons.atom[0]

    def _parse_c_next(self, c_current):
        neighbors = self.graph.neighbors(c_current)
        neighbors = self._df[self._df.atom.isin(neighbors)]
        neighbors = neighbors[neighbors.element == "C"]
        neighbors = neighbors[neighbors.label == "none"]
        neighbors = neighbors.sort_values("total", ascending=False)
        # neighbors = neighbors.sort_values(
        #     ["neighbors", "neighbor_element_sum"], ascending=False
        # )
        neighbors = neighbors.reset_index(drop=True)
        if len(neighbors) == 0:
            return None
        return neighbors.atom[0]

    def _parse_carbon_labels(self):
        carbons = self._df[self._df.element == "C"]
        c1 = self._parse_c1()
        self._df.loc[self._df.atom == c1, "label"] = "1"
        self._c1_row = self._df[self._df.atom == c1]
        idx = 2
        c_current = c1
        while (carbons.label == "none").any():
            c_next = self._parse_c_next(c_current)
            if c_next is None:
                c_next = self._parse_c1(carbons=carbons[carbons.label == "none"].atom)
            self._df.loc[self._df.atom == c_next, "label"] = f"{idx}"
            idx += 1
            c_current = c_next
            carbons = self._df[self._df.element == "C"]

    def _parse_hetero_labels(self):
        _neighbor_connect_dict = {}
        heteros = self._df[(self._df.element != "C") * (self._df.element != "H")]
        while (heteros.label == "none").any():
            for hetero in heteros.atom:
                neighbors = self.graph.neighbors(hetero)
                neighbors = self._df[self._df.atom.isin(neighbors)]
                # remove c1 row
                if len(neighbors[neighbors.element != "H"]) > 1:
                    neighbors = neighbors[neighbors.atom != self._c1_row.atom.iloc[0]]
                neighbors = neighbors[neighbors.label != "none"]
                neighbors = neighbors.sort_values("label", ascending=False)
                # neighbors = neighbors.sort_values(
                #     ["neighbors", "neighbor_element_sum"], ascending=False
                # )
                neighbors = neighbors.reset_index(drop=True)
                if len(neighbors) == 0:
                    continue

                neighbor = neighbors.atom.iloc[-1]
                label = neighbors.label.iloc[-1]
                if neighbor not in _neighbor_connect_dict:
                    _neighbor_connect_dict[neighbor] = [hetero]
                else:
                    _neighbor_connect_dict[neighbor].append(hetero)
                self._df.loc[self._df.atom == hetero, "label"] = label
            heteros = self._df[(self._df.element != "C") * (self._df.element != "H")]
        for _heteros in _neighbor_connect_dict.values():
            idx = 1
            if len(_heteros) > 1:
                for h in _heteros:
                    self._df.loc[self._df.atom == h, "label"] += str(idx)
                    idx += 1

    def _parse_hydrogen_labels(self):
        _neighbor_connect_dict = {}
        hydrogens = self._df[self._df.element == "H"]
        for hydrogen in hydrogens.atom:
            neighbors = self.graph.neighbors(hydrogen)
            neighbors = self._df[self._df.atom.isin(neighbors)]
            neighbors = neighbors[neighbors.label != "none"]
            neighbors = neighbors.sort_values("label", ascending=False)
            # neighbors = neighbors.sort_values(
            #     ["neighbors", "neighbor_element_sum"], ascending=False
            # )
            neighbors = neighbors.reset_index(drop=True)
            if len(neighbors) == 0:
                continue
            neighbor = neighbors.atom.iloc[-1]
            if neighbor not in _neighbor_connect_dict:
                _neighbor_connect_dict[neighbor] = [hydrogen]
            else:
                _neighbor_connect_dict[neighbor].append(hydrogen)
            element = neighbors.element.iloc[-1]
            if element == "C":
                element = ""
            label = element + neighbors.label.iloc[-1]
            self._df.loc[self._df.atom == hydrogen, "label"] = label
        for _hydrogens in _neighbor_connect_dict.values():
            idx = 1
            if len(_hydrogens) > 1:
                for h in _hydrogens:
                    self._df.loc[self._df.atom == h, "label"] += str(idx)
                    idx += 1


def autolabel(molecule):
    """
    Automatically relabel atoms in a structure to match the CHARMM naming scheme.
    Note, this function is not guaranteed to produce the correct labels in all cases,
    validation of the labels is recommended.

    Parameters
    ----------
    molecule : biobuild.core.Molecule
        The molecule that holds the atoms to be relabeled.
        This molecule needs to have bonds assigned or computed.
    """
    labeler = AutoLabel(molecule._AtomGraph)
    df = labeler.relabel()
    for idx in range(len(df)):
        atom = df.atom.iloc[idx]
        label = df.label.iloc[idx]
        atom.id = label

    return molecule


def relabel_hydrogens(molecule):
    """
    Relabel hydrogen atoms in a structure to match the CHARMM naming scheme.

    Parameters
    ----------
    molecule : biobuild.structural.base.Molecule
        The molecule that holds the atoms to be relabeled.
        This molecule needs to have bonds assigned or computed.
    """

    _neighbors_H_dict = {}
    for atom in molecule.get_atoms():
        if not atom.element == "H":
            continue

        _neighbors = molecule.get_neighbors(atom)
        if len(_neighbors) != 1:
            Warning(
                f"Atom {atom} (full_id: {atom.full_id}) has {len(_neighbors)} neighbors, but should have 1!"
            )
            continue

        _neighbor = _neighbors.pop()
        _neighbors_H_dict.setdefault(_neighbor, []).append(atom)

    # -------------------------   IMPORTANT   -------------------------
    # The code below assumes that we are only dealing with default
    # organic molecules whose atoms have single letter element symbols.
    # -------------------------   IMPORTANT   -------------------------

    for _neighbor, hydrogens in _neighbors_H_dict.items():
        if len(hydrogens) == 1:
            if _neighbor.element == "C":
                _n = _neighbor.id[1:]
            else:
                _n = _neighbor.id
            hydrogens[0].id = f"H{_n}"
        else:
            for i, hydrogen in enumerate(hydrogens):
                if _neighbor.element == "C":
                    _n = _neighbor.id[1:]
                else:
                    _n = _neighbor.id
                hydrogen.id = f"H{_n}{i+1}"

    return molecule


def compute_residue_radius(residue):
    """
    Compute the radius of a residue by computing the distance between its center of mass
    and the furthest atom.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        The residue to compute the radius for.

    Returns
    -------
    radius : float
        The radius of the residue.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    radius = np.max(distances)
    return radius


def compute_outlier_atoms(residue, f: float = 1.5):
    """
    Compute which atoms of a residue are especially far away from the residues center of mass.
    This function compute the distances between the center of mass of the residue and its atoms
    and returns all atoms are further than `f * p75` away from the center of mass, where `p75`
    is the 75th percentile of the distances.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        The residue to compute the outlier atoms for.
    f : float
        The factor to multiply the 75th percentile with.

    Returns
    -------
    outlier_atoms : list
        A list of atoms that are considered outliers.
    """
    atoms = list(residue.get_atoms())
    atom_coords = np.array([atom.get_coord() for atom in atoms])
    center = residue.center_of_mass()

    distances = np.linalg.norm(atom_coords - center, axis=1)
    p75 = np.percentile(distances, 75)
    f = f * p75

    outlier_atoms = [atom for atom, distance in zip(atoms, distances) if distance > f]
    return outlier_atoms


def infer_surface_residues(
    structure,
    cutoff: int = 75,
    fraction: float = None,
):
    """
    Infer residues that are likely to be on the surface of the structure
    using the Solvent accessible surface area (SASA) of the structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the surface residues from.
    n_points : int
        The number of points to sample on the surface of the structure.
    cutoff : int
        The cutoff to use for classifying residues as surface residues.
    fraction : float
        The fraction of residues to classify as surface residues. In this case,
        the cutoff is adjusted to match the fraction of residues.

    Returns
    -------
    surface_residues : list
        A list of residues that are likely to be on the surface of the structure.
    """

    sasa = defaults.get_default_instance("bioSASA")
    sasa.compute(structure, level="R")

    sasa_values = np.array([residue.sasa for residue in structure.get_residues()])
    sasa_values = sasa_values / sasa_values.max() * 100

    if fraction is not None:
        cutoff = np.percentile(sasa_values, 100 - fraction * 100)

    surface_residues = [
        residue
        for residue, sasa in zip(structure.get_residues(), sasa_values)
        if sasa > cutoff
    ]
    return surface_residues


def infer_residue_connections(
    structure,
    bond_length: float = None,
    triplet: bool = False,
):
    """
    Infer the connectivity graph of residues from the distances between atoms of residue pairs.
    This will establish only bonds between close-by atoms from non-identical residues.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host at least two Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
    bond_length : float
        The maximum distance between two atoms to be considered a bond.
    triplet : bool
        If True, bonds between atoms of the same residue are also added, if one
        of the atoms is considered bonded to another residue. Like this residue connections
        are described not by a single bond with a pair of atoms, but two bonds with a triplet of atoms.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms from different Residues that are bonded.
    """
    if not triplet:
        if bond_length is None:
            bond_length = defaults.DEFAULT_BOND_LENGTH / 2, defaults.DEFAULT_BOND_LENGTH
        elif isinstance(bond_length, (int, float)):
            bond_length = defaults.DEFAULT_BOND_LENGTH / 2, bond_length
        min_length, max_length = bond_length

        bonds = []
        _seen_residues = set()
        for residue1 in structure.get_residues():
            for residue2 in structure.get_residues():
                if residue1 == residue2:
                    continue
                elif residue2 in _seen_residues:
                    continue

                atoms = list(residue1.get_atoms())
                atoms.extend(residue2.get_atoms())

                _neighbors = NeighborSearch(atoms)
                _neighbors = _neighbors.search_all(radius=max_length)

                _neighbors = (
                    i
                    for i in _neighbors
                    if (i[0].element != "H" and i[1].element != "H")
                    and i[0].get_parent() != i[1].get_parent()
                    and np.linalg.norm(i[0].coord - i[1].coord) > min_length
                )

                bonds.extend(_neighbors)

            _seen_residues.add(residue1)
    else:
        connections = infer_residue_connections(
            structure, bond_length=bond_length, triplet=False
        )
        base_bonds = infer_bonds(
            structure, bond_length=bond_length, restrict_residues=True
        )
        _additional_bonds = []
        for atom1, atom2 in connections:
            _new = (bond for bond in base_bonds if atom1 in bond)
            _additional_bonds.extend(_new)
        bonds = connections + _additional_bonds

    return bonds


def infer_bonds(structure, bond_length: float = None, restrict_residues: bool = True):
    """
    Generate a connectivity graph by inferring bonds from the distance between atoms.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to infer the bonds from. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
    bond_length : float or tuple
        The maximum distance between two atoms to be considered a bond.
        If a tuple is provided, it specifies the minimal and maximal distances between atoms.
    restrict_residues : bool
        If set to `True`, only bonds between atoms of the same residue will be considered.

    Returns
    -------
    bonds : list
        The connectivity graph of the molecule, storing tuples of `Bio.PDB.Atom` objects.
    """
    if bond_length is None:
        bond_length = (defaults.DEFAULT_BOND_LENGTH / 2, defaults.DEFAULT_BOND_LENGTH)
    elif isinstance(bond_length, (int, float)):
        bond_length = (defaults.DEFAULT_BOND_LENGTH / 2, bond_length)
    min_length, max_length = bond_length

    if restrict_residues:
        bonds = []
        for residue in structure.get_residues():
            atoms = list(residue.get_atoms())
            _neighbors = NeighborSearch(atoms)
            bonds.extend(
                _neighbors.search_all(radius=max_length),
            )

    else:
        atoms = list(structure.get_atoms())
        _neighbors = NeighborSearch(atoms)
        bonds = _neighbors.search_all(radius=max_length)

    bonds = [
        i
        for i in bonds
        if not (i[0].element == "H" and i[1].element == "H")
        and np.linalg.norm(i[0].coord - i[1].coord) > min_length
    ]
    bonds = _prune_H_triplets(bonds)
    return bonds


def identify_glycans(structure):
    """
    Identify glycans in a structure.

    Parameters
    ----------
    structure : Bio.PDB.Structure.Structure
        The structure to identify glycans in. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`

    Returns
    -------
    glycans : list
        A list of glycans, each represented by a list of `Bio.PDB.Residue` objects.
    """
    glycans = []
    for residue in structure.get_residues():
        if residue.get_resname() in defaults.GLYCAN_RESIDUES:
            glycans.append(residue)
    return glycans


def apply_standard_bonds(structure, _topology=None):
    """
    Apply bonds according to loaded standard topologies. This will compute a list of tuples with bonded
    atoms from the same residue. It will not infer residue-residue connections! It is possible to provide a single
    atom as input structure, in which case all bonds that involve said atom are returned.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which hosts Atoms:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`
        - `Bio.PDB.Atom`

    _topology : Topology
        The topology to use for applying the bonds.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms that are bonded.
    """
    if structure.level == "R":
        residue = structure

        if _topology is None:
            _topology = resources.get_default_topology()

        if not _topology.has_residue(residue.resname):
            warnings.warn(
                f"[ignoring] No abstract residue found in Topology for {residue.resname}!"
            )
            return

        atoms = residue.child_dict
        abstract_residue = _topology.get_residue(residue.resname)
        bonds = [
            (atoms.get(bond.atom1.id), atoms.get(bond.atom2.id))
            for bond in abstract_residue.bonds
        ]
        bonds = [
            bond for bond in bonds if bond[0] and bond[1]
        ]  # make sure to have no None entries...
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_standard_bonds(residue, _topology)
        bonds = [
            bond
            for bond in bonds
            if bond[0].id == structure.id or bond[1].id == structure.id
        ]
        return bonds

    else:
        bonds = []
        for residue in structure.get_residues():
            bonds.extend(apply_standard_bonds(residue, _topology))
        return bonds


def fill_missing_atoms(structure, _topology=None, _compounds=None):
    """
    Fill missing atoms in residues either based on the loaded standard topology
    or based on the loaded compounds dictionary. By default, if a residue is not
    available in the topology, the compounds database will be used. Note, however,
    that for structures that contain residues in non-standard conformations (e.g. proteins)
    using the compounds database may lead to poor results and is not recommended. In such cases,
    it is recommended to use a molecular dynamics software to fill missing atoms in appropriate coordinates.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which host Residues:
        - `Bio.PDB.Structure`
        - `Bio.PDB.Model`
        - `Bio.PDB.Chain`
        - `Bio.PDB.Residue`

    _topology
        A specific topology to use for references.
    _compounds
        A specific compounds dictionary to use for references.
    """
    if structure.level in set(("S", "M", "C")):
        for residue in structure.get_residues():
            fill_missing_atoms(residue, _topology=_topology, _compounds=_compounds)
        return
    elif structure.level == "A":
        raise ValueError("Cannot fill missing atoms in an Atom object!")
    elif structure.level != "R":
        raise ValueError(f"Invalid level for structure: {structure.level}")

    if not _topology:
        _topology = resources.get_default_topology()
    if not _compounds:
        _compounds = resources.get_default_compounds()

    if _topology.has_residue(structure.resname):
        _infer_missing_via_topology(structure, _topology)
    elif _compounds.has_residue(structure.resname):
        _infer_missing_via_compounds(structure, _compounds)


def compute_internal_coordinates(bonds: list):
    """
    Compute internal coordinates for a structure.

    Parameters
    ----------
    bonds : list
        A list of tuples of atoms that are bonded.
        The atoms must be `Bio.PDB.Atom` objects with coordinates.

    Returns
    -------
    list
        A list of InternalCoordinates
    """
    quartets = neighbors.compute_quartets(bonds)
    ics = []
    for quartet in quartets:
        angle_123 = base.compute_angle(quartet[0], quartet[1], quartet[2])
        angle_234 = base.compute_angle(quartet[1], quartet[2], quartet[3])
        dihedral = base.compute_dihedral(quartet[0], quartet[1], quartet[2], quartet[3])
        l_12 = base.compute_distance(quartet[0], quartet[1])
        # l_23 = base.compute_distance(quartet[1], quartet[2])
        l_34 = base.compute_distance(quartet[2], quartet[3])

        ic = _ic.InternalCoordinates(
            quartet[0].id,
            quartet[1].id,
            quartet[2].id,
            quartet[3].id,
            l_12,
            l_34,
            angle_123,
            angle_234,
            dihedral,
            improper=quartet.improper,
        )
        ics.append(ic)
    return ics


def impute_atom_from_IC(atom, residue, abstract):
    """
    Compute the coordinates of an atom from internal coordinates from an AbstractResidue into a structure.

    Parameters
    ----------
    atom : AbstractAtom
        The atom to compute the coordinates for.
    residue : Bio.PDB.Residue
        The target residue wherein to add the atom.
    abstract : AbstractResidue
        The abstract residue to reference data from.

    Returns
    -------
    success: bool
        Whether the atom could successfully be integrated.
        It may be that a specific atom could not be computed due
        to lacking reference atoms. In this case it may be necessary
        to first impute other atoms.
    """
    success, _coord = compute_atom_coordinate_from_IC(atom, residue, abstract)
    if not success:
        return False

    # create the new atom
    _new_atom = atom.to_biopython()
    _new_atom.set_coord(_coord)

    # add the new atom to the residue
    residue.add(_new_atom)
    return True


def compute_atom_coordinate_from_IC(atom, residue, abstract):
    """
    Compute the coordinates of an atom from internal coordinates from an AbstractResidue into a structure.

    Parameters
    ----------
    atom : AbstractAtom
        The atom to compute the coordinates for.
    residue : Bio.PDB.Residue
        The target residue where other reference atoms are located.
    abstract : AbstractResidue
        The abstract residue to reference internal coordinates from.

    Returns
    -------
    success: bool
        Whether the atom could successfully be integrated.
    coords : np.ndarray
        The coordinates of the atom
    """
    # get the internal coordinates to impute the atom positions from
    # since we can either compute atoms at position 1 or 4, we will specifically
    # search for ICs with these properties...
    _coord_func = compute_atom1_from_others
    _ics = abstract.get_internal_coordinates(atom, None, None, None, mode="partial")
    if len(_ics) == 0:
        _coord_func = compute_atom4_from_others
        _ics = abstract.get_internal_coordinates(None, None, None, atom, mode="partial")
        if len(_ics) == 0:
            # warnings.warn(f"Cannot find suitable IC for {atom}!")
            return False, None

    # try finding an IC that has available reference atoms
    for ic in _ics:
        _ref = ic.get_reference_atoms(residue)
        if len(_ref) == 3:
            break
    else:
        # warnings.warn(f"Cannot find suitable IC for {atom}!")
        return False, None

    # get the coordinates of the reference atoms
    _ref = [i.coord for i in _ref]

    # now use the internal coordinates to compute the new atom coordinate
    _coord = _coord_func(*_ref, ic)

    return True, _coord


def compute_atom1_from_others(coords2, coords3, coords4, ic):
    """
    Compute the coordinates of the first atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : InternalCoordinates
        The internal coordinates

    Returns
    -------
    coords1 : array-like
        The coordinates of the first atom
    """
    if ic.is_proper:
        return base._IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=coords2,
            r=-ic.bond_length_12,
            theta=-np.radians(ic.bond_angle_123),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        _vec = base._IC_to_xyz(
            coords4,
            coords3,
            coords2,
            anchor=np.full(3, 0),
            r=1,
            theta=np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
        if ic.bond_length_13:
            _vec *= ic.bond_length_13
        else:
            BC = np.linalg.norm(coords2 - coords3)
            AC = ic.bond_length_12
            AB = np.sqrt(
                BC**2 + AC**2 - 2 * BC * AC * np.cos(np.radians(ic.bond_angle_123))
            )
            _vec *= AB
        final = coords3 + _vec
        return final


def compute_atom4_from_others(coords1, coords2, coords3, ic):
    """
    Compute the coordinates of the fourth atom from internal coordinates and the coordinates of the other three.

    Parameters
    ----------
    *coords: array-like
        The coordinates of the other three atoms

    ic : InternalCoordinates
        The internal coordinates

    Returns
    -------
    coords4 : array-like
        The coordinates of the fourth atom
    """
    if ic.is_proper:
        return base._IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )
    else:
        return base._IC_to_xyz(
            coords1,
            coords2,
            coords3,
            anchor=coords3,
            r=-ic.bond_length_34,
            theta=-np.radians(ic.bond_angle_234),
            dihedral=np.radians(ic.dihedral),
        )


def _infer_missing_via_topology(residue, _topology):
    """
    Infer missing atoms in a residue through the internal coordinates
    of an abstract residue in a Topology.

    Parameters
    ----------
    residue : Bio.PDB.Residue
        The residue to impute atoms into.
    _topology : Topology
        The Topology to reference abstract residues from.
    """
    if not _topology.has_residue(residue.resname):
        raise ValueError(
            f"No abstract residue found in Topology for {residue.resname}!"
        )
        return
    _abstract = _topology.get_residue(residue.resname)

    # first get all atoms that need to be added
    _atoms_to_add = _abstract.get_missing_atoms(residue)

    # now iterate through the atoms and impute them one by one
    # try to keep atoms that could not be imputed in a queue
    # in hopes of imputing them later, but abort if it's no good...
    _fail_timer = 0
    while len(_atoms_to_add) != 0:
        _atom = _atoms_to_add.pop(0)
        _success = impute_atom_from_IC(_atom, residue, _abstract)
        if not _success:
            _atoms_to_add.append(_atom)
            _fail_timer += 1
            if _fail_timer == len(_atoms_to_add) + 1:
                warnings.warn(
                    f"Could not impute atoms: {[i.id for i in _atoms_to_add]}! Aborting...",
                    RuntimeWarning,
                )
                return
            continue
        _fail_timer = 0


def _infer_missing_via_compounds(residue, _compounds):
    """
    Infer missing atoms in a residue through a reference residue from the
    compounds database.

    Parameters
    ----------
    residue : Bio.PDB.Residue
        The residue to impute atoms into.
    _compounds : PDBECompounds
        The Compounds database to reference abstract residues from.
    """
    if not _compounds.has_residue(residue.resname):
        raise ValueError(
            f"No abstract residue found in Compounds database for {residue.resname}!"
        )
    ref = _compounds.get(residue.resname, return_type="residue")

    # first get all atoms that need to be added
    _atoms_to_add = [i for i in ref.get_atoms() if i.id not in residue.child_dict]
    _ref_atoms = [i.id for i in ref.get_atoms() if i.id in residue.child_dict]

    # now overlay the abstract residue onto the real one
    imposer = SVDSuperimposer.SVDSuperimposer()

    residue_coords = np.array([residue.child_dict[i].coord for i in _ref_atoms])
    ref_coords = np.array([ref.child_dict[i].coord for i in _ref_atoms])

    imposer.set(residue_coords, ref_coords)
    imposer.run()

    rot, tran = imposer.get_rotran()
    for atom in _atoms_to_add:
        atom.coord = np.dot(atom.coord, rot) + tran
        residue.add(atom)


def _prune_H_triplets(bonds):
    """
    Remove and erroneous bonds that connect hydrogens to multiple other atoms.

    Parameters
    ----------
    bonds : list of tuples
        The bonds to prune.

    Returns
    -------
    bonds : list of tuples
        The pruned bonds.
    """
    bonds_with_H = [
        bond for bond in bonds if bond[0].element == "H" or bond[1].element == "H"
    ]
    triplets = neighbors.generate_triplets(bonds_with_H)
    bond_mappings = defaultdict(int)
    for a, b in bonds_with_H:
        bond_mappings[a] += 1
        bond_mappings[b] += 1

    for triplet in triplets:
        if triplet[1].element != "H":
            continue

        non_H1, H, non_H2 = triplet

        e_non_H1 = non_H1.element
        e_non_H2 = non_H2.element

        if bond_mappings[non_H1] > element_connectivity[e_non_H1]:
            bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1
        elif bond_mappings[non_H2] > element_connectivity[e_non_H2]:
            bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1
        else:
            raise ValueError("Could not prune H triplet!")
    return bonds

    # bonds_with_H = [
    #     set(bond) for bond in bonds if bond[0].element == "H" or bond[1].element == "H"
    # ]
    # bonds = [set(bond) for bond in bonds]

    # for bond in bonds_with_H:
    #     if bond not in bonds:
    #         continue
    #     partners = [
    #         i
    #         for i in bonds_with_H
    #         if i != bond
    #         and not i.isdisjoint(bond)
    #         and i.intersection(bond).pop().element == "H"
    #     ]
    #     if len(partners) == 0:
    #         continue
    #     for partner in partners:
    #         a, b = bond
    #         not_H_1 = b if a.element == "H" else a
    #         not_H_2 = partner.difference(bond).pop()

    #         not_H_1_allowed = element_connectivity[not_H_1.element]
    #         not_H_1_have = sum(1 for _b in bonds if not_H_1 in _b)
    #         if not_H_1_have > not_H_1_allowed:
    #             bonds.remove(bond)
    #             continue

    #         not_H_2_allowed = element_connectivity[not_H_2.element]
    #         not_H_2_have = sum(1 for _b in bonds if not_H_2 in _b)
    #         if not_H_2_have > not_H_2_allowed:
    #             bonds.remove(partner)
    #             continue

    # bonds = [tuple(i) for i in bonds]
    # return bonds


if __name__ == "__main__":
    import biobuild as bb
    from scipy.spatial import distance
    import time

    # def infer_bonds2(structure, bond_length):
    #     min_length, max_length = bond_length

    #     atoms = np.array([atom for atom in structure.get_atoms()])
    #     elements = np.array([atom.element for atom in atoms])
    #     coords = np.stack([atom.coord for atom in atoms])

    #     dists = distance.cdist(coords, coords, metric="euclidean")

    #     mask = np.logical_and(dists < max_length, dists > min_length)

    #     bonds = []
    #     bond_set = []
    #     for i in range(len(atoms)):
    #         for j in np.where(mask[i])[0]:
    #             if elements[i] == "H" and elements[j] == "H":
    #                 continue
    #             b = (atoms[i], atoms[j])
    #             _s = set(b)
    #             if _s in bond_set:
    #                 continue
    #             bonds.append(b)
    #             bond_set.append(_s)

    #     bonds = _prune_H_triplets(bonds)
    #     return bonds

    # repeats = 4
    # print("size, infer_bonds, infer_bonds2")
    # for size in (5, 10, 50, 100, 200):
    #     for i in range(repeats):
    #         glc = bb.Molecule.from_compound("GLC")
    #         glc.repeat(size, "14bb")
    #         glc.purge_bonds()

    #         start = time.time()
    #         bonds = infer_bonds(
    #             glc._base_struct, bond_length=(0.8, 1.6), restrict_residues=False
    #         )
    #         end = time.time()
    #         old_length = len(bonds)
    #         t1 = end - start

    #         glc.purge_bonds()
    #         start = time.time()
    #         bonds = infer_bonds2(glc._base_struct, (0.8, 1.6))
    #         # for b in bonds:
    #         #     glc.add_bond(*b)
    #         # glc.view()
    #         # glc.purge_bonds()
    #         # bonds = _prune_H_triplets(bonds)
    #         # for b in bonds:
    #         #     glc.add_bond(*b)
    #         # glc.view()
    #         end = time.time()
    #         new_length = len(bonds)
    #         t2 = end - start
    #         print(f"{size}, {t1}, {t2}")

    #         # print("old: ", old_length, "new: ", new_length)
    #         pass
