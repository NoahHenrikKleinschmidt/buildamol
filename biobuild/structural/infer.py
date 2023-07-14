"""
Functions to infer structural data such as missing atom coordinates, bond connectivity, or atom labels.
"""

import re
import pandas as pd
import networkx as nx
from collections import defaultdict
import warnings
import numpy as np

from Bio.PDB import NeighborSearch
import periodictable as pt

import biobuild.utils.ic as _ic
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
        self._cycles = nx.cycle_basis(self.graph)
        self._df_all = self._make_df()
        self._df = None

    @property
    def carbons(self):
        """
        All carbon atoms in the molecule.
        """
        return [n for n in self.graph.nodes if n.element == "C"]

    def autolabel(self):
        """
        Generate labels for the atoms in the molecule

        Returns
        -------
        pd.DataFrame
            A dataframe with the atom objects and their new labels.
        """
        for res_df in self._df_all.groupby("residue"):
            self._df = res_df[1]
            self._parse_carbon_labels()
            self._parse_hetero_labels()
            self._parse_hydrogen_labels()
            self._df["label"] = self._df["element"] + self._df["label"]
            self._final_vet()
            self._df_all.loc[self._df.index, "label"] = self._df["label"]

        return self._df_all[["atom", "label"]]

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
                n_num = pt.elements.symbol(neighbor.element.title()).number
                if n_num != 1 and n_num != 6:
                    n_num *= 10
                n_num *= self._bond_orders.get((atom, neighbor), 1)
                edx += n_num
            neighbors.append(ndx)
            neighbor_element_sum.append(edx)
        return neighbors, neighbor_element_sum

    def _in_cycle(self, atom):
        """
        Check if a list of atoms is in a cycle.
        """
        for cycle in self._cycles:
            if set(atom).issubset(cycle):
                return True
        return False

    def _make_df(self):
        """
        Make a dataframe of the molecule connectivity.
        """

        neighbors, neighbor_element_sum = self._neighbors(self.graph.nodes)
        in_cycle = [self._in_cycle([a]) for a in self.graph.nodes]
        self._df = pd.DataFrame(
            {
                "atom": list(self.graph.nodes),
                "element": [a.element.title() for a in self.graph.nodes],
                "neighbors": neighbors,
                "neighbor_element_sum": neighbor_element_sum,
                "in_cycle": in_cycle,
                "residue": [a.get_parent().id[1] for a in self.graph.nodes],
            }
        )
        self._df["total"] = (
            self._df.neighbors
            * self._df.neighbor_element_sum
            * (10 * self._df.in_cycle + 1)
        )
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
        c1 = self._parse_c1()
        self._df.loc[self._df.atom == c1, "label"] = "1"
        self._c1_row = self._df[self._df.atom == c1]
        idx = 2
        c_current = c1
        carbons = self._df[self._df.element == "C"]
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

    def _final_vet(self):
        _label_counts = self._df.label.value_counts()
        _label_counts = _label_counts[_label_counts > 1]
        _label_counts = _label_counts.sort_index()
        _label_counts = _label_counts.reset_index()
        _label_counts.columns = ["label", "count"]
        _label_counts = _label_counts.sort_values("count", ascending=False)
        _label_counts = _label_counts.reset_index(drop=True)
        for idx in range(len(_label_counts)):
            label = _label_counts.label.iloc[idx]
            count = _label_counts["count"].iloc[idx]
            if count == 1:
                continue
            atoms = self._df[self._df.label == label].atom
            for atom in atoms:
                self._df.loc[self._df.atom == atom, "label"] += str(idx)


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
    df = labeler.autolabel()
    # df.loc[:, "residue"] = df.atom.apply(lambda x: x.get_parent())
    # _ = df.residue.apply(lambda x: x.child_dict.clear())
    for idx in range(len(df)):
        atom = df.atom.iloc[idx]
        label = df.label.iloc[idx]
        # p = df.residue.iloc[idx]
        # p.child_dict[label] = atom
        atom.id = label
        atom.name = label

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


def vet_structure(
    molecule, clash_range: tuple = (0.6, 1.7), angle_range: tuple = (90, 180)
) -> bool:
    """
    Check for basic structure integrity.
    This will return False if there are clashes within the structure,
    or invalid angles.

    Parameters
    ----------
    molecule : Molecule
        A biobuild Molecule
    clash_range : tuple
        The minimal and maximal allowed distances between bonded atoms (in Angstrom).
    angle_range : tuple
        The minimal and maximal allowed angle between a triplet of adjacent bonded atoms (in degrees).

    Returns
    -------
    bool
        True if the structure is free from any obstacles,
        False otherwise.
    """
    for a, b in molecule.get_bonds():
        d = base.compute_distance(a, b)
        if not clash_range[0] <= d <= clash_range[1]:
            return False
    for angle in molecule.angles.values():
        if not angle_range[0] <= angle <= angle_range[1]:
            return False
    return True


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


def _atom_from_residue(id, residue):
    return next((atom for atom in residue.get_atoms() if atom.id == id), None)


def apply_reference_bonds(structure, _compounds=None):
    """
    Apply bonds according to loaded reference compounds. This will compute a list of tuples with bonded
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

    _compounds : PDBECompounds
        The reference compounds to use for bond inference. If not provided, the default compounds are used.

    Returns
    -------
    bonds : list
        A list of tuples of Atoms that are bonded.
    """
    if _compounds is None:
        _compounds = resources.get_default_compounds()

    if structure.level == "R":
        residue = structure

        if not _compounds.has_residue(residue.resname):
            warnings.warn(
                f"[ignoring] No reference residue found in Compounds for {residue.resname}!"
            )
            return []

        ref = _compounds.get(residue.resname)
        bonds = ((a.id, b.id) for a, b in ref.get_bonds())
        bonds = (
            (_atom_from_residue(a, residue) and _atom_from_residue(b, residue))
            for a, b in bonds
        )
        bonds = [
            bond for bond in bonds if bond[0] and bond[1]
        ]  # make sure to have no None entries...
        return bonds

    elif structure.level == "A":
        residue = structure.get_parent()
        bonds = apply_reference_bonds(residue, _compounds)
        bonds = [
            bond
            for bond in bonds
            if bond[0].id == structure.id or bond[1].id == structure.id
        ]
        return bonds

    else:
        bonds = []
        for residue in structure.get_residues():
            bonds.extend(apply_reference_bonds(residue, _compounds))
        return bonds


def atoms_in_area(structure, center, radius):
    """
    Get all atoms in a given area around a center point.

    Parameters
    ----------
    structure : Structure
        The structure to search for atoms.
    center : np.ndarray
        The center point of the area to search.
    radius : float
        The radius of the area to search.

    Yields
    ------
    Atom
        The atoms in the area.
    """
    for atom in structure.get_atoms():
        if np.linalg.norm(atom.coord - center) < radius:
            yield atom


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
        l_13 = base.compute_distance(quartet[0], quartet[2])
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
            l_13 if quartet.improper else None,
            improper=quartet.improper,
        )
        ics.append(ic)
    return ics


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
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif bond_mappings[non_H2] > element_connectivity[e_non_H2]:
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1

        elif _H_id_match(H, non_H1):
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif _H_id_match(H, non_H2):
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1

        elif _H_dist_match(H, non_H1):
            if triplet[:2] in bonds:
                bonds.remove(triplet[:2])
            bond_mappings[non_H1] -= 1

        elif _H_dist_match(H, non_H2):
            if triplet[1:] in bonds:
                bonds.remove(triplet[1:])
            bond_mappings[non_H2] -= 1
        else:
            warnings.warn(f"Could not prune H triplet! ({triplet=})", RuntimeWarning)
    return bonds


def _H_dist_match(H, non_H):
    """
    Check if the distance between a hydrogen and a non-hydrogen atom may suggest
    they are bonded.

    Parameters
    ----------
    H : Bio.PDB.Atom
        The hydrogen atom.
    non_H : Bio.PDB.Atom
        The non-hydrogen atom.

    Returns
    -------
    bool
        True if the hydrogen may be bonded to the non-hydrogen atom.
    """
    d = np.linalg.norm(H.coord - non_H.coord)
    if non_H.element == "C":
        return 0.94 < d < 1.1
    elif non_H.element == "S":
        return 1.15 < d < 1.42
    elif non_H.element in ("O", "N"):
        return 0.94 < d < 1.07
    else:
        return 0.7 < d < 1.4


def _H_id_match(H, non_H):
    """
    Check if the id of a hydrogen may suggest it belongs to a particular non-H atom.

    Parameters
    ----------
    H : Bio.PDB.Atom
        The hydrogen atom.
    non_H : Bio.PDB.Atom
        The non-hydrogen atom.

    Returns
    -------
    bool
        True if the hydrogen may belong to the non-hydrogen atom.
    """
    if (
        non_H.element == "C"
        and re.match("C\d.*", non_H.id.upper()) is not None
        and re.match("H\d.*", H.id.upper()) is None
    ):
        return False
    elif non_H.element != "C" and re.match("H\d.*", H.id.upper()) is not None:
        return False
    return re.search(non_H.id[1:].upper(), H.id[1:].upper()) is not None


if __name__ == "__main__":
    import biobuild as bb

    mol = bb.Molecule.from_pubchem(
        "3-Deoxy-D-Lyxo-Heptopyran-2-ularic Acid"
    )  # ("/Users/noahhk/GIT/biobuild/support/examples/man9.pdb")
    autolabel(mol)
    print(mol.atoms)
