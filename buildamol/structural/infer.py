"""
Functions to infer structural data such as missing atom coordinates, bond connectivity, or atom labels.
"""

import re
import pandas as pd
import networkx as nx
from collections import defaultdict
import warnings
import numpy as np
from scipy.spatial.distance import cdist

from Bio.PDB import NeighborSearch
import periodictable as pt

import buildamol.utils.ic as _ic
import buildamol.utils.defaults as defaults
import buildamol.resources as resources
import buildamol.structural.base as base
import buildamol.structural.neighbors as neighbors
import buildamol.structural.geometry as geometry

# might cause circular imports!
import buildamol.base_classes as base_classes

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

element_to_hydrogen_bond_lengths = {
    "C": 1.09,
    "N": 1.01,
    "O": 0.96,
    "S": 1.04,
    "P": 1.10,
    "F": 0.92,
    "Cl": 0.99,
    "Br": 1.14,
    "I": 1.33,
    "B": 1.19,
    "Si": 1.17,
    "Se": 1.17,
    "Zn": 1.31,
    "Ca": 1.74,
    "Mg": 1.36,
    "Fe": 1.24,
    "Cu": 1.28,
    "Mn": 1.20,
}

single_bond_lengths = {
    "C": {
        "C": 1.54,
        "N": 1.47,
        "O": 1.43,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 1.09,
    },
    "N": {
        "C": 1.47,
        "N": 1.45,
        "O": 1.43,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 1.01,
    },
    "O": {
        "C": 1.43,
        "N": 1.43,
        "O": 1.34,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 0.96,
    },
    "S": {
        "C": 1.81,
        "N": 1.81,
        "O": 1.81,
        "S": 2.05,
        "P": 2.05,
        "F": 1.81,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.04,
    },
    "P": {
        "C": 1.87,
        "N": 1.87,
        "O": 1.87,
        "S": 2.05,
        "P": 2.05,
        "F": 1.87,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.10,
    },
    "F": {
        "C": 1.35,
        "N": 1.35,
        "O": 1.35,
        "S": 1.81,
        "P": 1.87,
        "F": 1.35,
        "Cl": 1.77,
        "Br": 1.94,
        "I": 2.14,
        "H": 0.92,
    },
    "Cl": {
        "C": 1.77,
        "N": 1.77,
        "O": 1.77,
        "S": 2.05,
        "P": 2.05,
        "F": 1.77,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 0.99,
    },
    "Br": {
        "C": 1.94,
        "N": 1.94,
        "O": 1.94,
        "S": 2.05,
        "P": 2.05,
        "F": 1.94,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.14,
    },
    "I": {
        "C": 2.14,
        "N": 2.14,
        "O": 2.14,
        "S": 2.05,
        "P": 2.05,
        "F": 2.14,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
        "H": 1.33,
    },
    "H": {
        "C": 1.09,
        "N": 1.01,
        "O": 0.96,
        "S": 1.04,
        "P": 1.10,
        "F": 0.92,
        "Cl": 0.99,
        "Br": 1.14,
        "I": 1.33,
    },
}

double_bond_lengths = {
    "C": {
        "C": 1.34,
        "N": 1.30,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.34,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "N": {
        "C": 1.30,
        "N": 1.28,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.30,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "O": {
        "C": 1.21,
        "N": 1.21,
        "O": 1.20,
        "S": 2.05,
        "P": 2.05,
        "F": 1.21,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "S": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "P": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "F": {
        "C": 1.34,
        "N": 1.30,
        "O": 1.21,
        "S": 2.05,
        "P": 2.05,
        "F": 1.34,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "Cl": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "Br": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
    "I": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.28,
        "P": 2.28,
        "F": 2.05,
        "Cl": 2.28,
        "Br": 2.28,
        "I": 2.28,
    },
}

triple_bond_lengths = {
    "C": {
        "C": 1.20,
        "N": 1.16,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.20,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "N": {
        "C": 1.16,
        "N": 1.14,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.16,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "O": {
        "C": 1.13,
        "N": 1.13,
        "O": 1.12,
        "S": 2.05,
        "P": 2.05,
        "F": 1.13,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
    "S": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.51,
        "P": 2.51,
        "F": 2.05,
        "Cl": 2.51,
        "Br": 2.51,
        "I": 2.51,
    },
    "P": {
        "C": 2.05,
        "N": 2.05,
        "O": 2.05,
        "S": 2.51,
        "P": 2.51,
        "F": 2.05,
        "Cl": 2.51,
        "Br": 2.51,
        "I": 2.51,
    },
    "F": {
        "C": 1.20,
        "N": 1.16,
        "O": 1.13,
        "S": 2.05,
        "P": 2.05,
        "F": 1.20,
        "Cl": 2.05,
        "Br": 2.05,
        "I": 2.05,
    },
}

bond_length_by_order = {
    1: single_bond_lengths,
    2: double_bond_lengths,
    3: triple_bond_lengths,
}


def atomic_number(element: str):
    """
    Get the atomic number of an element.

    Parameters
    ----------
    element : str
        The element symbol.

    Returns
    -------
    int
        The atomic number of the element.
    """
    return pt.elements.symbol(element.title()).number


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
        self._element_counter = defaultdict(int)

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

    @staticmethod
    def hydrogen_neighbors(atom):
        """
        Generate automated labels for an atom's possible hydrogen neighbors.

        Parameters
        ----------
        atom : Atom
            The atom

        Returns
        -------
        list
            A list of possible hydrogen labels

        Examples
        --------
        >>> from buildamol import Atom
        >>> a = Atom(id="C1")
        >>> AutoLabel.hydrogen_neighbors(a)
        ['H1', 'H11', 'H12', 'H13', 'H14']
        >>> a = Atom(id="OXT")
        >>> AutoLabel.hydrogen_neighbors(a)
        ['HOXT', 'HOXT1', 'HOXT2']
        """
        connectivity = element_connectivity.get(atom.element, 0)
        if atom.element == "C":
            return [
                f"H{atom.id[1:]}",
                *(f"H{atom.id[1:]}{i+1}" for i in range(connectivity)),
            ]
        return [f"H{atom.id}", *(f"H{atom.id}{i+1}" for i in range(connectivity))]

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
                if (
                    len(neighbors[(neighbors.element != "H")]) > 1
                    and len(neighbors[neighbors.element == "C"]) > 1
                ):
                    neighbors = neighbors[neighbors.atom != self._c1_row.atom.iloc[0]]
                neighbors = neighbors[neighbors.label != "none"]

                if len(neighbors[neighbors.element == "C"]) > 0:
                    neighbors = neighbors[neighbors.element == "C"]
                    use_blank_label = True
                else:
                    use_blank_label = False

                # used to by: sort by "label", and ascending=False
                neighbors = neighbors.sort_values("total", ascending=True)

                neighbors = neighbors.reset_index(drop=True)
                if len(neighbors) == 0:
                    continue

                neighbor = neighbors.atom.iloc[-1]
                label = neighbors.label.iloc[-1]
                if neighbor.element not in _neighbor_connect_dict:
                    _neighbor_connect_dict[neighbor] = {hetero.element: [hetero]}
                else:
                    _neighbor_connect_dict[neighbor][hetero.element].append(hetero)

                if not use_blank_label:
                    self._element_counter[hetero.element] += 1
                    if not label[-1].isdigit():
                        label = label[:-1]
                    label += chr(self._element_counter[hetero.element] + 64)

                self._df.loc[self._df.atom == hetero, "label"] = label
            heteros = self._df[(self._df.element != "C") * (self._df.element != "H")]
        for _heteros in _neighbor_connect_dict.values():
            if len(_heteros) > 1:
                # idx = 1
                for h in _heteros.values():
                    for idx, atom in enumerate(h):
                        self._df.loc[self._df.atom == h, "label"] += str(idx + 1)
                    # idx += 1

    def _parse_hydrogen_labels(self):
        _neighbor_connect_dict = {}
        _hydrogen_unique_label_dict = {}
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
            mask = self._df.atom == hydrogen
            self._df.loc[mask, "label"] = label
        for neighbor, _hydrogens in _neighbor_connect_dict.items():
            if len(_hydrogens) == 1:
                continue
            idx = 1
            for h in _hydrogens:
                mask = self._df.atom == h
                if len(self._df.loc[mask, "label"].values[0]) >= 2:
                    ref = str(self._df.loc[mask, "label"].str[:2].values[0])
                    _hydrogen_unique_label_dict[ref] = (
                        _hydrogen_unique_label_dict.get(ref, 0) + 1
                    )
                    self._df.loc[mask, "label"] = ref + chr(
                        64 + _hydrogen_unique_label_dict[ref]
                    )
                else:
                    self._df.loc[mask, "label"] += str(idx)

                # new = self._df.loc[mask, "label"] + str(idx)
                # self._df.loc[mask, "label"] += str(idx)

                #     self._df.loc[mask, "label"] = self._df.loc[mask, "label"].str[
                #         :2
                #     ] + chr(64 + idx)
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
            atoms = self._df[self._df.label == label].atom
            for j, atom in enumerate(atoms):
                self._df.loc[self._df.atom == atom, "label"] += str(j + 1)


def autolabel(molecule):
    """
    Automatically relabel atoms in a structure to match the CHARMM naming scheme.
    Note, this function is not guaranteed to produce the correct labels in all cases,
    validation of the labels is recommended.

    Parameters
    ----------
    molecule : buildamol.core.Molecule
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


class Hydrogenator:
    """
    A class to automatically add hydrogen atoms to organic molecules.

    Note
    ----
    This is designed specifically to infer hydrogens for organic molecules.
    So it can infer hydrogen coordinates for CHNOPS atoms, but may not work reliably for non organic molecules.

    """

    tetrahedral = geometry.Tetrahedral()
    trigonal_planar = geometry.TrigonalPlanar()
    linear = geometry.Linear()

    __geometries__ = {
        (1, 4): tetrahedral,
        (1, 3): tetrahedral,
        (1, 2): trigonal_planar,
        (2, 3): trigonal_planar,
        (2, 4): trigonal_planar,
        (3, 4): linear,
        (2, 5): tetrahedral,  # an organic phosphate
    }

    # some default bond length for hydrogen atoms
    _bond_length = 1.05

    def infer_hydrogens(
        self, molecule: "Molecule", bond_length: float = 1
    ) -> "Molecule":
        """
        Add hydrogen atoms to a molecule.

        Parameters
        ----------
        molecule : buildamol.core.Molecule
            The molecule to add hydrogen atoms to.
        bond_length : float
            The bond length to use for the hydrogen atoms.

        Returns
        -------
        Molecule
            The molecule with hydrogen atoms added.
        """
        self._molecule = molecule
        self._bond_length = bond_length

        for atom in self._molecule._atoms:

            if atom.element == "H":
                continue

            connectivity = element_connectivity.get(atom.element, 0)
            if connectivity == 0:
                continue

            bonds = self._molecule.get_bonds(atom)
            free_slots = (
                connectivity - sum(b.order for b in bonds) - (atom.pqr_charge or 0)
            )

            if free_slots > 0:
                neighbors = set(j for i in bonds for j in i if j != atom)
                self._add_hydrogens(
                    atom=atom,
                    neighbors=neighbors,
                    free_slots=free_slots,
                    bond_order=max(i.order for i in bonds),
                    connectivity=connectivity,
                )

        return self._molecule

    def add_hydrogens(self, atom, _molecule=None):
        """
        Add hydrogens to one particular atom.
        """
        self._molecule = _molecule or self._molecule
        connectivity = element_connectivity.get(atom.element, 0)
        if connectivity == 0:
            return

        bonds = self._molecule.get_bonds(atom)
        free_slots = connectivity - sum(b.order for b in bonds) - (atom.pqr_charge or 0)
        if free_slots > 0:
            neighbors = set(j for i in bonds for j in i if j != atom)
            self._add_hydrogens(
                atom=atom,
                neighbors=neighbors,
                free_slots=free_slots,
                bond_order=max(i.order for i in bonds),
                connectivity=connectivity,
            )

    def _add_hydrogens(self, atom, neighbors, free_slots, bond_order, connectivity):
        """
        Add hydrogen atoms to a molecule.

        Parameters
        ----------
        atom : Atom
            The atom to add hydrogen atoms to.
        neighbors : set
            The eighbors the atom has.
        free_slots : int
            The number of free slots the atom has.
        bond_order : int
            The highest bond order of all bonds that connects the atom to a neighbor.
        connectivity: int
            The theoretically possible number of bonds the atom can form.

        Returns
        -------
        list
            A list of hydrogen atoms.
        list
            A list of bonds to the atom for each hydrogen.
        """
        geometry = Hydrogenator.__geometries__.get((bond_order, connectivity), None)
        if geometry is None:
            return

        _neighbors = list(neighbors)[: geometry.max_points - 1]
        out = geometry.make_coords(atom, *_neighbors, length=self._bond_length)[
            len(_neighbors) :
        ]

        labels = AutoLabel.hydrogen_neighbors(atom)
        if free_slots > 1:
            labels.pop(0)

        # technically this is convenient because
        # it ensures that we don't have multiple Hs
        # with the same id present
        # labels = list(set(labels) - set(i.id for i in neighbors))

        # now figure out which of the coordinates already belong to a neighbor
        # and filter for the free locations
        if len(out) > 1:
            neighbor_coords = np.array([i.coord for i in neighbors])
            d = cdist(out, neighbor_coords)
            out = [out[i] for i in range(len(out)) if not np.any(d[i] < 0.95)]

        if len(out) < free_slots:
            raise ValueError(
                f"Could not add the correct number of hydrogen atoms to {atom}! Expected {free_slots}, got {len(out)}."
            )

        Hs = [
            base_classes.Atom(labels[i], out[i], element="H") for i in range(free_slots)
        ]
        bonds = [base_classes.Bond(atom, H, 1) for H in Hs]

        # now adjust the bond lengths to match the elements
        length = element_to_hydrogen_bond_lengths.get(atom.element, self._bond_length)
        for b in bonds:
            base.adjust_bond_length(b, length)

        self._molecule.add_atoms(*Hs, residue=atom.parent)
        self._molecule.add_bonds(*bonds)


def adjust_protonation(molecule, atom, new_charge):
    """
    Adjust the protonation state of an atom in a molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to adjust the protonation state in.
    atom : Atom
        The atom to adjust the protonation state of.
    new_charge : int
        The new charge of the atom.
    """
    if new_charge == atom.pqr_charge:
        return

    connectivity = element_connectivity.get(atom.element, 0)
    if connectivity == 0:
        raise ValueError(
            f"Cannot adjust protonation state of {atom}. No connectivity information available for element {atom.element}."
        )

    if new_charge >= 0:
        H = Hydrogenator()
        if new_charge == 0:
            # remove all hydrogens
            hydrogens = molecule.get_hydrogens(atom)
            molecule.remove_atoms(*hydrogens)
            atom.pqr_charge = new_charge
            H.add_hydrogens(atom, molecule)
        else:
            # little hack here
            element_connectivity[atom.element] += new_charge
            H.add_hydrogens(atom, molecule)
            atom.pqr_charge = new_charge
            element_connectivity[atom.element] = connectivity
    else:
        hydrogens = tuple(molecule.get_hydrogens(atom))
        if len(hydrogens) < abs(new_charge):
            raise ValueError(
                f"Cannot adjust protonation state of {atom} to {new_charge}. Not enough hydrogens present."
            )
        molecule.remove_atoms(*hydrogens[: abs(new_charge)])
        atom.pqr_charge = new_charge

    return molecule


def relabel_hydrogens(molecule):
    """
    Relabel hydrogen atoms in a structure to match the CHARMM naming scheme.

    Parameters
    ----------
    molecule : buildamol.structural.base.Molecule
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


def change_bond_order(molecule, atom1, atom2, order):
    """
    Change the bond order between two atoms. This will also adjust the protonation state and add or remove hydrogens if necessary.
    This function will not, however, adjust the overall geometry of the molecule. It can currently only adjust the local geometry of the atoms involved and their hydrogen partners.

    Parameters
    ----------
    molecule : Molecule
        The molecule to change the bond order in.
    atom1 : Atom
        The first atom of the bond.
    atom2 : Atom
        The second atom of the bond.
    order : int
        The new bond order. This can be either 1, 2, or 3.
    """
    # check if the bond order is valid
    if order not in (1, 2, 3):
        raise ValueError(f"Invalid bond order {order}.")

    bond = molecule.get_bond(atom1, atom2)
    if bond is None:
        raise ValueError(f"No bond found between {atom1} and {atom2}.")

    # check if the bond order is the same
    if bond.order == order:
        return

    # now remove all hydrogens from the atoms
    hydrogens = []
    for atom in (atom1, atom2):
        for hydrogen in molecule.get_neighbors(atom):
            if hydrogen.element == "H":
                hydrogens.append(hydrogen)
    molecule.remove_atoms(*hydrogens)

    # set the new bond order
    bond.order = order

    # adjust the bond length
    length = bond_length_by_order[order][atom1.element][atom2.element]
    base.adjust_bond_length(bond, length)

    # now add the hydrogens back
    H = Hydrogenator()
    H.add_hydrogens(atom1, molecule)
    H.add_hydrogens(atom2, molecule)

    return molecule


def find_equatorial_substituents(molecule):
    """
    Find all atoms in a molecule that are equatorial to some ring within the molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in

    Returns
    -------
    list
        A list of atoms that are equatorial to some ring.
    """
    rings = molecule._AtomGraph.find_cycles()
    equatorial_atoms = []

    for ring in rings:
        ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
        ring_mean = np.array(ring_mean)

        for atom in ring:
            neighbors = molecule.get_neighbors(atom)
            if len(neighbors) < 3:
                continue

            neighbors = list(neighbors)
            vec = atom.get_coord() - ring_mean

            for neighbor in neighbors:
                if neighbor in ring:
                    continue

                neighbor_vec = neighbor.get_coord() - atom.get_coord()
                angle = np.arccos(
                    np.dot(vec, neighbor_vec)
                    / (np.linalg.norm(vec) * np.linalg.norm(neighbor_vec))
                )
                angle = np.degrees(angle)
                if angle < 45:
                    equatorial_atoms.append(neighbor)
                    break

    return equatorial_atoms


def find_axial_substituents(molecule):
    """
    Find all atoms in a molecule that are axial to some ring within the molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in

    Returns
    -------
    list
        A list of atoms that are axial to some ring.
    """
    rings = molecule._AtomGraph.find_cycles()
    axial_atoms = []

    for ring in rings:
        ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
        ring_mean = np.array(ring_mean)

        for atom in ring:
            neighbors = molecule.get_neighbors(atom)
            if len(neighbors) < 3:
                continue

            neighbors = list(neighbors)
            vec = atom.get_coord() - ring_mean

            for neighbor in neighbors:
                if neighbor in ring:
                    continue

                neighbor_vec = neighbor.get_coord() - atom.get_coord()
                angle = np.arccos(
                    np.dot(vec, neighbor_vec)
                    / (np.linalg.norm(vec) * np.linalg.norm(neighbor_vec))
                )
                angle = np.degrees(angle)
                if angle > 55:
                    axial_atoms.append(neighbor)
                    break

    return axial_atoms


def find_equatorial_hydrogens(molecule):
    """
    Find the equatorial hydrogen atoms in a molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in

    Returns
    -------
    list
        The equatorial hydrogen atoms.
    """
    rings = molecule._AtomGraph.find_cycles()
    equatorial_Hs = []

    # v = molecule.draw()

    for ring in rings:

        ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
        ring_mean = np.array(ring_mean)

        for atom in ring:

            neighbors = molecule.get_neighbors(atom)
            if len(neighbors) < 3:
                continue

            neighbors = list(neighbors)
            vec = atom.get_coord() - ring_mean

            Hs = [i for i in neighbors if i.element == "H"]
            if len(Hs) < 1:
                continue

            for H in Hs:
                H_vec = H.get_coord() - atom.get_coord()
                angle = np.arccos(
                    np.dot(vec, H_vec) / (np.linalg.norm(vec) * np.linalg.norm(H_vec))
                )
                angle = np.degrees(angle)
                if angle < 45:
                    equatorial_Hs.append(H)
                    break

    # v.draw_atoms(*equatorial_Hs, colors="red")
    # v.show()

    return equatorial_Hs


def find_axial_hydrogens(molecule):
    """
    Find the axial hydrogen atoms in a molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in

    Returns
    -------
    list
        The axial hydrogen atoms.
    """
    rings = molecule._AtomGraph.find_cycles()
    axial_Hs = []

    # v = molecule.draw()

    for ring in rings:

        ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
        ring_mean = np.array(ring_mean)

        for atom in ring:

            neighbors = molecule.get_neighbors(atom)
            if len(neighbors) < 3:
                continue

            neighbors = list(neighbors)
            vec = atom.get_coord() - ring_mean

            Hs = [i for i in neighbors if i.element == "H"]
            if len(Hs) < 1:
                continue

            for H in Hs:
                H_vec = H.get_coord() - atom.get_coord()
                angle = np.arccos(
                    np.dot(vec, H_vec) / (np.linalg.norm(vec) * np.linalg.norm(H_vec))
                )
                angle = np.degrees(angle)
                if angle > 55:
                    axial_Hs.append(H)
                    break

    # v.draw_atoms(*axial_Hs, colors="red")
    # v.show()

    return axial_Hs


def get_equatorial_hydrogen_neighbor(molecule, atom):
    """
    Get the equatorial hydrogen atom of a given atom.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The equatorial hydrogen atom (if present, or None if not present)
    """
    ring = next((i for i in molecule._AtomGraph.find_cycles() if atom in i), None)
    if ring is None:
        return None

    ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
    ring_mean = np.array(ring_mean)

    neighbors = molecule.get_neighbors(atom)
    if len(neighbors) < 3:
        return None

    neighbors = list(neighbors)
    vec = atom.get_coord() - ring_mean

    Hs = [i for i in neighbors if i.element == "H"]
    if len(Hs) < 1:
        return None

    for H in Hs:
        H_vec = H.get_coord() - atom.get_coord()
        angle = np.arccos(
            np.dot(vec, H_vec) / (np.linalg.norm(vec) * np.linalg.norm(H_vec))
        )
        angle = np.degrees(angle)
        if angle < 45:
            return H

    return None


def get_axial_hydrogen_neighbor(molecule, atom):
    """
    Get the axial hydrogen atom of a given atom.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The axial hydrogen atom (if present, or None if not present)
    """
    ring = next((i for i in molecule._AtomGraph.find_cycles() if atom in i), None)
    if ring is None:
        return None

    ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
    ring_mean = np.array(ring_mean)

    neighbors = molecule.get_neighbors(atom)
    if len(neighbors) < 3:
        return None

    neighbors = list(neighbors)
    vec = atom.get_coord() - ring_mean

    Hs = [i for i in neighbors if i.element == "H"]
    if len(Hs) < 1:
        return None

    for H in Hs:
        H_vec = H.get_coord() - atom.get_coord()
        angle = np.arccos(
            np.dot(vec, H_vec) / (np.linalg.norm(vec) * np.linalg.norm(H_vec))
        )
        angle = np.degrees(angle)
        if angle > 55:
            return H

    return None


def get_equatorial_neighbor(molecule, atom):
    """
    Get the equatorial atom of a given atom.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The equatorial atom (if present, or None if not present)
    """
    ring = next((i for i in molecule._AtomGraph.find_cycles() if atom in i), None)
    if ring is None:
        return None

    ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
    ring_mean = np.array(ring_mean)

    neighbors = molecule.get_neighbors(atom)
    if len(neighbors) < 3:
        return None

    neighbors = list(neighbors)
    vec = atom.get_coord() - ring_mean

    for neighbor in neighbors:
        if neighbor in ring:
            continue

        neighbor_vec = neighbor.get_coord() - atom.get_coord()
        angle = np.arccos(
            np.dot(vec, neighbor_vec)
            / (np.linalg.norm(vec) * np.linalg.norm(neighbor_vec))
        )
        angle = np.degrees(angle)
        if angle < 45:
            return neighbor

    return None


def get_axial_neighbor(molecule, atom):
    """
    Get the axial atom of a given atom.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The axial atom (if present, or None if not present)
    """
    ring = next((i for i in molecule._AtomGraph.find_cycles() if atom in i), None)
    if ring is None:
        return None

    ring_mean = np.mean([atom.get_coord() for atom in ring], axis=0)
    ring_mean = np.array(ring_mean)

    neighbors = molecule.get_neighbors(atom)
    if len(neighbors) < 3:
        return None

    neighbors = list(neighbors)
    vec = atom.get_coord() - ring_mean

    for neighbor in neighbors:
        if neighbor in ring:
            continue

        neighbor_vec = neighbor.get_coord() - atom.get_coord()
        angle = np.arccos(
            np.dot(vec, neighbor_vec)
            / (np.linalg.norm(vec) * np.linalg.norm(neighbor_vec))
        )
        angle = np.degrees(angle)
        if angle > 55:
            return neighbor

    return None


def get_left_hydrogen_neighbor(molecule, atom):
    """
    Get the "left-protruding" hydrogen neighbor of an atom
    that has two hydrogen neighbors.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The left hydrogen atom.
        If the nomenclature does not apply to the atom (e.g. because it is completely symmetrical on both sides) any of the two hydrogen atoms is returned.
        If the atom does not have two hydrogen neighbors, None may be returned if the only hydrogen neighbor is not in the correct orientation, otherwise the hydrogen neighbor is returned.

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
    Hs, vec, H_vecs = _core_left_right_hydrogens(molecule, atom)

    # get the cross product of the two vectors
    cross = np.cross(vec, H_vecs[0])

    # if the cross product is pointing in the same direction as the second hydrogen, return the first hydrogen
    if np.dot(cross, H_vecs[1]) > 0:
        return Hs[0]
    return Hs[1]


def get_right_hydrogen_neighbor(molecule, atom):
    """
    Get the "right-protruding" hydrogen neighbor of an atom
    that has two hydrogen neighbors.

    Parameters
    ----------
    molecule : Molecule
        The molecule to search in
    atom : Atom
        The atom to search for

    Returns
    -------
    Atom
        The right hydrogen atom.
        If the nomenclature does not apply to the atom (e.g. because it is completely symmetrical on both sides) any of the two hydrogen atoms is returned.
        If the atom does not have two hydrogen neighbors, None may be returned if the only hydrogen neighbor is not in the correct orientation, otherwise the hydrogen neighbor is returned.

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
    Hs, vec, H_vecs = _core_left_right_hydrogens(molecule, atom)

    # get the cross product of the two vectors
    cross = np.cross(vec, H_vecs[0])

    # if the cross product is pointing in the same direction as the second hydrogen, return the second hydrogen
    if np.dot(cross, H_vecs[1]) > 0:
        return Hs[1]
    return Hs[0]


def _core_left_right_hydrogens(molecule, atom):
    """The base function for finding left/right hydrogens"""
    neighbors = molecule.get_neighbors(atom)
    Hs = [i for i in neighbors if i.element == "H"]
    if len(Hs) != 2:
        return None

    non_Hs = [i for i in neighbors if i.element != "H"]
    if len(non_Hs) != 2:
        Hs = (Hs[0], None)

    # sort the non-H neighbors
    a_num_A = atomic_number(non_Hs[0].element)
    a_num_B = atomic_number(non_Hs[1].element)
    if a_num_A > a_num_B:
        non_Hs = (non_Hs[0], non_Hs[1])
    elif a_num_A < a_num_B:
        non_Hs = (non_Hs[1], non_Hs[0])
    else:
        non_Hs = sorted(non_Hs, key=lambda x: _neighbor_sort_key(molecule, x))

    # get the vector between the two non-H neighbors
    vec = non_Hs[0].get_coord() - non_Hs[1].get_coord()

    # get the vector between the central atom and the two hydrogens
    H_vecs = [i.get_coord() - atom.get_coord() for i in Hs]
    return Hs, vec, H_vecs


def _neighbor_sort_key(molecule, atom):
    neighbors = molecule.get_neighbors(atom)
    key = sum(
        atomic_number(i.element) * molecule.get_bond(atom, i).order ** 2
        for i in neighbors
    )
    return key


def vet_structure(
    molecule, clash_range: tuple = (0.6, 2.7), angle_range: tuple = (90, 180)
) -> bool:
    """
    Check for basic structure integrity.
    This will return False if there are clashes within the structure,
    or invalid angles.

    Parameters
    ----------
    molecule : Molecule
        A BuildAMol Molecule
    clash_range : tuple
        The minimal and maximal allowed distances between bonded atoms (in Angstrom).
        The lower limit is also used for non-bonded atoms.
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
    for angle in molecule.compute_angles().values():
        if not angle_range[0] <= angle <= angle_range[1]:
            return False
    for a in molecule.get_atoms():
        for b in molecule.get_atoms():
            if a is b:
                continue
            dist = base.compute_distance(a, b)
            if dist <= clash_range[0]:
                return False
    return True


def find_clashes(molecule, min_dist: float = 1.0, ignore_hydrogens: bool = False):
    """
    Find all clashing atoms in a molecule.

    Parameters
    ----------
    molecule : Molecule
        The molecule to check for clashes.
    min_dist : float
        The minimal allowed distance between atoms (in Angstrom).
    ignore_hydrogens : bool
        If set to True, hydrogen atoms are ignored.

    Yields
    ------
    tuple
        A tuple of clashing atoms.
    """
    warnings.deprecation.warn(
        "find_clashes is deprecated and will be removed in a future version. Use find_clashes_between instead.",
        DeprecationWarning,
    )
    yield from find_clashes_between(
        molecule, molecule, min_dist, ignore_hydrogens, False
    )


def find_clashes_between(
    mol_a,
    mol_b,
    min_dist: float = 1.0,
    ignore_hydrogens: bool = False,
    coarse_precheck: bool = True,
):
    """
    Find all clashing atoms between two sets of atoms.

    Parameters
    ----------
    mol_a, mol_b : Molecule
        Two molecules to check for clashes.
    min_dist : float
        The minimal allowed distance between atoms (in Angstrom).
    ignore_hydrogens : bool
        If set to True, hydrogen atoms are ignored.
    coarse_precheck : bool
        If set to True, a coarse-grain precheck is performed on residue-level before checking for clashes on atom-level.
        This will speed up the process for large molecules but may miss clashes if individual residues are particularly large (e.g. lipids with long tails).

    Yields
    ------
    tuple
        A tuple of clashing atoms.
    """
    # first make a rough distance matrix on residue-level
    residues_a = list(mol_a.get_residues())
    residues_b = list(mol_b.get_residues())
    r_a = np.empty((len(residues_a)), dtype=object)
    r_b = np.empty((len(residues_b)), dtype=object)
    r_a[:] = residues_a
    r_b[:] = residues_b
    residues_a = r_a
    residues_b = r_b

    if coarse_precheck:

        residue_coords_a = np.array([r.center_of_mass() for r in residues_a])
        residue_coords_b = np.array([r.center_of_mass() for r in residues_b])

        residue_dists = cdist(residue_coords_a, residue_coords_b)
        np.fill_diagonal(residue_dists, np.inf)
        residue_edge_mask = np.zeros(residue_dists.shape, dtype=bool)
        for i in range(len(residues_a)):
            for j in range(len(residues_b)):
                if residue_dists[i, j] < 12:
                    residue_edge_mask[i, j] = True
    else:

        residue_edge_mask = np.ones((len(residues_a), len(residues_b)), dtype=bool)

    for i in range(len(residues_a)):
        residue_a = residues_a[i]
        close_by_residues = residues_b[residue_edge_mask[i]]

        if ignore_hydrogens:
            atoms_a = [a for a in residue_a.get_atoms() if a.element != "H"]
            atoms_b = []
            for residue_b in close_by_residues:
                atoms_b.extend([a for a in residue_b.get_atoms() if a.element != "H"])
        else:
            atoms_a = list(residue_a.get_atoms())
            atoms_b = []
            for residue_b in close_by_residues:
                atoms_b.extend(residue_b.get_atoms())

        atoms_a = np.array(atoms_a, dtype=object)
        atoms_b = np.array(atoms_b, dtype=object)
        coords_a = np.array([a.get_coord() for a in atoms_a])
        coords_b = np.array([a.get_coord() for a in atoms_b])
        dists = cdist(coords_a, coords_b)
        np.fill_diagonal(dists, np.inf)

        xs, ys = np.where((0 < dists) * (dists < min_dist))
        for x, y in zip(xs, ys):
            yield atoms_a[x], atoms_b[y]


def sample_atoms_around_reference(
    reference_coord: np.ndarray,
    candidates: np.ndarray,
    num_samples: int,
    max_radius: float = 10.0,
):
    """
    Sample atoms around a reference coordinate. Such that they are spacially evenly distributed around the central reference coordinates.

    Parameters
    ----------
    reference_coord : np.ndarray
        The reference coordinate to sample around.
    candidates : np.ndarray
        The atoms to sample from. This must be an array of Atom objects.
    num_samples : int
        The number of samples to generate.
    max_radius : float
        The maximum radius to sample for.

    Returns
    -------
    samples : np.ndarray
        The sampled atoms.
    """
    # Get the coordinates of the atoms
    coordinates = np.array([i.coord for i in candidates])

    # Generate azimuth and elevation angles evenly spaced in a sphere
    phi = np.linspace(0, 2 * np.pi, num_samples)
    theta = np.linspace(0, np.pi, num_samples)

    # Create a grid of angles
    phi_grid, theta_grid = np.meshgrid(phi, theta)

    # Convert spherical coordinates to Cartesian coordinates
    x = max_radius * np.sin(theta_grid) * np.cos(phi_grid) + reference_coord[0]
    y = max_radius * np.sin(theta_grid) * np.sin(phi_grid) + reference_coord[1]
    z = max_radius * np.cos(theta_grid) + reference_coord[2]

    # Combine the x, y, z coordinates to get the sampled points
    sampled_coordinates = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

    # Find the closest coordinates from the possible_coordinates array
    closest_indices = np.argmin(
        np.linalg.norm(sampled_coordinates[:, None] - coordinates, axis=2),
        axis=1,
    )
    samples = candidates[closest_indices]

    return samples


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

                atoms = list(i for i in residue1.get_atoms() if i.element != "H")
                atoms.extend(i for i in residue2.get_atoms() if i.element != "H")

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


def infer_bond_orders(molecule):
    """
    Infer the bond orders using the registered higher order functional groups (i.e. functional groups with bonds of order > 1).

    Parameters
    ----------
    molecule : Molecule
        The molecule to infer the bond orders for.
    """
    import buildamol.structural.groups as groups

    group_matches = {}
    for atom in molecule.get_atoms():
        neighbors = molecule.get_neighbors(atom)
        atoms = (atom, *neighbors)
        connectivity = [(0, i + 1) for i in range(len(neighbors))]

        for group in groups.higher_order_groups:
            if group.matches(molecule, atoms):
                if atoms in group_matches:
                    if group_matches[atoms][0].rank < group.rank:
                        group_matches[atoms] = (group, group._assignment)
                else:
                    group_matches[atoms] = (group, group._assignment)

    for atoms, (group, assignment) in group_matches.items():
        group._assignment = assignment
        group.apply_connectivity(molecule, atoms)


def _atom_from_residue(id, residue):
    return next((atom for atom in residue.get_atoms() if atom.id == id), None)


def has_free_valence(atom, bonds, needed: int = 1):
    """
    Check if an atom has free valence.

    Parameters
    ----------
    atom : Atom
        The atom to check for free valence.
    bonds : list
        A list of Bond objects that connect the atom to other atoms.
    needed : int
        The number of free valences needed.

    Returns
    -------
    bool
        True if the atom has free valence, False otherwise.
    """
    degree = sum(bond.order for bond in bonds)
    return degree <= element_connectivity.get(atom.element, 0) - needed


def change_element(atom, new_element: str, _molecule):
    """
    Change the element of an atom and update its connectivity accordingly.
    This will remove hydrogen atoms if the new element has a lower connectivity,
    or add them.
    """
    old_element = atom.element
    atom.set_element(new_element, adjust_id=False)
    new_connectivity = element_connectivity.get(new_element, 0)
    reference_connectivity = element_connectivity.get(old_element, 0)

    if _molecule.get_degree(atom) > new_connectivity:
        neighbors = _molecule.get_neighbors(atom)
        to_remove = []
        loops = 0
        while _molecule.get_degree(atom) > new_connectivity:
            for neighbor in neighbors:
                if neighbor.element == "H":
                    to_remove.append(neighbor)
                    neighbors.remove(neighbor)
                    new_connectivity += 1
                    break
            loops += 1
            if loops > 100:
                raise ValueError(
                    f"Could not change element of atom {atom} to {new_element}. Not enough hydrogen atoms available to remove."
                )

        _molecule.remove_atoms(*to_remove)
    elif _molecule.get_degree(atom) < new_connectivity:
        hydrogenator = Hydrogenator()
        hydrogenator.add_hydrogens(atom, _molecule)


# def aromatic_to_double(bonds: list) -> list:
#     """
#     Set bonds between aromatic atoms to single and double bonds
#     This is used to convert aromatic bonds from RDKit which have a bond order of 1.5 to a double bond.

#     Parameters
#     ----------
#     bonds : list
#         A list of Bond objects.

#     Yields
#     -------
#     int
#         The bond order of the bond.
#     """

#     connectivity = defaultdict(int)
#     for bond in bonds:
#         atom1, atom2 = bond
#         connectivity[atom1] += bond.order
#         connectivity[atom2] += bond.order

#     for bond in bonds:
#         atom1, atom2 = bond

#         n_bonds_1 = connectivity[atom1]
#         n_bonds_2 = connectivity[atom2]

#         # aromatic carbons can have 3 or 4 bonds
#         # with 3 bonds they still have room left for one more bond
#         # with 4 bonds they are fully saturated
#         if atom1.element == "C" and atom2.element == "C":

#             if bond.order == int(bond.order):
#                 yield int(bond.order)
#                 continue

#             if n_bonds_1 == 3 or n_bonds_2 == 3:
#                 connectivity[atom1] += 1
#                 connectivity[atom2] += 1
#                 yield 2  # double bond
#                 continue

#             elif n_bonds_1 == 4 or n_bonds_2 == 4:
#                 connectivity[atom1] -= 1
#                 connectivity[atom2] -= 1
#                 yield 1
#                 continue

#             # aromatic bonds from RDKit have 1.5 bond order
#             elif bond.order == 1.5:
#                 yield 1  # single bond
#                 continue

#         # yield any other bond orders as they are
#         if bond.order != int(bond.order):
#             yield bond.order
#         else:
#             yield int(bond.order)


def apply_reference_bonds(structure, _compounds=None):
    """
    Apply bonds according to loaded reference compounds. This will compute a list of tuples with bonded
    atoms from the same residue. It will not infer residue-residue connections! It is possible to provide a single
    atom as input structure, in which case all bonds that involve said atom are returned.

    Parameters
    ----------
    structure
        The structure to apply the bonds to. This can be any of the following objects which hosts Atoms:
        - `buildamol.Structure`
        - `buildamol.Model`
        - `buildamol.Chain`
        - `buildamol.Residue`

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
        bonds = (
            (
                _atom_from_residue(bond.atom1.id, residue),
                _atom_from_residue(bond.atom2.id, residue),
                bond.order,
            )
            for bond in ref.get_bonds()
        )
        bonds = [
            bond for bond in bonds if bond[1] is not None and bond[0] is not None
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

        e_non_H1 = non_H1.element.title()
        e_non_H2 = non_H2.element.title()

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
    import buildamol as bam

    bam.load_sugars()

    mol = bam.Molecule.from_compound("GLC")
    mol.bonds = []
    bonds = apply_reference_bonds(mol.structure)
    mol._add_bonds(*bonds)
    mol.show()

    x = find_clashes(mol)

    # mol = bam.Molecule.from_pubchem("GlcNAc")
    autolabel(mol)
    mol.show()
