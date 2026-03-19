"""
Atom labeling inference utilities.
"""

from collections import defaultdict

import networkx as nx
import pandas as pd
import periodictable as pt

from .constants import element_connectivity


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
        carbons = carbons.reset_index(drop=True)
        return carbons.atom[0]

    def _parse_c_next(self, c_current):
        neighbors = self.graph.neighbors(c_current)
        neighbors = self._df[self._df.atom.isin(neighbors)]
        neighbors = neighbors[neighbors.element == "C"]
        neighbors = neighbors[neighbors.label == "none"]
        neighbors = neighbors.sort_values("total", ascending=False)
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
                for h in _heteros.values():
                    for idx, atom in enumerate(h):
                        self._df.loc[self._df.atom == h, "label"] += str(idx + 1)

    def _parse_hydrogen_labels(self):
        _neighbor_connect_dict = {}
        _hydrogen_unique_label_dict = {}
        hydrogens = self._df[self._df.element == "H"]
        for hydrogen in hydrogens.atom:
            neighbors = self.graph.neighbors(hydrogen)
            neighbors = self._df[self._df.atom.isin(neighbors)]
            neighbors = neighbors[neighbors.label != "none"]
            neighbors = neighbors.sort_values("label", ascending=False)
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
    for idx in range(len(df)):
        atom = df.atom.iloc[idx]
        label = df.label.iloc[idx]
        atom.id = label
        atom.name = label

    return molecule


def autolabel_atoms(bonds, to_label):
    """
    Automatically relabel a list of atoms to match the CHARMM naming scheme.
    Note, this function is not guaranteed to produce the correct labels in all cases,
    validation of the labels is recommended.

    Parameters
    ----------
    bonds : list
        A list of Bond objects that define the connectivity of atoms relevant for the labelling.
    to_label : list
        A list of Atom objects to relabel. This must be a subset of the atoms present in the bonds.

    Returns
    -------
    dict
        A dictionary mapping the original Atom objects to their new labels.
    """
    from buildamol.graphs import AtomGraph

    g = AtomGraph("autolabel_temp_graph", bonds=bonds)

    labeler = AutoLabel(g)
    df = labeler.autolabel()
    label_dict = {}
    for idx in range(len(df)):
        atom = df.atom.iloc[idx]
        label = df.label.iloc[idx]
        if atom in to_label:
            label_dict[atom] = label
    return label_dict
