"""
Hydrogen-related and protonation-related inference utilities.
"""

from typing import Union

import numpy as np
from scipy.spatial.distance import cdist

import buildamol.base_classes as base_classes
import buildamol.structural.base as base
import buildamol.structural.geometry as geometry
import buildamol.utils.auxiliary as aux

from .constants import (
    bond_length_by_order,
    element_connectivity,
    element_to_hydrogen_bond_lengths,
)
from .labeling import AutoLabel


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
        (2, 5): tetrahedral,
    }

    _bond_length = 1.05

    def infer_hydrogens(self, molecule, bond_length: float = 1):
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
        _geometry = Hydrogenator.__geometries__.get((bond_order, connectivity), None)
        if _geometry is None:
            return

        _neighbors = list(neighbors)[: _geometry.max_points - 1]
        out = _geometry.make_coords(atom, *_neighbors, length=self._bond_length)[
            len(_neighbors) :
        ]

        labels = AutoLabel.hydrogen_neighbors(atom)
        if free_slots > 1:
            labels.pop(0)

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
            hydrogens = molecule.get_hydrogens(atom)
            molecule.remove_atoms(*hydrogens)
            atom.pqr_charge = new_charge
            H.add_hydrogens(atom, molecule)
        else:
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


def adjust_to_ph(
    molecule, ph: Union[float, int, tuple], inplace: bool = True, **kwargs
):
    """
    Adjust the protonation state of a molecule to suit a specific pH range.
    """
    if not (
        aux.has_package("scrubber") or aux.has_package("molscrub")
    ) or not aux.has_package("rdkit"):
        raise ImportError(
            "The `molscrub` and `rdkit` packages are required for this function."
        )
    if not inplace:
        molecule = molecule.copy()
    mol = molecule.remove_hydrogens().to_rdkit()

    if isinstance(ph, (float, int)):
        pH_low = ph
        pH_high = None
    elif isinstance(ph, (tuple, list, np.ndarray)):
        pH_low, pH_high = ph
    else:
        raise ValueError(
            f"Invalid pH value or range provided. Expected float, int, or tuple, got {type(ph)}."
        )

    if aux.has_package("scrubber"):
        import importlib

        Scrub = importlib.import_module("scrubber").Scrub
    elif aux.has_package("molscrub"):
        import importlib

        Scrub = importlib.import_module("molscrub").Scrub
    else:
        raise ImportError("The `molscrub` package not found in the environment.")

    skip_tautomers = kwargs.pop("skip_tautomers", True)
    scrub = Scrub(
        ph_high=pH_high, ph_low=pH_low, skip_tautomers=skip_tautomers, **kwargs
    )
    out = scrub(mol)
    if len(out) == 0:
        raise ValueError(
            "Could not adjust the protonation state of the molecule. Try a different pH value or range or manually run the molecule through `molscrub` to investigate potential error sources."
        )

    if len(out) > 1:
        _out = []
        for i in out:
            _out.append(molecule.__class__.from_rdkit(i))
        out = _out
        for i in out:
            i.id = molecule.id
    else:
        out = molecule.__class__.from_rdkit(out[0])
        molecule._base_struct = out._base_struct
        molecule._model = out._model
        molecule._bonds = out._bonds
        molecule._AtomGraph = out._AtomGraph
        out = molecule
    return out


def relabel_hydrogens(molecule):
    """
    Relabel hydrogen atoms in a structure to match the CHARMM naming scheme.
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
    Change the bond order between two atoms.
    """
    if order not in (1, 2, 3):
        raise ValueError(f"Invalid bond order {order}.")

    bond = molecule.get_bond(atom1, atom2)
    if bond is None:
        raise ValueError(f"No bond found between {atom1} and {atom2}.")

    if bond.order == order:
        return

    hydrogens = []
    for atom in (atom1, atom2):
        for hydrogen in molecule.get_neighbors(atom):
            if hydrogen.element == "H":
                hydrogens.append(hydrogen)
    molecule.remove_atoms(*hydrogens)

    bond.order = order

    length = bond_length_by_order[order][atom1.element][atom2.element]
    base.adjust_bond_length(bond, length)

    H = Hydrogenator()
    H.add_hydrogens(atom1, molecule)
    H.add_hydrogens(atom2, molecule)

    return molecule
