"""
Ring-orientation related inference helpers (equatorial/axial and left/right hydrogens).
"""

import numpy as np
import networkx as nx

import buildamol.base_classes as base_classes

from .constants import atomic_number


def find_equatorial_substituents(molecule):
    """
    Find all atoms in a molecule that are equatorial to some ring within the molecule.
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
    """
    rings = molecule._AtomGraph.find_cycles()
    equatorial_Hs = []

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

    return equatorial_Hs


def find_axial_hydrogens(molecule):
    """
    Find the axial hydrogen atoms in a molecule.
    """
    rings = molecule._AtomGraph.find_cycles()
    axial_Hs = []

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

    return axial_Hs


def get_equatorial_hydrogen_neighbor(molecule, atom):
    """
    Get the equatorial hydrogen atom of a given atom.
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
    Get the "left-protruding" hydrogen neighbor of an atom that has two hydrogen neighbors.
    """
    out = _core_left_right_hydrogens(molecule, atom)
    if out is None:
        return None
    Hs, vec, H_vecs = out

    cross = np.cross(vec, H_vecs[0])

    if np.dot(cross, H_vecs[1]) > 0:
        return Hs[0]
    return Hs[1]


def split_into_contiguous_residues(molecule, target_residues: list = None):
    """
    Split disconnected atom groups into separate residues.
    """
    if not molecule.count_bonds():
        raise ValueError(
            "The molecule has no bonds. Cannot split into contiguous residues due to missing connectivity information. First add bonds, for example by running `molecule.infer_bonds()`."
        )

    atom_graph = molecule.get_atom_graph()

    if target_residues is not None:
        target_residues = molecule.get_residues(target_residues)
        other_residues = [
            res for res in molecule.get_residues() if res not in target_residues
        ]
        atom_nodes_to_drop = set()
        for res in other_residues:
            atom_nodes_to_drop.update(res.child_list)
        atom_graph.remove_nodes_from(atom_nodes_to_drop)

    contiguous_subgraphs = nx.connected_components(atom_graph)
    old_residues_to_remove = set()
    for i, subgraph in enumerate(contiguous_subgraphs):
        new_residue = base_classes.Residue(resname=f"UNL_{i+1}")
        molecule.add_residues(new_residue)
        for atom in subgraph:
            current_residue = atom.parent
            current_residue.detach_child(atom.get_id())
            new_residue.add(atom)
            old_residues_to_remove.add(current_residue)

    molecule.remove_residues(old_residues_to_remove)
    return molecule


def get_right_hydrogen_neighbor(molecule, atom):
    """
    Get the "right-protruding" hydrogen neighbor of an atom that has two hydrogen neighbors.
    """
    out = _core_left_right_hydrogens(molecule, atom)
    if out is None:
        return None
    Hs, vec, H_vecs = out

    cross = np.cross(vec, H_vecs[0])

    if np.dot(cross, H_vecs[1]) > 0:
        return Hs[1]
    return Hs[0]


def _core_left_right_hydrogens(molecule, atom):
    """The base function for finding left/right hydrogens."""
    neighbors = molecule.get_neighbors(atom)
    Hs = [i for i in neighbors if i.element == "H"]
    if len(Hs) != 2:
        return None

    non_Hs = [i for i in neighbors if i.element != "H"]
    if len(non_Hs) != 2:
        Hs = (Hs[0], None)

    a_num_A = atomic_number(non_Hs[0].element)
    a_num_B = atomic_number(non_Hs[1].element)
    if a_num_A > a_num_B:
        non_Hs = (non_Hs[0], non_Hs[1])
    elif a_num_A < a_num_B:
        non_Hs = (non_Hs[1], non_Hs[0])
    else:
        non_Hs = sorted(non_Hs, key=lambda x: _neighbor_sort_key(molecule, x))

    vec = non_Hs[0].get_coord() - non_Hs[1].get_coord()

    H_vecs = [i.get_coord() - atom.get_coord() for i in Hs]
    return Hs, vec, H_vecs


def _neighbor_sort_key(molecule, atom):
    neighbors = molecule.get_neighbors(atom)
    key = sum(
        atomic_number(i.element) * molecule.get_bond(atom, i).order ** 2
        for i in neighbors
    )
    return key
