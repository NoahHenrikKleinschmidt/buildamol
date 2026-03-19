"""
Template-based atom mapping and bond transfer inference.
"""

from collections import deque
import warnings

import numpy as np
from scipy.spatial.distance import cdist


def infer_mapping_from_template(
    target,
    template,
    anchors: dict,
    strict: bool = True,
    distance_tolerance: float = 0.3,
):
    """
    Match atoms of a target molecule to a template and transfer bond topology.
    """
    if not isinstance(anchors, dict) or len(anchors) < 3:
        raise ValueError(
            "Anchors must be a dictionary of atom mappings between the target (key) and template (value) molecules with at least 3 entries."
        )

    def _resolve_anchor(atom_id, molecule, which):
        atom = molecule.get_atom(atom_id)
        if atom is None:
            raise ValueError(
                f"{which} atom '{atom_id}' not found in {which.lower()} molecule."
            )
        return atom

    target_template_pairs = []
    for target_id, template_id in anchors.items():
        target_atom = _resolve_anchor(target_id, target, "Target")
        template_atom = _resolve_anchor(template_id, template, "Template")
        target_template_pairs.append((target_atom, template_atom))

    if strict:
        if target.count_atoms() != template.count_atoms():
            raise ValueError(
                "Target and template molecules must have the same number of atoms in strict mode."
            )
        element_hist = {}
        for atom in target.get_atoms():
            element_hist[atom.element] = element_hist.get(atom.element, 0) + 1
        for atom in template.get_atoms():
            element_hist[atom.element] = element_hist.get(atom.element, 0) - 1
        mismatching = [element for element, count in element_hist.items() if count != 0]
        if mismatching:
            raise ValueError(
                "Element counts do not match in strict mode for: "
                + ", ".join(sorted(mismatching))
            )

    template_atoms = list(template.get_atoms())
    target_atoms = list(target.get_atoms())
    template_indices = {atom: idx for idx, atom in enumerate(template_atoms)}
    target_indices = {atom: idx for idx, atom in enumerate(target_atoms)}

    template_coords = np.array([atom.get_coord() for atom in template_atoms])
    target_coords = np.array([atom.get_coord() for atom in target_atoms])
    template_pairwise_distances = cdist(template_coords, template_coords)
    target_pairwise_distances = cdist(target_coords, target_coords)

    template_to_target = {tpl: tgt for tgt, tpl in target_template_pairs}
    target_to_template = {tgt: tpl for tgt, tpl in target_template_pairs}
    mapped_template_indices = [template_indices[tpl] for tpl in template_to_target]
    mapped_target_indices = [
        target_indices[template_to_target[tpl]] for tpl in template_to_target
    ]

    min_anchor_matches = 2

    unmapped_targets = set(target_atoms) - set(target_to_template.keys())

    queue = deque()
    queued = set()

    def _enqueue_template_neighbors(template_atom):
        for neighbor in template_atom.get_neighbors():
            if neighbor in template_to_target or neighbor in queued:
                continue
            queue.append(neighbor)
            queued.add(neighbor)

    def _match_target(template_atom):
        template_idx = template_indices[template_atom]
        template_distances = template_pairwise_distances[template_idx][
            mapped_template_indices
        ]
        best_candidate = None
        best_score = (-1, np.inf)
        required_matches = min(min_anchor_matches, len(mapped_template_indices))

        for target_atom in list(unmapped_targets):
            if target_atom.element != template_atom.element:
                continue
            target_idx = target_indices[target_atom]
            target_distances = target_pairwise_distances[target_idx][
                mapped_target_indices
            ]
            diff = np.abs(target_distances - template_distances)
            matches = diff <= distance_tolerance
            match_count = int(matches.sum())
            if match_count < required_matches:
                continue
            score = diff[matches].sum()
            if match_count > best_score[0] or (
                match_count == best_score[0] and score < best_score[1]
            ):
                best_candidate = target_atom
                best_score = (match_count, score)
        return best_candidate

    for _, template_anchor in target_template_pairs:
        _enqueue_template_neighbors(template_anchor)

    template_bond_orders = [(bond, bond.order) for bond in template.get_bonds()]
    try:
        for bond in template.get_bonds():
            bond.order = 1

        while queue:
            progress_made = False
            level_size = len(queue)
            for _ in range(level_size):
                template_atom = queue.popleft()
                queued.discard(template_atom)
                candidate_target = _match_target(template_atom)
                if candidate_target is None:
                    queue.append(template_atom)
                    queued.add(template_atom)
                    continue

                progress_made = True
                template_to_target[template_atom] = candidate_target
                target_to_template[candidate_target] = template_atom
                unmapped_targets.discard(candidate_target)
                mapped_template_indices.append(template_indices[template_atom])
                mapped_target_indices.append(target_indices[candidate_target])

                _enqueue_template_neighbors(template_atom)

            if queue and not progress_made:
                unresolved = sorted(atom.id for atom in queue)
                msg = "Could not resolve mapping for template atoms: " + ", ".join(
                    unresolved
                )
                if strict:
                    raise ValueError(msg)
                else:
                    warnings.warn(msg)
                    break

        if len(template_to_target) != len(template_atoms):
            missing = [
                atom.id for atom in template_atoms if atom not in template_to_target
            ]
            if strict:
                raise ValueError(
                    "Failed to map all template atoms. Missing: "
                    + ", ".join(sorted(missing))
                )
            else:
                warnings.warn(
                    "Failed to map all template atoms. Missing: "
                    + ", ".join(sorted(missing))
                )

        mapped_bonds = [
            (template_to_target[bond.atom1], template_to_target[bond.atom2], order)
            for bond, order in template_bond_orders
            if bond.atom1 in template_to_target and bond.atom2 in template_to_target
        ]

    finally:
        for bond, order in template_bond_orders:
            bond.order = order

    target_to_template = {v: k for k, v in template_to_target.items()}
    return target_to_template, mapped_bonds
