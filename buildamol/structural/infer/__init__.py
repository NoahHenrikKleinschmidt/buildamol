"""
Functions to infer structural data such as missing atom coordinates,
bond connectivity, or atom labels.
"""

from .analysis import (
    compute_outlier_atoms,
    compute_residue_radius,
    find_clashes_between,
    infer_surface_residues,
    sample_atoms_around_reference,
    vet_structure,
)
from .bond_orders import (
    _infer_bond_orders_native,
    _infer_bond_orders_rdkit,
    infer_bond_orders,
)
from .chemistry import change_element, has_free_valence
from .connectivity import (
    _atom_from_residue,
    apply_reference_bonds,
    atoms_in_area,
    infer_bonds,
    infer_residue_connections,
)
from .constants import (
    _MIN_BOND_LENGTH,
    _VDW_DEFAULT_RADIUS,
    _VDW_PADDING,
    _VDW_SCALE,
    _bond_cutoff_vdw,
    _get_vdw_radius,
    _max_search_radius_vdw,
    acceptable_surplus_charge,
    atomic_number,
    bond_length_by_order,
    double_bond_lengths,
    element_connectivity,
    element_to_hydrogen_bond_lengths,
    element_vdw_radii,
    single_bond_lengths,
    triple_bond_lengths,
)
from .hydrogens import (
    Hydrogenator,
    adjust_protonation,
    adjust_to_ph,
    change_bond_order,
    relabel_hydrogens,
)
from .internal_coordinates import (
    _H_dist_match,
    _H_id_match,
    _prune_H_triplets,
    compute_atom1_from_others,
    compute_atom4_from_others,
    compute_internal_coordinates,
)
from .labeling import AutoLabel, autolabel, autolabel_atoms
from .mapping import infer_mapping_from_template
from .ring import (
    _core_left_right_hydrogens,
    _neighbor_sort_key,
    find_axial_hydrogens,
    find_axial_substituents,
    find_equatorial_hydrogens,
    find_equatorial_substituents,
    get_axial_hydrogen_neighbor,
    get_axial_neighbor,
    get_equatorial_hydrogen_neighbor,
    get_equatorial_neighbor,
    get_left_hydrogen_neighbor,
    get_right_hydrogen_neighbor,
    split_into_contiguous_residues,
)

# Backward-compatible aliases used in tests and user code.
get_left_hydrogen = get_left_hydrogen_neighbor
get_right_hydrogen = get_right_hydrogen_neighbor

# Export public API for autodoc and downstream imports.
__all__ = [name for name in globals() if not name.startswith("_")]
