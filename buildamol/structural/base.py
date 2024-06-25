"""
Basic structure related functions
"""

import numpy as np
import Bio.PDB as bio

import buildamol.utils.auxiliary as aux

origin = np.array([0, 0, 0], dtype=np.float64)
"""
The origin of the 3D coordinate system
"""

x_axis = np.array([1, 0, 0], dtype=np.float64)
"""
Unit vector along the x-axis
"""
y_axis = np.array([0, 1, 0], dtype=np.float64)
"""
Unit vector along the y-axis
"""

z_axis = np.array([0, 0, 1], dtype=np.float64)
"""
Unit vector along the z-axis
"""

xy_plane = np.array([0, 0, 1], dtype=np.float64)
"""
Unit vector normal to the xy-plane
"""

xz_plane = np.array([0, 1, 0], dtype=np.float64)
"""
Unit vector normal to the xz-plane
"""

yz_plane = np.array([1, 0, 0], dtype=np.float64)
"""
Unit vector normal to the yz-plane
"""


def atom_make_full_id(self):
    """
    A self-adjusting full_id for an Biopython Atom
    """
    p = self.get_parent()
    if p:
        return (*p.full_id, (self.id, self.altloc))
    else:
        return (self.id, self.altloc)


def residue_make_full_id(self):
    """
    A self-adjusting full_id for an Biopython Residue
    """
    p = self.get_parent()
    if p:
        return (*p.full_id, self._id)
    else:
        return self._id


def chain_make_full_id(self):
    """
    A self-adjusting full_id for an Biopython Chain
    """
    p = self.get_parent()
    if p:
        return (*p.full_id, self.id)
    return (self.id,)


def model_make_full_id(self):
    """
    A self-adjusting full_id for an Biopython Model
    """
    p = self.get_parent()
    if p:
        return (p.id, self.id)
    return (self.id,)


# --------------------------- POSSIBLE DELETE ---------------------------
# the whole set_full_id and whatever can probably be deleted since
# we are using our own wrapper around the biopython structure anyway
# --------------------------- POSSIBLE DELETE ---------------------------
def set_full_id(self, value):
    pass


bio.Atom.Atom.full_id = property(atom_make_full_id, set_full_id)
bio.Residue.Residue.full_id = property(residue_make_full_id, set_full_id)
bio.Chain.Chain.full_id = property(chain_make_full_id, set_full_id)
bio.Model.Model.full_id = property(model_make_full_id, set_full_id)

# --------------------------- POSSIBLE DELETE ---------------------------


def make_empty_structure(id: str = "empty"):
    """
    Make an empty PDB structure with a single model and chain.

    Returns
    -------
    structure : Bio.PDB.Structure
        The empty structure
    """
    s = bio.Structure.Structure(id)
    m = bio.Model.Model(0)
    s.add(m)
    c = bio.Chain.Chain("A")
    m.add(c)
    return s


def rename_residue(residue: "bio.Residue.Residue", new_name: str):
    """
    Rename a biopython residue object.
    This happens in-place.

    Parameters
    ----------
    residue : Bio.PDB.Residue.Residue
        The residue to rename.
    new_name : str
        The new name of the residue.

    Returns
    -------
    residue : Bio.PDB.Residue.Residue
        The renamed residue.
    """
    residue.resname = new_name
    residue.id = ("H_" + new_name, *residue.id[1:])
    for atom in residue.get_atoms():
        atom.full_id = (*residue.full_id, atom.full_id[-1])
    return residue


def vector_between(atom1, atom2):
    """
    Compute the vector between two atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom

    Returns
    -------
    vector : numpy.ndarray
        The vector between the two atoms
    """
    return atom2.coord - atom1.coord


def norm_vector(atom1, atom2):
    """
    Compute the normalized vector between two atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom

    Returns
    -------
    vector : numpy.ndarray
        The normalized vector between the two atoms
    """
    v = vector_between(atom1, atom2)
    return v / np.linalg.norm(v)


def plane_vector(vec1, vec2):
    """
    Compute the vector of the plane of two vectors.

    Parameters
    ----------
    vec1 : numpy.ndarray
        The first vector
    vec2 : numpy.ndarray
        The second vector

    Returns
    -------
    vector : numpy.ndarray
        The norm-vector of the plane of the two vectors
    """
    v = np.cross(vec1, vec2)
    return v / np.linalg.norm(v)


def bond_vector(bond):
    """
    Compute the vector between two atoms in a bond.

    Parameters
    ----------
    bond : Bio.PDB.Bond
        The bond

    Returns
    -------
    vector : numpy.ndarray
        The vector between the two atoms in the bond
    """
    return vector_between(*bond)


def compute_angle(atom1, atom2, atom3):
    """
    Compute the angle between three atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom

    Returns
    -------
    angle : float
        The angle between the three atoms in degrees
    """
    return angle_between(atom1.coord, atom2.coord, atom3.coord)


def angle_between(coords1, coords2, coords3):
    """
    Compute the angle between three atoms.

    Parameters
    ----------
    coords1 : numpy.ndarray
        The coordinates of the first atom
    coords2 : numpy.ndarray
        The coordinates of the second atom
    coords3 : numpy.ndarray
        The coordinates of the third atom

    Returns
    -------
    angle : float
        The angle between the three atoms in degrees
    """
    a = coords1 - coords2
    b = coords3 - coords2
    return np.degrees(np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))))


def bond_angle(bond1, bond2):
    """
    Compute the angle between two bonds.

    Parameters
    ----------
    bond1 : Bio.PDB.Bond
        The first bond
    bond2 : Bio.PDB.Bond
        The second bond

    Returns
    -------
    angle : float
        The angle between the two bonds in degrees
    """
    if isinstance(bond1[0], np.ndarray):
        return angle_between(bond1[0], bond1[1], bond2[1])
    else:
        return angle_between(bond1[0].coord, bond1[1].coord, bond2[1].coord)


def compute_dihedral(atom1, atom2, atom3, atom4):
    """
    Compute the dihedral angle between four atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom
    atom3 : Bio.PDB.Atom
        The third atom
    atom4 : Bio.PDB.Atom
        The fourth atom

    Returns
    -------
    dihedral : float
        The dihedral angle between the four atoms in degrees
    """
    return dihedral_between(atom1.coord, atom2.coord, atom3.coord, atom4.coord)


def dihedral_between(coords1, coords2, coords3, coords4):
    """
    Compute the dihedral angle between four points

    Parameters
    ----------
    coords1 : numpy.ndarray
        The coordinates of the first atom
    coords2 : numpy.ndarray
        The coordinates of the second atom
    coords3 : numpy.ndarray
        The coordinates of the third atom
    coords4 : numpy.ndarray
        The coordinates of the fourth atom

    Returns
    -------
    dihedral : float
        The dihedral angle between the four atoms in degrees
    """
    ab = coords1 - coords2
    bc = coords3 - coords2
    cd = coords4 - coords3

    # normalize bc so that it does not influence magnitude of vector
    # rejections that come next
    bc /= np.linalg.norm(bc)

    # vector rejections
    v = ab - np.dot(ab, bc) * bc
    w = cd - np.dot(cd, bc) * bc

    # angle between v and w in radians
    x = np.dot(v, w)
    y = np.dot(np.cross(bc, v), w)
    return np.degrees(np.arctan2(y, x))


compute_torsional = compute_dihedral


def distance_between(coord1, coord2):
    """
    Compute the distance between two 3d coordinates.

    Parameters
    ----------
    coord1 : numpy.ndarray
        The first coordinate
    coord2 : numpy.ndarray
        The second coordinate

    Returns
    -------
    distance : float
        The euclidean distance between the two coordinates
    """
    return np.linalg.norm(coord1 - coord2)


def compute_distance(atom1, atom2):
    """
    Compute the distance between two atoms.

    Parameters
    ----------
    atom1 : Bio.PDB.Atom
        The first atom
    atom2 : Bio.PDB.Atom
        The second atom

    Returns
    -------
    distance : float
        The distance between the two atoms
    """
    return np.linalg.norm(atom1.coord - atom2.coord)


def center_of_gravity(masses, coords):
    """
    Compute the center of gravity of a molecule.

    Parameters
    ----------
    masses : array-like
        The masses of the atoms as an nx1 vector
    coords : array-like
        The coordinates of the atoms as an nx3 array

    Returns
    -------
    cog : array-like
        The center of gravity
    """
    return np.average(coords, axis=0, weights=masses)


center_of_mass = center_of_gravity


def center_of_geometry(coords):
    """
    Compute the center of geometry of a molecule.

    Parameters
    ----------
    coords : array-like
        The coordinates of the atoms as an nx3 array

    Returns
    -------
    cog : array-like
        The center of geometry
    """
    return np.mean(coords, axis=0)


def adjust_bond_length(bond, new_length: float):
    """
    Adjust the bond length of a bond.

    Parameters
    ----------
    bond : buildamol.Bond
        The bond to adjust
    new_length : float
        The new bond length

    Returns
    -------
    adjusted_bond : Bio.PDB.Bond
        The adjusted bond
    """
    atom1, atom2 = bond
    vec = atom2.coord - atom1.coord
    vec /= np.linalg.norm(vec)
    atom2.coord = atom1.coord + new_length * vec
    return bond


def adjust_distance(coord1, coord2, new_length: float):
    """
    Adjust the distance between two points.

    Parameters
    ----------
    coord1 : array-like
        The first point
    coord2 : array-like
        The second point
    new_length : float
        The new distance

    Returns
    -------
    new_coord2 : array-like
        The new coordinates of the second point
    """
    vec = coord2 - coord1
    vec /= np.linalg.norm(vec)
    return coord1 + new_length * vec


def rotate_molecule(
    molecule, angle: float, axis: np.ndarray, center: np.ndarray = None
):
    """
    Rotate a molecule around an axis by a given angle.

    Parameters
    ----------
    molecule : Bio.PDB.Structure
        The molecule to rotate
    axis : array-like
        The axis to rotate around
    angle : float
        The angle to rotate by (in radians)

    Returns
    -------
    rotated_molecule : Bio.PDB.Structure
        The rotated molecule
    """
    atoms = list(molecule.get_atoms())
    coords = np.array([a.coord for a in atoms])
    if center is not None:
        coords -= center

    new_coords = rotate_coords(coords=coords, angle=angle, axis=axis)
    if center is not None:
        new_coords += center

    for a, c in zip(atoms, new_coords):
        a.set_coord(c)

    return molecule


def rotate_coords(
    coords: np.ndarray,
    angle: float,
    axis: np.ndarray,
):
    """
    Rotate a set of coordinates around an axis by a given angle.

    Parameters
    ----------
    coords : array-like
        The coordinates to rotate
    angle : float
        The angle to rotate by (in radians)
    axis : array-like
        The axis to rotate around

    Returns
    -------
    rotated_coords : array-like
        The rotated coordinates
    """
    rot = _rotation_matrix(axis, angle)
    return np.dot(np.asarray(coords), rot.T)


def superimpose_points(
    coords: np.ndarray,
    points1: tuple,
    points2: tuple,
):
    """
    Superimpose two structures by aligning two sets of points.
    This will place the first point of points1 on the first point of points2 and so on, while transforming the
    rest of the points accordingly.

    Parameters
    ----------
    coords : array-like
        The coordinates of the atoms to move and align. The coordinates of the participants of Bond1 must be
        part of the coordinates.
    points1 : tuple
        The first set of points. This must be a tuple of two or three coordinates (numpy.ndarray) which are part of the coords.
    points2 : tuple
        The second bond. This must be a tuple of two or three coordinates (numpy.ndarray), depending on how many were provided as points1 (the same number of points).

    Returns
    -------
    new_coords : array-like
        The new coordinates of the atoms
    """

    _old_coords = np.array(points1)
    _new_coords = np.array(points2)

    if len(_old_coords) != len(_new_coords):
        raise ValueError(
            "The number of points in points1 and points2 must be the same."
        )

    # compute translation vector
    old_centroid = _old_coords.mean(axis=0)
    new_centroid = _new_coords.mean(axis=0)

    _relative_old_coords = _old_coords - old_centroid
    _relative_new_coords = _new_coords - new_centroid

    H = (_relative_old_coords).T.dot(_relative_new_coords)
    U, S, VT = np.linalg.svd(H)
    R = VT.T @ U.T

    # Check for reflection
    if np.linalg.det(R) < 0:
        VT[-1, :] *= -1
        R = VT.T @ U.T

    new_coords = (R @ (coords - old_centroid).T).T + new_centroid

    return new_coords

    # get the bond vectors
    v1 = points1[1] - points1[0]
    v2 = points2[1] - points2[0]

    # normalize the bond vectors
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)

    # compute the rotation axis
    axis = np.cross(v1, v2)
    axis /= np.linalg.norm(axis)

    # compute the angle between the bond vectors
    angle = np.arctan(np.dot(v1, v2))

    # rotate the coordinates
    centroid = np.mean([points2[0], points2[1]], axis=0)
    new_coords = coords + centroid
    new_coords = rotate_coords(new_coords, angle, axis)
    new_coords -= centroid

    return new_coords

    # # get the bond vectors
    # bond1_atom1_idx = int(np.where(coords == bond1[0])[0])

    # v0 = bond2[0] - bond1[0]
    # v1 = bond1[1] - bond1[0] - v0
    # v2 = bond2[1] - bond2[0] - v0

    # # normalize the bond vectors
    # v1 /= np.linalg.norm(v1)
    # v2 /= np.linalg.norm(v2)

    # # compute the rotation axis
    # axis = np.cross(v1, v2)
    # axis /= np.linalg.norm(axis)

    # # compute the angle between the bond vectors
    # angle = np.arctan(np.dot(v1, v2))

    # # rotate the coordinates
    # new_coords = coords - v0 - bond1[0]
    # new_coords = rotate_coords(new_coords, angle, axis)
    # new_coords += bond1[0]

    # return new_coords


def plane_of_points(points) -> np.ndarray:
    """
    Compute the plane vector that best describes the plane defined by a set of points

    Parameters
    ----------
    points : array-like
        The points that define the plane

    Returns
    -------
    plane : np.ndarray
        The plane vector
    """
    points = np.array(points)
    centroid = np.mean(points, axis=0)
    points -= centroid
    _, _, v = np.linalg.svd(points)
    return v[2]


def principal_axis(coords: np.ndarray) -> np.ndarray:
    """
    Compute the principle axis of a set of coordinates
    The principle axis is the eigenvector of the covariance matrix with the largest eigenvalue

    Parameters
    ----------
    coords : array-like
        The coordinates of the atoms

    Returns
    -------
    axis : np.ndarray
        The principle axis
    """
    coords = coords - np.mean(coords, axis=0)
    _, _, v = np.linalg.svd(coords)
    return v[0]


def length_along_axis(coords, axis) -> float:
    """
    Compute the length of a set of coordinates along an axis.

    Parameters
    ----------
    coords : array-like
        The coordinates of the atoms
    axis : array-like
        The axis to compute the length along

    Returns
    -------
    length : float
        The length of the coordinates along the axis
    """
    return np.dot(coords, axis).ptp()


@aux.njit
def _numba_wrapper_rotate_coords(
    coords: np.ndarray,
    angle: float,
    axis: np.ndarray,
):
    rot = _numba_wrapper_rotation_matrix(axis, angle)
    return np.dot(np.asarray(coords), rot.T)


def _rotate_coords_base_classes(
    obj,
    angle: float,
    axis: np.ndarray,
    axis_is_absolute: bool = False,
):
    """
    Rotate a set of coordinates around an axis by a given angle.

    Parameters
    ----------
    obj
        The object to rotate
    angle : float
        The angle to rotate by (in radians)
    axis : array-like
        The axis to rotate around
    axis_is_absolute : bool
        Whether the axis is absolute or relative to the object.
        If True, the axis is absolute and the object is rotated around the axis
        which will also incur a translation of the object. If False, the axis
        is relative to the object and the object is rotated around the axis
        without translation.

    Returns
    -------
    rotated_coords : array-like
        The rotated coordinates
    """
    if hasattr(obj, "get_atoms"):
        coords = np.array([a.coord for a in obj.get_atoms()])
        if axis_is_absolute:
            center = np.zeros(3)
        else:
            center = coords.mean()

    else:
        if axis_is_absolute:
            center = np.zeros(3)
        else:
            # since we only have one point
            # and we are just rotating around it
            # there is not going to be any change
            # in the position of the point
            return obj
            center = obj.coord
        coords = obj.coord

    rot = _rotation_matrix(axis, angle)
    new = np.dot(coords - center, rot.T) + center
    if hasattr(obj, "get_atoms"):
        for a, c in zip(obj.get_atoms(), new):
            a.set_coord(c)
    else:
        obj.set_coord(new)
    return obj


bio.Atom.Atom.rotate = _rotate_coords_base_classes
bio.Residue.Residue.rotate = _rotate_coords_base_classes
bio.Chain.Chain.rotate = _rotate_coords_base_classes
bio.Model.Model.rotate = _rotate_coords_base_classes
bio.Structure.Structure.rotate = _rotate_coords_base_classes


def flip_molecule(mol, plane_vector: np.ndarray, center: np.ndarray = None):
    """
    Flip a molecule around an axis.

    Parameters
    ----------
    mol : Molecule
        The molecule to flip
    plane_vector : array-like
        The vector describing the plane to flip around

    Returns
    -------
    flipped_molecule : Bio.PDB.Structure
        The flipped molecule
    """
    atoms = list(mol.get_atoms())
    coords = np.array([a.coord for a in atoms])
    plane_vector = np.array(plane_vector)
    if center is not None:
        center = np.array(center)
        coords -= center

    new_coords = flip_coords(coords=coords, plane_vector=plane_vector)
    if center is not None:
        new_coords += center

    for a, c in zip(atoms, new_coords):
        a.set_coord(c)

    return mol


def flip_coords(coords: np.ndarray, plane_vector: np.ndarray):
    """
    Flip a set of coordinates around an axis.

    Parameters
    ----------
    coords : array-like
        The coordinates to flip
    plane_vector : array-like
        The vector describing the plane to flip around

    Returns
    -------
    flipped_coords : array-like
        The flipped coordinates
    """
    # project the coordinates onto the plane
    # and then subtract the projection from the
    # original coordinates to get the flipped coordinates
    return coords - 2 * np.dot(coords, plane_vector.T)[:, None] * plane_vector


def _rotation_matrix(axis, angle):
    """
    Compute the rotation matrix about an arbitrary axis in 3D

    Source
    ------
    Stackoverflow thread: http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    Parameters
    ----------
    axis : array-like
        The axis to rotate around
    angle : float
        The angle to rotate by (in radians)

    Returns
    -------
    rotation_matrix : array-like
        The rotation matrix
    """
    # axis = np.asarray(axis)
    # angle = np.asarray(angle)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array(
        [
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
        ]
    )


_numba_wrapper_rotation_matrix = aux.njit(_rotation_matrix)


def _euclidean_distances(X, Y):
    """
    Compute the euclidean distances between two sets of points.
    """
    result = np.empty((X.shape[0], Y.shape[0]), dtype=X.dtype)
    for i in range(X.shape[0]):
        for j in range(Y.shape[0]):
            result[i, j] = np.sqrt(np.sum((X[i] - Y[j]) ** 2))
    return result


_numba_wrapper_euclidean_distances = aux.njit(_euclidean_distances)


def _IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    """
    compute the coordinates of a fourth atom from a proper internal coordinate
    system and the other three atom coordinates.

    Parameters
    ----------
    a, b, c : np.ndarray
        coordinates of the other three atoms
    anchor : np.ndarray
        coordinates of the anchor atom relative to which the new coordinate should be calculated
    r : float
        bond length of the new atom relative to the anchor
    theta : float
        bond angle between the new atom and its plane partners
    dihedral : float
        dihedral angle of the internal coordinates
    """
    ab = b - a
    bc = c - b

    # compute normalized bond vectors for available atoms
    ab /= np.linalg.norm(ab)
    bc /= np.linalg.norm(bc)

    # compute plane vector for atoms 1-2-3
    plane_abc = np.cross(ab, bc)
    plane_abc /= np.linalg.norm(plane_abc)

    # rotate the plane vector around the middle bond (2-3) to get the plane 2-3-4
    _rot = _rotation_matrix(bc, dihedral)
    plane_bcd = np.dot(_rot, plane_abc)
    plane_bcd /= np.linalg.norm(plane_bcd)

    # rotate the middle bond around the new plane
    _rot = _rotation_matrix(plane_bcd, theta)
    cd = np.dot(_rot, bc)
    cd /= np.linalg.norm(cd)

    # compute the coordinates of the fourth atom
    d = anchor + r * cd
    return d


@aux.njit
def _numba_wrapper_IC_to_xyz(a, b, c, anchor, r, theta, dihedral):
    ab = b - a
    bc = c - b

    # compute normalized bond vectors for available atoms
    ab /= np.linalg.norm(ab)
    bc /= np.linalg.norm(bc)

    # compute plane vector for atoms 1-2-3
    plane_abc = np.cross(ab, bc)
    plane_abc /= np.linalg.norm(plane_abc)

    # rotate the plane vector around the middle bond (2-3) to get the plane 2-3-4
    _rot = _numba_wrapper_rotation_matrix(bc, dihedral)
    plane_bcd = np.dot(_rot, plane_abc)
    plane_bcd /= np.linalg.norm(plane_bcd)

    # rotate the middle bond around the new plane
    _rot = _numba_wrapper_rotation_matrix(plane_bcd, theta)
    cd = np.dot(_rot, bc)
    cd /= np.linalg.norm(cd)

    # compute the coordinates of the fourth atom
    d = anchor + r * cd
    return d
