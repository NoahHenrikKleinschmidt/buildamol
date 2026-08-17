import buildamol.core as core
import buildamol.structural as structural
from buildamol.optimizers import rdkit_optimize
import numpy as np
from typing import Union

from buildamol.extensions.bio.proteins import phi, psi, omega

__all__ = ["alpha_helix", "beta_sheet"]

# Ideal backbone dihedrals (phi, psi) for common secondary structures.
_ALPHA_HELIX = (-57.0, -47.0)
_BETA_SHEET_ANTIPARALLEL = (-139.0, 135.0)
_BETA_SHEET_PARALLEL = (-119.0, 113.0)


def _set_backbone_dihedrals(
    mol: core.Molecule,
    target_phi: float,
    target_psi: float,
    enforce_omega: bool = True,
) -> core.Molecule:
    """
    Set every residue's phi and psi dihedrals to the given target values.

    Residues are processed sequentially from the N-terminus. The current phi,
    psi, and omega are recomputed just before each rotation so that cascading
    changes introduced by earlier residues are accounted for automatically.

    Parameters
    ----------
    mol : Molecule
        The peptide chain. Modified in place.
    target_phi : float
        Target phi angle in degrees.
    target_psi : float
        Target psi angle in degrees.
    enforce_omega : bool
        If True (default), enforce a trans peptide bond (ω = 180°) between
        residue i and i+1 after setting psi_i.
    """
    n_res = mol.count_residues()

    for i in range(1, n_res + 1):
        res = mol.get_residue(i)
        if res is None:
            continue

        CA = res.get_atom("CA")
        N = res.get_atom("N")
        C = res.get_atom("C")
        if None in (CA, N, C):
            continue

        # phi: rotate around N–CA
        _phi = phi(mol, i)
        if not np.isnan(_phi):
            delta = target_phi - _phi
            mol._rotate_around_bond(N, CA, np.deg2rad(delta), descendants_only=True)

        # psi: rotate around CA–C
        _psi = psi(mol, i)
        if not np.isnan(_psi):
            delta = target_psi - _psi
            mol._rotate_around_bond(CA, C, np.deg2rad(delta), descendants_only=True)

        # omega: enforce trans peptide bond at C_i — N_{i+1}
        if enforce_omega:
            _next = mol.get_residue(i + 1)
            if _next is not None:
                N_next = _next.get_atom("N")
                if N_next is not None:
                    _omega = omega(mol, i)
                    if not np.isnan(_omega):
                        delta = 180.0 - _omega
                        mol._rotate_around_bond(
                            C, N_next, np.deg2rad(delta), descendants_only=True
                        )

    return mol


def alpha_helix(
    mol: core.Molecule,
    optimize: bool = False,
) -> core.Molecule:
    """
    Fold a peptide chain into an alpha-helix conformation.

    Sets each residue's backbone dihedrals to the ideal alpha-helical values
    (φ = −57°, ψ = −47°) and enforces trans (ω = 180°) peptide bonds.
    Residues are processed sequentially from the N-terminus; current angles are
    recomputed before every rotation so downstream geometry changes from
    earlier residues propagate correctly.

    Parameters
    ----------
    mol : Molecule
        The peptide chain. Modified in place.
    optimize : bool
        If True, run a short RDKit geometry optimisation after setting the
        dihedrals to relieve any residual local strain. Default False.

    Returns
    -------
    Molecule
        The modified molecule (same object as *mol*).

    Examples
    --------
    >>> from buildamol.extensions.bio.proteins import peptide
    >>> from buildamol.extensions.bio.proteins.secondary_structure import alpha_helix
    >>> pep = peptide("AAAAAAAAAAAA")
    >>> helix = alpha_helix(pep)
    """
    _set_backbone_dihedrals(mol, *_ALPHA_HELIX)
    if optimize:
        mol = rdkit_optimize(mol, steps=500)
    return mol


def beta_sheet(
    mol: core.Molecule,
    parallel: bool = False,
    optimize: bool = False,
) -> core.Molecule:
    """
    Fold a peptide chain into a beta-sheet conformation.

    Sets each residue's backbone dihedrals to ideal beta-sheet values and
    enforces trans (ω = 180°) peptide bonds.  Two geometries are supported:

    * **Antiparallel** (default): φ = −139°, ψ = +135°
    * **Parallel**: φ = −119°, ψ = +113°

    Note that this places a single strand into the extended beta conformation.
    To build a full beta sheet, create multiple strands and arrange them
    spatially with the appropriate hydrogen-bonding geometry.

    Parameters
    ----------
    mol : Molecule
        The peptide chain. Modified in place.
    parallel : bool
        Use parallel beta-sheet dihedrals instead of antiparallel.
        Default False (antiparallel).
    optimize : bool
        If True, run a short RDKit geometry optimisation after setting the
        dihedrals. Default False.

    Returns
    -------
    Molecule
        The modified molecule (same object as *mol*).

    Examples
    --------
    >>> from buildamol.extensions.bio.proteins import peptide
    >>> from buildamol.extensions.bio.proteins.secondary_structure import beta_sheet
    >>> strand = peptide("VVVVVVVV")
    >>> sheet_strand = beta_sheet(strand)
    """
    target = _BETA_SHEET_PARALLEL if parallel else _BETA_SHEET_ANTIPARALLEL
    _set_backbone_dihedrals(mol, *target)
    if optimize:
        mol = rdkit_optimize(mol, steps=500)
    return mol


if __name__ == "__main__":
    from buildamol.extensions.bio.proteins import peptide

    pep = peptide("AAAAAAAAAAAA")

    helix = alpha_helix(pep.copy())
    helix.to_pdb("helix.pdb")

    strand = beta_sheet(pep.copy())
    strand.to_pdb("beta_sheet.pdb")
