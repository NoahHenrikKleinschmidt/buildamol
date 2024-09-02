import buildamol.core as core
import buildamol.structural as structural
from buildamol.optimizers import rdkit_optimize
import numpy as np
from typing import Union

from buildamol.extensions.bio.proteins import phi, psi


# THIS DOES NOT WORK YET
# (I mean, it does but the conformation overall does not look good yet...)
def alpha_helix(mol: core.Molecule) -> core.Molecule:
    """
    Turn a peptide chain into an alpha helix

    Parameters
    ----------
    mol : Molecule
        The peptide chain

    Returns
    -------
    Molecule
        The alpha helix. The input molecule is modified in place.
    """
    # get the phi and psi angles
    phis = phi(mol)
    psis = psi(mol)

    for i in range(len(phis)):
        _phi = phis[i]
        if np.isnan(_phi):
            continue

        # rotate phi angle to -60
        res = mol.get_residue(i + 1)
        CA = res.get_atom("CA")
        N = res.get_atom("N")
        mol._rotate_around_bond(
            N, CA, np.pi - np.deg2rad(_phi) + np.deg2rad(-60), descendants_only=True
        )

        # rotate psi angle to -45
        _psi = psis[i]
        if np.isnan(_psi):
            continue
        C = res.get_atom("C")
        mol._rotate_around_bond(
            CA, C, np.pi - np.deg2rad(_psi) + np.deg2rad(-45), descendants_only=True
        )

    # energy minimize the structure
    mol = rdkit_optimize(mol, steps=5000)

    return mol


if __name__ == "__main__":

    # pep = core.read_pdb("/Users/noahhk/val.pdb")
    from buildamol.extensions.bio.proteins import peptide

    pep = peptide("MLQSMVSLLQSLVSLIIQ")

    helix = alpha_helix(pep.copy())
    helix.to_pdb("helix.pdb")
    # v = pep.draw(atoms=False)
    # v += helix.draw(atoms=False, line_color="red")
    # helix2 = helix.copy().optimize(algorithm="rdkit")
    # helix2.to_pdb("helix2.pdb")
    # v += helix2.draw(atoms=False, line_color="blue")

    # v.show()
