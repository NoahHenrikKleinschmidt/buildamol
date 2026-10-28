.. _molecular_dynamics:

Molecular Dynamics Extension
============================

The `molecular_dynamics` extension provides a set of tools to help with various tasks related to running molecular dynamics simulations from BuildAMol Molecules.
Currently, it provides the following sub-packages:


.. tab-set::

    .. tab-item:: OpenFF-MD
        
        The `openff_md` package provides an end-to-end pipeline to parametrize a molecule and run a molecular dynamics simulation using the OpenFF toolkit and OpenMM.

        .. dropdown:: OpenFF-MD Settings

            .. automodule:: buildamol.extensions.molecular_dynamics.openff_md
                :members:
                :undoc-members:
                :show-inheritance:

    .. tab-item:: Solvation

        The `solvate` package provides the `solvate` function that can be used to create a solvation box around a molecule.
        It can use the following backend softwares to create the solvation box (requires the respective software to be installed):

        - `biobb_amber <https://github.com/bioexcel/biobb_amber>`_
        
        - `pdbfixer <https://github.com/openmm/pdbfixer>`_

        .. dropdown:: Backend Settings

            .. automodule:: buildamol.extensions.molecular_dynamics.solvate
                :members:
                :undoc-members:
                :show-inheritance:
    
        .. dropdown:: Biobb Amber Backend

            .. automodule:: buildamol.extensions.molecular_dynamics.solvate.backend_biobb
                :members:
                :undoc-members:
                :show-inheritance:

        .. dropdown:: PDBFixer Backend

            .. automodule:: buildamol.extensions.molecular_dynamics.solvate.backend_pdbfixer
                :members:
                :undoc-members:
                :show-inheritance:

    .. tab-item:: Atom Types

        The `atom_typing` package provides a base-class `AtomTyper` class that can be used to assign atom types to a molecule.
        It also implements a `CHARMMTyper` class as an example class. 


        .. dropdown:: Base Atom Typer

            .. automodule:: buildamol.extensions.molecular_dynamics.atom_typing.atom_typer
                :members:
                :undoc-members:
                :show-inheritance:

        .. dropdown:: CHARMM Typer

            .. automodule:: buildamol.extensions.molecular_dynamics.atom_typing.charmm_typer
                :members:
                :undoc-members:
                :show-inheritance:
    

    .. tab-item:: PSF Generation

        The `psf` module works with the :ref:`atom_typing <atom_typing>` package to generate a PSF file for a molecule.
    
        .. automodule:: buildamol.extensions.molecular_dynamics.psf
            :members:
            :undoc-members:
            :show-inheritance: