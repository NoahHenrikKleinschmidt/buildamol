.. _docking:

Molecular Docking Extension
============================

The `docking` extension provides a set of tools to help with docking drug-like molecules to proteins.

The main function is the `dock` function that can be used to perform the docking using several backend softwares (requires the respective software to be installed):

Available backends are:

- `easydock <https://github.com/ci-lab-cz/easydock>`_

- `dockstring <https://github.com/dockstring/dockstring>`_


.. dropdown:: Backend Settings

    .. automodule:: buildamol.extensions.docking
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: EasyDock Backend

    .. automodule:: buildamol.extensions.docking.easydock
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: DockString Backend

    .. automodule:: buildamol.extensions.docking.dockstring
        :members:
        :undoc-members:
        :show-inheritance:
