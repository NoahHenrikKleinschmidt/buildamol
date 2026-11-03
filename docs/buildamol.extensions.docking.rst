.. _docking:

Molecular Docking Extension
============================

The `docking` extension provides a set of tools to help with docking drug-like molecules to proteins.

The main function is the `dock` function that can be used to perform the docking using several backend softwares (requires the respective software to be installed):

Available backends are:

- `vina <https://github.com/ccsb-scripps/AutoDock-Vina>`_ (default) - uses the AutoDock Vina Python bindings directly. This is the fastest backend and
  is designed for generative pipelines and optimization loops: use a `Docker` to prepare the receptor and affinity maps once (maps are only re-computed when a ligand introduces a new Vina atom type),
  locally refine (``mode="optimize"``) or score (``mode="score"``) candidates that were built inside the pocket in milliseconds, and dock batches of candidates in parallel.
  Optimizing and scoring require the whole ligand to lie inside the box; use ``box_around(ligand)`` to get a box that fits.

- `easydock <https://github.com/ci-lab-cz/easydock>`_

- `dockstring <https://github.com/dockstring/dockstring>`_

- `opendock <https://github.com/guyuehuo/opendock-open>`_

- `diffdock <https://github.com/gcorso/DiffDock>`_

The backend can be changed using `set_docking_backend`.


.. dropdown:: Backend Settings

    .. automodule:: buildamol.extensions.docking
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: Vina Backend

    .. automodule:: buildamol.extensions.docking.backend_vina
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: EasyDock Backend

    .. automodule:: buildamol.extensions.docking.backend_easydock
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: DockString Backend

    .. automodule:: buildamol.extensions.docking.backend_dockstring
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: Opendock Backend

    .. automodule:: buildamol.extensions.docking.backend_opendock
        :members:
        :undoc-members:
        :show-inheritance:

.. dropdown:: Diffdock Backend

    .. automodule:: buildamol.extensions.docking.backend_diffdock
        :members:
        :undoc-members:
        :show-inheritance:
