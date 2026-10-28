Installation
============

BuildAMol supports Python 3.8+ and is published on PyPI.

Install from PyPI
-----------------

.. code-block:: bash

   pip install buildamol

This installs BuildAMol together with its required dependencies as defined in `pyproject.toml`.

.. note::

   RDKit is now a required dependency of BuildAMol (not an optional extra).


Optional extras
---------------

BuildAMol provides optional extras for visualization and selected integrations:

.. code-block:: bash

   # visualization stack (py3Dmol, nglview, ipywidgets, ...)
   pip install "buildamol[visual]"

   # Open Babel integration
   pip install "buildamol[openbabel]"

   # OpenMM integration
   pip install "buildamol[openmm]"

   # install all optional extras defined by BuildAMol
   pip install "buildamol[full]"

Current extras in this release:

- `visual`
- `openbabel`
- `openmm`
- `full`


Conda/mamba installation (recommended for compiled dependencies)
-----------------------------------------------------------------

If you run into binary/wheel issues for compiled packages on your platform, use conda-forge:

.. code-block:: bash

   mamba create -n buildamol -c conda-forge python=3.11 buildamol
   mamba activate buildamol

Equivalent with conda:

.. code-block:: bash

   conda create -n buildamol -c conda-forge python=3.11 buildamol
   conda activate buildamol


Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/NoahHenrikKleinschmidt/buildamol
   cd buildamol
   pip install .

For development (editable install):

.. code-block:: bash

   pip install -e .


Verify installation
-------------------

.. code-block:: python

   import buildamol as bam
   print(bam.__version__)


Updating BuildAMol
------------------

.. warning::

   Before updating, export any custom default `PDBECompounds` or `CHARMMTopology`
   settings. Updating may overwrite stored pickle files for defaults.

We recommend exporting custom resources in a version-agnostic format (JSON/XML)
before updating:

.. code-block:: python

   import buildamol as bam

   bam.export_custom_resources("my_settings")

   # update BuildAMol ...

   bam.import_custom_resources("my_settings")


Notes on additional external tools
----------------------------------

Some extension backends require third-party tools that are not installed by
default with BuildAMol extras (for example backend-specific MD solvation tools).
Refer to the corresponding extension pages for backend-specific setup steps.

