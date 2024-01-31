The visual module
=================

Biobuild can generate a variety of visualizations for molecules using multiple backends.
The `visual` module is technically part of the `buildamol.utils` package but is documented
separately here.

Visualizing from the `Molecule` class
-------------------------------------

Each `Molecule` (and `AtomGraph`, and `ResidueGraph`) comes with a number of visualization methods that can be used to generate the below listed visualizations.
So, there is no need to set them up manually, but if the user wishes to customize the visualizations, they can do so by
using the visualization classes directly.


Visualizing with Plotly
-----------------------

Plotly is a web-based visualization library that can be used to generate interactive plots of any kind. 
Biobuild uses Plotly's 3D scatter plots to visualize molecules and their graph representations.

.. dropdown:: PlotlyViewer3D

   .. autoclass:: buildamol.utils.visual.PlotlyViewer3D
      :members:
      :undoc-members:
      :show-inheritance:
   
.. dropdown:: MoleculeViewer3D

   .. autoclass:: buildamol.utils.visual.MoleculeViewer3D
      :members:
      :undoc-members:
      :show-inheritance:


.. dropdown:: AtomGraphViewer3D

   .. autoclass:: buildamol.utils.visual.AtomGraphViewer3D
      :members:
      :undoc-members:
      :show-inheritance:


.. dropdown:: ResidueGraphViewer3D

   .. autoclass:: buildamol.utils.visual.ResidueGraphViewer3D
      :members:
      :undoc-members:
      :show-inheritance:



Visualizing with NglView   
------------------------

NglView is a Python package that can be used to visualize molecules in Jupyter notebooks.
If it is installed, Biobuild can use it to visualize molecules (but not their graph representations).

.. dropdown:: NglViewer3D

   .. autoclass:: buildamol.utils.visual.NglViewer
      :members:
      :undoc-members:
      :show-inheritance:
   

Visualizing with 3DMol.js
-------------------------

3DMol.js is a JavaScript library that can be used to visualize molecules in Jupyter notebooks.
It has a python wrapper `py3Dmol` that Biobuild can use (if installed) to visualize molecules (but not their graph representations).

.. dropdown:: Py3DmolViewer

   .. autoclass:: buildamol.utils.visual.Py3DmolViewer
      :members:
      :undoc-members:
      :show-inheritance:


Visualizing 2D molecules
------------------------

Biobuild can also visualize 2D molecules using the `rdkit` package (if installed).
RDKit is a cheminformatics package that can be used to generate high quality 2D molecule schematics.

.. dropdown:: Chem2DViewer

   .. autoclass:: buildamol.utils.visual.Chem2DViewer
      :members:
      :undoc-members:
      :show-inheritance:

