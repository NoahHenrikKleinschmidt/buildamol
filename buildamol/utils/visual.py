"""
Visualization auxiliary functions
"""

import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.colors as colors

import buildamol.utils as utils
import buildamol.utils.auxiliary as aux

import periodictable

Draw = aux.Draw
Chem = aux.Chem

DEFAULT_BACKEND = "plotly"
"""
The default visualization backend for bare `draw` and `show` methods on objects. 
"""


def set_backend(backend: str):
    """
    Set the default visualization backend, which will be used by objects when calling `draw` and `show` methods.

    Parameters
    ----------
    backend : str
        Available backends are:
            - plotly (default)
            - py3dmol
            - nglview
    """
    backend = backend.strip().lower()
    if backend not in ("plotly", "py3dmol", "nglview"):
        raise ValueError(f"Unsupported backend: {backend}")
    global DEFAULT_BACKEND
    DEFAULT_BACKEND = backend


default_plotly_opacity = 1.0
"""
The default opacity for plotly-based visualizations.
"""

default_plotly_marker_size = 5
"""
The default marker size for plotly-based visualizations.
"""

default_plotly_bond_color = "black"
"""
The default color for plotly-based bond visualizations.
"""

default_plotly_linewidth = 1.2
"""
The default linewidth for plotly-based bond visualizations.
"""


class Chem2DViewer:
    """
    View a molecule in 2D using the RDKit library.

    Parameters
    ----------
    molecule
        The molecule to view. This may be any object that holds
        a biopython structure e.g. a Molecule, AtomGraph, or ResidueGraph.
    drawer: str
        The 2D drawer to use. This can be any of:
            - png (default, uses `MolDraw2DCairo`, requires cairo to be installed)
            - svg (uses `MolDraw2DSVG`)
    highlight_color : str
        The color to use for highlighting atoms and bonds (deprecated, specify a color when calling `highlight_atoms` or `highlight_bonds` instead).
    linewidth : float
            The linewidth of the bonds (deprecated, specify a linewidth when calling the `draw` method instead).
    atoms : str
        The label to use for the atoms (deprecated, specify a label when calling the `label_atoms` method instead).
        This can be any of the following:
        - None (element, except for carbons)
        - "element" (elements, even for carbons)
        - "serial" (the atom serial number)
        - "id" (the atom id / name)
        - "resid" (atom id + parent residue)
        - "off" (no label)
        - any function that takes an (rdkit) atom and returns a string.
    """

    def __init__(
        self,
        molecule,
        drawer: str = "png",
        highlight_color: str = None,
        linewidth: float = None,
        atoms: str = None,
    ):
        if Chem is None:
            raise ImportError(
                "rdkit is not available. Please install it and be sure to use a compatible environment."
            )

        drawer = drawer.strip().lower()
        if drawer not in ("svg", "png"):
            raise ValueError(f"Unsupported drawer: {drawer}")

        self._drawer_type = drawer

        self._raw_molecule = None
        self._raw_is_rdkit = False
        if hasattr(molecule, "to_rdkit"):
            mol = molecule.to_rdkit()
            self._raw_molecule = molecule
        elif molecule.__class__.__name__ in ("AtomGraph", "ResidueGraph"):
            mol = molecule._molecule.to_rdkit()
            self._raw_molecule = molecule
        elif "Chem" in str(molecule.__class__.mro()[0]):
            mol = molecule
            self._raw_molecule = molecule
            self._raw_is_rdkit = True
        else:
            raise ValueError(
                f"Unsupported molecule type: {molecule.__class__.__name__}"
            )
        mol.RemoveAllConformers()
        self.mol = mol

        if atoms is not None:
            aux.deprecation_warning(
                "The `atoms` argument is deprecated and will be removed in future versions. Please use the `label_atoms` method instead."
            )
            if atoms == "element":
                atoms = lambda atom: atom.GetSymbol()
            elif atoms == "serial":
                atoms = lambda atom: str(atom.GetPDBResidueInfo().GetSerialNumber())
            elif atoms == "id":
                atoms = lambda atom: atom.GetPDBResidueInfo().GetName().strip()
            elif atoms == "off":
                atoms = lambda atom: ""
            elif atoms == "resid":

                def atoms(atom):
                    info = atom.GetPDBResidueInfo()
                    return f"{info.GetName().strip()}@{info.GetResidueName().strip()}[{info.GetResidueNumber()}]"

            elif callable(atoms):
                pass
            else:
                raise ValueError(f"Unsupported atom label: {atoms}")

            for a in mol.GetAtoms():
                a.SetProp("atomLabel", atoms(a))

        self._atoms_to_highlight = []
        self._atoms_highlight_colors = {}
        self._atoms_highlight_radii = {}
        self._bonds_to_highlight = {}

        if highlight_color is not None:
            aux.deprecation_warning(
                "The `highlight_color` argument is deprecated and will be removed in future versions. Please specify a color when calling `highlight_atoms` or `highlight_bonds` instead."
            )
        if linewidth is not None:
            aux.deprecation_warning(
                "The `linewidth` argument is deprecated and will be removed in future versions. Please specify a linewidth when calling the `draw` method instead."
            )
        self.highlight_color = highlight_color or "cyan"
        self.linewidth = linewidth or 1
        self.options = Draw.MolDrawOptions()
        self._custom_colors = {}

    def label_atoms(self, func_or_mapping, rdkit: bool = None):
        """
        Generate custom atom labels

        Parameters
        ----------
        func_or_mapping : str, callable or dict
            If a string is provided it has to be one of the following:
            - "element" (elements, even for carbons)
            - "serial" (the atom serial number)
            - "id" (the atom id / name)
            - "resid" (atom id + parent residue)
            - "off" (no label)
            Alternatively, either a function that takes an atom and returns a string. Or a dictionary mapping atoms to strings.
            Only one type of key can be included in the dictionary!
            Supported dictionary keys are:
            - BuildAMol Atoms
            - atom serial numbers (int)
            - atom ids (str) (will match all atoms with that id)
        rdkit : bool
            Whether the function takes an RDKit atom or a BuildAMol atom.
        """
        if isinstance(func_or_mapping, str):
            if func_or_mapping == "element":
                func_or_mapping = lambda atom: atom.GetSymbol()
            elif func_or_mapping == "serial":
                func_or_mapping = lambda atom: str(
                    atom.GetPDBResidueInfo().GetSerialNumber()
                )
            elif func_or_mapping == "id":
                func_or_mapping = (
                    lambda atom: atom.GetPDBResidueInfo().GetName().strip()
                )
            elif func_or_mapping == "off":
                func_or_mapping = lambda atom: ""
            elif func_or_mapping == "resid":
                func_or_mapping = (
                    lambda atom: f"{atom.GetPDBResidueInfo().GetName().strip()}@{atom.GetPDBResidueInfo().GetResidueName().strip()}[{atom.GetPDBResidueInfo().GetResidueNumber()}]"
                )
            else:
                raise ValueError(
                    f"Unsupported atom label for string identifier: '{func_or_mapping}'. Supported are 'element', 'serial', 'id', 'resid', and 'off'."
                )
            rdkit = True

        if rdkit is False and self._raw_is_rdkit:
            raise ValueError(
                "The underlying molecule is an RDKit molecule. Cannot perform BuildAMol operations on RDKit Atoms. Please set `rdkit=True`."
            )

        if callable(func_or_mapping):
            try:

                test_atom = next(iter(self.mol.GetAtoms()))
                func_or_mapping(test_atom)
                rdkit = True
            except Exception:
                if self._raw_is_rdkit:
                    raise ValueError(
                        "The underlying molecule is an RDKit molecule. Cannot perform BuildAMol operations on RDKit Atoms. Please set `rdkit=True`."
                    )
                rdkit = False

        elif isinstance(func_or_mapping, dict):
            # allow for Atoms as well under the hood
            first_key = next(iter(func_or_mapping.keys()))
            if hasattr(first_key, "GetPDBResidueInfo"):
                func_or_mapping = {
                    atom.GetPDBResidueInfo().GetSerialNumber(): label
                    for atom, label in func_or_mapping.items()
                }
            elif hasattr(first_key, "element") and hasattr(first_key, "serial_number"):
                func_or_mapping = {
                    atom.serial_number: label for atom, label in func_or_mapping.items()
                }
            elif isinstance(first_key, str):
                _func_or_mapping = {}
                for key, label in func_or_mapping.items():
                    matching_atoms = self._raw_molecule.get_atoms(key)
                    for atom in matching_atoms:
                        _func_or_mapping[atom.serial_number] = label
                func_or_mapping = _func_or_mapping

            elif not isinstance(first_key, int):
                raise ValueError(
                    "When providing a mapping, the keys must be either BuildAMol Atoms, RDKit Atoms, or atom serial numbers (ints)."
                )

            def default_label(atom):
                element = atom.GetSymbol()
                is_carbon = element == "C"
                if is_carbon:
                    return ""
                else:
                    return element

            func = lambda atom: func_or_mapping.get(
                atom.GetPDBResidueInfo().GetSerialNumber(), default_label(atom)
            )
            for a in self.mol.GetAtoms():
                a.SetProp("atomLabel", func(a))
            return self

        else:
            raise ValueError(
                f"func_or_mapping must be a valid string identifier, a callable or a dictionary, got {type(func_or_mapping)}."
            )

        if rdkit or (rdkit is None and self._raw_is_rdkit):
            for a in self.mol.GetAtoms():
                a.SetProp("atomLabel", str(func_or_mapping(a)))
        else:

            for a in self.mol.GetAtoms():
                serial = a.GetPDBResidueInfo().GetSerialNumber()
                bam_atom = self._raw_molecule.get_atom(serial, by="serial")
                a.SetProp("atomLabel", str(func_or_mapping(bam_atom)))
        return self

    def highlight_atoms(self, *atoms, color, radius=0.3):
        """
        Highlight atoms in the molecule.

        Parameters
        ----------
        atoms : list
            The Atoms to highlight.
        color
            The color to use for highlighting. This can be either a string, a tuple of RGB values, or a callable that takes an atom and returns a color.
        radius
            The radius to use for highlighting. This can be a float or a callable that takes an atom and returns a float.
        """
        if isinstance(atoms[0], (list, tuple, set)) and len(atoms) == 1:
            atoms = atoms[0]

        a = atoms[0]
        if isinstance(a, (str, int)):
            if self._raw_is_rdkit:
                raise ValueError(
                    "When providing atom ids or serial numbers the underlying molecule cannot be an RDKit molecule. Please provide RDKit or BuildAMol Atoms directly."
                )
            atoms = self._raw_molecule.get_atoms(atoms)
        elif hasattr(a, "element") and hasattr(a, "serial_number"):
            pass
        elif hasattr(a, "GetPDBResidueInfo"):
            if not self._raw_is_rdkit:
                atoms = [
                    self._raw_molecule.get_atom(a.GetPDBResidueInfo().GetSerialNumber())
                    for a in atoms
                ]

        self._atoms_to_highlight.extend(atoms)
        if callable(color):
            color = {atom: color(atom) for atom in atoms}
        else:
            color = {atom: color for atom in atoms}
        self._atoms_highlight_colors.update(color)

        if callable(radius):
            radius = {atom: radius(atom) for atom in atoms}
        else:
            radius = {atom: radius for atom in atoms}
        self._atoms_highlight_radii.update(radius)
        return self

    def highlight_bonds(self, *bonds, color=None):
        """
        Highlight bonds in the molecule.

        Parameters
        ----------
        bonds : list
            The bonds (tuples of BuildAMol Atoms) to highlight.
        color
            The color to use for highlighting. This can be either a string or a tuple of RGB values, or a callable that takes a bond and returns a color.
        """
        if isinstance(bonds[0], (list, tuple, set)) and len(bonds) == 1:
            bonds = bonds[0]


        if callable(color):
            bonds = {bond: color(bond) for bond in bonds}
        else:
            bonds = {bond: color or self.highlight_color for bond in bonds}

        self._bonds_to_highlight.update(bonds)
        return self

    def highlight_residues(self, *residues, color):
        """
        Highlight all bonds and atoms in the given residues.

        Parameters
        ----------
        residues : list
            The residues (BuildAMol Residue objects) whose bonds to highlight.
        color
            The color to use for highlighting. This can be either a string or a tuple of RGB values, or a callable that takes a bond and returns a color.
        """
        if self._raw_is_rdkit:
            raise ValueError(
                "When providing residues the underlying molecule cannot be an RDKit molecule. Please provide RDKit or BuildAMol Atoms directly."
            )
        if isinstance(residues[0], (list, tuple, set)) and len(residues) == 1:
            residues = residues[0]

        bonds = []
        atoms = []
        for residue in residues:
            residue = self._raw_molecule.get_residue(residue)
            atoms.extend(residue.get_atoms())
            bonds.extend(self._raw_molecule.get_bonds(residue))

        self.highlight_bonds(*bonds, color=color).highlight_atoms(*atoms, color=color)
        return self

    def set_colors(self, _element_colors: dict):
        """
        Set the colors for the atoms.

        Parameters
        ----------
        _element_colors : dict
            A dictionary mapping element symbols to colors.
            The keys should be elementy symbols (i.e. 6 for carbon, etc.)
            and the values should be RGB tuples.
        """
        _element_colors = {
            key: colors.to_rgba(value) for key, value in _element_colors.items()
        }
        keys = list(_element_colors.keys())
        for k in keys:
            if isinstance(k, str):
                _k = periodictable.elements.symbol(k).number
                _element_colors[_k] = _element_colors.pop(k)

        self._custom_colors.update(_element_colors)
        return self

    def set_options(self, **kwargs):
        """
        Set the default drawing options.

        Parameters
        ----------
        **kwargs
            Any additional arguments to pass to `MolDrawOptions`.
        """
        aux.deprecation_warning(
            "The `set_options` method is deprecated and will be removed in future versions. Provide drawing options directly as kwargs to the `draw` method instead."
        )
        for k, v in kwargs.items():
            if not k.startswith("_") and hasattr(self.options, k):
                setattr(self.options, k, v)
        return self

    def draw(
        self,
        draw_hydrogens: bool = False,
        linewidth: float = 1,
        fontsize: float = 20,
        width: int = 1000,
        height: int = 500,
        background: tuple = None,
        **kwargs,
    ):
        """
        Generate the 2D image.

        Parameters
        ----------
        draw_hydrogens : bool
            Whether to draw hydrogens.
        linewidth: float
            The linewidth of the bonds.
        fontsize : float
            The font size of the atom labels.
        width : int
            The width of the image in pixels.
        height : int
            The height of the image in pixels.
        background : tuple
            The background color to use. Use `None` for a transparent background.
        **kwargs
            Any additional arguments to pass to the `MolDrawOptions` of the RDKit drawer (either `MolDraw2DSVG` or `MolDraw2DCairo`).

        Returns
        -------
        str or PIL.Image.Image
            The SVG string (if drawer is "svg") or a PIL Image (if drawer is "png").
        """
        if not draw_hydrogens:
            mol = Chem.rdmolops.RemoveHs(self.mol)
        else:
            mol = self.mol

        drawer = (
            Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
            if self._drawer_type == "svg"
            else Draw.rdMolDraw2D.MolDraw2DCairo(width, height)
        )

        draw_options = drawer.drawOptions()

        for k, v in self.options.__dict__.items():
            setattr(draw_options, k, v)

        for k, v in self.options.__dict__.items():
            setattr(draw_options, k, v)
        for k, v in kwargs.items():
            if not k.startswith("_") and hasattr(draw_options, k):
                setattr(draw_options, k, v)

        draw_options.bondLineWidth = linewidth
        if self._custom_colors:
            draw_options.updateAtomPalette(self._custom_colors)

        draw_options.fixedFontSize = fontsize

        if background is None:
            draw_options.clearBackground = False
        else:
            if isinstance(background, str):
                background = colors.to_rgba(background)
            draw_options.setBackgroundColour(colors.to_rgba(background))

        kws = self._prepare_highlighting(mol, draw_hydrogens)

        drawer.DrawMoleculeWithHighlights(mol, legend="", **kws)
        drawer.FinishDrawing()
        if self._drawer_type == "svg":
            svg = drawer.GetDrawingText()
            return svg
        else:
            from PIL import Image
            from io import BytesIO

            img = drawer.GetDrawingText()
            img = Image.open(BytesIO(img))
            return img

    def show(self, draw_hydrogens: bool = False, **kwargs):
        """
        Show the molecule

        Parameters
        ----------
        draw_hydrogens : bool
            Whether to draw hydrogens.
        **kwargs
            Any additional keyword arguments to pass to `draw`.
        """
        if self._drawer_type == "svg":
            dpi = kwargs.pop("dpi", 300)
        out = self.draw(draw_hydrogens=draw_hydrogens, **kwargs)
        if self._drawer_type == "svg":
            # turn SVG string into a PIL image
            from io import BytesIO
            from svglib.svglib import svg2rlg
            from reportlab.graphics import renderPM

            img = svg2rlg(BytesIO(out.encode("utf-8")))
            img = renderPM.drawToPIL(img, dpi=dpi)
            img.show()

        else:
            out.show()

    def _prepare_highlighting(self, mol, include_hydrogens: bool):
        kws = {}
        has_atoms_to_highlight = len(self._atoms_to_highlight) > 0
        has_bonds_to_highlight = len(self._bonds_to_highlight) > 0

        kws["highlight_atom_map"] = {}
        kws["highlight_bond_map"] = {}
        kws["highlight_radii"] = {}
        kws["highlight_linewidth_multipliers"] = {}

        if not has_atoms_to_highlight and not has_bonds_to_highlight:
            return kws

        if not include_hydrogens:
            atom_filter = lambda atom: atom.element != "H"
            bond_filter = lambda bond: (
                bond[0].element != "H" and bond[1].element != "H"
            )
        else:
            atom_filter = lambda atom: True
            bond_filter = lambda bond: True

        if has_atoms_to_highlight:
            highlight_atom_map = {
                atom: [colors.to_rgba(color)]
                for atom, color in self._atoms_highlight_colors.items()
                if atom_filter(atom)
            }
            highlight_radii = {
                atom: radius
                for atom, radius in self._atoms_highlight_radii.items()
                if atom_filter(atom)
            }

            highlight_atom_map = {
                self._rdkit_atom_from_buildamol_atom(atom, mol).GetIdx(): color
                for atom, color in highlight_atom_map.items()
            }
            highlight_radii = {
                self._rdkit_atom_from_buildamol_atom(atom, mol).GetIdx(): radius
                for atom, radius in highlight_radii.items()
            }

            kws["highlight_atom_map"] = highlight_atom_map
            kws["highlight_radii"] = highlight_radii

        if has_bonds_to_highlight:

            highlight_bond_map = {
                idx: color
                for idx, color in self._bonds_to_highlight.items()
                if bond_filter(idx)
            }

            highlight_bond_map = {
                bond.GetIdx(): [colors.to_rgba(color)]
                for (atom1, atom2), color in highlight_bond_map.items()
                if (
                    bond := mol.GetBondBetweenAtoms(
                        self._rdkit_atom_from_buildamol_atom(atom1, mol).GetIdx(),
                        self._rdkit_atom_from_buildamol_atom(atom2, mol).GetIdx(),
                    )
                )
            }
            kws["highlight_bond_map"] = highlight_bond_map

        return kws

    def _rdkit_atom_from_buildamol_atom(self, atom, mol=None):
        mol = mol or self.mol
        if hasattr(atom, "GetPDBResidueInfo"):
            return atom
        elif hasattr(atom, "serial_number"):
            _atom = next(
                (
                    _atom
                    for _atom in mol.GetAtoms()
                    if _atom.GetPDBResidueInfo().GetSerialNumber() == atom.serial_number
                ),
                None,
            )
            if _atom is None:
                raise ValueError(f"Atom {atom} did not have an RDKit equivalent.")
        elif isinstance(atom, int):
            _atom = next(
                (_atom for _atom in mol.GetAtoms() if _atom.GetIdx() == atom),
                None,
            )
            if _atom is None:
                raise ValueError(
                    f"Atom with index {atom} did not have an RDKit equivalent."
                )
        else:
            raise ValueError(
                f"Unsupported atom type: {atom.__class__.__name__}. The input has to be a BuildAMol Atom, an RDKit Atom, or an atom index (int)."
            )
        return _atom


class Py3DmolViewer:
    """
    View a molecule in 3D using the py3Dmol library.

    Attributes
    ----------
    view : py3Dmol.view
        The py3Dmol view object.

    Parameters
    ----------
    molecule
        The molecule to view.
    width : int
        The width of the viewer in pixels.
    height : int
        The height of the viewer in pixels.
    style : dict
        The style to apply to the visualization.
    """

    default_style = {"stick": {}}

    def __init__(
        self, molecule, width: int = 500, height: int = 500, style: dict = None
    ) -> None:
        try:
            import py3Dmol
        except ImportError:
            py3Dmol = None

        if py3Dmol is None or Chem is None:
            raise ImportError(
                "py3Dmol and/or rdkit are not available. Please install them and be sure to use a compatible (Jupyter) environment."
            )
        if isinstance(molecule, (list, tuple, set)):
            molecule = aux.AtomIterator(molecule)
        if not hasattr(molecule, "get_atoms"):
            raise ValueError(
                f"Unsupported molecule type: {molecule.__class__.__name__}. The input has to be a Py3DmolViewer, Molecule, or any other class with an 'get_atoms' method that can be converted to PDB."
            )

        if hasattr(molecule, "to_pdb"):
            self.pdb = utils.pdb.encode_pdb(molecule)
        else:
            self.pdb = utils.pdb.make_atoms_table(molecule)

        self.style = dict(Py3DmolViewer.default_style)
        if style:
            self.style.update(style)

        self.view = py3Dmol.view(width=width, height=height)
        self.view.addModel(self.pdb, "pdb")
        self.n_models = 1
        self.view.setStyle(self.style)
        self.view.zoomTo()

    def set_style(self, style: dict, model=None) -> None:
        """
        Set the visualization style.

        Parameters
        ----------
        style : dict
            The style to add.
        model : int
            A specific model to apply the style to.
        """
        if model:
            if model > self.n_models:
                raise ValueError(
                    f"Model {model} does not exist. The viewer contains {self.n_models} models."
                )
            self.view.setStyle({"model": model}, style)
        else:
            self.view.setStyle(style)
        return self

    def add(self, other, style=None):
        """
        Add a second molecule to the viewer.

        Parameters
        ----------
        other
            This may either be another Py3DmolViewer, a molecule object that can be converted to an RDKit molecule.
        """
        if isinstance(other, Py3DmolViewer):
            self.view.addModel(other.pdb, "pdb")
            if style is None:
                style = other.style
        elif hasattr(other, "to_pdb"):
            pdb = utils.pdb.encode_pdb(other)
            self.view.addModel(pdb, "pdb")
            if style is None:
                style = self.style
        elif hasattr(other, "get_atoms"):
            pdb = utils.pdb.make_atoms_table(other)
            self.view.addModel(pdb, "pdb")
            if style is None:
                style = self.style
        elif isinstance(other, str):
            self.view.addModel(other, "pdb")
            if style is None:
                style = self.style
        else:
            raise ValueError(
                f"Unsupported molecule type: {other.__class__.__name__}. The input has to be a Py3DmolViewer or Molecule."
            )

        self.view.setStyle({"model": self.n_models}, style)
        self.n_models += 1
        return self

    def __iadd__(self, other):
        return self.add(other)

    def __add__(self, other):
        return self.add(other)

    def show(self):
        """
        Show the molecule in a Jupyter notebook
        """
        return self.view.show()


class NglViewer:
    """
    View a molecule in 3D using
    the NGLView library.

    Parameters
    ----------
    molecule
        The molecule to view. This may be any object that holds
        a biopython structure e.g. a Molecule, AtomGraph, or ResidueGraph.
    """

    def __init__(self, molecule):
        try:
            import nglview
        except ImportError:
            nglview = None

        if nglview is None:
            raise ImportError(
                "NGLView is not available. Please install it with `pip install nglview` and be sure to use a compatible environment."
            )
        if hasattr(molecule, "to_pdb"):
            self.pdb = utils.pdb.encode_pdb(molecule)
        elif molecule.__class__.__name__ in ("AtomGraph", "ResidueGraph"):
            self.pdb = utils.pdb.encode_pdb(molecule._molecule)
        else:
            raise ValueError(
                f"Unsupported molecule type: {molecule.__class__.__name__}"
            )

    def show(self):
        """
        Show the molecule in a Jupyter notebook
        """
        import nglview
        import io

        f = io.StringIO(self.pdb)
        f.seek(0)
        fig = nglview.show_file(f, ext="pdb")
        return fig


def rgba_to_hex(rgba: tuple) -> str:
    """
    Convert an rgba color to hex.

    Parameters
    ----------
    rgba : tuple
        The rgba color to convert.

    Returns
    -------
    str
        The hex color.
    """
    return "#" + "".join([hex(int(i * 255))[2:] for i in rgba])


class PlotlyViewer3D:
    __continuous_colors__ = [
        "navy",
        "blue",
        "teal",
        "green",
        "lightgreen",
        "yellow",
        "orange",
        "red",
        "crimson",
        "darkred",
        "brown",
        "purple",
        "pink",
    ]

    __atom_colors__ = {
        "C": "darkslategray",
        "O": "red",
        "H": "lightgray",
        "N": "blue",
        "S": "yellow",
        "P": "purple",
        "F": "green",
        "Cl": "green",
        "Br": "green",
        "I": "green",
    }

    def __init__(self) -> None:
        PlotlyViewer3D.reset(self)
        self._color_idx = 0
        self.opacity = default_plotly_opacity
        self.size = default_plotly_marker_size
        self.bond_color = default_plotly_bond_color
        self.bond_linewidth = default_plotly_linewidth

    def _get_color(self):
        color = self.__continuous_colors__[self._color_idx]
        self._color_idx = (self._color_idx + 1) % len(self.__continuous_colors__)
        return color

    def _get_atom_color(self, atom):
        return self.__atom_colors__.get(atom.element.title(), "black")

    def add(self, fig):
        """
        Add a plotly figure to the viewer.
        """
        if isinstance(fig, PlotlyViewer3D):
            data = fig.figure.data
        else:
            data = getattr(fig, "data", fig)
        self.figure.add_traces(data)

    def __add__(self, fig):
        self.add(fig)
        return self

    def show(self):
        self.figure.show()

    def write_html(self, path):
        self.figure.write_html(path)

    def reset(self, **kwargs):
        self.figure = go.Figure(
            layout=go.Layout(
                scene=dict(
                    xaxis=dict(
                        showgrid=False,
                        showline=False,
                        showticklabels=False,
                        range=kwargs.pop("xlim", None),
                    ),
                    yaxis=dict(
                        showgrid=False,
                        showline=False,
                        showticklabels=False,
                        range=kwargs.pop("ylim", None),
                    ),
                    zaxis=dict(
                        showgrid=False,
                        showline=False,
                        showticklabels=False,
                        range=kwargs.pop("zlim", None),
                    ),
                    # aspectmode="cube",
                ),
                template="simple_white",
            )
        )
        return self

    def viewbox(self, xlim=None, ylim=None, zlim=None):
        if isinstance(xlim, (int, float)):
            xlim = [-xlim, xlim]
        if isinstance(ylim, (int, float)):
            ylim = [-ylim, ylim]
        if isinstance(zlim, (int, float)):
            zlim = [-zlim, zlim]
        self.figure.update_layout(
            scene=dict(
                xaxis=dict(range=xlim),
                yaxis=dict(range=ylim),
                zaxis=dict(range=zlim),
            )
        )
        return self

    def update_layout(self, **kwargs):
        self.figure.update_layout(**kwargs)
        return self

    def draw_point(
        self,
        id: str,
        coord,
        color="black",
        opacity=1.0,
        size=5,
        showlegend=True,
        **kwargs,
    ):
        new = go.Scatter3d(
            x=[coord[0]],
            y=[coord[1]],
            z=[coord[2]],
            mode="markers",
            marker=dict(opacity=opacity, color=color, size=size),
            name=id,
            hoverinfo="name",
            showlegend=showlegend,
            **kwargs,
        )
        self.add(new)
        return self

    def draw_vector(
        self,
        id,
        coord_a,
        coord_b,
        color="black",
        linewidth=1.5,
        opacity=1.0,
        showlegend=True,
        hoverinfo: str = "skip",
        elongate: float = 1.0,
        legendgroup: str = None,
    ):
        new = go.Scatter3d(
            x=[coord_a[0], coord_a[0] + (coord_b[0] - coord_a[0]) * elongate],
            y=[coord_a[1], coord_a[1] + (coord_b[1] - coord_a[1]) * elongate],
            z=[coord_a[2], coord_a[2] + (coord_b[2] - coord_a[2]) * elongate],
            mode="lines",
            line=dict(color=color, width=linewidth),
            name=id,
            hoverinfo=hoverinfo,
            opacity=opacity,
            showlegend=showlegend,
            legendgroup=legendgroup,
        )
        self.add(new)
        return self

    def draw_edges(
        self,
        *edges,
        color="black",
        linewidth=1,
        opacity=1.0,
        elongate: float = 1.0,
        showlegend: bool = True,
        name: str = None,
    ):
        for edge in edges:
            self.draw_vector(
                f"{edge[0].id}-{edge[1].id}",
                edge[0].coord,
                edge[1].coord,
                color=color,
                linewidth=linewidth,
                opacity=opacity,
                elongate=elongate,
                showlegend=showlegend,
                legendgroup=name,
            )
        return self

    def draw_points(
        self,
        coords: list,
        ids: list = None,
        colors: list = None,
        opacities: list = None,
        sizes: list = None,
        showlegends: list = None,
        **kwargs,
    ):
        if ids is None:
            ids = [str(i) for i in range(len(coords))]
        if colors is None:
            colors = ["black" for _ in range(len(coords))]
        elif isinstance(colors, str):
            colors = [colors for _ in range(len(coords))]
        if opacities is None:
            opacities = [1.0 for _ in range(len(coords))]
        elif isinstance(opacities, (int, float)):
            opacities = [opacities for _ in range(len(coords))]
        if showlegends is None:
            showlegends = [True for _ in range(len(coords))]
        elif isinstance(showlegends, bool):
            showlegends = [showlegends for _ in range(len(coords))]
        if sizes is None:
            sizes = [self.size for _ in range(len(coords))]
        elif isinstance(sizes, (int, float)):
            sizes = [sizes for _ in range(len(coords))]

        for idx, coord in enumerate(coords):
            self.draw_point(
                ids[idx],
                coord,
                colors[idx],
                opacities[idx],
                sizes[idx],
                showlegends[idx],
                **kwargs,
            )
        return self

    def highlight_atoms(
        self,
        *atoms,
        names: list = None,
        colors: list = None,
        opacity: float = 1,
        size: int = 10,
        showlegend: bool = True,
        hoverinfo: str = "name",
    ):
        if colors is not None and not isinstance(colors, list):
            colors = [colors] * len(atoms)

        atom_scatter = []
        for idx, atom in enumerate(atoms):
            atom = self._src.get_atom(atom)
            if colors is None:
                color = self.__atom_colors__.get(atom.element.title(), "black")
            else:
                color = colors[idx]
            if names is None:
                name = repr(atom)
            else:
                name = names[idx]
            new = go.Scatter3d(
                x=[atom.coord[0]],
                y=[atom.coord[1]],
                z=[atom.coord[2]],
                mode="markers",
                marker=dict(color=color, opacity=opacity, size=size),
                hoverinfo=hoverinfo,
                showlegend=showlegend,
                name=name,
                legendgroup="Highlighted",
            )
            atom_scatter.append(new)
        self.add(atom_scatter)
        return self

    def highlight_residues(
        self,
        *residues,
        bond_colors: list = None,
        opacity: float = 0.6,
        linewidth: float = 2,
        draw_atoms: bool = False,
    ):
        if not isinstance(bond_colors, list):
            bond_colors = [bond_colors] * len(residues)

        residue_traces = []
        for idx, residue in enumerate(residues):
            residue = self._src.get_residue(residue)
            atoms = self._atom_df[self._atom_df["residue_serial"] == residue.id[1]]
            bonds = self._bond_df[
                self._bond_df["a"].isin(atoms.index)
                & self._bond_df["b"].isin(atoms.index)
            ]
            if bond_colors:
                bonds.loc[:, "bond_color"] = bond_colors[idx]
            bonds.loc[:, "bond_order"] = bonds["bond_order"] + linewidth
            _op = self.opacity
            self.opacity = opacity

            fig = self._setup_fig(atoms, bonds, draw_atoms=draw_atoms)

            residue_traces.extend(fig.data)
            self.opacity = _op
            bonds.loc[:, "bond_order"] = bonds["bond_order"] - linewidth
        self.add(residue_traces)
        return self

    def draw_atoms(
        self,
        *atoms,
        names: list = None,
        colors: list = None,
        opacity: float = None,
        showlegend: bool = True,
        hoverinfo: str = "name",
    ):
        if not opacity:
            opacity = self.opacity
        self.highlight_atoms(
            *atoms,
            names=names,
            colors=colors,
            opacity=opacity,
            showlegend=showlegend,
            hoverinfo=hoverinfo,
        )
        return self

    def draw_residues(
        self,
        *residues,
        bond_colors: list = None,
        opacity: float = None,
        linewidth: float = 2,
    ):
        if not opacity:
            opacity = self.opacity
        self.highlight_residues(
            *residues, bond_colors=bond_colors, opacity=opacity, linewidth=linewidth
        )
        return self

    def draw_atom(self, atom, id=None, color=None, opacity=None, size=None):
        if color is None:
            color = self.__atom_colors__.get(atom.element)
        if opacity is None:
            opacity = min(1, self.opacity * 2)
        if id is None:
            id = str(atom.id) + " " + str(atom.serial_number)
        self.draw_point(
            id,
            atom.coord,
            color,
            opacity,
            size,
        )
        return self

    def draw_bond(
        self,
        atom_a,
        atom_b,
        color="black",
        linewidth=1.5,
        showlegend=True,
        elongate: float = 1.0,
    ):
        self.draw_vector(
            f"{atom_a.id}-{atom_b.id}",
            atom_a.coord,
            atom_b.coord,
            color,
            linewidth,
            showlegend,
            elongate=elongate,
        )
        return self


class MoleculeViewer3D(PlotlyViewer3D):
    def make_df(self, mol) -> tuple:
        _atom_df = {
            "x": [atom.coord[0] for atom in mol.get_atoms()],
            "y": [atom.coord[1] for atom in mol.get_atoms()],
            "z": [atom.coord[2] for atom in mol.get_atoms()],
            "atom_id": [atom.id for atom in mol.get_atoms()],
            "atom_serial": [atom.serial_number for atom in mol.get_atoms()],
            "atom_element": [atom.element.title() for atom in mol.get_atoms()],
            "residue_serial": [atom.get_parent().id[1] for atom in mol.get_atoms()],
            "residue_name": [atom.get_parent().resname for atom in mol.get_atoms()],
            "chain_id": [atom.get_parent().get_parent().id for atom in mol.get_atoms()],
        }
        _atom_df = pd.DataFrame(_atom_df)
        _atom_df.set_index("atom_serial", drop=False, inplace=True)

        _bond_df = {
            "a": [i[0].serial_number for i in mol.get_bonds()],
            "b": [i[1].serial_number for i in mol.get_bonds()],
            "bond_color": [self.bond_color for i in mol.get_bonds()],
            "bond_order": [self.bond_linewidth * i.order for i in mol.get_bonds()],
        }
        _bond_df = pd.DataFrame(_bond_df)

        return _atom_df, _bond_df

    def link(self, mol: "Molecule"):
        """
        Link a source molecule to the viewer.
        """
        self._src = mol
        atom_df, bond_df = self.make_df(mol)
        self._atom_df = atom_df
        self._bond_df = bond_df

    def setup(self, draw_atoms=True):
        """
        Setup the viewer with the molecule.
        """
        self.add(self._setup_fig(self._atom_df, self._bond_df, draw_atoms=draw_atoms))

    def _setup_fig(self, atom_df, bond_df, draw_atoms=True):
        if not draw_atoms:
            fig = go.Figure()
        else:
            atom_df["__marker_size"] = self.size
            fig = px.scatter_3d(
                atom_df,
                x="x",
                y="y",
                z="z",
                color="atom_element",
                color_discrete_map=self.__atom_colors__,
                opacity=self.opacity,
                size="__marker_size",
                hover_data=[
                    "atom_id",
                    "atom_serial",
                    "residue_serial",
                    "residue_name",
                    "chain_id",
                ],
                template="none",
            )
        bonds = []
        for i, row in bond_df.iterrows():
            a1 = atom_df.loc[row["a"]]
            a2 = atom_df.loc[row["b"]]
            new = go.Scatter3d(
                x=[a1["x"], a2["x"]],
                y=[a1["y"], a2["y"]],
                z=[a1["z"], a2["z"]],
                mode="lines",
                line=dict(
                    color=row["bond_color"],
                    width=row["bond_order"] ** 2,
                    # opacity=min(1, self.opacity * 2),
                ),
                hoverinfo="skip",
                showlegend=False,
            )
            bonds.append(new)
        fig.add_traces(bonds)

        return fig

    def reset(self):
        self.figure = self._setup_fig(self._atom_df, self._bond_df)

    def rainbow(self):
        """
        Colorize the residues in rainbow colors
        """
        self.highlight_residues(
            *self._src.get_residues(),
            bond_colors=[self._get_color() for i in self._src.get_residues()],
        )


class AtomGraphViewer3D(PlotlyViewer3D):
    def link(self, graph):
        self._src = graph
        self._atom_df, self._bond_df = self.make_df(graph)
        self.add(self._setup_fig(self._atom_df, self._bond_df))

    def make_df(self, graph):
        _atom_df = {
            "x": [atom.coord[0] for atom in graph.nodes],
            "y": [atom.coord[1] for atom in graph.nodes],
            "z": [atom.coord[2] for atom in graph.nodes],
            "atom_id": [atom.id for atom in graph.nodes],
            "atom_serial": [atom.serial_number for atom in graph.nodes],
            "atom_element": [atom.element.title() for atom in graph.nodes],
            "residue_serial": [atom.get_parent().id[1] for atom in graph.nodes],
            "residue_name": [atom.get_parent().resname for atom in graph.nodes],
            "chain_id": [atom.get_parent().get_parent().id for atom in graph.nodes],
        }
        _atom_df = pd.DataFrame(_atom_df)
        _atom_df.set_index("atom_serial", drop=False, inplace=True)

        bond_orders = nx.get_edge_attributes(graph, "bond_order")
        _bond_df = {
            "a": [i[0].serial_number for i in bond_orders.keys()],
            "b": [i[1].serial_number for i in bond_orders.keys()],
            "bond_color": [self.bond_color for i in bond_orders.keys()],
            "bond_order": [self.bond_linewidth * i for i in bond_orders.values()],
        }

        _bond_df = pd.DataFrame(_bond_df)

        return _atom_df, _bond_df

    def _setup_fig(self, atom_df, bond_df):
        fig = px.scatter_3d(
            atom_df,
            x="x",
            y="y",
            z="z",
            color="atom_element",
            color_discrete_map=self.__atom_colors__,
            opacity=self.opacity,
            hover_data=[
                "atom_id",
                "atom_serial",
                "residue_serial",
                "residue_name",
                "chain_id",
            ],
            template="none",
        )
        bonds = []
        for i, row in bond_df.iterrows():
            a1 = atom_df.loc[row["a"]]
            a2 = atom_df.loc[row["b"]]
            new = go.Scatter3d(
                x=[a1["x"], a2["x"]],
                y=[a1["y"], a2["y"]],
                z=[a1["z"], a2["z"]],
                mode="lines",
                line=dict(
                    color=row["bond_color"],
                    width=row["bond_order"] ** 2,
                ),
                opacity=min(1, self.opacity * 2),
                hoverinfo="skip",
                showlegend=False,
            )
            bonds.append(new)
        fig.add_traces(bonds)

        return fig


class ResidueGraphViewer3D(PlotlyViewer3D):
    def link(self, graph):
        for node in graph.nodes:
            if getattr(node, "coord", None) is None:
                node.coord = node.center_of_mass()
        self._src = graph
        self._atom_df, self._bond_df = self.make_df(graph)
        self.add(self._setup_fig(self._atom_df, self._bond_df))

    def make_df(self, graph):
        _atom_df = {
            "_id": [atom.get_id() for atom in graph.nodes],
            "x": [atom.coord[0] for atom in graph.nodes],
            "y": [atom.coord[1] for atom in graph.nodes],
            "z": [atom.coord[2] for atom in graph.nodes],
            "id": [str(atom.id) for atom in graph.nodes],
            "serial": [atom.serial_number for atom in graph.nodes],
            "element_or_resname": [
                getattr(atom, "element", getattr(atom, "resname", "")).title()
                for atom in graph.nodes
            ],
            "parent_id": [str(atom.get_parent().id) for atom in graph.nodes],
            "parent_serial": [
                getattr(atom.get_parent(), "serial_number", -1) for atom in graph.nodes
            ],
        }

        _atom_df = pd.DataFrame(_atom_df)
        _atom_df.set_index("_id", drop=False, inplace=True)

        _bond_df = {
            "a": [i[0].get_id() for i in graph.edges],
            "b": [i[1].get_id() for i in graph.edges],
            "bond_color": [self.bond_color for i in graph.edges],
            "bond_order": [self.bond_linewidth for i in graph.edges],
        }

        _bond_df = pd.DataFrame(_bond_df)

        return _atom_df, _bond_df

    def _setup_fig(self, atom_df, bond_df):
        fig = px.scatter_3d(
            atom_df,
            x="x",
            y="y",
            z="z",
            color="element_or_resname",
            color_discrete_map=self.__atom_colors__,
            opacity=self.opacity,
            hover_data=[
                "id",
                "serial",
                "parent_serial",
                "parent_id",
            ],
            template="none",
        )
        bonds = []
        for i, row in bond_df.iterrows():
            a1 = atom_df.loc[row["a"]]
            a2 = atom_df.loc[row["b"]]
            new = go.Scatter3d(
                x=[a1["x"], a2["x"]],
                y=[a1["y"], a2["y"]],
                z=[a1["z"], a2["z"]],
                mode="lines",
                line=dict(
                    color=row["bond_color"],
                    width=row["bond_order"] ** 2,
                ),
                opacity=min(1, self.opacity * 2),
                hoverinfo="skip",
                showlegend=False,
            )
            bonds.append(new)
        fig.add_traces(bonds)

        return fig

    def rainbow(self):
        for node in self._src.nodes:
            if getattr(node, "element", None) is not None:
                continue
            self.draw_atom(node, color=self._get_color())


if __name__ == "__main__":
    import buildamol as bam

    bam.load_sugars()
    man = bam.molecule("MAN")
    man = man % "14bb" * 2
    man.change_element(1, "Au")
    v = Chem2DViewer(man, drawer="png")
    v.label_atoms(lambda a: a.GetSymbol())
    v.highlight_atoms(
        man.atoms,
        color=lambda a: (0, 0, 1, ((a.mass or 0) / man.mass) ** 0.5),
        # radius=lambda a: 0.2 + 0.2 * (a.mass or 0),
    )

    v.highlight_bonds(man.get_bonds(man.get_residue(1)), color="blue")
    v.highlight_residues(1, color=(1, 1, 0, 0.3))
    v.set_colors(
        {
            6: "pink",
        }
    )
    v.show(draw_hydrogens=False, linewidth=5, background="white")
    pass
    # v = MoleculeViewer3D()
    # v.link(man)
    # v.setup()
    # v.show()
    # man.repeat(5, "14bb")
    # # v = Chem2DViewer(man)
    # # v.show()

    # v = MoleculeViewer3D()
    # v.link(man)
    # v.setup()
    # v.show()
    # atoms = man.atoms[:10]
    # v.draw_atoms(*atoms)
    # v.show()

    # manv = ResidueGraphViewer3D()
    # manv.link(man.make_residue_graph(detailed=True))
    # # manv.highlight_residues(1, bond_colors=["red"])
    # manv.rainbow()
    # manv.show()
    # pass
