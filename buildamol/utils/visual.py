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
    """

    def __init__(self, molecule, highlight_color: str = "cyan"):
        if Chem is None:
            raise ImportError(
                "rdkit is not available. Please install it and be sure to use a compatible environment."
            )
        if hasattr(molecule, "to_rdkit"):
            mol = molecule.to_rdkit()
        elif molecule.__class__.__name__ in ("AtomGraph", "ResidueGraph"):
            mol = molecule._molecule.to_rdkit()
        elif "Chem" in str(molecule.__class__.mro()[0]):
            mol = molecule
        else:
            raise ValueError(
                f"Unsupported molecule type: {molecule.__class__.__name__}"
            )
        mol.RemoveAllConformers()
        self.mol = mol
        self._atoms_to_highlight = []
        self._bonds_to_highlight = []
        self.highlight_color = highlight_color

    def draw(self, draw_hydrogens: bool = False, width: int = 1000, height: int = 500):
        """
        Generate the 2D image.

        Parameters
        ----------
        draw_hydrogens : bool
            Whether to draw hydrogens.
        width : int
            The width of the image in pixels.
        height : int
            The height of the image in pixels.
        """
        if not draw_hydrogens:
            mol = Chem.rdmolops.RemoveHs(self.mol)
        else:
            mol = self.mol
        return Draw.MolToImage(
            mol,
            size=(width, height),
            highlightAtoms=self._atoms_to_highlight,
            highlightBonds=self._bonds_to_highlight,
            highlightColor=colors.to_rgb(self.highlight_color),
        )

    def show(self, draw_hydrogens: bool = False):
        """
        Show the molecule

        Parameters
        ----------
        draw_hydrogens : bool
            Whether to draw hydrogens.
        """
        return self.draw(draw_hydrogens=draw_hydrogens).show()

    def highlight_atoms(self, *atoms):
        """
        Highlight atoms in the molecule.

        Parameters
        ----------
        atoms : list
            The BuildAMol Atoms to highlight.
        """
        self._atoms_to_highlight.extend(atom.serial_number for atom in atoms)

    def highlight_bonds(self, *bonds):
        """
        Highlight bonds in the molecule.

        Parameters
        ----------
        bonds : list
            The bonds (tuples of BuildAMol Atoms) to highlight.
        """

        self._bonds_to_highlight.extend(
            self.mol.GetBondBetweenAtoms(a.serial_number, b.serial_number).GetIdx()
            for a, b in bonds
        )


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

    def update_layout(self, **kwargs):
        self.figure.update_layout(**kwargs)

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
            atom_b.coord - atom_a.coord,
            color,
            linewidth,
            showlegend,
            elongate=elongate,
        )


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
            "a": [i[0].serial_number for i in mol.bonds],
            "b": [i[1].serial_number for i in mol.bonds],
            "bond_color": [self.bond_color for i in mol.bonds],
            "bond_order": [self.bond_linewidth * i.order for i in mol.bonds],
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

    v = MoleculeViewer3D()
    v.link(man)
    v.setup()
    v.show()
    man.repeat(5, "14bb")
    # v = Chem2DViewer(man)
    # v.show()

    v = MoleculeViewer3D()
    v.link(man)
    v.setup()
    v.show()
    # atoms = man.atoms[:10]
    # v.draw_atoms(*atoms)
    # v.show()

    # manv = ResidueGraphViewer3D()
    # manv.link(man.make_residue_graph(detailed=True))
    # # manv.highlight_residues(1, bond_colors=["red"])
    # manv.rainbow()
    # manv.show()
    # pass
