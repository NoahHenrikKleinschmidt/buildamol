"""
Visualization auxiliary functions
"""

import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from copy import deepcopy
import pandas as pd
import networkx as nx

try:
    import nglview
except:
    nglview = None


class NglViewer:
    """
    View a molecule or graph object in 3D using
    the NGLView library.

    Parameters
    ----------
    molecule
        The molecule to view. This may be any object that holds
        a biopython structure e.g. a Molecule, AtomGraph, or ResidueGraph.
    """

    def __init__(self, molecule):
        if nglview is None:
            raise ImportError(
                "NGLView is not available. Please install it with `pip install nglview` and be sure to use a compatible environment."
            )
        if molecule.__class__.__name__ in ("Molecule", "AtomGraph", "ResidueGraph"):
            self.structure = molecule.structure.to_biopython()
        else:
            self.structure = molecule.to_biopython()

    def show(self):
        """
        Show the molecule in a Jupyter notebook
        """
        return nglview.show_biopython(self.structure)


class MoleculeViewer3D:
    """
    View a molecule or graph object in 3D

    Parameters
    ----------
    molecule
        The molecule to view. This may be any object that holds
        nodes to draw in an attribute such as "atoms" or "nodes",
        e.g. a Molecule, AtomGraph, or ResidueGraph.
    bonds
        The bonds to draw. This may be a list of tuples or equivalent
        iterable, or some other object that holds such data in an attribute
        such as "bonds", or "edges".
    """

    __atom_colors__ = {
        "C": "black",
        "O": "red",
        "H": "lightgray",
        "N": "blue",
        "S": "yellow",
        "P": "purple",
        "F": "green",
        "Cl": "green",
        "CL": "green",
        "Br": "green",
        "BR": "green",
        "I": "green",
    }

    def __init__(self, molecule, bonds=None):
        self.mol = molecule
        self.opacity = 0.3
        self._bonds_obj = bonds if bonds else molecule

        # preprocess to make sure a residue graph can
        # be drawn without issue as residues by default
        # do not have coordinates
        if self.mol.__class__.__name__ == "ResidueGraph":
            for residue in self.mol.residues:
                if hasattr(residue, "coord"):
                    continue
                residue.coord = residue.center_of_mass()

        self._fig = self.setup()
        self._backup_fig = deepcopy(self._fig)

    @property
    def figure(self):
        return self._fig

    def get_color(self, atom):
        """
        Get the color of an atom

        Parameters
        ----------
        atom : Bio.PDB.Atom
            The atom to get the color of

        Returns
        -------
        color : str
            The color of the atom
        """
        return self.__atom_colors__.get(atom.element)

    def update_colors(self, color_dict):
        """
        Update the colors of the atoms

        Parameters
        ----------
        color_dict : dict
            A dictionary of the form {atom: color}
        """
        self.__atom_colors__.update(color_dict)

    def highlight(self, ids):
        """
        Highlight a set of atoms or residues from their ids to draw them at full opacity.

        Parameters
        ----------
        ids : list
            The ids of the atoms to highlight
        """
        atoms = {i.id: i for i in self._get_atoms(self.mol)}
        for i in ids:
            atom = atoms.get(i)
            if not atom:
                continue
            _ = self._fig.add_trace(
                go.Scatter3d(
                    x=[atom.coord[0]],
                    y=[atom.coord[1]],
                    z=[atom.coord[2]],
                    mode="markers",
                    marker=dict(color=self.__atom_colors__.get(i[0])),
                    name=i,
                    showlegend=False,
                ),
            )

    def draw_point(self, id, coords, color="black", opacity=1.0, showlegend=True):
        """
        Draw a point on the figure

        Parameters
        ----------
        id : str
            The id of the point
        coords : list
            The coordinates of the point
        color : str
            The color of the point
        opacity : float
            The opacity of the point
        showlegend : bool
            Whether to show the legend or not
        """
        new = go.Scatter3d(
            x=[coords[0]],
            y=[coords[1]],
            z=[coords[2]],
            mode="markers",
            marker=dict(opacity=opacity, color=color),
            name=id,
            hoverinfo="name",
            showlegend=showlegend,
        )
        _ = self._fig.add_trace(new)

    def draw_vector(self, id, point, vec, color="black", showlegend=True):
        """
        Draw a vector on the figure

        Parameters
        ----------
        id : str
            The id of the vector
        point : list
            The coordinates of the point from which to start the vector
        vec : list
            The vector to draw
        color : str
            The color of the vector
        showlegend : bool
            Whether to show the legend or not
        """
        new = go.Scatter3d(
            x=[point[0], point[0] + vec[0]],
            y=[point[1], point[1] + vec[1]],
            z=[point[2], point[2] + vec[2]],
            mode="lines",
            line=dict(color=color, width=10),
            name=id,
            hoverinfo="skip",
            showlegend=showlegend,
        )
        _ = self._fig.add_trace(new)

    def save(self, filename=None):
        """
        Save the current state of the figure, to revert to this
        when calling the `reset()` method. If also a filename is given,
        the figure is saved to a file.
        """
        self._backup_fig = deepcopy(self._fig)
        if filename:
            self._fig.write_html(filename)

    def reset(self):
        """
        Reset the figure to the last saved state
        """
        self._fig = deepcopy(self._backup_fig)

    def show(self):
        self._fig.show()

    def draw_edges(self, edges=None, color="black", linewidth=1, opacity=1.0):
        """
        Draw edges on the figure

        Parameters
        ----------
        edges : list
            The edges to draw
        color : str
            The color of the edges
        linewidth : int
            The width of the edges
        opacity : float
            The opacity of the edges
        """
        if not edges:
            edges = self._get_bonds(self._bonds_obj)
        for edge in edges:
            a1, a2 = edge
            new = go.Scatter3d(
                x=[a1.coord[0], a2.coord[0]],
                y=[a1.coord[1], a2.coord[1]],
                z=[a1.coord[2], a2.coord[2]],
                mode="lines",
                line=dict(color=color, width=linewidth),
                name=f"{a1.id}-{a2.id}",
                hoverinfo="skip",
                showlegend=False,
                opacity=opacity,
            )
            _ = self._fig.add_trace(new)

    def setup(self):
        """
        Draw the basic molecule setup with all atoms and bonds.
        """
        # draw all atoms
        atoms = self._get_atoms(self.mol)
        self._fig = px.scatter_3d(
            x=[i.coord[0] for i in atoms],
            y=[i.coord[1] for i in atoms],
            z=[i.coord[2] for i in atoms],
            color=[i.id[0] for i in atoms],
            color_discrete_map=self.__atom_colors__,
            opacity=self.opacity,
            hover_name=[i.id for i in atoms],
            template="none",
        )

        # draw all bonds
        bonds = self._get_bonds(self._bonds_obj)
        atoms = {i.id: i for i in atoms}
        for bond in bonds:
            a1 = atoms.get(bond[0], bond[0])
            a2 = atoms.get(bond[1], bond[1])
            if not a1 or not a2:
                continue
            new = go.Scatter3d(
                x=[a1.coord[0], a2.coord[0]],
                y=[a1.coord[1], a2.coord[1]],
                z=[a1.coord[2], a2.coord[2]],
                mode="lines",
                line=dict(color="black", width=1),
                hoverinfo="skip",
                showlegend=False,
            )
            _ = self._fig.add_trace(new)

        self._fig.update_scenes(
            xaxis_showgrid=False,
            xaxis_showline=False,
            xaxis_showticklabels=False,
            yaxis_showgrid=False,
            yaxis_showline=False,
            yaxis_showticklabels=False,
            zaxis_showgrid=False,
            zaxis_showline=False,
            zaxis_showticklabels=False,
        )
        return self._fig

    @staticmethod
    def _get_atoms(obj):
        if hasattr(obj, "nodes"):
            return obj.nodes
        elif hasattr(obj, "atoms"):
            return obj.atoms
        elif hasattr(obj, "get_atoms"):
            return obj.get_atoms()
        elif isinstance(obj, (tuple, list, np.ndarray)):
            return obj
        else:
            raise TypeError("Unknown object type")

    @staticmethod
    def _get_bonds(obj):
        if hasattr(obj, "edges"):
            return obj.edges
        elif hasattr(obj, "bonds"):
            return obj.bonds
        elif hasattr(obj, "get_bonds"):
            return obj.get_bonds()
        elif isinstance(obj, (tuple, list, np.ndarray)):
            return obj
        else:
            raise TypeError("Unknown object type")


class MoleculeViewer3D_v2:
    __atom_colors__ = {
        "C": "black",
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
        self.opacity = 0.3
        self.bond_color = "black"

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
        _atom_df.set_index("atom_serial", inplace=True)

        bonds = nx.get_edge_attributes(mol._AtomGraph, "bond_order")
        _bond_df = {
            "a": [i[0].serial_number for i in bonds.keys()],
            "b": [i[1].serial_number for i in bonds.keys()],
            "bond_color": [self.bond_color for i in bonds.keys()],
            "bond_order": [i for i in bonds.values()],
        }
        _bond_df = pd.DataFrame(_bond_df)

        return _atom_df, _bond_df

    def link(self, mol: "Molecule"):
        self._mol = mol
        atom_df, bond_df = self.make_df(mol)
        self._atom_df = atom_df
        self._bond_df = bond_df
        self._fig = self._setup_fig(atom_df, bond_df)
        self._fig.update_scenes(
            xaxis_showgrid=False,
            xaxis_showline=False,
            xaxis_showticklabels=False,
            yaxis_showgrid=False,
            yaxis_showline=False,
            yaxis_showticklabels=False,
            zaxis_showgrid=False,
            zaxis_showline=False,
            zaxis_showticklabels=False,
        )

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
                line=dict(color=row["bond_color"], width=row["bond_order"] + 1),
                hoverinfo="skip",
                showlegend=False,
            )
            bonds.append(new)
        fig.add_traces(bonds)

        return fig

    def show(self):
        self._fig.show()

    def write_html(self, path):
        self._fig.write_html(path)

    def reset(self):
        self._fig = self._setup_fig(self._atom_df, self._bond_df)

    def highlight_atoms(
        self,
        *atoms,
        names: list = None,
        colors: list = None,
        opacity: float = 1,
        showlegend: bool = True,
    ):
        atom_scatter = []
        for idx, atom in enumerate(atoms):
            atom = self._mol.get_atom(atom)
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
                marker=dict(color=color, opacity=opacity, size=10),
                hoverinfo="skip",
                showlegend=showlegend,
                name=name,
            )
            atom_scatter.append(new)
        self._fig.add_traces(atom_scatter)

    def highlight_residues(
        self, *residues, bond_colors: list = None, opacity: float = 1
    ):
        residue_traces = []
        for idx, residue in enumerate(residues):
            residue = self._mol.get_residue(residue)
            atoms = self._atom_df[self._atom_df["residue_serial"] == residue.id[1]]
            bonds = self._bond_df[
                self._bond_df["a"].isin(atoms.index)
                & self._bond_df["b"].isin(atoms.index)
            ]
            if bond_colors:
                bonds["bond_color"] = bond_colors[idx]
            bonds["bond_order"] += 3
            _op = self.opacity
            self.opacity = opacity
            fig = self._setup_fig(atoms, bonds)
            residue_traces.extend(fig.data)
            self.opacity = _op
        self._fig.add_traces(residue_traces)


if __name__ == "__main__":
    import biobuild as bb

    man = bb.molecule("MAN")
    manv = MoleculeViewer3D_v2()
    manv.link(man)
    manv.highlight_residues(1, bond_colors=["red"])
    manv.show()
