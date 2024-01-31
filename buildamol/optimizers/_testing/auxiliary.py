"""
Auxiliary functions
"""
import buildamol.structural as structural
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist


def rotatron_factory(rotatron_class, graph, rotatable_edges, **params):
    """
    Factory function for Rotatron environments

    Parameters
    ----------
    rotatron_class : class
        The Rotatron class to instantiate
    graph : nx.Graph
        The graph to use
    rotatable_edges : list
        The rotatable edges
    params : dict
        Additional parameters

    Returns
    -------
    Rotatron
        The Rotatron environment
    """
    return rotatron_class(graph, rotatable_edges, **params)


def graph_factory(molecule, **params):
    """
    Factory function for graphs

    Parameters
    ----------
    molecule : Molecule
        The molecule to use
    params : dict
        Additional parameters

    Returns
    -------
    graph, edges
        The graph and rotatable edges
    """
    graph = molecule.get_atom_graph()
    edges = graph.edges
    return graph, edges


def count_clashes(graph, clash_cutoff=0.5):
    gen = getattr(graph, "nodes", None)
    gen = gen or getattr(graph, "get_atoms")
    coords = np.array([a.coord for a in gen()])
    dists = cdist(coords, coords)
    np.fill_diagonal(dists, np.inf)
    clashes = np.sum(dists < clash_cutoff)
    return clashes // 2


def compute_clashes(dists, clash_cutoff=0.5):
    clashes = np.sum(dists < clash_cutoff)
    return clashes // 2


def apply_solution(sol: "np.ndarray", env, molecule):
    for angle, edge in zip(sol, env.rotatable_edges):
        molecule.rotate_around_bond(
            edge[0].serial_number,
            edge[1].serial_number,
            angle,
            descendants_only=True,
            angle_is_degrees=False,
        )
    return molecule


def transform_to_df(dict, other_dict=None, name_a=None, name_b=None):
    name_a = name_a or 0
    name_b = name_b or 1

    df = pd.DataFrame.from_dict(dict)
    df = df.melt()
    df = df.rename(columns={"variable": "key", "value": name_a})
    if other_dict is not None:
        other_df = pd.DataFrame.from_dict(other_dict)
        other_df = other_df.melt()
        other_df = other_df.rename(columns={"variable": "key", "value": name_b})
        df[name_b] = other_df[name_b]
    return df
