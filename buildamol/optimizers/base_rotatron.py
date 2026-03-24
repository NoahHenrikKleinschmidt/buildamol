"""
This is the basic Rotatron environment. It provides the basic functionality for preprocessing a graph into numpy arrays, masking rotatable edges, and evaluating a possible solution.
All other Rotatron environments inherit from this class.
"""

import gymnasium as gym
import numpy as np

from scipy.spatial.distance import cdist

import buildamol.utils.auxiliary as aux
import buildamol.graphs.base_graph as base_graph
import buildamol.structural.base as structural
from copy import deepcopy

from multiprocessing import Pool

__all__ = ["Rotatron"]


class Rotatron(gym.Env):
    """
    The base class for rotational optimization environments.

    Parameters
    ----------
    graph : AtomGraph or ResidueGraph
        The graph to optimize
    rotatable_edges : list
        A list of edges that can be rotated during optimization.
        If None, all non-locked edges are used.
    n_processes : int
        The number of processes to use to speed up the computation of edge masks and lengths
    setup : bool
        Whether to set up the edge masks and lengths during initialization
    numba : bool
        Whether to use numba to speed up the rotation function.
    """

    def __init__(
        self,
        graph: "base_graph.BaseGraph",
        rotatable_edges: list = None,
        n_processes: int = 1,
        setup: bool = True,
        numba: bool = False,
        backend: str = None,
        **kwargs,
    ):
        self.graph = graph
        self.rotatable_edges = self._get_rotatable_edges(graph, rotatable_edges)

        self.node_dict = {n: i for i, n in enumerate(self.graph.nodes)}
        self.n_nodes = len(self.node_dict)
        self.n_edges = len(self.rotatable_edges)

        self.state = self._make_state_from_graph(self.graph).astype(np.float64)
        self._backup_state = self.state.copy()

        self.rotation_unit_masks = np.ones(
            (len(graph.nodes), len(graph.nodes)), dtype=bool
        )

        self.edge_lengths = np.zeros(self.n_edges)
        self.edge_masks = np.zeros((self.n_edges, self.n_nodes), dtype=bool)

        self.n_processes = n_processes

        if setup:
            self._generate_edge_masks(n_processes=n_processes)
            self._generate_edge_lengths()

        self._edge_node_coords = np.array(
            [[self.node_dict[e[0]], self.node_dict[e[1]]] for e in self.rotatable_edges]
        )

        if backend is None:
            if (
                numba
                or aux.USE_ALL_NUMBA
                or (self.n_edges * self.n_nodes > 10000 and aux.USE_NUMBA)
            ):
                backend = "numba"
            else:
                backend = aux.resolve_compute_backend(None)

        self.backend = aux.resolve_compute_backend(backend)
        self._backend_contexts = {}
        self._backend_states = {}
        self._setup_backend_contexts()
        self.set_backend(self.backend)

    def eval(self, state):
        """
        Calculate the evaluation score for a given state

        Parameters
        ----------
        state : np.ndarray
            The state of the environment

        Returns
        -------
        float
            The evaluation for the state
        """
        return np.inf

    def step(self, action, backend: str = None):
        """
        Take a step in the environment

        Parameters
        ----------
        action : np.ndarray
            The action to take
        backend : str, optional
            The backend to use for computation. If None, uses the
            environment's configured backend.

        Returns
        -------
        np.ndarray
            The new state of the environment
        float
            The evaluation for the new state
        bool
            Whether the environment is done
        dict
            Additional information
        """
        if backend is None:
            return self._step_handle(action)

        step_handle = self._resolve_backend_handle("step", backend)
        return step_handle(action)

    def _step_numpy(self, action):
        action = np.asarray(action)
        new_state = self.state
        for edge in range(self.n_edges):
            new_state = self._normal_rotate(new_state, edge, action[edge])

        self._backend_states["numpy"] = new_state
        e = self.eval(new_state)
        done = self.is_done(new_state)
        return new_state, e, done, {}

    def _step_numba(self, action):
        action = np.asarray(action)
        new_state = self.state
        for edge in range(self.n_edges):
            new_state = self._numba_rotate(new_state, edge, action[edge])

        self._backend_states["numpy"] = new_state
        e = self.eval(new_state)
        done = self.is_done(new_state)
        return new_state, e, done, {}

    def _step_jax(self, action):
        context = self._get_backend_context("jax")
        state = self._backend_states.get("jax", None)
        if state is None:
            state = self.to_backend_state("jax", self.state, cache=True)

        action_j = self.to_backend_state("jax", action, cache=False)
        for edge in range(self.n_edges):
            state = self._jax_rotate(state, edge, action_j[edge], context)

        self._backend_states["jax"] = state

        # Current environment eval/is_done implementations are NumPy-first.
        # Convert explicitly here until backend-native eval paths are implemented.
        eval_state = self.to_numpy_state(state)
        e = self.eval(eval_state)
        done = self.is_done(eval_state)
        return state, e, done, {}

    def is_done(self, state):
        """
        Check whether the environment is done

        Parameters
        ----------
        state : np.ndarray
            The state of the environment

        Returns
        -------
        bool
            Whether the environment is done
        """
        return False

    def reset(self, *args, **kwargs):
        """
        Reset the environment
        """
        self.state[:, :] = self._backup_state
        self._backend_states["numpy"] = self.state
        for _backend in list(self._backend_states.keys()):
            if _backend == "numpy":
                continue
            self._backend_states[_backend] = self.to_backend_state(
                _backend, self._backup_state, cache=False
            )

    def blank(self):
        """
        A blank action
        """
        return np.zeros(len(self.rotatable_edges))

    def copy(self):
        """
        Make a deep copy of the environment
        """
        return deepcopy(self)

    def _setup_backend_contexts(self):
        """
        Set up backend contexts and cached backend state containers.
        """
        self._backend_contexts["numpy"] = {
            "edge_masks": self.edge_masks,
            "edge_lengths": self.edge_lengths,
            "edge_node_coords": self._edge_node_coords,
        }
        self._backend_states["numpy"] = self.state

    def _resolve_backend_handle(self, handle: str, backend: str = None):
        """
        Resolve a backend-specific handle.

        For a handle name ``step`` and backend ``jax``, this looks for
        ``_step_jax``. If missing, it falls back to ``_step_numpy``.
        """
        backend = aux.resolve_compute_backend(backend or self.backend)
        backend_name = f"_{handle}_{backend}"
        default_name = f"_{handle}_numpy"

        fn = getattr(self, backend_name, None)
        if fn is not None:
            return fn

        fn = getattr(self, default_name, None)
        if fn is not None:
            return fn

        raise AttributeError(
            f"No handle found for '{handle}' with backend '{backend}'."
        )

    def set_backend(self, backend: str = None):
        """
        Set the environment backend and bind backend-specific method handles.
        """
        self.backend = aux.resolve_compute_backend(backend or self.backend)
        self._step_handle = self._resolve_backend_handle("step", self.backend)
        self._rotate = self._resolve_backend_handle("rotate", self.backend)
        self._state_to_backend_handle = self._resolve_backend_handle(
            "to_backend_state", self.backend
        )
        self._build_context_handle = self._resolve_backend_handle(
            "build_backend_context", self.backend
        )

    def _build_backend_context(self, backend: str):
        """
        Backward-compatible wrapper for backend context builders.
        """
        return self._resolve_backend_handle("build_backend_context", backend)(backend)

    def _build_backend_context_numpy(self, backend: str = "numpy"):
        """
        Build a backend-specific immutable context from NumPy source arrays.
        """
        return self._backend_contexts["numpy"]

    def _build_backend_context_numba(self, backend: str = "numba"):
        return self._backend_contexts["numpy"]

    def _build_backend_context_jax(self, backend: str = "jax"):
        jnp = aux.get_jax_numpy()
        return {
            "edge_masks": jnp.asarray(self.edge_masks),
            "edge_lengths": jnp.asarray(self.edge_lengths),
            "edge_node_coords": jnp.asarray(self._edge_node_coords),
        }

    def _get_backend_context(self, backend: str):
        """
        Get or lazily initialize a backend context.
        """
        backend = aux.resolve_compute_backend(backend)
        ctx = self._backend_contexts.get(backend, None)
        if ctx is None:
            ctx = self._build_backend_context(backend)
            self._backend_contexts[backend] = ctx
        return ctx

    def to_backend_state(self, backend: str = None, state=None, cache: bool = True):
        """
        Port a NumPy state/action array to a backend-native representation.
        """
        backend = aux.resolve_compute_backend(backend or self.backend)
        if state is None:
            state = self.state

        out = self._resolve_backend_handle("to_backend_state", backend)(state)

        if cache:
            self._backend_states[backend] = out
        return out

    def _to_backend_state_numpy(self, state):
        return np.asarray(state)

    def _to_backend_state_numba(self, state):
        return np.asarray(state)

    def _to_backend_state_jax(self, state):
        jnp = aux.get_jax_numpy()
        return jnp.asarray(state)

    def to_numpy_state(self, state=None, update_state: bool = False):
        """
        Convert a backend-native state back to a NumPy array.

        This is the explicit back-conversion hook intended to be called at the
        end of optimization runs that used non-NumPy backends.
        """
        if state is None:
            state = self._backend_states.get(self.backend, self.state)
        out = np.asarray(state)
        if update_state:
            self.state[:, :] = out
            self._backend_states["numpy"] = self.state
        return out

    def _make_state_from_graph(self, graph):
        """
        Set up the state of the environment
        """
        state = np.array([i.coord for i in graph.nodes])
        return state

    def _make_state_from_dict(self, dict):
        """
        Set up the state of the environment
        """
        state = np.array([v for v in dict.values()])
        return state

    def _get_rotatable_edges(self, graph, rotatable_edges):
        """
        Get the rotatable edges

        Parameters
        ----------
        graph : AtomGraph or ResidueGraph
            The graph to optimize
        rotatable_edges : list
            A list of edges that can be rotated during optimization.
            If None, all non-locked edges are used.

        Returns
        -------
        list
            The rotatable edges
        """
        if rotatable_edges is None:
            _circulars = graph.nodes_in_cycles
            rotatable_edges = [
                e
                for e in graph.edges
                if e not in graph._locked_edges
                and graph.edges[e].get("bond_order", 1) == 1
                and not "Residue" in type(e[0]).__name__
                and not "Residue" in type(e[1]).__name__
                and not (e[0] in _circulars and e[1] in _circulars)
                and len(graph.get_descendants(*e)) > 1
            ]
        return rotatable_edges

    def _generate_rotation_unit_masks(self):
        """
        Generate a boolean mask (n_nodes, n_nodes) where all nodes
        that are part of the same rotation unit are set to False
        """
        dists1 = cdist(self.state, self.state)
        for i, angle in enumerate(np.random.random(self.n_edges)):
            state2 = self._rotate_numpy(self.state, i, angle)
        dists2 = cdist(state2, state2)
        for i, angle in enumerate(np.random.random(self.n_edges)):
            state3 = self._rotate_numpy(state2, i, angle)
        dists3 = cdist(state3, state3)

        d12 = np.abs(dists1 - dists2) < 1e-4
        d13 = np.abs(dists1 - dists3) < 1e-4
        d23 = np.abs(dists2 - dists3) < 1e-4

        dists = np.sum([d12, d13, d23], axis=0) == 3
        self.rotation_unit_masks = ~dists
        self.reset()

    def _find_rotation_units(self):
        self.rotation_units = {}
        patterns = []
        rdx = 0
        for edx, mask in enumerate(self.rotation_unit_masks):
            pattern = ~mask
            if any(np.all(i == pattern) for i in patterns):
                pdx = next(
                    idx for idx, i in enumerate(patterns) if np.all(i == pattern)
                )
                self.rotation_units[pdx].add(edx)
                continue
            patterns.append(pattern)
            self.rotation_units[rdx] = {edx}
            rdx += 1
        self.rotation_units = {
            r: np.array(list(v)) for r, v in self.rotation_units.items()
        }

    def _generate_edge_lengths(self):
        """
        Compute the lengths of the edges
        """
        self.edge_lengths = np.array(
            [
                np.linalg.norm(
                    self.state[self.node_dict[e[0]]] - self.state[self.node_dict[e[1]]]
                )
                for e in self.rotatable_edges
            ]
        )

    def _generate_edge_masks(self, n_processes):
        """
        Compute the edge masks of downstream nodes
        """
        if n_processes > 1:
            p = Pool(n_processes)
            self.edge_masks = np.array(
                p.map(self._generate_edge_mask, [e for e in self.rotatable_edges]),
                dtype=bool,
            )
            p.close()
            p.join()
        else:
            self.edge_masks = np.array(
                [self._generate_edge_mask(e) for e in self.rotatable_edges],
                dtype=bool,
            )

    def _generate_edge_mask(self, edge):
        return np.array(
            [
                1 if i in self.graph.get_descendants(*edge) else 0
                for i in self.graph.nodes
            ]
        )

    def _rotate_numpy(self, state, edx, angle):
        if -1e-3 < angle < 1e-3:
            return state

        mask = self.edge_masks[edx]

        # vec = self._get_edge_vector(edx)
        adx, bdx = self._edge_node_coords[edx]
        vec = state[bdx] - state[adx]
        vec /= self.edge_lengths[edx]

        # ref_coord = self._get_edge_ref_coord(edx)
        ref_coord = state[adx]

        state[mask] = (
            structural.rotate_coords(state[mask] - ref_coord, angle, vec) + ref_coord
        )
        return state

    # Backward-compatible alias
    _normal_rotate = _rotate_numpy

    def _rotate_numba(self, state, edx, angle):
        if -1e-3 < angle < 1e-3:
            return state

        return _numba_wrapper_rotate(
            state,
            edx,
            angle,
            self._edge_node_coords,
            self.edge_lengths,
            self.edge_masks,
        )

    # Backward-compatible alias
    _numba_rotate = _rotate_numba

    def _rotate_jax(self, state, edx, angle, context=None):
        if -1e-3 < angle < 1e-3:
            return state

        if context is None:
            context = self._get_backend_context("jax")

        jnp = aux.get_jax_numpy()
        mask = context["edge_masks"][edx]
        adx, bdx = context["edge_node_coords"][edx]

        vec = state[bdx] - state[adx]
        vec = vec / context["edge_lengths"][edx]
        ref_coord = state[adx]

        rotated = structural.rotate_coords(
            state[mask] - ref_coord,
            angle,
            vec,
            backend="jax",
        )
        return state.at[mask].set(jnp.asarray(rotated) + ref_coord)

    # Backward-compatible alias
    _jax_rotate = _rotate_jax

    # ============================================================
    # The setup helper functions can be used by other environments
    # that inherit from this base class
    # ============================================================

    def _setup_helpers_crop_faraway_nodes(self, radius, graph=None, edges=None):
        """
        This is a helper function to remove nodes that are too far away from the rotatable edges.
        """
        if graph and edges:
            rotatable_edges = self._get_rotatable_edges(graph, edges)
        else:
            rotatable_edges = self.rotatable_edges

        edge_coords = np.array([(a.coord + b.coord) / 2 for a, b in rotatable_edges])

        nodes = list(graph.nodes)
        node_coords = np.array([node.coord for node in nodes])

        dists = cdist(edge_coords, node_coords)
        dists = dists > radius

        dists = np.apply_along_axis(np.all, 0, dists)

        if np.max(dists) != 0:
            nodes_to_drop = [nodes[i] for i, d in enumerate(dists) if d]
            graph.remove_nodes_from(nodes_to_drop)

        return graph, rotatable_edges


# ============================================================
# The numba functions are used to speed up things.
# For each function there must be a _normal_ and a _numba_
# version. The _normal_ version is used in the setup if numba is not installed
# The _numba_ version is used in the step function if numba is installed
# ============================================================


@aux.njit
def _numba_wrapper_rotate(
    state, edx, angle, edge_node_coords, edge_lengths, edge_masks
):
    """
    Rotate the graph around an edge. This is the version that is used in the step function.

    Parameters
    ----------
    edx : int
        The edge index to rotate around
    angle : float
        The angle to rotate by

    Returns
    -------
    np.ndarray
        The new state of the environment
    """
    mask = edge_masks[edx]
    adx, bdx = edge_node_coords[edx]
    vec = state[bdx] - state[adx]
    vec /= edge_lengths[edx]

    ref_coord = state[adx]

    rot = structural._numba_wrapper_rotation_matrix(vec, angle)
    rot = np.transpose(rot)

    state[mask] -= ref_coord
    _c = state[mask]
    _c = np.dot(_c, rot, out=_c)
    state[mask] = _c
    state[mask] += ref_coord

    return state


# if __name__ == "__main__":
#     import buildamol as bam

#     bam.load_sugars()
#     mol = bam.molecule("GLC") % "14bb"
#     mol *= 4

#     rot = Rotatron(mol.get_atom_graph(), n_processes=4)
#     rot._generate_rotation_unit_masks()
#     rot._find_rotation_units()


# if __name__ == "__main__":
#     import buildamol as bam

#     bam.load_sugars()
#     mol = bam.molecule("GLC") % "14bb"
#     mol *= 4

#     rot = Rotatron(mol.get_atom_graph(), n_processes=4)
#     rot._generate_rotation_unit_masks()
#     rot._find_rotation_units()
