from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Hashable, Iterable, List, Optional, Tuple
import re
import warnings

import networkx as nx

from ._graph_utils import iter_subgraph_candidates


def _normalize_label(label: Optional[str]) -> str:
    if not label:
        return ""
    return re.sub(r"[^A-Za-z0-9]", "", str(label)).upper()


def _element_symbol_from_name(name: str) -> str:
    token = _normalize_label(name)
    if not token:
        return ""
    return token[0]


def _is_hydrogen_element(element: Optional[str]) -> bool:
    return _normalize_label(element) == "H"


@dataclass
class AtomNameTemplate:
    """Stores one atom naming template as a connectivity graph."""

    residue_name: str
    graph: nx.Graph
    atom_data: Dict[str, dict] = field(default_factory=dict)
    source: Optional[str] = None


@dataclass
class AtomNameMatch:
    """Result of matching a target residue graph to a naming template."""

    residue_name: str
    template_name: str
    atom_name_by_node: Dict[Hashable, str]
    score: float


class AtomNameGraphEngine:
    """
    Generic graph-based atom naming engine.

    The engine stores reference templates as graphs and performs graph
    isomorphism against residue-internal connectivity to infer atom names.
    """

    def __init__(
        self,
        allow_partial_matching: bool = True,
        min_node_coverage: float = 0.6,
        min_edge_coverage: float = 0.5,
    ):
        self._templates_by_residue: Dict[str, List[AtomNameTemplate]] = defaultdict(
            list
        )
        self.allow_partial_matching = allow_partial_matching
        self.min_node_coverage = min_node_coverage
        self.min_edge_coverage = min_edge_coverage

    @staticmethod
    def _validate_coverage_thresholds(
        min_node_coverage: float,
        min_edge_coverage: float,
    ):
        if not (0.0 <= min_node_coverage <= 1.0):
            raise ValueError(
                f"min_node_coverage must be in [0, 1], got {min_node_coverage}."
            )
        if not (0.0 <= min_edge_coverage <= 1.0):
            raise ValueError(
                f"min_edge_coverage must be in [0, 1], got {min_edge_coverage}."
            )

    @staticmethod
    def _element_counts_contained(template_graph: nx.Graph, target_graph: nx.Graph) -> bool:
        template_counts = Counter(
            template_graph.nodes[node].get("element", "") for node in template_graph.nodes
        )
        target_counts = Counter(
            target_graph.nodes[node].get("element", "") for node in target_graph.nodes
        )
        return all(target_counts[element] <= template_counts.get(element, 0) for element in target_counts)

    @classmethod
    def load(self, filename: str) -> AtomNameGraphEngine:
        """
        Load a pickled AtomNameEngine from file.
        """
        import pickle

        with open(filename, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, AtomNameGraphEngine):
            raise ValueError(f"File {filename} does not contain an AtomNameEngine")
        return obj

    def save(self, filename: str):
        """
        Save this AtomNameEngine to file using pickle.
        """
        import pickle

        with open(filename, "wb") as f:
            pickle.dump(self, f)

    @property
    def residue_names(self) -> List[str]:
        return sorted(self._templates_by_residue.keys())

    def has_templates(self, residue_name: str) -> bool:
        return residue_name.upper() in self._templates_by_residue

    def register_molecule(self, mol):
        """
        Register template residues from a molecule. The molecule should have residues with correct connectivity.

        Parameters
        ----------
        mol
            A molecule object with residues and connectivity information.
        """
        templates = []
        for residue in mol.get_residues():
            residue_name = residue.resname
            atoms = []
            for atom in residue.get_atoms():
                atom_info = {"name": str(getattr(atom, "id", ""))}
                element = getattr(atom, "element", None)
                if element:
                    atom_info["element"] = element
                atoms.append(atom_info)

            bonds = []
            for bond in residue.get_bonds(residue_internal=True):
                atom1, atom2 = bond
                name1 = str(getattr(atom1, "id", ""))
                name2 = str(getattr(atom2, "id", ""))
                bonds.append((name1, name2))

            template = self.register_template(
                residue_name=residue_name,
                atoms=atoms,
                bonds=bonds,
                source="molecule",
            )
            templates.append(template)
        return templates

    def register_template(
        self,
        residue_name: str,
        atoms: Iterable[dict],
        bonds: Iterable[Tuple[str, str]],
        source: Optional[str] = None,
    ) -> AtomNameTemplate:
        """
        Register one template.

        Parameters
        ----------
        residue_name
            Template residue name.
        atoms
            Iterable of dictionaries with at least keys `name` and optional `element`.
        bonds
            Iterable of 2-tuples with atom names.
        source
            Optional string describing the source of this template (e.g. file name).
        """
        graph = nx.Graph()
        atom_data: Dict[str, dict] = {}

        for atom in atoms:
            atom_name = str(atom["name"])
            element = atom.get("element") or _element_symbol_from_name(atom_name)
            if _is_hydrogen_element(element):
                continue
            graph.add_node(
                atom_name,
                element=_normalize_label(element),
                atom_name=atom_name,
            )
            atom_data[atom_name] = dict(atom)

        for a1, a2 in bonds:
            if a1 not in graph or a2 not in graph:
                continue
            graph.add_edge(a1, a2)

        template = AtomNameTemplate(
            residue_name=residue_name.upper(),
            graph=graph,
            atom_data=atom_data,
            source=source,
        )
        self._templates_by_residue[template.residue_name].append(template)
        return template

    def build_graph_from_residue(self, residue) -> nx.Graph:
        """Build a residue-internal graph using atom elements and current names."""
        graph = nx.Graph()

        atoms = list(residue.get_atoms())
        for atom in atoms:
            element = getattr(atom, "element", None) or _element_symbol_from_name(
                atom.id
            )
            if _is_hydrogen_element(element):
                continue
            graph.add_node(
                atom,
                element=_normalize_label(element),
                current_name=str(atom.id),
            )

        for bond in residue.get_bonds(residue_internal=True):
            atom1, atom2 = bond
            if atom1 in graph and atom2 in graph:
                graph.add_edge(atom1, atom2)

        return graph

    def _element_counts(self, graph: nx.Graph) -> Counter:
        return Counter(graph.nodes[node].get("element", "") for node in graph.nodes)

    def _score_mapping(
        self, target_graph: nx.Graph, mapping: Dict[str, Hashable]
    ) -> float:
        """Heuristic score used to disambiguate symmetric mappings."""
        score = 0.0
        for template_name, node in mapping.items():
            current_name = target_graph.nodes[node].get("current_name", "")
            t_name = _normalize_label(template_name)
            c_name = _normalize_label(current_name)

            if t_name == c_name:
                score += 4.0
            if re.sub(r"\d+$", "", t_name) == re.sub(r"\d+$", "", c_name):
                score += 2.0
            if t_name[:1] and c_name[:1] and t_name[0] == c_name[0]:
                score += 0.5
        return score

    def _match_template(
        self,
        template: AtomNameTemplate,
        target_graph: nx.Graph,
        max_isomorphisms: int = 256,
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Optional[AtomNameMatch]:
        if allow_partial_matching is None:
            allow_partial_matching = self.allow_partial_matching
        if min_node_coverage is None:
            min_node_coverage = self.min_node_coverage
        if min_edge_coverage is None:
            min_edge_coverage = self.min_edge_coverage
        self._validate_coverage_thresholds(min_node_coverage, min_edge_coverage)

        template_nodes = template.graph.number_of_nodes()
        target_nodes = target_graph.number_of_nodes()
        if target_nodes > template_nodes:
            return None

        if not self._element_counts_contained(template.graph, target_graph):
            return None

        if not allow_partial_matching and target_nodes != template_nodes:
            return None

        # Early exit if node coverage is already below threshold before any pruning.
        node_cov_ratio = target_nodes / template_nodes if template_nodes else 1.0
        if node_cov_ratio < min_node_coverage:
            return None

        target_element_counts = Counter(
            target_graph.nodes[n].get("element", "") for n in target_graph.nodes
        )
        template_edge_count = template.graph.number_of_edges()

        best_mapping: Optional[Dict[str, Hashable]] = None
        best_score = float("-inf")

        # Iterate over candidate subgraphs produced by peripheral pruning
        # (leaf-first peeling), then run fast exact isomorphism on each one.
        # This replaces the open-ended subgraph_isomorphisms_iter() call which
        # can enumerate millions of embeddings for realistic residue sizes.
        for candidate, node_cov in iter_subgraph_candidates(
            template.graph, target_nodes, target_element_counts
        ):
            if node_cov < min_node_coverage:
                continue

            edge_cov = (
                candidate.number_of_edges() / template_edge_count
                if template_edge_count
                else 1.0
            )
            if edge_cov < min_edge_coverage:
                continue

            matcher = nx.algorithms.isomorphism.GraphMatcher(
                candidate,
                target_graph,
                node_match=lambda a, b: a.get("element") == b.get("element"),
            )

            for idx, mapping in enumerate(matcher.isomorphisms_iter()):
                score = (
                    self._score_mapping(target_graph, mapping)
                    + node_cov
                    + 0.5 * edge_cov
                )
                if score > best_score:
                    best_score = score
                    best_mapping = mapping
                if idx + 1 >= max_isomorphisms:
                    break

            if best_mapping is not None:
                break  # greedy: first successful candidate wins

        if best_mapping is None:
            return None

        atom_name_by_node = {
            node: template_name for template_name, node in best_mapping.items()
        }
        return AtomNameMatch(
            residue_name=template.residue_name,
            template_name=template.residue_name,
            atom_name_by_node=atom_name_by_node,
            score=best_score,
        )

    def match_residue(
        self,
        residue,
        residue_name: Optional[str] = None,
        fallback_residue_names: Iterable[str] = (),
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Optional[AtomNameMatch]:
        """Infer atom-name mapping for one residue."""
        name_candidates = []
        if residue_name:
            name_candidates.append(residue_name.upper())
        if getattr(residue, "resname", None):
            name_candidates.append(str(residue.resname).upper())
        for candidate in fallback_residue_names:
            if candidate:
                name_candidates.append(str(candidate).upper())

        # Keep order while removing duplicates.
        unique_candidates = list(dict.fromkeys(name_candidates))
        if not unique_candidates:
            return None

        target_graph = self.build_graph_from_residue(residue)
        best_match = None
        for candidate in unique_candidates:
            for template in self._templates_by_residue.get(candidate, []):
                match = self._match_template(
                    template,
                    target_graph,
                    allow_partial_matching=allow_partial_matching,
                    min_node_coverage=min_node_coverage,
                    min_edge_coverage=min_edge_coverage,
                )
                if match is None:
                    continue
                if best_match is None or match.score > best_match.score:
                    best_match = match
        return best_match

    def infer_atom_names(
        self,
        residue,
        residue_name: Optional[str] = None,
        fallback_residue_names: Iterable[str] = (),
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Dict[Hashable, str]:
        match = self.match_residue(
            residue,
            residue_name=residue_name,
            fallback_residue_names=fallback_residue_names,
            allow_partial_matching=allow_partial_matching,
            min_node_coverage=min_node_coverage,
            min_edge_coverage=min_edge_coverage,
        )
        if match is None:
            return {}
        return match.atom_name_by_node

    def rename(
        self,
        residue,
        residue_name: Optional[str] = None,
        fallback_residue_names: Iterable[str] = (),
        on_missing: str = "ignore",
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ):
        """
        Rename atoms of a residue in-place.

        Parameters
        ----------
        residue
            Residue-like object with `get_atoms` and `get_bonds`.
        residue_name
            Optional residue name override for matching.
        fallback_residue_names
            Additional residue-name candidates used during template lookup.
        on_missing
            Behavior when no naming template match is found.
            One of: `"ignore"`, `"warn"`, `"error"`.
        """
        if on_missing not in ("ignore", "warn", "error"):
            raise ValueError(
                f"Invalid on_missing={on_missing!r}; expected 'ignore', 'warn', or 'error'."
            )

        if hasattr(residue, "get_residues"):
            for _residue in residue.get_residues():
                self.rename(
                    _residue,
                    residue_name=residue_name,
                    fallback_residue_names=fallback_residue_names,
                    on_missing=on_missing,
                    allow_partial_matching=allow_partial_matching,
                    min_node_coverage=min_node_coverage,
                    min_edge_coverage=min_edge_coverage,
                )
            return residue

        atom_name_by_node = self.infer_atom_names(
            residue,
            residue_name=residue_name,
            fallback_residue_names=fallback_residue_names,
            allow_partial_matching=allow_partial_matching,
            min_node_coverage=min_node_coverage,
            min_edge_coverage=min_edge_coverage,
        )

        if len(atom_name_by_node) == 0:
            message = (
                "No atom naming template could be matched for residue "
                f"'{getattr(residue, 'resname', '<unknown>')}'."
            )
            if on_missing == "warn":
                warnings.warn(message, UserWarning)
            elif on_missing == "error":
                raise KeyError(message)
            return residue

        for node, new_name in atom_name_by_node.items():
            setattr(node, "id", new_name)
            if hasattr(node, "name"):
                setattr(node, "name", new_name)
        return residue

    def infer_residue_atom_name_map(
        self,
        residue,
        residue_name: Optional[str] = None,
        fallback_residue_names: Iterable[str] = (),
    ) -> Dict[str, str]:
        """Return current_name -> inferred_template_name for a residue."""
        by_node = self.infer_atom_names(
            residue,
            residue_name=residue_name,
            fallback_residue_names=fallback_residue_names,
        )
        out = {}
        for node, inferred in by_node.items():
            out[str(getattr(node, "id", ""))] = inferred
        return out
