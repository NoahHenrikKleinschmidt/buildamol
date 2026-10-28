from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
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
class ResidueNameTemplate:
    residue_name: str
    graph: nx.Graph
    source: Optional[str] = None


@dataclass
class ResidueNameMatch:
    input_name: Optional[str]
    matched_name: str
    score: float


class ResidueNameLookupEngine:
    """Simple dictionary-based residue naming engine."""

    def __init__(self, mapping: Dict[str, object], value_key: Optional[str] = None):
        self._mapping = {}
        for key, value in mapping.items():
            if value_key and isinstance(value, dict):
                value = value.get(value_key)
            self._mapping[str(key).upper()] = value

    def translate(
        self,
        residue_name: str,
        explicit_form: bool = False,
        default=None,
    ):
        value = self._mapping.get(str(residue_name).upper())
        if value is None:
            return default
        if explicit_form:
            return value if isinstance(value, list) else [value]
        return value[0] if isinstance(value, list) else value

    def candidates(self, residue_name: str) -> List[str]:
        names = self.translate(residue_name, explicit_form=True, default=[])
        return [str(i) for i in names]

    def bulk_translate(
        self,
        residue_names: Iterable[str],
        keep_unknown: bool = True,
    ) -> List[Optional[str]]:
        out = []
        for name in residue_names:
            mapped = self.translate(name)
            if mapped is None and keep_unknown:
                out.append(name)
            else:
                out.append(mapped)
        return out

    def rename(self, mol, on_missing: str = "ignore"):
        """
        Rename residues in a molecule according to the mapping.

        Parameters
        ----------
        mol
            Molecule-like object exposing `get_residues`.
        on_missing
            Behavior when no mapping is available for a residue name.
            One of: `"ignore"`, `"warn"`, `"error"`.
        """
        if on_missing not in ("ignore", "warn", "error"):
            raise ValueError(
                f"Invalid on_missing={on_missing!r}; expected 'ignore', 'warn', or 'error'."
            )

        for residue in mol.get_residues():
            new_name = self.translate(residue.resname)
            if new_name is not None:
                residue.resname = new_name
                continue

            message = (
                "No residue-name mapping could be found for residue "
                f"'{getattr(residue, 'resname', '<unknown>')}'."
            )
            if on_missing == "warn":
                warnings.warn(message, UserWarning)
            elif on_missing == "error":
                raise KeyError(message)
        return mol

    def __call__(self, *args, **kwds):
        return self.rename(*args, **kwds)


class ResidueNameGraphEngine:
    """
    Graph-based residue naming engine.

    Templates are stored as residue-internal connectivity graphs and matched by
    graph isomorphism and element labels.
    """

    def __init__(
        self,
        allow_partial_matching: bool = True,
        min_node_coverage: float = 0.6,
        min_edge_coverage: float = 0.5,
    ):
        self._templates_by_name: Dict[str, List[ResidueNameTemplate]] = defaultdict(
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
    def _element_counts_contained(
        template_graph: nx.Graph, target_graph: nx.Graph
    ) -> bool:
        template_counts = Counter(
            template_graph.nodes[node].get("element", "")
            for node in template_graph.nodes
        )
        target_counts = Counter(
            target_graph.nodes[node].get("element", "") for node in target_graph.nodes
        )
        return all(
            target_counts[element] <= template_counts.get(element, 0)
            for element in target_counts
        )

    @classmethod
    def from_molecule(cls, mol, **kwargs) -> "ResidueNameGraphEngine":
        engine = cls(**kwargs)
        engine.register_molecule(mol)
        return engine

    @classmethod
    def from_molecules(cls, molecules: Iterable, **kwargs) -> "ResidueNameGraphEngine":
        engine = cls(**kwargs)
        for mol in molecules:
            engine.register_molecule(mol)
        return engine

    @classmethod
    def from_compounds(
        cls,
        compounds,
        compound_ids: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> "ResidueNameGraphEngine":
        engine = cls(**kwargs)
        if compound_ids is None:
            molecules = compounds.iter_molecules()
        else:
            molecules = (
                compounds.get(compound_id, by="id", return_type="molecule")
                for compound_id in compound_ids
            )

        for mol in molecules:
            if mol is None:
                continue
            engine.register_molecule(mol)
        return engine

    @classmethod
    def from_json(
        cls,
        filename: str,
        compound_ids: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> "ResidueNameGraphEngine":
        from buildamol.resources.pdbe_compounds import PDBECompounds

        return cls.from_compounds(
            PDBECompounds.from_json(filename),
            compound_ids=compound_ids,
            **kwargs,
        )

    @classmethod
    def from_xml(
        cls,
        filename: str,
        compound_ids: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> "ResidueNameGraphEngine":
        from buildamol.resources.pdbe_compounds import PDBECompounds

        return cls.from_compounds(
            PDBECompounds.from_xml(filename),
            compound_ids=compound_ids,
            **kwargs,
        )

    @classmethod
    def from_pickle(
        cls,
        filename: str,
        compound_ids: Optional[Iterable[str]] = None,
        **kwargs,
    ) -> "ResidueNameGraphEngine":
        import pickle

        with open(filename, "rb") as handle:
            obj = pickle.load(handle)

        if isinstance(obj, cls):
            return obj
        if obj.__class__.__name__ == "PDBECompounds":
            return cls.from_compounds(obj, compound_ids=compound_ids, **kwargs)
        raise ValueError(
            f"File {filename} does not contain an {cls.__name__} or PDBECompounds object"
        )

    @classmethod
    def load(cls, filename: str) -> "ResidueNameGraphEngine":
        """Load a pickled ResidueNameGraphEngine from file."""
        import pickle

        with open(filename, "rb") as handle:
            obj = pickle.load(handle)
        if not isinstance(obj, cls):
            raise ValueError(f"File {filename} does not contain a {cls.__name__}")
        return obj

    def save(self, filename: str):
        """Save this ResidueNameGraphEngine to file using pickle."""
        import pickle

        with open(filename, "wb") as handle:
            pickle.dump(self, handle)

    @property
    def residue_names(self) -> List[str]:
        return sorted(self._templates_by_name.keys())

    def register_template(
        self,
        residue_name: str,
        atoms: Iterable[dict],
        bonds: Iterable[Tuple[str, str]],
        source: Optional[str] = None,
    ) -> ResidueNameTemplate:
        graph = nx.Graph()
        for atom in atoms:
            atom_name = str(atom["name"])
            element = atom.get("element") or _element_symbol_from_name(atom_name)
            if _is_hydrogen_element(element):
                continue
            graph.add_node(atom_name, element=_normalize_label(element))

        for atom1, atom2 in bonds:
            if atom1 in graph and atom2 in graph:
                graph.add_edge(atom1, atom2)

        template = ResidueNameTemplate(
            residue_name=str(residue_name).upper(),
            graph=graph,
            source=source,
        )
        self._templates_by_name[template.residue_name].append(template)
        return template

    def register_molecule(self, mol):
        for residue in mol.get_residues():
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
                bonds.append(
                    (str(getattr(atom1, "id", "")), str(getattr(atom2, "id", "")))
                )

            self.register_template(residue.resname, atoms, bonds, source="molecule")

    def build_graph_from_residue(self, residue) -> nx.Graph:
        graph = nx.Graph()
        for atom in residue.get_atoms():
            element = getattr(atom, "element", None) or _element_symbol_from_name(
                atom.id
            )
            if _is_hydrogen_element(element):
                continue
            graph.add_node(atom, element=_normalize_label(element))

        for bond in residue.get_bonds(residue_internal=True):
            atom1, atom2 = bond
            if atom1 in graph and atom2 in graph:
                graph.add_edge(atom1, atom2)

        return graph

    def _element_counts(self, graph: nx.Graph) -> Counter:
        return Counter(graph.nodes[node].get("element", "") for node in graph.nodes)

    def _score(
        self,
        template_name: str,
        residue,
        preferred_names: Iterable[str],
    ) -> float:
        score = 0.0
        current_name = str(getattr(residue, "resname", "")).upper()
        if template_name == current_name:
            score += 3.0

        preferred = [str(i).upper() for i in preferred_names if i]
        if template_name in preferred:
            score += 2.0

        if current_name and template_name[:1] == current_name[:1]:
            score += 0.1
        return score

    def match_residue(
        self,
        residue,
        preferred_names: Iterable[str] = (),
        allowed_names: Optional[Iterable[str]] = None,
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Optional[ResidueNameMatch]:
        if allow_partial_matching is None:
            allow_partial_matching = self.allow_partial_matching
        if min_node_coverage is None:
            min_node_coverage = self.min_node_coverage
        if min_edge_coverage is None:
            min_edge_coverage = self.min_edge_coverage
        self._validate_coverage_thresholds(min_node_coverage, min_edge_coverage)

        target = self.build_graph_from_residue(residue)
        allowed = None
        if allowed_names is not None:
            allowed = {str(i).upper() for i in allowed_names}

        best = None
        for name, templates in self._templates_by_name.items():
            if allowed is not None and name not in allowed:
                continue
            for template in templates:
                template_nodes = template.graph.number_of_nodes()
                target_nodes = target.number_of_nodes()
                if target_nodes > template_nodes:
                    continue
                if not self._element_counts_contained(template.graph, target):
                    continue

                if not allow_partial_matching and target_nodes != template_nodes:
                    continue

                node_cov_ratio = (
                    target_nodes / template_nodes if template_nodes else 1.0
                )
                if node_cov_ratio < min_node_coverage:
                    continue

                target_element_counts = Counter(
                    target.nodes[v].get("element", "") for v in target.nodes
                )
                template_edge_count = template.graph.number_of_edges()

                # Use peripheral pruning + exact isomorphism instead of
                # subgraph_isomorphisms_iter() which enumerates all embeddings.
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
                        target,
                        node_match=lambda a, b: a.get("element") == b.get("element"),
                    )

                    # is_isomorphic() is sufficient — we only need the name, not a mapping.
                    if matcher.is_isomorphic():
                        score = (
                            self._score(name, residue, preferred_names)
                            + node_cov
                            + 0.5 * edge_cov
                        )
                        if best is None or score > best.score:
                            best = ResidueNameMatch(
                                input_name=getattr(residue, "resname", None),
                                matched_name=name,
                                score=score,
                            )
                        break  # greedy: first matching candidate for this template wins

        return best

    def infer_residue_name(
        self,
        residue,
        preferred_names: Iterable[str] = (),
        allowed_names: Optional[Iterable[str]] = None,
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Optional[str]:
        match = self.match_residue(
            residue,
            preferred_names=preferred_names,
            allowed_names=allowed_names,
            allow_partial_matching=allow_partial_matching,
            min_node_coverage=min_node_coverage,
            min_edge_coverage=min_edge_coverage,
        )
        if match is None:
            return None
        return match.matched_name

    def rename_residue(
        self,
        residue,
        preferred_names: Iterable[str] = (),
        allowed_names: Optional[Iterable[str]] = None,
        on_missing: str = "ignore",
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ) -> Optional[str]:
        if on_missing not in ("ignore", "warn", "error"):
            raise ValueError(
                f"Invalid on_missing={on_missing!r}; expected 'ignore', 'warn', or 'error'."
            )

        new_name = self.infer_residue_name(
            residue,
            preferred_names=preferred_names,
            allowed_names=allowed_names,
            allow_partial_matching=allow_partial_matching,
            min_node_coverage=min_node_coverage,
            min_edge_coverage=min_edge_coverage,
        )
        if new_name is not None:
            residue.resname = new_name
            return new_name

        message = (
            "No residue naming template could be matched for residue "
            f"'{getattr(residue, 'resname', '<unknown>')}'."
        )
        if on_missing == "warn":
            warnings.warn(message, UserWarning)
        elif on_missing == "error":
            raise KeyError(message)
        return new_name

    def rename(
        self,
        mol,
        preferred_names: Iterable[str] = (),
        allowed_names: Optional[Iterable[str]] = None,
        on_missing: str = "ignore",
        allow_partial_matching: Optional[bool] = None,
        min_node_coverage: Optional[float] = None,
        min_edge_coverage: Optional[float] = None,
    ):
        for residue in mol.get_residues():
            self.rename_residue(
                residue,
                preferred_names=preferred_names,
                allowed_names=allowed_names,
                on_missing=on_missing,
                allow_partial_matching=allow_partial_matching,
                min_node_coverage=min_node_coverage,
                min_edge_coverage=min_edge_coverage,
            )
        return mol

    def __call__(self, *args, **kwds):
        return self.rename(*args, **kwds)
