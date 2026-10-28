"""
Generic in-memory library for collections of ``Molecule`` objects with
optional custom metadata.

``MoleculeLibrary`` fills the gap between raw Python lists and the
PDBE-specific ``PDBECompounds`` class: it works with any set of
``Molecule`` objects and lets you attach arbitrary key/value metadata to
each entry, query by id/name/SMILES or any metadata key, and persist
the whole library to a pickle or portable JSON file.

Typical use cases
-----------------

* Storing a custom fragment library for use with the ``Assembler`` /
  ``Generator`` classes.
* Maintaining a curated set of reference molecules with associated
  experimental or computed properties.
* Bulk-loading / bulk-saving a heterogeneous compound collection that
  was not sourced from PDBE.

Quick start
-----------

.. code-block:: python

    import buildamol as bam
    from buildamol.resources import MoleculeLibrary

    # --- build from a list of molecules --------------------------------
    frags = [bam.Molecule.from_smiles(s, id=name) for s, name in [
        ("c1ccccc1",  "BNZ"),
        ("CCO",       "EOH"),
        ("c1ccncc1",  "PYR"),
    ]]
    lib = MoleculeLibrary.from_molecules(frags)

    # --- attach metadata at creation time ------------------------------
    lib = MoleculeLibrary.from_molecules(
        frags,
        metadata=[
            {"source": "manual", "mw": 78},
            {"source": "manual", "mw": 46},
            {"source": "manual", "mw": 79},
        ],
    )

    # --- add / remove individual entries -------------------------------
    lib.add(bam.Molecule.from_smiles("C1CCNCC1", id="PIP"), mw=85)
    lib.remove("EOH")

    # --- query ----------------------------------------------------------
    mol = lib.get("BNZ")                              # by id (default)
    mol = lib.get("c1ccccc1", by="smiles")            # by canonical SMILES
    mol = lib.get("benzene",  by="name")              # by name metadata field
    hits = lib.get(78, by="mw")                       # by any metadata key

    # --- persist --------------------------------------------------------
    lib.save("my_fragments.pkl")                      # pickle
    lib.to_json("my_fragments.json")                  # portable JSON
    lib2 = MoleculeLibrary.load("my_fragments.pkl")
    lib3 = MoleculeLibrary.from_json("my_fragments.json")

    # --- iterate --------------------------------------------------------
    for mol, meta in lib.iter_items():
        print(mol.id, meta.get("mw"))

    # --- filter to a sub-library ---------------------------------------
    heavy = lib.filter(lambda mol, meta: meta.get("mw", 0) > 70)
"""

import pickle
import warnings

import buildamol.utils.auxiliary as aux


class MoleculeLibrary:
    """
    A generic, metadata-aware collection of ``Molecule`` objects.

    Parameters
    ----------
    id : str, optional
        A name/identifier for this library instance.

    Notes
    -----
    Internally molecules are stored in two parallel dicts keyed by the
    molecule's ``id`` attribute:

    * ``_molecules`` — the ``Molecule`` objects themselves.
    * ``_metadata``  — one plain ``dict`` of user-defined key/value pairs
      per molecule.  The ``"smiles"`` key is populated automatically on
      ``add()`` when the molecule can produce a SMILES string.
    """

    def __init__(self, id: str = None) -> None:
        self.id = id
        self._molecules = {}
        self._metadata = {}
        self._filename = None

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_molecules(
        cls,
        molecules,
        metadata=None,
        id: str = None,
    ) -> "MoleculeLibrary":
        """
        Build a library from an iterable of ``Molecule`` objects.

        Parameters
        ----------
        molecules : iterable of Molecule
            The molecules to add.
        metadata : list of dict, optional
            One metadata dict per molecule (same order as *molecules*).
            Missing entries default to ``{}``.
        id : str, optional
            Library identifier.

        Returns
        -------
        MoleculeLibrary
        """
        lib = cls(id=id)
        if metadata is None:
            metadata = [{}] * len(list(molecules)) if hasattr(molecules, "__len__") else []
        for i, mol in enumerate(molecules):
            meta = metadata[i] if i < len(metadata) else {}
            lib.add(mol, **meta)
        return lib

    @classmethod
    def load(cls, filename: str) -> "MoleculeLibrary":
        """
        Load a ``MoleculeLibrary`` from a pickle file.

        Parameters
        ----------
        filename : str
            Path to a ``.pkl`` file previously created by :meth:`save`.

        Returns
        -------
        MoleculeLibrary
        """
        with open(filename, "rb") as f:
            obj = pickle.load(f)
        if not isinstance(obj, cls):
            raise TypeError(f"File does not contain a MoleculeLibrary (got {type(obj)}).")
        return obj

    @classmethod
    def from_xml(cls, filename: str) -> "MoleculeLibrary":
        """
        Load a ``MoleculeLibrary`` from an XML file written by :meth:`to_xml`.

        Parameters
        ----------
        filename : str

        Returns
        -------
        MoleculeLibrary
        """
        import buildamol.utils.xml as _xml

        root = _xml.read_xml(filename)
        lib = _xml.decode_molecule_library(root)
        lib._filename = filename
        return lib

    @classmethod
    def from_json(cls, filename: str) -> "MoleculeLibrary":
        """
        Load a ``MoleculeLibrary`` from a JSON file previously written by
        :meth:`to_json`.

        Parameters
        ----------
        filename : str
            Path to the JSON file.

        Returns
        -------
        MoleculeLibrary
        """
        import json
        from buildamol.core import Molecule

        with open(filename) as f:
            data = json.load(f)

        lib = cls(id=data.get("id"))
        lib._filename = filename
        for entry in data.get("molecules", []):
            mol_id = entry["id"]
            smiles = entry.get("smiles")
            meta = entry.get("metadata", {})
            if smiles:
                try:
                    mol = Molecule.from_smiles(smiles, id=mol_id)
                except Exception:
                    warnings.warn(f"Could not reconstruct molecule '{mol_id}' from SMILES '{smiles}'. Skipping.")
                    continue
            else:
                warnings.warn(f"Entry '{mol_id}' has no SMILES; cannot reconstruct molecule. Skipping.")
                continue
            lib._molecules[mol_id] = mol
            lib._metadata[mol_id] = {**meta, "smiles": smiles}
        return lib

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, filename: str = None) -> None:
        """
        Save the library to a pickle file.

        Parameters
        ----------
        filename : str, optional
            Destination path.  If omitted, reuses the path from the last
            :meth:`load` or :meth:`save` call.  Raises ``ValueError`` if
            neither is available.
        """
        if filename is None:
            if self._filename is None:
                raise ValueError("No filename specified.")
            filename = self._filename
        if not filename.endswith(".pkl"):
            filename = aux.change_suffix(filename, ".pkl")
        self._filename = filename
        with open(filename, "wb") as f:
            pickle.dump(self, f)

    def to_xml(self, filename: str = None) -> None:
        """
        Export the library to an XML file.

        Each molecule is stored as its canonical SMILES string together
        with the associated metadata.  The XML file can be reloaded with
        :meth:`from_xml` and is version-portable (no pickle dependency).

        Parameters
        ----------
        filename : str, optional
            Destination path.  Defaults to the path used at load time.
        """
        import buildamol.utils.xml as _xml

        if filename is None:
            if self._filename is None:
                raise ValueError("No filename specified.")
            filename = self._filename
        if not filename.endswith(".xml"):
            filename = aux.change_suffix(filename, ".xml")
        self._filename = filename

        root = _xml.encode_molecule_library(self)
        _xml.write_xml(filename, root)

    def to_json(self, filename: str = None) -> None:
        """
        Export the library to a portable JSON file.

        Each molecule is serialised as its canonical SMILES string plus
        the associated metadata dict.  The JSON file can be shared across
        BuildAMol versions without pickle compatibility issues.

        Parameters
        ----------
        filename : str, optional
            Destination path.  Defaults to the path used at load time.
        """
        import json

        if filename is None:
            if self._filename is None:
                raise ValueError("No filename specified.")
            filename = self._filename
        if not filename.endswith(".json"):
            filename = aux.change_suffix(filename, ".json")
        self._filename = filename

        entries = []
        for mol_id, mol in self._molecules.items():
            meta = dict(self._metadata.get(mol_id, {}))
            smiles = meta.pop("smiles", None)
            if smiles is None:
                try:
                    smiles = mol.to_smiles()
                except Exception:
                    smiles = None
            entries.append({"id": mol_id, "smiles": smiles, "metadata": meta})

        with open(filename, "w") as f:
            json.dump({"id": self.id, "molecules": entries}, f, indent=2)

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def add(self, mol, **metadata) -> None:
        """
        Add a molecule to the library.

        Parameters
        ----------
        mol : Molecule
            The molecule to store.  ``mol.id`` is used as the key; if an
            entry with that id already exists it will be overwritten with
            a warning.
        **metadata
            Arbitrary key/value metadata to associate with this molecule.
            ``"smiles"`` is auto-populated if not provided and the molecule
            supports ``to_smiles()``.
        """
        mol_id = mol.id
        if mol_id in self._molecules:
            warnings.warn(f"Molecule '{mol_id}' already present; it will be overwritten.")
        self._molecules[mol_id] = mol
        if "smiles" not in metadata:
            try:
                metadata["smiles"] = mol.to_smiles()
            except Exception:
                pass
        self._metadata[mol_id] = metadata

    def remove(self, id: str) -> None:
        """
        Remove a molecule from the library.

        Parameters
        ----------
        id : str
            The molecule id to remove.  Silently ignored if not found.
        """
        self._molecules.pop(id, None)
        self._metadata.pop(id, None)

    def merge(self, other: "MoleculeLibrary", overwrite: bool = False) -> None:
        """
        Merge another ``MoleculeLibrary`` into this one.

        Parameters
        ----------
        other : MoleculeLibrary
            The library to pull molecules from.
        overwrite : bool
            If ``False`` (default), existing entries are kept and a
            warning is issued for any collision.  If ``True``, colliding
            entries are silently replaced.
        """
        for mol_id, mol in other._molecules.items():
            if mol_id in self._molecules and not overwrite:
                warnings.warn(
                    f"Molecule '{mol_id}' already present; skipping (pass overwrite=True to replace)."
                )
                continue
            self._molecules[mol_id] = mol
            self._metadata[mol_id] = dict(other._metadata.get(mol_id, {}))

    # ------------------------------------------------------------------
    # Querying
    # ------------------------------------------------------------------

    def get(self, query, by: str = "id"):
        """
        Retrieve one or more molecules matching a query.

        Parameters
        ----------
        query
            The value to search for.
        by : str
            What to search against.  Built-in options:

            * ``"id"``    — the molecule's id (exact match, O(1)).
            * ``"name"``  — the ``"name"`` key in the molecule's metadata,
              case-insensitive; also matches the molecule's id field.
            * ``"smiles"`` — the canonical SMILES string stored in metadata.

            Any other string is treated as a metadata key and searched
            for an exact string match of ``str(metadata[by]) == str(query)``.

        Returns
        -------
        Molecule or list of Molecule or None
            A single ``Molecule`` when exactly one match is found,
            a list when multiple matches are found, and ``None`` when
            no match is found.
        """
        ids = self._match(query, by)
        if not ids:
            return None
        mols = [self._molecules[i] for i in ids if i in self._molecules]
        return mols[0] if len(mols) == 1 else mols

    def has(self, query, by: str = "id") -> bool:
        """
        Check whether the library contains a molecule matching *query*.

        Parameters
        ----------
        query
            The value to search for.
        by : str
            Same options as :meth:`get`.

        Returns
        -------
        bool
        """
        return len(self._match(query, by)) > 0

    def filter(self, fn) -> "MoleculeLibrary":
        """
        Return a new ``MoleculeLibrary`` containing only molecules for
        which ``fn(mol, metadata)`` returns a truthy value.

        Parameters
        ----------
        fn : callable
            ``fn(mol: Molecule, metadata: dict) -> bool``

        Returns
        -------
        MoleculeLibrary
        """
        sub = MoleculeLibrary(id=self.id)
        for mol_id, mol in self._molecules.items():
            meta = self._metadata.get(mol_id, {})
            if fn(mol, meta):
                sub._molecules[mol_id] = mol
                sub._metadata[mol_id] = dict(meta)
        return sub

    def get_metadata(self, id: str) -> dict:
        """
        Return the metadata dict for a given molecule id.

        Parameters
        ----------
        id : str

        Returns
        -------
        dict or None
        """
        return self._metadata.get(id)

    def set_metadata(self, id: str, **metadata) -> None:
        """
        Update (merge) metadata for a molecule already in the library.

        Parameters
        ----------
        id : str
            The molecule id.
        **metadata
            Key/value pairs to set or update.
        """
        if id not in self._molecules:
            raise KeyError(f"Molecule '{id}' not in library.")
        self._metadata[id].update(metadata)

    # ------------------------------------------------------------------
    # Iteration / properties
    # ------------------------------------------------------------------

    @property
    def ids(self) -> list:
        """All molecule ids in insertion order."""
        return list(self._molecules.keys())

    @property
    def molecules(self) -> list:
        """All molecules as a list."""
        return list(self._molecules.values())

    def iter_molecules(self):
        """Iterate over all ``Molecule`` objects."""
        return iter(self._molecules.values())

    def iter_items(self):
        """
        Iterate over ``(Molecule, metadata_dict)`` pairs.
        """
        for mol_id, mol in self._molecules.items():
            yield mol, self._metadata.get(mol_id, {})

    # ------------------------------------------------------------------
    # Dunder
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return len(self._molecules)

    def __iter__(self):
        """Yield ``(id, Molecule, metadata)`` triples."""
        for mol_id, mol in self._molecules.items():
            yield mol_id, mol, self._metadata.get(mol_id, {})

    def __getitem__(self, key):
        return self._molecules[key]

    def __contains__(self, key) -> bool:
        return key in self._molecules

    def __repr__(self) -> str:
        name = f" '{self.id}'" if self.id else ""
        return f"MoleculeLibrary{name} ({len(self)} molecules)"

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _match(self, query, by: str) -> list:
        """Return a list of ids that match *query* under the given *by* mode."""
        if by == "id":
            return [str(query)] if str(query) in self._molecules else []
        elif by == "name":
            q = str(query).lower()
            return [
                k for k, meta in self._metadata.items()
                if q == k.lower() or q == str(meta.get("name", "")).lower()
            ]
        elif by == "smiles":
            return [k for k, meta in self._metadata.items() if meta.get("smiles") == str(query)]
        else:
            return [
                k for k, meta in self._metadata.items()
                if str(meta.get(by, "")) == str(query)
            ]


__all__ = ["MoleculeLibrary"]


if __name__ == "__main__":
    import buildamol as bam

    mols = [
        bam.Molecule.from_smiles("c1ccccc1", id="BNZ").autolabel(),
        bam.Molecule.from_smiles("CCO",       id="EOH").autolabel(),
        bam.Molecule.from_smiles("c1ccncc1",  id="PYR").autolabel(),
    ]

    lib = MoleculeLibrary.from_molecules(
        mols,
        metadata=[
            {"name": "benzene",  "mw": 78},
            {"name": "ethanol",  "mw": 46},
            {"name": "pyridine", "mw": 79},
        ],
        id="demo_lib",
    )

    print(lib)
    print("ids:", lib.ids)
    print("get BNZ:", lib.get("BNZ"))
    print("get by name:", lib.get("benzene", by="name"))
    print("get by smiles:", lib.get("c1ccccc1", by="smiles"))
    print("get by mw=78:", lib.get(78, by="mw"))
    print("has EOH:", lib.has("EOH"))
    print("has XYZ:", lib.has("XYZ"))

    heavy = lib.filter(lambda m, meta: meta.get("mw", 0) > 70)
    print("filtered (mw>70):", heavy)

    lib.to_json("/tmp/demo_lib.json")
    lib2 = MoleculeLibrary.from_json("/tmp/demo_lib.json")
    print("round-trip:", lib2)
    for mol, meta in lib2.iter_items():
        print(f"  {mol.id}: {meta}")
