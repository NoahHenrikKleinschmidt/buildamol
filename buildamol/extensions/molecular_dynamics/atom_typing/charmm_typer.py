import warnings
from typing import Optional

from buildamol.extensions.molecular_dynamics.atom_typing import AtomTyper
from buildamol.extensions.naming.charmm import (
    CHARMMAtomNameGraphEngine,
    CHARMMResidueNameLookupEngine,
)


def _iter_residues(atom_or_higher):
    if hasattr(atom_or_higher, "get_residues"):
        yield from atom_or_higher.get_residues()
        return
    if hasattr(atom_or_higher, "get_atoms"):
        yield atom_or_higher
        return
    residue = getattr(atom_or_higher, "parent", None)
    if residue is not None:
        yield residue


def _handle_missing_residue_name(residue, on_missing: str):
    message = (
        "No CHARMM residue-name mapping could be found for residue "
        f"'{getattr(residue, 'resname', '<unknown>')}'."
    )
    if on_missing == "warn":
        warnings.warn(message, UserWarning)
    elif on_missing == "error":
        raise KeyError(message)


def rename_for_charmm_typing(
    atom_or_higher,
    typer: "CHARMMTyper",
    residue_name_engine: Optional[CHARMMResidueNameLookupEngine] = None,
    atom_name_engine: Optional[CHARMMAtomNameGraphEngine] = None,
    filename: Optional[str] = None,
    rename_residues: bool = True,
    rename_atoms: bool = True,
    on_missing: str = "error",
):
    """
    Rename residues and atoms in-place so a pure lookup-only CHARMMTyper can type them.

    Parameters
    ----------
    atom_or_higher
        Atom, residue, molecule, or any object exposing `get_atoms`/`get_residues`.
    typer
        A lookup-only CHARMMTyper.
    residue_name_engine
        Optional prebuilt residue-name translator. Defaults to CHARMMResidueNameLookupEngine.
    atom_name_engine
        Optional prebuilt atom-name graph engine. If omitted, `filename` is used to build one.
    filename
        CHARMM RTF file used to build `atom_name_engine` when not provided.
    rename_residues
        Whether to normalize residue names to CHARMM names before atom renaming/typing.
    rename_atoms
        Whether to normalize atom names using the CHARMM graph naming engine.
    on_missing
        Missing-name behavior: `ignore`, `warn`, or `error`.
    """
    if on_missing not in ("ignore", "warn", "error"):
        raise ValueError(
            f"Invalid on_missing={on_missing!r}; expected 'ignore', 'warn', or 'error'."
        )

    if rename_residues:
        residue_name_engine = residue_name_engine or CHARMMResidueNameLookupEngine()
        valid_names = set(typer.residue_names)
        for residue in _iter_residues(atom_or_higher):
            current_name = str(getattr(residue, "resname", "")).upper()
            if current_name in valid_names:
                continue

            translated = residue_name_engine.translate(current_name)
            if translated is not None:
                residue.resname = translated
                continue

            if on_missing != "ignore":
                _handle_missing_residue_name(residue, on_missing)

    if rename_atoms:
        if atom_name_engine is None:
            if filename is None:
                raise ValueError(
                    "filename is required when atom_name_engine is not provided."
                )
            atom_name_engine = CHARMMAtomNameGraphEngine.from_file(
                filename,
                residue_whitelist=typer.residue_names,
            )
        atom_name_engine.rename(atom_or_higher, on_missing=on_missing)

    if hasattr(atom_or_higher, "relabel_hydrogens"):
        atom_or_higher.relabel_hydrogens()
    return atom_or_higher


def type_with_charmm(
    atom_or_higher,
    filename: Optional[str] = None,
    typer: Optional["CHARMMTyper"] = None,
    residue_name_engine: Optional[CHARMMResidueNameLookupEngine] = None,
    atom_name_engine: Optional[CHARMMAtomNameGraphEngine] = None,
    rename_residues: bool = True,
    rename_atoms: bool = True,
    on_missing: str = "error",
    assign_types: bool = True,
    keep_existing_type: bool = False,
    keep_existing_charge: bool = False,
    keep_existing: Optional[bool] = None,
):
    """
    Prepare CHARMM residue/atom names, then perform pure lookup-based atom typing.

    Returns a single atom type for atoms and an atom->type mapping for residues/molecules.
    Types are assigned to atoms in-place when `assign_types=True`.

    Parameters (additional)
    -----------------------
    keep_existing_type : bool
        If True, atom types already present on atoms take precedence over CHARMM lookup.
    keep_existing_charge : bool
        If True, atom charges already present on atoms take precedence over CHARMM lookup.
    keep_existing : bool, optional
        Backward-compatible alias that enables both ``keep_existing_type`` and
        ``keep_existing_charge`` when True.
    """
    if typer is None:
        if filename is None:
            raise ValueError("filename is required when typer is not provided.")
        typer = CHARMMTyper.from_file(filename)

    if keep_existing is not None:
        keep_existing_type = keep_existing_type or keep_existing
        keep_existing_charge = keep_existing_charge or keep_existing

    rename_for_charmm_typing(
        atom_or_higher,
        typer=typer,
        residue_name_engine=residue_name_engine,
        atom_name_engine=atom_name_engine,
        filename=filename,
        rename_residues=rename_residues,
        rename_atoms=rename_atoms,
        on_missing=on_missing,
    )

    if assign_types:
        typer.assign_types(
            atom_or_higher,
            keep_existing_type=keep_existing_type,
        )

    if hasattr(atom_or_higher, "get_atoms"):
        return typer.get_types(
            atom_or_higher,
            keep_existing_type=keep_existing_type,
        )
    return typer.get_type(
        atom_or_higher,
        keep_existing_type=keep_existing_type,
    )


class CHARMMTyper(AtomTyper):
    """
    A class to assign CHARMM atom types to atoms in a molecule based on a CHARMM Topology file
    This requires that the molecule has residue names and atom names match with those in the CHARMM Topology file
    """

    def __init__(self, _dict=None, _pres_dict=None):
        super().__init__(_dict=_dict)
        self._pres_dict = _pres_dict or {}
        self._pres_fallback_reported = set()
        self._keep_existing_reported = set()

    def _report_keep_existing(self, atom, key: str, attribute: str, value):
        sig = (id(atom), key, attribute)
        if sig in self._keep_existing_reported:
            return
        if attribute == "type":
            print(
                "[INFO] CHARMMTyper: keeping existing type for "
                f"{key}; using atom.type='{value}'."
            )
        elif attribute == "charge":
            print(
                "[INFO] CHARMMTyper: keeping existing charge for "
                f"{key}; using atom.pqr_charge={value}."
            )
        self._keep_existing_reported.add(sig)

    def atom_key(self, atom):
        return f"{atom.parent.resname}:{atom.id}"

    def read(self, filename: str):
        """
        Read atom type data from a CHARMM Topology file and add to this instance.

        Parameters
        ----------
        filename : str
            The filename of the CHARMM Topology file

        Returns
        -------
        self
            Returns self for method chaining
        """
        _type_masses = {}
        residue = None
        patch = None
        with open(filename, "r") as f:
            for line in f:
                line = line.split("!", 1)[0].strip()
                if not line:
                    continue
                if line.startswith("MASS"):
                    _, idx, atom_type, mass, *_ = line.split()
                    _type_masses[atom_type] = float(mass)
                    continue
                if line.startswith("RESI"):
                    residue = line.split()[1]
                    patch = None
                    continue
                if line.startswith("PRES"):
                    patch = line.split()[1]
                    residue = None
                    continue
                if line.startswith("ATOM"):
                    _, atom_name, atom_type, charge, *_ = line.split()
                    data = {
                        "type": atom_type,
                        "charge": float(charge),
                    }
                    if residue is not None:
                        self._dict[f"{residue}:{atom_name}"] = data
                    elif patch is not None and atom_name not in self._pres_dict:
                        # PRES entries are stored without residue key for soft fallback lookup.
                        self._pres_dict[atom_name] = data
        return self

    @classmethod
    def from_file(
        cls,
        filename,
    ):
        """
        Create a CHARMMTyper object from one or more CHARMM Topology files

        Parameters
        ----------
        filename : str or list of str
            The filename(s) of the CHARMM Topology file(s)

        Returns
        -------
        CHARMMTyper
            The CHARMMTyper object
        """
        typer = cls()

        # Handle both single file and list of files
        filenames = filename if isinstance(filename, (list, tuple)) else [filename]

        for fname in filenames:
            typer.read(fname)

        return typer

    def get_data(
        self,
        atom,
        keep_existing_type: bool = False,
        keep_existing_charge: bool = False,
        keep_existing: Optional[bool] = None,
        requested_attribute: Optional[str] = None,
    ) -> dict:
        """
        Get the atom type data for a given atom.

        Parameters
        ----------
        atom : Atom
            The atom for which the type data is requested.
        keep_existing_type : bool
            If True and atom.type is set, prefer it over lookup data.
        keep_existing_charge : bool
            If True and atom.pqr_charge is set, prefer it over lookup data.
        keep_existing : bool, optional
            Backward-compatible alias that enables both keep-existing flags.
        requested_attribute : str, optional
            Internal hint ("type" or "charge") used to make lookup diagnostics
            more explicit.
        """
        if keep_existing is not None:
            keep_existing_type = keep_existing_type or keep_existing
            keep_existing_charge = keep_existing_charge or keep_existing

        existing_type = getattr(atom, "type", None) if keep_existing_type else None
        existing_charge = (
            getattr(atom, "pqr_charge", None) if keep_existing_charge else None
        )

        key = self.atom_key(atom)

        if existing_type is not None:
            self._report_keep_existing(atom, key, "type", existing_type)

        if existing_charge is not None:
            self._report_keep_existing(atom, key, "charge", existing_charge)

        data = None
        if existing_type is None or existing_charge is None:
            data = self._dict.get(key, None)
            if data is None:
                data = self._pres_dict.get(atom.id, None)
            if data is not None and key not in self._dict:
                # Print info once per atom object to keep feedback useful but not too noisy.
                context = requested_attribute or "type/charge"
                atom_sig = (id(atom), key, context)
                if atom_sig not in self._pres_fallback_reported:
                    print(
                        "[INFO] CHARMMTyper: standard lookup failed for "
                        f"{key}; using PRES fallback for {context} of atom '{atom.id}'."
                    )
                    self._pres_fallback_reported.add(atom_sig)

        if data is None and existing_type is None and existing_charge is None:
            raise KeyError(f"No data could be found for {atom} (key={key})!")

        return {
            "type": (
                existing_type if existing_type is not None else (data or {}).get("type")
            ),
            "charge": (
                existing_charge
                if existing_charge is not None
                else (data or {}).get("charge")
            ),
        }

    def get_type(
        self,
        atom,
        keep_existing_type: bool = False,
        keep_existing: Optional[bool] = None,
    ) -> str:
        if keep_existing is not None:
            keep_existing_type = keep_existing_type or keep_existing

        key = self.atom_key(atom)
        if keep_existing_type:
            existing_type = getattr(atom, "type", None)
            if existing_type is not None:
                self._report_keep_existing(atom, key, "type", existing_type)
                return existing_type

        data = self.get_data(
            atom,
            keep_existing_type=keep_existing_type,
            requested_attribute="type",
        )
        if data["type"] is None:
            raise KeyError(
                f"No type could be found for {atom} (key={self.atom_key(atom)})!"
            )
        return data["type"]

    def get_charge(
        self,
        atom,
        keep_existing_charge: bool = False,
        keep_existing: Optional[bool] = None,
    ) -> float:
        if keep_existing is not None:
            keep_existing_charge = keep_existing_charge or keep_existing

        key = self.atom_key(atom)
        if keep_existing_charge:
            existing_charge = getattr(atom, "pqr_charge", None)
            if existing_charge is not None:
                self._report_keep_existing(atom, key, "charge", existing_charge)
                return existing_charge

        data = self.get_data(
            atom,
            keep_existing_charge=keep_existing_charge,
            requested_attribute="charge",
        )
        if data["charge"] is None:
            raise KeyError(
                f"No charge could be found for {atom} (key={self.atom_key(atom)})!"
            )
        return data["charge"]

    def get_types(
        self,
        residue_or_higher,
        keep_existing_type: bool = False,
        keep_existing: Optional[bool] = None,
    ) -> dict:
        if keep_existing is not None:
            keep_existing_type = keep_existing_type or keep_existing

        return {
            atom: self.get_type(
                atom,
                keep_existing_type=keep_existing_type,
            )
            for atom in residue_or_higher.get_atoms()
        }

    def assign_types(
        self,
        atom_or_higher,
        keep_existing_type: bool = False,
        keep_existing: Optional[bool] = None,
    ):
        """
        Assign a ``type`` attribute on one or more atoms.

        Parameters
        ----------
        atom_or_higher : Atom or Residue or Molecule or ...
            Any atom or object with a ``get_atoms`` method.
        keep_existing_type : bool
            If True, atoms that already have a type keep that value.
        keep_existing : bool, optional
            Backward-compatible alias for ``keep_existing_type``.
        """
        if keep_existing is not None:
            keep_existing_type = keep_existing_type or keep_existing

        if hasattr(atom_or_higher, "get_atoms"):
            for atom in atom_or_higher.get_atoms():
                atom.type = self.get_type(
                    atom,
                    keep_existing_type=keep_existing_type,
                )
        else:
            atom_or_higher.type = self.get_type(
                atom_or_higher,
                keep_existing_type=keep_existing_type,
            )

    def assign_charges(
        self,
        atom_or_higher,
        keep_existing_charge: bool = False,
        keep_existing: Optional[bool] = None,
    ):
        if keep_existing is not None:
            keep_existing_charge = keep_existing_charge or keep_existing

        if hasattr(atom_or_higher, "get_atoms"):
            for atom in atom_or_higher.get_atoms():
                atom.pqr_charge = self.get_charge(
                    atom,
                    keep_existing_charge=keep_existing_charge,
                )
        else:
            atom_or_higher.pqr_charge = self.get_charge(
                atom_or_higher,
                keep_existing_charge=keep_existing_charge,
            )

    @property
    def residue_names(self) -> list:
        """
        Get the residue names in the CHARMM Topology file

        Returns
        -------
        list
            The residue names
        """
        return list(sorted(set([key.split(":")[0] for key in self._dict.keys()])))


if __name__ == "__main__":
    import buildamol as bam

    bam.load_sugars()
    man = bam.get_compound("MAN")
    typer = CHARMMTyper.from_file(
        "/Users/noahhk/GIT/glycosylator/support/toppar_charmm/carbohydrates.rtf"
    )

    out = typer.get_types(man)
    print(out)
    print(typer._dict)
