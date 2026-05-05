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
):
    """
    Prepare CHARMM residue/atom names, then perform pure lookup-based atom typing.

    Returns a single atom type for atoms and an atom->type mapping for residues/molecules.
    Types are assigned to atoms in-place when `assign_types=True`.
    """
    if typer is None:
        if filename is None:
            raise ValueError("filename is required when typer is not provided.")
        typer = CHARMMTyper.from_file(filename)

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
        typer.assign_types(atom_or_higher)

    if hasattr(atom_or_higher, "get_atoms"):
        return typer.get_types(atom_or_higher)
    return typer.get_type(atom_or_higher)


class CHARMMTyper(AtomTyper):
    """
    A class to assign CHARMM atom types to atoms in a molecule based on a CHARMM Topology file
    This requires that the molecule has residue names and atom names match with those in the CHARMM Topology file
    """

    def __init__(self, _dict=None):
        super().__init__(_dict=_dict)

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
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("!"):
                    continue
                if line.startswith("MASS"):
                    _, idx, atom_type, mass, *_ = line.split()
                    _type_masses[atom_type] = float(mass)
                    continue
                if line.startswith("RESI"):
                    residue = line.split()[1]
                    continue
                if line.startswith("ATOM"):
                    _, atom_name, atom_type, charge, *_ = line.split()
                    self._dict[f"{residue}:{atom_name}"] = {
                        "type": atom_type,
                        "charge": float(charge),
                    }
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

    def get_data(self, atom) -> dict:
        return super().get_data(atom)

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
