from buildamol.extensions.molecular_dynamics.atom_typing import AtomTyper


class CHARMMTyper(AtomTyper):
    """
    A class to assign CHARMM atom types to atoms in a molecule based on a CHARMM Topology file
    This requires that the molecule has residue names and atom names match with those in the CHARMM Topology file
    """

    def atom_key(self, atom):
        return f"{atom.parent.resname}:{atom.id}"

    @classmethod
    def from_file(cls, filename: str):
        """
        Create a CHARMMTyper object from a CHARMM Topology file

        Parameters
        ----------
        filename : str
            The filename of the CHARMM Topology file

        Returns
        -------
        CHARMMTyper
            The CHARMMTyper object
        """
        typer = cls()
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
                    typer._dict[f"{residue}:{atom_name}"] = {
                        "type": atom_type,
                        "charge": float(charge),
                    }
        return typer

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
