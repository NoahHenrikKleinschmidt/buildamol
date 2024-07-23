import buildamol.base_classes as base_classes


class AtomTyper:
    """
    A base class to assign atom types to atoms in a molecule based on reference data
    """

    def __init__(self, _dict=None):
        self._dict = _dict or {}

    @classmethod
    def from_dict(cls, _dict):
        return cls(_dict)

    @classmethod
    def from_file(cls, filename: str):
        raise NotImplementedError("This method needs to be implemented in a subclass")

    def update(self, _dict_or_typer):
        if isinstance(_dict_or_typer, AtomTyper):
            _dict_or_typer = _dict_or_typer._dict
        self._dict.update(_dict_or_typer)

    def get(self, key):
        return self._dict.get(key, None)

    def atom_key(self, atom: base_classes.Atom):
        """
        The atom key to locate the reference data for a given atom in the dictionary
        Each subclass should implement this method to return an appropriate key for the dictionary
        """
        raise NotImplementedError("This method needs to be implemented in a subclass")

    def set_data(self, atom: base_classes.Atom, data: dict):
        """
        Set the atom type data for a given atom

        Parameters
        ----------
        atom : Atom
            The atom for which the type data is set
        data : dict
            The atom type data
        """
        self._dict[self.atom_key(atom)] = data

    def set_type(self, atom: base_classes.Atom, atom_type: str):
        """
        Set the atom type for a given atom

        Parameters
        ----------
        atom : Atom
            The atom for which the type is set
        atom_type : str
            The atom type
        """
        data = self.get_data(atom)
        data["type"] = atom_type

    def get_data(self, atom: base_classes.Atom) -> dict:
        """
        Get the atom type data for a given atom

        Parameters
        ----------
        atom : Atom
            The atom for which the type data is requested

        Returns
        -------
        dict
            The atom type data
        """
        key = self.atom_key(atom)
        data = self._dict.get(key, None)
        if data is None:
            raise KeyError(f"No data could be found for {atom} (key={key})!")
        return data

    def get_type(self, atom: base_classes.Atom) -> str:
        """
        Get the atom type for a given atom

        Parameters
        ----------
        atom : Atom
            The atom for which the type is requested

        Returns
        -------
        str
            The atom type
        """
        return self.get_data(atom)["type"]

    def get_types(self, residue_or_higher) -> dict:
        """
        Get the atom types for all atoms in a residue

        Parameters
        ----------
        residue_or_higher : Molecule or Residue or Chain or Model or Structure
            Any object that has a `get_atoms` method

        Returns
        -------
        dict
            A dictionary with atoms as keys and atom types as values
        """
        return {atom: self.get_type(atom) for atom in residue_or_higher.get_atoms()}

    def assign_types(self, atom_or_higher):
        """
        Assign a "type" attribute on one or more atoms based on the atom type data

        Parameters
        ----------
        atom_or_higher : Atom or Residue or Chain or Model or Structure
            Any object that is an Atom or has a `get_atoms` method
        """
        if hasattr(atom_or_higher, "get_atoms"):
            for atom in atom_or_higher.get_atoms():
                atom.type = self.get_type(atom)
        else:
            atom_or_higher.type = self.get_type(atom_or_higher)

    def assign_charges(self, atom_or_higher):
        """
        Assign the "pqr_charge" attribute on one or more atoms based on the atom type data

        Parameters
        ----------
        atom_or_higher : Atom or Residue or Chain or Model or Structure
            Any object that is an Atom or has a `get_atoms` method
        """
        if hasattr(atom_or_higher, "get_atoms"):
            for atom in atom_or_higher.get_atoms():
                atom.pqr_charge = self.get_data(atom)["charge"]
        else:
            atom_or_higher.pqr_charge = self.get_data(atom_or_higher)["charge"]

    def assign_attribute(self, atom_or_higher, attribute: str):
        """
        Assign an attribute on one or more atoms based on the atom type data

        Parameters
        ----------
        atom_or_higher : Atom or Residue or Chain or Model or Structure
            Any object that is an Atom or has a `get_atoms` method
        attribute : str
            The attribute to assign
        """
        if hasattr(atom_or_higher, "get_atoms"):
            for atom in atom_or_higher.get_atoms():
                setattr(atom, attribute, self.get_data(atom)[attribute])
        else:
            setattr(atom_or_higher, attribute, self.get_data(atom_or_higher)[attribute])

    def prepare_dict_for(self, atom_or_higher) -> dict:
        """
        Prepare a dictionary with atom types for one or more atoms which can be added to the atom typer dictionary

        Returns
        -------
        dict
        """
        if getattr(atom_or_higher, "get_atoms", None) is None:
            return {
                self.atom_key(atom_or_higher): {
                    "type": getattr(atom_or_higher, "type"),
                    "charge": (getattr(atom_or_higher, "pqr_charge")),
                }
            }
        return {
            atom: {
                self.atom_key(atom): {
                    "type": getattr(atom, "type"),
                    "charge": (getattr(atom, "pqr_charge")),
                }
            }
            for atom in atom_or_higher.get_atoms()
        }
