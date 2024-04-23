from typing import Union
import buildamol.core as core
import buildamol.structural as structural
import buildamol.optimizers as optimizers

_geometry_hydrogen_mapping = {
    "Octahedral": {
        "planar": ["H1", "H2", "H3", "H4"],
        "axial": ["H5", "H6"],
    },
    "Tetrahedral": {
        "planar": ["H1", "H2", "H3"],
        "axial": ["H4"],
    },
    "TrigonalBipyramidal": {
        "planar": ["H1", "H2", "H3"],
        "axial": ["H4", "H5"],
    },
    "TrigonalPlanar": {
        "planar": ["H1", "H2", "H3"],
        "axial": ["H1", "H2", "H3"],
    },
    "Linear": {
        "planar": ["H1", "H2"],
        "axial": ["H1", "H2"],
    },
}


class MetalComplexer:
    """
    The MetalComplexer can be used to make metal complexes by
    aligning ligands to a metal core.

    Parameters
    ----------
    metal : str or Atom
        The metal atom or its element symbol.
    geometry : Geometry
        The geometry of the metal complex.
    """

    def __init__(
        self, metal: Union[str, "core.Atom"], geometry: "structural.geometry.Geometry"
    ):
        self.metal = core.Atom.new(metal) if isinstance(metal, str) else metal
        self.geometry = geometry
        self._core = None
        self._backup = None
        self._hydrogen_mapping = dict(
            **_geometry_hydrogen_mapping[geometry.__class__.__name__]
        )
        self._bond_length = 2.0

    def store(self):
        """
        Store the current state of the complexer.
        """
        self._backup = self._core.copy()

    def reset(self):
        """
        Reset the complexer to its stored state.
        """
        self._core = self._backup.copy()

    def get_complex(self) -> "core.Molecule":
        """
        Get the current metal complex.

        Returns
        -------
        Molecule
            The metal complex.
        """
        return self._core

    def make_core(self, bond_length: float = 2.0) -> "core.Molecule":
        """
        Make the metal core.

        Parameters
        ----------
        bond_length : float, optional
            The bond length of the metal core, by default 2.0

        Returns
        -------
        Molecule
            The metal core.
        """
        mol = core.Molecule.from_geometry(
            self.geometry, atoms=[self.metal], id="complex", resname="MOL"
        )
        for bond in mol.get_bonds():
            mol.adjust_bond_length(*bond, bond_length)

        self._core = mol
        self._backup = mol.copy()
        self._bond_length = bond_length
        return mol

    def add_ligand(
        self,
        ligand: "core.Molecule",
        binders: list,
        delete: list = None,
        direction: str = "planar",
        acceptors: list = None,
        optimize: bool = None,
        optimize_kwargs: dict = None,
        copy: bool = True,
    ) -> "core.Molecule":
        """
        Add a ligand to the metal core.

        Parameters
        ----------
        ligand : Molecule
            The ligand to add.
        binders : list
            A list of atom objects or other atom identifiers to get atoms from the ligand,
            that will bind to the metal core.
        delete: list, optional
            A list of atom objects or other atom identifiers to get atoms from the ligand,
            that will be deleted after alignment. Each of these must be a direct neighbor of a binder atom.
            They must be specified in the same order as the binders, by default None, in which case a random hydrogen neighbor from each binder is deleted.
        direction : str
            The direction in which to align the ligand to the metal core, if multiple binders are specified. This can be either
            "planar" or "axial". To achieve a mixed alignment or specific alignment to hydrogens use the "acceptors" argument instead.
        acceptors : list, optional
            A list of atom objects or other atom identifiers to get hydrogen atoms from the core. These will be used to align the ligand to the core.
            If specified, the "direction" argument will be ignored, by default None. The acceptors must be specified in the same order as the binders to which they should be aligned.
        optimize : bool, optional
            If True, optimize the ligand after alignment, by default optimization will be performed if multiple binders are specified.
        optimize_kwargs : dict, optional
            The keyword arguments to pass to the optimizer, by default None
        copy : bool, optional
            If True, make a copy of the ligand, by default True
        """
        if self._core is None:
            self.make_core()

        if optimize_kwargs is None:
            optimize_kwargs = {
                "rotatron": optimizers.DistanceRotatron,
                "algorithm": "scipy",
            }
        else:
            optimize_kwargs.setdefault("rotatron", optimizers.DistanceRotatron)
            optimize_kwargs.setdefault("algorithm", "scipy")

        if optimize is None:
            optimize = len(binders) > 1

        ligand = ligand.copy() if copy else ligand

        if len(binders) == 1:

            b = ligand.get_atom(binders[0])
            core_h = self._core.get_atoms("H", by="element")[0]

            if delete is None:
                delete = ligand.get_neighbors(
                    b, filter=lambda x: x.element == "H"
                ).pop()
            else:
                delete = ligand.get_atom(delete)

            ligand.superimpose_to_bond(
                (b, delete),
                (core_h, self.metal),
            ).move(core_h.coord - b.coord)

            ligand._remove_atoms(delete)
            self._core._remove_atoms(core_h)
            self._core.merge(ligand)
            self._core._add_bond(b, self.metal)

        if acceptors is None:
            acceptors = self._core.get_atoms(*self._hydrogen_mapping[direction])
        else:
            acceptors = [self._core.get_atom(a) for a in acceptors]

        binders = [ligand.get_atom(b) for b in binders]

        if len(binders) > len(acceptors):
            raise ValueError("Not enough acceptor-hydrogen atoms to align the ligand")

        if delete is None:
            delete = [
                ligand.get_neighbors(b, filter=lambda x: x.element == "H").pop()
                for b in binders
            ]
        else:
            delete = [ligand.get_atom(d) for d in delete]

        b, binders = binders[0], binders[1:]
        core_h, core_Hs = acceptors[0], acceptors[1:]

        core_Hs = core_Hs[: len(binders)]
        reference_points = [H.coord for H in core_Hs]

        ligand.superimpose_to_bond(
            (b, delete[0]),
            (core_h, self.metal),
        ).move(core_h.coord - b.coord)

        if optimize:

            rotatron = self._make_constraint_rotatron_for_one(
                ligand,
                b,
                delete[0],
                binders,
                reference_points,
                optimize_kwargs.pop("rotatron"),
            )
            ligand = optimizers.optimize(ligand, rotatron, **optimize_kwargs)

        ligand._remove_atoms(*delete)
        self._core._remove_atoms(core_h, *core_Hs)
        self._core.merge(ligand)
        self._core._add_bonds(
            (b, self.metal), *((binder, self.metal) for binder in binders)
        )

        return self._core

    def add_ligands(
        self,
        ligands: list,
        binders: list,
        delete: list = None,
        acceptors: list = None,
        optimize: bool = None,
        optimize_kwargs: dict = None,
        copy: bool = True,
    ) -> "core.Molecule":
        """
        Add multiple ligands to the metal core. This method will only perform one optimization at the end of the alignment and should be faster than adding ligands one by one.

        Parameters
        ----------
        ligands: list
            A list of ligand molecules to add.
        binders : list
            Lists of atom objects or other atom identifiers to get atoms from the ligands,
            that will bind to the metal core. This is a list of lists where each sublist corresponds to a ligand's binder atoms.
        delete: list, optional
            Lists of atom objects or other atom identifiers to get atoms from the ligand,
            that will be deleted after alignment. Each of these must be a direct neighbor of a binder atom.
            They must be specified in the same order as the binders, by default None, in which case a random hydrogen neighbor from each binder is deleted.
            Similarly, this is a list of lists where each sublist corresponds to a ligand's delete atoms.
        acceptors : list, optional
            Lists of atom objects or other atom identifiers to get hydrogen atoms from the core. These will be used to align the ligands to the core.
            This is a list of lists where each sublist corresponds to a ligand's acceptor atoms.
        optimize : bool, optional
            If True, optimize the ligands after alignment, by default optimization will be performed if multiple binders or ligands are specified.
        optimize_kwargs : dict, optional
            The keyword arguments to pass to the optimizer, by default None
        copy : bool, optional
            If True, make a copy of the ligands before adding, by default True
        """
        if self._core is None:
            self.make_core()

        if optimize_kwargs is None:
            optimize_kwargs = {
                "rotatron": optimizers.DistanceRotatron,
                "algorithm": "scipy",
            }
        else:
            optimize_kwargs.setdefault("rotatron", optimizers.DistanceRotatron)
            optimize_kwargs.setdefault("algorithm", "scipy")

        if len(ligands) != len(binders):
            raise ValueError("The number of ligands and binders must match.")

        if acceptors is None:
            acceptors = self._core.get_atoms("H", by="element")
            if len(acceptors) < sum(len(b) for b in binders):
                raise ValueError(
                    "Not enough acceptor-hydrogen atoms to align the ligands"
                )

            _a = []
            for b in binders:
                _a.append([])
                for _ in b:
                    _a[-1].append(acceptors.pop(0))
            acceptors = _a
        else:
            acceptors = [[self._core.get_atom(_a) for _a in a] for a in acceptors]

        if optimize is None:
            optimize = len(binders[0]) > 1 or len(binders) > 1

        if copy:
            ligands = [ligand.copy() for ligand in ligands]

        _ligand_bonds = []
        _reference_points = {}
        _bonds_to_make = []
        _atoms_to_delete = []

        for idx, ligand in enumerate(ligands):

            _binders = binders[idx]
            if len(_binders) == 1:

                b = ligand.get_atom(_binders[0])
                core_h = acceptors[idx][0]

                if delete is None:
                    _delete = ligand.get_neighbors(
                        b, filter=lambda x: x.element == "H"
                    ).pop()
                else:
                    _delete = ligand.get_atom(delete[idx][0])

                ligand.superimpose_to_bond(
                    (b, _delete),
                    (core_h, self.metal),
                ).move(core_h.coord - b.coord)

                ligand._remove_atoms(_delete)
                self._core._remove_atoms(core_h)

                if optimize:
                    _ligand_bonds.extend(
                        ligand.get_atom_graph().find_edges(
                            root_node=b, exclude_cycles=True, exclude_locked=True
                        )
                    )
                    _ligand_bonds.append((self.metal, b))

                self._core.merge(ligand)
                self._core._add_bond(b, self.metal)

            else:

                _binders = [ligand.get_atom(b) for b in _binders]
                _acceptors = acceptors[idx]

                if len(_binders) > len(_acceptors):
                    raise ValueError(
                        "Not enough acceptor-hydrogen atoms to align the ligand"
                    )

                if delete is None:
                    _delete = [
                        ligand.get_neighbors(b, filter=lambda x: x.element == "H").pop()
                        for b in _binders
                    ]
                else:
                    _delete = [ligand.get_atom(d) for d in _delete]

                b, _binders = _binders[0], _binders[1:]
                core_h, core_Hs = _acceptors[0], _acceptors[1:]
                core_Hs = core_Hs[: len(_binders)]

                ligand.superimpose_to_bond(
                    (b, _delete[0]),
                    (core_h, self.metal),
                ).move(core_h.coord - b.coord)

                if optimize:
                    _ligand_bonds.extend(
                        ligand.get_atom_graph().find_edges(
                            root_node=b, exclude_cycles=True, exclude_locked=True
                        )
                    )
                    _ligand_bonds.append((self.metal, b))
                    _reference_points.update(
                        {b: h.coord for b, h in zip(_binders, core_Hs)}
                    )

                _bonds_to_make.extend(((self.metal, _b) for _b in _binders))
                _atoms_to_delete.extend(core_Hs + _delete[1:])

                ligand._remove_atoms(_delete[0])
                self._core._remove_atoms(core_h)
                self._core.merge(ligand)
                self._core._add_bond(b, self.metal)

        if optimize:
            graph = self._core.get_atom_graph()
            base_rotatron = optimize_kwargs.pop("rotatron")
            rotatron = base_rotatron(graph, _ligand_bonds)
            # if we don't have any reference points that we need to overlap
            # we can just use a normal rotatron and don't need a constraint
            if len(_reference_points) == 0:
                self._core = optimizers.optimize(
                    self._core, rotatron, **optimize_kwargs
                )
            else:
                nodes = list(graph.nodes)
                _reference_points = {
                    nodes.index(b): ref for b, ref in _reference_points.items()
                }

                def constraint(rotatron, state, **kwargs):
                    dist = 0
                    for idx, ref in _reference_points.items():
                        s = state[idx]
                        dist += (s - ref) ** 2

                    dist = dist.sum()
                    rotatron._target_dist = dist
                    return dist

                def finisher(rotatron, state, **kwargs):
                    return rotatron._target_dist < 0.1

                rotatron = optimizers.ConstraintRotatron(
                    rotatron,
                    constraint,
                    finisher,
                )

                self._core = optimizers.optimize(
                    self._core, rotatron, **optimize_kwargs
                )

        self._core._remove_atoms(*_atoms_to_delete)
        self._core._add_bonds(*_bonds_to_make)
        return self._core

    def _make_constraint_rotatron_for_one(
        self, ligand, root_binder, root_delete, binders, reference_points, base_rotatron
    ):

        def constraint(rotatron, state, binder_indices, reference_points, **kwargs):
            dist = 0
            for idx, ref in zip(binder_indices, reference_points):
                s = state[idx]
                dist += (s - ref) ** 2

            dist = dist.sum()
            rotatron._target_dist = dist
            return dist

        def finisher(rotatron, state, **kwargs):
            return rotatron._target_dist < 0.1

        graph = ligand.get_atom_graph()
        edges = graph.find_edges(
            root_node=root_binder, exclude_cycles=True, exclude_locked=True
        )
        edges.append((root_delete, root_binder))
        r = base_rotatron(graph, edges)

        binder_indices = [list(r.graph.nodes).index(b) for b in binders]
        out = optimizers.ConstraintRotatron(
            r,
            constraint,
            finisher,
            reference_points=reference_points,
            binder_indices=binder_indices,
        )
        return out


if __name__ == "__main__":

    iron = core.Atom.new("FE", pqr_charge=2)

    complexer = MetalComplexer(iron, structural.geometry.Octahedral())

    ligand = core.Molecule.from_smiles(
        "C=[N+](CC=[N+](Br)[H])[H]", id="LIG"
    ).autolabel()

    complexer.make_core(2)

    complex = complexer.add_ligands(
        ligands=[ligand] * 3,
        binders=[["N1", "N3"] for i in range(3)],
        acceptors=[("H1", "H5"), ("H6", "H4"), ("H2", "H3")],
        optimize_kwargs={
            "algorithm": "genetic",
            "rotatron": optimizers.DistanceRotatron,
        },
    )

    complex.show()
    pass
