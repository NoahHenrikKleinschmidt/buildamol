"""
The ConstraintRotatron allows for the optimization of a molecule's conformation while also accepting an additional constraint function that will also contribute to the evaluation. 
"""

import buildamol.optimizers.base_rotatron as Rotatron


__all__ = ["ConstraintRotatron"]


class ConstraintRotatron(Rotatron.Rotatron):
    """
    The ConstraintRotatron is a Meta-Rotatron environment that uses one of the other Rotatron environments to optimize the conformation of a molecule while also accepting an additional freely definable constraint function that will also contribute to the evaluation.

    Parameters
    ----------
    rotatron : Rotatron
        The rotatron object to use for basic conformer evaluation.
        This needs to be already set up and ready to use.
    constraint : callable
        A function that evaluates additional constraints and returnes a scalar value that will be added to the evaluation.
        This function will receive both the base-rotatron object, as well as the current state to evaluate as arguments,
        and can receive any additional arguments that are passes as keyword arguments during initialization of the ConstraintRotatron.
    finisher : callable
        The function that evaluates if the constraints are met and returns a boolean. This function will receive the base-rotatron object as first argument and the current state as second state.
        Also, all kwargs that are passed during initialization will be passed to this function (just like with the constraint function).
    **kwargs
        Additional keyword arguments to pass to the constraint function.
    """

    def __init__(
        self,
        rotatron: Rotatron,
        constraint: callable,
        finisher: callable = None,
        **kwargs
    ):
        self.rotatron = rotatron
        self.constraint = constraint
        self.finisher = finisher
        self.kwargs = kwargs

        self.action_space = self.rotatron.action_space
        self._bounds_tuple = self.rotatron._bounds_tuple
        self.rotatable_edges = self.rotatron.rotatable_edges

        if finisher is None:
            self.__step__ = self._step_without_finish
        else:
            self.__step__ = self._step_with_finish

    def eval(self, state):
        """
        Evaluate the state of the rotatron.

        Parameters
        ----------
        state : np.ndarray
            The state of the rotatron.

        Returns
        -------
        float
            The evaluation of the state.
        """
        return self.rotatron.eval(state) + self.constraint(
            self.rotatron, state, **self.kwargs
        )

    def reset(self):
        """
        Reset the rotatron to its initial state.
        """
        self.rotatron.reset()

    def step(self, action):
        """
        Perform a step in the rotatron.

        Parameters
        ----------
        action : np.ndarray
            The action to perform.

        Returns
        -------
        np.ndarray
            The new state of the rotatron.
        float
            The evaluation of the new state.
        bool
            Whether the rotatron is done.
        dict
            Additional information.
        """
        return self.__step__(action)

    def _step_with_finish(self, action):
        new_state, _eval, done, info = self.rotatron.step(action)
        _eval += self.constraint(self.rotatron, new_state, **self.kwargs)
        done = done and self.finisher(self.rotatron, new_state, **self.kwargs)
        return new_state, _eval, done, info

    def _step_without_finish(self, action):
        new_state, _eval, done, info = self.rotatron.step(action)
        _eval += self.constraint(self.rotatron, new_state, **self.kwargs)
        return new_state, _eval, done, info

    def done(self, state):
        """
        Check if the state is done.

        Parameters
        ----------
        state : np.ndarray
            The state to check.

        Returns
        -------
        bool
            Whether the state is done.
        """
        _done = self.rotatron.done(state)
        if self.finisher is not None:
            _done = _done and self.finisher(self.rotatron, state, **self.kwargs)
        return _done


if __name__ == "__main__":
    import buildamol as bam

    iron = bam.Atom("FE", coord=[0, 0, 0], pqr_charge=2)
    atoms, bonds = bam.structural.geometry.octahedral.fill_hydrogens(iron)
    mol = bam.Molecule.new(atoms=atoms, bonds=bonds)

    # the 1.7 bond length is a little cheating, but it's the only way of making the
    # octahedral geometry thing work with the current implementation
    for bond in mol.get_bonds():
        mol.adjust_bond_length(*bond, 1.7)

    # now add the ligands
    ligand = bam.Molecule.from_smiles("C=[N+](CC=[N+](Br)[H])[H]", id="LIG").autolabel()

    ligand.superimpose_to_bond(
        ligand.get_atoms("N1", "HN1"), mol.get_atoms("H1", "FE", keep_order=True)
    )

    mol.remove_atoms("H1")
    ligand.remove_atoms("HN1")
    mol.merge(ligand).add_bond("N1", "FE")

    H5_ref_coord = mol.get_atom("H5").coord

    def constraint(rotatron, state, **kwargs):
        dist = state[rotatron.HN3_idx] - H5_ref_coord
        dist = (dist**2).sum()
        rotatron._target_dist = dist
        return dist

    def finisher(rotatron, state, **kwargs):
        return rotatron._target_dist < 0.1

    graph = mol.get_atom_graph()
    edges = graph.find_edges(root_node=mol.get_atom("FE"), min_descendants=2)
    rotatron = bam.optimizers.DistanceRotatron(graph, edges)
    rotatron.HN3_idx = list(rotatron.graph.nodes).index(mol.get_atom("HN3"))

    constraint_rotatron = ConstraintRotatron(rotatron, constraint, finisher)

    v = mol.draw(atoms=False)
    out = bam.optimize(mol.copy(), constraint_rotatron, algorithm="scipy")
    v += out.draw(atoms=False, line_color="red")
    v.show()
    pass
