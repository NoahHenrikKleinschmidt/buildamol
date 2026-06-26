"""
OpenFF/OpenMM parameterization and simulation pipeline for BuildAMol molecules.

This module provides tools to convert BuildAMol molecules into fully parameterized
OpenMM simulation systems using the OpenFF SMIRNOFF force field family.

Requires
--------
- openff-toolkit  (``conda install -c conda-forge openff-toolkit``)
- openff-units    (installed with openff-toolkit)
- openmm          (``conda install -c conda-forge openmm``)

Usage
-----

First parametrize a molecule with the OpenFF force field

.. code-block:: python

    import buildamol as bam
    from buildamol.extensions.molecular_dynamics import openff_md as offmd

    aspirin = bam.molecule("aspirin")

    # Parametrize and run MD in one go
    system = offmd.parametrize(aspirin)

    # Or use the class for more control
    p = offmd.OpenFFParameterizer(force_field="openff-2.1.0.offxml")
    system = p.parametrize(aspirin)

Then run energy minimization and/or molecular dynamics:

.. code-block:: python

    system.minimize()
    system.run(10_000, reporter="md_log.csv")
    aspirin_relaxed = system.to_molecule()


Notes
-----
Partial charge assignment uses ``mmff94`` by default, which requires only RDKit
(always available alongside BuildAMol). For production MD simulations the
OpenFF-recommended method is ``am1bccelf10``, which requires AmberTools or the
OpenEye toolkit to be installed. Use ``charge_method='gasteiger'`` for the
fastest possible option with no additional dependencies.
"""

import sys
from io import StringIO

# openff.interchange checks for foyer by actually importing it. In some environments
# foyer is installed but broken (e.g. gmso/pydantic v1 conflict), which prevents
# interchange from loading. Marking it as unavailable here is safe because foyer
# is only used for non-SMIRNOFF force fields; we never need it.
if "foyer" not in sys.modules:
    sys.modules["foyer"] = None  # type: ignore[assignment]

try:
    from openff.toolkit import Molecule as OFFMolecule
    from openff.toolkit import Topology as OFFTopology
    from openff.toolkit.typing.engines.smirnoff import ForceField
    from openff.units import unit as off_unit

    HAS_OPENFF = True
except ImportError:
    OFFMolecule = None
    OFFTopology = None
    ForceField = None
    off_unit = None
    HAS_OPENFF = False

try:
    import openmm
    from openmm import unit as mm_unit
    from openmm.app import DCDReporter, PDBFile, Simulation, StateDataReporter

    HAS_OPENMM = True
except ImportError:
    openmm = None
    mm_unit = None
    DCDReporter = None
    PDBFile = None
    Simulation = None
    StateDataReporter = None
    HAS_OPENMM = False


DEFAULT_FORCEFIELD = "openff-2.1.0.offxml"
DEFAULT_CHARGE_METHOD = "mmff94"
DEFAULT_TEMPERATURE = 300  # K
DEFAULT_FRICTION = 1.0  # 1/ps
DEFAULT_TIMESTEP = 0.002  # ps


def _require_openff():
    if not HAS_OPENFF:
        raise ImportError(
            "openff-toolkit is required for this feature. "
            "Install with: conda install -c conda-forge openff-toolkit"
        )


def _require_openmm():
    if not HAS_OPENMM:
        raise ImportError(
            "OpenMM is required for this feature. "
            "Install with: conda install -c conda-forge openmm"
        )


def to_openff(mol, charge_method: str = DEFAULT_CHARGE_METHOD) -> "OFFMolecule":
    """
    Convert a BuildAMol molecule to an OpenFF Molecule.

    The conversion preserves 3D coordinates. Connectivity and bond orders are
    read from the molecule's RDKit representation; if that fails (e.g. for
    molecules with unusual aromaticity), the function falls back to deriving
    topology from the SMILES string and generating a fresh ETKDG conformer.

    Parameters
    ----------
    mol : Molecule
        The BuildAMol molecule. Must have 3D coordinates.
    charge_method : str
        Partial charge method passed to
        ``openff.toolkit.Molecule.assign_partial_charges``. Common options:

        - ``'mmff94'`` (default) — RDKit-based, no extra dependencies needed
        - ``'am1bccelf10'`` — OpenFF-recommended for production, requires AmberTools or OpenEye
        - ``'am1bcc'`` — classic AM1-BCC, requires AmberTools
        - ``'gasteiger'`` — fastest option, dependency-free, lowest accuracy

    Returns
    -------
    openff.toolkit.Molecule
        An OpenFF Molecule with a 3D conformer and assigned partial charges.
    """
    _require_openff()

    off_mol = _build_openff_molecule(mol)
    off_mol.assign_partial_charges(charge_method)
    return off_mol


def _build_openff_molecule(mol) -> "OFFMolecule":
    """
    Build an OpenFF Molecule from a BuildAMol molecule, preserving 3D coordinates.

    Primary path: RDKit mol from to_rdkit() → OFFMolecule.from_rdkit().
    This keeps the atom ordering and 3D geometry intact.

    Fallback: SMILES-based topology + ETKDG conformer (loses BuildAMol geometry
    but is more tolerant of unusual bond perception from the PDB intermediate).
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem

    # --- Primary path: via RDKit (preserves 3D coords and atom ordering) ---
    try:
        rdmol = mol.to_rdkit()
        Chem.SanitizeMol(rdmol)
        off_mol = OFFMolecule.from_rdkit(rdmol, allow_undefined_stereo=True)
        # Verify a conformer came through
        if off_mol.n_conformers == 0:
            _embed_conformer_from_rdkit(off_mol, rdmol)
        return off_mol
    except Exception:
        pass

    # --- Fallback: SMILES topology + ETKDG geometry ---
    smiles = mol.to_smiles(isomeric=True, write_hydrogens=False)
    off_mol = OFFMolecule.from_smiles(smiles, allow_undefined_stereo=True)

    # Try to embed BuildAMol coordinates onto the SMILES-ordered mol
    try:
        _embed_conformer_from_bam(off_mol, mol)
    except Exception:
        # Last resort: generate fresh ETKDG conformer
        rdmol_ff = off_mol.to_rdkit()
        AllChem.EmbedMolecule(rdmol_ff, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(rdmol_ff)
        positions = rdmol_ff.GetConformer().GetPositions()
        off_mol.add_conformer(off_unit.Quantity(positions, off_unit.angstrom))

    return off_mol


def _embed_conformer_from_rdkit(off_mol, rdmol):
    """Add the conformer from an RDKit mol to an OpenFF molecule."""
    conf = rdmol.GetConformer()
    positions = conf.GetPositions()
    off_mol.add_conformer(off_unit.Quantity(positions, off_unit.angstrom))


def _embed_conformer_from_bam(off_mol, bam_mol):
    """
    Map BuildAMol heavy-atom coordinates onto the SMILES-ordered OpenFF molecule
    using MCS atom matching, then refine H positions with MMFF.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS

    rdmol_off = off_mol.to_rdkit()  # OpenFF atom ordering, no conformer
    rdmol_bam = bam_mol.to_rdkit()  # BuildAMol atom ordering, has conformer

    # Find maximum common substructure to map atom indices
    mcs = rdFMCS.FindMCS(
        [rdmol_off, rdmol_bam],
        completeRingsOnly=True,
        bondCompare=rdFMCS.BondCompare.CompareAny,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        timeout=5,
    )
    query = Chem.MolFromSmarts(mcs.smartsString)
    off_match = rdmol_off.GetSubstructMatch(query)
    bam_match = rdmol_bam.GetSubstructMatch(query)

    off_to_bam = dict(zip(off_match, bam_match))
    bam_conf = rdmol_bam.GetConformer()

    rwmol = Chem.RWMol(rdmol_off)
    conf = Chem.Conformer(rwmol.GetNumAtoms())

    # Place MCS-matched atoms at BuildAMol positions
    for off_idx in range(rwmol.GetNumAtoms()):
        if off_idx in off_to_bam:
            pos = bam_conf.GetAtomPosition(off_to_bam[off_idx])
        else:
            # Unmapped atom (usually H added by OpenFF) — use centroid as placeholder
            pos = bam_conf.GetAtomPosition(next(iter(off_to_bam.values())))
        conf.SetAtomPosition(off_idx, pos)

    rwmol.AddConformer(conf, assignId=True)

    # MMFF refinement to place unmatched atoms (H) correctly
    AllChem.MMFFOptimizeMolecule(rwmol, maxIters=500, mmffVariant="MMFF94s")

    positions = rwmol.GetConformer().GetPositions()
    off_mol.add_conformer(off_unit.Quantity(positions, off_unit.angstrom))


class OpenFFParameterizer:
    """
    Parametrize BuildAMol molecules using the OpenFF SMIRNOFF force field.

    Parameters
    ----------
    force_field : str
        The OpenFF force field to use. Defaults to ``'openff-2.1.0.offxml'``.
        Any force field file resolvable by the openff-toolkit can be used
        (e.g. ``'openff_unconstrained-2.1.0.offxml'``).
    charge_method : str
        Partial charge assignment method. Defaults to ``'am1bccelf10'``.

    Examples
    --------
    .. code-block:: python

        import buildamol as bam
        from buildamol.extensions.molecular_dynamics import openff as offmd

        p = offmd.OpenFFParameterizer()

        aspirin = bam.molecule("aspirin")
        system = p.parametrize(aspirin)
        system.minimize().save_pdb("aspirin_min.pdb")

        # Reuse the same parametrizer for multiple molecules
        ibuprofen = bam.molecule("ibuprofen")
        system2 = p.parametrize(ibuprofen)
    """

    def __init__(
        self,
        force_field: str = DEFAULT_FORCEFIELD,
        charge_method: str = DEFAULT_CHARGE_METHOD,
    ):
        _require_openff()
        _require_openmm()
        self.force_field = force_field
        self.charge_method = charge_method
        self._ff = None

    @property
    def ff(self) -> "ForceField":
        """The underlying OpenFF ForceField object (lazy-loaded)."""
        if self._ff is None:
            self._ff = ForceField(self.force_field)
        return self._ff

    def parametrize(self, mol) -> "OpenMMSystem":
        """
        Parametrize a BuildAMol molecule and return an :class:`OpenMMSystem`.

        Parameters
        ----------
        mol : Molecule
            The BuildAMol molecule to parametrize. Must have 3D coordinates.

        Returns
        -------
        OpenMMSystem
            A ready-to-simulate system with the OpenMM System, Topology,
            and initial coordinates.
        """
        off_mol = to_openff(mol, charge_method=self.charge_method)
        topology = off_mol.to_topology()
        # Pass charge_from_molecules so the pre-computed partial charges are used
        # rather than being recalculated inside create_openmm_system.
        system = self.ff.create_openmm_system(topology, charge_from_molecules=[off_mol])
        return OpenMMSystem(system, topology, off_mol, source_mol=mol)

    def __call__(self, mol) -> "OpenMMSystem":
        return self.parametrize(mol)


class OpenMMSystem:
    """
    A parameterized molecular system ready for OpenMM simulation.

    Instances are returned by :func:`parametrize` or
    :meth:`OpenFFParameterizer.parametrize`. The simulation is set up lazily
    on the first call to :meth:`minimize` or :meth:`run`.

    Parameters
    ----------
    system : openmm.System
        The OpenMM system with force field parameters.
    topology : openff.toolkit.Topology
        The OpenFF topology (used to build the OpenMM topology).
    off_mol : openff.toolkit.Molecule
        The OpenFF molecule carrying the 3D conformer used as initial positions.
    source_mol : Molecule, optional
        The original BuildAMol molecule. When provided, :meth:`to_molecule`
        returns a copy of this molecule with MD coordinates applied, preserving
        bond orders and residue structure. If ``None``, :meth:`to_molecule`
        falls back to a PDB round-trip (bond orders are lost).
    """

    def __init__(self, system, topology, off_mol, source_mol=None):
        _require_openmm()
        self._system = system
        self._topology = topology
        self._off_mol = off_mol
        self._source_mol = source_mol
        self._simulation = None

    # ------------------------------------------------------------------
    # Simulation setup
    # ------------------------------------------------------------------

    def _build_simulation(
        self,
        temperature: float = DEFAULT_TEMPERATURE,
        friction: float = DEFAULT_FRICTION,
        timestep: float = DEFAULT_TIMESTEP,
    ) -> "Simulation":
        """
        Construct the OpenMM Simulation with a Langevin integrator.

        Parameters
        ----------
        temperature : float
            Temperature in Kelvin (default 300).
        friction : float
            Friction coefficient in 1/ps (default 1.0).
        timestep : float
            Integration timestep in ps (default 0.002).
        """
        integrator = openmm.LangevinMiddleIntegrator(
            temperature * mm_unit.kelvin,
            friction / mm_unit.picosecond,
            timestep * mm_unit.picoseconds,
        )
        simulation = Simulation(
            self._topology.to_openmm(),
            self._system,
            integrator,
        )
        simulation.context.setPositions(self._off_mol.conformers[0].to_openmm())
        self._simulation = simulation
        return simulation

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def minimize(
        self,
        tolerance: float = 10.0,
        max_iterations: int = 0,
        **simulation_kwargs,
    ) -> "OpenMMSystem":
        """
        Run energy minimization.

        Parameters
        ----------
        tolerance : float
            Energy convergence tolerance in kJ/mol (default 10.0).
        max_iterations : int
            Maximum minimization iterations. 0 means run until convergence.
        **simulation_kwargs
            Forwarded to :meth:`_build_simulation` if the simulation has not
            been set up yet (``temperature``, ``friction``, ``timestep``).

        Returns
        -------
        OpenMMSystem
            Returns ``self`` to allow method chaining.
        """
        if self._simulation is None:
            self._build_simulation(**simulation_kwargs)
        self._simulation.minimizeEnergy(
            tolerance=tolerance * mm_unit.kilojoules_per_mole,
            maxIterations=max_iterations,
        )
        return self

    def run(
        self,
        n_steps: int,
        reporter: str = None,
        report_interval: int = 1000,
        dcd_reporter: str = None,
        dcd_interval: int = None,
        **simulation_kwargs,
    ) -> "OpenMMSystem":
        """
        Run MD simulation.

        Parameters
        ----------
        n_steps : int
            Number of integration steps to perform.
        reporter : str, optional
            Path to a CSV file for energy/temperature logging. If ``None``,
            no CSV reporter is attached.
        report_interval : int
            Steps between CSV reporter outputs (default 1000).
        dcd_reporter : str, optional
            Path to a DCD trajectory file. DCD is a compact binary format
            widely supported by trajectory analysis tools (MDAnalysis, VMD,
            GROMACS, etc.). If ``None``, no DCD reporter is attached.
        dcd_interval : int, optional
            Steps between DCD frames. Defaults to ``report_interval`` if not
            specified.
        **simulation_kwargs
            Forwarded to :meth:`_build_simulation` if the simulation has not
            been set up yet (``temperature``, ``friction``, ``timestep``).

        Returns
        -------
        OpenMMSystem
            Returns ``self`` to allow method chaining.

        Examples
        --------
        .. code-block:: python

            # CSV log + DCD trajectory, both every 500 steps
            system.run(50_000, reporter="md.csv", dcd_reporter="traj.dcd", report_interval=500)

            # DCD only, frame every 100 steps
            system.run(10_000, dcd_reporter="traj.dcd", dcd_interval=100)
        """
        if self._simulation is None:
            self._build_simulation(**simulation_kwargs)
        if reporter is not None:
            self._simulation.reporters.append(
                StateDataReporter(
                    reporter,
                    report_interval,
                    step=True,
                    potentialEnergy=True,
                    temperature=True,
                    progress=True,
                    totalSteps=n_steps,
                )
            )
        if dcd_reporter is not None:
            self._simulation.reporters.append(
                DCDReporter(
                    dcd_reporter,
                    dcd_interval if dcd_interval is not None else report_interval,
                )
            )
        self._simulation.step(n_steps)
        return self

    def get_energy(self) -> float:
        """
        Return the current potential energy in kJ/mol.

        Returns
        -------
        float
            Potential energy in kJ/mol.

        Raises
        ------
        RuntimeError
            If called before :meth:`minimize` or :meth:`run`.
        """
        if self._simulation is None:
            raise RuntimeError(
                "No simulation has been set up yet. Call minimize() or run() first."
            )
        state = self._simulation.context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(mm_unit.kilojoules_per_mole)

    def save_pdb(self, path: str) -> None:
        """
        Write the current atomic positions to a PDB file.

        Parameters
        ----------
        path : str
            Output file path.

        Raises
        ------
        RuntimeError
            If called before :meth:`minimize` or :meth:`run`.
        """
        if self._simulation is None:
            raise RuntimeError(
                "No simulation has been set up yet. Call minimize() or run() first."
            )
        state = self._simulation.context.getState(getPositions=True)
        with open(path, "w") as fh:
            PDBFile.writeFile(
                self._simulation.topology,
                state.getPositions(),
                fh,
            )

    def to_molecule(self):
        """
        Return the current atomic positions as a BuildAMol Molecule.

        When the system was created via :func:`parametrize` or
        :meth:`OpenFFParameterizer.parametrize`, bond orders, residue structure,
        and atom names are fully preserved — only coordinates are updated from
        the MD state. This is the recommended path.

        If no source molecule is available (system was constructed manually),
        the method falls back to a PDB round-trip, which loses bond orders.

        Returns
        -------
        Molecule
            A BuildAMol molecule with coordinates from the current simulation state.

        Raises
        ------
        RuntimeError
            If called before :meth:`minimize` or :meth:`run`.
        """
        if self._simulation is None:
            raise RuntimeError(
                "No simulation has been set up yet. Call minimize() or run() first."
            )

        state = self._simulation.context.getState(getPositions=True)
        # Positions come back in nanometres; BuildAMol uses Angstroms.
        positions_ang = state.getPositions(asNumpy=True).value_in_unit(mm_unit.angstrom)

        if self._source_mol is not None:
            # Fast path: copy the original molecule and update only coordinates.
            # Atom ordering is preserved through the BuildAMol→RDKit→OpenFF→OpenMM
            # chain, so positional index i in OpenMM matches atom i in get_atoms().
            mol = self._source_mol.copy()
            mol.set_coords(positions_ang)
            return mol

        # Fallback when no source molecule is stored: PDB round-trip.
        # Bond orders will be lost (PDB format has no bond-order field).
        import buildamol as bam

        buf = StringIO()
        PDBFile.writeFile(
            self._simulation.topology,
            state.getPositions(),
            buf,
        )
        buf.seek(0)
        return bam.Molecule._from_pdb_string(buf.read())

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def system(self):
        """The underlying ``openmm.System`` object."""
        return self._system

    @property
    def topology(self):
        """The underlying ``openff.toolkit.Topology`` object."""
        return self._topology

    @property
    def simulation(self):
        """
        The ``openmm.app.Simulation`` object, or ``None`` if not yet built.
        Set up by the first call to :meth:`minimize` or :meth:`run`.
        """
        return self._simulation


# ------------------------------------------------------------------
# Convenience function
# ------------------------------------------------------------------


def parametrize(
    mol,
    force_field: str = DEFAULT_FORCEFIELD,
    charge_method: str = DEFAULT_CHARGE_METHOD,
) -> "OpenMMSystem":
    """
    Parametrize a BuildAMol molecule with the OpenFF SMIRNOFF force field.

    This is a convenience wrapper around :class:`OpenFFParameterizer`. For
    repeated use on multiple molecules, instantiate the parametrizer directly
    to avoid reloading the force field each time.

    Parameters
    ----------
    mol : Molecule
        The BuildAMol molecule. Must have 3D coordinates.
    force_field : str
        OpenFF force field file name (default ``'openff-2.1.0.offxml'``).
    charge_method : str
        Partial charge method (default ``'am1bccelf10'``).

    Returns
    -------
    OpenMMSystem
        A ready-to-simulate system.

    Examples
    --------
    .. code-block:: python

        import buildamol as bam
        from buildamol.extensions.molecular_dynamics import openff as offmd

        aspirin = bam.molecule("aspirin")
        system = offmd.parametrize(aspirin)

        # Energy minimization
        system.minimize()
        print(f"Energy after minimization: {system.get_energy():.2f} kJ/mol")

        # Short MD run with logging
        system.run(10_000, reporter="md_log.csv", report_interval=500)

        # Save final structure and bring it back into BuildAMol
        system.save_pdb("aspirin_md.pdb")
        mol_relaxed = system.to_molecule()
    """
    return OpenFFParameterizer(force_field, charge_method).parametrize(mol)


# ------------------------------------------------------------------
# Module-level self-test
# ------------------------------------------------------------------

if __name__ == "__main__":
    import buildamol as bam

    bam.load_small_molecules()
    aspirin = bam.molecule("aspirin")

    print("Parametrizing aspirin with OpenFF...")
    system = parametrize(aspirin)  # uses mmff94 charges by default

    print("Running energy minimization...")
    system.minimize()
    print(f"Potential energy: {system.get_energy():.2f} kJ/mol")

    system.save_pdb("/tmp/aspirin_minimized.pdb")
    print("Saved to /tmp/aspirin_minimized.pdb")

    mol_out = system.to_molecule()
    print(f"Got back molecule with {mol_out.count_atoms()} atoms")
