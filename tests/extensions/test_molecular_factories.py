"""
Tests for buildamol.extensions.molecular_factories pipeline building blocks.

Coverage:
  - Context dataclass
  - ChainableBlock / ChainedBlock (| operator)
  - Source blocks: Choice, Compound
  - Transform blocks: Modify, SetAttachResidue, SetLinkerAtoms
  - Annotation blocks: FindLinkerAtoms
  - Combinator blocks: Connect
  - Integration: multi-step pipelines including the oxazolidinone tutorial pattern
"""

import pytest
import buildamol as bam

from buildamol.extensions.molecular_factories.base import (
    ChainableBlock,
    ChainedBlock,
    Context,
)
from buildamol.extensions.molecular_factories.sources import Choice, Compound
from buildamol.extensions.molecular_factories.modifiers import (
    Connect,
    FindLinkerAtoms,
    Modify,
    SetAttachResidue,
    SetLinkerAtoms,
    Random,
    Optimize,
)
from buildamol.extensions.molecular_factories.optimizer import Optimizer

# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def ethane():
    return bam.read_smiles("CC", id="eth").autolabel()


@pytest.fixture(scope="module")
def propane():
    return bam.read_smiles("CCC", id="pro").autolabel()


@pytest.fixture(scope="module")
def benzene():
    return bam.read_smiles("c1ccccc1", id="ben").autolabel()


@pytest.fixture(scope="module")
def oxazolidinone_core():
    # O1C(=O)NCC1  — the 5-membered oxazolidinone ring used in the tutorial
    return bam.read_smiles("O1C(=O)NCC1", id="ox").autolabel()


# ── Context ───────────────────────────────────────────────────────────────────


def test_context_defaults():
    ctx = Context()
    assert ctx.molecule is None
    assert ctx.attach_residue is None
    assert ctx.linker_atom is None
    assert ctx.deleter_atom is None


# ── ChainableBlock / ChainedBlock ─────────────────────────────────────────────


def test_pipe_operator_creates_chained_block(ethane):
    a = Choice([ethane])
    b = Modify(lambda mol: mol)
    chain = a | b
    assert isinstance(chain, ChainedBlock)


def test_chained_block_calls_blocks_in_order(ethane):
    log = []

    class First(ChainableBlock):
        def __call__(self, *args, **kwargs):
            log.append("first")
            ctx = Context()
            ctx.molecule = ethane
            return ctx

    class Second(ChainableBlock):
        def __call__(self, ctx, *args, **kwargs):
            log.append("second")
            return ctx

    (First() | Second())()
    assert log == ["first", "second"]


def test_triple_chain_order(ethane):
    log = []

    class Tag(ChainableBlock):
        def __init__(self, tag):
            self.tag = tag

        def __call__(self, ctx=None, *args, **kwargs):
            log.append(self.tag)
            if ctx is None:
                ctx = Context()
                ctx.molecule = ethane
            return ctx

    (Tag("a") | Tag("b") | Tag("c"))()
    assert log == ["a", "b", "c"]


def test_pipe_type_error(ethane):
    with pytest.raises(TypeError):
        Choice([ethane]) | "not_a_block"


# ── Choice ────────────────────────────────────────────────────────────────────


def test_choice_returns_context(ethane, propane):
    ctx = Choice([ethane, propane])()
    assert isinstance(ctx, Context)
    assert ctx.molecule is not None


def test_choice_picks_from_provided_library(ethane, propane):
    block = Choice([ethane, propane])
    for _ in range(30):
        ctx = block()
        assert ctx.molecule is ethane or ctx.molecule is propane


def test_choice_single_molecule_always_returns_it(ethane):
    for _ in range(5):
        ctx = Choice([ethane])()
        assert ctx.molecule is ethane


def test_choice_with_probabilities(ethane, propane):
    # p=[1, 0] must always return ethane
    block = Choice([ethane, propane], p=[1.0, 0.0])
    for _ in range(10):
        assert block().molecule is ethane


# ── Compound ──────────────────────────────────────────────────────────────────


def test_compound_returns_context(ethane):
    ctx = Compound(ethane)()
    assert isinstance(ctx, Context)
    assert ctx.molecule is not None


def test_compound_copies_by_default(ethane):
    ctx = Compound(ethane)()
    assert ctx.molecule is not ethane


def test_compound_no_copy(ethane):
    ctx = Compound(ethane, copy=False)()
    assert ctx.molecule is ethane


def test_compound_repeated_calls_return_independent_copies(ethane):
    block = Compound(ethane)
    ctx1 = block()
    ctx2 = block()
    assert ctx1.molecule is not ctx2.molecule


# ── Modify ────────────────────────────────────────────────────────────────────


def test_modify_calls_function_with_molecule(ethane):
    seen = []

    def fn(mol):
        seen.append(mol)
        return mol

    ctx = Context()
    ctx.molecule = ethane
    Modify(fn)(ctx)
    assert len(seen) == 1 and seen[0] is ethane


def test_modify_replaces_molecule_with_return_value(ethane, propane):
    ctx = Context()
    ctx.molecule = ethane
    out = Modify(lambda _: propane)(ctx)
    assert out.molecule is propane


def test_modify_passes_context_through(ethane):
    ctx = Context()
    ctx.molecule = ethane
    ctx.linker_atom = "sentinel"
    out = Modify(lambda mol: mol)(ctx)
    assert out.linker_atom == "sentinel"


# ── SetAttachResidue ──────────────────────────────────────────────────────────


def test_set_attach_residue_by_seqid(ethane):
    ctx = Context()
    ctx.molecule = ethane.copy()
    out = SetAttachResidue(1)(ctx)
    assert out.attach_residue is not None
    assert out.attach_residue is out.molecule.get_attach_residue()


def test_set_attach_residue_callable(ethane):
    mol = ethane.copy()
    res = mol.residues[0]

    ctx = Context()
    ctx.molecule = mol
    out = SetAttachResidue(lambda m: m.residues[0])(ctx)
    assert out.attach_residue is res


# ── SetLinkerAtoms ────────────────────────────────────────────────────────────


def test_set_linker_atoms_by_string_id(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = SetLinkerAtoms(link="C1", delete="H11")(ctx)
    assert out.linker_atom is not None
    assert out.deleter_atom is not None
    assert out.linker_atom.id == "C1"
    assert out.deleter_atom.id == "H11"


def test_set_linker_atoms_callable_link(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = SetLinkerAtoms(link=lambda m: m.get_atom("C1"))(ctx)
    assert out.linker_atom is not None
    assert out.linker_atom.id == "C1"


def test_set_linker_atoms_only_link_set(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.deleter_atom = "previous"

    out = SetLinkerAtoms(link="C1")(ctx)
    assert out.linker_atom.id == "C1"
    assert out.deleter_atom == "previous"  # unchanged


def test_set_linker_atoms_only_delete_set(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = "previous"

    out = SetLinkerAtoms(delete="H11")(ctx)
    assert out.linker_atom == "previous"  # unchanged
    assert out.deleter_atom.id == "H11"


# ── FindLinkerAtoms ───────────────────────────────────────────────────────────


def test_find_linker_atoms_invalid_how():
    with pytest.raises(ValueError):
        FindLinkerAtoms(how="nonsense")


def test_find_linker_atoms_random_sets_linker_and_deleter(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms(how="random")(ctx)
    assert out.linker_atom is not None
    assert out.deleter_atom is not None


def test_find_linker_atoms_linker_is_non_hydrogen(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms()(ctx)
    assert out.linker_atom.element != "H"


def test_find_linker_atoms_deleter_is_hydrogen(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms()(ctx)
    assert out.deleter_atom.element == "H"


def test_find_linker_atoms_deleter_is_neighbor_of_linker(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms()(ctx)
    nbrs = out.linker_atom.get_neighbors()
    assert out.deleter_atom in nbrs


def test_find_linker_atoms_callable_how(ethane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])
    c1 = mol.get_atom("C1")
    h11 = mol.get_atom("H11")

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms(how=lambda m: (c1, h11))(ctx)
    assert out.linker_atom is c1
    assert out.deleter_atom is h11


def test_find_linker_atoms_falls_back_to_first_residue_when_no_attach_residue(ethane):
    """FindLinkerAtoms must not crash when get_attach_residue() returns None."""
    mol = ethane.copy()
    # deliberately do NOT call set_attach_residue

    ctx = Context()
    ctx.molecule = mol

    out = FindLinkerAtoms()(ctx)
    assert out.linker_atom is not None
    assert out.deleter_atom is not None


# ── Connect ───────────────────────────────────────────────────────────────────


def test_connect_raises_when_context_has_no_linker(ethane, propane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = None
    ctx.deleter_atom = mol.get_atom("H11")

    R = Compound(propane) | FindLinkerAtoms()
    with pytest.raises(ValueError, match="[Ll]inker atom"):
        Connect(R)(ctx)


def test_connect_without_explicit_deleter_auto_detects_hydrogen(ethane, propane):
    """
    When deleter_atom is None the linkage auto-detects a H on the linker atom,
    so Connect should succeed and produce a valid molecule.
    This is the behaviour exploited by SetLinkerAtoms(link=_N) without delete.
    """
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = mol.get_atom("C1")
    ctx.deleter_atom = None  # auto-detect

    R = Compound(propane) | FindLinkerAtoms()
    out = Connect(R)(ctx)

    assert out.molecule is not None
    assert out.molecule.count_residues() == 2


def test_connect_raises_when_source_has_no_linker(ethane):
    """If the source pipeline doesn't set a linker, Connect must raise."""

    class EmptySource(ChainableBlock):
        def __call__(self, *args, **kwargs):
            ctx = Context()
            ctx.molecule = ethane.copy()
            # linker_atom stays None
            return ctx

    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = mol.get_atom("C1")
    ctx.deleter_atom = mol.get_atom("H11")

    with pytest.raises(ValueError, match="linker atom"):
        Connect(EmptySource())(ctx)


def test_connect_produces_context_with_molecule(ethane, propane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = mol.get_atom("C1")
    ctx.deleter_atom = mol.get_atom("H11")

    source = Compound(propane) | FindLinkerAtoms()
    out = Connect(source)(ctx)

    assert isinstance(out, Context)
    assert out.molecule is not None


def test_connect_result_has_more_atoms_than_either_fragment(ethane, propane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])
    n_eth = len(mol.atoms)
    n_pro = len(propane.atoms)

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = mol.get_atom("C1")
    ctx.deleter_atom = mol.get_atom("H11")

    source = Compound(propane) | FindLinkerAtoms()
    out = Connect(source)(ctx)

    # each connection removes one H from each side → n_eth + n_pro - 2
    assert len(out.molecule.atoms) == n_eth + n_pro - 2


def test_connect_result_has_two_residues(ethane, propane):
    mol = ethane.copy()
    mol.set_attach_residue(mol.residues[0])

    ctx = Context()
    ctx.molecule = mol
    ctx.linker_atom = mol.get_atom("C1")
    ctx.deleter_atom = mol.get_atom("H11")

    source = Compound(propane) | FindLinkerAtoms()
    out = Connect(source)(ctx)

    assert out.molecule.count_residues() == 2


# ── Pipeline integration ──────────────────────────────────────────────────────


def test_source_pipeline_sets_linker_atoms(ethane, propane):
    """Choice | FindLinkerAtoms gives a context with both linker and deleter set."""
    pipeline = Choice([ethane, propane]) | FindLinkerAtoms()
    ctx = pipeline()
    assert ctx.linker_atom is not None
    assert ctx.deleter_atom is not None


def test_two_fragment_pipeline(ethane, propane):
    """Compound | SetAttachResidue | SetLinkerAtoms | Connect(source) end-to-end."""
    source = Compound(propane) | FindLinkerAtoms()
    pipeline = (
        Compound(ethane)
        | SetAttachResidue(1)
        | SetLinkerAtoms(link="C1", delete="H11")
        | Connect(source)
    )
    ctx = pipeline()
    assert ctx.molecule is not None
    assert ctx.molecule.count_residues() == 2


def test_pipeline_can_run_multiple_times_independently(ethane, propane):
    """Each pipeline() call produces a fresh independent molecule."""
    source = Compound(propane) | FindLinkerAtoms()
    pipeline = (
        Compound(ethane)
        | SetAttachResidue(1)
        | SetLinkerAtoms(link="C1", delete="H11")
        | Connect(source)
    )
    results = [pipeline() for _ in range(3)]
    molecules = [r.molecule for r in results]
    # all molecules should be distinct objects
    assert molecules[0] is not molecules[1]
    assert molecules[1] is not molecules[2]


def test_oxazolidinone_core_connects_aryl_to_nitrogen(oxazolidinone_core, benzene):
    """
    Reproduce the key structural requirement from the tutorial:
    the oxazolidinone N should be connectable to an aryl fragment via a pipeline.
    """
    core = oxazolidinone_core.copy()

    # Locate N by element and pick one of its H neighbours as deleter
    N = core.get_atom("N", by="element")
    N_h = next(iter(N.get_hydrogens()))

    aryl_source = Compound(benzene) | FindLinkerAtoms()

    # Set context manually to avoid relying on SetLinkerAtoms Atom-object handling
    ctx = Context()
    ctx.molecule = core
    ctx.linker_atom = N
    ctx.deleter_atom = N_h

    out = Connect(aryl_source)(ctx)

    assert out.molecule is not None
    assert out.molecule.count_residues() == 2


def test_two_connect_steps_produces_three_residue_molecule(ethane, propane, benzene):
    """
    A pipeline with two successive Connect steps must produce a 3-residue molecule.

    After Connect(R1) the propane residue becomes residue 2; FindLinkerAtoms
    picks a free H-bearing atom from that residue for the second connection.
    """
    R1 = Compound(propane) | FindLinkerAtoms()
    R2 = Compound(benzene) | FindLinkerAtoms()

    pipeline = (
        Compound(ethane)
        | FindLinkerAtoms()
        | Connect(R1)
        | SetAttachResidue(2)
        | FindLinkerAtoms()
        | Connect(R2)
    )

    ctx = pipeline()
    assert ctx.molecule is not None
    assert ctx.molecule.count_residues() == 3


# ── Optimization ──────────────────────────────────────────────────────


def test_choice_param_bounds(ethane, propane):
    block = Choice([ethane, propane])
    bounds = block.param_bounds
    assert isinstance(bounds, list)
    assert len(bounds) == 1
    lo, hi = bounds[0]
    assert lo == 0 and hi == 1


def test_compound_param_bounds(ethane):
    block = Compound(ethane)
    bounds = block.param_bounds
    assert isinstance(bounds, list)
    assert len(bounds) == 0


def test_chained_block_param_bounds(ethane, propane):
    block1 = (
        Choice([ethane, propane]) | FindLinkerAtoms() | Random(Modify(bam.acetylate))
    )
    bounds = block1.param_bounds
    assert isinstance(bounds, list)
    assert len(bounds) == 3
    lo, hi = bounds[0]
    assert lo == 0 and hi == 1
    lo, hi = bounds[1]
    assert lo == 0 and hi == len(FindLinkerAtoms._deterministic_how) - 1
    lo, hi = bounds[2]
    assert lo == 0 and hi == 1


def test_can_optimize_choice(ethane, propane):

    def score(mol):
        return mol.mass

    pipeline = Choice([ethane, propane])
    optimizer = Optimizer(pipeline, score)
    optimizer.run(strategy="random", steps=40)
    top_results = optimizer.top(n=1)
    assert len(top_results) == 1
    assert top_results[0].to_smiles() == propane.to_smiles()  # heavier than ethane

    def score2(mol):
        return -mol.mass

    optimizer = Optimizer(pipeline, score2)
    optimizer.run(strategy="random", steps=40)
    top_results = optimizer.top(n=1)
    assert len(top_results) == 1
    assert top_results[0].to_smiles() == ethane.to_smiles()


def test_can_optimize_findlinkeratoms():

    A = bam.read_smiles("CC(O)CC(N)CC")
    B = bam.read_smiles("CC(O)(O)CC(O)CC")

    def score(mol):
        N = mol.get_atoms("N", by="element")[0]
        Os_next_to_N = N.get_neighbors(filter=lambda a: a.element == "O")
        Os = mol.get_atoms("O", by="element")
        g = mol.get_atom_graph()
        d = sum(g.distance(N, o) for o in Os)
        return len(Os_next_to_N) - 0.1 * d

    pipeline = (
        Compound(A) | FindLinkerAtoms() | Connect(Compound(B) | FindLinkerAtoms())
    )
    optimizer = Optimizer(pipeline, score)
    optimizer.run(strategy="swarm", steps=30, population=3)
    top_results = optimizer.top(n=10)

    smi = "CC[C@H](O)C[C@@](C)(O)ON[C@@H](CC)C[C@H](C)O"
    assert any(mol.to_smiles() == smi for mol in top_results)

    # top = top_results[0]
    # assert (
    #     score(top) == -0.3
    # )  # corresponds to the optimal solution with 1 O next to N and a total distance of 3 to all Os

    # import matplotlib.pyplot as plt

    # fig, axs = plt.subplots(2, 5, figsize=(12, 6))
    # axs = axs.flatten()
    # for i, mol in enumerate(top_results):
    #     axs[i].imshow(mol.draw2d().draw())
    #     axs[i].set_title(f"Score: {score(mol):.2f}")
    #     axs[i].axis("off")

    # plt.tight_layout()
    # plt.show()


def test_can_optimize_random_modify(ethane):

    def score(mol):
        return mol.mass

    pipeline = Compound(ethane) | FindLinkerAtoms() | Random(Modify(bam.acetylate))

    Op = Optimizer(pipeline, score)
    Op.run(strategy="random", steps=10)
    top_results = Op.top(n=5)

    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 5, figsize=(12, 3))
    for i, mol in enumerate(top_results):
        axs[i].imshow(mol.draw2d().draw())
        axs[i].set_title(f"Score: {score(mol):.2f}")
        axs[i].axis("off")

    fig.tight_layout()
    plt.show()
