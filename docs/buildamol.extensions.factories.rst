.. _factories:

Molecular Factories
===================

The `molecular_factories` package provides a collection of tools for the automatic generation of molecular structures.

Specifically implemented are:

.. tab-set::

    .. tab-item:: Fragment Libraries
        

        This module provides functions to split existing molecules into `BRICS fragments <https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cmdc.200800178>`_ or to query `ChEMBL <https://www.ebi.ac.uk/chembl/>`_ for common fragments in drug-like molecules.

        .. automodule:: buildamol.extensions.molecular_factories.fragment_library
            :members:
            :undoc-members:
            :show-inheritance:


    .. tab-item:: Factories

        The new factories provide a `Langchain <https://github.com/langchain-ai/langchain>`-inspired syntax for molecular generation. They offer freely combinable building blocks to direct a generation process,
        which can be coupled to a scoring function to optimize the generated molecules for a specific property.

        Usage
        -----
        For instance we could draft a generation process that will 
        graft some fragments onto a core molecule like so:

        .. code-block:: python
            
            import buildamol as bam
            from buildamol.extensions import molecular_factories as mf

            fragments = mf.get_chembl_fragments(n=10, max_heavy_atoms=5)

            core = bam.Molecule.from_smiles("c1ccccc1")  # benzene
            core.change_element("H1", "F")
            core.change_element("C3", "N")


            # let's say we want to graft one fragment opposite the fluorine and next to the nitrogen
            # and another smaller fragment next to the fluorine but opposite the nitrogen 
            def get_linker_a(mol):
                # first atom is linker the second defines the deleter(s)
                return mol.get_atom("C4"), None # could be anything more robust


            def get_linker_b(mol):
                return mol.get_atom("C6"), None # could be anything more robust
            

            # make a "larger" fragment by using two fragments and connecting them
            frag_a1 = mf.Choice(fragments) | mf.FindLinkerAtoms()
            frag_a2 = mf.Choice(fragments) | mf.FindLinkerAtoms()
            frag_a = frag_a2 | mf.Connect(frag_a1) | mf.FindLinkerAtoms()

            frag_b = mf.Choice(fragments) | mf.FindLinkerAtoms()

            factory = (
                mf.Compound(core) |
                mf.FindLinkerAtoms(get_linker_a) |
                mf.Connect(frag_a) |
                mf.SetAttachResidue(1) |
                mf.FindLinkerAtoms(get_linker_b) |
                mf.Connect(frag_b) |
                mf.Forge() # the final step that will compose the whole pipeline
            )


        To then sample from this factory we can simply call it like so:

        .. code-block:: python

            from buildamol.utils.visual import gallery_grid
                        
            fig, axs = gallery_grid(10, draw_molecules=False)
            for i, ax in enumerate(axs.flat):
                
                ctx = factory() # <-- this will run the factory and return a context object with the generated molecule
                
                mol = ctx.molecule
                mol.id = f"sample_{i}"
                mol.draw2d(ax=ax)
                ax.set_title(mol.id)
                
        
        .. image:: examples/files/factories_output_example1.png


        This so far blindly samples from the factors. We can also use the ``Optimizer`` class to steer the factory
        to optimise for certain properties. Let's make a simple case of optimising for QED, and making sure our "fragment A slightly bigger than fragment B" actually holds true, since we're currently not enforcing that in the pipeline. 


        .. code-block:: python

            from rdkit.Chem import QED
   
            def qed_score(mol):
                return QED.qed(mol.to_rdkit())


            def size_violation_penalty(mol):
                frag_a = len(mol.get_residue(2).atoms) + len(mol.get_residue(3).atoms)
                frag_b = len(mol.get_residue(4).atoms)
                return frag_b - frag_a  # negative if frag_a is bigger than frag_b, positive otherwise

            def score(mol):
                return qed_score(mol) - size_violation_penalty(mol)
            

            # now we can set up an optimiser to navigate our sampling space
            opt = mf.Optimizer(factory, scoring_fn=score)
            opt.run("random", n=100)
            opt.run("pso", n=50)

            topk_df = opt.to_dataframe().head(n=10)
            topk = opt.top(10)

            fig, axs = gallery_grid(topk)

            for i, ax in enumerate(axs.flat):
                if i >= len(topk):
                    break
                score = topk_df.iloc[i]["score"]
                ax.set_title(f"top {i+1} score: {score:.3f}")

        .. image:: examples/files/factories_output_example2.png

        
        .. dropdown:: Inputs
            
            .. automodule:: buildamol.extensions.molecular_factories.sources
                :members:
                :undoc-members:
                :show-inheritance:
        
        .. dropdown:: Actions & Modifiers

            .. automodule:: buildamol.extensions.molecular_factories.modifiers
                :members:
                :undoc-members:
                :show-inheritance:




        .. dropdown:: Optimizer

            .. automodule:: buildamol.extensions.molecular_factories.optimizer
                :members:
                :undoc-members:
                :show-inheritance:

        .. dropdown:: Base Blocks

            .. automodule:: buildamol.extensions.molecular_factories.basE
                :members:
                :undoc-members:
                :show-inheritance:
                
    .. tab-item:: Assembler

        .. automodule:: buildamol.extensions.molecular_factories.assembler
            :members:
            :undoc-members:
            :show-inheritance:

    .. tab-item:: Derivator
        
        .. automodule:: buildamol.extensions.molecular_factories.derivator
            :members:
            :undoc-members:
            :show-inheritance:

    .. tab-item:: Generator
        
        .. automodule:: buildamol.extensions.molecular_factories.generator
            :members:
            :undoc-members:
            :show-inheritance:
