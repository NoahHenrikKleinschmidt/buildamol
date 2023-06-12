"""
The `Molecule` class is a wrapper around a biopython structure and a core part 
of biobuild functionality. It provides a convenient interface to molecular structures
and their properties, such as atoms, bonds, residues, chains, etc.


Making Molecules
================

To make a new molecule, the easiest way is to use the `Molecule.from_compound` classmethod.
This method is linked to the PDBECompounds database, which contains a large number of reference
compounds (currently restricted to sugars and amino acids due to the nature of "biobuild"). The database
is queried using the `by` parameter, which can be one of the following:
- "id" for the PDB id (default)
- "name" for the name of the compound (must match any known synonym of the iupac name)
- "formula" for the chemical formula (usually ambiguous and will therefore often raise an error)
- "smiles" for the SMILES string (also accepts InChI)

.. code-block:: python

    from biobuild import Molecule

    # create a new galactose molecule
    gal = Molecule.from_compound("GAL") # use the PDB id

    # create a new glucose molecule
    glc = Molecule.from_compound("alpha-d-glucose", by="name") # use the name

Alternatively, an existing structure can be loaded to a molecule from a PDB file using the `Molecule.from_pdb` classmethod. 
In this case the molecule may be as large as the PDB format supports (currently 99999 atoms).

.. code-block:: python

    my_glycan = Molecule.from_pdb("my_glycan.pdb")

Unless loaded from a PDB file, each `Molecule` starts with a single chain and residue,
specifying a single compound, such as Galactose or Glucose. Starting from there, multiple molecules
can be assembled together to form complex structures. 

Of course, an existing `biopython.Structure` object can also be used to create a `Molecule`:

.. code-block:: python

    from Bio.PDB import PDBParser

    parser = PDBParser()
    structure = parser.get_structure("my_structure", "my_structure.pdb")

    # do any modifications to the structure here
    # ...
    
    my_molecule = Molecule(structure)

    
If a molecule is created from a large PDB file and has undergone a lot of preprocessing, it may be useful to `save` the molecule to a pickle file
which can be loaded again later to avoid repeated preprocessing. This can be done using the `save` and `load` methods.

.. code-block:: python

    # save the molecule to a pickle file
    my_molecule.save("my_molecule.pkl")

    # load the molecule from the pickle file
    my_molecule = Molecule.load("my_molecule.pkl")

    

Modifying Molecules
===================

Once a molecule is created, it can be modified in a number of ways. The most common modifications are
- adding bonds (because when loading from a PDB file, bonds are not inferred automatically!)
- adding additional residues and atoms
- removing residues and atoms
- adjusting labelling (e.g. changing the chain names and residue seqids)

Adding bonds
------------

To add a bond between two atoms, we use the `add_bond` method, which expects two arguments for the two connected atoms.
The atoms can be specified by their `full_id` tuple, their `id` string, their `serial_number` (always starting at 1) or directly (the `biopython.Atom` object).

.. code-block:: python

    glc = Molecule.from_compound("GLC")

    # add a bond between the first and second atom
    glc.add_bond(1, 2)

    # and also add a bond between "O1" and "HO1" atoms
    glc.add_bond("O1", "HO1")


Already for small molecules such as bbucose with only 24 atoms, it would be very tedious to add all bonds manually. Good thing that the molecules
created using `from_compound` already contain all the default bonds!

However, for molecules that were loaded from a PDB file , the `Molecule` class offers two methods: `apply_standard_bonds` and `infer_bonds`. 
The former uses reference connectivity information from the PDBECompounds database to add all bonds that are known for the compound (if it exists in the database).
The latter will use a simple distance-based approach to infer bonds between atoms that are closer than a specified threshold (default: 1.6 Å), which can be restricted
further to a min-max window.

.. code-block:: python

    # add all standard bonds for Glucose
    glc.apply_standard_bonds()

    # add all bonds that are closer than 1.6 Å
    glc.infer_bonds(bond_length=1.6)

    # add all bonds that are closer than 1.6 Å, but not closer than 1.0 Å
    glc.infer_bonds(bond_length=(1.0, 1.6))


.. note::

    By default `infer_bonds` will not attempt to add bonds between atoms that belong to different residues.
    This is because in cases of suboptimal conformations or in very large structures atoms that are close in space may not be connected in the structure.
    To override this behaviour, set the `restrict_residues` parameter to `False`.

    .. code-block:: python

        # add all bonds that are closer than 1.6 Å, even if they belong to different residues
        glc.infer_bonds(bond_length=1.6, restrict_residues=False)


    Alternatively, instead of `infer_bonds` one may use `infer_residue_connections` to get bonds that connect different residues.
    This method will infer all bonds between atoms from different residues based on the distance between them. The inferred bonds are
    saved in the `Molecule` and also returned in a `list`. If the optional argument `triplet` is set to `True`, the methodd will also
    return the bonds immediately adjacent to the inferred bonds.

    Take the following example of a molecule with two residues `A` and `B` that are connected by a bond between `OA` and `(1)CB`:
    ```     
             connection -->  OA     OB --- H
                            /  \\   /
              (1)CA --- (2)CA   (1)CB 
               /         \\        \\
           (6)CA         (3)CA    (2)CB --- (3)CB
               \\         /
             (5)CA --- (4)CA
    ``` 
    
    If `triplet=False` the method will only return the bond between `OA` and `(1)CB`. However, if `triplet=True` it will also return the bond
    between `(2)CA` and `OA` - thus forming a triplet of atoms `(2)CA`, `OA` and `(1)CB` that connect the two residues A and B.


    .. code-block:: python

        # infer all bonds between atoms from different residues
        inferred_bonds = glc.infer_residue_connections()
        


Adding residues and atoms
-------------------------

To add one or more new residue(s), we use the `add_residues` method, which expects a number `biopython.Residue` objects as unnamed arguments.
Similarly, to add one or more new atom(s), we use the `add_atoms` method, which expects a number of `biopython.Atom` objects as unnamed arguments.
Both methods allow to specify the parent object (`chain` or `residue`) via an optional argument and will automatically choose the last
chain or residue if none is specified.

.. code-block:: python

    from Bio.PDB import Residue, Atom

    new_residue = Residue("XYZ", " ", 1, " ")
    
    # do things with the residue here
    # ...

    # add the residue to the molecule
    # (add it to the last chain, whichever that may be)
    glc.add_residues(new_residue)

    new_atom = Atom("X", [0, 0, 0], 0, 0, "X", "X", 0, "X")

    # add the atom to the first residue in the `glc` molecule
    glc.add_atoms(new_atom, residue=1)

Removing residues and atoms
---------------------------

In order to remove residues or atoms or bonds, we can use the `remove_residues`, `remove_atoms` and `remove_bond`(yes, singluar!) mathods.
They work exactly like their `add_` counterparts, but instead of adding, they remove the specified objects.

.. code-block:: python

    # remove the first residue
    glc.remove_residues(1)

    # remove the first atom
    glc.remove_atoms(1)

    # remove the bond between the first and second atom
    glc.remove_bond(1, 2)


.. warning::

    When adding and removing atoms and residues, the `Molecule` object will not automatically update the `biopython.Structure` object as well as its internal
    connectivity graph (the `AtomGraph`). However, if the user chooses to edit the biopython structure directly, the `AtomGraph` will **not** be updated automatically!
    In this case, the user must call the `update_atom_graph` method to update the `AtomGraph` manually. 

    .. code-block:: python

        # add a new residue to the molecule
        glc.add_residues(new_residue)

        # now add atoms into the residue directly instead of via `add_atoms`
        new_residue.add(new_atom)
    
        # now update the graph
        glc.update_atom_graph()


Adjusting labelling
-------------------

Single-residue molecules that were loaded from a PDB file may not use the same atom labelling as the `PDBE` and `CHARMM` databases.
In order to quickly adjust the labelling, a method `autolabel` exists. However, it is restricted to default-conformation residues at the moment
as it uses a structure overlay to determine the correct labelling. If your molecules follow another labelling scheme it may be more efficient to simply
define your own linkage recipies or patches (see the documentation of `linkages`).

.. code-block:: python

    # load a molecule from a PDB file
    glc = Molecule.from_pdb("glucose.pdb")

    # adjust the labelling
    glc.autolabel()

    # save the molecule to a new PDB file
    glc.to_pdb("glucose_autolabelled.pdb")


A more common operation is the adjustment of chain names and residue seqids. This can be done using the `reindex` method.
This method accepts three starting values for the chain name, residue seqid and atom serial number and will then reindex all chains, residues and atoms
to ensure they are continuously numbered and labelled. Some internal methods used when connecting different molecules are reliant on a continuous
numbering scheme, so this method should be called before connecting molecules.

.. code-block:: python

    # load a molecule from a PDB file
    glc = Molecule.from_pdb("glucose.pdb")

    # reindex the molecule
    glc.reindex()

We can also use one molecule as a "reference" for reindexing another molecule to make sure there are now labelling conflicts between them in case we want to connect them together later (this is usually done internally by _biobuild_ automatically).

.. code-block:: python  

    # load a molecule from a PDB file
    glc = Molecule.from_pdb("glucose.pdb")

    # load another molecule from a PDB file
    cel = Molecule.from_pdb("cellulose.pdb")
    cel.reindex() # make sure we have a continuous numbering scheme

    # reindex the glucose molecule using the cellulose molecule as a reference
    cel.adjust_indexing(glc)



Connecting Molecules
====================

Since most modifications are not simply single residues but rather complex structures, the second main purpose of a `Molecule` is to be easily 
connected to other Molecules to form a larger structure. To this end, the `Molecule` class provides a number of methods to easily assemble complex structures from small single residue molecules.


Forming Polymers
----------------

The simplest way to generate a large structure is probably the `repeat` method, which will repeat the given molecule `n` times
to form a homo-polymer.

.. code-block:: python

    # create a glucose molecule
    glc = Molecule.from_compound("GLC")

    # create cellulose from glucose
    # using a 1-4 beta-beta glycosidic linkage
    glc.repeat(10, "14bb")

    # Now we have a cellulose of 10 glucoses

In the above example we used the `repeat` method explicitly, but we could also achieve the same with the short-hand `*=`. For this to work, we need to specify 
the linkage type beforehand. We do this by setting the `patch` attribute before using any operator.

.. code-block:: python

    # specify the "default" linkage type for connecting 
    # other molecules to this glucose
    glc.patch = "14bb"

    # now make a cellulose by multiplying glucoses
    glc *= 20

    # Now we have a cellulose of 20 glucoses

    
If we wish to keep `glc` as a single residue Glucose and still get our desired cellulose, 
we can set `inplace=False` when calling `repeat` or simply use the `*` operator, both of which 
will have the same effect of creating a new copy that houses the appropriate residues.

.. code-block:: python

    cel = glc.repeat(10, "14bb", inplace=False)

    # or (equivalently)
    glc.patch = "14bb"
    cel = glc * 10


    
Connecting different Molecules
------------------------------

What if we want to connect different molecules? For example, we may want to connect a Galactose to a Glucose to form Lactose.
This can be achieved using the `attach` method, which will attach a given molecule to to another molecule.

.. code-block:: python

    glc = Molecule.from_compound("GLC")
    gal = Molecule.from_compound("GAL")

    # attach the galactose to the glucose
    # (we want a copy, so we set inplace=False just like with 'repeat')
    lac = glc.attach(gal, "14bb", inplace=False)

    # Now we have a lactose molecule


In the above example, the `attach` method is used to attach the galactose molecule to the glucose, but for those
among us who prefer a more shorty syntax, the `+` operator will do the same thing. 

.. code-block:: python

    # specify that incoming molecules shall be
    # attached using a 1-4 beta linkage
    glc.patch = "14bb"

    # now attach the galactose to the glucose
    lac = glc + gal


Of course, if there is a `+` operator there should also be a `+=` operator, which is simply the equivalent of `attach` with `inplace=True`.

.. code-block:: python

    glc.patch = "14bb"
    glc += gal


    # Now 'glc' is a lactose molecule


Setting default Modifiers
-------------------------

So far, we have always worked with a 1-4 beta-beta glycosidic linkage, which we apparently could select using the string `"14bb"`. 
But what if we want to use a different linkage type? For example, a 1-4 alpha-beta glycosidic linkage?  You of course noticed, 
that `attach` and `repeat` accept an argument `recipe_or_patch` which allows you to specify the linkage type, and that if you leave it blank
the default linkage type is used. But how do we set the default linkage type?

Let's first check what linkage types are available by default anyway. Have you noticed an argument named `_topology`
at the end of the `attach` or `repeat` methods? The `topology` refers to the underlying _CHARMM_ topology which hosts the linkage type information.
By default a topology is already loaded in _biobuild_'s framework so it is not necessary for the user to specify anything here, but we can check which linkage types are available by:

.. code-block:: python

    import biobuild as bb

    # get the default topology
    topology = bb.get_default_topology()

    print(topology.patches)


Any of these linkages can be referenced by their name, e.g. `"14bb"` or `"14ab"`. 

Wait a second, my desired linkage is not in the list! What now?! Well, you can always define your own linkage type by creating either a new `Patch` or a new `Recipe`.
Check out the documentation on :ref:`linkages` for more information on how to do this. If you have your desired patch or recipe
ready to go, set it as the default patch or recipe by:

.. code-block:: python

    my_molecule.patch = my_patch
    # or
    my_molecule.recipe = my_recipe

    # or if you feel "super correct"
    my_molecule.set_ linkage(my_ linkage)

    # or if you feel "extra cocky"
    my_molecule % my_ linkage # <- the modulo operator assignes the "modifier" to the molecule


Now any call to `attach`, `repeat`, or any of its operator proxies will use your defined linkage by default.

Setting the default Residue for attachment
------------------------------------------

When specifying a `Patch` or `Recipe` we specify which atoms are supposed to be connected and removed, but we do not specify which residues
these atoms belong to. We specify this as arguments inside the `attach` method for instance, but we can also leave this blank, in which case the
last residue in the molecule is used by default. This is obviously not always what we want, however! Hence, if we do not want to specify the residue for attachment
at every `attach` call or if we want to use the `+` operator, we can set the default residue for attachment by setting the `attach_residue` attribute:

.. code-block:: python

    # set the default attachment residue to the first residue
    my_molecule.attach_residue = 1

    # or
    my_molecule.set_attach_residue(1)

    # or (if you feel "extra cocky")
    my_molecule @ 1 # <- the 'at' operator sets the attachment residue


    # serial number indexing also works in reverse
    # (set the second last residue as the default attachment residue)
    my_molecule.attach_residue = -2

Attaching Molecules to a Scaffold
=================================

Ultimately, the fate of a `Molecule` is to be attached to a `Scaffold`.
The `Scaffold` class is very similar to the `Molecule` class but instead of representing modifying molecules such as bbycans, it 
represents the structures which are getting modified, such as proteins or membranes. A `Molecule` is always attached to a `Scaffold` 
at its (the molecule's) `root_atom`. This can be set using the `root_atom` argument when creating a `Molecule` or by setting the `root_atom` attribute later on.

.. code-block:: python

    # set the root atom to the atom with serial number 1
    my_molecule.root_atom = 1

    # or
    my_molecule.set_root(1)

    # or (if you fell "extra cocky")
    my_molecule ^ 1 # <- the caret operator sets the root atom 

    # now we can attach the molecule to a scaffold
    # ... some scaffold creation code ...

    my_scaffold.attach(my_molecule)
    
Now, we can use the `Molecule` to attach it to a `Scaffold` (described in the :ref:`Scaffold` documentation).
"""

from copy import deepcopy

from typing import Union

import numpy as np
import Bio.PDB as bio

import biobuild.core.entity as entity
import biobuild.core.Linkage as Linkage
import biobuild.utils as utils
import biobuild.structural as structural
import biobuild.resources as resources


class Molecule(entity.BaseEntity):
    """
    A molecule to add onto a scaffold.
    A molecule consists of a single chain.

    Parameters
    ----------
    structure : bio.PDB.Structure
        A biopython structure object
    root_atom : str or int or bio.PDB.Atom
        The id or the serial number of the root atom
        at which the molecule would be attached to a another
        structure such as protein scaffold or another Molecule.
    model : int
        The model to use from the structure. Defaults to 0. This may be any
        valid identifier for a model in the structure, such as an integer or string.
    chain : str
        The chain to use from the structure. Defaults to the first chain in the structure.
    """

    def __init__(
        self,
        structure,
        root_atom: Union[str, int, bio.Atom.Atom] = None,
        model: int = 0,
        chain: str = None,
    ):
        super().__init__(structure, model)

        if not chain or len(self._model.child_list) == 1:
            self._working_chain = self._model.child_list[0]
        else:
            self._working_chain = self._model.child_dict.get(chain)
            if not self._working_chain:
                raise ValueError("The chain {} is not in the structure".format(chain))

        if root_atom:
            self.set_root(root_atom)
            if not root_atom in self._working_chain.get_atoms():
                raise ValueError("The root atom is not in the structure")
        self._root_atom = root_atom

    @classmethod
    def empty(cls, id: str = None) -> "Molecule":
        """
        Create an empty Molecule object

        Parameters
        ----------
        id : str
            The id of the Molecule. By default an id is inferred from the filename.

        Returns
        -------
        Molecule
            An empty Molecule object
        """
        structure = structural.make_empty_structure(id)
        return cls(structure)

    @classmethod
    def from_pdb(
        cls,
        filename: str,
        root_atom: Union[str, int] = None,
        id: str = None,
        model: int = 0,
        chain: str = None,
    ) -> "Molecule":
        """
        Read a Molecule from a PDB file

        Parameters
        ----------
        filename : str
            Path to the PDB file
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        model : int
            The model to use from the structure. Defaults to 0. This may be any
            valid identifier for a model in the structure, such as an integer or string.
        chain : str
            The chain to use from the structure. Defaults to the first chain in the structure.

        Returns
        -------
        Molecule
            The Molecule object
        """
        if id is None:
            id = utils.filename_to_id(filename)
        struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
        new = cls(struct, root_atom, model=model, chain=chain)
        bonds = utils.pdb.parse_connect_lines(filename)
        if len(bonds) != 0:
            for b in bonds:
                new.add_bond(*b)
        return new

    @classmethod
    def from_compound(
        cls,
        compound: str,
        by: str = "id",
        root_atom: Union[str, int] = None,
    ) -> "Molecule":
        """
        Create a Molecule from a reference compound from the PDBECompounds database

        Parameters
        ----------
        compound : str
            The compound to search for
        by : str
            The field to search by. This can be
            - "id" for the PDB id
            - "name" for the name of the compound (must match any known synonym of the iupac name)
            - "formula" for the chemical formula
            - "smiles" for the SMILES string (also accepts InChI)
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        """
        mol = resources.get_default_compounds().get(compound, by=by)
        if isinstance(mol, list):
            raise ValueError(
                f"Multiple compounds found using '{by}={compound}', choose any of these ids specifically {[i.id for i in mol]}"
            )
        mol.set_root(root_atom)
        return mol

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        id: str = None,
        root_atom: Union[str, int] = None,
        add_hydrogens: bool = True,
    ) -> "Molecule":
        """
        Read a Molecule from a SMILES string

        Parameters
        ----------
        smiles : str
            The SMILES string
        id : str
            The id of the Molecule. By default the provided smiles string is used.
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        add_hydrogens : bool
            Whether to add hydrogens to the molecule

        Returns
        -------
        Molecule
            The Molecule object
        """
        struct = structural.read_smiles(smiles, add_hydrogens)
        struct.id = id if id else smiles
        return cls(struct, root_atom)

    @classmethod
    def from_pubchem(
        cls,
        query: str,
        root_atom: Union[str, int] = None,
        by: str = "name",
        idx: int = 0,
    ) -> "Molecule":
        """
        Create a Molecule from PubChem

        Note
        ----
        PubChem follows a different atom labelling scheme than the CHARMM force field!
        This means that atom names may not match the names required by the default patches that are
        integrated in biobuild. It is advisable to run `autolabel` or `relabel_hydrogens` on the molecule.
        Naturally, custom patches or recipies working with adjusted atom names will always work.

        Parameters
        ----------
        query : str
            The query to search for in the PubChem database
        root_atom : str or int
            The id or the serial number of the root atom (optional)
        by : str
            The method to search by. This can be any of the following:
            - cid
            - name
            - smiles
            - sdf
            - inchi
            - inchikey
            - formula
        idx : int
            The index of the result to use if multiple are found. By default, the first result is used.

        Returns
        -------
        Molecule
            The Molecule object
        """
        _compound_2d, _compound_3d = resources.pubchem.query(query, by=by, idx=idx)
        new = _molecule_from_pubchem(_compound_2d.iupac_name, _compound_3d)
        new.set_root(root_atom)
        return new

    def get_residue_connections(self, triplet: bool = True, direct_by: str = "resid"):
        """
        Get bonds between atoms that connect different residues in the structure
        This method is different from `infer_residue_connections` in that it works
        with the already present bonds in the molecule instead of computing new ones.

        Parameters
        ----------
        triplet : bool
            Whether to include bonds between atoms that are in the same residue
            but neighboring a bond that connects different residues. This is useful
            for residues that have a side chain that is connected to the main chain.
            This is mostly useful if you intend to use the returned list for some purpose,
            because the additionally returned bonds are already present in the structure
            from inference or standard-bond applying and therefore do not actually add any
            particular information to the Molecule object itself.
        direct_by : str
            The attribute to sort by. Can be either "serial", "resid" or "root".
            In the case of "serial", the bonds are sorted by the serial number of the first atom.
            In the case of "resid", the bonds are sorted by the residue id of the first atom.
            In the case of "root", the bonds are sorted by the root atom of the first atom.
            Set to None to not sort the bonds.

        Returns
        -------
        set
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = super().get_residue_connections(triplet)
        if direct_by is not None:
            direct_connections = None
            if direct_by == "resid":
                direct_connections = self.get_residue_connections(
                    triplet=False, direct_by=None
                )
                direct_connections1, direct_connections2 = set(
                    a for a, b in direct_connections
                ), set(b for a, b in direct_connections)
                direct_connections = direct_connections1.union(direct_connections2)
            bonds = self._direct_bonds(bonds, direct_by, direct_connections)
        return set(bonds)

    def repeat(self, n: int, link=None, inplace: bool = True):
        """
        Repeat the molecule n times into a homo-polymer.

        Parameters
        ----------
        n : int
            The number or units of the final polymer.
        link : str or Patch or Recipe
            The patch or recipe to use when patching individual units together.
            If noe is given, the default patch or recipe is used (if defined).
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.

        Returns
        -------
        molecule
            The modified molecule (either the original object or a copy)
        """
        if not isinstance(n, int):
            raise TypeError("Can only multiply a molecule by an integer")
        if n <= 0:
            raise ValueError("Can only multiply a molecule by a positive integer")
        if not link and not self._linkage:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        _other = deepcopy(self)
        if not inplace:
            obj = deepcopy(self)
        else:
            obj = self

        _patch = obj._linkage
        if link:
            obj % link

        for i in range(n - 1):
            obj += _other

        if _patch:
            obj % _patch

        return obj

    def attach(
        self,
        other: "Molecule",
        link: Union[str, "Linkage"] = None,
        at_residue: Union[int, bio.Residue.Residue] = None,
        other_residue: Union[int, bio.Residue.Residue] = None,
        inplace: bool = True,
        _topology=None,
    ):
        """
        Attach another structure to this one using a Patch or a Recipe.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        link : str or Patch or Recipe
            Either a Patch to apply when attaching or a Recipe to use when stitching.
            If None is defined, the default patch or recipe that was set earlier on the molecule is used.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the defined `attach_residue` is used.
        other_residue : int or Residue
            The residue in the other molecule to attach this molecule to. If None, the defined `attach_residue` of the other molecule is used.
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.
        _topology : Topology
            The topology to use when attaching. If None, the topology of the molecule is used. Only used if the patch is a string.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only attach a Molecule to another Molecule")

        if not inplace:
            obj = deepcopy(self)
        else:
            obj = self

        if not link:
            link = obj._linkage
            if not link:
                raise ValueError("Cannot attach a molecule without a patch defined")

        if isinstance(link, str):
            if not _topology:
                _topology = resources.get_default_topology()
            link = _topology.get_patch(link)

        if link.has_IC:
            obj.patch_attach(
                other,
                link,
                at_residue=at_residue,
                other_residue=other_residue,
                _topology=_topology,
            )
        else:
            obj.stitch_attach(
                other,
                link,
                at_residue=at_residue,
                other_residue=other_residue,
            )
        return obj

    def patch_attach(
        self,
        other: "Molecule",
        patch: Union["Linkage", str] = None,
        at_residue: Union[int, bio.Residue.Residue] = None,
        other_residue: Union[int, bio.Residue.Residue] = None,
        _topology=None,
    ):
        """
        Attach another structure to this one using a Patch.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        patch : str or Linkage
            A linkage to apply when attaching. If none is given, the default link that was set earlier
            on the molecule is used. If no patch was set, an AttributeError is raised. If a string is
            given, it is interpreted as the name of a patch in the topology.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the last residue of the molecule.
        other_residue : int or Residue
            The residue of the other molecule to attach. If None, the first residue of the other molecule.
        _topology
            A specific topology to use for referencing.
            If None, the default CHARMM topology is used.
        """

        if not patch:
            patch = self._linkage
        if not patch:
            raise AttributeError(
                "No patch was set for this molecule. Either set a default patch or provide a patch when attaching."
            )
        if isinstance(patch, str):
            if not _topology:
                _topology = resources.get_default_topology()
            patch = _topology.get_patch(patch)
        if not patch:
            raise ValueError(
                "No patch was found with the given name. Either set a default patch or provide a patch when attaching."
            )

        p = structural.__default_keep_copy_patcher__
        p.apply(patch, self, other, at_residue, other_residue)
        p.merge()
        return self

    def stitch_attach(
        self,
        other: "Molecule",
        recipe: Linkage = None,
        remove_atoms=None,
        other_remove_atoms=None,
        at_atom=None,
        other_at_atom=None,
        at_residue=None,
        other_residue=None,
    ):
        """
        Stitch two molecules together by removing atoms and connecting them with a bond. This works without a pre-defined patch.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        recipe : Recipe
            The recipe to use when stitching. If None, the default recipe that was set earlier on the molecule is used (if defined).
        remove_atoms : list of int
            The atoms to remove from this molecule. Only used if no recipe is provided.
        other_remove_atoms : list of int
            The atoms to remove from the other molecule. Only used if no recipe is provided.
        at_atom : int or str or Bio.PDB.Atom
            The atom forming the bond in this molecule. If a string is provided, an `at_residue` needs to be defined from which to get the atom. If None is provided, the root atom is used (if defined). Only used if no recipe is provided.
        other_at_atom : int or str or Bio.PDB.Atom
            The atom to attach to in the other molecule. If a string is provided, an `other_residue` needs to be defined from which to get the atom. If None is provided, the root atom is used (if defined). Only used if no recipe is provided.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the `attach_residue` is used. Only used if a recipe is provided and the atoms
        """
        if not recipe and not remove_atoms:
            if not self._linkage:
                raise AttributeError(
                    "No recipe was set for this molecule and no manual instructions were found. Either set a default recipe, provide a recipe when stitching, or provide the information about removed and bonded atoms directly."
                )
            recipe = self._linkage

        if recipe:
            return self.stitch_attach(
                other,
                remove_atoms=recipe.deletes[0],
                other_remove_atoms=recipe.deletes[1],
                at_atom=recipe.bonds[0][0],
                other_at_atom=recipe.bonds[0][1],
                at_residue=at_residue,
                other_residue=other_residue,
            )

        if not isinstance(other, Molecule):
            raise TypeError("Can only stitch two molecules together")

        if not at_atom:
            at_atom = self._root_atom
            if not at_atom:
                raise ValueError(
                    "No atom to attach to was provided and no root atom was defined in this molecule"
                )

        if not other_at_atom:
            other_at_atom = other._root_atom
            if not other_at_atom:
                raise ValueError(
                    "No atom to attach to was provided and no root atom was defined in the other molecule"
                )

        # at_atom = self.get_atom(at_atom, residue=at_residue)
        # other_at_atom = other.get_atom(other_at_atom, residue=other_residue)

        p = structural.__default_keep_copy_stitcher__
        p.apply(self, other, remove_atoms, other_remove_atoms, at_atom, other_at_atom)
        self = p.merge()
        return self

    def __add__(self, other) -> "Molecule":
        """
        Add two molecules together. This will return a new molecule.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only add two molecules together")

        patch = self._linkage
        if not patch:
            patch = other._linkage
        if not patch:
            raise RuntimeError(
                "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 + mol2)"
            )

        if patch.has_IC:
            p = structural.__default_copy_copy_patcher__
            p.apply(patch, self, other)
        else:
            p = structural.__default_copy_copy_stitcher__
            p.apply(
                self,
                other,
                patch.deletes[0],
                patch.deletes[1],
                patch.bonds[0][0],
                patch.bonds[0][1],
            )
        new = p.merge()
        return new

    def __iadd__(self, other) -> "Molecule":
        """
        Attach another molecule to this one
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only add two molecules together")

        patch = self._linkage
        if not patch:
            patch = other._linkage
        if not patch:
            raise RuntimeError(
                "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 += mol2)"
            )

        self.attach(other, patch)
        return self

    def __mul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the * operator (i.e. mol * 3)
        This requires that the molecule has a patch defined
        """
        if not isinstance(n, int):
            raise TypeError("Can only multiply a molecule by an integer")
        if n <= 0:
            raise ValueError("Can only multiply a molecule by a positive integer")
        if not self._linkage:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        new = self + self
        for i in range(n - 2):
            new += self

        return new

    def __imul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the *= operator (i.e. mol *= 3)
        This requires that the molecule has a patch defined
        """
        if not self._linkage:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        self.repeat(n)

        return self

    def __repr__(self):
        return f"Molecule({self.id})"


def _molecule_from_pubchem(id, comp):
    """
    Convert a PubChem compound to a Molecule.

    Parameters
    ----------
    id : str
        The id of the molecule
    comp : pcp.Compound
        The PubChem compound

    Returns
    -------
    Molecule
        The molecule
    """
    bonds = comp.bonds
    structure = bio.Structure.Structure(id)
    model = bio.Model.Model(0)
    structure.add(model)
    chain = bio.Chain.Chain("A")
    model.add(chain)
    residue = bio.Residue.Residue((" ", 1, " "), "UNL", 1)
    chain.add(residue)

    element_counts = {}
    for atom in comp.atoms:
        element = atom.element
        element_counts[element] = element_counts.get(element, 0) + 1

        _atom = bio.Atom.Atom(
            f"{element}{element_counts[element]}",
            np.array((atom.x, atom.y, atom.z)),
            1.0,
            1.0,
            0,
            element,
            atom.aid,
            element=element,
        )
        residue.add(_atom)

    mol = Molecule(structure)
    for bond in bonds:
        mol.add_bond(bond.aid1, bond.aid2)

    return mol


if __name__ == "__main__":
    recipe = Linkage()
    recipe.add_delete("O1", "target")
    recipe.add_delete("HO1", "target")
    recipe.add_delete("HO4", "source")
    recipe.add_bond(("C1", "O4"))

    glc = Molecule.from_compound("GLC")
    glc.set_linkage(recipe)

    glc2 = deepcopy(glc)

    _current_residues = len(glc.residues)
    glc3 = glc + glc2
    glc3.show()
    # glc3 = glc % "14bb" + glc
    pass
    # import pickle

    # tmp = pickle.load(open("/Users/noahhk/GIT/biobuild/tmp.pickle", "rb"))
    # tmp.get_residue_connections()
    # exit()
    # from timeit import timeit

    # # man = Molecule.from_compound("MAN")
    # glc = Molecule.from_compound("GLC")

    # # man.adjust_indexing(glc)
    # # assert glc.atoms[0].serial_number == len(man.atoms) + 1

    # t1 = timeit()

    # glc.repeat(5, "14bb")
    # # glc % "14bb"
    # # glc *= 10

    # t2 = timeit()

    # print(t2 - t1)

    from biobuild.utils import visual

    # v = visual.MoleculeViewer3D(glc)
    # v.show()

    # man = Molecule.from_pdb("support/examples/membrane.pdb", model=4)
    # man.infer_bonds()
    # man.

    man = Molecule.from_pdb("/Users/noahhk/GIT/biobuild/support/examples/MAN9.pdb")
    man.infer_bonds(restrict_residues=False)

    g = man.make_residue_graph(True)
    v = visual.MoleculeViewer3D(g)

    for c in man.get_residue_connections():
        v.draw_vector(
            None, c[0].coord, 1.3 * (c[1].coord - c[0].coord), color="magenta"
        )
    v.show()
