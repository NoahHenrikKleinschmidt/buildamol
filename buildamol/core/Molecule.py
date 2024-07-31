"""
The `Molecule` class is a wrapper around a biopython structure and a core part 
of BuildAMol functionality. It provides a convenient interface to molecular structures
and their properties, such as atoms, bonds, residues, chains, etc. 

.. note:: 

    To help with identifying individual atoms, residues, etc. BuildAMol uses a different identification scheme than biopython does.
    Therefore BuildAMol comes with its own child classes of the biopython classes
    that are used to represent the structure. These classes are called `Atom`, `Residue`, `Chain`, etc.
    and can be used as drop-in replacements for the biopython classes and should not break any existing code.
    However, in case any incompatibility is observed anyway, the classes are equipped with a `to_biopython` method
    that will remove the BuildAMol overhead and return pure biopython objects (this is not an in-place operation, however,
    and will return a new object). 


Making Molecules
================

The easiest way to make a new molecule is to use the toplevel `molecule` function, which will automatically
try to detect the type of user provided input and generate a molecule from it. Currently supported inputs are:
- A biopython structure
- An RDKit molecule
- An OpenBabel molecule
- An OpenMM topology
- An STK molecule
- A PDB id
- A PDB file
- A CIF file
- A MOL file
- A JSON file
- An XML file
- A SMILES string
- An InChI string
- An IUPAC name or abbreviation, or any name that matches a known compound synonym that is associated with the PubChem database

.. code-block:: python

    from buildamol import molecule

    my_glucose = molecule("GLC") # use the PDB id

    # or
    my_glucose = molecule("GLC.pdb") # use the PDB file

    # or
    my_glucose = molecule("alpha-d-glucose") # use the name

    # ...

Since the `molecule` function is a try-and-error function it is convenient but not the most efficient. Hence, the `Molecule` class
offers already a number of convenient methods to easily generate molecules directly from specified data sources. Available methods are:

- `Molecule.from_pdb` to generate a molecule from a PDB file
- `Molecule.from_cif` to generate a molecule from a CIF file
- `Molecule.from_smiles` to generate a molecule from a SMILES string
- `Molecule.from_pubchem` to generate a molecule from a PubChem entry
- `Molecule.from_compound` to generate a molecule from a PDBECompounds entry
- `Molecule.from_rdkit` to generate a molecule from an RDKit molecule object
- `Molecule.from_openmm` to generate a molecule from an OpenMM topology object
- `Molecule.from_stk` to generate a molecule from an STK molecule object
- `Molecule.from_molfile` to generate a molecule from a MOL file
- `Molecule.from_json` to generate a molecule from a JSON file
- `Molecule.from_xml` to generate a molecule from an XML file
- `Molecule.empty` to generate an empty molecule


Hence, if we know that "glucose" is already available in our local PDBECompounds database, we can generate the molecule also as follows:

.. code-block:: python

    from buildamol import Molecule

    my_glucose = Molecule.from_compound("GLC") # use the PDB id


The quickest way to query the local PDBECompounds database is to use the `PDB Id` of the desired compounds.
However, the `from_compound` accepts other inputs as well. The database is queried using the `by` parameter, which can be one of the following:
- "id" for the PDB id (default)
- "name" for the name of the compound (must match any known synonym of the iupac name)
- "formula" for the chemical formula (usually ambiguous and will therefore often raise an error)
- "smiles" for the SMILES string (also accepts InChI)

.. code-block:: python

    # create a new glucose molecule
    glc = Molecule.from_compound("alpha-d-glucose", by="name") # use the name


Saving Molecules
----------------
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


Already for small molecules such as glucose with only 24 atoms, it would be very tedious to add all bonds manually. Good thing that the molecules
created using `from_compound` or `from_pdb` already contain all the default bonds!

However, in case the bonds are missing, or the PDB file did not specify any to begin with, the `Molecule` class offers two methods: `apply_standard_bonds` and `infer_bonds`. 
The former uses reference connectivity information from the PDBECompounds database or CHARMMTopology to add all bonds that are known for the compound (if it exists in the database).
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

        
Residue Connections
-------------------

Alternatively, instead of `infer_bonds` one may use `infer_residue_connections` to get bonds that connect different residues.
This method will infer all bonds between atoms from different residues based on the distance between them. The inferred bonds are
saved in the `Molecule` and also returned in a `list`. If the optional argument `triplet` is set to `True`, the methodd will also
return the bonds immediately adjacent to the inferred bonds.

Take the following example of a molecule with two residues `A` and `B` that are connected by a bond between `OA` and `(1)CB`:

.. code-block:: text    

     connection -->  OA     OB --- H
                    /  \\   /
     (1)CA --- (2)CA   (1)CB 
        /         \\        \\
    (6)CA         (3)CA    (2)CB --- (3)CB
        \\         /
     (5)CA --- (4)CA
 
     
If `triplet=False` the method will only return the bond between `OA` and `(1)CB`. However, if `triplet=True` it will also return the bond
between `(2)CA` and `OA` - thus forming a triplet of atoms `(2)CA`, `OA` and `(1)CB` that connect the two residues A and B.

.. code-block:: python

    # infer all bonds between atoms from different residues
    >>> glc.infer_residue_connections(triplet=False)
    [(OA, (1)CB)]
    >>> glc.infer_residue_connections(triplet=True)
    [(OA, (1)CB), ((2)CA, OA)]
    


Adding residues and atoms
-------------------------

To add one or more new residue(s), we use the `add_residues` method, which expects a number `buildamol.Residue` objects as unnamed arguments.
Similarly, to add one or more new atom(s), we use the `add_atoms` method, which expects a number of `buildamol.Atom` objects as unnamed arguments.
Both methods allow to specify the parent object (`chain` or `residue`) via an optional argument and will automatically choose the last
chain or residue if none is specified.

.. code-block:: python

    from buildamol import Residue, Atom

    new_residue = Residue("XYZ", 1, " ")
    
    # do things with the residue here
    # ...

    # add the residue to the molecule
    # (add it to the last chain, whichever that may be)
    glc.add_residues(new_residue)

    new_atom = Atom("X", [0, 0, 0])

    # add the atom to the first residue in the `glc` molecule
    glc.add_atoms(new_atom, residue=1)

Removing residues and atoms
---------------------------

In order to remove residues or atoms or bonds, we can use the `remove_residues`, `remove_atoms` and `remove_bond`(yes, singluar!) methods.
They work exactly like their `add_` counterparts, but instead of adding, they remove the specified objects.

.. code-block:: python

    # remove the first residue
    glc.remove_residues(1)

    # remove the first atom
    glc.remove_atoms(1)

    # remove the bond between the first and second atom
    glc.remove_bond(1, 2)


Adjusting labelling
-------------------

Single-residue molecules that were loaded from a PDB file may not use the same atom labelling as the `PDBE` and `CHARMM` databases.
In order to quickly adjust the labelling, a method `autolabel` exists. `autolabel` uses the atom connectivity and a rule-based algorithm to infer 
the most likely atom labels. However, since this method is not perfect, it is recommended to check the labelling afterwards and adjust it manually if necessary.

If working with large molecules that follow another labelling scheme it may be more efficient to simply
define your own linkage recipies or patches (see the documentation of `linkages`) that use the the appropriate labelling scheme.

.. code-block:: python

    # load a molecule from a PDB file
    glc = Molecule.from_pdb("glucose.pdb")

    # adjust the labelling
    glc.autolabel()

    # save the molecule to a new PDB file
    glc.to_pdb("glucose_autolabelled.pdb")


Another common operation is the adjustment of chain names and residue seqids. This can be done using the `reindex` method.
This method accepts three starting values for the chain name, residue seqid and atom serial number and will then reindex all chains, residues and atoms
to ensure they are continuously numbered and labelled. Some internal methods used when connecting different molecules are reliant on a continuous
numbering scheme, so this method should be called before connecting molecules that were loaded from PDB files.

.. code-block:: python

    # load a molecule from a PDB file
    glc = Molecule.from_pdb("glucose.pdb")

    # reindex the molecule
    glc.reindex()

We can also use one molecule as a "reference" for reindexing another molecule to make sure there are now labelling conflicts between them in case we want to connect them together later (this is usually done internally by BuildAMol automatically).

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
    # which is pre-defined in the default CHARMMTopology
    glc.repeat(10, "14bb")

    # Now we have a cellulose of 10 glucoses

In the above example we used the `repeat` method explicitly, but we could also achieve the same with the short-hand `*=`. For this to work, we need to specify 
the linkage type beforehand. We do this by setting the `patch` attribute before using any operator.

.. code-block:: python

    # specify the "default" linkage type for connecting 
    # other molecules to this glucose
    glc.linkage = "14bb"

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
    glc.linkage = "14bb"

    # now attach the galactose to the glucose
    lac = glc + gal


Of course, if there is a `+` operator there should also be a `+=` operator, which is simply the equivalent of `attach` with `inplace=True`.

.. code-block:: python

    glc.linkage = "14bb"
    glc += gal


    # Now 'glc' is a lactose molecule


Setting default Modifiers
-------------------------

So far, we have always worked with a 1-4 beta-beta glycosidic linkage, which we apparently could select using the string `"14bb"`. 
But what if we want to use a different linkage type? For example, a 1-4 alpha-beta glycosidic linkage?  You of course noticed, 
that `attach` and `repeat` accept an argument `link` which allows you to specify the linkage type, and that if you leave it blank
the default linkage type is used. But how do we set the default linkage type?

Let's first check what linkage types are available by default anyway. Have you noticed an argument named `_topology`
at the end of the `attach` or `repeat` methods? The `topology` refers to the underlying _CHARMM_ topology which hosts the linkage type information.
By default a topology is already loaded in BuildAMol's framework so it is not necessary for the user to specify anything here, but we can check which linkage types are available by:

.. code-block:: python

    import buildamol as bam

    # get the default topology
    topology = bam.get_default_topology()

    print(topology.patches)


Any of these linkages can be referenced by their name, e.g. `"14bb"` or `"14ab"`. 

Wait a second, my desired linkage is not in the list! What now?! Well, you can always define a new `Linkage` to suit your needs.
Check out the documentation on :ref:`linkages` for more information on how to do this. If you have your desired linkage 
ready to go, set it as the default by:

.. code-block:: python

    my_molecule.linkage = my_linkage
    
    # or if you feel "super correct"
    my_molecule.set_linkage(my_linkage)

    # or if you feel "extra cocky"
    my_molecule % my_linkage # <- the modulo operator assignes the "modifier" to the molecule


Now any call to `attach`, `repeat`, or any of its operator proxies will use your defined linkage by default.

Setting the default Residue for attachment
------------------------------------------

When defining a `Linkage` we specify which atoms are supposed to be connected and removed, but we do not specify which residues
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
"""

from copy import deepcopy
import os
from typing import Union

import numpy as np
import Bio.PDB as bio

import buildamol.core.entity as entity
import buildamol.core.Linkage as Linkage
import buildamol.utils as utils
import buildamol.structural as structural
import buildamol.resources as resources

__all__ = [
    "Molecule",
    "read_pdb",
    "write_pdb",
    "read_molfile",
    "write_molfile",
    "read_cif",
    "write_cif",
    "read_smiles",
    "molecule",
    "connect",
    "react",
    "polymerize",
    "make_smiles",
    "query_pubchem",
    "phosphorylate",
    "methylate",
    "acetylate",
    "amidate",
    "hydroxylate",
    "carboxylate",
    "phenolate",
    "benzylate",
    "thiolate",
]


def read_pdb(
    filename: str, id: str = None, multimodel: bool = False, has_atom_ids: bool = True
) -> "Molecule":
    """
    Read a PDB file and return a molecule.

    Parameters
    ----------
    filename : str
        The path to the PDB file
    id : str
        The id of the molecule
    multimodel : bool
        Whether to read all models and return a list of molecules.
    has_atom_ids : bool
        Whether the PDB file contains atom ids. If the file does not, the atom ids can be auto-generated if this is set to false.

    Returns
    -------
    molecule : Molecule or list
        The molecule or a list of molecules if multimodel is True
    """
    if multimodel:
        models = utils.pdb.find_models(filename)
        molecules = []
        for model in models:
            new = Molecule.from_pdb(
                filename, id=id, model=model, has_atom_ids=has_atom_ids
            )
            molecules.append(new)
        return molecules
    return Molecule.from_pdb(filename, id=id)


def write_pdb(mol: "Molecule", filename: str) -> None:
    """
    Write a molecule to a PDB file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The path to the PDB file
    """
    mol.to_pdb(filename)


def read_cif(filename: str, id: str = None) -> "Molecule":
    """
    Read a CIF file and return a molecule.

    Parameters
    ----------
    filename : str
        The path to the CIF file
    id : str
        The id of the molecule

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_cif(filename, id=id)


def write_cif(mol: "Molecule", filename: str) -> None:
    """
    Write a molecule to a CIF file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The path to the CIF file
    """
    mol.to_cif(filename)


def read_molfile(filename: str, id: str = None) -> "Molecule":
    """
    Read a MOL file and return a molecule.

    Parameters
    ----------
    filename : str
        The path to the MOL file
    id : str
        The id of the molecule

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_molfile(filename, id=id)


def write_molfile(mol: "Molecule", filename: str) -> None:
    """
    Write a molecule to a MOL file.

    Parameters
    ----------
    mol : Molecule
        The molecule to write
    filename : str
        The path to the MOL file
    """
    mol.to_molfile(filename)


def read_smiles(smiles: str, id: str = None) -> "Molecule":
    """
    Read a SMILES string and return a molecule.

    Parameters
    ----------
    smiles : str
        The SMILES string

    Returns
    -------
    molecule : Molecule
        The molecule
    """
    return Molecule.from_smiles(smiles, id=id)


def make_smiles(
    mol: "Molecule", isomeric: bool = True, write_hydrogens: bool = False
) -> str:
    """
    Generate a SMILES string from a molecule.

    Parameters
    ----------
    mol : Molecule
        The molecule
    isomeric : bool
        Whether to include stereochemistry information
    write_hydrogens : bool
        Whether to include hydrogens in the SMILES string

    Returns
    -------
    smiles : str
        The SMILES string
    """
    return mol.to_smiles(isomeric, write_hydrogens)


def query_pubchem(query: str, by: str = "name") -> "Molecule":
    """
    Query the PubChem database for a given
    query string to obtain a Molecule object.

    Parameters
    ----------
    query : str
        The query string
    by : str
        The type of query to perform. Can be one of:
        he method to search by. This can be any of the following:

        - cid
        - name
        - smiles
        - sdf
        - inchi
        - inchikey
        - formula

    Returns
    -------
    Molecule or None
        The molecule or None if no match was found
    """
    try:
        return Molecule.from_pubchem(query, by=by)
    except:
        return None


def molecule(mol=None) -> "Molecule":
    """
    Generate a molecule from an input. If the input is a string, the string can be a PDB id, some filename, SMILES or InChI string, IUPAC name or abbreviation.
    This function will try its best to automatically generate the molecule with minimal user effort. However, using a dedicated classmethod is
    recommended for more efficient and predictable results.

    Parameters
    ----------
    mol : str or structure-like object
        An input string or structure-like object such as a BioPython Structure or RDKit Molecule, etc.
        If nothing is provided, a new empty molecule is generated.

    Returns
    -------
    molecule : Molecule
        The generated molecule

    Examples
    --------
    >>> from buildamol import molecule
    >>> mol = molecule("GLC")
    >>> mol = molecule("GLC.pdb")
    >>> mol = molecule("alpha-d-glucose")
    """
    if mol is None:
        return Molecule.empty()

    if isinstance(mol, bio.Structure.Structure) or isinstance(
        mol, entity.base_classes.Structure
    ):
        return Molecule(mol)
    elif isinstance(mol, Molecule):
        return mol.copy()
    elif "openbabel" in str(type(mol).__mro__[0]):
        return Molecule.from_pybel(mol)
    elif "rdkit" in str(type(mol).__mro__[0]):
        return Molecule.from_rdkit(mol)
    elif (
        "openmm" in str(type(mol).__mro__[0])
        and hasattr(mol, "topology")
        and hasattr(mol, "positions")
    ):
        return Molecule.from_openmm(mol.topology, mol.positions)
    elif "stk" in str(type(mol).__mro__[0]):
        return Molecule.from_stk(mol)

    if not isinstance(mol, str):
        raise ValueError("input must be a structure object or a string")

    # ------------------
    # mol may be a PDB id
    # ------------------
    if resources.has_compound(mol):
        return resources.get_compound(mol)

    if os.path.isfile(mol):
        _mol = mol.lower().strip()
        if _mol.endswith(".pdb"):
            return Molecule.from_pdb(mol)
        elif _mol.endswith(".cif"):
            return Molecule.from_cif(mol)
        elif _mol.endswith(".pkl"):
            return Molecule.load(mol)
        elif _mol.endswith(".json"):
            return Molecule.from_json(mol)
        elif _mol.endswith(".xml"):
            return Molecule.from_xml(mol)
        elif (
            _mol.endswith(".mol")
            or _mol.endswith(".mol2")
            or _mol.endswith(".sdf")
            or _mol.endswith(".sd")
        ):
            return Molecule.from_molfile(mol)

    if " " not in mol:
        try:
            return Molecule.from_smiles(mol)
        except:
            pass

    try:
        return Molecule.from_pubchem(mol)
    except:
        pass

    raise ValueError(f"Could not generate molecule from input: {mol}")


def polymerize(
    mol: "Molecule",
    n: int,
    link: Union[str, "Linkage.Linkage"] = None,
    inplace: bool = False,
) -> "Molecule":
    """
    Polymerize a molecule

    Parameters
    ----------
    mol : Molecule
        The molecule to polymerize
    n : int
        The number of monomers to add
    link : str or Linkage
        The linkage to use for polymerization. If None, the default linkage of the molecule is used.
    inplace : bool
        Whether to polymerize the molecule in place or return a new molecule

    Returns
    -------
    Molecule
        The polymerized molecule
    """
    if link is None and mol._linkage is None:
        raise ValueError(
            "No patch or recipe provided and no default is set on the molecule"
        )
    return mol.repeat(n, link, inplace=inplace)


def connect(
    mol_a: "Molecule",
    mol_b: "Molecule",
    link: Union[str, "Linkage.Linkage"],
    at_residue_a: Union[int, "bio.Residue.Residue"] = None,
    at_residue_b: Union[int, "bio.Residue.Residue"] = None,
    copy_a: bool = True,
    copy_b: bool = True,
    _topology=None,
    use_patch: bool = True,
) -> "Molecule":
    """
    Connect two molecules together

    Parameters
    ----------
    mol_a : Molecule
        The first (target) molecule
    mol_b : Molecule
        The second (source) molecule
    link : Linkage or str
        The linkage to use for connection. This can be either an instance of the Linkage class or a string identifier
        of a pre-defined patch in the (currently loaded default or specified) CHARMMTopology.
    at_residue_a : int or bio.PDB.Residue
        The residue of the first molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    at_residue_b : int or bio.PDB.Residue
        The residue of the second molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    copy_a : bool
        Whether to copy the first molecule before connecting
    copy_b : bool
        Whether to copy the second molecule before connecting.
        If False, all atoms of the second molecule will be added to the first molecule.
    _topology : CHARMMTopology
        A specific topology to use in case a pre-existing patch is used as link and only the string identifier
        is supplied.
    use_patch : bool
        If the linkage has internal coordinates available (i.e. is a "patch") these are used by default. Set this to False
        to force-use stitching and its associated conformational optimization instead.

    Returns
    -------
    Molecule
        The connected molecule
    """
    if copy_b:
        mol_b = mol_b.copy()
    new = mol_a.attach(
        mol_b,
        link,
        at_residue=at_residue_a,
        other_residue=at_residue_b,
        inplace=not copy_a,
        _topology=_topology,
        use_patch=use_patch,
    )
    return new


def react(
    mol_a: "Molecule",
    mol_b: "Molecule",
    egroup: "FunctionalGroup",
    ngroup: "FunctionalGroup",
    a_is_electrophile: bool = True,
    at_residue_a: Union[int, "bio.Residue.Residue"] = None,
    at_residue_b: Union[int, "bio.Residue.Residue"] = None,
    copy_a: bool = True,
    copy_b: bool = True,
) -> "Molecule":
    """
    Connect two molecules together by imitating a chemical reaction based on functional groups.

    Parameters
    ----------
    mol_a : Molecule
        The first (target) molecule
    mol_b : Molecule
        The second (source) molecule
    egroup : FunctionalGroup
        The functional group of the first molecule to connect to.
    ngroup : FunctionalGroup
        The functional group of the second molecule to connect to.
    a_is_electrophile : bool
        Whether the first molecule is the electrophile (True) or nucleophile (False)
    at_residue_a : int or bio.PDB.Residue
        The residue of the first molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    at_residue_b : int or bio.PDB.Residue
        The residue of the second molecule to connect to. If an integer is provided, the seqid must be used, starting at 1.
    copy_a : bool
        Whether to copy the first molecule before connecting
    copy_b : bool
        Whether to copy the second molecule before connecting.
        If False, all atoms of the second molecule will be added to the first molecule.

    Returns
    -------
    Molecule
        The connected molecule
    """
    new = mol_a.react_with(
        mol_b,
        egroup,
        ngroup,
        a_is_electrophile,
        at_residue_a,
        at_residue_b,
        not copy_a,
        not copy_b,
    )
    return new


def _modify(
    mol: "Molecule",
    modifier: "Molecule",
    at_atom,
    delete,
    modifier_at_atom: "str",
    modifier_deletes: "List[str]",
    inplace: bool = True,
):
    """
    The core function of all modification functions. This function is used to modify a molecule by attaching a modifier (e.g. a functional group) molecule to it.

    Parameters
    ----------
    mol : Molecule
        The molecule to modify
    modifier : Molecule
        The modifier molecule
    at_atom
        The atom of the molecule to modify. This can be any identifier that will allow to get an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete
        The atom to delete in the molecule. This can be any identifier that will allow to get an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well. If a list is provided as at_atom, use either None to specify no deleters
        or a list of equal length to specify the deleters for each atom in at_atom. None is allowed for individual entries in the list but the list
        needs the same length as at_atom.
    modifier_at_atom : str
        The atom of the modifier molecule to attach to the molecule
    modifier_deletes : List[str]
        The atoms to delete in the modifier molecule. The first entry needs to be a direct neighbor of the modifier_at_atom.
    inplace : bool
        Whether to modify the molecule in place or return a new molecule
    """
    if isinstance(at_atom, (list, tuple)):
        if delete is None:
            delete = (None for i in at_atom)
        elif isinstance(delete, (list, tuple)) and len(delete) != len(at_atom):
            raise ValueError(
                f"delete must be a list of same length if at_atom is a list. at_atom has length {len(at_atom)} but delete has length {len(delete)}"
            )
        elif not isinstance(delete, (list, tuple)):
            raise ValueError(
                f"delete must be a list if at_atom is a list, got {type(delete)} instead"
            )

        if not inplace:
            mol = mol.copy()
        for a, d in zip(at_atom, delete):
            _modify(
                mol,
                modifier,
                a,
                d,
                modifier_at_atom,
                modifier_deletes,
                inplace=True,
            )
        return mol

    at_atom = mol.get_atom(at_atom)
    at_residue = at_atom.get_parent()
    if delete:
        delete = mol.get_atom(delete)
    else:
        delete = mol.get_hydrogen(at_atom)

    modifier_at_atom = modifier.get_atom(modifier_at_atom)

    dist = structural.single_bond_lengths[modifier_at_atom.element].get(
        at_atom.element, None
    )
    if dist:
        mol.adjust_bond_length(at_atom, delete, dist)
        modifier.adjust_bond_length(modifier_at_atom, modifier_deletes[0], dist)

    l = Linkage.linkage(
        at_atom.id, modifier_at_atom.id, [delete.id], modifier_deletes, id="MODIFY"
    )
    mol = mol.attach(modifier, l, at_residue=at_residue, inplace=inplace)
    return mol


def phosphorylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Phosphorylate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to phosphorylate
    at_atom : int or str or Atom
        The atom to phosphorylate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to phosphorylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to phosphorylate the molecule in place or return a new molecule

    Returns
    -------
    Molecule
        The phosphorylated molecule
    """
    resources.load_small_molecules()
    phos = Molecule.from_compound("PO4")
    return _modify(mol, phos, at_atom, delete, "P", ["O2"], inplace)


def methylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Methylate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to methylate
    at_atom : int or str or Atom
        The atom to methylate.This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
       The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to methylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to methylate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    methyl = Molecule.from_compound("CH3")
    return _modify(mol, methyl, at_atom, delete, "C", ["HC1"], inplace)


def acetylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Acetylate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to acetylate
    at_atom : int or str or Atom
        The atom to acetylate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to acetylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to acetylate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    acetyl = Molecule.from_compound("ACE")
    return _modify(mol, acetyl, at_atom, delete, "C", ["H"], inplace)


def hydroxylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Hydroxylate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to hydroxylate
    at_atom : int or str or Atom
        The atom to hydroxylate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to hydroxylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to hydroxylate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    hydroxyl = Molecule.from_compound("HOH")
    return _modify(mol, hydroxyl, at_atom, delete, "O", ["H1"], inplace)


def amidate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Amidate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to amidate
    at_atom : int or str or Atom
        The atom to amidate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to amidate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to amidate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    amide = Molecule.from_compound("NH3")
    return _modify(mol, amide, at_atom, delete, "N", ["HN1"], inplace)


def carboxylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Carboxylate a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to carboxylate
    at_atom : int or str or Atom
        The atom to carboxylate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to carboxylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to carboxylate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    carboxyl = Molecule.from_compound("CBX")
    return _modify(mol, carboxyl, at_atom, delete, "C", ["H"], inplace)


def benzylate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Add a benzyl group to a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to benzylate
    at_atom : int or str or Atom
        The atom to benzylate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to benzylate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to benzylate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    benzyl = Molecule.from_compound("BNZ")
    return _modify(mol, benzyl, at_atom, delete, "C1", ["H1"], inplace)


def phenolate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    how: str = "para",
    inplace: bool = True,
) -> "Molecule":
    """
    Add a phenol group to a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to phenolate
    at_atom : int or str or Atom
        The atom to phenolate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to phenolate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    how: str
        The position of the hydroxyl group on the phenol. Can be one of "ortho", "meta", or "para".
    inplace : bool
        Whether to phenolate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    phenol = Molecule.from_compound("MBN")
    if how == "ortho":
        a, d = "C2", "H2"
    elif how == "meta":
        a, d = "C3", "H3"
    elif how == "para":
        a, d = "C4", "H4"
    else:
        raise ValueError(
            "Invalid value for 'how'. Must be one of 'ortho', 'meta', or 'para'"
        )
    return _modify(mol, phenol, at_atom, delete, a, [d], inplace)


def thiolate(
    mol: "Molecule",
    at_atom: Union[int, str, entity.base_classes.Atom],
    delete: Union[int, str, entity.base_classes.Atom] = None,
    inplace: bool = True,
) -> "Molecule":
    """
    Add a thiol group to a molecule at one or more specific atoms

    Parameters
    ----------
    mol : Molecule
        The molecule to thiolate
    at_atom : int or str or Atom
        The atom to thiolate. This can be any input that will allow to obtain an Atom object from the molecule.
        Alternatively, a list of such inputs can be provided as well.
    delete : int or str or Atom
        The atom to delete. This can be any input that will allow to obtain an Atom object from the molecule.
        This atom needs to be in the same residue as the atom to thiolate. If not provided, any Hydrogen atom attached to the at_atom will be deleted.
        If at_atom is a list, delete can be a list of the same length or None.
    inplace : bool
        Whether to thiolate the molecule in place or return a new molecule
    """
    resources.load_small_molecules()
    thiol = Molecule.from_compound("HOH")
    thiol.get_atom("O").set_element("S")  # change the oxygen to sulfur
    return _modify(mol, thiol, at_atom, delete, "S", ["H1"], inplace)


class Molecule(entity.BaseEntity):
    """
    A molecule to add onto a scaffold.
    A molecule consists of a single chain.

    Parameters
    ----------
    structure : bio.PDB.Structure
        A biopython structure object
    root_atom : str or int or Atom
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
        root_atom: Union[str, int, entity.base_classes.Atom] = None,
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
        Create an Molecule without any atoms in it. This will have a single Model and Chain by default, however.

        Parameters
        ----------
        id : str
            An id of the Molecule.

        Returns
        -------
        Molecule
            An empty Molecule object
        """
        structure = structural.make_empty_structure(id)
        new = cls(structure)
        return new

    @classmethod
    def new(
        self,
        id: str = None,
        resname: str = "UNK",
        atoms: list = None,
        bonds: list = None,
    ) -> "Molecule":
        """
        Create a new Molecule with a single residue

        Parameters
        ----------
        id : str
            The id of the Molecule. By default an id is inferred from the filename.
        resname : str
            The resname of the residue to add
        atoms : list
            A list of Atom objects to add to the residue
        bonds : list
            A list of bonds to add to the residue

        Returns
        -------
        Molecule
            A new Molecule object
        """
        new = self.empty(id)
        new.add_residues(entity.base_classes.Residue(resname))
        if atoms:
            new.add_atoms(atoms)
        if bonds:
            new.add_bonds(bonds)
        return new

    # @classmethod
    # def from_pdb(
    #     cls,
    #     filename: str,
    #     root_atom: Union[str, int] = None,
    #     id: str = None,
    #     model: int = 0,
    #     chain: str = None,
    # ) -> "Molecule":
    #     """
    #     Read a Molecule from a PDB file

    #     Parameters
    #     ----------
    #     filename : str
    #         Path to the PDB file
    #     root_atom : str or int
    #         The id or the serial number of the root atom (optional)
    #     id : str
    #         The id of the Molecule. By default an id is inferred from the filename.
    #     model : int
    #         The model to use from the structure. Defaults to 0. This may be any
    #         valid identifier for a model in the structure, such as an integer or string.
    #     chain : str
    #         The chain to use from the structure. Defaults to the first chain in the structure.

    #     Returns
    #     -------
    #     Molecule
    #         The Molecule object
    #     """
    #     if id is None:
    #         id = utils.filename_to_id(filename)
    #     struct = utils.defaults.__bioPDBParser__.get_structure(id, filename)
    #     new = cls(struct, root_atom, model=model, chain=chain)
    #     bonds = utils.pdb.parse_connect_lines(filename)
    #     if len(bonds) != 0:
    #         for b in bonds:
    #             new.add_bond(*b)
    #     return new

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
        if mol is None:
            raise ValueError(f"Could not find compound '{compound}' using '{by}'")
        elif isinstance(mol, list):
            raise ValueError(
                f"Multiple compounds found using '{by}={compound}', choose any of these ids specifically {[i.id for i in mol]}"
            )
        if root_atom:
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
        obj = structural.read_smiles(smiles, add_hydrogens)
        if structural.smiles.use_rdkit:
            new = cls.from_rdkit(obj, id=id)
        elif structural.smiles.use_openbabel:
            obj.title = id
            new = cls.from_pybel(obj)
        if root_atom:
            new.set_root(root_atom)
        if id is not None:
            new.rename_residue(1, id[:3])
        return new

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
        integrated in buildamol. It is advisable to run `autolabel` or `relabel_hydrogens` on the molecule.
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
        _new = cls(new.structure)
        _new.add_bonds(*(i.to_tuple() for i in new._bonds))
        new = _new
        new.id = _compound_2d.iupac_name
        if root_atom:
            new.set_root(root_atom)
        return new

    @classmethod
    def from_geometry(
        cls,
        geometry: "structural.geometry.Geometry",
        atoms: list,
        id: str = None,
        resname: str = "UNK",
        direction: str = None,
    ):
        """
        Create a new Molecule using a molecular geometry and a list of starting atoms.
        This will place the atoms in the right spacial coordinates and fill up the provided starting atoms with hydrogens to match the geometry's defined number of atoms.

        Parameters
        ----------
        geometry : Geometry
            The geometry object
        atoms : list
            A list of Atom objects
        id : str
            The id of the Molecule.
        resname : str
            The resname of the residue to add
        direction : str
            The direction of the atoms (in case of a geometry that has planar and axial directions).
            This can be either "planar" or "axial".

        Returns
        -------
        Molecule
            The Molecule object
        """
        if len(atoms) > geometry.max_points:
            n = sum(1 for idx, atom in enumerate(atoms) if atom.coord.sum() != 0)
            n = n or 1
            _atoms = atoms[:n]
            geometry.make_and_apply(_atoms, atoms, direction=direction)

        if len(atoms) == geometry.size:
            bonds = [entity.base_classes.Bond(atoms[0], a) for a in atoms[1:]]
        else:
            atoms, bonds = geometry.fill_hydrogens(*atoms, direction=direction)
        new = cls.new(id=id, atoms=atoms, bonds=bonds, resname=resname)
        return new

    def to_smiles(self, isomeric: bool = True, write_hydrogens: bool = False) -> str:
        """
        Convert the molecule to a SMILES string

        Parameters
        ----------
        isomeric : bool
            Whether to include stereochemistry information in the SMILES string
        write_hydrogens : bool
            Whether to include hydrogens in the SMILES string

        Returns
        -------
        str
            The SMILES string
        """
        return structural.make_smiles(
            self, isomeric=isomeric, add_hydrogens=write_hydrogens
        )

    def get_residue_connections(
        self,
        residue_a=None,
        residue_b=None,
        triplet: bool = True,
        direct_by: str = None,
    ):
        """
        Get bonds between atoms that connect different residues in the structure
        This method is different from `infer_residue_connections` in that it works
        with the already present bonds in the molecule instead of computing new ones.

        Parameters
        ----------
        residue_a, residue_b : Union[int, str, tuple, bio.Residue.Residue]
            The residues to consider. If None, all residues are considered.
            Otherwise, only between the specified residues are considered.
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
            In the case of "root", the bonds are sorted so that the first atom is graph-closer to the root than the second atom.
            Set to None to not sort the bonds.

        Returns
        -------
        set
            A set of tuples of atom pairs that are bonded and connect different residues
        """
        bonds = super().get_residue_connections(
            residue_a=residue_a, residue_b=residue_b, triplet=triplet
        )
        if direct_by is not None:
            direct_connections = None
            if direct_by == "resid":
                direct_connections = super().get_residue_connections(
                    residue_a=residue_a, residue_b=residue_b, triplet=False
                )
                direct_connections1, direct_connections2 = set(
                    a for a, b in direct_connections
                ), set(b for a, b in direct_connections)
                direct_connections = direct_connections1.union(direct_connections2)
            bonds = self._direct_bonds(bonds, direct_by, direct_connections)
        return bonds

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

        _other = self.copy()
        if not inplace:
            obj = self.copy()
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
        link: Union[str, "Linkage.Linkage"] = None,
        at_residue: Union[int, "entity.base_classes.Residue"] = None,
        other_residue: Union[int, "entity.base_classes.Residue"] = None,
        use_patch: bool = None,
        inplace: bool = True,
        other_inplace: bool = False,
        _topology=None,
    ):
        """
        Attach another structure to this one using a Patch or a Recipe.

        Parameters
        ----------
        other : Molecule
            The other molecule to attach to this one
        link : str or Linkage
            Either a Patch to apply when attaching or a Recipe to use when stitching.
            If None is defined, the default patch or recipe that was set earlier on the molecule is used.
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the defined `attach_residue` is used.
        other_residue : int or Residue
            The residue in the other molecule to attach this molecule to. If None, the defined `attach_residue` of the other molecule is used.
        use_patch : bool
            If the specified linkage is a patch (has internal coordinates) it can and is by default applied as a patch. However, it can also be used as a recipe.
            Set this to false if you want to use the patch as a recipe.
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.
        other_inplace : bool
            All atoms from the other molecule are integrated into this one. Hence, the other molecule is left empty. If False, a copy of the other molecule is used.
            Thus leaving the original molecule intact.
        _topology : Topology
            The topology to use when attaching. If None, the topology of the molecule is used. Only used if the patch is a string.
        """
        if not isinstance(other, Molecule):
            raise TypeError("Can only attach a Molecule to another Molecule")

        if use_patch is None:
            use_patch = utils.auxiliary.USE_IC

        if not inplace:
            obj = self.copy()
        else:
            obj = self

        if not link:
            link = obj._linkage
            if not link:
                raise ValueError("Cannot attach a molecule without a patch defined")

        if not other_inplace:
            _other = other.copy()
        else:
            _other = other

        if isinstance(link, str):
            if not _topology:
                _topology = resources.get_default_topology()
            link = _topology.get_patch(link)

        if link.has_IC and use_patch:
            obj.patch_attach(
                _other,
                link,
                at_residue=at_residue,
                other_residue=other_residue,
                _topology=_topology,
            )
        else:
            obj.stitch_attach(
                _other,
                link,
                at_residue=at_residue,
                other_residue=other_residue,
            )
        return obj

    def patch_attach(
        self,
        other: "Molecule",
        patch: Union["Linkage.Linkage", str] = None,
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

        p = structural.__default_keep_keep_patcher__
        p.apply(patch, self, other, at_residue, other_residue)
        p.merge()
        return self

    def stitch_attach(
        self,
        other: "Molecule",
        recipe: "Linkage.Linkage" = None,
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
        if not recipe and remove_atoms is None:
            if not self._linkage:
                raise AttributeError(
                    "No recipe was set for this molecule and no manual instructions were found. Either set a default recipe, provide a recipe when stitching, or provide the information about removed and bonded atoms directly."
                )
            recipe = self._linkage

        if recipe:
            target_atom, source_atom = recipe._stitch_ref_atoms
            return self.stitch_attach(
                other,
                remove_atoms=recipe.deletes[0],
                other_remove_atoms=recipe.deletes[1],
                at_atom=target_atom,
                other_at_atom=source_atom,
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

        if not at_residue:
            at_residue = self.attach_residue

        if not other_residue:
            other_residue = other.attach_residue

        # since we are already copying before we can use the keep-keep stitcher, actually...
        p = structural.__default_keep_keep_stitcher__
        p.apply(
            self,
            other,
            remove_atoms,
            other_remove_atoms,
            at_atom,
            other_at_atom,
            at_residue,
            other_residue,
        )
        self = p.merge()
        return self

    def react_with(
        self,
        other: "Molecule",
        egroup: "FunctionalGroup",
        ngroup: "FunctionalGroup",
        as_electrophile: bool = True,
        at_residue: Union[int, "entity.base_classes.Residue"] = None,
        other_residue: Union[int, "entity.base_classes.Residue"] = None,
        inplace: bool = True,
        other_inplace: bool = False,
    ) -> "Molecule":
        """
        React this molecule with another molecule using functional groups to
        automatically create a linkage.

        Parameters
        ----------
        other : Molecule
            The other molecule to react with
        egroup : FunctionalGroup
            The electrophilic functional group to use
        ngroup : FunctionalGroup
            The nucleophilic functional group to use
        as_electrophile : bool
            Whether to use this molecule as the electrophile or the nucleophile
        at_residue : int or Residue
            The residue to attach the other molecule to. If None, the last residue of the molecule.
        other_residue : int or Residue
            The residue of the other molecule to attach. If None, the first residue of the other molecule.
        inplace : bool
            If True the molecule is directly modified, otherwise a copy of the molecule is returned.
        other_inplace : bool
            All atoms from the other molecule are integrated into this one. Hence, the other molecule is left empty. If False, a copy of the other molecule is used.
            Thus leaving the original molecule intact.
        Returns
        -------
        molecule
            The modified molecule (either the original object or a copy)
        """
        if not inplace:
            obj = self.copy()
        else:
            obj = self

        if not other_inplace:
            _other = other.copy()
        else:
            _other = other

        _backup_attach_residue = obj.attach_residue
        if at_residue is not None:
            obj.set_attach_residue(at_residue)
        if other_residue is not None:
            _other.set_attach_residue(other_residue)

        if as_electrophile:
            link = Linkage.Linkage.from_functional_groups(obj, egroup, _other, ngroup)
        else:
            link = Linkage.Linkage.from_functional_groups(_other, egroup, obj, ngroup)

        obj.stitch_attach(_other, link)
        obj.set_attach_residue(_backup_attach_residue)
        return obj

    def optimize(
        self,
        residue_graph: bool = None,
        algorithm: str = None,
        rotatron: str = None,
        rotatron_kws: dict = None,
        algorithm_kws: dict = None,
        inplace: bool = True,
    ):
        """
        Optimize the molecule's conformation. This is a convenience method with less
        customizability than a manual optimization using the `optimizers` module.

        Parameters
        ----------
        residue_graph : bool
            Whether to use the residue graph or the full atom graph for optimization.
            The residue graph is faster but less accurate. If the molecule is larger than 100 atoms,
            the residue graph is used by default.
        algorithm : str
            The optimization algorithm to use. If not provided, an algorithm is determined based on the molecule's size.
            This can be one of the following:
            - "genetic" for a genetic algorithm
            - "scipy" for a scipy-implemented gradient-based optimization
            - "swarm" for a particle swarm optimization
            - "anneal" for a simulated annealing optimization
            - "rdkit" for an RDKit-implemented force-field-based optimization (if RDKit is installed)
        rotatron : str
            The rotatron to use. This can be one of the following:
            - "distance" for a distance-based rotatron (default)
            - "overlap" for an overlap-based rotatron
            - "forcefield" for a force-field-based rotatron
        algorithm_kws : dict
            Keyword arguments to pass to the optimization algorithm
        rotatron_kws : dict
            Keyword arguments to pass to the rotatron
        inplace : bool
            Whether to optimize the molecule in place or return a copy.

        Returns
        -------
        molecule
            The optimized molecule (either the original object or a copy)
        """
        import buildamol.optimizers as optimizers

        if residue_graph is None:
            residue_graph = self.count_atoms() > 300

        if not algorithm_kws:
            algorithm_kws = {}

        if not rotatron_kws:
            rotatron_kws = {}

        algorithm = algorithm or optimizers.auto_algorithm(self)

        if algorithm == "rdkit":
            new = self.copy() if not inplace else self
            out = optimizers.rdkit_optimize(new)
            out.id = self.id
            return out
        elif algorithm not in ("scipy", "genetic", "swarm", "anneal"):
            raise ValueError(f"Unknown algorithm {algorithm}")

        rotatron = rotatron or "distance"
        if rotatron == "distance":
            rotatron = optimizers.DistanceRotatron
        elif rotatron == "overlap":
            rotatron = optimizers.OverlapRotatron
        elif rotatron == "forcefield":
            rotatron = optimizers.ForceFieldRotatron
        else:
            raise ValueError(f"Unknown rotatron {rotatron}")

        if not inplace:
            new = self.copy()
        else:
            new = self

        if residue_graph:
            graph = new.make_residue_graph(detailed=True)
        else:
            graph = new.make_atom_graph(_copy=False)

        edges = graph.find_rotatable_edges(
            graph.central_node, min_ancestors=1, min_descendants=1
        )
        if len(edges) > 150:
            edges = graph.sample_edges(edges, edges // 5)

        env = rotatron(graph, edges, **rotatron_kws)

        return optimizers.optimize(new, env, algorithm, **algorithm_kws)

    def __add__(self, other) -> "Molecule":
        """
        Add two molecules together. This will return a new molecule.
        """
        if isinstance(other, (tuple, list)):
            new = self + other[0]
            for o in other[1:]:
                new += o
            return new

        elif isinstance(other, Molecule):
            patch = self._linkage
            if not patch:
                patch = other._linkage
            if not patch:
                raise RuntimeError(
                    "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 + mol2)"
                )

            # if linkages are functional groups instead, we use the react_with method
            if (
                self._linkage.__class__.__name__ == "FunctionalGroup"
                and other._linkage.__class__.__name__ == "FunctionalGroup"
            ):
                new = self.react_with(
                    other, self._linkage, other._linkage, inplace=False
                )
                return new

            # otherwise we use the attach method
            new = self.attach(other, patch, inplace=False, other_inplace=False)
            return new

        new = self.copy()
        if isinstance(other, entity.base_classes.Atom):
            new.add_atoms(other)

        elif isinstance(other, entity.base_classes.Residue):
            new.add_residues(other)

        elif isinstance(other, entity.base_classes.Chain):
            new._model.add(other)

        elif isinstance(other, entity.base_classes.Model):
            new._base_struct.add(other)

        return new

    def __iadd__(self, other) -> "Molecule":
        """
        Attach another molecule to this one
        """
        if isinstance(other, (tuple, list)):
            for o in other:
                self += o
            return self

        elif isinstance(other, Molecule):
            patch = self._linkage
            if not patch:
                patch = other._linkage
            if not patch:
                raise RuntimeError(
                    "Cannot add two molecules together without a patch, set a default patch on either of them (preferably on the one you are adding to, i.e. mol1 in mol1 += mol2)"
                )

            # if linkages are functional groups instead, we use the react_with method
            if (
                self._linkage.__class__.__name__ == "FunctionalGroup"
                and other._linkage.__class__.__name__ == "FunctionalGroup"
            ):
                self.react_with(other, self._linkage, other._linkage)
                return self

            # otherwise we use the attach method
            self.attach(other, patch)
            return self

        elif isinstance(other, entity.base_classes.Atom):
            self.add_atoms(other)

        elif isinstance(other, entity.base_classes.Residue):
            self.add_residues(other)

        elif isinstance(other, entity.base_classes.Chain):
            self._model.add(other)

        elif isinstance(other, entity.base_classes.Model):
            self._base_struct.add(other)

        return self

    def __sub__(self, other) -> "Molecule":
        """
        Remove atoms, residues, chains or models from the molecule
        """
        if isinstance(other, (tuple, list)):
            new = self.copy()
            for o in other:
                new -= o
            return new

        new = self.copy()
        if isinstance(other, entity.base_classes.Atom):
            new.remove_atoms(other)

        elif isinstance(other, entity.base_classes.Residue):
            new.remove_residues(other)

        elif isinstance(other, entity.base_classes.Chain):
            new._model.detach_child(other.get_id())

        elif isinstance(other, entity.base_classes.Model):
            new._base_struct.detach_child(other.get_id())

        elif isinstance(other, Molecule):
            new.remove_residues(*other.residues)

        return new

    def __isub__(self, other) -> "Molecule":
        """
        Remove atoms, residues, chains or models from the molecule
        """
        if isinstance(other, (tuple, list)):
            for o in other:
                self -= o
            return self

        elif isinstance(other, entity.base_classes.Atom):
            self.remove_atoms(other)

        elif isinstance(other, entity.base_classes.Residue):
            self.remove_residues(other)

        elif isinstance(other, entity.base_classes.Chain):
            self._model.detach_child(other.get_id())

        elif isinstance(other, entity.base_classes.Model):
            self._base_struct.detach_child(other.get_id())

        elif isinstance(other, Molecule):
            self.remove_residues(*other.residues)

        return self

    def __mul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the * operator (i.e. mol * 3)
        This requires that the molecule has a patch defined
        """
        if not self._linkage:
            raise RuntimeError("Cannot multiply a molecule without a patch defined")

        if not isinstance(n, int):
            raise TypeError("Can only multiply a molecule by an integer")
        if n <= 0:
            raise ValueError("Can only multiply a molecule by a positive integer")

        new = self.attach(self, self._linkage, inplace=False, other_inplace=False)
        for i in range(n - 2):
            new = new.attach(self, self._linkage, inplace=True, other_inplace=False)

        return new

    def __imul__(self, n) -> "Molecule":
        """
        Add multiple identical molecules together using the *= operator (i.e. mol *= 3)
        This requires that the molecule has a patch defined
        """
        self.repeat(n)
        return self

    def __repr__(self):
        return f"Molecule({self.id})"

    def __str__(self):
        string = f"Molecule {self.id}:\n"
        if len(self.residues) > 1:
            string += str(self.make_residue_graph())
        else:
            string += str(self.residues[0])
        return string


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
    mol = Molecule.empty(id)
    residue = entity.base_classes.Residue("UNK", " ", 1)
    mol.add_residues(residue)
    # structure = bio.Structure.Structure(id)
    # model = bio.Model.Model(0)
    # structure.add(model)
    # chain = bio.Chain.Chain("A")
    # model.add(chain)
    # residue = bio.Residue.Residue((" ", 1, " "), "UNK", 1)
    # chain.add(residue)

    element_counts = {}
    adx = 1
    for atom in comp.atoms:
        element = atom.element
        element_counts[element] = element_counts.get(element, 0) + 1

        id = f"{element}{element_counts[element]}"
        _atom = entity.base_classes.Atom(
            id,
            coord=np.array((atom.x, atom.y, atom.z)),
            serial_number=adx,
            bfactor=0.0,
            occupancy=0.0,
            element=element.upper(),
            fullname=id,
            altloc=" ",
            pqr_charge=0.0,
        )
        residue.add(_atom)
        adx += 1

    for bond in bonds:
        mol.add_bond(bond.aid1, bond.aid2, bond.order)

    return mol

    # def _molecule_from_ctfile(id, ctfile):
    #     """
    #     Convert a CTFile object to a Molecule.

    #     Parameters
    #     ----------
    #     id : str
    #         The id of the molecule
    #     ctfile : CTFile
    #         The CTFile object

    #     Returns
    #     -------
    #     Molecule
    #         The molecule
    #     """
    #     new = Molecule.empty(id)
    #     residue = entity.base_classes.Residue("UNK", " ", 1)
    #     new.add_residues(residue)

    # element_counts = {}
    # for atom in ctfile.atoms:
    #     element_counts.setdefault(atom.atom_symbol, 1)
    #     count = element_counts[atom.atom_symbol]
    #     element_counts[atom.atom_symbol] += 1

    #     x = float(atom._ctab_data["x"])
    #     y = float(atom._ctab_data["y"])
    #     z = float(atom._ctab_data["z"])

    #     _atom = entity.base_classes.Atom(
    #         id=f"{atom.atom_symbol}{count}",
    #         coord=np.array((x, y, z)),
    #         serial_number=int(atom.atom_number),
    #         pqr_charge=float(atom.charge),
    #     )
    #     residue.add(_atom)

    # for bond in ctfile.bonds:
    #     a, b = bond.first_atom, bond.second_atom
    #     a = int(a.atom_number)
    #     b = int(b.atom_number)
    #     order = int(bond._ctab_data["bond_type"])
    #     for i in range(order):
    #         new.add_bond(a, b)

    # return new


if __name__ == "__main__":
    # ser = Molecule.from_json("ser.json")
    # link = Linkage.from_json("peptide_linkage.json")

    # pep = ser.attach(ser, link, use_patch=False)
    # pep.show()

    # # ser = Molecule.from_pubchem("SER")
    # ser.to_cif("ser.cif")
    # ser.to_pdb("ser.pdb")
    # recipe = Linkage()
    # recipe.add_delete("O1", "target")
    # recipe.add_delete("HO1", "target")
    # recipe.add_delete("HO4", "source")
    # recipe.add_bond(("C1", "O4"))

    # glc = Molecule.from_compound("GLC")
    # glc.set_linkage(recipe)

    # glc2 = glc.copy()

    # _current_residues = len(glc.residues)
    # glc3 = glc + glc2
    # glc3.show()
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

    from buildamol.utils import visual

    # v = visual.MoleculeViewer3D(glc)
    # v.show()

    # man = Molecule.from_pdb("support/examples/membrane.pdb", model=4)
    # man.infer_bonds()
    # man.

    man = Molecule.from_pubchem("benzene")
    man.infer_bonds(restrict_residues=False)

    g = man.make_residue_graph(True)
    v = visual.MoleculeViewer3D(g)

    for c in man.get_residue_connections():
        v.draw_vector(
            None, c[0].coord, 1.3 * (c[1].coord - c[0].coord), color="magenta"
        )
    v.show()
