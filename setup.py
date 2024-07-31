import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="buildamol",
    version="1.2.6",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@unibe.ch",
    description="A fragment-based molecular assembly toolkit for python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={
        "buildamol.resources": ["*.pkl", "*.json", "*.cif", "*.xml"],
    },
    url="https://github.com/NoahHenrikKleinschmidt/buildamol",
    packages=[
        "buildamol",
        "buildamol.core",
        "buildamol.resources",
        "buildamol.utils",
        "buildamol.structural",
        "buildamol.graphs",
        "buildamol.optimizers",
        "buildamol.extensions",
        "buildamol.extensions.complexes",
        "buildamol.extensions.polymers",
        "buildamol.extensions.bio",
        "buildamol.extensions.bio.proteins",
        "buildamol.extensions.bio.lipids",
        "buildamol.extensions.bio.glycans",
        "buildamol.extensions.molecular_dynamics",
        "buildamol.extensions.molecular_dynamics.atom_typing",
    ],
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "seaborn",
        "networkx",
        "biopython",
        "pdbecif",
        "periodictable",
        "plotly",
        "gym",
        "pubchempy",
        "tabulate",
        "scikit-learn",
        "ipywidgets",  # "ipywidgets<8.0.1", # I don't remember why the 8.0.1 was there
        "attrs",
    ],
    optional_requires={
        "visual": ["py3Dmol", "nglview"],
        "openbabel": ["openbabel"],
        "rdkit": ["rdkit"],
        "openmm": ["openmm"],
        "full": ["rdkit", "openbabel", "openmm", "py3Dmol", "nglview"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
)
