import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="biobuild",
    version="3.3.51-a",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A python package for building organic molecules",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={
        "biobuild.resources": ["*.pkl", "*.json", "*.cif"],
    },
    url="TBA",
    packages=[
        "biobuild",
        "biobuild.core",
        "biobuild.resources",
        "biobuild.utils",
        "biobuild.structural",
        "biobuild.graphs",
        "biobuild.optimizers",
        "biobuild.optimizers.environments",
        "biobuild.optimizers.agents",
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
        "alive_progress",
        "pyswarms",
        "gym",
        "pubchempy",
        "tabulate",
        "scikit-learn",
        "ipywidgets<8.0.1...",
    ],
    optional_requires={
        "openbabel": ["openbabel"],
        "nglview": ["nglview"],
        "rdkit": ["rdkit"],
        "openmm": ["openmm"],
        "md": ["openmm", "openbabel"],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
)
