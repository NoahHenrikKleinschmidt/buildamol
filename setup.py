import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="biobuild",
    version="3.10.10",
    author="Noah H. Kleinschmidt",
    author_email="noah.kleinschmidt@students.unibe.ch",
    description="A fragment-based molecular assembly toolkit for python.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    package_data={
        "biobuild.resources": ["*.pkl", "*.json", "*.cif"],
    },
    url="https://github.com/NoahHenrikKleinschmidt/biobuild",
    packages=[
        "biobuild",
        "biobuild.core",
        "biobuild.resources",
        "biobuild.utils",
        "biobuild.structural",
        "biobuild.graphs",
        "biobuild.optimizers",
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
        "ipywidgets<8.0.1",
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
    python_requires=">=3.10",
)
