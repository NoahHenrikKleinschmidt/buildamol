import os

BASE = os.path.dirname(os.path.dirname(__file__))
"""
The base directory of the package
"""

RESOURCES = os.path.join(BASE, "resources")
"""
The resources directory of the package
"""

DEFAULT_BONDS_FILE = os.path.join(RESOURCES, "default.bonds.pkl")
"""
The path to the default bonds file
"""
