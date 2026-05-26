"""
NCIP - NCI-Profiler
====================

A Python tool for detecting and analyzing non-covalent interactions
between proteins and small molecule ligands from PDB files.

Supported interaction types:
- Hydrogen Bonds
- Salt Bridges
- Hydrophobic Contacts
- π-π Stacking
- Cation-π Interactions
- Halogen Bonds
- Water Bridges
- XWH Bridges (Halogen-Water-Hydrogen)
- Metal Coordination

Quick Start:
    >>> from ncip import analyze, print_results
    >>> results, *_ = analyze("protein.pdb")
    >>> print_results(results)

Or from command line:
    $ ncip protein.pdb
    $ ncip protein.pdb --json results.json
    $ ncip --fetch 1ABC -l LIG --pml output.pml
"""

__version__ = "1.0.0"
__author__ = "Rahul"

# Main analysis functions
from .core import (
    analyze,
    analyze_single_ligand,
    print_results,
    save_tabular_data
)

# Data structures
from .parser import (
    Atom,
    Ring,
    Interaction,
    parse_pdb,
    distance,
    angle,
    centroid
)

# Individual interaction detection functions
from .interactions import (
    find_hydrogen_bonds,
    find_salt_bridges,
    find_hydrophobic_contacts,
    find_pi_stacking,
    find_cation_pi,
    find_halogen_bonds,
    find_water_bridges,
    detect_metal_coordination
)

# Visualization
from .visualizer import (
    generate_pymol_script,
    save_script_only
)

# Molecule detection
from .parser import detect_molecules

__all__ = [
    "__version__",
    "__author__",
    "analyze",
    "analyze_single_ligand",
    "print_results",
    "save_tabular_data",
    "Atom",
    "Ring",
    "Interaction",
    "parse_pdb",
    "distance",
    "angle",
    "centroid",
    "find_hydrogen_bonds",
    "find_salt_bridges",
    "find_hydrophobic_contacts",
    "find_pi_stacking",
    "find_cation_pi",
    "find_halogen_bonds",
    "find_water_bridges",
    "detect_metal_coordination",
    "generate_pymol_script",
    "save_script_only",
    "detect_molecules",
]
