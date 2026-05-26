"""
Interaction detection subpackage.

Provides all interaction detection functions for protein-ligand analysis.
Split into submodules by interaction type for maintainability.
"""

from .hbonds import find_hydrogen_bonds
from .salt_bridges import find_salt_bridges
from .hydrophobic import find_hydrophobic_contacts
from .pi_stacking import find_pi_stacking
from .cation_pi import find_cation_pi
from .halogen import find_halogen_bonds
from .water_bridges import find_water_bridges, find_xwh_bridges
from .metal import detect_metal_coordination
from .validation import validate_interactions, get_interaction_summary

__all__ = [
    'find_hydrogen_bonds',
    'find_salt_bridges',
    'find_hydrophobic_contacts',
    'find_pi_stacking',
    'find_cation_pi',
    'find_halogen_bonds',
    'find_water_bridges',
    'find_xwh_bridges',
    'detect_metal_coordination',
    'validate_interactions',
    'get_interaction_summary',
]
