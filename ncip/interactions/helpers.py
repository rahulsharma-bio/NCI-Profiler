"""
Helper functions for interaction detection.

Provides donor/acceptor classification for protein and ligand atoms.
"""

from typing import List
from ..parser import Atom
from ..constants import PROTEIN_DONORS, PROTEIN_ACCEPTORS


def is_protein_donor(atom: Atom) -> bool:
    """Check if a protein atom can act as H-bond donor."""
    if atom.name == 'N':
        return True
    if atom.res_name in PROTEIN_DONORS:
        if atom.name in PROTEIN_DONORS[atom.res_name]:
            return True
    return False


def is_protein_acceptor(atom: Atom) -> bool:
    """Check if a protein atom can act as H-bond acceptor."""
    if atom.name == 'O':
        return True
    if atom.res_name in PROTEIN_ACCEPTORS:
        if atom.name in PROTEIN_ACCEPTORS[atom.res_name]:
            return True
    return False


def is_ligand_donor(atom: Atom) -> bool:
    """Check if a ligand atom can act as H-bond donor."""
    return atom.element in ('N', 'O', 'S')


def is_ligand_acceptor(atom: Atom) -> bool:
    """Check if a ligand atom can act as H-bond acceptor."""
    return atom.element in ('N', 'O', 'S', 'F')
