"""
PDB file parser and data structures for protein-ligand interaction analysis.

This module provides classes and functions for reading PDB files and extracting
atom, ring, and residue information.
"""

import math
import re
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set, Optional
from collections import defaultdict

from .constants import (
    STANDARD_AMINO_ACIDS, COMMON_IONS, COFACTORS, SOLVENTS,
    AROMATIC_RESIDUES, ATOM_ELEMENT_MAP, NUCLEIC_ACID_RESIDUES,
    TRUE_COFACTORS, BUFFERS_SALTS, GLYCANS
)

# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass
class Atom:
    """Represents a single atom from a PDB file."""
    serial: int
    name: str
    res_name: str
    chain_id: str
    res_seq: int
    x: float
    y: float
    z: float
    element: str
    occupancy: float = 1.0
    b_factor: float = 0.0
    alt_loc: str = ''
    
    @property
    def coords(self) -> Tuple[float, float, float]:
        """Return coordinates as a tuple."""
        return (self.x, self.y, self.z)


@dataclass
class Ring:
    """
    Represents an aromatic ring.
    
    For fused ring systems (e.g., TRP, porphyrins), multiple Ring objects
    share the same ring_system_id, allowing selection of the closest ring
    per protein residue interaction (PInteract approach).
    
    Literature: Li et al. (2025) Biomolecules - PInteract algorithm
    """
    atoms: List[Atom]
    centroid: Tuple[float, float, float]
    normal: Tuple[float, float, float]
    res_name: str
    res_seq: int
    chain_id: str
    ring_size: int = 0                    # Number of atoms in ring (5 or 6)
    ring_system_id: Optional[int] = None  # ID for grouping fused rings
    is_aromatic: bool = True              # True if aromatic character verified
    
    @property
    def residue_id(self) -> str:
        """Return residue identifier string."""
        return f"{self.chain_id}:{self.res_name}{self.res_seq}"


@dataclass
class ChargedGroup:
    """
    Represents a charged functional group for salt bridge detection.
    
    Based on proper chemical identification rather than element-only detection.
    Literature: Barlow & Thornton (1983) J Mol Biol
    """
    group_type: str          # e.g., 'carboxylate', 'primary_amine', 'guanidinium'
    charge: int              # +1 or -1
    atoms: List[Atom]        # All atoms in the group
    charge_center: Tuple[float, float, float]  # Geometric centroid for distance calc
    parent_residue: str = "" # Residue identifier
    confidence: float = 1.0  # How confident we are in this assignment
    
    @property
    def is_positive(self) -> bool:
        return self.charge > 0
    
    @property
    def is_negative(self) -> bool:
        return self.charge < 0


@dataclass
class Interaction:
    """
    Represents a detected interaction between protein and ligand.
    
    The confidence score indicates how well the interaction matches
    ideal geometric criteria (0.0 = distance-only, 1.0 = ideal geometry).
    """
    interaction_type: str
    protein_residue: str
    distance: float
    protein_atoms: List[Atom] = field(default_factory=list)
    ligand_atoms: List[Atom] = field(default_factory=list)
    extra_info: Dict = field(default_factory=dict)
    confidence: float = 1.0  # 0.0 to 1.0, higher = more confident


# ============================================================================
# GEOMETRY FUNCTIONS
# ============================================================================

def distance(a1: Atom, a2: Atom) -> float:
    """Calculate Euclidean distance between two atoms."""
    return math.sqrt(
        (a1.x - a2.x)**2 + 
        (a1.y - a2.y)**2 + 
        (a1.z - a2.z)**2
    )


def distance_coords(c1: Tuple[float, float, float], 
                    c2: Tuple[float, float, float]) -> float:
    """Calculate distance between two coordinate tuples."""
    return math.sqrt(
        (c1[0] - c2[0])**2 + 
        (c1[1] - c2[1])**2 + 
        (c1[2] - c2[2])**2
    )


def angle(a1: Atom, a2: Atom, a3: Atom) -> float:
    """Calculate angle in degrees formed by three atoms (a1-a2-a3)."""
    v1 = (a1.x - a2.x, a1.y - a2.y, a1.z - a2.z)
    v2 = (a3.x - a2.x, a3.y - a2.y, a3.z - a2.z)
    
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    
    if mag1 == 0 or mag2 == 0:
        return 0.0
    
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))


def centroid(atoms: List[Atom]) -> Tuple[float, float, float]:
    """Calculate the centroid of a list of atoms."""
    if not atoms:
        return (0.0, 0.0, 0.0)
    n = len(atoms)
    return (
        sum(a.x for a in atoms) / n,
        sum(a.y for a in atoms) / n,
        sum(a.z for a in atoms) / n
    )


def normal_vector(atoms: List[Atom]) -> Tuple[float, float, float]:
    """
    Calculate normal vector to the plane defined by ring atoms.
    Uses first 3 atoms to define the plane.
    """
    if len(atoms) < 3:
        return (0.0, 0.0, 1.0)
    
    # Vectors from atom 0 to atoms 1 and 2
    v1 = (atoms[1].x - atoms[0].x, 
          atoms[1].y - atoms[0].y, 
          atoms[1].z - atoms[0].z)
    v2 = (atoms[2].x - atoms[0].x, 
          atoms[2].y - atoms[0].y, 
          atoms[2].z - atoms[0].z)
    
    # Cross product
    normal = (
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    )
    
    # Normalize
    mag = math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
    if mag == 0:
        return (0.0, 0.0, 1.0)
    
    return (normal[0]/mag, normal[1]/mag, normal[2]/mag)


def angle_between_normals(n1: Tuple[float, float, float], 
                          n2: Tuple[float, float, float]) -> float:
    """Calculate angle between two normal vectors in degrees."""
    dot = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
    # Take absolute value since direction of normal is arbitrary
    cos_angle = max(-1.0, min(1.0, abs(dot)))
    return math.degrees(math.acos(cos_angle))


def angle_coords(c1: Tuple[float, float, float],
                 c2: Tuple[float, float, float],
                 c3: Tuple[float, float, float]) -> float:
    """
    Calculate angle in degrees formed by three coordinate tuples (c1-c2-c3).
    
    Args:
        c1, c2, c3: Coordinate tuples (x, y, z)
        
    Returns:
        Angle in degrees at c2
    """
    v1 = (c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2])
    v2 = (c3[0] - c2[0], c3[1] - c2[1], c3[2] - c2[2])
    
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    
    if mag1 == 0 or mag2 == 0:
        return 0.0
    
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))


# ============================================================================
# CHARGED GROUP DETECTION
# ============================================================================

def find_bonded_atoms(atom: Atom, all_atoms: List[Atom], 
                      bond_threshold: float = 1.9) -> List[Atom]:
    """
    Find all atoms bonded to a given atom based on distance.
    
    Args:
        atom: The central atom
        all_atoms: All atoms to search
        bond_threshold: Maximum distance to consider bonded (Å)
        
    Returns:
        List of bonded atoms
    """
    bonded = []
    for other in all_atoms:
        if other.serial != atom.serial:
            d = distance(atom, other)
            if d < bond_threshold:
                bonded.append(other)
    return bonded


def detect_charged_groups(atoms: List[Atom]) -> Tuple[List['ChargedGroup'], List['ChargedGroup']]:
    """
    Detect charged functional groups in a molecule for salt bridge detection.
    
    This function identifies chemically correct charged groups rather than
    simply using element identity (which leads to false positives).
    
    Literature references:
    - Barlow & Thornton (1983) J Mol Biol - Ion pairs in proteins
    - PLIP (Salentin et al. 2015) - Charged group patterns
    
    Positive groups detected:
    - Primary amines (R-NH3+): N with 3+ H or 4 connections
    - Secondary amines (R2-NH2+): N with 2H, 2C
    - Tertiary amines (R3-NH+): N with 1H, 3C
    - Quaternary ammonium (R4-N+): N with 4C, no H
    - Guanidinium: Central C bonded to 3 N
    - Amidinium: C(=NR)(NR2) pattern
    - Imidazolium: 5-membered ring with 2 N (protonated)
    
    Negative groups detected:
    - Carboxylate (R-COO-): C bonded to 2 O with equal C-O distances (~1.25Å)
    - Sulfonate (R-SO3-): S bonded to 3 O
    - Phosphate (R-PO4-/R-PO3-): P bonded to 3-4 O
    - Sulfate (R-O-SO3-): Similar to sulfonate
    
    Args:
        atoms: List of atoms in the molecule
        
    Returns:
        Tuple of (positive_groups, negative_groups)
    """
    positive_groups = []
    negative_groups = []
    
    # Group atoms by element for faster lookup
    atoms_by_element = defaultdict(list)
    for a in atoms:
        atoms_by_element[a.element].append(a)
    
    hydrogens = atoms_by_element.get('H', [])
    carbons = atoms_by_element.get('C', [])
    nitrogens = atoms_by_element.get('N', [])
    oxygens = atoms_by_element.get('O', [])
    sulfurs = atoms_by_element.get('S', [])
    phosphorus = atoms_by_element.get('P', [])
    
    # Build residue identifier if available
    def get_residue_id(atom_list: List[Atom]) -> str:
        if atom_list:
            a = atom_list[0]
            return f"{a.chain_id}:{a.res_name}{a.res_seq}"
        return ""
    
    # =========================================================================
    # DETECT NEGATIVE GROUPS
    # =========================================================================
    
    # 1. Carboxylate detection (R-COO-)
    # A true carboxylate has C bonded to 2 O with nearly equal distances (~1.25Å)
    # A protonated COOH has C=O (~1.2Å) and C-OH (~1.35Å)
    for c in carbons:
        bonded = find_bonded_atoms(c, atoms, bond_threshold=1.7)
        bonded_oxygens = [a for a in bonded if a.element == 'O']
        
        if len(bonded_oxygens) >= 2:
            # Check if it's a carboxylate (deprotonated) vs COOH (protonated)
            # Carboxylate: both C-O distances ~1.25Å
            # COOH: one C=O ~1.2Å, one C-OH ~1.35Å
            o1, o2 = bonded_oxygens[0], bonded_oxygens[1]
            d1 = distance(c, o1)
            d2 = distance(c, o2)
            
            # Carboxylate criterion: |d1 - d2| < 0.15 Å and both ~1.25Å
            if abs(d1 - d2) < 0.15 and 1.15 < d1 < 1.35 and 1.15 < d2 < 1.35:
                charge_center = centroid([o1, o2])
                group = ChargedGroup(
                    group_type='carboxylate',
                    charge=-1,
                    atoms=[c, o1, o2],
                    charge_center=charge_center,
                    parent_residue=get_residue_id([c]),
                    confidence=0.95 if abs(d1 - d2) < 0.10 else 0.8
                )
                negative_groups.append(group)
    
    # 2. Sulfonate detection (R-SO3-)
    for s in sulfurs:
        bonded = find_bonded_atoms(s, atoms, bond_threshold=1.9)
        bonded_oxygens = [a for a in bonded if a.element == 'O']
        
        if len(bonded_oxygens) >= 3:
            charge_center = centroid(bonded_oxygens)
            group = ChargedGroup(
                group_type='sulfonate',
                charge=-1,
                atoms=[s] + bonded_oxygens[:3],
                charge_center=charge_center,
                parent_residue=get_residue_id([s]),
                confidence=0.9
            )
            negative_groups.append(group)
    
    # 3. Phosphate detection (R-PO4- / R-PO3-)
    for p in phosphorus:
        bonded = find_bonded_atoms(p, atoms, bond_threshold=2.0)
        bonded_oxygens = [a for a in bonded if a.element == 'O']
        
        if len(bonded_oxygens) >= 3:
            charge_center = centroid(bonded_oxygens)
            group = ChargedGroup(
                group_type='phosphate',
                charge=-1,
                atoms=[p] + bonded_oxygens,
                charge_center=charge_center,
                parent_residue=get_residue_id([p]),
                confidence=0.9
            )
            negative_groups.append(group)
    
    # =========================================================================
    # DETECT POSITIVE GROUPS
    # =========================================================================
    
    # Check if molecule has explicit hydrogen atoms — most PDB files do NOT.
    # When H atoms are absent, we infer protonation from heavy-atom geometry
    # (bond distances, bond angles, hybridization) instead of H-counting.
    has_explicit_h = len(hydrogens) > 0
    
    # Track nitrogens already assigned to groups
    assigned_nitrogens = set()
    
    # 1. Guanidinium detection (like ARG)
    # Central C bonded to 3 N atoms — does NOT require H atoms
    for c in carbons:
        bonded = find_bonded_atoms(c, atoms, bond_threshold=1.6)
        bonded_nitrogens = [a for a in bonded if a.element == 'N']
        
        if len(bonded_nitrogens) >= 3:
            # This is a guanidinium group
            charge_center = centroid([c] + bonded_nitrogens)
            group = ChargedGroup(
                group_type='guanidinium',
                charge=+1,
                atoms=[c] + bonded_nitrogens,
                charge_center=charge_center,
                parent_residue=get_residue_id([c]),
                confidence=0.95
            )
            positive_groups.append(group)
            for n in bonded_nitrogens:
                assigned_nitrogens.add(n.serial)
    
    # 2. Amidinium detection C(=NR)(NR2)
    # Carbon bonded to exactly 2 N atoms with resonance-stabilized C-N bonds.
    # Works with AND without explicit H atoms.
    # Examples: benzamidinium (MID/NAPAP), formamidinium
    for c in carbons:
        bonded = find_bonded_atoms(c, atoms, bond_threshold=1.6)
        bonded_nitrogens = [a for a in bonded if a.element == 'N' 
                           and a.serial not in assigned_nitrogens]
        
        if len(bonded_nitrogens) == 2:
            n1, n2 = bonded_nitrogens
            d1 = distance(c, n1)
            d2 = distance(c, n2)
            
            # Amidinium C-N bonds: ~1.28-1.40Å, resonance-stabilized (similar lengths)
            if (1.25 < d1 < 1.45) and (1.25 < d2 < 1.45):
                # Exclude urea/carbonyl: C must NOT be bonded to O
                bonded_oxygens = [a for a in bonded if a.element == 'O']
                if bonded_oxygens:
                    continue  # This is likely urea C(=O)(NR)(NR) or similar
                
                # Check that this C is sp2-like (3 connections including the 2 N's)
                bonded_heavy = [a for a in bonded if a.element != 'H']
                
                n1_bonded = find_bonded_atoms(n1, atoms, bond_threshold=1.5)
                n2_bonded = find_bonded_atoms(n2, atoms, bond_threshold=1.5)
                
                # Reject if either N is an aromatic ring nitrogen.
                # Aromatic ring N: exactly 2 C neighbors at aromatic bond distance
                # (< 1.40 Å) with sp2-like angle (> 115°). Such N's are part of
                # heteroaromatic rings (pyridine, pyrimidine, etc.), not amidinium.
                has_ring_n = False
                for n_check, n_bonded_check in [(n1, n1_bonded), (n2, n2_bonded)]:
                    c_neighbors = [a for a in n_bonded_check if a.element == 'C']
                    if len(c_neighbors) == 2:
                        d_a = distance(n_check, c_neighbors[0])
                        d_b = distance(n_check, c_neighbors[1])
                        if d_a < 1.40 and d_b < 1.40:
                            ang = angle(c_neighbors[0], n_check, c_neighbors[1])
                            if ang > 115:
                                has_ring_n = True
                                break
                if has_ring_n:
                    continue
                
                # Determine if we can confirm amidinium
                is_amidinium = False
                
                if has_explicit_h:
                    # With H: at least one N must have H attached
                    has_h1 = any(a.element == 'H' for a in n1_bonded)
                    has_h2 = any(a.element == 'H' for a in n2_bonded)
                    if has_h1 or has_h2:
                        is_amidinium = True
                else:
                    # Without H: infer from heavy-atom connectivity.
                    # Amidinium N's are terminal or near-terminal.
                    # Key distinction from ring N: in a true amidinium, at least one N
                    # has only 1 heavy-atom neighbor (the central C), leaving room for 
                    # implicit H's. Ring N's (pyrimidine, etc.) have ≥2 heavy neighbors.
                    n1_heavy = [a for a in n1_bonded if a.element != 'H']
                    n2_heavy = [a for a in n2_bonded if a.element != 'H']
                    # Exclude the central C itself to count OTHER heavy connections
                    n1_other_heavy = [a for a in n1_heavy if a.serial != c.serial]
                    n2_other_heavy = [a for a in n2_heavy if a.serial != c.serial]
                    # At least one N must have ≤1 other heavy connection (terminal-ish)
                    # If both have ≥2 other connections, they're deeply embedded → not amidinium
                    if len(n1_other_heavy) <= 1 or len(n2_other_heavy) <= 1:
                        # Additional check: both C-N distances should be similar
                        # (resonance-stabilized), |d1-d2| < 0.12Å
                        if abs(d1 - d2) < 0.12:
                            # Ring closure check: if both N's have 1 other heavy neighbor,
                            # verify those "other" atoms don't connect back to each other
                            # (which would indicate a ring like pyrimidine, imidazole, etc.)
                            is_ring = False
                            if len(n1_other_heavy) == 1 and len(n2_other_heavy) == 1:
                                a1 = n1_other_heavy[0]
                                a2 = n2_other_heavy[0]
                                # Direct bond between them → 5-membered ring (C-N1-A1-A2-N2)
                                if distance(a1, a2) < 1.7:
                                    is_ring = True
                                else:
                                    # Shared neighbor → 6-membered ring (C-N1-A1-X-A2-N2)
                                    a1_nb = find_bonded_atoms(a1, atoms, bond_threshold=1.6)
                                    a2_nb = find_bonded_atoms(a2, atoms, bond_threshold=1.6)
                                    a1_serials = {a.serial for a in a1_nb}
                                    a2_serials = {a.serial for a in a2_nb}
                                    excluded = {c.serial, n1.serial, n2.serial}
                                    if (a1_serials & a2_serials) - excluded:
                                        is_ring = True
                            if not is_ring:
                                is_amidinium = True
                            
                if is_amidinium:
                    charge_center = centroid([c, n1, n2])
                    group = ChargedGroup(
                        group_type='amidinium',
                        charge=+1,
                        atoms=[c, n1, n2],
                        charge_center=charge_center,
                        parent_residue=get_residue_id([c]),
                        confidence=0.80 if has_explicit_h else 0.72
                    )
                    positive_groups.append(group)
                    assigned_nitrogens.add(n1.serial)
                    assigned_nitrogens.add(n2.serial)
    
    # 3. Amine detection — handles both with-H and without-H PDB files.
    # When H atoms are absent (standard for PDB), uses heavy-atom connectivity
    # and bond geometry to infer sp3 hybridization (protonatable).
    for n in nitrogens:
        if n.serial in assigned_nitrogens:
            continue
        
        bonded = find_bonded_atoms(n, atoms, bond_threshold=1.6)
        bonded_h = [a for a in bonded if a.element == 'H']
        bonded_c = [a for a in bonded if a.element == 'C']
        bonded_heavy = [a for a in bonded if a.element != 'H']
        bonded_other = [a for a in bonded if a.element not in ('H', 'C')]
        
        # Skip if bonded to electron-withdrawing groups (amides, etc.)
        # Amide N: bonded to C that is also bonded to O (C=O)
        is_amide = False
        for c in bonded_c:
            c_bonded = find_bonded_atoms(c, atoms, bond_threshold=1.7)
            c_oxygens = [a for a in c_bonded if a.element == 'O']
            if c_oxygens:
                # Check for C=O (carbonyl)
                for o in c_oxygens:
                    d_co = distance(c, o)
                    if d_co < 1.30:  # Double bond distance
                        is_amide = True
                        break
            if is_amide:
                break
        
        if is_amide:
            continue
        
        # Skip nitrile N (triple bonded to C, very short C≡N ~1.15Å)
        is_nitrile = False
        for c in bonded_c:
            d_cn = distance(n, c)
            if d_cn < 1.20:  # Triple bond distance
                is_nitrile = True
                break
        
        if is_nitrile:
            continue
        
        # Skip N bonded to electron-withdrawing atoms other than C/H
        # (e.g., N-O in nitro groups, N-S in sulfonamides)
        # These N's are poor bases and unlikely to be protonated
        is_ew_bonded = False
        for other in bonded_other:
            if other.element in ('O', 'S', 'P'):
                is_ew_bonded = True
                break
        if is_ew_bonded:
            continue
        
        total_connections = len(bonded)
        num_h = len(bonded_h)
        num_c = len(bonded_c)
        num_heavy = len(bonded_heavy)
        
        # Classify amine type
        amine_type = None
        confidence = 0.7
        
        if has_explicit_h:
            # ---- Classification with explicit H atoms ----
            if num_h >= 3 or (num_h >= 2 and total_connections >= 4):
                # Primary amine (R-NH3+)
                amine_type = 'primary_amine'
                confidence = 0.85
            elif num_h == 2 and num_c >= 1:
                # Secondary amine (R2-NH2+)
                amine_type = 'secondary_amine'
                confidence = 0.80
            elif num_h == 1 and num_c >= 2:
                # Tertiary amine (R3-NH+) - protonatable
                amine_type = 'tertiary_amine'
                confidence = 0.70
            elif num_h == 0 and num_c >= 4:
                # Quaternary ammonium (R4-N+)
                amine_type = 'quaternary_ammonium'
                confidence = 0.90
            elif num_h >= 1 and total_connections >= 3:
                # Generic protonatable nitrogen
                amine_type = 'protonatable_nitrogen'
                confidence = 0.60
        else:
            # ---- H-free classification (standard PDB files) ----
            # Infer hybridization from geometry: sp3 ~109.5° vs sp2 ~120°
            
            # Skip sp2 aromatic ring N (pyridine, pyrimidine-like):
            # N with 2 C neighbors, both with aromatic-length C-N bonds, sp2 angles
            if num_heavy == 2 and num_c == 2:
                d1 = distance(n, bonded_c[0])
                d2 = distance(n, bonded_c[1])
                bond_angle = angle(bonded_c[0], n, bonded_c[1])
                # Aromatic C-N: both < 1.40Å and angle > 115° (sp2)
                if d1 < 1.40 and d2 < 1.40 and bond_angle > 115:
                    continue  # Aromatic ring N — pKa too low, skip
            
            if num_heavy == 1 and num_c == 1:
                # Terminal N with single C bond → primary amine (R-NH₂/NH₃⁺)
                # Excludes nitrile (caught above) and amide (caught above)
                # Check C-N distance: single bond ~1.47Å; if very short (<1.32Å)
                # it could be an imine (C=N) — weaker base, lower confidence
                d_cn = distance(n, bonded_c[0])
                if d_cn > 1.35:
                    amine_type = 'primary_amine'
                    confidence = 0.75
                elif d_cn > 1.25:
                    # Could be imine (C=N) — still potentially protonatable
                    amine_type = 'primary_amine'
                    confidence = 0.55
                    
            elif num_heavy >= 3:
                # N with ≥3 heavy connections → tertiary/quaternary amine
                # Compute average bond angle to determine sp3 vs sp2
                angles_list = []
                for i in range(len(bonded_heavy)):
                    for j in range(i + 1, len(bonded_heavy)):
                        a_val = angle(bonded_heavy[i], n, bonded_heavy[j])
                        angles_list.append(a_val)
                avg_angle = sum(angles_list) / len(angles_list) if angles_list else 120.0
                
                if avg_angle < 115:  # sp3 geometry (~109.5°)
                    if num_heavy >= 4:
                        amine_type = 'quaternary_ammonium'
                        confidence = 0.85
                    else:
                        # Tertiary amine (3 heavy connections, sp3)
                        # This catches dimethylamino groups like in 4JRV KJV
                        amine_type = 'tertiary_amine'
                        confidence = 0.65
                        
            elif num_heavy == 2 and num_c >= 1:
                # N with 2 heavy connections — not caught by aromatic check above
                # Could be secondary amine (R₂-NH)
                # Check sp3 from bond angle
                a_val = angle(bonded_heavy[0], n, bonded_heavy[1])
                if a_val < 115:  # sp3
                    amine_type = 'secondary_amine'
                    confidence = 0.60
        
        if amine_type:
            group = ChargedGroup(
                group_type=amine_type,
                charge=+1,
                atoms=[n] + bonded_h,
                charge_center=(n.x, n.y, n.z),
                parent_residue=get_residue_id([n]),
                confidence=confidence
            )
            positive_groups.append(group)
            assigned_nitrogens.add(n.serial)
    
    return positive_groups, negative_groups


def is_valid_carboxylate(carbon: Atom, atoms: List[Atom]) -> Tuple[bool, List[Atom]]:
    """
    Check if a carbon atom is part of a deprotonated carboxylate group.
    
    A true carboxylate (COO-) has nearly equal C-O distances (~1.25Å each).
    A protonated carboxylic acid (COOH) has C=O (~1.2Å) and C-OH (~1.35Å).
    
    Args:
        carbon: The carbon atom to check
        atoms: All atoms in the molecule
        
    Returns:
        Tuple of (is_carboxylate, [oxygen1, oxygen2])
    """
    bonded = find_bonded_atoms(carbon, atoms, bond_threshold=1.7)
    bonded_oxygens = [a for a in bonded if a.element == 'O']
    
    if len(bonded_oxygens) < 2:
        return False, []
    
    o1, o2 = bonded_oxygens[0], bonded_oxygens[1]
    d1 = distance(carbon, o1)
    d2 = distance(carbon, o2)
    
    # Carboxylate: |d1 - d2| < 0.15Å and both around 1.25Å
    is_carboxylate = (abs(d1 - d2) < 0.15 and 
                     1.15 < d1 < 1.35 and 
                     1.15 < d2 < 1.35)
    
    return is_carboxylate, [o1, o2] if is_carboxylate else []


# ============================================================================
# ELEMENT INFERENCE
# ============================================================================

def infer_element(atom_name: str, res_name: str = '') -> str:
    """
    Infer element symbol from atom name.
    
    Args:
        atom_name: PDB atom name (e.g., 'CA', 'CL1', 'FE')
        res_name: Residue name for context
        
    Returns:
        Element symbol (1-2 characters)
    """
    name = atom_name.strip().upper()
    res = res_name.strip().upper()
    
    # IMPORTANT: For standard amino acids, handle common atom names first
    # to avoid confusing CA (C-alpha) with CA (Calcium)
    if res in STANDARD_AMINO_ACIDS:
        # Standard protein atom names - extract element from first letter
        # Common: N, CA, C, O, CB, CG, CD, CE, CZ, NZ, OG, etc.
        element = re.sub(r'[0-9\*\']', '', name)
        if element and element[0] in 'CNOSH':
            return element[0]
    
    # Special case: 'CA' atom name is ambiguous
    # - In proteins/ligands: C-alpha (Carbon)
    # - As ion (residue name 'CA'): Calcium
    if name == 'CA':
        # Only Calcium if the residue itself is named 'CA' (the metal ion)
        if res == 'CA':
            return 'CA'
        else:
            return 'C'  # C-alpha or carbon in ligand
    
    # For metals and halogens, check explicit mappings
    if name in ATOM_ELEMENT_MAP:
        return ATOM_ELEMENT_MAP[name]
    
    # Check 2-letter elements (before stripping)
    two_letter = name[:2] if len(name) >= 2 else name
    if two_letter in ATOM_ELEMENT_MAP:
        return ATOM_ELEMENT_MAP[two_letter]
    
    # Strip numbers and get first letter(s)
    element = re.sub(r'[0-9\*\']', '', name)
    
    if not element:
        return 'C'  # Default to carbon
    
    # Common cases
    if element[0] in 'CNOS':
        return element[0]
    
    # Handle special cases
    if element.startswith('CL'):
        return 'CL'
    if element.startswith('BR'):
        return 'BR'
    if element.startswith('FE'):
        return 'FE'
    if element.startswith('ZN'):
        return 'ZN'
    if element.startswith('MG'):
        return 'MG'
    if element.startswith('SE'):
        return 'SE'
    
    # Default: first character
    return element[0]


# ============================================================================
# RING DETECTION
# ============================================================================

def detect_rings(atoms: List[Atom], residue_name: str) -> List[Ring]:
    """
    Detect aromatic rings in a protein residue.
    
    For fused ring systems (TRP), both rings are returned with the same
    ring_system_id. The π-π detection algorithm will select the closest
    ring per interaction (PInteract approach).
    
    Literature: 
    - Li et al. (2025) Biomolecules - PInteract: for fused rings, use closest ring only
    
    Args:
        atoms: List of atoms in the residue
        residue_name: Three-letter residue code
        
    Returns:
        List of Ring objects found (with ring_system_id for fused systems)
    """
    rings = []
    
    # Check if this is an aromatic residue
    if residue_name in AROMATIC_RESIDUES:
        ring_atoms_names = AROMATIC_RESIDUES[residue_name]
        ring_atoms = [a for a in atoms if a.name in ring_atoms_names]
        
        if len(ring_atoms) >= 4:  # Need at least 4 atoms for a ring
            cent = centroid(ring_atoms)
            norm = normal_vector(ring_atoms)
            
            # TRP has two fused rings - create separate Ring objects with same system ID
            if residue_name == 'TRP':
                # 6-membered ring (benzene portion)
                six_ring_names = {'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'}
                six_ring = [a for a in atoms if a.name in six_ring_names]
                
                # 5-membered ring (pyrrole portion)  
                five_ring_names = {'CG', 'CD1', 'NE1', 'CE2', 'CD2'}
                five_ring = [a for a in atoms if a.name in five_ring_names]
                
                # Use same ring_system_id for both fused rings
                system_id = hash((ring_atoms[0].chain_id, ring_atoms[0].res_seq, 'TRP'))
                
                if len(six_ring) >= 5:
                    cent = centroid(six_ring)
                    norm = normal_vector(six_ring)
                    ring = Ring(
                        atoms=six_ring,
                        centroid=cent,
                        normal=norm,
                        res_name=residue_name,
                        res_seq=six_ring[0].res_seq,
                        chain_id=six_ring[0].chain_id,
                        ring_size=6,
                        ring_system_id=system_id,
                        is_aromatic=True
                    )
                    rings.append(ring)
                
                if len(five_ring) >= 4:
                    cent = centroid(five_ring)
                    norm = normal_vector(five_ring)
                    ring = Ring(
                        atoms=five_ring,
                        centroid=cent,
                        normal=norm,
                        res_name=residue_name,
                        res_seq=five_ring[0].res_seq,
                        chain_id=five_ring[0].chain_id,
                        ring_size=5,
                        ring_system_id=system_id,
                        is_aromatic=True
                    )
                    rings.append(ring)
            else:
                # Single ring residues (PHE, TYR, HIS)
                ring_size = len(ring_atoms)
                ring = Ring(
                    atoms=ring_atoms,
                    centroid=cent,
                    normal=norm,
                    res_name=residue_name,
                    res_seq=ring_atoms[0].res_seq,
                    chain_id=ring_atoms[0].chain_id,
                    ring_size=ring_size,
                    ring_system_id=None,  # No fused system
                    is_aromatic=True
                )
                rings.append(ring)
    
    return rings


def detect_ligand_rings(atoms: List[Atom]) -> List[Ring]:
    """
    Detect aromatic rings in a ligand with improved handling for fused ring systems.
    
    Enhanced algorithm based on literature:
    - Li et al. (2025) Biomolecules - PInteract: for fused rings, group by connectivity
    - Calinsky & Levy (2024) J. Phys. Chem. B - proper aromatic ring validation
    
    Key improvements:
    1. Groups fused rings into ring systems with shared ring_system_id
    2. Validates aromaticity via planarity AND bond length consistency
    3. Only accepts 5-7 membered rings (excludes 4-membered)
    4. Returns ring_system_id for downstream "closest ring" selection
    
    Args:
        atoms: List of ligand atoms
        
    Returns:
        List of Ring objects found, with ring_system_id for fused systems
    """
    rings = []
    
    # Get atoms that can participate in aromatic systems
    aromatic_candidates = [a for a in atoms if a.element in ('C', 'N', 'O', 'S')]
    
    if len(aromatic_candidates) < 5:
        return rings
    
    # Build connectivity based on distance (aromatic bond lengths: 1.25-1.55 Å)
    connectivity = defaultdict(list)
    bond_lengths = {}  # Track bond lengths for aromaticity check
    
    for i, a1 in enumerate(aromatic_candidates):
        for j, a2 in enumerate(aromatic_candidates[i+1:], i+1):
            d = distance(a1, a2)
            if 1.25 < d < 1.55:  # Aromatic/conjugated bond length range
                connectivity[i].append(j)
                connectivity[j].append(i)
                bond_lengths[(min(i,j), max(i,j))] = d
    
    # Find cycles of length 5, 6, or 7 (valid aromatic ring sizes)
    def find_cycles(start: int, current: int, path: List[int], 
                    visited: Set[int], depth: int) -> List[List[int]]:
        """DFS to find cycles of valid sizes."""
        if depth > 7:  # Max 7-membered rings
            return []
        
        cycles = []
        for neighbor in connectivity[current]:
            if neighbor == start and depth >= 5:  # Min 5-membered rings
                cycles.append(path + [current])
            elif neighbor not in visited:
                new_visited = visited | {neighbor}
                cycles.extend(find_cycles(start, neighbor, path + [current], 
                                         new_visited, depth + 1))
        return cycles
    
    # Find all valid rings
    found_rings = set()
    for start in range(len(aromatic_candidates)):
        cycles = find_cycles(start, start, [], {start}, 1)
        for cycle in cycles:
            # Normalize cycle for comparison (start from min index)
            min_idx = cycle.index(min(cycle))
            normalized = tuple(cycle[min_idx:] + cycle[:min_idx])
            if normalized not in found_rings:
                found_rings.add(normalized)
    
    # Helper to check aromaticity (planarity + bond length uniformity)
    def check_aromaticity(ring_atoms: List[Atom], ring_indices: tuple) -> Tuple[bool, bool]:
        """
        Check if a ring is aromatic based on:
        1. Planarity (atoms deviate < 0.5 Å from best-fit plane)
        2. Bond length uniformity (variance indicates resonance)
        
        Returns:
            (is_planar, is_aromatic)
        """
        if len(ring_atoms) < 5:
            return False, False
            
        cent = centroid(ring_atoms)
        norm = normal_vector(ring_atoms)
        
        # Check planarity
        max_deviation = 0
        for atom in ring_atoms:
            v = (atom.x - cent[0], atom.y - cent[1], atom.z - cent[2])
            dist_to_plane = abs(v[0]*norm[0] + v[1]*norm[1] + v[2]*norm[2])
            max_deviation = max(max_deviation, dist_to_plane)
        
        is_planar = max_deviation <= 0.5  # 0.5 Å tolerance
        
        # Check bond length uniformity (aromatic bonds are ~1.38-1.40 Å)
        # Non-aromatic would have alternating single (~1.54 Å) and double (~1.34 Å)
        ring_bonds = []
        indices_list = list(ring_indices)
        for k in range(len(indices_list)):
            i, j = indices_list[k], indices_list[(k+1) % len(indices_list)]
            key = (min(i,j), max(i,j))
            if key in bond_lengths:
                ring_bonds.append(bond_lengths[key])
        
        if ring_bonds:
            avg_bond = sum(ring_bonds) / len(ring_bonds)
            variance = sum((b - avg_bond)**2 for b in ring_bonds) / len(ring_bonds)
            # Aromatic rings have uniform bonds (low variance)
            # Typical aromatic: variance < 0.01, non-aromatic: variance > 0.02
            is_aromatic = variance < 0.015 and 1.35 < avg_bond < 1.45
        else:
            # If we can't measure bonds, assume planar rings are aromatic
            is_aromatic = is_planar
        
        return is_planar, is_aromatic
    
    # Group rings by shared atoms (fused ring detection)
    ring_list = list(found_rings)
    ring_groups = []  # List of sets of ring indices that share atoms
    ring_to_group = {}
    
    for i, ring_i in enumerate(ring_list):
        set_i = set(ring_i)
        found_group = None
        
        for group_idx, group in enumerate(ring_groups):
            for j in group:
                set_j = set(ring_list[j])
                # Rings are fused if they share 2+ atoms (shared edge)
                if len(set_i & set_j) >= 2:
                    found_group = group_idx
                    break
            if found_group is not None:
                break
        
        if found_group is not None:
            ring_groups[found_group].add(i)
            ring_to_group[i] = found_group
        else:
            new_group_idx = len(ring_groups)
            ring_groups.append({i})
            ring_to_group[i] = new_group_idx
    
    # Convert to Ring objects with ring_system_id
    for i, ring_indices in enumerate(ring_list):
        ring_atoms = [aromatic_candidates[idx] for idx in ring_indices]
        ring_size = len(ring_atoms)
        
        # Skip rings outside valid size range
        if ring_size < 5 or ring_size > 7:
            continue
        
        # Check aromaticity
        is_planar, is_aromatic = check_aromaticity(ring_atoms, ring_indices)
        
        if is_planar:  # Only include planar rings
            cent = centroid(ring_atoms)
            norm = normal_vector(ring_atoms)
            
            # Generate ring_system_id from group membership
            group_idx = ring_to_group.get(i, i)
            system_id = hash((ring_atoms[0].chain_id, ring_atoms[0].res_seq, 
                            ring_atoms[0].res_name, group_idx))
            
            ring = Ring(
                atoms=ring_atoms,
                centroid=cent,
                normal=norm,
                res_name=ring_atoms[0].res_name,
                res_seq=ring_atoms[0].res_seq,
                chain_id=ring_atoms[0].chain_id,
                ring_size=ring_size,
                ring_system_id=system_id if len(ring_groups[group_idx]) > 1 else None,
                is_aromatic=is_aromatic
            )
            rings.append(ring)
    
    return rings


# ============================================================================
# COVALENT BOND PARSING
# ============================================================================

def parse_covalent_links(pdb_file: str) -> List[Dict]:
    """
    Parse LINK records from PDB file to detect covalent bonds between
    protein and ligand.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        List of dictionaries containing covalent bond information
    """
    covalent_bonds = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('LINK'):
                # LINK record format:
                # Columns 13-16: Atom name 1
                # Columns 17: Alt loc 1
                # Columns 18-20: Residue name 1
                # Column 22: Chain ID 1
                # Columns 23-26: Residue sequence number 1
                # Columns 43-46: Atom name 2
                # Column 47: Alt loc 2
                # Columns 48-50: Residue name 2
                # Column 52: Chain ID 2
                # Columns 53-56: Residue sequence number 2
                # Columns 74-78: Link distance (if present)
                
                try:
                    atom1_name = line[12:16].strip()
                    res1_name = line[17:20].strip()
                    chain1 = line[21].strip() or 'A'
                    res1_num = int(line[22:26].strip())
                    
                    atom2_name = line[42:46].strip()
                    res2_name = line[47:50].strip()
                    chain2 = line[51].strip() or 'A'
                    res2_num = int(line[52:56].strip())
                    
                    # Try to get distance if present
                    dist = None
                    if len(line) >= 78:
                        try:
                            dist = float(line[73:78].strip())
                        except (ValueError, IndexError):
                            pass
                    
                    # Determine which is protein and which is ligand
                    is_res1_protein = res1_name in STANDARD_AMINO_ACIDS
                    is_res2_protein = res2_name in STANDARD_AMINO_ACIDS
                    
                    if is_res1_protein and not is_res2_protein:
                        covalent_bonds.append({
                            'prot_atom': atom1_name,
                            'prot_res': res1_name,
                            'prot_chain': chain1,
                            'prot_resnum': res1_num,
                            'lig_atom': atom2_name,
                            'lig_res': res2_name,
                            'lig_chain': chain2,
                            'lig_resnum': res2_num,
                            'distance': dist
                        })
                    elif is_res2_protein and not is_res1_protein:
                        covalent_bonds.append({
                            'prot_atom': atom2_name,
                            'prot_res': res2_name,
                            'prot_chain': chain2,
                            'prot_resnum': res2_num,
                            'lig_atom': atom1_name,
                            'lig_res': res1_name,
                            'lig_chain': chain1,
                            'lig_resnum': res1_num,
                            'distance': dist
                        })
                        
                except (ValueError, IndexError):
                    continue
    
    return covalent_bonds


# ============================================================================
# MOLECULE DETECTION
# ============================================================================

def detect_molecules(pdb_file: str) -> list:
    """
    Scan a PDB file and return a categorized table of all non-protein,
    non-solvent molecules found.
    
    Returns:
        List of dicts with keys: chain, resname, resnum, category, atoms
        where category is one of: 'Ligand', 'Cofactor', 'Glycan', 'Buffer/Salt'
    """
    molecules = []
    seen = set()
    atom_counts = {}
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            try:
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_seq = int(line[22:26].strip())
                
                # Skip protein, nucleic acids, solvents, ions
                if res_name in STANDARD_AMINO_ACIDS:
                    continue
                if res_name in NUCLEIC_ACID_RESIDUES:
                    continue
                if res_name in SOLVENTS:
                    continue
                if res_name in COMMON_IONS:
                    continue
                element = line[76:78].strip().upper() if len(line) > 76 else ''
                if element in COMMON_IONS:
                    continue
                
                key = (chain_id, res_name, res_seq)
                atom_counts[key] = atom_counts.get(key, 0) + 1
                
                if key in seen:
                    continue
                seen.add(key)
                
                # Categorize
                if res_name in TRUE_COFACTORS:
                    category = 'Cofactor'
                elif res_name in GLYCANS:
                    category = 'Glycan'
                elif res_name in BUFFERS_SALTS:
                    category = 'Buffer/Salt'
                else:
                    category = 'Ligand'
                
                molecules.append({
                    'chain': chain_id or '-',
                    'resname': res_name,
                    'resnum': res_seq,
                    'category': category,
                })
            except (ValueError, IndexError):
                continue
    
    # Add atom counts
    for mol in molecules:
        key = (mol['chain'] if mol['chain'] != '-' else '', mol['resname'], mol['resnum'])
        mol['atoms'] = atom_counts.get(key, 0)
    
    return molecules


# ============================================================================
# PDB PARSING
# ============================================================================

def parse_pdb(pdb_file: str, 
              ligand_name: str = None,
              ligand_chain: str = None, 
              ligand_resnum: int = None,
              include_cofactors: bool = True,
              exclude_cofactors: bool = False,
              exclude_glycans: bool = False) -> Tuple:
    """
    Parse a PDB file and extract protein, ligand, water, and metal atoms.
    
    Args:
        pdb_file: Path to PDB file
        ligand_name: Optional specific ligand residue name to filter
        ligand_chain: Optional specific ligand chain ID
        ligand_resnum: Optional specific ligand residue number
        include_cofactors: If True, include cofactors as ligands (default True)
        exclude_cofactors: If True, explicitly exclude cofactors
        exclude_glycans: If True, explicitly exclude glycans
        
    Returns:
        Tuple of:
        - protein_atoms: List of Atom
        - ligand_atoms: List of Atom
        - water_atoms: List of Atom
        - metal_atoms: List of Atom
        - ligand_names: Set of ligand residue names found
        - ligand_residues: List of (chain, resname, resnum) tuples
        - has_chains: Boolean indicating if file has chain IDs
    """
    protein_atoms = []
    ligand_atoms = []
    water_atoms = []
    metal_atoms = []
    ligand_names = set()
    ligand_residues = []
    has_chains = False
    
    seen_residues = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            
            try:
                # Parse atom record
                record_type = line[:6].strip()
                # Handle hybrid-36 serial numbers (alphanumeric for atoms > 99999)
                serial_str = line[6:11].strip()
                try:
                    serial = int(serial_str)
                except ValueError:
                    # Hybrid-36 encoding: convert to unique integer > 99999
                    serial = hash(serial_str) % (10**9) + 100000
                name = line[12:16].strip()
                alt_loc = line[16].strip()
                res_name = line[17:20].strip()
                raw_chain_id = line[21].strip()  # Keep original for has_chains check
                chain_id = raw_chain_id or 'A'   # Default to 'A' for internal use
                res_seq = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                # Try to get occupancy and B-factor
                try:
                    occupancy = float(line[54:60].strip())
                except (ValueError, IndexError):
                    occupancy = 1.0
                    
                try:
                    b_factor = float(line[60:66].strip())
                except (ValueError, IndexError):
                    b_factor = 0.0
                
                # Try to get element from columns 77-78
                try:
                    element = line[76:78].strip().upper()
                except IndexError:
                    element = ''
                
                if not element:
                    element = infer_element(name, res_name)
                
                # Check if file uses chain IDs (using original value before default)
                if raw_chain_id:
                    has_chains = True
                
                # Create atom object
                atom = Atom(
                    serial=serial,
                    name=name,
                    res_name=res_name,
                    chain_id=chain_id,
                    res_seq=res_seq,
                    x=x, y=y, z=z,
                    element=element,
                    occupancy=occupancy,
                    b_factor=b_factor,
                    alt_loc=alt_loc
                )
                
                # Classify atom
                if res_name in SOLVENTS:
                    water_atoms.append(atom)
                elif res_name in COMMON_IONS or element in COMMON_IONS:
                    metal_atoms.append(atom)
                elif res_name in STANDARD_AMINO_ACIDS:
                    protein_atoms.append(atom)
                elif res_name in NUCLEIC_ACID_RESIDUES:
                    # Skip DNA/RNA residues - not supported for interaction analysis
                    continue
                else:
                    # Potential ligand - check classification
                    is_buffer = res_name in BUFFERS_SALTS
                    is_cofactor = res_name in TRUE_COFACTORS
                    is_glycan = res_name in GLYCANS
                    
                    # Always exclude buffers/salts unless explicitly requested by name
                    if is_buffer and not ligand_name:
                        continue
                    
                    # Handle cofactor exclusion (included by default)
                    if is_cofactor and (exclude_cofactors or not include_cofactors):
                        continue
                    
                    # Handle glycan exclusion (included by default)
                    if is_glycan and exclude_glycans:
                        continue
                    
                    # Apply ligand filters if specified
                    if ligand_name and res_name != ligand_name:
                        continue
                    if ligand_chain and chain_id != ligand_chain:
                        continue
                    if ligand_resnum and res_seq != ligand_resnum:
                        continue
                    
                    ligand_atoms.append(atom)
                    ligand_names.add(res_name)
                    
                    # Track unique residue instances
                    res_key = (chain_id, res_name, res_seq)
                    if res_key not in seen_residues:
                        seen_residues.add(res_key)
                        ligand_residues.append(res_key)
                        
            except (ValueError, IndexError) as e:
                continue
    
    return (protein_atoms, ligand_atoms, water_atoms, metal_atoms, 
            ligand_names, ligand_residues, has_chains)
