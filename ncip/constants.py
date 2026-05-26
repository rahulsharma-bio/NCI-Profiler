"""
Constants and parameters for protein-ligand interaction analysis.

This module contains all the cutoff values, atom classifications, and 
interaction parameters used throughout the analysis.
"""

from typing import Set, Dict

# ============================================================================
# RESIDUE CLASSIFICATIONS
# ============================================================================

STANDARD_AMINO_ACIDS: Set[str] = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Modified residues
    'MSE', 'SEC', 'PYL', 'HYP', 'SEP', 'TPO', 'PTR',
    # Protonation states
    'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN',
}

# Common metal ions and their typical oxidation states
COMMON_IONS: Set[str] = {
    'ZN', 'MG', 'CA', 'FE', 'MN', 'CU', 'NI', 'CO', 'NA', 'K',
    'CD', 'HG', 'PB', 'SR', 'BA', 'AL', 'MO', 'W', 'V', 'CR',
}

# True biological cofactors (enzyme prosthetic groups, coenzymes)
TRUE_COFACTORS: Set[str] = {
    # Nucleotide cofactors
    'NAD', 'NAP', 'NDP', 'FAD', 'FMN', 'ADP', 'ATP', 'GTP', 'GDP',
    'CTP', 'UTP', 'NAI', 'NAJ', 'AMP', 'GMP', 'SAM', 'SAH',
    # Heme and related
    'HEM', 'HEC', 'HEA', 'HEB', 'HEO', 'HAS', 'HDD', 'HDM',
    # Coenzyme A and related
    'COA', 'ACO', 'CAA',
    # Biotin
    'BTN', 'BIO',
    # Thiamine
    'TPP', 'TDP',
    # Pyridoxal phosphate
    'PLP', 'PMP', 'PXP',
    # Flavins
    'RBF', 'FLA',
    # Folate
    'FOL', 'THF', 'DHF',
}

# Buffers, salts, and crystallization additives (not biologically relevant)
BUFFERS_SALTS: Set[str] = {
    # Common buffers
    'GOL', 'PEG', 'PG4', 'P6G', 'P4G', '1PE', '2PE',  # Polyethylene glycols
    'EDO', 'EOH',  # Ethylene glycol, ethanol
    'DMS', 'DMSO',  # DMSO
    'MPD', 'MRD',  # MPD
    'ACT', 'ACE', 'ACY',  # Acetate
    'TRS', 'TAM', 'MES', 'EPE', 'HEZ',  # Buffer molecules
    'IMD', 'IPA', 'IPH',  # Imidazole
    # Common salts and ions
    'SO4', 'PO4', 'PO3',  # Sulfate, Phosphate
    'CL', 'BR', 'IOD', 'FLC',  # Halides
    'NO3', 'NO2',  # Nitrate/Nitrite
    'SCN', 'AZI',  # Thiocyanate, Azide
    'CIT', 'FLC', 'FMT', 'MLI',  # Citrate, Formate, Malonate
    'BOG', 'LDA', 'LMT', 'OLC',  # Detergents
    'BME', 'DTT', 'TCE',  # Reducing agents
    'NH4', 'NH2',  # Ammonium
}

# Combined set for backward compatibility (includes both cofactors and buffers)
COFACTORS: Set[str] = TRUE_COFACTORS | BUFFERS_SALTS

# Glycans and sugar residues commonly found in PDB structures
GLYCANS: Set[str] = {
    # N-linked glycan core
    'NAG', 'NDG',  # N-Acetyl-D-glucosamine
    'BMA', 'MAN',  # Mannose (beta and alpha)
    'FUC', 'FUL',  # Fucose
    # Common sugars
    'GAL', 'GLA',  # Galactose
    'GLC', 'BGC',  # Glucose
    'SIA', 'SLB',  # Sialic acid
    'NGA', 'A2G',  # N-Acetyl-D-galactosamine
    # Other sugars
    'XYS', 'XYP',  # Xylose
    'GCU', 'BDP',  # Glucuronic acid
    'IDR', 'IDS',  # Iduronic acid
    'RAM', 'RHA',  # Rhamnose
    'ARA', 'ARB',  # Arabinose
    'RIB',         # Ribose
    'FRU',         # Fructose
    'LAT',         # Lactose
    'MAL',         # Maltose
    'SUC', 'TRE',  # Sucrose, Trehalose
}

# Solvents and crystallization artifacts to always exclude
SOLVENTS: Set[str] = {'HOH', 'WAT', 'H2O', 'DOD', 'D2O'}

# DNA and RNA nucleotide residues (exclude from ligand detection)
NUCLEIC_ACID_RESIDUES: Set[str] = {
    # Standard DNA nucleotides
    'DA', 'DT', 'DG', 'DC', 'DU',
    # Standard RNA nucleotides  
    'A', 'U', 'G', 'C', 'T',
    # Alternative naming conventions
    'ADE', 'THY', 'GUA', 'CYT', 'URA',
    # Modified nucleotides (common)
    'PSU', '5MC', '1MA', '2MG', '7MG', 'M2G', 'OMC', 'OMG', 'YG',
    'H2U', 'DHU',  # Dihydrouridine
    'I', 'DI',      # Inosine
    # Nucleotide monophosphates as part of nucleic acids
    'AMP', 'GMP', 'CMP', 'UMP', 'TMP', 'DMP',
    # 5' terminal phosphate forms
    '5GP', '5CM', '5MU',
    # DNA/RNA backbone and related
    'N', 'DN',
}

# ============================================================================
# METAL COORDINATION PARAMETERS
# ============================================================================

# Typical metal-ligand bond distances (in Angstroms)
METAL_COORD_DISTANCES: Dict[str, tuple] = {
    # Metal: (min_dist, max_dist, ideal_dist)
    'ZN': (1.9, 2.6, 2.1),
    'MG': (1.9, 2.5, 2.1),
    'CA': (2.2, 2.8, 2.4),
    'FE': (1.8, 2.6, 2.0),
    'MN': (1.9, 2.6, 2.2),
    'CU': (1.8, 2.5, 2.0),
    'NI': (1.8, 2.4, 2.0),
    'CO': (1.8, 2.4, 2.1),
    'NA': (2.2, 3.0, 2.4),
    'K':  (2.5, 3.3, 2.8),
    'CD': (2.1, 2.8, 2.4),
    'HG': (2.0, 2.8, 2.4),
    'MO': (1.9, 2.6, 2.1),
    'W':  (1.9, 2.6, 2.1),
    'V':  (1.8, 2.4, 2.0),
    'CR': (1.8, 2.4, 2.0),
}

# Default coordination distance for unknown metals
DEFAULT_METAL_COORD_DIST: tuple = (1.8, 2.8, 2.2)

# Coordination number to geometry mapping
COORDINATION_GEOMETRIES: Dict[int, str] = {
    2: 'linear',
    3: 'trigonal planar',
    4: 'tetrahedral/square planar',
    5: 'trigonal bipyramidal/square pyramidal',
    6: 'octahedral',
    7: 'pentagonal bipyramidal',
    8: 'square antiprismatic',
}

# Atoms that can coordinate metals (by element or atom type)
METAL_COORDINATING_ATOMS: Set[str] = {
    'N', 'O', 'S', 'SE',  # Elements
    'NE2', 'ND1', 'NZ', 'NE', 'NH1', 'NH2',  # Nitrogen atoms
    'OD1', 'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH',  # Oxygen atoms
    'SG', 'SD',  # Sulfur atoms
}

# Residues known to coordinate metals
METAL_COORDINATING_RESIDUES: Set[str] = {
    'HIS', 'CYS', 'ASP', 'GLU', 'MET', 'SER', 'THR', 'TYR', 'ASN', 'GLN',
    'HID', 'HIE', 'HIP',  # Histidine protonation states
}

# Protein backbone atoms that can coordinate metals
BACKBONE_METAL_COORDS: Set[str] = {'O', 'N'}

# ============================================================================
# HYDROGEN BOND PARAMETERS
# ============================================================================

# Protein donor atoms (by residue type)
PROTEIN_DONORS: Dict[str, Set[str]] = {
    # Backbone - all residues
    'ALL': {'N'},
    # Sidechains
    'ARG': {'NE', 'NH1', 'NH2'},
    'ASN': {'ND2'},
    'CYS': {'SG'},  # Weak donor
    'GLN': {'NE2'},
    'HIS': {'ND1', 'NE2'},
    'HID': {'ND1'},
    'HIE': {'NE2'},
    'HIP': {'ND1', 'NE2'},
    'LYS': {'NZ'},
    'SER': {'OG'},
    'THR': {'OG1'},
    'TRP': {'NE1'},
    'TYR': {'OH'},
}

# Protein acceptor atoms (by residue type)
PROTEIN_ACCEPTORS: Dict[str, Set[str]] = {
    # Backbone - all residues
    'ALL': {'O'},
    # Sidechains
    'ASN': {'OD1'},
    'ASP': {'OD1', 'OD2'},
    'GLN': {'OE1'},
    'GLU': {'OE1', 'OE2'},
    'HIS': {'ND1', 'NE2'},
    'HID': {'NE2'},
    'HIE': {'ND1'},
    'MET': {'SD'},  # Weak acceptor
    'SER': {'OG'},
    'THR': {'OG1'},
    'TYR': {'OH'},
}

# Hydrogen bond distance cutoffs (Jeffrey 1997; operational thresholds)
HBOND_DISTANCE_MAX: float = 3.5  # Maximum D-A distance (Å)
HBOND_DISTANCE_MIN: float = 2.5  # Minimum D-A distance (avoid clashes)
HBOND_ANGLE_MIN: float = 120.0   # Minimum D-H-A angle (°)

# ============================================================================
# SALT BRIDGE PARAMETERS
# ============================================================================

# Positively charged atoms in proteins
PROTEIN_POSITIVE: Dict[str, Set[str]] = {
    'ARG': {'NE', 'NH1', 'NH2', 'CZ'},  # Guanidinium
    'LYS': {'NZ'},                       # Amino
    'HIS': {'ND1', 'NE2'},               # Imidazole (protonated)
    'HIP': {'ND1', 'NE2'},
}

# Negatively charged atoms in proteins
PROTEIN_NEGATIVE: Dict[str, Set[str]] = {
    'ASP': {'OD1', 'OD2', 'CG'},  # Carboxylate
    'GLU': {'OE1', 'OE2', 'CD'},  # Carboxylate
}

# Salt bridge cutoffs
# Barlow & Thornton (1983) + 1.5 Å buffer for charge-centre distances
SALT_BRIDGE_DIST_MAX: float = 5.5  # Maximum charge-center distance (Å)

# ============================================================================
# HYDROPHOBIC CONTACT PARAMETERS
# ============================================================================

# Hydrophobic atoms in proteins (by residue)
HYDROPHOBIC_RESIDUES: Dict[str, Set[str]] = {
    'ALA': {'CB'},
    'VAL': {'CB', 'CG1', 'CG2'},
    'LEU': {'CB', 'CG', 'CD1', 'CD2'},
    'ILE': {'CB', 'CG1', 'CG2', 'CD1'},
    'MET': {'CB', 'CG', 'SD', 'CE'},
    'PHE': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'TYR': {'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2'},  # Ring only
    'TRP': {'CB', 'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'HIS': {'CB', 'CG', 'CD2'},  # Hydrophobic carbons only (CD2 is purely C-bonded)
    'PRO': {'CB', 'CG', 'CD'},
}

# Hydrophobic contact cutoffs
HYDROPHOBIC_DIST_MIN: float = 3.3  # Minimum distance (avoid clashes)
HYDROPHOBIC_DIST_MAX: float = 4.5  # Maximum distance (Bissantz et al. 2010)

# ============================================================================
# π-π STACKING PARAMETERS
# ============================================================================

# Aromatic ring atoms in proteins
AROMATIC_RESIDUES: Dict[str, Set[str]] = {
    'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'TRP': {'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'HIS': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HID': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HIE': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'HIP': {'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
}

# Pi stacking cutoffs (Calinsky & Levy 2024; Li et al. 2025; McGaughey et al. 1998)
PI_STACK_DIST_MAX: float = 5.0      # Maximum centroid-centroid distance (Å)
PI_STACK_ANGLE_PARALLEL: float = 30.0  # Max interplane angle for parallel stacking (P < 30°)
PI_STACK_ANGLE_TSHAPE: float = 60.0    # Min interplane angle for T-shaped (P > 60°)
PI_STACK_ELEVATION_MIN: float = 65.0   # Min elevation angle (Tθ) for true stacking (Calinsky & Levy)
PI_STACK_VERT_DIST_MIN: float = 3.0    # Min vertical distance for parallel stacking
PI_STACK_VERT_DIST_MAX: float = 4.2    # Max vertical distance for parallel stacking

# ============================================================================
# CATION-π PARAMETERS
# ============================================================================

# Cationic groups in proteins
PROTEIN_CATIONS: Dict[str, Set[str]] = {
    'ARG': {'CZ', 'NE', 'NH1', 'NH2'},  # Guanidinium
    'LYS': {'NZ'},                       # Amino
}

# Cation-π cutoffs
CATION_PI_DIST_MAX: float = 6.0  # Maximum distance to ring centroid

# ============================================================================
# HALOGEN BOND PARAMETERS
# ============================================================================

# Halogens that form halogen bonds
HALOGENS: Set[str] = {'CL', 'BR', 'I', 'F'}

# Halogen bond cutoffs (Auffinger et al. 2004; Jiang et al. 2005)
HALOGEN_BOND_DIST_MAX: float = 4.0  # Maximum X...A distance (Å)
HALOGEN_BOND_ANGLE_MIN: float = 140.0  # Minimum C-X...A angle (°)

# ============================================================================
# WATER BRIDGE PARAMETERS
# ============================================================================

WATER_BRIDGE_DIST_MAX: float = 3.5  # Maximum water-protein/ligand distance (Å)
WATER_BRIDGE_ANGLE_MIN: float = 80.0  # Min D-W-A angle (Jiang et al. 2005)
WATER_BRIDGE_ANGLE_MAX: float = 140.0  # Max D-W-A angle (Jiang et al. 2005)

# ============================================================================
# XWH (HALOGEN-WATER-HYDROGEN) BRIDGE PARAMETERS
# ============================================================================

XWH_XBOND_ANGLE_MIN: float = 130.0    # Min C-X···Ow angle (Zhou et al. 2010)
XWH_HBOND_DIST_MAX: float = 3.9       # Max Ow···D/A distance for N,O
XWH_HBOND_ANGLE_MIN: float = 90.0     # Min Ow···D/A-C angle

# ============================================================================
# METAL COORDINATION SIMPLE CUTOFFS
# ============================================================================
# These are the default (unknown-metal) min/max coordination distances.
# Per-metal distances in METAL_COORD_DISTANCES take precedence when available.
METAL_COORD_DIST_MIN: float = 1.8     # Default min coordination distance (Å)
METAL_COORD_DIST_MAX: float = 2.8     # Default max coordination distance (Å)

# ============================================================================
# ELEMENT INFERENCE
# ============================================================================

# Common atom name to element mapping
ATOM_ELEMENT_MAP: Dict[str, str] = {
    'CL': 'CL', 'BR': 'BR', 'FE': 'FE', 'ZN': 'ZN', 'MG': 'MG',
    'CA': 'CA', 'MN': 'MN', 'CU': 'CU', 'NI': 'NI', 'CO': 'CO',
    'SE': 'SE', 'NA': 'NA', 'K': 'K', 'CD': 'CD', 'HG': 'HG',
    'PB': 'PB', 'SR': 'SR', 'BA': 'BA', 'AL': 'AL', 'MO': 'MO',
    'W': 'W', 'V': 'V', 'CR': 'CR',
}

# ============================================================================
# VISUALIZATION COLORS (for PyMOL)
# ============================================================================

# Colors matched to 2D diagram (INTERACTION_LINE_COLORS in diagram.py)
# Using PyMOL color names that match the hex values
INTERACTION_COLORS: Dict[str, str] = {
    'Hydrogen Bonds': 'forest',      # #228B22 Forest Green
    'Salt Bridges': 'magenta',       # #FF00FF Magenta
    'Hydrophobic': 'gray50',         # #808080 Gray
    'π-π Stacking': 'teal',          # #008080 Teal
    'Cation-π': 'orange',            # #FF8C00 Dark Orange
    'π-Cation': 'chocolate',         # #D2691E Chocolate
    'Halogen Bonds': 'purple',       # #9400D3 Dark Violet
    'Water Bridges': 'marine',       # #1E90FF Dodger Blue
    'XWH Bridges': 'slate',          # #7B68EE Medium Slate Blue
    'Metal Coordination': 'copper',  # #B87333 Copper
}

# Line widths for visualization
INTERACTION_LINE_WIDTH: float = 2.0
DASHED_GAP: float = 0.3
