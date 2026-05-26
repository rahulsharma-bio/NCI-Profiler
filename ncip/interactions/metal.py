"""
Metal coordination detection with RMS-based geometry analysis.

Supported geometries:
- CN=2: linear, bent
- CN=3: trigonal planar, trigonal pyramidal, T-shaped
- CN=4: tetrahedral, square planar, seesaw (+ τ₄ parameter)
- CN=5: trigonal bipyramidal, square pyramidal (+ τ₅ parameter)
- CN=6: octahedral, trigonal prismatic
- CN=7: pentagonal bipyramidal, capped octahedral
- CN=8: square antiprismatic, cubic

References:
- Harding (2001) Acta Cryst D57, 401-411
- Yang et al. (2007) Dalton Trans. 955-964 - τ₄
- Addison et al. (1984) J. Chem. Soc. Dalton Trans. 1349-1356 - τ₅
"""

from typing import List, Dict, Tuple, Optional

from ..parser import Atom, Interaction, distance, angle_coords
from ..constants import (
    METAL_COORD_DISTANCES, DEFAULT_METAL_COORD_DIST,
    COORDINATION_GEOMETRIES, METAL_COORDINATING_ATOMS,
    METAL_COORDINATING_RESIDUES, BACKBONE_METAL_COORDS,
    METAL_COORD_DIST_MIN, METAL_COORD_DIST_MAX,
)


# ============================================================================
# IDEAL ANGLE SIGNATURES
# ============================================================================

IDEAL_GEOMETRIES = {
    # CN = 2
    'linear': {
        'coordination_number': 2,
        'ideal_angles': [180.0],
        'description': 'Two ligands opposite each other (180°)'
    },
    'bent': {
        'coordination_number': 2,
        'ideal_angles': [109.5],
        'description': 'Two ligands at angle < 180°'
    },

    # CN = 3
    'trigonal_planar': {
        'coordination_number': 3,
        'ideal_angles': [120.0, 120.0, 120.0],
        'description': 'Three ligands in a plane, 120° apart'
    },
    'trigonal_pyramidal': {
        'coordination_number': 3,
        'ideal_angles': [109.5, 109.5, 109.5],
        'description': 'Three ligands in pyramid arrangement'
    },
    'T-shaped': {
        'coordination_number': 3,
        'ideal_angles': [90.0, 90.0, 180.0],
        'description': 'Three ligands in T arrangement'
    },

    # CN = 4
    'tetrahedral': {
        'coordination_number': 4,
        'ideal_angles': [109.5, 109.5, 109.5, 109.5, 109.5, 109.5],
        'description': 'Four ligands at tetrahedral vertices (Td symmetry)'
    },
    'square_planar': {
        'coordination_number': 4,
        'ideal_angles': [90.0, 90.0, 90.0, 90.0, 180.0, 180.0],
        'description': 'Four ligands in a plane, 90° apart (D4h symmetry)'
    },
    'seesaw': {
        'coordination_number': 4,
        'ideal_angles': [90.0, 90.0, 90.0, 90.0, 120.0, 180.0],
        'description': 'Equatorial-axial distortion from TBP (C2v symmetry)'
    },

    # CN = 5
    'trigonal_bipyramidal': {
        'coordination_number': 5,
        'ideal_angles': [90.0, 90.0, 90.0, 90.0, 90.0, 90.0,
                         120.0, 120.0, 120.0, 180.0],
        'description': 'Trigonal bipyramidal geometry (D3h symmetry)'
    },
    'square_pyramidal': {
        'coordination_number': 5,
        'ideal_angles': [90.0, 90.0, 90.0, 90.0,
                         105.0, 105.0, 105.0, 105.0,
                         150.0, 150.0],
        'description': 'Square pyramidal geometry (C4v symmetry)'
    },

    # CN = 6
    'octahedral': {
        'coordination_number': 6,
        'ideal_angles': [90.0, 90.0, 90.0, 90.0, 90.0, 90.0,
                         90.0, 90.0, 90.0, 90.0, 90.0, 90.0,
                         180.0, 180.0, 180.0],
        'description': 'Octahedral geometry (Oh symmetry)'
    },
    'trigonal_prismatic': {
        'coordination_number': 6,
        'ideal_angles': [82.0, 82.0, 82.0, 82.0, 82.0, 82.0,
                         82.0, 82.0, 82.0, 136.0, 136.0, 136.0,
                         136.0, 136.0, 136.0],
        'description': 'Trigonal prismatic geometry (D3h symmetry)'
    },

    # CN = 7
    'pentagonal_bipyramidal': {
        'coordination_number': 7,
        'ideal_angles': [72.0, 72.0, 72.0, 72.0, 72.0,
                         90.0, 90.0, 90.0, 90.0, 90.0,
                         90.0, 90.0, 90.0, 90.0, 90.0,
                         144.0, 144.0, 144.0, 144.0, 144.0,
                         180.0],
        'description': 'Pentagonal bipyramidal geometry (D5h symmetry)'
    },
    'capped_octahedral': {
        'coordination_number': 7,
        'ideal_angles': [70.5, 70.5, 70.5, 90.0, 90.0, 90.0,
                         90.0, 90.0, 90.0, 90.0, 90.0, 90.0,
                         109.5, 109.5, 109.5, 131.8, 131.8, 131.8,
                         180.0, 180.0, 180.0],
        'description': 'Capped octahedral geometry (C3v symmetry)'
    },

    # CN = 8
    'square_antiprismatic': {
        'coordination_number': 8,
        'ideal_angles': [70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5,
                         74.9, 74.9, 74.9, 74.9, 74.9, 74.9, 74.9, 74.9,
                         118.5, 118.5, 118.5, 118.5, 118.5, 118.5, 118.5, 118.5,
                         141.1, 141.1, 141.1, 141.1],
        'description': 'Square antiprismatic geometry (D4d symmetry)'
    },
    'cubic': {
        'coordination_number': 8,
        'ideal_angles': [70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5, 70.5,
                         70.5, 70.5, 70.5, 70.5,
                         109.5, 109.5, 109.5, 109.5, 109.5, 109.5,
                         109.5, 109.5, 109.5, 109.5, 109.5, 109.5,
                         180.0, 180.0, 180.0, 180.0],
        'description': 'Cubic geometry (Oh symmetry)'
    },
}


# ============================================================================
# RMS-BASED GEOMETRY FITTING
# ============================================================================

def calculate_rms_deviation(observed_angles: List[float],
                            ideal_angles: List[float]) -> Tuple[float, List[Tuple[int, int]]]:
    """
    Calculate RMS deviation between observed and ideal angles using optimal assignment.

    Uses sorted alignment: optimal for minimizing sum of squared differences.
    """
    if len(observed_angles) != len(ideal_angles):
        return float('inf'), []

    if len(observed_angles) == 0:
        return 0.0, []

    obs_sorted = sorted(observed_angles)
    ideal_sorted = sorted(ideal_angles)

    squared_diffs = [(o - i) ** 2 for o, i in zip(obs_sorted, ideal_sorted)]
    rms = (sum(squared_diffs) / len(squared_diffs)) ** 0.5

    obs_indices = sorted(range(len(observed_angles)), key=lambda i: observed_angles[i])
    ideal_indices = sorted(range(len(ideal_angles)), key=lambda i: ideal_angles[i])
    mapping = list(zip(obs_indices, ideal_indices))

    return rms, mapping


def calculate_geometry_confidence(best_rms: float,
                                  next_best_rms: Optional[float],
                                  cn: int) -> float:
    """
    Calculate confidence score for geometry assignment.
    """
    if best_rms < 5:
        base_confidence = 0.95
    elif best_rms < 10:
        base_confidence = 0.85
    elif best_rms < 15:
        base_confidence = 0.70
    elif best_rms < 20:
        base_confidence = 0.55
    elif best_rms < 30:
        base_confidence = 0.40
    else:
        base_confidence = 0.25

    if next_best_rms is not None:
        rms_gap = next_best_rms - best_rms
        if rms_gap > 15:
            discrimination_bonus = 0.10
        elif rms_gap > 8:
            discrimination_bonus = 0.05
        elif rms_gap > 3:
            discrimination_bonus = 0.0
        elif rms_gap > 1:
            discrimination_bonus = -0.10
        else:
            discrimination_bonus = -0.20
    else:
        discrimination_bonus = 0.0

    cn_penalty = 0.0
    if cn >= 7:
        cn_penalty = -0.15
    elif cn >= 5:
        cn_penalty = -0.05

    confidence = base_confidence + discrimination_bonus + cn_penalty
    return max(0.1, min(0.99, confidence))


# ============================================================================
# TAU PARAMETERS
# ============================================================================

def calculate_tau_4(angles: List[float]) -> Optional[float]:
    """
    Calculate τ₄ parameter for 4-coordinate metal centers.

    τ₄ = (360° - (α + β)) / 141°

    Reference: Yang et al. (2007) Dalton Trans., 955-964.
    """
    if len(angles) != 6:
        return None

    sorted_angles = sorted(angles, reverse=True)
    alpha = sorted_angles[0]
    beta = sorted_angles[1]

    tau_4 = (360.0 - (alpha + beta)) / 141.0
    return round(tau_4, 3)


def calculate_tau_5(angles: List[float]) -> Optional[float]:
    """
    Calculate τ₅ parameter for 5-coordinate metal centers.

    τ₅ = (β - α) / 60°

    Reference: Addison et al. (1984) J. Chem. Soc., Dalton Trans., 1349-1356.
    """
    if len(angles) != 10:
        return None

    sorted_angles = sorted(angles, reverse=True)
    beta = sorted_angles[0]
    alpha = sorted_angles[1]

    tau_5 = (beta - alpha) / 60.0
    return round(tau_5, 3)


def interpret_tau_parameter(tau_value: Optional[float], tau_type: str) -> str:
    """Provide human-readable interpretation of τ parameter."""
    if tau_value is None:
        return "unknown"

    if tau_type == 'tau_4':
        if tau_value >= 0.85:
            return "tetrahedral"
        elif tau_value >= 0.65:
            return "distorted tetrahedral / trigonal pyramidal"
        elif tau_value >= 0.35:
            return "seesaw / intermediate"
        elif tau_value >= 0.15:
            return "distorted square planar"
        else:
            return "square planar"

    elif tau_type == 'tau_5':
        if tau_value >= 0.85:
            return "trigonal bipyramidal"
        elif tau_value >= 0.65:
            return "distorted trigonal bipyramidal"
        elif tau_value >= 0.35:
            return "intermediate 5-coordinate"
        elif tau_value >= 0.15:
            return "distorted square pyramidal"
        else:
            return "square pyramidal"

    return "unknown"


# ============================================================================
# MAIN RMS GEOMETRY DETERMINATION
# ============================================================================

def determine_coordination_geometry_rms(
    metal_pos: Tuple[float, float, float],
    coord_positions: List[Tuple[float, float, float]],
    allowed_geometries: Optional[List[str]] = None
) -> Dict:
    """
    Determine metal coordination geometry using RMS fitting to ideal angles.
    """
    cn = len(coord_positions)

    result = {
        'geometry': 'unknown',
        'rms': float('inf'),
        'confidence': 0.0,
        'coordination_number': cn,
        'observed_angles': [],
        'all_fits': {},
        'tau_parameter': None,
        'tau_type': None,
        'tau_interpretation': None,
        'next_best': None
    }

    if cn < 2:
        result['geometry'] = f'coordination_{cn}'
        result['confidence'] = 0.5
        return result

    if cn > 8:
        result['geometry'] = f'high_coordination_{cn}'
        result['confidence'] = 0.3
        return result

    observed_angles = []
    for i in range(len(coord_positions)):
        for j in range(i + 1, len(coord_positions)):
            ang = angle_coords(coord_positions[i], metal_pos, coord_positions[j])
            observed_angles.append(ang)

    result['observed_angles'] = [round(a, 1) for a in observed_angles]

    candidates = {}
    for geom_name, geom_data in IDEAL_GEOMETRIES.items():
        if geom_data['coordination_number'] != cn:
            continue
        if geom_data['ideal_angles'] is None:
            continue
        if allowed_geometries and geom_name not in allowed_geometries:
            continue
        candidates[geom_name] = geom_data

    if not candidates:
        result['geometry'] = f'coordination_{cn}'
        result['confidence'] = 0.4
        return result

    all_fits = {}
    for geom_name, geom_data in candidates.items():
        rms, mapping = calculate_rms_deviation(observed_angles, geom_data['ideal_angles'])
        all_fits[geom_name] = {
            'rms': rms,
            'mapping': mapping,
            'description': geom_data['description']
        }

    result['all_fits'] = {k: round(v['rms'], 2) for k, v in all_fits.items()}

    sorted_fits = sorted(all_fits.items(), key=lambda x: x[1]['rms'])
    best_geom, best_data = sorted_fits[0]

    result['geometry'] = best_geom
    result['rms'] = round(best_data['rms'], 2)

    if len(sorted_fits) > 1:
        next_geom, next_data = sorted_fits[1]
        result['next_best'] = {
            'geometry': next_geom,
            'rms': round(next_data['rms'], 2)
        }

    result['confidence'] = calculate_geometry_confidence(
        best_rms=best_data['rms'],
        next_best_rms=sorted_fits[1][1]['rms'] if len(sorted_fits) > 1 else None,
        cn=cn
    )

    if cn == 4:
        tau_4 = calculate_tau_4(observed_angles)
        result['tau_parameter'] = tau_4
        result['tau_type'] = 'tau_4'
        result['tau_interpretation'] = interpret_tau_parameter(tau_4, 'tau_4')
    elif cn == 5:
        tau_5 = calculate_tau_5(observed_angles)
        result['tau_parameter'] = tau_5
        result['tau_type'] = 'tau_5'
        result['tau_interpretation'] = interpret_tau_parameter(tau_5, 'tau_5')

    return result


def determine_coordination_geometry(metal_pos: Tuple[float, float, float],
                                    coord_positions: List[Tuple[float, float, float]]) -> Tuple[str, float]:
    """Backward-compatible wrapper around RMS-based geometry determination."""
    result = determine_coordination_geometry_rms(metal_pos, coord_positions)
    return result['geometry'], result['confidence']


# ============================================================================
# METAL COORDINATION DETECTION
# ============================================================================

def detect_metal_coordination(metal_atoms: List[Atom],
                               protein_atoms: List[Atom],
                               ligand_atoms: List[Atom],
                               water_atoms: List[Atom]) -> Dict[str, List[Interaction]]:
    """
    Detect metal coordination interactions with RMS-based geometry analysis.

    Args:
        metal_atoms: List of metal ion atoms
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms
        water_atoms: List of water atoms

    Returns:
        Dictionary mapping metal IDs to lists of coordination Interactions
    """
    results = {}

    for metal in metal_atoms:
        metal_id = f"{metal.chain_id}:{metal.element}{metal.res_seq}"
        coordinations = []
        coord_positions = []

        element = metal.element.upper()
        if element in METAL_COORD_DISTANCES:
            min_dist, max_dist, ideal_dist = METAL_COORD_DISTANCES[element]
        else:
            min_dist = METAL_COORD_DIST_MIN
            max_dist = METAL_COORD_DIST_MAX
            ideal_dist = (min_dist + max_dist) / 2

        metal_pos = (metal.x, metal.y, metal.z)

        # Check protein coordination
        for pa in protein_atoms:
            can_coord = (pa.element in METAL_COORDINATING_ATOMS or
                        pa.name in METAL_COORDINATING_ATOMS)
            if not can_coord:
                continue

            d = distance(metal, pa)
            if min_dist <= d <= max_dist:
                is_backbone = pa.name in BACKBONE_METAL_COORDS
                coord_type = "backbone" if is_backbone else "sidechain"
                dist_quality = 1.0 - abs(d - ideal_dist) / (max_dist - ideal_dist)

                inter = Interaction(
                    interaction_type="Metal Coordination",
                    protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                    distance=d,
                    protein_atoms=[pa],
                    ligand_atoms=[metal],
                    confidence=max(0.5, dist_quality),
                    extra_info={
                        'metal': element,
                        'metal_id': metal_id,
                        'metal_atom_name': metal.name,
                        'metal_atom_serial': getattr(metal, 'serial', None),
                        'coord_atom': pa.name,
                        'coord_type': coord_type,
                        'source': 'protein_backbone' if is_backbone else 'protein',
                        'ideal_dist': f"{ideal_dist:.2f}"
                    }
                )
                coordinations.append(inter)
                coord_positions.append((pa.x, pa.y, pa.z))

        # Check ligand coordination
        for la in ligand_atoms:
            can_coord = la.element in ('N', 'O', 'S', 'SE')
            if not can_coord:
                continue

            d = distance(metal, la)
            if min_dist <= d <= max_dist:
                dist_quality = 1.0 - abs(d - ideal_dist) / (max_dist - ideal_dist)

                inter = Interaction(
                    interaction_type="Metal Coordination",
                    protein_residue=f"{la.chain_id}:{la.res_name}{la.res_seq}",
                    distance=d,
                    protein_atoms=[],
                    ligand_atoms=[la],
                    confidence=max(0.5, dist_quality),
                    extra_info={
                        'metal': element,
                        'metal_id': metal_id,
                        'metal_atom_name': metal.name,
                        'metal_atom_serial': getattr(metal, 'serial', None),
                        'coord_atom': la.name,
                        'coord_atom_serial': getattr(la, 'serial', None),
                        'coord_type': la.element,
                        'source': 'ligand',
                        'ideal_dist': f"{ideal_dist:.2f}"
                    }
                )
                coordinations.append(inter)
                coord_positions.append((la.x, la.y, la.z))

        # Check water coordination
        for wa in water_atoms:
            d = distance(metal, wa)
            if min_dist <= d <= max_dist:
                dist_quality = 1.0 - abs(d - ideal_dist) / (max_dist - ideal_dist)

                inter = Interaction(
                    interaction_type="Metal Coordination",
                    protein_residue=f"HOH{wa.res_seq}",
                    distance=d,
                    protein_atoms=[],
                    ligand_atoms=[metal],
                    confidence=max(0.5, dist_quality),
                    extra_info={
                        'metal': element,
                        'metal_id': metal_id,
                        'metal_atom_name': metal.name,
                        'metal_atom_serial': getattr(metal, 'serial', None),
                        'coord_atom': 'O',
                        'coord_type': 'water',
                        'source': 'water',
                        'water_atom': wa,
                        'ideal_dist': f"{ideal_dist:.2f}"
                    }
                )
                coordinations.append(inter)
                coord_positions.append((wa.x, wa.y, wa.z))

        # Determine coordination geometry
        cn = len(coordinations)
        if cn >= 2:
            geom_result = determine_coordination_geometry_rms(metal_pos, coord_positions)
            geometry = geom_result['geometry']
            geom_confidence = geom_result['confidence']
            rms_deviation = geom_result['rms']
            observed_angles = geom_result['observed_angles']
            tau_param = geom_result['tau_parameter']
            tau_type = geom_result['tau_type']
            tau_interp = geom_result['tau_interpretation']
            all_fits = geom_result['all_fits']
            next_best = geom_result['next_best']
        else:
            geometry = COORDINATION_GEOMETRIES.get(cn, f"coordination_{cn}")
            geom_confidence = 0.5
            rms_deviation = None
            observed_angles = []
            tau_param = None
            tau_type = None
            tau_interp = None
            all_fits = {}
            next_best = None

        for inter in coordinations:
            inter.extra_info['coordination_number'] = cn
            inter.extra_info['geometry'] = geometry
            inter.extra_info['geometry_confidence'] = f"{geom_confidence:.2f}"

            if rms_deviation is not None:
                inter.extra_info['geometry_rms'] = f"{rms_deviation:.1f}°"

            if tau_param is not None:
                inter.extra_info['tau_parameter'] = f"{tau_param:.3f}"
                inter.extra_info['tau_type'] = tau_type
                inter.extra_info['tau_interpretation'] = tau_interp

            if all_fits:
                inter.extra_info['alternative_geometries'] = all_fits

            if next_best:
                inter.extra_info['next_best_geometry'] = next_best['geometry']
                inter.extra_info['next_best_rms'] = f"{next_best['rms']:.1f}°"

            if observed_angles:
                inter.extra_info['observed_angles'] = observed_angles

            inter.confidence = (inter.confidence + geom_confidence) / 2

        if coordinations:
            results[metal_id] = coordinations

    return results
