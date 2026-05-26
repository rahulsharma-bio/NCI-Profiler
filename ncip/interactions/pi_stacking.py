"""
π-π stacking detection.

Enhanced algorithm based on:
- Calinsky & Levy (2024) - 4-parameter geometric model
- Li et al. (2025) - Fused ring handling
- McGaughey et al. (1998) - Original geometry definitions
"""

import math
from typing import List, Tuple
from collections import defaultdict

from ..parser import (
    Atom, Ring, Interaction,
    distance_coords, angle_between_normals,
    detect_rings, detect_ligand_rings
)
from ..constants import PI_STACK_DIST_MAX


def calculate_elevation_angle(ring1_centroid: Tuple[float, float, float],
                              ring2_centroid: Tuple[float, float, float],
                              ring2_normal: Tuple[float, float, float]) -> float:
    """
    Calculate elevation angle (Tθ) of ring1 centroid relative to ring2 plane.

    Per Calinsky & Levy (2024):
    - Tθ = arcsin(vertical_distance / centroid_distance)
    - For true π-stacking: Tθ should be > 65°

    Returns:
        Elevation angle in degrees (0-90°)
    """
    v = (ring1_centroid[0] - ring2_centroid[0],
         ring1_centroid[1] - ring2_centroid[1],
         ring1_centroid[2] - ring2_centroid[2])

    d = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    if d < 0.001:
        return 90.0

    vert_dist = abs(v[0]*ring2_normal[0] + v[1]*ring2_normal[1] + v[2]*ring2_normal[2])
    sin_theta = min(1.0, vert_dist / d)
    theta = math.degrees(math.asin(sin_theta))

    return theta


def classify_pi_stacking_enhanced(ring_angle: float, vert_dist: float,
                                  horiz_offset: float, elevation_theta1: float,
                                  elevation_theta2: float) -> Tuple[str, float, str]:
    """
    Enhanced π-π stacking classification using 4-parameter model.

    Based on Calinsky & Levy (2024) J. Phys. Chem. B.

    Returns:
        Tuple of (stacking_type, confidence, classification_note)
    """
    min_elevation = min(elevation_theta1, elevation_theta2)
    max_elevation = max(elevation_theta1, elevation_theta2)

    confidence = 0.5
    note = ""

    if ring_angle >= 60:
        stack_type = "T-shaped"
        note = "CH-π dominant (edge-to-face)"
        angle_quality = 1.0 - abs(90 - ring_angle) / 30
        elevation_quality = 1.0 - min_elevation / 90
        confidence = (angle_quality * 0.6 + elevation_quality * 0.4)

    elif ring_angle <= 30:
        if min_elevation >= 65:
            if horiz_offset <= 2.0:
                stack_type = "π-stacked"
                note = "face-centered stacking"
            else:
                stack_type = "offset-π"
                note = "parallel-displaced stacking"

            angle_quality = 1.0 - ring_angle / 30
            elevation_quality = (min_elevation - 65) / 25
            dist_quality = 1.0 - abs(vert_dist - 3.5) / 1.5
            confidence = (angle_quality + elevation_quality + dist_quality) / 3

        elif min_elevation >= 45:
            stack_type = "offset-parallel"
            note = "tilted parallel (partial overlap)"
            angle_quality = 1.0 - ring_angle / 30
            elevation_quality = (min_elevation - 45) / 20
            confidence = (angle_quality + elevation_quality) / 2 * 0.8

        else:
            stack_type = "edge-parallel"
            note = "edge approach despite parallel planes"
            confidence = 0.35

    else:
        stack_type = "intermediate"
        note = "tilted geometry"
        confidence = 0.4

    return stack_type, max(0.2, min(1.0, confidence)), note


def select_closest_rings_per_residue(protein_rings: List[Ring],
                                     ligand_rings: List[Ring]) -> List[Tuple[Ring, Ring, float]]:
    """
    For fused ring systems, select only the closest ring per protein residue.

    Per Li et al. (2025): prevents over-counting for porphyrins
    and other polycyclic aromatics.

    Returns:
        List of (protein_ring, ligand_ring, distance) tuples
    """
    pairs = []

    ligand_systems = defaultdict(list)
    for lring in ligand_rings:
        if lring.ring_system_id is not None:
            ligand_systems[lring.ring_system_id].append(lring)
        else:
            ligand_systems[id(lring)].append(lring)

    protein_by_residue = defaultdict(list)
    for pring in protein_rings:
        res_key = (pring.chain_id, pring.res_seq, pring.res_name)
        protein_by_residue[res_key].append(pring)

    for res_key, prot_rings_in_res in protein_by_residue.items():
        for system_id, lig_rings_in_system in ligand_systems.items():

            best_pair = None
            best_dist = float('inf')

            for pring in prot_rings_in_res:
                for lring in lig_rings_in_system:
                    d = distance_coords(pring.centroid, lring.centroid)
                    if d < best_dist:
                        best_dist = d
                        best_pair = (pring, lring)

            if best_pair is not None and best_dist <= PI_STACK_DIST_MAX:
                pairs.append((best_pair[0], best_pair[1], best_dist))

    return pairs


def find_pi_stacking(protein_atoms: List[Atom],
                     ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect π-π stacking interactions between aromatic rings.

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected π-π interactions
    """
    interactions = []

    prot_by_res = defaultdict(list)
    for a in protein_atoms:
        key = (a.chain_id, a.res_name, a.res_seq)
        prot_by_res[key].append(a)

    protein_rings = []
    for (chain, res_name, res_seq), atoms in prot_by_res.items():
        rings = detect_rings(atoms, res_name)
        protein_rings.extend(rings)

    ligand_rings = detect_ligand_rings(ligand_atoms)

    ring_pairs = select_closest_rings_per_residue(protein_rings, ligand_rings)

    for pring, lring, d in ring_pairs:
        ring_angle = angle_between_normals(pring.normal, lring.normal)

        v = (lring.centroid[0] - pring.centroid[0],
             lring.centroid[1] - pring.centroid[1],
             lring.centroid[2] - pring.centroid[2])

        vert_dist = abs(v[0]*pring.normal[0] +
                       v[1]*pring.normal[1] +
                       v[2]*pring.normal[2])

        horiz_offset = (d**2 - vert_dist**2)**0.5 if d > vert_dist else 0

        theta1 = calculate_elevation_angle(pring.centroid, lring.centroid, lring.normal)
        theta2 = calculate_elevation_angle(lring.centroid, pring.centroid, pring.normal)

        stack_type, confidence, note = classify_pi_stacking_enhanced(
            ring_angle, vert_dist, horiz_offset, theta1, theta2
        )

        valid = True

        if stack_type in ("π-stacked", "offset-π"):
            if vert_dist < 3.0 or vert_dist > 4.2:
                valid = False
            if min(theta1, theta2) < 55:
                valid = False

        elif stack_type == "offset-parallel":
            if vert_dist < 2.8 or vert_dist > 4.5:
                valid = False

        elif stack_type == "T-shaped":
            if d < 4.0 or d > 6.5:
                valid = False
            if min(theta1, theta2) > 45:
                valid = False

        elif stack_type == "edge-parallel":
            confidence *= 0.7

        elif stack_type == "intermediate":
            if min(theta1, theta2) < 40:
                valid = False

        if not valid:
            continue

        inter = Interaction(
            interaction_type="π-π Stacking",
            protein_residue=pring.residue_id,
            distance=d,
            protein_atoms=pring.atoms,
            ligand_atoms=lring.atoms,
            confidence=confidence,
            extra_info={
                'type': stack_type,
                'note': note,
                'angle': f"{ring_angle:.1f}°",
                'elevation_theta1': f"{theta1:.1f}°",
                'elevation_theta2': f"{theta2:.1f}°",
                'vertical_dist': f"{vert_dist:.2f}",
                'horizontal_offset': f"{horiz_offset:.2f}",
                'prot_ring_size': getattr(pring, 'ring_size', 'N/A'),
                'lig_ring_size': getattr(lring, 'ring_size', 'N/A'),
                'geometry_quality': 'good' if confidence > 0.6 else 'moderate' if confidence > 0.4 else 'weak'
            }
        )
        interactions.append(inter)

    return interactions
