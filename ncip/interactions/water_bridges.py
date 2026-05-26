"""
Water bridge and XWH (halogen-water-hydrogen) bridge detection.

Water bridges: Jiang et al. (2005)
XWH bridges: Zhou et al. (2010), Lu et al. (2010), Verteramo et al. (2024)
"""

from typing import List, Optional
from ..parser import Atom, Interaction, distance, angle
from ..constants import (
    WATER_BRIDGE_DIST_MAX, WATER_BRIDGE_ANGLE_MIN, WATER_BRIDGE_ANGLE_MAX,
    XWH_XBOND_ANGLE_MIN, XWH_HBOND_DIST_MAX, XWH_HBOND_ANGLE_MIN,
)
from .helpers import is_protein_donor, is_protein_acceptor, is_ligand_donor, is_ligand_acceptor


def find_water_bridges(protein_atoms: List[Atom],
                       ligand_atoms: List[Atom],
                       water_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect water-mediated hydrogen bonds between protein and ligand.

    Geometry requirements:
    - Protein-Water distance: ≤3.5Å
    - Water-Ligand distance: ≤3.5Å
    - Water angle (D-W-A): 80-140° (Jiang et al. 2005)

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms
        water_atoms: List of water oxygen atoms

    Returns:
        List of Interaction objects for detected water bridges
    """
    interactions = []

    if not water_atoms:
        return interactions

    WATER_DIST_MAX = WATER_BRIDGE_DIST_MAX
    WATER_ANGLE_MIN_VAL = WATER_BRIDGE_ANGLE_MIN
    WATER_ANGLE_MAX_VAL = WATER_BRIDGE_ANGLE_MAX

    prot_hbond = [a for a in protein_atoms
                  if is_protein_donor(a) or is_protein_acceptor(a)]
    lig_hbond = [a for a in ligand_atoms
                 if is_ligand_donor(a) or is_ligand_acceptor(a)]

    for water in water_atoms:
        near_prot = []
        for pa in prot_hbond:
            d = distance(water, pa)
            if d <= WATER_DIST_MAX:
                near_prot.append((pa, d))

        near_lig = []
        for la in lig_hbond:
            d = distance(water, la)
            if d <= WATER_DIST_MAX:
                near_lig.append((la, d))

        if near_prot and near_lig:
            for pa, d_prot in near_prot:
                for la, d_lig in near_lig:
                    water_angle = angle(pa, water, la)

                    if water_angle < WATER_ANGLE_MIN_VAL or water_angle > WATER_ANGLE_MAX_VAL:
                        continue

                    optimal_angle = 105.0
                    angle_deviation = abs(water_angle - optimal_angle)
                    angle_quality = max(0, 1.0 - angle_deviation / 35)

                    dist_quality = 1.0 - ((d_prot + d_lig) / 2) / WATER_DIST_MAX * 0.3

                    confidence = (angle_quality + dist_quality) / 2

                    inter = Interaction(
                        interaction_type="Water Bridge",
                        protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                        distance=(d_prot + d_lig) / 2,
                        protein_atoms=[pa, water],
                        ligand_atoms=[la],
                        confidence=confidence,
                        extra_info={
                            'water_id': water.serial,
                            'water_res_seq': water.res_seq,
                            'prot_water_dist': f"{d_prot:.2f}",
                            'lig_water_dist': f"{d_lig:.2f}",
                            'water_angle': f"{water_angle:.1f}°",
                            'geometry_quality': 'good' if confidence > 0.7 else 'moderate'
                        }
                    )
                    interactions.append(inter)

    return interactions


def find_xwh_bridges(protein_atoms: List[Atom],
                     ligand_atoms: List[Atom],
                     water_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect halogen-water-hydrogen (XWH) bridges between ligand and protein.

    XWH bridges occur when a ligand halogen (Cl, Br, I) forms a halogen bond
    to a water oxygen, and that same water simultaneously H-bonds to a protein
    donor/acceptor atom.

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms
        water_atoms: List of water oxygen atoms

    Returns:
        List of Interaction objects for detected XWH bridges
    """
    interactions = []

    if not water_atoms:
        return interactions

    xwh_halogen_elements = ('CL', 'BR', 'I')

    xwh_xbond_dist_cutoffs = {
        'CL': 3.27,
        'BR': 3.37,
        'I':  3.50,
    }

    _xwh_xbond_angle_min = XWH_XBOND_ANGLE_MIN
    _xwh_hbond_dist_max_no = XWH_HBOND_DIST_MAX
    _xwh_hbond_dist_max_s  = XWH_HBOND_DIST_MAX + 0.2
    _xwh_hbond_angle_min   = XWH_HBOND_ANGLE_MIN

    lig_halogens = [a for a in ligand_atoms
                    if a.element.upper() in xwh_halogen_elements]

    if not lig_halogens:
        return interactions

    lig_carbons = [a for a in ligand_atoms if a.element == 'C']

    prot_hbond = [a for a in protein_atoms
                  if a.element in ('O', 'N', 'S')]

    def find_attached_carbon(hal: Atom) -> Optional[Atom]:
        best_c = None
        best_d = 2.2
        for c in lig_carbons:
            d = distance(hal, c)
            if d < best_d:
                best_d = d
                best_c = c
        return best_c

    def find_bonded_heavy(atom: Atom, max_dist: float = 1.8) -> Optional[Atom]:
        best = None
        best_d = max_dist
        for other in protein_atoms:
            if other.serial != atom.serial and other.element != 'H':
                d = distance(atom, other)
                if d < best_d:
                    best_d = d
                    best = other
        return best

    for hal in lig_halogens:
        hal_element = hal.element.upper()
        xbond_cutoff = xwh_xbond_dist_cutoffs.get(hal_element, 3.50)

        attached_c = find_attached_carbon(hal)
        if not attached_c:
            continue

        for water in water_atoms:
            d1 = distance(hal, water)
            if d1 > xbond_cutoff:
                continue

            theta = angle(attached_c, hal, water)
            if theta < _xwh_xbond_angle_min:
                continue

            for pa in prot_hbond:
                if pa.element == 'S':
                    hbond_cutoff = _xwh_hbond_dist_max_s
                else:
                    hbond_cutoff = _xwh_hbond_dist_max_no

                d2 = distance(water, pa)
                if d2 > hbond_cutoff:
                    continue

                bonded = find_bonded_heavy(pa)
                if bonded:
                    phi = angle(water, pa, bonded)
                    if phi < _xwh_hbond_angle_min:
                        continue
                else:
                    phi = None

                xbond_penetration = 1.0 - (d1 / xbond_cutoff)
                xbond_angle_quality = (theta - _xwh_xbond_angle_min) / (180.0 - _xwh_xbond_angle_min)

                hbond_dist_quality = 1.0 - (d2 / hbond_cutoff) * 0.3
                hbond_angle_quality = ((phi - XWH_HBOND_ANGLE_MIN) / (180.0 - XWH_HBOND_ANGLE_MIN)) if phi else 0.5

                confidence = (xbond_penetration + xbond_angle_quality + hbond_dist_quality + hbond_angle_quality) / 4
                confidence = max(0.1, min(1.0, confidence))

                pct_shortening = (1.0 - d1 / xbond_cutoff) * 100

                inter = Interaction(
                    interaction_type="XWH Bridge",
                    protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                    distance=(d1 + d2) / 2,
                    protein_atoms=[pa, water],
                    ligand_atoms=[hal],
                    confidence=confidence,
                    extra_info={
                        'halogen': hal_element,
                        'water_id': water.serial,
                        'water_res_seq': water.res_seq,
                        'xbond_dist': f"{d1:.2f}",
                        'xbond_angle': f"{theta:.1f}°",
                        'hbond_dist': f"{d2:.2f}",
                        'hbond_angle': f"{phi:.1f}°" if phi else "N/A",
                        'pct_shortening': f"{pct_shortening:.1f}%",
                        'protein_atom': f"{pa.res_name}{pa.res_seq} {pa.name}",
                        'carbon_atom': attached_c.name,
                        'geometry_quality': 'good' if confidence > 0.5 else 'moderate'
                    }
                )
                interactions.append(inter)

    return interactions
