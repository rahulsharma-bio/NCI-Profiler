"""
Hydrogen bond detection.

Enhanced algorithm with strength classification based on Jeffrey (1997).
"""

from typing import List, Tuple, Optional
from ..parser import Atom, Interaction, distance, angle
from ..constants import HBOND_DISTANCE_MAX, HBOND_DISTANCE_MIN, HBOND_ANGLE_MIN
from .helpers import is_protein_donor, is_protein_acceptor, is_ligand_donor, is_ligand_acceptor


def classify_hbond_strength(d_a_dist: float, h_a_dist: Optional[float],
                            dha_angle: Optional[float]) -> Tuple[str, float]:
    """
    Classify hydrogen bond strength based on geometry.

    Based on Jeffrey (1997) "An Introduction to Hydrogen Bonding":

    Strong H-bonds:
    - D-A distance: 2.2-2.5 Å
    - H-A distance: 1.2-1.5 Å
    - D-H-A angle: 170-180°

    Moderate H-bonds:
    - D-A distance: 2.5-3.2 Å
    - H-A distance: 1.5-2.2 Å
    - D-H-A angle: 130-170°

    Weak H-bonds:
    - D-A distance: 3.2-4.0 Å
    - H-A distance: 2.2-3.2 Å
    - D-H-A angle: 90-130°

    Returns:
        Tuple of (strength_category, confidence_score)
    """
    if h_a_dist is None:
        h_a_dist = d_a_dist - 1.0
    if dha_angle is None:
        dha_angle = 150.0

    dist_score = 0.0
    angle_score = 0.0
    strength = "weak"

    if d_a_dist <= 2.5:
        dist_score = 1.0
    elif d_a_dist <= 3.2:
        dist_score = 0.7 - (d_a_dist - 2.5) * 0.3
    else:
        dist_score = max(0, 0.4 - (d_a_dist - 3.2) * 0.2)

    if dha_angle >= 170:
        angle_score = 1.0
    elif dha_angle >= 130:
        angle_score = 0.6 + (dha_angle - 130) * 0.01
    else:
        angle_score = max(0, (dha_angle - 90) * 0.015)

    combined_score = (dist_score + angle_score) / 2

    if d_a_dist <= 2.5 and dha_angle >= 170:
        strength = "strong"
    elif d_a_dist <= 3.2 and dha_angle >= 130:
        strength = "moderate"
    else:
        strength = "weak"

    return strength, combined_score


def find_hydrogen_bonds(protein_atoms: List[Atom],
                        ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect hydrogen bonds between protein and ligand.

    Enhanced algorithm with strength classification based on Jeffrey (1997).
    Both protein-donor/ligand-acceptor and ligand-donor/protein-acceptor
    interactions are detected.

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected H-bonds
    """
    interactions = []
    seen = set()

    prot_hydrogens = [a for a in protein_atoms if a.element == 'H']
    lig_hydrogens = [a for a in ligand_atoms if a.element == 'H']

    def find_bonded_hydrogen(donor, hydrogen_list, max_dist=1.3):
        """Find hydrogen atom bonded to a donor atom."""
        best_h = None
        best_dist = max_dist
        for h in hydrogen_list:
            d = distance(donor, h)
            if d < best_dist:
                best_dist = d
                best_h = h
        return best_h

    prot_donors = [a for a in protein_atoms if is_protein_donor(a)]
    prot_acceptors = [a for a in protein_atoms if is_protein_acceptor(a)]
    lig_donors = [a for a in ligand_atoms if is_ligand_donor(a)]
    lig_acceptors = [a for a in ligand_atoms if is_ligand_acceptor(a)]

    # Protein donor -> Ligand acceptor
    for pd in prot_donors:
        for la in lig_acceptors:
            d = distance(pd, la)
            if HBOND_DISTANCE_MIN <= d <= HBOND_DISTANCE_MAX:
                key = (pd.serial, la.serial)
                if key not in seen:
                    seen.add(key)

                    extra = {
                        'role': f"Protein donor ({pd.name}) -> Ligand acceptor ({la.name})",
                        'D-A': f"{d:.2f}"
                    }

                    h_a_dist = None
                    dha_angle = None

                    h_atom = find_bonded_hydrogen(pd, prot_hydrogens)
                    if h_atom:
                        h_a_dist = distance(h_atom, la)
                        extra['H-A'] = f"{h_a_dist:.2f}"
                        dha_angle = angle(pd, h_atom, la)
                        extra['angle'] = f"{dha_angle:.1f}°"

                    strength, confidence = classify_hbond_strength(d, h_a_dist, dha_angle)
                    extra['strength'] = strength

                    inter = Interaction(
                        interaction_type="Hydrogen Bond",
                        protein_residue=f"{pd.chain_id}:{pd.res_name}{pd.res_seq}",
                        distance=d,
                        protein_atoms=[pd],
                        ligand_atoms=[la],
                        confidence=confidence,
                        extra_info=extra
                    )
                    interactions.append(inter)

    # Ligand donor -> Protein acceptor
    for ld in lig_donors:
        for pa in prot_acceptors:
            d = distance(pa, ld)
            if HBOND_DISTANCE_MIN <= d <= HBOND_DISTANCE_MAX:
                key = (pa.serial, ld.serial)
                if key not in seen:
                    seen.add(key)

                    extra = {
                        'role': f"Ligand donor ({ld.name}) -> Protein acceptor ({pa.name})",
                        'D-A': f"{d:.2f}"
                    }

                    h_a_dist = None
                    dha_angle = None

                    h_atom = find_bonded_hydrogen(ld, lig_hydrogens)
                    if h_atom:
                        h_a_dist = distance(h_atom, pa)
                        extra['H-A'] = f"{h_a_dist:.2f}"
                        dha_angle = angle(ld, h_atom, pa)
                        extra['angle'] = f"{dha_angle:.1f}°"

                    strength, confidence = classify_hbond_strength(d, h_a_dist, dha_angle)
                    extra['strength'] = strength

                    inter = Interaction(
                        interaction_type="Hydrogen Bond",
                        protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                        distance=d,
                        protein_atoms=[pa],
                        ligand_atoms=[ld],
                        confidence=confidence,
                        extra_info=extra
                    )
                    interactions.append(inter)

    return interactions
