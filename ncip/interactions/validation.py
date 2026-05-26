"""
Post-detection interaction validation and filtering.
"""

from typing import List, Dict
from collections import defaultdict
from ..parser import Interaction


def validate_interactions(interactions: List[Interaction],
                          min_confidence: float = 0.3,
                          remove_duplicates: bool = True,
                          check_physical_constraints: bool = True) -> List[Interaction]:
    """
    Post-detection validation and filtering of interactions.

    Validation criteria:
    1. Confidence threshold
    2. Physical constraints (impossible geometries)
    3. Duplicate removal
    4. Mutual exclusion (same atoms, different types)
    """
    if not interactions:
        return []

    validated = []
    seen_keys = set()

    for inter in interactions:
        if inter.confidence < min_confidence:
            continue

        if check_physical_constraints:
            valid = True

            if inter.distance < 0.5:
                valid = False

            if inter.interaction_type == "Hydrogen Bond":
                if inter.protein_atoms:
                    pa = inter.protein_atoms[0]
                    if pa.element not in ('N', 'O', 'S', 'F'):
                        valid = False
                if inter.ligand_atoms:
                    la = inter.ligand_atoms[0]
                    if la.element not in ('N', 'O', 'S', 'F', 'H'):
                        valid = False

            elif inter.interaction_type == "Halogen Bond":
                if inter.ligand_atoms:
                    la = inter.ligand_atoms[0]
                    if la.element.upper() not in ('CL', 'BR', 'I'):
                        valid = False

            elif inter.interaction_type == "Hydrophobic":
                if inter.ligand_atoms:
                    la = inter.ligand_atoms[0]
                    if la.element != 'C':
                        valid = False

            if not valid:
                continue

        if remove_duplicates:
            prot_atoms = tuple(sorted(a.serial for a in inter.protein_atoms))
            lig_atoms = tuple(sorted(a.serial for a in inter.ligand_atoms))
            key = (inter.interaction_type, prot_atoms, lig_atoms)

            if key in seen_keys:
                continue
            seen_keys.add(key)

        validated.append(inter)

    # Mutual exclusion checks
    final = []
    validated_by_atoms = defaultdict(list)

    for inter in validated:
        prot_atoms = tuple(sorted(a.serial for a in inter.protein_atoms))
        lig_atoms = tuple(sorted(a.serial for a in inter.ligand_atoms))
        key = (prot_atoms, lig_atoms)
        validated_by_atoms[key].append(inter)

    for key, inters in validated_by_atoms.items():
        if len(inters) == 1:
            final.extend(inters)
        else:
            priority = {
                'Salt Bridge': 4,
                'Metal Coordination': 5,
                'Halogen Bond': 3,
                'Hydrogen Bond': 2,
                'π-π Stacking': 2,
                'Cation-π': 2,
                'Water Bridge': 1,
                'Hydrophobic': 0
            }

            inters_sorted = sorted(inters, key=lambda x: priority.get(x.interaction_type, 0), reverse=True)
            top_priority = priority.get(inters_sorted[0].interaction_type, 0)
            for inter in inters_sorted:
                if priority.get(inter.interaction_type, 0) == top_priority:
                    final.append(inter)
                    break

    return final


def get_interaction_summary(interactions: List[Interaction]) -> Dict:
    """Generate a summary of detected interactions."""
    summary = {
        'total': len(interactions),
        'by_type': defaultdict(int),
        'by_residue': defaultdict(list),
        'avg_confidence': 0.0,
        'high_confidence_count': 0,
        'medium_confidence_count': 0,
        'low_confidence_count': 0
    }

    if not interactions:
        return summary

    confidences = []
    for inter in interactions:
        summary['by_type'][inter.interaction_type] += 1
        summary['by_residue'][inter.protein_residue].append(inter.interaction_type)
        confidences.append(inter.confidence)

        if inter.confidence >= 0.7:
            summary['high_confidence_count'] += 1
        elif inter.confidence >= 0.4:
            summary['medium_confidence_count'] += 1
        else:
            summary['low_confidence_count'] += 1

    summary['avg_confidence'] = sum(confidences) / len(confidences)
    summary['by_type'] = dict(summary['by_type'])
    summary['by_residue'] = dict(summary['by_residue'])

    return summary
