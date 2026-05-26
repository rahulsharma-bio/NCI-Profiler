"""
Cation-π and π-cation interaction detection.

Detects both directions:
- Cation-π: Ligand cation -> Protein aromatic ring
- π-Cation: Protein cation -> Ligand aromatic ring
"""

from typing import List
from collections import defaultdict

from ..parser import (
    Atom, Ring, Interaction,
    distance_coords, centroid,
    detect_rings, detect_ligand_rings, detect_charged_groups
)
from ..constants import PROTEIN_CATIONS, CATION_PI_DIST_MAX


def find_cation_pi(protein_atoms: List[Atom],
                   ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect cation-π and π-cation interactions between protein and ligand.

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected cation-π and π-cation interactions
    """
    interactions = []

    prot_by_res = defaultdict(list)
    for a in protein_atoms:
        key = (a.chain_id, a.res_name, a.res_seq)
        prot_by_res[key].append(a)

    protein_cations = []
    for (chain, res_name, res_seq), atoms in prot_by_res.items():
        if res_name in PROTEIN_CATIONS:
            cat_atoms = [a for a in atoms if a.name in PROTEIN_CATIONS[res_name]]
            if cat_atoms:
                protein_cations.append({
                    'atoms': cat_atoms,
                    'centroid': centroid(cat_atoms),
                    'chain': chain,
                    'res_name': res_name,
                    'res_seq': res_seq
                })

    protein_rings = []
    for (chain, res_name, res_seq), atoms in prot_by_res.items():
        rings = detect_rings(atoms, res_name)
        protein_rings.extend(rings)

    ligand_rings = detect_ligand_rings(ligand_atoms)

    lig_positive, _ = detect_charged_groups(ligand_atoms)

    def calc_cation_pi_confidence(d: float, ring) -> float:
        """Calculate confidence based on distance and alignment."""
        optimal_dist = 4.5
        dist_quality = max(0, 1.0 - abs(d - optimal_dist) / 2.5)
        if d < 3.0 or d > 7.0:
            dist_quality *= 0.5
        return max(0.3, dist_quality)

    # Protein cation -> Ligand aromatic = π-Cation
    for cat in protein_cations:
        for lring in ligand_rings:
            d = distance_coords(cat['centroid'], lring.centroid)
            if d <= CATION_PI_DIST_MAX:
                confidence = calc_cation_pi_confidence(d, lring)

                inter = Interaction(
                    interaction_type="π-Cation",
                    protein_residue=f"{cat['chain']}:{cat['res_name']}{cat['res_seq']}",
                    distance=d,
                    protein_atoms=cat['atoms'],
                    ligand_atoms=lring.atoms,
                    confidence=confidence,
                    extra_info={
                        'type': f"Ligand π -> Protein cation ({cat['res_name']})",
                        'cation_type': cat['res_name'],
                        'geometry_quality': 'good' if confidence > 0.6 else 'moderate'
                    }
                )
                interactions.append(inter)

    # Ligand cation -> Protein aromatic = Cation-π
    for lig_group in lig_positive:
        for pring in protein_rings:
            d = distance_coords(lig_group.charge_center, pring.centroid)
            if d <= CATION_PI_DIST_MAX:
                confidence = calc_cation_pi_confidence(d, pring) * lig_group.confidence

                inter = Interaction(
                    interaction_type="Cation-π",
                    protein_residue=pring.residue_id,
                    distance=d,
                    protein_atoms=pring.atoms,
                    ligand_atoms=lig_group.atoms,
                    confidence=confidence,
                    extra_info={
                        'type': f"Ligand cation ({lig_group.group_type}) -> Protein π ({pring.res_name})",
                        'cation_type': lig_group.group_type,
                        'geometry_quality': 'good' if confidence > 0.6 else 'moderate'
                    }
                )
                interactions.append(inter)

    return interactions
