"""
Salt bridge (ionic interaction) detection.

Enhanced algorithm using proper charged group detection with charge center
distances (Barlow & Thornton 1983).
"""

from typing import List
from collections import defaultdict

from ..parser import (
    Atom, Interaction, distance_coords, centroid,
    detect_charged_groups, is_valid_carboxylate
)
from ..constants import SALT_BRIDGE_DIST_MAX


def find_salt_bridges(protein_atoms: List[Atom],
                      ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect salt bridges (ionic interactions) between protein and ligand.

    Detects:
    - Protein positive (ARG guanidinium, LYS ε-amino, HIS imidazolium)
      with ligand negative (carboxylate, sulfonate, phosphate)
    - Protein negative (ASP/GLU carboxylate)
      with ligand positive (amines, guanidinium, amidinium)

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected salt bridges
    """
    interactions = []

    lig_positive, lig_negative = detect_charged_groups(ligand_atoms)

    prot_by_res = defaultdict(list)
    for a in protein_atoms:
        key = (a.chain_id, a.res_name, a.res_seq)
        prot_by_res[key].append(a)

    prot_positive = []
    prot_negative = []

    for (chain, res_name, res_seq), atoms in prot_by_res.items():
        residue_id = f"{chain}:{res_name}{res_seq}"

        if res_name == 'ARG':
            guanidinium = [a for a in atoms if a.name in ('CZ', 'NH1', 'NH2', 'NE')]
            if len(guanidinium) >= 3:
                center = centroid(guanidinium)
                prot_positive.append((residue_id, center, guanidinium, 0.95, 'guanidinium'))

        elif res_name == 'LYS':
            nz = [a for a in atoms if a.name == 'NZ']
            if nz:
                prot_positive.append((residue_id, (nz[0].x, nz[0].y, nz[0].z), nz, 0.90, 'primary_amine'))

        elif res_name == 'HIS':
            his_n = [a for a in atoms if a.name in ('ND1', 'NE2')]
            if len(his_n) == 2:
                center = centroid(his_n)
                prot_positive.append((residue_id, center, his_n, 0.50, 'imidazolium'))

        elif res_name == 'ASP':
            cg = [a for a in atoms if a.name == 'CG']
            od = [a for a in atoms if a.name in ('OD1', 'OD2')]
            if cg and len(od) == 2:
                is_carboxylate, _ = is_valid_carboxylate(cg[0], atoms)
                if is_carboxylate:
                    center = centroid(od)
                    prot_negative.append((residue_id, center, od, 0.95, 'carboxylate'))
                else:
                    center = centroid(od)
                    prot_negative.append((residue_id, center, od, 0.60, 'carboxylate'))

        elif res_name == 'GLU':
            cd = [a for a in atoms if a.name == 'CD']
            oe = [a for a in atoms if a.name in ('OE1', 'OE2')]
            if cd and len(oe) == 2:
                is_carboxylate, _ = is_valid_carboxylate(cd[0], atoms)
                if is_carboxylate:
                    center = centroid(oe)
                    prot_negative.append((residue_id, center, oe, 0.95, 'carboxylate'))
                else:
                    center = centroid(oe)
                    prot_negative.append((residue_id, center, oe, 0.60, 'carboxylate'))

    charge_center_dist_max = SALT_BRIDGE_DIST_MAX

    # Protein positive -> Ligand negative
    for res_id, prot_center, prot_atoms, prot_conf, prot_type in prot_positive:
        for lig_group in lig_negative:
            d = distance_coords(prot_center, lig_group.charge_center)
            if d <= charge_center_dist_max:
                combined_conf = prot_conf * lig_group.confidence

                inter = Interaction(
                    interaction_type="Salt Bridge",
                    protein_residue=res_id,
                    distance=d,
                    protein_atoms=prot_atoms,
                    ligand_atoms=lig_group.atoms,
                    confidence=combined_conf,
                    extra_info={
                        'type': f"Protein+ ({prot_type}) - Ligand- ({lig_group.group_type})",
                        'prot_group': prot_type,
                        'lig_group': lig_group.group_type,
                        'charge_center_dist': f"{d:.2f}"
                    }
                )
                interactions.append(inter)

    # Protein negative -> Ligand positive
    for res_id, prot_center, prot_atoms, prot_conf, prot_type in prot_negative:
        for lig_group in lig_positive:
            d = distance_coords(prot_center, lig_group.charge_center)
            if d <= charge_center_dist_max:
                combined_conf = prot_conf * lig_group.confidence

                inter = Interaction(
                    interaction_type="Salt Bridge",
                    protein_residue=res_id,
                    distance=d,
                    protein_atoms=prot_atoms,
                    ligand_atoms=lig_group.atoms,
                    confidence=combined_conf,
                    extra_info={
                        'type': f"Protein- ({prot_type}) - Ligand+ ({lig_group.group_type})",
                        'prot_group': prot_type,
                        'lig_group': lig_group.group_type,
                        'charge_center_dist': f"{d:.2f}"
                    }
                )
                interactions.append(inter)

    return interactions
