"""
Hydrophobic contact detection.

Enhanced algorithm that properly validates hydrophobic carbons by
excluding polar-adjacent carbons and classifying aliphatic vs aromatic.
"""

from typing import List, Tuple
from ..parser import Atom, Interaction, distance
from ..constants import HYDROPHOBIC_RESIDUES, HYDROPHOBIC_DIST_MIN, HYDROPHOBIC_DIST_MAX


def is_hydrophobic_carbon(carbon: Atom, all_atoms: List[Atom]) -> Tuple[bool, str]:
    """
    Determine if a carbon atom is truly hydrophobic.

    A carbon is NOT hydrophobic if:
    - It is directly bonded to N, O, S (polar atoms)
    - It is a carbonyl carbon (C=O)
    - It is part of a carboxyl group (COOH/COO-)

    Args:
        carbon: The carbon atom to check
        all_atoms: All atoms in the molecule

    Returns:
        Tuple of (is_hydrophobic, classification)
        classification: 'aliphatic', 'aromatic', or 'polar'
    """
    if carbon.element != 'C':
        return False, 'not_carbon'

    polar_elements = ('N', 'O', 'S', 'F')
    bonded_polar = []
    bonded_carbon = []
    bonded_hydrogen = []

    for other in all_atoms:
        if other.serial == carbon.serial:
            continue
        d = distance(carbon, other)

        if d < 1.8:
            if other.element in polar_elements:
                bonded_polar.append((other, d))
            elif other.element == 'C':
                bonded_carbon.append((other, d))
            elif other.element == 'H':
                bonded_hydrogen.append((other, d))

    if bonded_polar:
        for atom, d in bonded_polar:
            if atom.element == 'O' and d < 1.3:
                return False, 'carbonyl'
        only_f_or_s = all(a.element in ('F', 'S') for a, _ in bonded_polar)
        num_ch_neighbors = len(bonded_carbon) + len(bonded_hydrogen)
        if only_f_or_s and num_ch_neighbors >= 2:
            pass
        else:
            return False, 'polar'

    num_c_neighbors = len(bonded_carbon)
    num_h_neighbors = len(bonded_hydrogen)
    total_bonds = num_c_neighbors + num_h_neighbors

    if total_bonds == 3 and num_c_neighbors >= 2:
        return True, 'aromatic'
    else:
        return True, 'aliphatic'


def find_hydrophobic_contacts(protein_atoms: List[Atom],
                               ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect hydrophobic contacts between protein and ligand.

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected hydrophobic contacts
    """
    interactions = []

    hydrophobic_prot = []
    for a in protein_atoms:
        if a.res_name in HYDROPHOBIC_RESIDUES:
            if a.name in HYDROPHOBIC_RESIDUES[a.res_name]:
                hydrophobic_prot.append(a)

    hydrophobic_lig = []
    for a in ligand_atoms:
        if a.element == 'C':
            is_hydrophobic, classification = is_hydrophobic_carbon(a, ligand_atoms)
            if is_hydrophobic:
                hydrophobic_lig.append((a, classification))

    seen = set()
    for pa in hydrophobic_prot:
        for la, la_type in hydrophobic_lig:
            d = distance(pa, la)
            if HYDROPHOBIC_DIST_MIN <= d <= HYDROPHOBIC_DIST_MAX:
                key = (pa.chain_id, pa.res_name, pa.res_seq, la.serial)
                if key not in seen:
                    seen.add(key)

                    optimal_dist = 3.8
                    dist_deviation = abs(d - optimal_dist)
                    confidence = max(0.3, 1.0 - dist_deviation / 2.0)

                    inter = Interaction(
                        interaction_type="Hydrophobic",
                        protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                        distance=d,
                        protein_atoms=[pa],
                        ligand_atoms=[la],
                        confidence=confidence,
                        extra_info={
                            'ligand_carbon_type': la_type,
                            'protein_atom': pa.name
                        }
                    )
                    interactions.append(inter)

    return interactions
