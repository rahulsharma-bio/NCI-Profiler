"""
Halogen bond detection.

Enhanced algorithm based on Auffinger et al. (2004) and Jiang et al. (2005).
Fluorine is excluded — only Cl, Br, I form sigma holes.
"""

from typing import List, Optional
from ..parser import Atom, Interaction, distance, angle
from ..constants import HALOGEN_BOND_DIST_MAX


def find_halogen_bonds(protein_atoms: List[Atom],
                       ligand_atoms: List[Atom]) -> List[Interaction]:
    """
    Detect halogen bonds between ligand halogens and protein acceptors.

    Geometry requirements:
    - C-X...A angle (donor angle): 140-180° (sigma-hole alignment)
    - X...A-Y angle (acceptor angle): 90-140° (lone pair alignment)
    - Distance: varies by halogen (Cl: 3.3Å, Br: 3.4Å, I: 3.5Å)

    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms

    Returns:
        List of Interaction objects for detected halogen bonds
    """
    interactions = []

    halogen_elements = ('CL', 'BR', 'I')

    lig_halogens = [a for a in ligand_atoms
                   if a.element.upper() in halogen_elements]

    if not lig_halogens:
        return interactions

    halogen_dist_cutoffs = {
        'CL': 3.3,
        'BR': 3.4,
        'I': 3.5,
    }

    prot_acceptors = [a for a in protein_atoms
                      if a.element in ('O', 'N', 'S')]

    def find_bonded_heavy(atom: Atom, atoms: List[Atom], max_dist: float = 1.8) -> List[Atom]:
        """Find heavy atoms bonded to the given atom."""
        bonded = []
        for other in atoms:
            if other.serial != atom.serial and other.element != 'H':
                d = distance(atom, other)
                if d < max_dist:
                    bonded.append(other)
        return bonded

    lig_carbons = [a for a in ligand_atoms if a.element == 'C']

    for hal in lig_halogens:
        hal_element = hal.element.upper()
        dist_cutoff = halogen_dist_cutoffs.get(hal_element, HALOGEN_BOND_DIST_MAX)

        attached_c = None
        min_dist = 2.2
        for c in lig_carbons:
            d = distance(hal, c)
            if d < min_dist:
                min_dist = d
                attached_c = c

        if not attached_c:
            continue

        for pa in prot_acceptors:
            d = distance(hal, pa)
            if d > dist_cutoff:
                continue

            donor_angle = angle(attached_c, hal, pa)

            if donor_angle < 140:
                continue

            bonded_to_acceptor = find_bonded_heavy(pa, protein_atoms)

            acceptor_angle = None
            valid_acceptor_geometry = True

            if bonded_to_acceptor:
                y_atom = bonded_to_acceptor[0]
                acceptor_angle = angle(hal, pa, y_atom)
                if not (90 <= acceptor_angle <= 140):
                    valid_acceptor_geometry = False

            if not valid_acceptor_geometry:
                continue

            donor_deviation = abs(170 - donor_angle)
            dist_quality = 1.0 - (d / dist_cutoff) * 0.3
            angle_quality = max(0, 1.0 - donor_deviation / 30)

            confidence = (dist_quality + angle_quality) / 2

            if acceptor_angle:
                acceptor_deviation = abs(120 - acceptor_angle)
                acceptor_quality = max(0, 1.0 - acceptor_deviation / 30)
                confidence = (confidence + acceptor_quality) / 2

            inter = Interaction(
                interaction_type="Halogen Bond",
                protein_residue=f"{pa.chain_id}:{pa.res_name}{pa.res_seq}",
                distance=d,
                protein_atoms=[pa],
                ligand_atoms=[hal],
                confidence=confidence,
                extra_info={
                    'type': f"{hal.element}...{pa.element}",
                    'donor_angle': f"{donor_angle:.1f}°",
                    'acceptor_angle': f"{acceptor_angle:.1f}°" if acceptor_angle else "N/A",
                    'halogen': hal_element,
                    'geometry_quality': 'good' if confidence > 0.7 else 'moderate'
                }
            )
            interactions.append(inter)

    return interactions
