"""
Core analysis functions for protein-ligand interaction analysis.

This module contains the main orchestration functions for running
complete interaction analyses.
"""

import os
import math
from datetime import datetime
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict

from .parser import (
    Atom, Ring, Interaction,
    parse_pdb
)
from .interactions import (
    find_hydrogen_bonds, find_salt_bridges, find_hydrophobic_contacts,
    find_pi_stacking, find_cation_pi, find_halogen_bonds,
    find_water_bridges, find_xwh_bridges, detect_metal_coordination
)
from .filters import apply_filters_to_results


# ============================================================================
# DISTANCE FILTER FOR MULTI-LIGAND SCENARIOS
# ============================================================================

# Maximum distance (Å) from ligand centroid for valid interactions
MAX_INTERACTION_DISTANCE = 12.0  # Generous but catches ~50Å spurious interactions

def _calculate_centroid(atoms: List[Atom]) -> Tuple[float, float, float]:
    """Calculate the centroid (center of mass) of a list of atoms."""
    if not atoms:
        return (0.0, 0.0, 0.0)
    x = sum(a.x for a in atoms) / len(atoms)
    y = sum(a.y for a in atoms) / len(atoms)
    z = sum(a.z for a in atoms) / len(atoms)
    return (x, y, z)

def _distance_to_centroid(atom: Atom, centroid: Tuple[float, float, float]) -> float:
    """Calculate distance from an atom to a centroid point."""
    dx = atom.x - centroid[0]
    dy = atom.y - centroid[1]
    dz = atom.z - centroid[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def _filter_interactions_by_distance(interactions: List[Interaction], 
                                     ligand_centroid: Tuple[float, float, float],
                                     max_distance: float = MAX_INTERACTION_DISTANCE) -> List[Interaction]:
    """
    Filter interactions to only keep those where the protein atom
    is within max_distance of the ligand centroid.
    
    This prevents spurious long-distance interactions when multiple
    similar ligands are present in the structure.
    """
    filtered = []
    for inter in interactions:
        # Check if any protein atom is within distance of ligand centroid
        if inter.protein_atoms:
            min_dist = min(_distance_to_centroid(pa, ligand_centroid) 
                          for pa in inter.protein_atoms)
            if min_dist <= max_distance:
                filtered.append(inter)
        else:
            # No protein atoms to check (e.g., some special cases)
            filtered.append(inter)
    return filtered


# ============================================================================
# SINGLE LIGAND ANALYSIS
# ============================================================================

def analyze_single_ligand(protein_atoms: List[Atom],
                          ligand_atoms: List[Atom],
                          water_atoms: List[Atom],
                          ligand_id: str,
                          covalent_bonds: List[Dict] = None,
                          filter_interactions: bool = True) -> Dict[str, List[Interaction]]:
    """
    Analyze non-covalent interactions for a single ligand instance.
    
    Args:
        protein_atoms: List of protein atoms
        ligand_atoms: List of ligand atoms for this specific ligand
        water_atoms: List of water atoms
        ligand_id: Identifier string for this ligand instance
        covalent_bonds: Deprecated, ignored (kept for API compatibility)
        filter_interactions: If True, apply interaction priority filtering
        
    Returns:
        Dictionary mapping interaction types to lists of Interaction objects
    """
    results = {}
    
    # Note: covalent_bonds parameter is ignored - PLI focuses on non-covalent interactions
    
    # Calculate ligand centroid for distance filtering
    ligand_centroid = _calculate_centroid(ligand_atoms)
    
    # Detect all non-covalent interaction types
    hbonds = find_hydrogen_bonds(protein_atoms, ligand_atoms)
    if hbonds:
        # Filter to remove spurious long-distance interactions (multi-ligand scenarios)
        hbonds = _filter_interactions_by_distance(hbonds, ligand_centroid)
        if hbonds:
            results['Hydrogen Bonds'] = hbonds
    
    salt = find_salt_bridges(protein_atoms, ligand_atoms)
    if salt:
        salt = _filter_interactions_by_distance(salt, ligand_centroid)
        if salt:
            results['Salt Bridges'] = salt
    
    hydrophobic = find_hydrophobic_contacts(protein_atoms, ligand_atoms)
    if hydrophobic:
        hydrophobic = _filter_interactions_by_distance(hydrophobic, ligand_centroid)
        if hydrophobic:
            results['Hydrophobic'] = hydrophobic
    
    pi_stack = find_pi_stacking(protein_atoms, ligand_atoms)
    if pi_stack:
        pi_stack = _filter_interactions_by_distance(pi_stack, ligand_centroid)
        if pi_stack:
            results['π-π Stacking'] = pi_stack
    
    # find_cation_pi returns BOTH Cation-π AND π-Cation interactions
    # We need to separate them by their actual interaction_type
    cat_pi_all = find_cation_pi(protein_atoms, ligand_atoms)
    if cat_pi_all:
        cat_pi_all = _filter_interactions_by_distance(cat_pi_all, ligand_centroid)
        if cat_pi_all:
            # Separate by actual interaction type
            cation_pi = [i for i in cat_pi_all if i.interaction_type == 'Cation-π']
            pi_cation = [i for i in cat_pi_all if i.interaction_type == 'π-Cation']
            if cation_pi:
                results['Cation-π'] = cation_pi
            if pi_cation:
                results['π-Cation'] = pi_cation
    
    halogen = find_halogen_bonds(protein_atoms, ligand_atoms)
    if halogen:
        halogen = _filter_interactions_by_distance(halogen, ligand_centroid)
        if halogen:
            results['Halogen Bonds'] = halogen
    
    # Water bridges
    water_bridges = find_water_bridges(protein_atoms, ligand_atoms, water_atoms)
    if water_bridges:
        water_bridges = _filter_interactions_by_distance(water_bridges, ligand_centroid)
        if water_bridges:
            results['Water Bridges'] = water_bridges
    
    # XWH bridges (halogen-water-hydrogen bridges)
    xwh = find_xwh_bridges(protein_atoms, ligand_atoms, water_atoms)
    if xwh:
        xwh = _filter_interactions_by_distance(xwh, ligand_centroid)
        if xwh:
            results['XWH Bridges'] = xwh
    
    # Apply interaction priority filtering if enabled
    if filter_interactions and results:
        results = apply_filters_to_results(results)
    
    return results


# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================

def analyze(pdb_file: str, 
            ligand_name: str = None,
            ligand_chain: str = None, 
            ligand_resnum: int = None,
            include_cofactors: bool = True,
            exclude_cofactors: bool = False,
            exclude_glycans: bool = False,
            verbose: bool = True,
            filter_interactions: bool = True) -> Tuple[Dict, Set[str], List, bool, List[Atom]]:
    """
    Run complete interaction analysis on a PDB file.
    
    Args:
        pdb_file: Path to PDB file
        ligand_name: Optional specific ligand residue name
        ligand_chain: Optional specific ligand chain ID
        ligand_resnum: Optional specific ligand residue number
        include_cofactors: If True, include cofactors (default True)
        exclude_cofactors: If True, explicitly exclude cofactors
        exclude_glycans: If True, explicitly exclude glycans
        verbose: If True, print progress messages
        filter_interactions: If True, apply interaction priority filtering (default True)
        
    Returns:
        Tuple of:
        - all_results: Dictionary keyed by ligand ID containing interaction results
        - lig_names: Set of ligand residue names
        - lig_residues: List of (chain, resname, resnum) tuples
        - has_chains: Boolean indicating if file uses chain IDs
        - metal_atoms: List of metal atom objects
    """
    def _print(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    
    _print(f"\n{'='*60}")
    _print(" NCIP - NON-COVALENT INTERACTION PROFILER v1.0")
    _print(f"{'='*60}")
    
    # Parse PDB file
    protein, all_ligand_atoms, water, metal_atoms, lig_names, lig_residues, has_chains = parse_pdb(
        pdb_file, ligand_name, ligand_chain, ligand_resnum, 
        include_cofactors=include_cofactors,
        exclude_cofactors=exclude_cofactors,
        exclude_glycans=exclude_glycans
    )
    
    _print(f"\nFile: {pdb_file}")
    _print(f"Protein atoms: {len(protein)}")
    _print(f"Total ligand atoms: {len(all_ligand_atoms)}")
    _print(f"Metal ions: {len(metal_atoms)}")
    _print(f"Ligand types: {', '.join(lig_names) if lig_names else 'None'}")
    _print(f"Ligand instances: {len(lig_residues)}")
    
    if not all_ligand_atoms:
        _print("\nNo ligand found!")
        return {}, set(), [], False, []
    
    # Group ligand atoms by instance
    ligand_by_instance = defaultdict(list)
    for a in all_ligand_atoms:
        key = (a.chain_id, a.res_name, a.res_seq)
        ligand_by_instance[key].append(a)
    
    # Covalent bonds - currently not detected from PDB, placeholder for future
    covalent_bonds = []
    
    # Group covalent bonds by ligand instance
    covalent_by_ligand = defaultdict(list)
    for cb in covalent_bonds:
        key = (cb['lig_chain'], cb['lig_res'], cb['lig_resnum'])
        covalent_by_ligand[key].append(cb)
    
    # Analyze each ligand instance separately
    _print(f"\n{'-'*60}")
    _print(" Detecting interactions for each ligand...")
    _print(f"{'-'*60}")
    
    all_results = {}  # Keyed by ligand instance ID
    
    for (chain, resname, resnum), lig_atoms in ligand_by_instance.items():
        ligand_id = f"{chain}:{resname}{resnum}"
        
        # Get covalent bonds for this ligand
        lig_covalent = covalent_by_ligand.get((chain, resname, resnum), [])
        
        # Diagnostics
        lig_elements = defaultdict(int)
        for a in lig_atoms:
            lig_elements[a.element] += 1
        
        _print(f"\n>>> Ligand: {ligand_id}")
        _print(f"    Atoms: {len(lig_atoms)}, Composition: {dict(lig_elements)}")
        if lig_covalent:
            _print(f"    Covalent attachment: {len(lig_covalent)} bond(s)")
        
        # Analyze this ligand instance (pass filter_interactions parameter)
        results = analyze_single_ligand(protein, lig_atoms, water, ligand_id, lig_covalent,
                                        filter_interactions=filter_interactions)
        
        total = sum(len(v) for v in results.values())
        _print(f"    Interactions found: {total}")
        
        if results:
            all_results[ligand_id] = {
                'results': results,
                'residue_info': (chain, resname, resnum),
                'has_chains': has_chains
            }
    
    # Analyze metal coordination
    if metal_atoms:
        _print(f"\n{'-'*60}")
        _print(" Detecting metal coordination...")
        _print(f"{'-'*60}")
        
        metal_results = detect_metal_coordination(metal_atoms, protein, all_ligand_atoms, water)
        
        for metal_id, coordinations in metal_results.items():
            _print(f"\n>>> Metal: {metal_id}")
            if coordinations:
                cn = coordinations[0].extra_info.get('coordination_number', 0)
                geom = coordinations[0].extra_info.get('geometry', 'unknown')
                _print(f"    Coordination number: {cn} ({geom})")
                
                # Group by source type
                by_source = {'protein': 0, 'water': 0, 'ligand': 0, 'protein_backbone': 0}
                for coord in coordinations:
                    src = coord.extra_info.get('source', 'unknown')
                    if src in by_source:
                        by_source[src] += 1
                sources = [f"{v} {k}" for k, v in by_source.items() if v > 0]
                _print(f"    Coordinating: {', '.join(sources)}")
            
            # Store metal results separately
            metal_chain = metal_id.split(':')[0]
            metal_res = metal_id.split(':')[1]
            metal_resnum = int(''.join(filter(str.isdigit, metal_res)))
            
            all_results[metal_id] = {
                'results': {'Metal Coordination': coordinations},
                'residue_info': (metal_chain, metal_res[:2], metal_resnum),
                'has_chains': has_chains,
                'is_metal': True
            }
        
        # Merge bridging metals into ligand results (creates synthetic interactions
        # with ligand_atoms pointing to actual ligand atoms, not metal atoms)
        _merge_bridging_metals(all_results, metal_results, metal_atoms)
    
    return all_results, lig_names, lig_residues, has_chains, metal_atoms


# ============================================================================
# BRIDGING METAL MERGE
# ============================================================================

def _merge_bridging_metals(all_results: Dict, metal_results: Dict, 
                           metal_atoms: List[Atom]) -> None:
    """
    Merge bridging metal coordination into ligand results.
    
    A bridging metal (e.g., ZN) coordinates BOTH protein residues and ligand atoms.
    For the 2D diagram, we want these shown as metal-mediated bonds on the ligand's
    diagram — similar to how water bridges show a water dot between residue and ligand.
    
    For each bridging metal:
    1. Find which ligand atom(s) it coordinates
    2. For each protein-side coordination, create a synthetic Interaction with:
       - ligand_atoms = [the bridging ligand atom] (so diagram can map to ligand)
       - extra_info includes metal symbol, metal-protein distance, metal-ligand distance
    3. Merge these into the ligand's results under 'Metal Coordination'
    4. Remove the standalone metal entry from all_results
    
    Args:
        all_results: Main results dict (modified in-place)
        metal_results: Dict from detect_metal_coordination {metal_id: [Interaction, ...]}
        metal_atoms: List of metal Atom objects
    """
    metals_to_remove = []
    
    for metal_id, coordinations in metal_results.items():
        if not coordinations:
            continue
        
        # Separate protein-side vs ligand-side coordinations
        protein_coords = []
        ligand_coords = []
        water_coords = []
        
        for coord in coordinations:
            source = coord.extra_info.get('source', '')
            if source in ('protein', 'protein_backbone'):
                protein_coords.append(coord)
            elif source == 'ligand':
                ligand_coords.append(coord)
            elif source == 'water':
                water_coords.append(coord)
        
        # Only merge if metal bridges protein AND ligand
        if not protein_coords or not ligand_coords:
            continue
        
        # Find which ligand entry this metal bridges to
        # Use the ligand-side coordination to identify the target ligand
        target_ligand_id = None
        bridging_lig_atom = None
        metal_lig_dist = 0.0
        
        for lig_coord in ligand_coords:
            # The ligand atom coordinating the metal
            if lig_coord.ligand_atoms:
                bridging_lig_atom = lig_coord.ligand_atoms[0]
                metal_lig_dist = lig_coord.distance
                # Find matching ligand in all_results
                for lid, ldata in all_results.items():
                    if ldata.get('is_metal'):
                        continue
                    res_info = ldata.get('residue_info', ())
                    if len(res_info) == 3:
                        chain, resname, resnum = res_info
                        if (bridging_lig_atom.chain_id == chain and 
                            bridging_lig_atom.res_name == resname and 
                            bridging_lig_atom.res_seq == resnum):
                            target_ligand_id = lid
                            break
            if target_ligand_id:
                break
        
        if not target_ligand_id or not bridging_lig_atom:
            continue
        
        # Get metal element symbol  
        metal_element = coordinations[0].extra_info.get('metal', 'M')
        metal_symbol = metal_element.capitalize()  # e.g., 'ZN' -> 'Zn'
        
        # Create synthetic interactions for each protein-side coordination
        # These will appear as metal-mediated bonds on the ligand's diagram
        target_results = all_results[target_ligand_id]['results']
        if 'Metal Coordination' not in target_results:
            target_results['Metal Coordination'] = []
        
        print(f"    Bridging metal {metal_id}: coordinates {len(ligand_coords)} ligand atom(s) + {len(protein_coords)} protein + {len(water_coords)} water")
        print(f"    → Merged {metal_id} coordination into {target_ligand_id}")
        
        for prot_coord in protein_coords:
            # Create synthetic interaction: protein residue ↔ ligand atom, mediated by metal
            synthetic = Interaction(
                interaction_type="Metal Coordination",
                protein_residue=prot_coord.protein_residue,
                distance=prot_coord.distance,  # Metal-protein distance (shown on residue side)
                protein_atoms=prot_coord.protein_atoms,
                ligand_atoms=[bridging_lig_atom],  # KEY: point to the actual ligand atom
                confidence=prot_coord.confidence,
                extra_info={
                    **prot_coord.extra_info,
                    'metal_bridged': True,
                    'metal_symbol': metal_symbol,
                    'metal_id': metal_id,
                    'metal_prot_dist': prot_coord.distance,
                    'metal_lig_dist': metal_lig_dist,
                }
            )
            target_results['Metal Coordination'].append(synthetic)
        
        # No separate "direct" metal-to-ligand interaction needed —
        # the metal dot on each bridging bond line already represents it
        # (same pattern as water bridges: no separate HOH sphere)
        
        metals_to_remove.append(metal_id)
    
    # Remove merged metal entries
    for mid in metals_to_remove:
        if mid in all_results:
            del all_results[mid]


# ============================================================================
# OUTPUT FUNCTIONS
# ============================================================================

def print_results(all_results: Dict) -> None:
    """
    Print formatted results for all ligands and metals.
    
    Args:
        all_results: Dictionary of analysis results
    """
    if not all_results:
        print("\nNo interactions detected.")
        return
    
    # Separate metals from ligands
    ligand_results = {k: v for k, v in all_results.items() if not v.get('is_metal', False)}
    metal_results = {k: v for k, v in all_results.items() if v.get('is_metal', False)}
    
    # Calculate grand total
    grand_total = 0
    for ligand_id, data in all_results.items():
        grand_total += sum(len(v) for v in data['results'].values())
    
    print(f"\n{'='*60}")
    print(f" RESULTS")
    print(f"{'='*60}")
    print(f"\nTotal ligands analyzed: {len(ligand_results)}")
    print(f"Total metal ions analyzed: {len(metal_results)}")
    print(f"Total interactions: {grand_total}")
    
    # Print results for each ligand
    for ligand_id, data in ligand_results.items():
        results = data['results']
        total = sum(len(v) for v in results.values())
        
        print(f"\n{'='*60}")
        print(f" LIGAND: {ligand_id}")
        print(f"{'='*60}")
        print(f"\nInteractions: {total}")
        
        for itype, interactions in results.items():
            print(f"\n{'-'*60}")
            print(f" {itype} ({len(interactions)})")
            print(f"{'-'*60}")
            
            for i, inter in enumerate(interactions, 1):
                prot_str = inter.protein_residue
                lig_str = ", ".join(a.name for a in inter.ligand_atoms[:3])
                if len(inter.ligand_atoms) > 3:
                    lig_str += "..."
                
                print(f"\n  [{i}] {prot_str} <-> {lig_str}")
                print(f"      Distance: {inter.distance:.2f} Å")
                for k, v in inter.extra_info.items():
                    print(f"      {k}: {v}")
    
    # Print metal coordination results
    if metal_results:
        print(f"\n{'='*60}")
        print(f" METAL COORDINATION")
        print(f"{'='*60}")
        
        for metal_id, data in metal_results.items():
            coordinations = data['results'].get('Metal Coordination', [])
            if not coordinations:
                continue
            
            cn = coordinations[0].extra_info.get('coordination_number', 0)
            geom = coordinations[0].extra_info.get('geometry', 'unknown')
            
            print(f"\n{'-'*60}")
            print(f" {metal_id} - CN: {cn} ({geom})")
            print(f"{'-'*60}")
            
            for i, coord in enumerate(coordinations, 1):
                coord_atom = coord.extra_info.get('coord_atom', '?')
                coord_type = coord.extra_info.get('coord_type', '?')
                res_id = coord.protein_residue
                
                print(f"  [{i}] {res_id}:{coord_atom} ({coord_type}) - {coord.distance:.2f} Å")
    
    # Grand Summary
    print(f"\n{'='*60}")
    print(f" GRAND SUMMARY")
    print(f"{'='*60}")
    print(f"\n  {'Entity':<25} {'Interactions':>12}")
    print(f"  {'-'*37}")
    for ligand_id, data in all_results.items():
        total = sum(len(v) for v in data['results'].values())
        is_metal = "🔷" if data.get('is_metal', False) else ""
        print(f"  {is_metal}{ligand_id:<24} {total:>12}")
    print(f"  {'-'*37}")
    print(f"  {'TOTAL':<25} {grand_total:>12}\n")


def save_tabular_data(pdb_file: str, 
                      all_results: Dict, 
                      output_file: str) -> None:
    """
    Save interaction data in tabular text format.
    
    Args:
        pdb_file: Path to source PDB file
        all_results: Dictionary of analysis results
        output_file: Path to output text file
    """
    if not all_results:
        print("\n✗ No interactions to save.")
        return
    
    # Flatten results for single-file output
    results = {}
    for ligand_id, data in all_results.items():
        for itype, interactions in data['results'].items():
            if itype not in results:
                results[itype] = []
            results[itype].extend(interactions)
    
    with open(output_file, 'w') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("NON-COVALENT INTERACTION ANALYSIS - TABULAR DATA\n")
        f.write("="*80 + "\n\n")
        f.write(f"PDB File: {os.path.basename(pdb_file)}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Interactions: {sum(len(v) for v in results.values())}\n\n")
        
        # Summary table
        f.write("-"*80 + "\n")
        f.write("SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Interaction Type':<25} {'Count':>10}\n")
        f.write("-"*35 + "\n")
        for itype, interactions in results.items():
            f.write(f"{itype:<25} {len(interactions):>10}\n")
        f.write("-"*35 + "\n")
        f.write(f"{'TOTAL':<25} {sum(len(v) for v in results.values()):>10}\n\n")
        
        # Detailed tables for each interaction type
        for itype, interactions in results.items():
            f.write("\n" + "="*80 + "\n")
            f.write(f"{itype.upper()} ({len(interactions)})\n")
            f.write("="*80 + "\n\n")
            
            if itype == 'Hydrogen Bonds':
                f.write(f"{'No.':<5} {'Residue':<15} {'Prot Atom':<12} {'Lig Atom':<12} {'Dist(Å)':<10} {'Role':<30}\n")
                f.write("-"*84 + "\n")
                for i, inter in enumerate(interactions, 1):
                    prot_atoms = ", ".join(a.name for a in inter.protein_atoms[:2]) if inter.protein_atoms else '-'
                    lig_atoms = ", ".join(a.name for a in inter.ligand_atoms[:2])
                    role = inter.extra_info.get('role', '-')[:28]
                    f.write(f"{i:<5} {inter.protein_residue:<15} {prot_atoms:<12} {lig_atoms:<12} {inter.distance:<10.2f} {role:<30}\n")
                    
            elif itype == 'Salt Bridges':
                f.write(f"{'No.':<5} {'Residue':<15} {'Prot Atom':<12} {'Lig Atom':<12} {'Dist(Å)':<10} {'Type':<25}\n")
                f.write("-"*79 + "\n")
                for i, inter in enumerate(interactions, 1):
                    prot_atoms = ", ".join(a.name for a in inter.protein_atoms[:2]) if inter.protein_atoms else '-'
                    lig_atoms = ", ".join(a.name for a in inter.ligand_atoms[:2])
                    sb_type = inter.extra_info.get('type', '-')
                    f.write(f"{i:<5} {inter.protein_residue:<15} {prot_atoms:<12} {lig_atoms:<12} {inter.distance:<10.2f} {sb_type:<25}\n")
                    
            elif itype == 'Hydrophobic':
                f.write(f"{'No.':<5} {'Residue':<15} {'Prot Atom':<12} {'Lig Atom':<12} {'Dist(Å)':<10}\n")
                f.write("-"*54 + "\n")
                for i, inter in enumerate(interactions, 1):
                    prot_atoms = ", ".join(a.name for a in inter.protein_atoms[:2]) if inter.protein_atoms else '-'
                    lig_atoms = ", ".join(a.name for a in inter.ligand_atoms[:2])
                    f.write(f"{i:<5} {inter.protein_residue:<15} {prot_atoms:<12} {lig_atoms:<12} {inter.distance:<10.2f}\n")
                    
            elif itype == 'π-π Stacking':
                f.write(f"{'No.':<5} {'Residue':<15} {'Type':<16} {'Dist(Å)':<10} {'Angle(°)':<10}\n")
                f.write("-"*56 + "\n")
                for i, inter in enumerate(interactions, 1):
                    pi_type = inter.extra_info.get('type', '-')
                    ang = inter.extra_info.get('angle', '-')
                    f.write(f"{i:<5} {inter.protein_residue:<15} {pi_type:<16} {inter.distance:<10.2f} {ang:<10}\n")
                    
            elif itype == 'Cation-π':
                f.write(f"{'No.':<5} {'Residue':<15} {'Dist(Å)':<10} {'Type':<40}\n")
                f.write("-"*70 + "\n")
                for i, inter in enumerate(interactions, 1):
                    cat_type = inter.extra_info.get('type', '-')
                    f.write(f"{i:<5} {inter.protein_residue:<15} {inter.distance:<10.2f} {cat_type:<40}\n")
                    
            elif itype == 'Halogen Bonds':
                f.write(f"{'No.':<5} {'Residue':<15} {'Halogen':<12} {'Dist(Å)':<10} {'Angle(°)':<10} {'Type':<20}\n")
                f.write("-"*72 + "\n")
                for i, inter in enumerate(interactions, 1):
                    halogen = inter.ligand_atoms[0].name if inter.ligand_atoms else '-'
                    ang = inter.extra_info.get('angle', '-')
                    hal_type = inter.extra_info.get('type', '-')
                    f.write(f"{i:<5} {inter.protein_residue:<15} {halogen:<12} {inter.distance:<10.2f} {ang:<10} {hal_type:<20}\n")
                    
            else:
                # Generic format
                f.write(f"{'No.':<5} {'Residue':<15} {'Dist(Å)':<10} {'Details':<50}\n")
                f.write("-"*80 + "\n")
                for i, inter in enumerate(interactions, 1):
                    details = "; ".join(f"{k}={v}" for k, v in list(inter.extra_info.items())[:3])
                    f.write(f"{i:<5} {inter.protein_residue:<15} {inter.distance:<10.2f} {details:<50}\n")
            
            f.write("\n")
        
        # Footer
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"\n✓ Tabular data saved: {output_file}")
