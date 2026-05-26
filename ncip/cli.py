#!/usr/bin/env python3
"""
NCIP - Non-Covalent Interaction Profiler — Command-Line Interface
==================================================================

Analyze non-covalent interactions between proteins and ligands from PDB files.

Usage:
    ncip structure.pdb [options]

Examples:
    ncip 1abc.pdb                         # Analyze all ligands (interactive)
    ncip 1abc.pdb -l ATP                  # Specific ligand (no prompt)
    ncip 1abc.pdb --pml output.pml        # PyMOL script
    ncip 1abc.pdb --json results.json     # JSON output
    ncip 1abc.pdb --csv results.csv       # CSV output
    ncip --fetch 1abc                     # Download from RCSB and analyze
    ncip --fetch 1abc --list              # Show available molecules
    ncip --batch pdbs/ -o results/        # Batch processing
"""

import argparse
import csv
import json
import os
import sys
import tempfile
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict

from . import __version__
from .core import analyze, print_results, save_tabular_data
from .visualizer import save_script_only
from .parser import parse_pdb, detect_molecules


def fetch_pdb(pdb_id: str, output_dir: str = None) -> str:
    """
    Download a PDB file from RCSB.

    Args:
        pdb_id: 4-character PDB ID
        output_dir: Directory to save the file (default: temp dir)

    Returns:
        Path to downloaded PDB file
    """
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"Invalid PDB ID: {pdb_id} (must be 4 characters)")

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    else:
        output_path = os.path.join(os.getcwd(), f"{pdb_id}.pdb")

    try:
        urllib.request.urlretrieve(url, output_path)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise FileNotFoundError(f"PDB ID {pdb_id} not found on RCSB") from None
        raise
    except urllib.error.URLError as e:
        raise ConnectionError(f"Failed to download {pdb_id}: {e}") from None

    return output_path


def list_molecules(pdb_file: str, exclude_cofactors: bool = False,
                   exclude_glycans: bool = False) -> list:
    """Detect and display all non-protein molecules in a PDB file.
    
    Returns the list of molecule dicts from detect_molecules(), after
    applying exclusion filters.
    """
    molecules = detect_molecules(pdb_file)
    
    # Apply exclusions
    if exclude_cofactors:
        molecules = [m for m in molecules if m['category'] != 'Cofactor']
    if exclude_glycans:
        molecules = [m for m in molecules if m['category'] != 'Glycan']
    # Always exclude Buffer/Salt from the interactive table
    molecules = [m for m in molecules if m['category'] != 'Buffer/Salt']
    
    if not molecules:
        print("No molecules found.")
        return []
    
    # Count categories
    n_lig = sum(1 for m in molecules if m['category'] == 'Ligand')
    n_cof = sum(1 for m in molecules if m['category'] == 'Cofactor')
    n_gly = sum(1 for m in molecules if m['category'] == 'Glycan')
    
    print(f"\n{'='*60}")
    print(f" NCIP v{__version__} — Molecule Detection")
    print(f"{'='*60}")
    print(f"\n  {'#':<4} {'Chain':<7} {'Name':<7} {'ResNum':<8} {'Atoms':<7} {'Category'}")
    print(f" {'─'*50}")
    
    for i, mol in enumerate(molecules, 1):
        print(f"  {i:<4} {mol['chain']:<7} {mol['resname']:<7} "
              f"{mol['resnum']:<8} {mol['atoms']:<7} {mol['category'].lower()}")
    
    # Summary line
    parts = [f"Ligands: {n_lig}", f"Cofactors: {n_cof}"]
    if exclude_cofactors:
        parts[1] += " (excluded)"
    parts.append(f"Glycans: {n_gly}")
    if exclude_glycans:
        parts[-1] += " (excluded)"
    print(f"\n  {' | '.join(parts)}")
    
    return molecules


def interactive_select(molecules: list) -> list:
    """Prompt user to select molecules for analysis.
    
    Returns filtered list of molecule dicts to analyze.
    """
    if not sys.stdin.isatty():
        # Non-interactive: analyze all
        return molecules
    
    n = len(molecules)
    if n == 0:
        return []
    
    prompt = f"\nAnalyze [{', '.join(str(i) for i in range(1, min(n+1, 4)))}{'...' if n > 3 else ''}, all, or q]: "
    
    while True:
        try:
            choice = input(prompt).strip().lower()
        except (EOFError, KeyboardInterrupt):
            print()
            sys.exit(0)
        
        if choice in ('q', 'quit', 'exit'):
            sys.exit(0)
        if choice in ('all', 'a', ''):
            return molecules
        
        # Parse comma-separated numbers and ranges (e.g., "1,3" or "1-3")
        selected = []
        try:
            for part in choice.split(','):
                part = part.strip()
                if '-' in part:
                    start, end = part.split('-', 1)
                    for idx in range(int(start), int(end) + 1):
                        if 1 <= idx <= n:
                            selected.append(molecules[idx - 1])
                else:
                    idx = int(part)
                    if 1 <= idx <= n:
                        selected.append(molecules[idx - 1])
        except ValueError:
            print(f"  Invalid input. Enter numbers (1-{n}), 'all', or 'q'.")
            continue
        
        if selected:
            return selected
        print(f"  No valid selection. Enter numbers (1-{n}), 'all', or 'q'.")


def load_threshold_config(config_path: str) -> dict:
    """Load threshold overrides from a YAML or JSON config file.
    
    Returns a flat dict mapping constants module variable names to values.
    """
    import json as _json
    
    with open(config_path, 'r') as f:
        content = f.read()
    
    # Try JSON first
    try:
        raw = _json.loads(content)
    except _json.JSONDecodeError:
        # Try YAML
        try:
            import yaml
            raw = yaml.safe_load(content)
        except ImportError:
            print("Error: PyYAML is required for YAML config files. "
                  "Install with: pip install pyyaml", file=sys.stderr)
            print("Alternatively, use JSON format.", file=sys.stderr)
            sys.exit(1)
    
    if not isinstance(raw, dict):
        print(f"Error: Config file must contain a mapping, got {type(raw).__name__}", 
              file=sys.stderr)
        sys.exit(1)
    
    # Map nested config keys to constants module variable names
    KEY_MAP = {
        ('hbond', 'dist_max'): 'HBOND_DISTANCE_MAX',
        ('hbond', 'don_angle_min'): 'HBOND_ANGLE_MIN',
        ('saltbridge', 'dist_max'): 'SALT_BRIDGE_DIST_MAX',
        ('hydrophobic', 'dist_min'): 'HYDROPHOBIC_DIST_MIN',
        ('hydrophobic', 'dist_max'): 'HYDROPHOBIC_DIST_MAX',
        ('pistack', 'centroid_dist_max'): 'PI_STACK_DIST_MAX',
        ('pistack', 'interplane_angle_parallel'): 'PI_STACK_ANGLE_PARALLEL',
        ('pistack', 'interplane_angle_tshaped'): 'PI_STACK_ANGLE_TSHAPE',
        ('pistack', 'elevation_angle_min'): 'PI_STACK_ELEVATION_MIN',
        ('pistack', 'vertical_dist_min'): 'PI_STACK_VERT_DIST_MIN',
        ('pistack', 'vertical_dist_max'): 'PI_STACK_VERT_DIST_MAX',
        ('cation_pi', 'dist_max'): 'CATION_PI_DIST_MAX',
        ('halogen', 'dist_cl'): 'HALOGEN_BOND_DIST_MAX',  # Use largest as max
        ('halogen', 'donor_angle_min'): 'HALOGEN_BOND_ANGLE_MIN',
        ('water_bridge', 'dist_max'): 'WATER_BRIDGE_DIST_MAX',
        ('water_bridge', 'omega_min'): 'WATER_BRIDGE_ANGLE_MIN',
        ('water_bridge', 'omega_max'): 'WATER_BRIDGE_ANGLE_MAX',
        ('xwh', 'cx_angle_min'): 'XWH_XBOND_ANGLE_MIN',
        ('metal', 'dist_max'): 'METAL_COORD_DIST_MAX',
    }
    
    overrides = {}
    for section_key, section_val in raw.items():
        if isinstance(section_val, dict):
            for param_key, param_val in section_val.items():
                lookup = (section_key, param_key)
                if lookup in KEY_MAP:
                    overrides[KEY_MAP[lookup]] = float(param_val)
        else:
            # Flat key format: e.g. "hbond_dist_max"
            pass  # Not supported, use nested format
    
    return overrides


def apply_thresholds(overrides: dict) -> None:
    """Apply threshold overrides to the constants module."""
    from . import constants
    for key, value in overrides.items():
        if hasattr(constants, key):
            setattr(constants, key, value)


def results_to_json(all_results: dict) -> list:
    """Convert results dict to JSON-serializable list."""
    output = []
    for ligand_id, data in all_results.items():
        entry = {
            'ligand_id': ligand_id,
            'is_metal': data.get('is_metal', False),
            'interactions': []
        }
        for itype, interactions in data['results'].items():
            for inter in interactions:
                record = {
                    'type': itype,
                    'protein_residue': inter.protein_residue,
                    'distance': round(inter.distance, 3),
                    'confidence': round(inter.confidence, 3),
                    'protein_atoms': [
                        {'name': a.name, 'element': a.element,
                         'chain': a.chain_id, 'res_name': a.res_name,
                         'res_seq': a.res_seq}
                        for a in inter.protein_atoms
                    ],
                    'ligand_atoms': [
                        {'name': a.name, 'element': a.element,
                         'chain': a.chain_id, 'res_name': a.res_name,
                         'res_seq': a.res_seq}
                        for a in inter.ligand_atoms
                    ],
                    'extra_info': {
                        k: v for k, v in inter.extra_info.items()
                        if not isinstance(v, (list, dict)) or k in ('observed_angles',)
                    }
                }
                entry['interactions'].append(record)
        output.append(entry)
    return output


def results_to_csv_rows(all_results: dict) -> list:
    """Convert results dict to flat CSV rows."""
    rows = []
    for ligand_id, data in all_results.items():
        for itype, interactions in data['results'].items():
            for inter in interactions:
                prot_atoms = "; ".join(
                    f"{a.chain_id}:{a.res_name}{a.res_seq}:{a.name}"
                    for a in inter.protein_atoms
                )
                lig_atoms = "; ".join(
                    f"{a.chain_id}:{a.res_name}{a.res_seq}:{a.name}"
                    for a in inter.ligand_atoms
                )
                row = {
                    'ligand_id': ligand_id,
                    'interaction_type': itype,
                    'protein_residue': inter.protein_residue,
                    'distance': f"{inter.distance:.3f}",
                    'confidence': f"{inter.confidence:.3f}",
                    'protein_atoms': prot_atoms,
                    'ligand_atoms': lig_atoms,
                }
                # Add selected extra_info fields
                for k in ('role', 'type', 'strength', 'angle', 'halogen',
                           'geometry', 'coordination_number', 'tau_parameter'):
                    if k in inter.extra_info:
                        row[k] = str(inter.extra_info[k])
                rows.append(row)
    return rows


def process_single(pdb_file: str, args, ligand_name: str = None,
                   ligand_chain: str = None, ligand_resnum: int = None) -> dict:
    """Run analysis on a single PDB file. Returns all_results dict.
    
    ligand_name/chain/resnum can override args for per-molecule analysis.
    """
    lig = ligand_name or args.ligand
    chain = ligand_chain or args.chain
    resnum = ligand_resnum or args.resnum
    
    all_results, detected_ligands, lig_residues, has_chains, metal_atoms = analyze(
        pdb_file,
        ligand_name=lig,
        ligand_chain=chain,
        ligand_resnum=resnum,
        exclude_cofactors=args.exclude_cofactors,
        exclude_glycans=args.exclude_glycans,
        verbose=not args.quiet,
        filter_interactions=not args.no_filter
    )

    if not args.quiet:
        print_results(all_results)

    # Determine ligand name for PyMOL
    ligand_for_pymol = lig
    if not ligand_for_pymol and detected_ligands:
        ligand_for_pymol = '+'.join(detected_ligands)

    ligand_info = {
        'residues': lig_residues,
        'has_chains': has_chains,
        'all_results': all_results,
        'metal_atoms': metal_atoms
    }

    # Save tabular text
    if args.txt:
        save_tabular_data(pdb_file, all_results, args.txt)

    # Save JSON
    if args.json:
        json_data = results_to_json(all_results)
        with open(args.json, 'w') as f:
            json.dump(json_data, f, indent=2)
        if not args.quiet:
            print(f"\n✓ JSON saved: {args.json}")

    # Save CSV
    if args.csv:
        rows = results_to_csv_rows(all_results)
        if rows:
            fieldnames = list(rows[0].keys())
            for row in rows:
                for k in row:
                    if k not in fieldnames:
                        fieldnames.append(k)
            with open(args.csv, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                writer.writeheader()
                writer.writerows(rows)
        else:
            with open(args.csv, 'w') as f:
                f.write("No interactions found\n")
        if not args.quiet:
            print(f"\n✓ CSV saved: {args.csv}")

    # PyMOL script generation
    if args.pml:
        pml_file = args.pml
        if not pml_file.endswith('.pml'):
            pml_file += '.pml'
        save_script_only(
            pdb_file, all_results, ligand_for_pymol,
            pml_file, ligand_info
        )

    return all_results


def batch_process(input_dir: str, output_dir: str, args) -> None:
    """Process all PDB files in a directory."""
    input_path = Path(input_dir)
    pdb_files = sorted(list(input_path.glob('*.pdb')) + list(input_path.glob('*.ent')))

    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"\nBatch processing {len(pdb_files)} files...")
    print(f"Output directory: {output_dir}\n")

    summary = []
    for i, pdb_file in enumerate(pdb_files, 1):
        stem = pdb_file.stem
        print(f"[{i}/{len(pdb_files)}] {pdb_file.name}...", end=" ", flush=True)

        try:
            # Temporarily override output paths
            orig_json = args.json
            orig_csv = args.csv
            orig_txt = args.txt
            orig_pml = args.pml

            if args.json:
                args.json = os.path.join(output_dir, f"{stem}.json")
            if args.csv:
                args.csv = os.path.join(output_dir, f"{stem}.csv")
            if args.txt:
                args.txt = os.path.join(output_dir, f"{stem}.txt")
            if args.pml:
                args.pml = os.path.join(output_dir, f"{stem}.pml")

            # Force quiet for batch
            orig_quiet = args.quiet
            args.quiet = True

            all_results = process_single(str(pdb_file), args)

            total = sum(
                len(v) for data in all_results.values()
                for v in data['results'].values()
            )
            print(f"{total} interactions")
            summary.append({'file': pdb_file.name, 'interactions': total})

            # Restore args
            args.json = orig_json
            args.csv = orig_csv
            args.txt = orig_txt
            args.pml = orig_pml
            args.quiet = orig_quiet

        except Exception as e:
            print(f"ERROR: {e}")
            summary.append({'file': pdb_file.name, 'interactions': -1, 'error': str(e)})

    # Write batch summary
    summary_path = os.path.join(output_dir, "batch_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*50}")
    print(f"Batch complete: {len(pdb_files)} files processed")
    print(f"Summary: {summary_path}")


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        prog='ncip',
        usage='%(prog)s [-h] [-v] [<.pdb>] [--fetch PDB_ID] [--batch DIR]\n'
              '            [-l NAME] [-c ID] [-r NUM]\n'
              '            [--exclude-cofactors] [--exclude-glycans]\n'
              '            [--no-filter] [--list] [--config FILE]\n'
              '            [--json [FILE]] [--csv [FILE]] [--txt [FILE]]\n'
              '            [--pml [FILE]] [-o DIR] [-q]',
        description=f"""
                 NCI-Profiler (NCIP) v{__version__}

SYNOPSIS

  %(prog)s structure.pdb [-l NAME] [--json FILE] [--csv FILE] [--pml FILE]
  %(prog)s --fetch PDB_ID [PDB_ID ...] [-l NAME] [--json FILE]
  %(prog)s --batch DIR/ -o results/ [--json] [--csv]
  %(prog)s structure.pdb --list
  %(prog)s --fetch PDB_ID --list

DESCRIPTION

  %(prog)s reads a PDB coordinate file (or fetches one from RCSB by its
  4-character accession code), identifies all small-molecule ligand instances
  in the structure, and detects non-covalent interactions between each ligand
  and surrounding protein residues, water molecules, and metal ions.

  Nine interaction types are detected using literature-validated geometric
  criteria: hydrogen bonds (Jeffrey 1997), salt bridges (Barlow & Thornton
  1983), hydrophobic contacts (Bissantz et al. 2010), pi-pi stacking
  (Calinsky & Levy 2024), cation-pi (Gallivan & Dougherty 1999), halogen bonds
  (Auffinger et al. 2004), water-mediated hydrogen bonds (Jiang et al.
  2005), halogen-water-hydrogen bridges (Zhou et al. 2010), and metal
  coordination with geometry fitting (Harding 2001).

  By default, all detected molecules (ligands, cofactors, glycans) are
  shown in an interactive selection table. Crystallization artifacts (PEG,
  GOL, SO4, etc.) are always excluded. Use --exclude-cofactors and/or
  --exclude-glycans to skip those categories, or -l to target a specific
  molecule directly.

  An 8-stage hierarchical priority filter removes redundant lower-priority
  interactions. Disable with --no-filter for raw detection output.

  Results are printed to the console by default. Structured output can be
  saved as JSON (--json), CSV (--csv), or formatted text (--txt). For 3D
  visualization, %(prog)s generates PyMOL scripts (.pml) with color-coded
  interaction dashes, distance labels, and toggleable groups per ligand.

  For high-throughput work, --batch processes all PDB files in a directory
  and writes per-structure output files plus a batch_summary.json.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=None
    )

    # Input source
    input_group = parser.add_argument_group('input')
    input_group.add_argument(
        'pdb_file',
        nargs='?',
        help='Input PDB file'
    )
    input_group.add_argument(
        '--fetch',
        nargs='+',
        metavar='PDB_ID',
        help='Download PDB(s) from RCSB (one or more 4-letter IDs)'
    )
    input_group.add_argument(
        '--batch',
        metavar='DIR',
        help='Process all PDB files in a directory'
    )

    # Molecule selection
    mol_group = parser.add_argument_group('molecule selection')
    mol_group.add_argument(
        '-l', '--ligand',
        metavar='NAME',
        help='Target a specific molecule by residue name'
    )
    mol_group.add_argument(
        '-c', '--chain',
        metavar='ID',
        help='Specify chain ID'
    )
    mol_group.add_argument(
        '-r', '--resnum',
        type=int,
        metavar='NUM',
        help='Specify residue number'
    )
    mol_group.add_argument(
        '--list',
        action='store_true',
        help='List detected molecules and exit'
    )
    mol_group.add_argument(
        '--exclude-cofactors',
        action='store_true',
        help='Skip cofactors (NAD, FAD, HEM, ATP, etc.)'
    )
    mol_group.add_argument(
        '--exclude-glycans',
        action='store_true',
        help='Skip glycans (NAG, MAN, GAL, etc.)'
    )

    # Analysis options
    analysis_group = parser.add_argument_group('analysis options')
    analysis_group.add_argument(
        '--no-filter',
        action='store_true',
        help='Disable 8-step interaction priority filter'
    )
    analysis_group.add_argument(
        '--config',
        metavar='FILE',
        help='Load custom thresholds from YAML/JSON file'
    )

    # Threshold overrides
    thresh_group = parser.add_argument_group('threshold overrides')
    thresh_group.add_argument('--hbond-dist-max', type=float, metavar='FLOAT',
                              help='H-bond D···A distance (default: 3.5 Å)')
    thresh_group.add_argument('--hbond-angle-min', type=float, metavar='FLOAT',
                              help='H-bond donor angle (default: 120°)')
    thresh_group.add_argument('--salt-dist-max', type=float, metavar='FLOAT',
                              help='Salt bridge distance (default: 5.5 Å)')
    thresh_group.add_argument('--hydro-dist-min', type=float, metavar='FLOAT',
                              help='Hydrophobic min distance (default: 3.3 Å)')
    thresh_group.add_argument('--hydro-dist-max', type=float, metavar='FLOAT',
                              help='Hydrophobic max distance (default: 4.5 Å)')
    thresh_group.add_argument('--pistack-dist-max', type=float, metavar='FLOAT',
                              help='π–π centroid distance (default: 5.0 Å)')
    thresh_group.add_argument('--catpi-dist-max', type=float, metavar='FLOAT',
                              help='Cation–π distance (default: 6.0 Å)')
    thresh_group.add_argument('--halogen-dist-max', type=float, metavar='FLOAT',
                              help='Halogen bond distance (default: 4.0 Å)')
    thresh_group.add_argument('--halogen-angle-min', type=float, metavar='FLOAT',
                              help='Halogen donor angle (default: 140°)')
    thresh_group.add_argument('--water-dist-max', type=float, metavar='FLOAT',
                              help='Water bridge distance (default: 3.5 Å)')
    thresh_group.add_argument('--water-omega-min', type=float, metavar='FLOAT',
                              help='Water bridge ω min (default: 80°)')
    thresh_group.add_argument('--water-omega-max', type=float, metavar='FLOAT',
                              help='Water bridge ω max (default: 140°)')
    thresh_group.add_argument('--xwh-angle-min', type=float, metavar='FLOAT',
                              help='XWH C–X···Ow angle (default: 130°)')
    thresh_group.add_argument('--metal-dist-max', type=float, metavar='FLOAT',
                              help='Metal coordination distance (default: 2.8 Å)')

    # Output options
    out_group = parser.add_argument_group('output options')
    out_group.add_argument(
        '--json',
        metavar='FILE',
        nargs='?',
        const='interactions.json',
        help='Save results as JSON (default: interactions.json)'
    )
    out_group.add_argument(
        '--csv',
        metavar='FILE',
        nargs='?',
        const='interactions.csv',
        help='Save results as CSV (default: interactions.csv)'
    )
    out_group.add_argument(
        '--txt',
        metavar='FILE',
        nargs='?',
        const='interactions.txt',
        help='Save detailed text report (default: interactions.txt)'
    )
    out_group.add_argument(
        '--pml',
        metavar='FILE',
        nargs='?',
        const='interactions.pml',
        help='Generate PyMOL script (.pml)'
    )
    out_group.add_argument(
        '-o', '--output-dir',
        metavar='DIR',
        default='.',
        help='Output directory for batch mode (default: current dir)'
    )

    # General
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='Suppress progress messages'
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )

    args = parser.parse_args()

    # ── Apply threshold overrides ──────────────────────────────────────
    # Priority: inline flags > config file > defaults
    if args.config:
        config_overrides = load_threshold_config(args.config)
        apply_thresholds(config_overrides)
    
    # Inline flags override config file values
    inline_map = {
        'hbond_dist_max': 'HBOND_DISTANCE_MAX',
        'hbond_angle_min': 'HBOND_ANGLE_MIN',
        'salt_dist_max': 'SALT_BRIDGE_DIST_MAX',
        'hydro_dist_min': 'HYDROPHOBIC_DIST_MIN',
        'hydro_dist_max': 'HYDROPHOBIC_DIST_MAX',
        'pistack_dist_max': 'PI_STACK_DIST_MAX',
        'catpi_dist_max': 'CATION_PI_DIST_MAX',
        'halogen_dist_max': 'HALOGEN_BOND_DIST_MAX',
        'halogen_angle_min': 'HALOGEN_BOND_ANGLE_MIN',
        'water_dist_max': 'WATER_BRIDGE_DIST_MAX',
        'water_omega_min': 'WATER_BRIDGE_ANGLE_MIN',
        'water_omega_max': 'WATER_BRIDGE_ANGLE_MAX',
        'xwh_angle_min': 'XWH_XBOND_ANGLE_MIN',
        'metal_dist_max': 'METAL_COORD_DIST_MAX',
    }
    inline_overrides = {}
    for arg_name, const_name in inline_map.items():
        val = getattr(args, arg_name, None)
        if val is not None:
            inline_overrides[const_name] = val
    if inline_overrides:
        apply_thresholds(inline_overrides)

    # ── Determine input PDB file ──────────────────────────────────────
    pdb_file = None

    if args.batch:
        batch_process(args.batch, args.output_dir, args)
        sys.exit(0)

    if args.fetch:
        pdb_ids = args.fetch
        fetched_files = []
        for pid in pdb_ids:
            try:
                if not args.quiet:
                    print(f"Downloading {pid.upper()} from RCSB...")
                path = fetch_pdb(pid)
                if not args.quiet:
                    print(f"Saved to: {path}")
                fetched_files.append(path)
            except (FileNotFoundError, ConnectionError, ValueError) as e:
                print(f"Error fetching {pid}: {e}", file=sys.stderr)

        if not fetched_files:
            print("Error: No PDB files were downloaded.", file=sys.stderr)
            sys.exit(1)

        for pdb_file in fetched_files:
            if not args.quiet and len(fetched_files) > 1:
                print(f"\n{'='*60}")
                print(f" Processing: {os.path.basename(pdb_file)}")
                print(f"{'='*60}")

            _run_analysis(pdb_file, args)
        sys.exit(0)

    elif args.pdb_file:
        pdb_file = args.pdb_file
    else:
        parser.print_help()
        sys.exit(1)

    if not Path(pdb_file).exists():
        print(f"Error: File not found: {pdb_file}", file=sys.stderr)
        sys.exit(1)

    _run_analysis(pdb_file, args)


def _run_analysis(pdb_file: str, args) -> None:
    """Detect molecules and run analysis with interactive selection."""
    
    # List mode: show molecules and exit
    if args.list:
        list_molecules(pdb_file, args.exclude_cofactors, args.exclude_glycans)
        return

    # If a specific ligand is given, skip the interactive table
    if args.ligand:
        try:
            process_single(pdb_file, args)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            if not args.quiet:
                import traceback
                traceback.print_exc()
            sys.exit(1)
        return
    
    # Interactive molecule detection
    if not args.quiet:
        molecules = list_molecules(pdb_file, args.exclude_cofactors, args.exclude_glycans)
    else:
        molecules = detect_molecules(pdb_file)
        if args.exclude_cofactors:
            molecules = [m for m in molecules if m['category'] != 'Cofactor']
        if args.exclude_glycans:
            molecules = [m for m in molecules if m['category'] != 'Glycan']
        molecules = [m for m in molecules if m['category'] != 'Buffer/Salt']
    
    if not molecules:
        if not args.quiet:
            print("No molecules found for analysis.")
        return
    
    # Select molecules
    if args.quiet:
        # Quiet mode: analyze all without prompting
        selected = molecules
    else:
        selected = interactive_select(molecules)
    
    if not selected:
        return
    
    # Analyze each selected molecule
    for mol in selected:
        try:
            process_single(pdb_file, args,
                          ligand_name=mol['resname'],
                          ligand_chain=mol['chain'] if mol['chain'] != '-' else None,
                          ligand_resnum=mol['resnum'])
        except Exception as e:
            print(f"Error analyzing {mol['resname']}: {e}", file=sys.stderr)
            if not args.quiet:
                import traceback
                traceback.print_exc()


if __name__ == "__main__":
    main()
