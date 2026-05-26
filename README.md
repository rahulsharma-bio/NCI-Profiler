# NCIP — NCI-Profiler CLI

Standalone command-line tool for detecting non-covalent interactions between
proteins and small molecule ligands from PDB files. Detects nine interaction
types (hydrogen bonds, salt bridges, hydrophobic contacts, π-π stacking,
cation-π, halogen bonds, water bridges, XWH bridges, and metal coordination)
using literature-validated geometric criteria. Requires Python ≥ 3.8 and no
third-party dependencies.

## Installation

**Option 1 — Virtual environment (recommended)**

Isolated install, works on all systems with no PATH issues:

```bash
git clone https://github.com/your-repo/ncip-cli.git
cd ncip-cli
python3 -m venv venv
source venv/bin/activate        # On Windows: venv\Scripts\activate
pip install -e .
ncip --help
```

**Option 2 — User-wide install (global command)**

Installs `ncip` as a global command without a virtual environment:

```bash
pip install -e . --break-system-packages
```

Then add `~/.local/bin` to your PATH if not already present (one-time setup):
```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

After that, `ncip` is available globally in any terminal.

## Usage

```bash
# Basic analysis
ncip protein.pdb

# Specific ligand
ncip protein.pdb -l ATP

# Download from RCSB and analyze
ncip --fetch 1ABC

# List available ligands
ncip --fetch 1ABC --list

# JSON output
ncip protein.pdb --json results.json

# CSV output
ncip protein.pdb --csv results.csv

# PyMOL session
ncip protein.pdb --pml output.pml

# Batch processing
ncip --batch pdbs/ -o results/ --json --csv

# Exclude cofactors (NAD, FAD, heme, etc. are included by default)
ncip protein.pdb --exclude-cofactors

# Skip interaction priority filtering
ncip protein.pdb --no-filter
```

## Python API

```python
from ncip import analyze, print_results

results, ligands, residues, has_chains, metals = analyze("protein.pdb")
print_results(results)

# Analyze a specific ligand
results, *_ = analyze("structure.pdb", ligand_name="ATP")

# Access individual interactions
for ligand_id, data in results.items():
    for itype, interactions in data['results'].items():
        for inter in interactions:
            print(f"{itype}: {inter.protein_residue} ({inter.distance:.2f} Å)")
```

## Output Formats

NCIP supports five output formats: console (default), JSON (`--json`), CSV
(`--csv`), formatted text (`--txt`), and PyMOL script (`--pml`).

### PyMOL Script (`--pml`)

Color-coded visualization with dashed and solid interaction lines, distance labels, and
toggleable groups per ligand and interaction type. Requires PyMOL and NumPy
to be installed. If PyMOL is not found, falls back to saving a `.pml` script.
To open a `.pml` script manually in PyMOL:

```bash
pymol output.pml
```
> **Note:** The `.pml` script references the PDB file by path. For the
> visualization to load correctly, the PDB file must be in the same folder
> as the `.pml` script, or the path inside the script must be updated.

| Color      | Interaction type     |
|------------|----------------------|
| Yellow     | Hydrogen bonds       |
| Magenta    | Salt bridges         |
| Gray       | Hydrophobic contacts |
| Cyan       | π–π stacking         |
| Green      | Cation–π             |
| Orange     | Halogen bonds        |
| Light blue | Water / XWH bridges  |
| Purple     | Metal coordination   |

## Supported Interaction Types

| Interaction Type     | Key Geometric Criteria             | Reference                    |
|----------------------|------------------------------------|------------------------------|
| Hydrogen Bonds       | D–A ≤ 3.5 Å, angle ≥ 120°         | Jeffrey (1997)               |
| Salt Bridges         | Distance ≤ 4.0 Å                   | Barlow & Thornton (1983)     |
| Hydrophobic Contacts | 3.3–4.0 Å (C···C)                  | Bissantz et al. (2010)       |
| π–π Stacking         | Centroid ≤ 5.5 Å, angle-dependent  | Calinsky & Levy (2024)       |
| Cation–π             | Cation–centroid ≤ 6.0 Å            | Dougherty (1996)             |
| Halogen Bonds        | D ≤ 3.5 Å, C–X···Y ≥ 140°         | Scholfield et al. (2013)     |
| Water Bridges        | D ≤ 3.5 Å via water                | Barillari et al. (2007)      |
| XWH Bridges          | X···W + W···H bridge               | Zhou et al. (2010)           |
| Metal Coordination   | M–L ≤ 2.8 Å + geometry fitting     | Harding (2001)               |

## License

NCIP is released under the [GNU Affero General Public License v3.0](LICENSE).
You are free to use, modify, and distribute this software under the terms of
the AGPL-3.0. If you use NCIP in a network service, the source code must also
be made available to users of that service.
