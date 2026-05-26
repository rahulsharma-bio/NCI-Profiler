"""
Download all benchmark PDB structures from RCSB for validation tests.
Run once from the tests/ directory: python3 download_pdbs.py
"""

import urllib.request
import os
import time
import sys

RCSB_URL = "https://files.rcsb.org/download/{}.pdb"

# All 85 benchmark structures organized by interaction type
STRUCTURES = {
    "cation_pi": [
        "2JKH", "2BOK", "1ACJ", "2Y5H", "1UV6",
        "1UW6", "4MQS", "3GWV", "2B4L",
    ],
    "halogen_bonds": [
        "5U00", "1P5E", "3SHZ", "3SIE", "3DN3",
        "4LBS", "4OEW", "4X21",
    ],
    "hydrogen_bonds": [
        "2ZC9", "1PXP", "3EHX", "2C3L", "2C3K",
        "1OWE", "1HVR", "1FKB", "3ERT", "4HJO", "4EKL",
    ],
    "hydrophobic": [
        "5C5P", "2FRA", "1J4H", "2ONZ", "3KWJ", "1M17", "2YXJ",
    ],
    "metal_coordination": [
        "4HHB", "2BB6", "1A42", "1AZM", "3HS4", "1JK3", "1O86",
        "1T64", "1W22", "3OYA", "3THH", "5YXC", "6DZQ",
    ],
    "pi_stacking": [
        "1EVE", "1L8B", "2W26", "3I1Y", "4TVJ", "4BFQ",
        "1N5R", "18GS", "1DMW", "1TLL", "1F20",
    ],
    "salt_bridges": [
        "4JRV", "2Z7K", "1BJI", "1PPC", "3PBL",
        "7Q2N", "7Q19", "5MU6", "3ODU",
    ],
    "water_bridges": [
        "1HPX", "1Z6F", "5J1X", "1B56", "1EY3", "1KT4", "1TT2",
    ],
    "xwh_bridges": [
        "2UY3", "2VIV", "1G3M", "1G5F", "1UV5",
        "1N8U", "1SN5", "3E3D", "2B1P", "1WCC",
    ],
}

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def download_pdb(pdb_id, dest_dir):
    pdb_id_lower = pdb_id.lower()
    dest = os.path.join(dest_dir, f"{pdb_id_lower}.pdb")
    if os.path.exists(dest):
        print(f"  [skip] {pdb_id_lower}.pdb already exists")
        return True
    url = RCSB_URL.format(pdb_id.upper())
    try:
        urllib.request.urlretrieve(url, dest)
        size_kb = os.path.getsize(dest) // 1024
        print(f"  [ok]   {pdb_id_lower}.pdb ({size_kb} KB)")
        return True
    except Exception as e:
        print(f"  [FAIL] {pdb_id_lower}: {e}", file=sys.stderr)
        return False

def main():
    total = sum(len(v) for v in STRUCTURES.values())
    done = 0
    failed = []

    for folder, pdb_ids in STRUCTURES.items():
        dest_dir = os.path.join(BASE_DIR, folder)
        os.makedirs(dest_dir, exist_ok=True)
        print(f"\n{folder.upper()} ({len(pdb_ids)} structures)")
        for pdb_id in pdb_ids:
            ok = download_pdb(pdb_id, dest_dir)
            if ok:
                done += 1
            else:
                failed.append(pdb_id)
            time.sleep(0.3)  # polite delay for RCSB

    print(f"\n{'='*50}")
    print(f"Downloaded: {done}/{total}")
    if failed:
        print(f"Failed:     {', '.join(failed)}")
    else:
        print("All downloads successful.")

if __name__ == "__main__":
    main()
