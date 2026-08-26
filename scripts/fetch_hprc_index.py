#!/usr/bin/env python3
"""
fetch_hprc_index.py — Fetch HPRC Release 2 assembly index.

TODO (Michael/Khoi): Verify URL still works for HPRC Release 2.
TODO (Ali): Add retry logic and checksum validation.

Usage:
    python3 scripts/fetch_hprc_index.py
"""

import csv
import os
import sys
import urllib.request

HPRC_INDEX_URL = (
    "https://raw.githubusercontent.com/human-pangenomics/"
    "hprc_intermediate_assembly/main/data_tables/"
    "assembly_index/HPRC_Release2_v1.0.1_assembly_index.csv"
)

REQUESTED_ASSEMBLIES = [
    "HG00673_mat_hprc_r2_v1.0.1",
    "HG00673_pat_hprc_r2_v1.0.1",
    "HG00733_mat_hprc_r2_v1.0.1",
    "HG00733_pat_hprc_r2_v1.0.1",
]


def fetch_index(url: str = HPRC_INDEX_URL) -> list:
    print(f"Fetching HPRC Release 2 index from:\n  {url}")
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode("utf-8")
        reader = csv.DictReader(content.splitlines())
        rows = list(reader)
        print(f"Downloaded {len(rows)} rows from index")
        return rows
    except Exception as e:
        print(f"WARNING: Could not fetch HPRC index: {e}")
        print("This is expected if offline. Creating empty manifest.")
        return []


def main():
    rows = fetch_index()
    os.makedirs("work/manifests", exist_ok=True)

    if not rows:
        with open("work/manifests/hprc_selected.csv", "w", newline="") as f:
            f.write("sample,assembly_name,haplotype,assembly_uri,checksum_uri,local_filename\n")
        print("Created empty manifest (offline mode).")
        return

    selected = [r for r in rows if r.get("assembly_name") in REQUESTED_ASSEMBLIES]
    found_names = {r.get("assembly_name") for r in selected}
    for name in REQUESTED_ASSEMBLIES:
        if name not in found_names:
            print(f"  WARNING: {name} not found in HPRC index")

    with open("work/manifests/hprc_selected.csv", "w", newline="") as f:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for s in selected:
            writer.writerow(s)

    print(f"Selected {len(selected)}/{len(REQUESTED_ASSEMBLIES)} requested assemblies")
    print(f"Manifest: work/manifests/hprc_selected.csv")


if __name__ == "__main__":
    main()