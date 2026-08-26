#!/usr/bin/env python3
"""
download_hprc.py — Download HPRC assemblies from the manifest.

REAL-DATA STEP: Requires S3 access. Each assembly is ~3GB.
TODO (Ali/Khoi): Add wget fallback if aws CLI unavailable.
TODO (Michael): Verify assemblies match expected checksums.

Usage:
    python3 scripts/download_hprc.py                      # list only
    aws s3 cp --no-sign-request <uri> work/downloads/     # actual download
"""

import csv
import os


def main():
    manifest_path = "work/manifests/hprc_selected.csv"
    if not os.path.exists(manifest_path):
        print(f"Manifest not found: {manifest_path}")
        print("Run: python3 scripts/fetch_hprc_index.py")
        return

    with open(manifest_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("Manifest is empty. No assemblies to download.")
        return

    os.makedirs("work/downloads", exist_ok=True)
    print(f"Found {len(rows)} assemblies in manifest:\n")

    for row in rows:
        name = row.get("assembly_name", "unknown")
        uri = row.get("assembly_uri", "N/A")
        print(f"  {name}")
        print(f"    URI: {uri}")
        print(f"    Local: work/downloads/{name}.fa")
        print()

    print("To download manually:")
    for row in rows:
        uri = row.get("assembly_uri", "")
        name = row.get("assembly_name", "unknown")
        print(f'  aws s3 cp --no-sign-request {uri} work/downloads/{name}.fa')

    print()
    print("REAL-DATA STEP: This downloads ~3GB per assembly.")
    print("Run on DNAnexus or a machine with sufficient bandwidth/disk.")


if __name__ == "__main__":
    main()