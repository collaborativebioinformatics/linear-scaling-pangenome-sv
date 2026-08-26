#!/usr/bin/env python3
"""
fetch_hprc_index.py — Select HG00673/HG00733 assemblies from the official
HPRC Release 2 index and write a compact manifest with exact S3 URIs.

The official index lives at:
  https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/
  refs/heads/main/data_tables/assemblies_release2_v1.0.index.csv

Official columns: sample_id, haplotype, phasing, assembly_method,
assembly_method_version, assembly_date, assembly_name, source,
genbank_accession, assembly_md5, assembly_fai, assembly_gzi, assembly

Usage:
    python3 scripts/fetch_hprc_index.py
    # Writes work/manifests/hprc_selected.csv
    # Exits non-zero if any of the 4 requested assemblies are missing.
"""

import csv
import os
import sys
import urllib.request

HPRC_INDEX_URL = (
    "https://raw.githubusercontent.com/human-pangenomics/"
    "hprc_intermediate_assembly/refs/heads/main/data_tables/"
    "assemblies_release2_v1.0.index.csv"
)

# The 4 requested assemblies: sample + haplotype
REQUESTED = [
    ("HG00673", "maternal"),
    ("HG00673", "paternal"),
    ("HG00733", "maternal"),
    ("HG00733", "paternal"),
]


def fetch_index(url: str) -> list[dict]:
    """Download and parse the official HPRC Release 2 index CSV."""
    print(f"Fetching HPRC Release 2 index from:\n  {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "pangenome-parallel/0.1"})
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            content = resp.read().decode("utf-8")
    except Exception as e:
        print(f"FATAL: Cannot fetch HPRC index: {e}", file=sys.stderr)
        sys.exit(1)

    reader = csv.DictReader(content.splitlines())
    rows = list(reader)
    print(f"  Downloaded {len(rows)} records from index")
    return rows


def main():
    rows = fetch_index(HPRC_INDEX_URL)

    # Build lookup keyed by (sample_id, haplotype)
    lookup: dict[tuple[str, str], dict] = {}
    for r in rows:
        sid = r.get("sample_id", "").strip()
        hap = r.get("haplotype", "").strip()
        if sid and hap:
            lookup[(sid, hap)] = r

    # Select the 4 we need
    selected: list[dict] = []
    missing: list[tuple[str, str]] = []
    for sample, haplotype in REQUESTED:
        r = lookup.get((sample, haplotype))
        if r is None:
            missing.append((sample, haplotype))
        else:
            selected.append(r)

    # Fail hard if any assembly is missing — no silent empty manifest
    if missing:
        print("FATAL: The following requested assemblies were NOT found "
              "in the official HPRC Release 2 index:", file=sys.stderr)
        for s, h in missing:
            print(f"  {s} ({h})", file=sys.stderr)
        print(f"\nIndex URL: {HPRC_INDEX_URL}", file=sys.stderr)
        sys.exit(1)

    # Write compact manifest
    os.makedirs("work/manifests", exist_ok=True)
    out_path = "work/manifests/hprc_selected.csv"

    fieldnames = [
        "sample_id", "haplotype", "assembly_name", "assembly_md5",
        "assembly_fai", "assembly_gzi", "assembly",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in selected:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    print(f"\nSelected {len(selected)}/{len(REQUESTED)} assemblies — all found.")
    print(f"Manifest: {out_path}")
    for r in selected:
        print(f"  {r['sample_id']} ({r['haplotype']}): {r['assembly_name']}")
        print(f"    S3: {r.get('assembly', 'N/A')}")

    # Also emit a human-readable summary
    print("\nNext step:")
    print("  python3 scripts/download_hprc.py")


if __name__ == "__main__":
    main()