#!/usr/bin/env python3
"""
fetch_hprc_index.py — Select the 4 requested HPRC assemblies by exact
assembly_name from the official Release 2 index, and write a compact
manifest with the numeric haplotype column preserved.

The official index uses NUMERIC haplotype values (1, 2).
Do NOT key on "maternal"/"paternal" strings — the official CSV does not
use them.  Derive human-readable labels from the assembly_name pattern.

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

# Canonical assembly names from config/samples.yaml.
# These are the exact strings in the official "assembly_name" column.
REQUESTED_ASSEMBLY_NAMES = [
    "HG00673_mat_hprc_r2_v1.0.1",
    "HG00673_pat_hprc_r2_v1.0.1",
    "HG00733_mat_hprc_r2_v1.0.1",
    "HG00733_pat_hprc_r2_v1.0.1",
]


def _haplotype_label(assembly_name: str) -> str:
    """Derive human-readable label from the canonical assembly name."""
    if "_mat_" in assembly_name:
        return "maternal"
    if "_pat_" in assembly_name:
        return "paternal"
    return "unknown"


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

    # Build a lookup keyed by assembly_name (exact match).
    # The official "assembly_name" column contains strings like
    # "HG00673_mat_hprc_r2_v1.0.1".
    lookup: dict[str, dict] = {}
    for r in rows:
        aname = r.get("assembly_name", "").strip()
        if aname:
            if aname in lookup:
                print(f"FATAL: Duplicate assembly_name in index: {aname}",
                      file=sys.stderr)
                sys.exit(1)
            lookup[aname] = r

    # Select the 4 we need by exact assembly_name match.
    selected: list[dict] = []
    missing: list[str] = []
    for aname in REQUESTED_ASSEMBLY_NAMES:
        r = lookup.get(aname)
        if r is None:
            missing.append(aname)
        else:
            selected.append(r)

    # Fail hard if any assembly is missing.
    if missing:
        print("FATAL: The following requested assemblies were NOT found "
              "in the official HPRC Release 2 index:", file=sys.stderr)
        for a in missing:
            print(f"  {a}", file=sys.stderr)
        print(f"\nIndex URL: {HPRC_INDEX_URL}", file=sys.stderr)
        sys.exit(1)

    # Fail on duplicates.
    if len(selected) != len(REQUESTED_ASSEMBLY_NAMES):
        print(f"FATAL: Expected {len(REQUESTED_ASSEMBLY_NAMES)} records, "
              f"got {len(selected)}. Possible duplicate.", file=sys.stderr)
        sys.exit(1)

    # Write compact manifest with both numeric haplotype and human-readable label.
    os.makedirs("work/manifests", exist_ok=True)
    out_path = "work/manifests/hprc_selected.csv"

    fieldnames = [
        "sample_id", "haplotype", "haplotype_label", "assembly_name",
        "assembly_md5", "assembly_fai", "assembly_gzi", "assembly",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in selected:
            aname = r.get("assembly_name", "")
            row_out = {
                "sample_id": r.get("sample_id", ""),
                "haplotype": r.get("haplotype", ""),
                "haplotype_label": _haplotype_label(aname),
                "assembly_name": aname,
                "assembly_md5": r.get("assembly_md5", ""),
                "assembly_fai": r.get("assembly_fai", ""),
                "assembly_gzi": r.get("assembly_gzi", ""),
                "assembly": r.get("assembly", ""),
            }
            w.writerow(row_out)

    print(f"\nSelected {len(selected)}/{len(REQUESTED_ASSEMBLY_NAMES)} assemblies — all found.")
    print(f"Manifest: {out_path}")
    for r in selected:
        aname = r.get("assembly_name", "")
        label = _haplotype_label(aname)
        sample = r.get("sample_id", "")
        numeric_hap = r.get("haplotype", "")
        print(f"  {sample} haplotype={numeric_hap} ({label}): {aname}")
        print(f"    S3: {r.get('assembly', 'N/A')}")

    print("\nNext step:")
    print("  python3 scripts/download_hprc.py")


if __name__ == "__main__":
    main()