#!/usr/bin/env python3
"""
download_hprc.py — Download HPRC assemblies from the manifest using S3 URIs.

REAL-DATA STEP: Each assembly is ~3GB. Requires network access to HPRC public S3.

Usage:
    python3 scripts/download_hprc.py              # print download commands
    python3 scripts/download_hprc.py --execute     # actually download (requires aws CLI)

The column 'assembly' in the manifest contains the exact S3 URI from the
official HPRC Release 2 index (e.g.
  s3://human-pangenomics/working/HPRC/HG00673/...)
"""

import argparse
import csv
import os
import subprocess
import sys


def load_manifest(path: str = "work/manifests/hprc_selected.csv") -> list[dict]:
    if not os.path.exists(path):
        print(f"Manifest not found: {path}", file=sys.stderr)
        print("Run: python3 scripts/fetch_hprc_index.py", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        print(f"Manifest is empty: {path}", file=sys.stderr)
        sys.exit(1)
    return rows


def main():
    parser = argparse.ArgumentParser(description="Download HPRC assemblies")
    parser.add_argument("--execute", action="store_true",
                        help="Actually run aws s3 cp commands")
    parser.add_argument("--dest", default="work/downloads",
                        help="Destination directory (default: work/downloads)")
    args = parser.parse_args()

    rows = load_manifest()
    os.makedirs(args.dest, exist_ok=True)

    print(f"Found {len(rows)} assemblies in manifest\n")

    for r in rows:
        sample = r.get("sample_id", "?")
        hap = r.get("haplotype", "?")
        name = r.get("assembly_name", f"{sample}_{hap}")
        uri = r.get("assembly", "").strip()
        dest = os.path.join(args.dest, f"{name}.fa")

        if not uri:
            print(f"  SKIP {name}: no S3 URI in manifest (column 'assembly')")
            continue

        if os.path.exists(dest):
            size_mb = os.path.getsize(dest) / (1024 * 1024)
            print(f"  EXISTS {name}: {dest} ({size_mb:.0f} MB)")
            continue

        # Build the aws cp command
        cmd = ["aws", "s3", "cp", "--no-sign-request", uri, dest]
        cmd_str = " ".join(cmd)

        if args.execute:
            print(f"  DOWNLOAD {name} ({sample} {hap})...")
            sys.stdout.flush()
            try:
                subprocess.run(cmd, check=True)
                size_mb = os.path.getsize(dest) / (1024 * 1024)
                print(f"    -> {dest} ({size_mb:.0f} MB)")
            except subprocess.CalledProcessError as e:
                print(f"    FAILED: {e}", file=sys.stderr)
            except FileNotFoundError:
                print("    FAILED: 'aws' CLI not found. Install with: "
                      "pip install awscli", file=sys.stderr)
                sys.exit(1)
        else:
            print(f"  QUEUED {name} ({sample} {hap})")
            print(f"    URI: {uri}")
            print(f"    ->   {dest}")
            print(f"    cmd: {cmd_str}")
            print()

    if not args.execute:
        print(f"\nTo actually download, run: {sys.argv[0]} --execute")
        print("Or pipe individual commands:")
        for r in rows:
            uri = r.get("assembly", "").strip()
            name = r.get("assembly_name", f"{r.get('sample_id','?')}_{r.get('haplotype','?')}")
            if uri:
                print(f"  aws s3 cp --no-sign-request {uri} {args.dest}/{name}.fa")


if __name__ == "__main__":
    main()