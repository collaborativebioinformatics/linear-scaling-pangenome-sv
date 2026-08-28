#!/usr/bin/env python3
"""
sync_web_results.py — Sync latest pipeline results to web public dir.

TODO (Quang): Add support for syncing multiple run IDs.
TODO (Ali): Integrate with DNAnexus download step.

Usage:
    python3 scripts/sync_web_results.py
"""

import json
import os
import shutil


def main():
    # Only compact JSON files go to web/public/data/ - no GFA/FASTA/VCF.
    sources = [
        ("results/benchmark/report.json", "web/public/data/latest.json"),
        ("results/web/manifest.json", "web/public/data/manifest.json"),
        ("results/web/overview.json", "web/public/data/overview.json"),
    ]

    for src, dst in sources:
        if not os.path.exists(src):
            print(f"  SKIP {src} (not found)")
            continue
        if not guard_no_genomic(src):
            continue
        if not guard_file_size(src, int(os.environ.get("WEB_MAX_FILE_MB", "10"))):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        print(f"Synced {src} -> {dst}")

    # Sync the per-sample graph directory (results/web/graphs -> web/public/data/graphs)
    src_graphs = "results/web/graphs"
    dst_graphs = "web/public/data/graphs"
    if os.path.isdir(src_graphs):
        for root, _dirs, files in os.walk(src_graphs):
            for fn in files:
                if not fn.endswith(".json"):
                    continue
                fp = os.path.join(root, fn)
                if not guard_no_genomic(fp):
                    continue
                if not guard_file_size(fp, int(os.environ.get("WEB_MAX_FILE_MB", "10"))):
                    continue
                rel = os.path.relpath(fp, src_graphs)
                dst = os.path.join(dst_graphs, rel)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy2(fp, dst)
        print(f"Synced {src_graphs} -> {dst_graphs}")

    # Scan web/public for forbidden files
    for root, dirs, files in os.walk("web/public"):
        for fn in files:
            fp = os.path.join(root, fn)
            guard_no_genomic(fp)

    # Fallback check
    demo_json = "web/public/data/latest.json"
    if not os.path.exists(demo_json) or os.path.getsize(demo_json) < 10:
        print("WARNING: No pipeline results yet. Run `make demo` for synthetic data.")



def guard_no_genomic(path):
    """Refuse GFA/FASTA/VCF/BAM/CRAM files in web/public."""
    forbidden = [".gfa", ".fa", ".fasta", ".vcf", ".vcf.gz", ".fa.gz",
                 ".fasta.gz", ".gfa.gz", ".bam", ".cram"]
    ext = os.path.splitext(path)[1].lower()
    if ext in forbidden or any(path.endswith(e) for e in forbidden):
        print(f"  BLOCKED: {path} - genomic files not allowed in web/public/")
        return False
    return True

def guard_file_size(filepath, max_mb=10):
    """Refuse files larger than max_mb."""
    if os.path.exists(filepath) and os.path.getsize(filepath) > max_mb * 1024 * 1024:
        print(f"  SKIP {filepath}: {os.path.getsize(filepath)/1024/1024:.1f}MB > {max_mb}MB limit")
        return False
    return True

if __name__ == "__main__":
    main()