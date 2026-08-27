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
    # Only compact JSON files go to web/public/data/ - no GFA
    sources = [
        ("results/benchmark/report.json", "web/public/data/latest.json"),
    ]

    for src, dst in sources:
        if os.path.exists(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)
            print(f"Synced {src} -> {dst}")
        else:
            print(f"  SKIP {src} (not found)")

    # Copy the demo JSON as fallback if real results don't exist
    demo_json = "web/public/data/latest.json"
    if not os.path.exists(demo_json) or os.path.getsize(demo_json) < 10:
        print("WARNING: No pipeline results yet. Run `make demo` for synthetic data.")


if __name__ == "__main__":
    main()