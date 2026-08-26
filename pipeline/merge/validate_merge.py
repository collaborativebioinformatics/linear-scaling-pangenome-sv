"""
validate_merge.py — Validate merged graph properties.
TODO (Quang): Check all haplotypes preserved, reference continuity, no orphaned nodes.
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def validate(baseline_path: str, merged_path: str) -> list:
    issues = []
    if not os.path.exists(merged_path):
        issues.append({"severity": "ERROR", "message": "Merged graph not found"})
        return issues

    merged = GfaGraph.parse_file(merged_path)
    if merged.node_count() == 0:
        issues.append({"severity": "ERROR", "message": "Merged graph has zero nodes"})
    if not merged.paths and not merged.walks:
        issues.append({"severity": "WARNING", "message": "No paths/walks in merged graph"})

    if os.path.exists(baseline_path):
        baseline = GfaGraph.parse_file(baseline_path)
        missing = baseline.get_sample_names() - merged.get_sample_names()
        if missing:
            issues.append({"severity": "WARNING",
                           "message": f"Samples lost in merge: {missing}"})

    return issues


def main():
    issues = validate("results/baseline/baseline.gfa", "results/merge/merged.gfa")
    if not issues:
        print("Validation: ALL CHECKS PASSED")
    else:
        for i in issues:
            print(f"[{i['severity']}] {i['message']}")


if __name__ == "__main__":
    main()