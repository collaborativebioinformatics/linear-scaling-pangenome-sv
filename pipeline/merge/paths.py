"""
paths.py — Path utilities for merged pangenome graphs.
TODO (Quang): Add path traversal across chunk boundaries.
TODO (Michael): Add sequence-level validation of merged paths.
"""
from pipeline.merge.gfa import GfaGraph


def compare_path_lengths(baseline: GfaGraph, merged: GfaGraph,
                         path_name: str) -> dict:
    bp = baseline.paths.get(path_name)
    mp = merged.paths.get(path_name)
    return {
        "path": path_name,
        "baseline_segments": len(bp.segment_names) if bp else 0,
        "merged_segments": len(mp.segment_names) if mp else 0,
        "in_baseline": bp is not None,
        "in_merged": mp is not None,
    }