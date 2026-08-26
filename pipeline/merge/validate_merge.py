"""
validate_merge.py — Validate merged graph properties.

Checks are independent, so a merged graph produced without a baseline still
gets everything except the comparison checks. The load-bearing one is
`_check_sequences`: a merged graph whose paths spell exactly what the
whole-region build spells is correct regardless of how its nodes ended up
numbered.
"""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def _issue(severity, message):
    return {"severity": severity, "message": message}


def _check_structure(merged):
    issues = []
    if merged.node_count() == 0:
        issues.append(_issue("ERROR", "Merged graph has zero nodes"))
    if not merged.paths and not merged.walks:
        issues.append(_issue("WARNING", "No paths/walks in merged graph"))
    dangling = merged.dangling_links()
    if dangling:
        ex = ", ".join(f"{l.from_node}->{l.to_node}" for l in dangling[:3])
        issues.append(_issue(
            "ERROR", f"{len(dangling)} link(s) point at missing segments: {ex}"))
    orphans = merged.orphan_segments()
    if orphans:
        ex = ", ".join(sorted(orphans)[:3])
        issues.append(_issue(
            "WARNING", f"{len(orphans)} orphaned segment(s) on no path: {ex}"))
    return issues


def _check_path_connectivity(merged):
    """Every consecutive pair of steps on a path must have a link behind it.

    A stitch that concatenated two chunks without joining them shows up here as
    a break, even though the path itself still spells the right sequence.
    """
    have = {(l.from_node, l.from_orient, l.to_node, l.to_orient)
            for l in merged.links}
    have |= {(l.to_node, "-" if l.to_orient == "+" else "+",
              l.from_node, "-" if l.from_orient == "+" else "+")
             for l in merged.links}
    issues = []
    for pn in sorted(merged.paths):
        steps = merged.path_steps(pn)
        breaks = [(a, b) for a, b in zip(steps, steps[1:])
                  if (a[0], a[1], b[0], b[1]) not in have]
        if breaks:
            a, b = breaks[0]
            issues.append(_issue(
                "ERROR", f"path {pn} has {len(breaks)} unlinked step(s), "
                         f"first at {a[0]}{a[1]}->{b[0]}{b[1]}"))
    return issues


def _check_reference(merged, ref_name="GRCh38"):
    refs = [pn for pn in merged.paths if pn.startswith(ref_name)]
    if not refs:
        return [_issue("ERROR", f"No {ref_name} reference path in merged graph")]
    if len(refs) > 1:
        return [_issue("WARNING",
                       f"{len(refs)} reference paths, expected 1: {sorted(refs)}")]
    return []


def _check_samples(baseline, merged):
    missing = baseline.get_sample_names() - merged.get_sample_names()
    if missing:
        return [_issue("WARNING", f"Samples lost in merge: {sorted(missing)}")]
    return []


def _check_sequences(baseline, merged):
    """Each baseline path must spell identically in the merged graph."""
    issues = []
    for pn in sorted(baseline.paths):
        if pn not in merged.paths:
            issues.append(_issue("ERROR", f"path {pn} missing from merged graph"))
            continue
        want, got = baseline.get_path_sequence(pn), merged.get_path_sequence(pn)
        if want != got:
            where = next((i for i, (a, b) in enumerate(zip(want, got)) if a != b),
                         min(len(want), len(got)))
            issues.append(_issue(
                "ERROR", f"path {pn} sequence differs: baseline {len(want)}bp vs "
                         f"merged {len(got)}bp, first difference at {where}"))
    extra = set(merged.paths) - set(baseline.paths)
    if extra:
        issues.append(_issue("WARNING",
                             f"paths only in merged: {sorted(extra)[:3]}"))
    return issues


def validate(baseline_path: str, merged_path: str, ref_name="GRCh38") -> list:
    if not os.path.exists(merged_path):
        return [_issue("ERROR", "Merged graph not found")]

    merged = GfaGraph.parse_file(merged_path)
    issues = _check_structure(merged)
    issues += _check_path_connectivity(merged)
    issues += _check_reference(merged, ref_name)

    if os.path.exists(baseline_path):
        baseline = GfaGraph.parse_file(baseline_path)
        issues += _check_samples(baseline, merged)
        issues += _check_sequences(baseline, merged)
    else:
        issues.append(_issue(
            "INFO", "No baseline graph: sequence equivalence not checked"))
    return issues


def main():
    issues = validate("results/baseline/baseline.gfa", "results/merge/merged.gfa")
    errors = [i for i in issues if i["severity"] == "ERROR"]
    if not issues:
        print("Validation: ALL CHECKS PASSED")
    else:
        for i in issues:
            print(f"[{i['severity']}] {i['message']}")
    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
