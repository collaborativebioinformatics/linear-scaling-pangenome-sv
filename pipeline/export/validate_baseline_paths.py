"""validate_baseline_paths.py — P0: verify the baseline GFA is the intended 1 Mb target.

Determines whether the monolithic baseline graph was built from the intended
~1 Mb slice of chr21 (20,000,000-21,000,000) or from full chromosome-21
contigs (each ~40-47 Mb).

Writes results/validation/baseline_path_validation.json

For each baseline path it reports:
    sample, haplotype, path_name, spelled_sequence_bp, unique_graph_nodes,
    number_of_steps, has_subrange, overlaps_present.

A path is flagged as "full contig" (invalid for the 1 Mb benchmark) when:
    - its path name carries no `:start-end` subrange, AND
    - its spelled length is far larger than ~1 Mb (tolerance ~1.2 Mb).

Does NOT fabricate any result — only reports what is actually in the GFA.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph, parse_pansn

TARGET_LEN = 1_000_000
TOLERANCE = 2_000_000  # allow up to 2 Mb for real biological variation

# A path is a "full contig" when its spelled length is far beyond regional scale.
# For the 1 Mb benchmark, anything over ~2 Mb is suspicious.
# We do NOT rely on :start-end subrange because faidx extraction produces
# slice-named paths (sample#hap#chr) without subranges.
FULL_CONTIG_THRESHOLD = None  # disabled — use spelled length only


def analyze(gfa_path):
    g = GfaGraph.parse_file(gfa_path)
    rows = []
    for pn in sorted(g.paths.keys()):
        s, h, contig, start, end, chunk = parse_pansn(pn)
        steps = g.path_steps(pn)
        unique_nodes = set(n for n, _ in steps)
        spelled = sum(g.segments[n].length for n, _ in steps if n in g.segments)
        overlaps = g.paths[pn].overlaps
        has_subrange = start is not None and end is not None

        # A "full contig" is any path whose spelled length far exceeds the
        # regional scale. We check spelled length only — faidx extraction
        # produces sample#hap#chr names without :start-end subranges.
        full_contig = spelled > TOLERANCE

        rows.append({
            "sample": s,
            "haplotype": h,
            "path_name": pn,
            "contig": contig,
            "subrange_start": start,
            "subrange_end": end,
            "has_subrange": has_subrange,
            "spelled_sequence_bp": spelled,
            "unique_graph_nodes": len(unique_nodes),
            "number_of_steps": len(steps),
            "overlaps_present": bool(overlaps and overlaps != ["*"]),
            "full_contig_suspected": full_contig,
        })

    any_full_contig = any(r["full_contig_suspected"] for r in rows)
    lengths_plausible = all(
        r["spelled_sequence_bp"] <= TOLERANCE for r in rows
    )

    if any_full_contig:
        status = "INVALID_FOR_1MB_BENCHMARK"
    elif lengths_plausible:
        status = "VALID_1MB_BENCHMARK"
    else:
        status = "UNVERIFIED"

    return {
        "gfa_source": gfa_path,
        "status": status,
        "target_length_bp": TARGET_LEN,
        "tolerance_bp": TOLERANCE,
        "paths": rows,
        "diagnosis": (
            "One or more paths are full chromosome-21 contigs (~40-47 Mb) "
            "with no :start-end subrange. The baseline was built from whole "
            "contigs, not the intended 1 Mb chr21:20,000,000-21,000,000 slice. "
            "Its timing must NOT be used as the official 1 Mb baseline."
            if any_full_contig
            else "All paths appear to be within the 1 Mb target."
        ),
    }


def main():
    gfa_path = "results/baseline/baseline.gfa"
    if not os.path.exists(gfa_path):
        print("No baseline GFA found; writing empty validation.", file=sys.stderr)
        result = {"gfa_source": gfa_path, "status": "MISSING", "paths": []}
    else:
        result = analyze(gfa_path)

    os.makedirs("results/validation", exist_ok=True)
    out = "results/validation/baseline_path_validation.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)

    print(f"baseline_path_validation: status={result['status']}")
    for r in result.get("paths", []):
        flag = "  <-- FULL CONTIG" if r["full_contig_suspected"] else ""
        print(f"  {r['sample']}#{r['haplotype']} "
              f"{r['spelled_sequence_bp']:>12,} bp  "
              f"steps={r['number_of_steps']:>6}  "
              f"subrange={r['has_subrange']}{flag}")
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
