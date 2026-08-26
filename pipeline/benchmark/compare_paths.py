"""Compare paths between baseline and merged graphs."""
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def main():
    rd = "results"
    bp, mp = f"{rd}/baseline/baseline.gfa", f"{rd}/merge/merged.gfa"
    if not (os.path.exists(bp) and os.path.exists(mp)):
        print("Both graphs required.")
        return
    b, m = GfaGraph.parse_file(bp), GfaGraph.parse_file(mp)
    all_p = sorted(set(b.paths) | set(m.paths))
    T, N = chr(9), chr(10)
    rows = []
    for pn in all_p:
        bs = len(b.paths[pn].segment_names) if pn in b.paths else 0
        ms = len(m.paths[pn].segment_names) if pn in m.paths else 0
        ok = "true" if (pn in b.paths and pn in m.paths) else "false"
        st = "OK" if ok == "true" else ("MISSING" if pn in b.paths else "EXTRA")
        sm = pn.split("#")[0] if "#" in pn else pn
        rows.append(f"{sm}{T}{pn}{T}{bs}{T}{ms}{T}{ok}{T}{st}")
    op = f"{rd}/benchmark/path_comparison.tsv"
    os.makedirs(os.path.dirname(op), exist_ok=True)
    with open(op, "w") as f:
        f.write(f"sample{T}haplotype{T}baseline_segments{T}")
        f.write(f"merged_segments{T}exact_match{T}status{N}")
        f.write(N.join(rows) + N)
    print(f"{len(rows)} paths -> {op}")


if __name__ == "__main__":
    main()