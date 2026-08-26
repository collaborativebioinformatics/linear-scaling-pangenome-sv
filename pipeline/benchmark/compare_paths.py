"""Compare paths between baseline and merged graphs using sequence-level comparison.

Supports BOTH P (Path) and W (Walk) records.
Spells each graph path/walk sequence and calculates:
    sequence length, SHA256, exact sequence equality, base-level identity.

Does NOT call path-name presence "exact matching".
"""
import hashlib, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline.merge.gfa import GfaGraph


def spell_path_sequence(g, path_name):
    """Reconstruct the nucleotide sequence for a path or walk by
    following segment references through the graph."""
    seq = []
    if path_name in g.paths:
        for sn in g.paths[path_name].segment_names:
            name, orient = (sn[:-1], sn[-1]) if sn[-1] in "+-" else (sn, "+")
            if name in g.segments:
                s = g.segments[name].sequence
                seq.append(s if orient == "+" else _revcomp(s))
    for w in g.walks:
        wk = f"WALK:{w.sample}#{w.haplotype}#{w.contig}"
        if path_name == wk or path_name == w.sample or path_name == f"{w.sample}#{w.haplotype}":
            for step in w.path:
                name, orient = (step[:-1], step[-1]) if step[-1] in "+-" else (step, "+")
                if name in g.segments:
                    s = g.segments[name].sequence
                    seq.append(s if orient == "+" else _revcomp(s))
    return "".join(seq)


def _revcomp(s):
    t = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(t)[::-1]


def main():
    rd = "results"
    bp, mp = f"{rd}/baseline/baseline.gfa", f"{rd}/merge/merged.gfa"
    if not (os.path.exists(bp) and os.path.exists(mp)):
        print("Both graphs required.")
        return
    b, m = GfaGraph.parse_file(bp), GfaGraph.parse_file(mp)
    T = chr(9)

    # Collect all path/walk names from both graphs
    all_names = sorted(set(b.paths) | set(m.paths))
    for w in b.walks:
        all_names.append(f"WALK:{w.sample}#{w.haplotype}#{w.contig}")
    for w in m.walks:
        all_names.append(f"WALK:{w.sample}#{w.haplotype}#{w.contig}")
    all_names = sorted(set(all_names))

    rows = []
    for pn in all_names:
        bs = spell_path_sequence(b, pn)
        ms = spell_path_sequence(m, pn)

        blen = len(bs)
        mlen = len(ms)
        bsha = hashlib.sha256(bs.encode()).hexdigest()[:16] if bs else ""
        msha = hashlib.sha256(ms.encode()).hexdigest()[:16] if ms else ""
        exact = "true" if bs == ms else "false"

        if bs and ms:
            mismatches = sum(1 for a, b in zip(bs, ms) if a != b)
            tail_diff = abs(blen - mlen)
            diffs = mismatches + tail_diff
            max_len = max(blen, mlen)
            identity = f"{(max_len - diffs) / max_len:.6f}" if max_len > 0 else "1.0"
        else:
            diffs = 0
            identity = "N/A"

        sm = pn.replace("WALK:", "").split("#")[0] if "#" in pn else pn
        st = "OK"
        if not bs and not ms:
            st = "BOTH_EMPTY"
        elif not bs:
            st = "MISSING_IN_BASELINE"
        elif not ms:
            st = "MISSING_IN_MERGED"
        elif bs != ms:
            st = "DIFFERENT"

        rows.append(f"{sm}{T}{pn}{T}{blen}{T}{mlen}{T}{bsha}{T}{msha}{T}"
                    f"{exact}{T}{identity}{T}{diffs}{T}{st}")

    op = f"{rd}/benchmark/path_comparison.tsv"
    os.makedirs(os.path.dirname(op), exist_ok=True)
    with open(op, "w") as f:
        f.write(f"sample{T}identifier{T}baseline_len{T}merged_len{T}"
                f"baseline_sha256{T}merged_sha256{T}"
                f"exact_match{T}identity{T}diffs{T}status\n")
        f.write("\n".join(rows) + "\n")
    print(f"{len(rows)} paths -> {op}")


if __name__ == "__main__":
    main()
