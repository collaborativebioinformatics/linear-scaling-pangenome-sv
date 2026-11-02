#!/usr/bin/env python3
"""Stitch the same chunk graphs twice and write both results side by side.

  results/merge_midpoint/merged.gfa   original cut (each haplotype at the
                                      midpoint of its own overlap)
  results/merge_anchor/merged.gfa     all haplotypes cut at one shared node

Run from a layout directory, the one that holds work/chunks/ and
results/preparation/:

    python3 scripts/merge_ab.py [--overlap-bp N] [--min-anchor-bp N]

Both runs read the same chunk graphs, so no PGGB run is repeated. The midpoint
output should reproduce the merged.gfa you already have; if it does not, stop
and find out why before trusting the anchor output.
"""
import argparse
import hashlib
import os
import sys
from collections import Counter

sys.path.insert(0, os.getcwd())
from pipeline.merge.merge_graphs import (  # noqa: E402
    _load_chunks, _load_chunk_mapping, _write_boundary_report,
    overlap_aware_stitch)


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()[:12]


def _overlap_from_config(default=100000):
    try:
        import yaml
        with open("config/pipeline.yaml") as f:
            return int(yaml.safe_load(f).get("parallel", {}).get("overlap_bp", default))
    except Exception:
        return default


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--manifest", default="work/chunks/chunk_manifest.tsv")
    ap.add_argument("--mapping", default="results/preparation/chunk_mapping.tsv")
    ap.add_argument("--overlap-bp", type=int, default=None)
    ap.add_argument("--min-anchor-bp", type=int, default=1)
    ap.add_argument("--out", default="results")
    a = ap.parse_args()

    chunks, rows = _load_chunks(a.manifest)
    if not chunks:
        sys.exit(f"no chunk graphs found next to {a.manifest}")
    mapping = _load_chunk_mapping(a.mapping)
    obp = a.overlap_bp if a.overlap_bp is not None else _overlap_from_config()
    print(f"{len(chunks)} chunks, {len(mapping)} chunk_mapping rows, "
          f"overlap_bp={obp}, min_anchor_bp={a.min_anchor_bp}")

    for mode in ("midpoint", "anchor"):
        merged, br = overlap_aware_stitch(
            chunks, overlap_bp=obp, chunk_rows=rows, chunk_mapping=mapping,
            seam_mode=mode, min_anchor_bp=a.min_anchor_bp)
        d = os.path.join(a.out, f"merge_{mode}")
        os.makedirs(d, exist_ok=True)
        gfa = os.path.join(d, "merged.gfa")
        merged.write_gfa(gfa)
        _write_boundary_report(br, os.path.join(d, "boundary_report.tsv"))
        status = Counter(b["status"] for b in br)
        seams = Counter(b["seam_mode"] for b in br)
        print(f"[{mode}] nodes={merged.node_count()} edges={merged.edge_count()} "
              f"bp={merged.total_sequence_bp()} sha256={_sha(gfa)}")
        print(f"[{mode}] boundaries={len(br)} status={dict(status)} "
              f"seam_mode={dict(seams)}")
        print(f"[{mode}] wrote {gfa}")


if __name__ == "__main__":
    main()
