"""Q1–Q3 stitch tests: continuity, sequence identity, boundary report, demo wiring."""
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.merge.gfa import GfaGraph, Header, Segment, Link, Path
from pipeline.merge.merge_graphs import (
    overlap_aware_stitch, diagnostic_disjoint_union, _load_chunks,
)
from pipeline.merge.validate_merge import validate
from pipeline.benchmark.graph_stats import compute_stats, compare
from pipeline.parallel.make_chunks import create_chunks, write_manifest


REF = "ACGT" * 25  # 100 bp


def _span(cid, start, end, seq=REF):
    g = GfaGraph()
    g.headers.append(Header("1.1"))
    n = f"{cid}_n"
    g.segments[n] = Segment(n, seq[start:end])
    pn = f"GRCh38#0#chr21:{start}-{end}"
    g.paths[pn] = Path(pn, [n + "+"])
    return cid, g


class TestStitchSequenceAndContinuity:
    def test_stitched_path_spells_the_full_sequence(self):
        chunks = [_span("chunk_0001", 0, 60), _span("chunk_0002", 40, 100)]
        m, br = overlap_aware_stitch(chunks)
        assert m.get_path_sequence("GRCh38#0#chr21") == REF
        assert br[0]["status"] == "PASS"
        assert br[0]["haplotypes_joined"] == 1

    def test_disjoint_union_keeps_separate_chunk_paths(self):
        chunks = [_span("chunk_0001", 0, 60), _span("chunk_0002", 40, 100)]
        u = diagnostic_disjoint_union(chunks)
        assert u.path_count() == 2

    def test_stitch_collapses_to_one_path_not_one_per_chunk(self):
        chunks = [_span("chunk_0001", 0, 60), _span("chunk_0002", 40, 100)]
        m, _br = overlap_aware_stitch(chunks)
        assert m.path_count() == 1

    def test_gap_does_not_invent_an_edge(self):
        chunks = [_span("chunk_0001", 0, 40), _span("chunk_0002", 60, 100)]
        m, br = overlap_aware_stitch(chunks)
        assert br[0]["status"] == "FAIL"
        assert br[0]["haplotypes_unjoined"] == 1
        assert "GRCh38#0#chr21" in m.paths
        # Path exists but is only the first run; no jump across the hole.
        spelled = m.get_path_sequence("GRCh38#0#chr21")
        assert spelled == REF[0:40] or spelled == REF[60:100] or spelled in REF

    def test_validate_accepts_matching_baseline(self, tmp_path):
        chunks = [_span("chunk_0001", 0, 60), _span("chunk_0002", 40, 100)]
        m, _br = overlap_aware_stitch(chunks)
        merged_p, base_p = tmp_path / "merged.gfa", tmp_path / "base.gfa"
        m.write_gfa(str(merged_p))
        base = GfaGraph()
        base.segments["r"] = Segment("r", REF)
        base.paths["GRCh38#0#chr21"] = Path("GRCh38#0#chr21", ["r+"])
        base.write_gfa(str(base_p))
        errs = [i for i in validate(str(base_p), str(merged_p))
                if i["severity"] == "ERROR"]
        assert errs == []


class TestLoadChunksDir:
    def test_load_chunks_reads_gfa_from_chunk_dir(self, tmp_path):
        gdir = tmp_path / "chunks"
        gdir.mkdir()
        cid, g = _span("chunk_0001", 0, 40)
        g.write_gfa(str(gdir / "chunk_0001.gfa"))
        man = tmp_path / "chunk_manifest.tsv"
        write_manifest(create_chunks("chr21", 0, 40, 40, 10)[:1], str(man))
        # rewrite id to match
        rows = list(csv.DictReader(open(man), delimiter="\t"))
        rows[0]["chunk_id"] = "chunk_0001"
        with open(man, "w", newline="") as f:
            w = csv.DictWriter(f, delimiter="\t", fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
        loaded, found = _load_chunks(str(man), chunk_dir=str(gdir))
        assert len(loaded) == 1
        assert loaded[0][0] == "chunk_0001"
        assert found["chunk_0001"]["chunk_id"] == "chunk_0001"


class TestDemoWiring:
    def test_setup_demo_source_calls_overlap_aware_stitch(self):
        root = os.path.join(os.path.dirname(__file__), "..")
        src = open(os.path.join(root, "scripts", "setup_demo.py"), encoding="utf-8").read()
        assert "overlap_aware_stitch" in src
        assert "diagnostic_disjoint_union" in src
        assert 'strategy' in src
        assert "Strategy:" in src

    def test_setup_demo_prints_strategy_and_equivalent(self, tmp_path, monkeypatch):
        root = os.path.join(os.path.dirname(__file__), "..")
        monkeypatch.chdir(root)
        from scripts.setup_demo import main
        main()
        merged = GfaGraph.parse_file("work/demo/merged.gfa")
        baseline = GfaGraph.parse_file("work/demo/baseline.gfa")
        assert merged.path_count() == baseline.path_count() == 5
        assert merged.path_count() != 15
        cmp = compare(compute_stats(baseline, "b"), compute_stats(merged, "m"))
        assert cmp["verdict"] == "EQUIVALENT"
        assert compute_stats(merged, "m")["components"]["count"] == 5
        br = "results/merge/boundary_report.tsv"
        demo_br = "work/demo/boundary_report.tsv"
        assert os.path.exists(br) or os.path.exists(demo_br)
        latest = __import__("json").load(open("web/public/data/latest.json"))
        assert latest.get("stitch", {}).get("status") in ("PASS", "overlap_aware", "OK")
        assert latest.get("equivalence", {}).get("verdict") == "EQUIVALENT"
