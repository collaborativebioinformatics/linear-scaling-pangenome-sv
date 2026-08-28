"""Test static results architecture + P0 baseline validation."""
import json
import os
import sys
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.export.validate_baseline_paths import analyze, TARGET_LEN, TOLERANCE
from pipeline.export import freeze_web_results


def _full_contig_gfa(tmp_path):
    """A GFA with full-contig paths (no :start-end subrange, >1Mb)."""
    gfa = tmp_path / "full.gfa"
    gfa.write_text(
        "H\tVN:Z:1.0\n"
        "S\t1\t" + "A" * 46_700_000 + "\n"
        "S\t2\t" + "C" * 40_000_000 + "\n"
        "P\tGRCh38#0#chr21\t1+\t*\n"
        "P\tHG00673#1#JAHBBZ020000061.1\t2+\t*\n"
    )
    return str(gfa)


def _valid_1mb_gfa(tmp_path):
    """A GFA with 1Mb subrange paths."""
    gfa = tmp_path / "valid.gfa"
    seq = "A" * 500_000
    gfa.write_text(
        "H\tVN:Z:1.0\n"
        "S\t1\t" + seq + "\n"
        "S\t2\t" + "C" * 500_000 + "\n"
        "P\tGRCh38#0#chr21:20000000-21000000\t1+\t*\n"
        "P\tHG00673#1#JAHBBZ020000061.1:20000000-21000000\t2+\t*\n"
    )
    return str(gfa)


def test_full_contig_flagged_invalid(tmp_path):
    r = analyze(_full_contig_gfa(tmp_path))
    assert r["status"] == "INVALID_FOR_1MB_BENCHMARK"
    assert all(p["full_contig_suspected"] for p in r["paths"])
    assert any(p["spelled_sequence_bp"] > TOLERANCE for p in r["paths"])


def test_1mb_subrange_is_valid(tmp_path):
    r = analyze(_valid_1mb_gfa(tmp_path))
    assert r["status"] == "VALID_1MB_BENCHMARK"
    assert all(p["spelled_sequence_bp"] <= TOLERANCE for p in r["paths"])


def test_placeholder_stitch_never_real_validated(monkeypatch, tmp_path):
    """freeze_web_results must never label a linear placeholder as real."""
    # Force the merged GFA to be a linear placeholder and no baseline.
    monkeypatch.chdir(tmp_path)
    os.makedirs("results/merge", exist_ok=True)
    with open("results/merge/merged.gfa", "w") as f:
        f.write(
            "H\tVN:Z:1.0\n"
            "S\t1\tA\nS\t2\tC\nS\t3\tG\n"
            "L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n"
            "P\tGRCh38#0#chr21\t1+,2+,3+\t*\n"
        )
    # No baseline validation -> NOT_RUN baseline
    files = freeze_web_results._validate_and_build()
    v = files["manifest.json"]["validation"]
    assert v["stitch"] != "REAL_HPRC_VALIDATED"
    assert v["stitch"] in ("PLACEHOLDER_LINEAR_ONLY", "PRESENT_BUT_NOT_VALIDATED")
    assert v["variant_equivalence"] == "NOT_RUN"


def test_freeze_never_writes_genomic_to_web(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    os.makedirs("web/public/data/final", exist_ok=True)
    os.makedirs("results/final_run", exist_ok=True)
    # Create a fake genomic file that must NOT be copied
    with open("results/final_run/leak.gfa", "w") as f:
        f.write("S\t1\tAAAA\n")
    assert not freeze_web_results._guard_no_genomic("leak.gfa")
    assert freeze_web_results._guard_no_genomic("manifest.json")


def test_visualization_selector_deterministic(tmp_path):
    """select_visualization_region returns identical results on repeat calls."""
    from pipeline.export.select_visualization_region import select_region
    gfa = tmp_path / "g.gfa"
    # A small graph with one clear branch node
    gfa.write_text(
        "H\tVN:Z:1.0\n"
        "S\t1\t" + "A" * 100 + "\n"
        "S\t2\t" + "C" * 100 + "\n"
        "S\t3\t" + "G" * 100 + "\n"
        "S\t4\t" + "T" * 100 + "\n"
        "L\t1\t+\t2\t+\t0M\n"
        "L\t1\t+\t3\t+\t0M\n"
        "L\t2\t+\t4\t+\t0M\n"
        "L\t3\t+\t4\t+\t0M\n"
        "P\tGRCh38#0#chr21:0-1000\t1+,2+,4+\t*\n"
        "P\tHG00673#1#c:0-1000\t1+,3+,4+\t*\n"
    )
    r1 = select_region(str(gfa), window_bp=50)
    r2 = select_region(str(gfa), window_bp=50)
    assert r1["score"] == r2["score"]
    assert r1["node_start"] == r2["node_start"]


def test_target_length_is_1mb(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # No config file -> default 1 Mb
    t = freeze_web_results._load_config_target()
    assert t["start"] == 20_000_000
    assert t["end"] == 21_000_000
    assert t["end"] - t["start"] == 1_000_000
