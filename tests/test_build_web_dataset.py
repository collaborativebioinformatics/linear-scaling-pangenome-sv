"""Behavioral tests for pipeline/export/build_web_dataset.py.

Verifies: sample discovery, bounded export caps, truncation metadata, and
that no genomic sequence is ever written into a web JSON file.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pipeline", "export"))

from pipeline.merge.gfa import GfaGraph
import build_web_dataset as bwd

REPO = os.path.join(os.path.dirname(__file__), "..")

LINEAR_GFA = """H\tVN:Z:1.0
S\ts1\tACGTACGTAC
S\ts2\tGGGGTTTTAA
S\ts3\tCCCCAAAATT
L\ts1\t+\ts2\t+\t0M
L\ts2\t+\ts3\t+\t0M
P\tGRCh38#0#chr21\ts1+,s2+,s3+\t*
W\tHG00673\t1\tchr21\t0\t30\t>s1>s2>s3
"""


def test_samples_from_graph_p_and_w():
    g = GfaGraph.parse(LINEAR_GFA)
    samples = bwd._samples_from_graph(g)
    assert "GRCh38" in samples
    assert samples["GRCh38"] == {"0"}
    assert "HG00673" in samples
    assert samples["HG00673"] == {"1"}


def test_build_selected_path_p_line():
    g = GfaGraph.parse(LINEAR_GFA)
    pn, steps = bwd._build_selected_path(g, "GRCh38", "0")
    assert pn == "GRCh38#0#chr21"
    assert [s for s, _ in steps] == ["s1", "s2", "s3"]


def test_build_selected_path_w_line():
    g = GfaGraph.parse(LINEAR_GFA)
    pn, steps = bwd._build_selected_path(g, "HG00673", "1")
    assert pn is not None
    assert [s for s, _ in steps] == ["s1", "s2", "s3"]


def test_extract_sample_graph_no_sequence_leak():
    g = GfaGraph.parse(LINEAR_GFA)
    out = bwd._extract_sample_graph(g, "GRCh38", "0", "baseline")
    blob = json.dumps(out)
    for seq in ["ACGTACGTAC", "GGGGTTTTAA", "CCCCAAAATT"]:
        assert seq not in blob, f"sequence leaked into web JSON: {seq}"


def test_extract_sample_graph_structure():
    g = GfaGraph.parse(LINEAR_GFA)
    out = bwd._extract_sample_graph(g, "HG00673", "1", "merged")
    assert out["sample"] == "HG00673"
    assert out["haplotype"] == "1"
    assert len(out["nodes"]) == 3
    assert len(out["edges"]) == 2
    assert out["path"]["length_bp"] == 30
    assert out["truncated"] is False


def test_node_cap_truncates():
    g = GfaGraph.parse(LINEAR_GFA)
    old = bwd.WEB_GRAPH_MAX_NODES
    try:
        bwd.WEB_GRAPH_MAX_NODES = 2
        out = bwd._extract_sample_graph(g, "HG00673", "1", "merged")
        assert out["truncated"] is True
        assert out["original_counts"]["nodes"] == 3
        assert len(out["nodes"]) == 2
    finally:
        bwd.WEB_GRAPH_MAX_NODES = old


def test_path_step_cap():
    g = GfaGraph.parse(LINEAR_GFA)
    old = bwd.WEB_GRAPH_MAX_PATH_STEPS
    try:
        bwd.WEB_GRAPH_MAX_PATH_STEPS = 2
        out = bwd._extract_sample_graph(g, "HG00673", "1", "merged")
        assert len(out["path"]["steps"]) == 2
    finally:
        bwd.WEB_GRAPH_MAX_PATH_STEPS = old


def test_hap_label():
    assert bwd._hap_label("0") == "reference"
    assert bwd._hap_label("1") == "paternal"
    assert bwd._hap_label("2") == "maternal"
    assert bwd._hap_label("3") == "3"


def test_manifest_samples_merge():
    samples_by_graph = {
        "baseline": {"GRCh38": {"0"}, "HG00673": {"1"}},
        "merged": {"HG00673": {"1", "2"}},
    }
    samples = bwd._manifest_samples(samples_by_graph)
    by_name = {s["sample"]: s for s in samples}
    assert by_name["GRCh38"]["haplotypes"] == ["0"]
    assert sorted(by_name["HG00673"]["haplotypes"]) == ["1", "2"]
