"""Tests for pipeline/benchmark/graph_stats.py."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.benchmark.graph_stats import (  # noqa: E402
    _DSU, _n50, compute_stats, compare)
from pipeline.merge.gfa import GfaGraph  # noqa: E402


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------
LINEAR_GFA = """H\tVN:Z:1.0
S\ts1\tACGTACGTAC
S\ts2\tGGGGTTTTAA
S\ts3\tCCCCAAAATT
L\ts1\t+\ts2\t+\t0M
L\ts2\t+\ts3\t+\t0M
P\tGRCh38#0#chr21\ts1+,s2+,s3+\t*
"""

# One bubble: s2a / s2b are alternate alleles between s1 and s3.
BUBBLE_GFA = """H\tVN:Z:1.0
S\ts1\tACGTACGTAC
S\ts2a\tGGGGTTTTAA
S\ts2b\tGG
S\ts3\tCCCCAAAATT
L\ts1\t+\ts2a\t+\t0M
L\ts1\t+\ts2b\t+\t0M
L\ts2a\t+\ts3\t+\t0M
L\ts2b\t+\ts3\t+\t0M
P\tGRCh38#0#chr21\ts1+,s2a+,s3+\t*
P\tHG00673#1#chr21\ts1+,s2b+,s3+\t*
"""

# Two disconnected linear graphs concatenated — the disjoint-union signature.
DISJOINT_GFA = """H\tVN:Z:1.0
S\ta1\tACGTACGTAC
S\ta2\tGGGGTTTTAA
S\tb1\tCCCCAAAATT
S\tb2\tTTTTGGGGCC
L\ta1\t+\ta2\t+\t0M
L\tb1\t+\tb2\t+\t0M
P\tchunk1#GRCh38#0#chr21\ta1+,a2+\t*
P\tchunk2#GRCh38#0#chr21\tb1+,b2+\t*
"""


def _graph(text):
    return GfaGraph.parse(text)


# --------------------------------------------------------------------------
# Union-Find
# --------------------------------------------------------------------------
def test_dsu_single_component():
    d = _DSU(["a", "b", "c"])
    d.union("a", "b")
    d.union("b", "c")
    assert d.components() == [3]


def test_dsu_two_components():
    d = _DSU(["a", "b", "c", "d"])
    d.union("a", "b")
    d.union("c", "d")
    assert d.components() == [2, 2]


def test_dsu_singletons():
    d = _DSU(["a", "b", "c"])
    assert d.components() == [1, 1, 1]


def test_dsu_idempotent_union():
    d = _DSU(["a", "b"])
    d.union("a", "b")
    d.union("a", "b")
    assert d.components() == [2]


# --------------------------------------------------------------------------
# N50
# --------------------------------------------------------------------------
def test_n50_empty():
    assert _n50([]) == 0


def test_n50_uniform():
    assert _n50([10, 10, 10, 10]) == 10


def test_n50_known_case():
    # total=100, half=50; sorted desc 50,30,20 -> cum 50 >= 50 at L=50
    assert _n50([50, 30, 20]) == 50


def test_n50_skewed():
    # total=110, half=55; 50 -> 50 (<55), then 30 -> 80 (>=55) => 30
    assert _n50([50, 30, 20, 10]) == 30


# --------------------------------------------------------------------------
# Topology metrics
# --------------------------------------------------------------------------
def test_linear_graph_single_component():
    s = compute_stats(_graph(LINEAR_GFA), "linear")
    assert s["nodes"] == 3
    assert s["edges"] == 2
    assert s["components"]["count"] == 1
    assert s["components"]["largest"] == 3
    assert s["components"]["largest_fraction"] == 1.0


def test_linear_graph_has_no_branching():
    s = compute_stats(_graph(LINEAR_GFA), "linear")
    assert s["degrees"]["branching_nodes"] == 0


def test_bubble_graph_detects_branching():
    s = compute_stats(_graph(BUBBLE_GFA), "bubble")
    # s1 has out-degree 2, s3 has in-degree 2 -> both branching
    assert s["degrees"]["branching_nodes"] == 2
    assert s["degrees"]["max_out_degree"] == 2
    assert s["degrees"]["max_in_degree"] == 2


def test_bubble_graph_still_one_component():
    s = compute_stats(_graph(BUBBLE_GFA), "bubble")
    assert s["components"]["count"] == 1


def test_disjoint_union_inflates_component_count():
    """The core diagnostic: a disjoint union is detectable as N components."""
    s = compute_stats(_graph(DISJOINT_GFA), "disjoint")
    assert s["components"]["count"] == 2
    assert s["components"]["largest_fraction"] == 0.5


def test_total_bp_matches_sequence_length():
    s = compute_stats(_graph(LINEAR_GFA), "linear")
    assert s["total_bp"] == 30
    assert s["node_lengths"]["total_bp"] == 30


def test_tip_nodes_counted():
    s = compute_stats(_graph(LINEAR_GFA), "linear")
    # s1 has in-degree 0, s3 has out-degree 0
    assert s["degrees"]["tip_nodes"] == 2


def test_duplicate_links_not_double_counted():
    dup = LINEAR_GFA + "L\ts1\t+\ts2\t+\t0M\n"
    s = compute_stats(_graph(dup), "dup")
    assert s["degrees"]["unique_edges"] == 2


def test_path_steps_recorded():
    s = compute_stats(_graph(BUBBLE_GFA), "bubble")
    assert s["path_steps"]["count"] == 2
    assert s["path_steps"]["total_steps"] == 6


def test_complexity_fields_present():
    s = compute_stats(_graph(BUBBLE_GFA), "bubble")
    cx = s["complexity"]
    assert cx["edge_node_ratio"] == 1.0  # 4 edges / 4 nodes
    assert 0.0 <= cx["branching_fraction"] <= 1.0


# --------------------------------------------------------------------------
# Comparison verdicts
# --------------------------------------------------------------------------
def test_identical_graphs_are_equivalent():
    a = compute_stats(_graph(BUBBLE_GFA), "baseline")
    b = compute_stats(_graph(BUBBLE_GFA), "merged")
    result = compare(a, b)
    assert result["verdict"] == "EQUIVALENT"
    assert result["n_fail"] == 0


def test_disjoint_merge_flagged_divergent():
    baseline = compute_stats(_graph(LINEAR_GFA), "baseline")
    merged = compute_stats(_graph(DISJOINT_GFA), "merged")
    result = compare(baseline, merged)
    assert result["verdict"] == "DIVERGENT"
    comp = [c for c in result["checks"] if c["metric"] == "component_count"][0]
    assert comp["status"] == "FAIL"


def test_path_count_check_is_strict():
    baseline = compute_stats(_graph(BUBBLE_GFA), "baseline")
    merged = compute_stats(_graph(LINEAR_GFA), "merged")  # 1 path vs 2
    result = compare(baseline, merged)
    paths = [c for c in result["checks"]
             if c["metric"] == "paths_plus_walks"][0]
    assert paths["strict"] is True
    assert paths["status"] == "FAIL"


def test_tolerance_band_allows_small_node_drift():
    """A 2% node-count difference should pass at the default 5% tolerance."""
    baseline = compute_stats(_graph(BUBBLE_GFA), "baseline")
    merged = compute_stats(_graph(BUBBLE_GFA), "merged")
    merged["nodes"] = int(baseline["nodes"] * 1.02) or baseline["nodes"]
    result = compare(baseline, merged, tolerance_pct=5.0)
    nodes = [c for c in result["checks"] if c["metric"] == "nodes"][0]
    assert nodes["status"] == "PASS"


def test_every_check_has_required_fields():
    a = compute_stats(_graph(BUBBLE_GFA), "baseline")
    b = compute_stats(_graph(BUBBLE_GFA), "merged")
    for c in compare(a, b)["checks"]:
        for field in ("metric", "baseline", "merged", "delta_pct",
                      "strict", "status"):
            assert field in c


# --------------------------------------------------------------------------
# Script availability
# --------------------------------------------------------------------------
@pytest.mark.parametrize("script", [
    "pipeline/linear/run_dipcall.sh",
    "pipeline/linear/run_svim_asm.sh",
    "pipeline/benchmark/benchmark_variants.sh",
])
def test_shell_scripts_exist_and_are_executable(script):
    root = os.path.join(os.path.dirname(__file__), "..")
    path = os.path.join(root, script)
    assert os.path.exists(path), f"{script} missing"
    assert os.access(path, os.X_OK), f"{script} not executable"


# --------------------------------------------------------------------------
# Regressions found while testing against team branches
# --------------------------------------------------------------------------
W_LINE_GFA = """H\tVN:Z:1.1
S\ts1\tACGTACGTAC
S\ts2a\tGGGGTTTTAA
S\ts2b\tGG
S\ts3\tCCCCAAAATT
L\ts1\t+\ts2a\t+\t0M
L\ts1\t+\ts2b\t+\t0M
L\ts2a\t+\ts3\t+\t0M
L\ts2b\t+\ts3\t+\t0M
P\tGRCh38#0#chr21\ts1+,s2a+,s3+\t*
W\tHG00673\t1\tchr21\t0\t22\t>s1>s2b>s3
W\tHG00673\t2\tchr21\t0\t30\t>s1>s2a>s3
"""

DANGLING_GFA = "H\tVN:Z:1.0\nS\ta\tACGT\nL\ta\t+\tGHOST\t+\t0M\n"


def test_walks_counted_alongside_paths():
    """Real PGGB/Minigraph-Cactus output uses W-lines, not just P-lines."""
    s = compute_stats(_graph(W_LINE_GFA), "w")
    assert s["paths"] == 1
    assert s["walks"] == 2
    assert s["path_steps"]["count"] == 3


def test_walk_names_use_pansn_form():
    s = compute_stats(_graph(W_LINE_GFA), "w")
    assert "HG00673#1#chr21" in s["path_steps"]["per_path"]


def test_empty_graph_does_not_divide_by_zero():
    s = compute_stats(_graph("H\tVN:Z:1.0\n"), "empty")
    assert s["nodes"] == 0
    assert s["complexity"]["edge_node_ratio"] == 0.0
    assert s["node_lengths"]["n50"] == 0


def test_dangling_link_is_reported_not_crashed():
    """A truncated chunk GFA yields edges pointing at absent segments."""
    s = compute_stats(_graph(DANGLING_GFA), "dangle")
    assert s["degrees"]["dangling_links"] == 1
    assert s["components"]["count"] == 1


def test_zero_baseline_delta_is_none_not_infinity():
    """REGRESSION: float('inf') serializes as bare `Infinity`, which is not
    valid strict JSON and makes JavaScript's JSON.parse() throw."""
    import json
    a = compute_stats(_graph("H\tVN:Z:1.0\nS\tn1\tACGT\n"), "baseline")
    b = compute_stats(_graph(LINEAR_GFA), "merged")
    result = compare(a, b)
    for c in result["checks"]:
        assert c["delta_pct"] != float("inf")
    json.dumps(result, allow_nan=False)  # raises if Infinity/NaN present


def test_undefined_ratio_counts_as_failure():
    """A zero baseline with a nonzero merged value is undefined, not a pass."""
    a = compute_stats(_graph("H\tVN:Z:1.0\nS\tn1\tACGT\n"), "baseline")
    b = compute_stats(_graph(LINEAR_GFA), "merged")
    for c in compare(a, b)["checks"]:
        if c["delta_pct"] is None:
            assert c["status"] == "FAIL"


def test_lex_tsv_row_schema():
    """The TSV row must match lex_testing/metrics/runs.tsv column order."""
    import tempfile
    from pipeline.benchmark.graph_stats import (
        append_lex_tsv_row, LEX_TSV_COLUMNS)
    d = tempfile.mkdtemp()
    tsv = os.path.join(d, "runs.tsv")
    s = compute_stats(_graph(W_LINE_GFA), "demo")
    append_lex_tsv_row(s, tsv, wall_s=12, peak_mib=100, threads=8)
    lines = open(tsv).read().strip().split("\n")
    assert lines[0].split("\t") == LEX_TSV_COLUMNS
    row = dict(zip(LEX_TSV_COLUMNS, lines[1].split("\t")))
    assert row["label"] == "demo"
    assert row["sequences"] == "3"      # 1 path + 2 walks
    assert row["status"] == "OK"


def test_lex_tsv_appends_without_duplicating_header():
    import tempfile
    from pipeline.benchmark.graph_stats import append_lex_tsv_row
    d = tempfile.mkdtemp()
    tsv = os.path.join(d, "runs.tsv")
    s = compute_stats(_graph(LINEAR_GFA), "a")
    append_lex_tsv_row(s, tsv)
    append_lex_tsv_row(s, tsv)
    lines = open(tsv).read().strip().split("\n")
    assert len(lines) == 3  # header + 2 rows


# --------------------------------------------------------------------------
# W-line spec conformance (GFA 1.1)
#
# Four of five team branches shipped a Walk parser with a phantom
# `step_count` column and a comma-separated walk field. Real PGGB and vg
# output uses `W Sample Hap Seq Start End >n1<n2` with NO step_count.
# These tests pin the correct behaviour so the bug cannot come back
# through a branch merge.
# --------------------------------------------------------------------------
W_ORIENT_GFA = """H\tVN:Z:1.1
S\ts1\tACGT
S\ts2\tGGGG
S\ts3\tTTTT
L\ts1\t+\ts2\t+\t0M
L\ts1\t+\ts3\t+\t0M
P\tGRCh38#0#chr21\ts1+,s2+\t*
W\tHG00673\t1\tchr21\t0\t8\t>s1<s2
W\tHG00673\t2\tchr21\t0\t8\t>s1>s3
"""


def test_wline_reverse_orientation_parsed():
    """`<` means the walk traverses that segment in reverse."""
    g = _graph(W_ORIENT_GFA)
    assert g.walks[0].path == ["s1+", "s2-"]


def test_wline_has_seven_columns_not_eight():
    """A step_count column would make it eight and break every parser."""
    g = _graph(W_ORIENT_GFA)
    assert len(g.walks[0].to_gfa().split("\t")) == 7


def test_wline_roundtrips_orientation():
    g = _graph(W_ORIENT_GFA)
    assert g.walks[0].to_gfa().split("\t")[6] == ">s1<s2"


def test_walk_counted_in_stats_with_orientation():
    s = compute_stats(_graph(W_ORIENT_GFA), "w")
    assert s["walks"] == 2
    assert s["paths"] == 1
    assert s["path_steps"]["count"] == 3


# --------------------------------------------------------------------------
# Naming / config alignment with the rest of the pipeline
# --------------------------------------------------------------------------
def test_web_metrics_block_matches_dashboard_schema():
    """web/app/page.tsx reads data.metrics.<label>.{nodes,edges,paths}."""
    from pipeline.benchmark.graph_stats import web_metrics_block
    stats = {
        "baseline": compute_stats(_graph(LINEAR_GFA), "baseline"),
        "merged": compute_stats(_graph(BUBBLE_GFA), "merged"),
    }
    block = web_metrics_block(stats)
    for label in ("baseline", "merged"):
        for key in ("nodes", "edges", "paths"):
            assert key in block[label], "%s.%s missing" % (label, key)


def test_web_metrics_block_skips_comparison_key():
    """`comparison` is not a graph and must not become a dashboard card."""
    from pipeline.benchmark.graph_stats import web_metrics_block
    a = compute_stats(_graph(LINEAR_GFA), "baseline")
    stats = {"baseline": a, "comparison": compare(a, a)}
    assert "comparison" not in web_metrics_block(stats)


def test_web_metrics_paths_counts_walks_too():
    """The dashboard's Paths card should include W-line haplotypes."""
    from pipeline.benchmark.graph_stats import web_metrics_block
    stats = {"g": compute_stats(_graph(W_LINE_GFA), "g")}
    assert web_metrics_block(stats)["g"]["paths"] == 3   # 1 P + 2 W


def test_results_dir_read_from_config():
    """merge_graphs.py honours output.results_dir; this must not diverge."""
    import tempfile
    from pipeline.benchmark.graph_stats import _results_dir_from_config
    d = tempfile.mkdtemp()
    cfg = os.path.join(d, "pipeline.yaml")
    with open(cfg, "w") as f:
        f.write("output:\n  results_dir: custom_results\n")
    assert _results_dir_from_config(cfg) == "custom_results"


def test_results_dir_falls_back_when_config_absent():
    from pipeline.benchmark.graph_stats import _results_dir_from_config
    assert _results_dir_from_config("/nonexistent/pipeline.yaml") == "results"
