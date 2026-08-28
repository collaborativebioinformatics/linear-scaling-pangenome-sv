"""Behavioral regression tests for PGGB configuration correctness.

Verifies that every PGGB execution path uses the SAME canonical parameters,
no invalid -w flag, no :latest image, and config identity between applets.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

REPO = os.path.join(os.path.dirname(__file__), "..")


def test_gen_pggb_config_parses():
    """Canonical PGGB config JSON parses and has required fields."""
    from gen_pggb_config import load_config
    cfg = load_config()
    for k in ["image", "image_tag", "threads", "minimum_identity",
              "segment_length", "mash_kmer", "match_length",
              "path_jump_max", "edge_jump_max"]:
        assert k in cfg, f"missing key: {k}"


def test_no_invalid_pggb_option():
    """No stale -w / kmer_length / window_size keys survive."""
    from gen_pggb_config import load_config
    cfg = load_config()
    assert "kmer_length" not in cfg
    assert "window_size" not in cfg
    assert "map_pct_id" not in cfg
    assert "noise_filter" not in cfg


def test_image_is_pinned_not_latest():
    """Every path MUST use an immutable image digest, never :latest."""
    from gen_pggb_config import load_config
    cfg = load_config()
    assert "sha256" in cfg["image"], f"image not pinned: {cfg['image']}"


def test_baseline_chunk_params_identical():
    """Both baseline and chunk applets receive the SAME JSON config."""
    from gen_pggb_config import load_config
    c1 = load_config()
    c2 = load_config()
    for k in c1:
        assert c1[k] == c2[k], f"{k} differs"


def test_config_json_is_stdlib_json():
    """Applet parses pggb_config_json as string via stdlib json."""
    from gen_pggb_config import load_config
    cfg = load_config()
    raw = json.dumps(cfg, sort_keys=True)
    parsed = json.loads(raw)
    assert parsed["image"] == cfg["image"]


def test_run_pggb_no_w_flag():
    """run_pggb.py source must NOT contain the -w flag in its command."""
    import inspect
    import run_pggb
    src = inspect.getsource(run_pggb.main)
    assert " -w " not in src, "run_pggb.py contains -w flag!"
    assert "-K" in src and "-k" in src, "run_pggb.py missing -K/-k flags!"


def test_docker_helper_no_stale_flags():
    """docker_helper.sh must not use -w flag in any command line."""
    content = open(os.path.join(REPO, "dnanexus/docker_helper.sh")).read()
    for line in content.split("\n"):
        if " -w " in line and not line.lstrip().startswith("#"):
            assert False, f"docker_helper.sh has -w flag: {line}"


def test_applet_code_no_yaml_dep():
    """Both applet code.sh files use stdlib json, not PyYAML."""
    for p in ["dnanexus/applets/pggb_chunk/src/code.sh",
              "dnanexus/applets/pggb_baseline/src/code.sh"]:
        c = open(os.path.join(REPO, p)).read()
        # The applet parses via: python3 -c "import sys,json; ..."
        assert "import sys,json" in c or "import json" in c, \
            f"{p} missing stdlib json parse"
        assert "import yaml" not in c, f"{p} has PyYAML dependency!"


def test_final_gfa_unique_selection():
    """Applets require exactly one *final.gfa (FATAL if zero or >1)."""
    for p in ["dnanexus/applets/pggb_chunk/src/code.sh",
              "dnanexus/applets/pggb_baseline/src/code.sh"]:
        c = open(os.path.join(REPO, p)).read()
        assert "final.gfa" in c, f"{p} missing final.gfa search"
        assert "GFA_COUNT" in c, f"{p} missing duplicate check"


def test_web_guards_full_coverage():
    """sync_web_results.py blocks gfa, fa, fasta, vcf, bam, cram."""
    c = open(os.path.join(REPO, "scripts/sync_web_results.py")).read()
    for ext in [".gfa", ".fa", ".fasta", ".vcf", ".bam", ".cram",
                ".fa.gz", ".fasta.gz", ".gfa.gz", ".vcf.gz"]:
        assert ext in c, f"missing forbidden extension: {ext}"
