"""freeze_web_results.py — Freeze scientific results into static web JSON.

This is an EXPORT operation, not analysis. It NEVER calls dx, PGGB, vg, or
Truvari. It reads existing artifacts and writes compact deterministic JSON.

    results/baseline/baseline.gfa  ─┐
    results/merge/merged.gfa        ─┤
    results/validation/*.json       ─┤  read once
    results/benchmark/*.json        ─┘
                                     │
                              results/final_run/
                                     │
                              web/public/data/final/

Genomic files are NEVER copied to web/public.

Usage:
    python3 pipeline/export/freeze_web_results.py
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

FINAL_DIR = "results/final_run"
WEB_FINAL_DIR = "web/public/data/final"

FORBIDDEN_EXT = (".gfa", ".fa", ".fasta", ".fa.gz", ".fasta.gz",
                 ".vcf", ".vcf.gz", ".bam", ".cram", ".gfa.gz")
WEB_MAX_FILE_MB = 10


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def _read_json(path, default):
    if path and os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except json.JSONDecodeError:
            return default
    return default


def _file_sha256(path):
    if not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_config_target():
    try:
        import yaml
        with open("config/pipeline.yaml") as f:
            cfg = yaml.safe_load(f)
        t = cfg.get("target", {})
        return {
            "reference": t.get("reference", "GRCh38"),
            "chromosome": t.get("chromosome", "chr21"),
            "start": int(t.get("start", 0)),
            "end": int(t.get("end", 0)),
        }
    except Exception:
        return {"reference": "GRCh38", "chromosome": "chr21",
                "start": 20000000, "end": 21000000}


def _validate_and_build():
    """Build the final_run manifest + status, honestly reporting what exists."""
    target = _load_config_target()
    target["length_bp"] = target["end"] - target["start"]

    baseline_validation = _read_json(
        "results/validation/baseline_path_validation.json",
        {"status": "MISSING"},
    )
    graph_comparison = _read_json(
        "results/benchmark/graph_comparison.json", {"status": "NOT_RUN"})
    path_comparison = _read_json(
        "results/benchmark/path_comparison.json", {"status": "NOT_RUN"})
    boundary_comparison = _read_json(
        "results/benchmark/boundary_comparison.json", {"status": "NOT_RUN"})
    variant_comparisons = _read_json(
        "results/benchmark/variant_comparisons.json", {"status": "NOT_RUN"})
    timing = _read_json("results/benchmark/timing.json", {"status": "NOT_RUN"})
    vis_region = _read_json("results/web/visualization_region.json", None)

    baseline_status = baseline_validation.get("status", "NOT_RUN")
    if baseline_status == "INVALID_FOR_1MB_BENCHMARK":
        baseline_level = "INVALID_FOR_1MB_BENCHMARK"
    elif baseline_status == "VALID_1MB_BENCHMARK":
        baseline_level = "REAL_HPRC_VALIDATED"
    elif baseline_status == "MISSING":
        baseline_level = "NOT_RUN"
    else:
        baseline_level = "UNVERIFIED"

    stitched_status = "NOT_RUN"
    merged_gfa = "results/merge/merged.gfa"
    if os.path.exists(merged_gfa):
        from pipeline.merge.gfa import GfaGraph
        g = GfaGraph.parse_file(merged_gfa)
        stitched_status = ("PRESENT_BUT_NOT_VALIDATED"
                           if g.edge_count() > 0 and g.path_count() > 0
                           else "PLACEHOLDER_LINEAR_ONLY")

    validation = {
        "baseline": baseline_level,
        "parallel": "NOT_RUN",
        "stitch": stitched_status,
        "path_equivalence": "NOT_RUN",
        "variant_equivalence": "NOT_RUN",
    }

    manifest = {
        "schema_version": "final-1",
        "dataset": target,
        "samples": [
            {"sample": "GRCh38", "haplotype": "0"},
            {"sample": "HG00673", "haplotype": "1"},
            {"sample": "HG00673", "haplotype": "2"},
            {"sample": "HG00733", "haplotype": "1"},
            {"sample": "HG00733", "haplotype": "2"},
        ],
        "validation": validation,
        "files": {
            "benchmark": "benchmark.json",
            "timing": "timing.json",
            "graph_comparison": "graph_comparison.json",
            "path_comparison": "path_comparison.json",
            "boundary_comparison": "boundary_comparison.json",
            "variant_comparisons": "variant_comparisons.json",
            "visualization_region": "visualization_region.json",
        },
        "provenance": {
            "git_sha": _git_sha(),
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "baseline_gfa_sha256": _file_sha256("results/baseline/baseline.gfa"),
            "merged_gfa_sha256": _file_sha256("results/merge/merged.gfa"),
        },
    }

    return {
        "manifest.json": manifest,
        "timing.json": timing,
        "graph_comparison.json": graph_comparison,
        "path_comparison.json": path_comparison,
        "boundary_comparison.json": boundary_comparison,
        "variant_comparisons.json": variant_comparisons,
        "visualization_region.json": vis_region or {
            "status": "NOT_SELECTED",
            "selection_method": "highest branch/path-diversity score",
        },
    }


def _guard_no_genomic(path):
    ext = os.path.splitext(path)[1].lower()
    if ext in FORBIDDEN_EXT or path.endswith(FORBIDDEN_EXT):
        print(f"  BLOCKED: {path} — genomic files not allowed in web/public")
        return False
    return True


def _guard_file_size(path, max_mb=WEB_MAX_FILE_MB):
    if os.path.exists(path) and os.path.getsize(path) > max_mb * 1024 * 1024:
        print(f"  SKIP {path}: > {max_mb}MB")
        return False
    return True


def main():
    files = _validate_and_build()

    os.makedirs(FINAL_DIR, exist_ok=True)
    for name, data in files.items():
        with open(os.path.join(FINAL_DIR, name), "w") as f:
            json.dump(data, f, indent=2, sort_keys=True)

    os.makedirs(WEB_FINAL_DIR, exist_ok=True)
    copied = 0
    for name, data in files.items():
        src = os.path.join(FINAL_DIR, name)
        if not _guard_no_genomic(name) or not _guard_file_size(src):
            continue
        shutil.copy2(src, os.path.join(WEB_FINAL_DIR, name))
        copied += 1

    print(f"freeze_web_results: wrote {len(files)} files to {FINAL_DIR}, "
          f"copied {copied} to {WEB_FINAL_DIR}")
    print(f"  baseline validation: "
          f"{files['manifest.json']['validation']['baseline']}")
    print(f"  stitched status:     "
          f"{files['manifest.json']['validation']['stitch']}")
    print(f"  variant comparison:  "
          f"{files['variant_comparisons.json'].get('status', 'NOT_RUN')}")


if __name__ == "__main__":
    main()