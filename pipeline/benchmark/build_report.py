"""Aggregate benchmark results into report JSON."""
import json
import os
import glob
from datetime import datetime


def main():
    rd = "results"
    dm = "synthetic"
    if os.path.exists("work/manifests/hprc_selected.csv"):
        with open("work/manifests/hprc_selected.csv") as f:
            if len(f.read()) > 50:
                dm = "real"
    report = dict(
        data_mode=dm, status="PARTIAL",
        merge_status="NOT_IMPLEMENTED",
        message="Overlap-aware stitching not implemented. Disjoint union only.",
        generated_at=datetime.now().isoformat(), components={},
    )
    for pattern, key in [
        (f"{rd}/benchmark/graph_metrics.json", "graph_metrics"),
        (f"{rd}/baseline/run_metadata.json", "run_metadata"),
    ]:
        if os.path.exists(pattern):
            report["components"][key] = json.load(open(pattern))
    for pattern, key in [
        (f"{rd}/baseline/*", "baseline_files"),
        (f"{rd}/merge/*", "merge_files"),
        (f"{rd}/benchmark/*", "benchmark_files"),
    ]:
        report["components"][key] = len(
            [f for f in glob.glob(pattern) if os.path.isfile(f)]
        )
    op = f"{rd}/benchmark/report.json"
    os.makedirs(os.path.dirname(op), exist_ok=True)
    json.dump(report, open(op, "w"), indent=2)
    print(f"Report: {op}")


if __name__ == "__main__":
    main()