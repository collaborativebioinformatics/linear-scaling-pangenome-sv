"""Save real DNAnexus timing for the final 1 Mb benchmark."""
import json
import subprocess
from pathlib import Path

JOBS = {
    "baseline": "job-JB8p1y00ZQvP86x4p5x7v8Fv",
    "c1": "job-JB8p2400ZQvPGXKyGJz78gVz",
    "c2": "job-JB8p2480ZQvG91xZY4bfFZPG",
    "c3": "job-JB8p24Q0ZQv1xgJYJpGkXF4b",
}

def describe(jid):
    raw = subprocess.check_output(["dx", "describe", jid, "--json"], text=True)
    d = json.loads(raw)
    return {
        "job_id": jid,
        "wall_seconds": (d["stoppedRunning"] - d["startedRunning"]) / 1000.0,
        "started_running_ms": int(d["startedRunning"]),
        "stopped_running_ms": int(d["stoppedRunning"]),
        "instance_type": "mem3_ssd1_v2_x16",
    }

bl = describe(JOBS["baseline"])
chunks = [describe(JOBS[k]) for k in ["c1", "c2", "c3"]]

par_wall = (max(x["stopped_running_ms"] for x in chunks) -
            min(x["started_running_ms"] for x in chunks)) / 1000.0
sum_w = sum(x["wall_seconds"] for x in chunks)
sp = round(bl["wall_seconds"] / par_wall, 2) if par_wall > 0 else None

out = {
    "benchmark_region": "GRCh38 chr21:20000000-21000000",
    "benchmark_length_bp": 1_000_000,
    "baseline": bl,
    "parallel": {
        "chunks": chunks,
        "graph_parallel_wall_seconds": round(par_wall, 1),
        "sum_worker_seconds": round(sum_w, 1),
    },
    "wall_clock_speedup": sp,
}

Path("results/final_run").mkdir(parents=True, exist_ok=True)
with open("results/final_run/timing.json", "w") as f:
    json.dump(out, f, indent=2)

print(f"Baseline: {bl['wall_seconds']:.0f}s")
print(f"Parallel wall: {par_wall:.0f}s")
print(f"Sum worker: {sum_w:.0f}s")
print(f"Speedup: {sp}x" if sp else "Speedup: N/A")
print("Saved: results/final_run/timing.json")