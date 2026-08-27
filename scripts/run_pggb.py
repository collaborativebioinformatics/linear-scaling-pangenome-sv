"""run_pggb.py — Run PGGB with canonical configuration from config/pipeline.yaml.

Reads pipeline.yaml for ALL PGGB parameters (single source of truth).
Records image, digest, full command, timestamps, instance_type, vcpus, and wall time in metadata.
Selects *final.gfa output exclusively (exactly one required).
Usage:
    python3 scripts/run_pggb.py input.fa output_dir [threads]
"""
import json, os, subprocess, sys, time, yaml
from pathlib import Path


def load_config():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    pggb_cfg = cfg.get("pggb", {})
    params = pggb_cfg.get("params", {})
    return {
        "image": pggb_cfg.get("image", "ghcr.io/pangenome/pggb:latest"),
        "image_tag": pggb_cfg.get("image_tag", "latest"),
        "threads": pggb_cfg.get("threads", 8) if len(sys.argv) < 4 else int(sys.argv[3]),
        "minimum_identity": params.get("minimum_identity", 90),
        "segment_length": params.get("segment_length", 5000),
        "kmer_length": params.get("kmer_length", 29),
        "window_size": params.get("window_size", 50000),
        "map_pct_id": params.get("map_pct_id", 0),
        "noise_filter": params.get("noise_filter", 0),
    }


def find_final_gfa(outdir):
    """Find exactly one *final.gfa in outdir. Zero or >1 = FATAL."""
    gfaf = list(Path(outdir).glob("*final.gfa"))
    if len(gfaf) == 0:
        # PGGB v0.6 produces subdir/ with {prefix}.final.gfa
        gfaf = list(Path(outdir).rglob("*final.gfa"))
    if len(gfaf) == 0:
        sys.exit(f"FATAL: No *final.gfa found in {outdir}")
    if len(gfaf) > 1:
        sys.exit(f"FATAL: Multiple *final.gfa found in {outdir}: {gfaf}")
    return str(gfaf[0])


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 scripts/run_pggb.py input.fa output_dir [threads]")
        sys.exit(1)

    input_fa = os.path.abspath(sys.argv[1])
    outdir = os.path.abspath(sys.argv[2])
    os.makedirs(outdir, exist_ok=True)

    cfg = load_config()
    num_paths = sum(1 for line in open(input_fa) if line.startswith(">"))

    # Build canonical PGGB command
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{os.path.dirname(input_fa)}:/data/input:ro",
        "-v", f"{outdir}:/data/output",
        cfg["image"],
        "pggb",
        "-i", f"/data/input/{os.path.basename(input_fa)}",
        "-o", "/data/output",
        "-t", str(cfg["threads"]),
        "-n", str(num_paths),
        "-p", str(cfg["minimum_identity"]),
        "-s", str(cfg["segment_length"]),
        "-k", str(cfg["kmer_length"]),
        "-w", str(cfg["window_size"]),
        "-j", str(cfg["map_pct_id"]),
        "-e", str(cfg["noise_filter"]),
    ]

    # Capture environment metadata
    instance_type = os.environ.get("DX_INSTANCE_TYPE", "unknown")
    vcpus_avail = os.cpu_count() or 0

    print(f"=== PGGB ({cfg['image_tag']}) ===")
    print(f"  Input: {input_fa} ({num_paths} paths)")
    print(f"  Output: {outdir}")
    print(f"  Threads: {cfg['threads']}")
    print(f"  Instance: {instance_type}")
    print(f"  VCPUs: {vcpus_avail}")
    print(f"  Command: {' '.join(cmd)}")

    start_ts = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_ts = time.time()
    wall_seconds = end_ts - start_ts

    # Write log
    log_path = os.path.join(outdir, "pggb_run.log")
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    if result.returncode != 0:
        print(f"FATAL: PGGB failed (exit {result.returncode})")
        print(result.stderr[-2000:])
        sys.exit(1)

    print(f"PGGB finished in {wall_seconds:.1f}s")

    # Locate exactly one *final.gfa
    final_gfa = find_final_gfa(outdir)
    canonical = os.path.join(outdir, "final.gfa")
    os.rename(final_gfa, canonical)
    gfa_size = os.path.getsize(canonical)
    print(f"  Final GFA: {canonical} ({gfa_size} bytes)")

    # Write metadata
    metadata = {
        "method": "pggb",
        "container": cfg["image"],
        "image_tag": cfg["image_tag"],
        "num_paths": num_paths,
        "threads": cfg["threads"],
        "instance_type": instance_type,
        "vcpus": vcpus_avail,
        "wall_seconds": round(wall_seconds, 1),
        "start_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(start_ts)),
        "stop_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(end_ts)),
        "status": "completed",
        "final_gfa": canonical,
        "gfa_size_bytes": gfa_size,
        "full_command": " ".join(cmd),
        "config_params": {k: v for k, v in cfg.items() if k != "image"},
    }
    meta_path = os.path.join(outdir, "run_metadata.json")
    json.dump(metadata, open(meta_path, "w"), indent=2)
    print(f"  Metadata: {meta_path}")

    return canonical


if __name__ == "__main__":
    main_path = main()
    print(f"OK: {main_path}")