"""run_pggb.py — Run PGGB locally with the canonical configuration.

Reads the SAME config/pipeline.yaml parameters as gen_pggb_config.py and the
DNAnexus applets. There must be exactly one PGGB parameter definition, so this
script mirrors gen_pggb_config.load_config() and never invents its own flags.

PGGB v0.6.0 flags (verified against pinned digest):
    -p = mapping identity (minimum_identity)
    -s = segment length (segment_length)
    -K = mash k-mer size (mash_kmer)
    -k = seqwish minimum match length (match_length)
    -j = smoothxg path-jump max (path_jump_max)
    -e = smoothxg edge-jump max (edge_jump_max)

NOTE: -w is NOT a valid PGGB v0.6.0 option.

Usage:
    python3 scripts/run_pggb.py input.fa output_dir [threads]
"""
import json, os, subprocess, sys, time, yaml
from pathlib import Path

# Import the canonical loader so there is a single source of truth.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_pggb_config import load_config


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
    if len(sys.argv) >= 4:
        cfg["threads"] = int(sys.argv[3])

    num_paths = sum(1 for line in open(input_fa) if line.startswith(">"))

    # Build canonical PGGB command — identical parameters to the applets.
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{os.path.dirname(input_fa)}:/data/input",
        "-v", f"{outdir}:/data/output",
        cfg["image"],
        "bash", "-lc",
        "samtools faidx /data/input/{input_file} && pggb "
        "-i /data/input/{input_file} "
        "-o /data/output "
        "-t {threads} "
        "-n {num_paths} "
        "-p {minimum_identity} "
        "-s {segment_length} "
        "-K {mash_kmer} "
        "-k {match_length} "
        "-j {path_jump_max} "
        "-e {edge_jump_max}".format(
            input_file=os.path.basename(input_fa),
            threads=cfg["threads"],
            num_paths=num_paths,
            minimum_identity=cfg["minimum_identity"],
            segment_length=cfg["segment_length"],
            mash_kmer=cfg["mash_kmer"],
            match_length=cfg["match_length"],
            path_jump_max=cfg["path_jump_max"],
            edge_jump_max=cfg["edge_jump_max"],
        ),
    ]

    instance_type = os.environ.get("DX_INSTANCE_TYPE", "local")
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

    log_path = os.path.join(outdir, "pggb_run.log")
    with open(log_path, "w") as f:
        f.write(result.stdout)
        f.write(result.stderr)

    if result.returncode != 0:
        print(f"FATAL: PGGB failed (exit {result.returncode})")
        print(result.stderr[-2000:])
        sys.exit(1)

    print(f"PGGB finished in {wall_seconds:.1f}s")

    final_gfa = find_final_gfa(outdir)
    canonical = os.path.join(outdir, "final.gfa")
    os.rename(final_gfa, canonical)
    gfa_size = os.path.getsize(canonical)
    print(f"  Final GFA: {canonical} ({gfa_size} bytes)")

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
        "config_sha256": cfg.get("config_sha256", ""),
        "pggb_params": {
            "minimum_identity": cfg["minimum_identity"],
            "segment_length": cfg["segment_length"],
            "mash_kmer": cfg["mash_kmer"],
            "match_length": cfg["match_length"],
            "path_jump_max": cfg["path_jump_max"],
            "edge_jump_max": cfg["edge_jump_max"],
        },
    }
    meta_path = os.path.join(outdir, "run_metadata.json")
    json.dump(metadata, open(meta_path, "w"), indent=2)
    print(f"  Metadata: {meta_path}")

    return canonical


if __name__ == "__main__":
    main_path = main()
    print(f"OK: {main_path}")