#!/usr/bin/env python3
"""Generate PGGB configuration JSON from config/pipeline.yaml.

Output is passed as pggb_config_json string input to ALL applets.
No YAML dependency needed inside applet workers - uses stdlib json.
"""
import hashlib, json, os, sys, yaml


def load_config():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    pggb = cfg.get("pggb", {})
    params = pggb.get("params", {})
    # NOTE: PGGB v0.6.0 does NOT have a -w option.
    # -K = mash k-mer (mash_kmer), -k = seqwish minimum match length (match_length)
    return {
        "image": pggb.get("image", "ghcr.io/pangenome/pggb:latest"),
        "image_tag": pggb.get("image_tag", "latest"),
        "threads": pggb.get("threads", 8),
        "minimum_identity": params.get("minimum_identity", 90),
        "segment_length": params.get("segment_length", 5000),
        "mash_kmer": params.get("mash_kmer", 31),
        "match_length": params.get("match_length", 29),
        "path_jump_max": params.get("path_jump_max", 0),
        "edge_jump_max": params.get("edge_jump_max", 0),
    }


def main():
    config = load_config()
    raw = json.dumps(config, sort_keys=True)
    config["config_sha256"] = hashlib.sha256(raw.encode()).hexdigest()[:16]
    output_path = sys.argv[1] if len(sys.argv) > 1 else "work/pggb_config.json"
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    json.dump(config, open(output_path, "w"), indent=2)
    # Print to stdout for shell callers
    print(json.dumps(config))


if __name__ == "__main__":
    main()