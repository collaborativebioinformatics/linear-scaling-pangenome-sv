"""Create overlapping chunk manifest for parallel graph building."""
import csv
import os
import yaml


def create_chunks(chrom, ref_start, ref_end, chunk_size=5000000, overlap=100000):
    chunks, cid, pos = [], 1, ref_start
    while pos < ref_end:
        ce = min(pos + chunk_size, ref_end)
        ol = overlap if pos > ref_start else 0
        or_ = overlap if ce < ref_end else 0
        chunks.append({
            "chunk_id": f"chunk_{cid:04d}", "chrom": chrom,
            "reference_start": max(ref_start, pos - ol),
            "reference_end": min(ref_end, ce + or_),
            "core_start": pos, "core_end": ce,
            "overlap_left": ol, "overlap_right": or_,
            "status": "pending",
        })
        pos = ce
        cid += 1
    return chunks


def write_manifest(chunks, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter=chr(9), fieldnames=[
            "chunk_id", "chrom", "reference_start", "reference_end",
            "core_start", "core_end", "overlap_left", "overlap_right", "status",
        ])
        w.writeheader()
        w.writerows(chunks)


def main():
    cfg = yaml.safe_load(open("config/pipeline.yaml"))
    t = cfg["target"]
    p = cfg["parallel"]
    chunks = create_chunks(t["chromosome"], t["start"], t["end"],
                           p["chunk_size_bp"], p["overlap_bp"])
    os.makedirs("work/chunks", exist_ok=True)
    write_manifest(chunks, "work/chunks/chunk_manifest.tsv")
    print(f"{len(chunks)} chunks -> work/chunks/chunk_manifest.tsv")


if __name__ == "__main__":
    main()
