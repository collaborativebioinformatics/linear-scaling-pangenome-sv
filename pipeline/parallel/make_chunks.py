"""Create overlapping chunk manifest. Pairwise overlap = overlap_bp."""
import csv, os, yaml


def create_chunks(chrom, ref_start, ref_end, chunk_size=5000000, overlap=100000):
    """Create chunks with pairwise overlap = exactly overlap_bp.

    Each chunk extends overlap/2 on each side (except at edges).
    Adjacent chunks share exactly overlap_bp of overlapping sequence.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if overlap >= chunk_size:
        raise ValueError(
            f"overlap {overlap} too large for chunk_size {chunk_size}: "
            f"needs overlap < chunk_size")
    chunks, cid, pos = [], 1, ref_start
    while pos < ref_end:
        ce = min(pos + chunk_size, ref_end)
        ext_left = overlap // 2 if pos > ref_start else 0
        ext_right = overlap // 2 if ce < ref_end else 0
        chunks.append({
            "chunk_id": f"chunk_{cid:04d}", "chrom": chrom,
            "reference_start": max(ref_start, pos - ext_left),
            "reference_end": min(ref_end, ce + ext_right),
            "core_start": pos, "core_end": ce,
            "overlap_left": ext_left, "overlap_right": ext_right,
            "pairwise_overlap_bp": overlap,
            "status": "pending",
        })
        pos = ce
        cid += 1
    return chunks


def write_manifest(chunks, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter=chr(9), fieldnames=[
            "chunk_id", "chrom", "reference_start", "reference_end",
            "core_start", "core_end", "overlap_left", "overlap_right",
            "pairwise_overlap_bp", "status"])
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
