"""Create overlapping chunk manifest. Pairwise overlap = overlap_bp.

Two strategies (config: parallel.strategy):
  fixed       - fixed stride of chunk_size_bp (original behaviour)
  haploblocks - core boundaries snapped to haploblock edges (BED), merged
                to [min_chunk_bp, max_chunk_bp]; overlap logic unchanged.
"""
import csv, os, math, yaml


def _chunks_from_cores(chrom, ref_start, ref_end, cores, overlap):
    """Turn core intervals (contiguous, tiling [ref_start, ref_end)) into
    manifest rows. Each chunk extends overlap/2 on each side (except at
    edges), so adjacent chunks share exactly overlap_bp."""
    chunks = []
    for cid, (pos, ce) in enumerate(cores, start=1):
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
    return chunks


def create_chunks(chrom, ref_start, ref_end, chunk_size=5000000, overlap=100000):
    """Fixed-stride chunks with pairwise overlap = exactly overlap_bp."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if overlap >= chunk_size:
        raise ValueError(
            f"overlap {overlap} too large for chunk_size {chunk_size}: "
            f"needs overlap < chunk_size")
    cores, pos = [], ref_start
    while pos < ref_end:
        ce = min(pos + chunk_size, ref_end)
        cores.append((pos, ce))
        pos = ce
    return _chunks_from_cores(chrom, ref_start, ref_end, cores, overlap)


def read_blocks_bed(path, chrom):
    """Read haploblock intervals (chrom<TAB>start<TAB>end, 0-based half-open)."""
    blocks = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split()[:3]
            if c == chrom:
                blocks.append((int(s), int(e)))
    blocks.sort()
    for (s0, e0), (s1, e1) in zip(blocks, blocks[1:]):
        if s1 < e0:
            raise ValueError(f"overlapping haploblocks: {s0}-{e0} / {s1}-{e1}")
    return blocks


def merge_blocks_to_cores(blocks, ref_start, ref_end, min_size, max_size):
    """Clip blocks to the target, fill gaps, merge small neighbours to reach
    min_size without exceeding max_size, split any single block > max_size
    into equal parts. Every core edge is a haploblock edge, except (a) the
    target's own ends and (b) cuts inside an oversized block."""
    # 1. clip + fill gaps so cores tile [ref_start, ref_end) exactly
    tiles, pos = [], ref_start
    for s, e in blocks:
        s, e = max(s, ref_start), min(e, ref_end)
        if e <= pos:
            continue
        if s > pos:                      # gap between blocks -> keep as its own tile
            tiles.append((pos, s))
        tiles.append((max(s, pos), e))
        pos = e
    if pos < ref_end:
        tiles.append((pos, ref_end))
    # 2. split oversized tiles
    split = []
    for s, e in tiles:
        n = math.ceil((e - s) / max_size)
        step = math.ceil((e - s) / n)
        split += [(s + i * step, min(e, s + (i + 1) * step)) for i in range(n)]
    # 3. greedy merge up to >= min_size, never exceeding max_size
    cores, cs, ce = [], split[0][0], split[0][1]
    for s, e in split[1:]:
        if ce - cs < min_size and e - cs <= max_size:
            ce = e
        else:
            cores.append((cs, ce)); cs, ce = s, e
    cores.append((cs, ce))
    # 4. fold an undersized tail into its neighbour (allowed to exceed max by 25%)
    if len(cores) > 1 and cores[-1][1] - cores[-1][0] < min_size \
            and cores[-1][1] - cores[-2][0] <= max_size * 1.25:
        cores = cores[:-2] + [(cores[-2][0], cores[-1][1])]
    return cores


def create_chunks_from_blocks(chrom, ref_start, ref_end, blocks, min_size,
                              max_size, overlap):
    if overlap >= min_size:
        raise ValueError(f"overlap {overlap} must be < min_chunk_bp {min_size}")
    if not (0 < min_size <= max_size):
        raise ValueError("need 0 < min_chunk_bp <= max_chunk_bp")
    cores = merge_blocks_to_cores(blocks, ref_start, ref_end, min_size, max_size)
    return _chunks_from_cores(chrom, ref_start, ref_end, cores, overlap)


def write_manifest(chunks, path):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter=chr(9), fieldnames=[
            "chunk_id", "chrom", "reference_start", "reference_end",
            "core_start", "core_end", "overlap_left", "overlap_right",
            "pairwise_overlap_bp", "status"])
        w.writeheader()
        w.writerows(chunks)


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default="config/pipeline.yaml")
    ap.add_argument("--out", default="work/chunks/chunk_manifest.tsv")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    t = cfg["target"]
    p = cfg["parallel"]
    if p.get("strategy", "fixed") == "haploblocks":
        blocks = read_blocks_bed(p["blocks_file"], t["chromosome"])
        chunks = create_chunks_from_blocks(
            t["chromosome"], t["start"], t["end"], blocks,
            p["min_chunk_bp"], p["max_chunk_bp"], p["overlap_bp"])
    else:
        chunks = create_chunks(t["chromosome"], t["start"], t["end"],
                               p["chunk_size_bp"], p["overlap_bp"])
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    write_manifest(chunks, a.out)
    print(f"{len(chunks)} chunks ({p.get('strategy', 'fixed')}) -> {a.out}")


if __name__ == "__main__":
    main()
