"""Generate chunk_mapping.tsv from the real PAF data.
Maps each haplotype to its HPRC source contig coordinates.
"""
import csv
import os

def main():
    paf_path = "work/chunks/chr21_per_target.paf"
    out_path = "results/preparation/chunk_mapping.tsv"

    seen = set()
    rows = []

    # GRCh38 reference (only once)
    ref = {
        "sample": "GRCh38", "haplotype": "0", "contig": "chr21",
        "source_start": "20000000", "source_end": "21000000", "strand": "+",
    }
    rows.append(ref)
    seen.add("GRCh38#0")

    with open(paf_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            tname = cols[5]  # HG00673#1#JAHBBZ020000061.1
            ts = cols[7]
            te = cols[8]
            strand = cols[4]
            parts = tname.split("#")
            if len(parts) < 3:
                continue
            key = f"{parts[0]}#{parts[1]}"
            if key in seen or parts[0] == "GRCh38":
                continue
            seen.add(key)
            rows.append({
                "sample": parts[0],
                "haplotype": parts[1],
                "contig": "#".join(parts[2:]),
                "source_start": ts,
                "source_end": te,
                "strand": strand,
            })

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "sample", "haplotype", "contig",
            "source_start", "source_end", "strand"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")
    for r in rows:
        print(f"  {r['sample']}#{r['haplotype']}  "
              f"{r['contig']}:{r['source_start']}-{r['source_end']}")


if __name__ == "__main__":
    main()