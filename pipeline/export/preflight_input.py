"""preflight_input.py — Strict input gate for the final 1 Mb DNAnexus run.

Verifies the regional FASTA before ANY PGGB submission:

    GRCh38 must be exactly 1,000,000 bp
    each HPRC haplotype ~regional scale (500 kb .. 2 Mb)
    GRCh38 must not be all-N
    total <= 10 Mb

Writes results/validation/final_input_preflight.json + .tsv
Exits non-zero on any failure.
"""
import csv
import hashlib
import json
import os
import sys

GRCH38_EXACT = 1_000_000
MAX_HAP = 2_000_000
MIN_HAP = 500_000
MAX_TOTAL = 10_000_000


def analyze(fasta_path):
    seqs = {}
    cur = None
    order = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = ""
                order.append(cur)
            elif cur is not None:
                seqs[cur] += line

    rows = []
    failed = []
    total = 0
    for name in order:
        seq = seqs[name]
        n = len(seq)
        total += n
        upper = seq.upper()
        counts = {b: upper.count(b) for b in "ACGTN"}
        n_pct = counts["N"] / n * 100 if n else 0
        gc_pct = (counts["G"] + counts["C"]) / n * 100 if n else 0
        sha = hashlib.sha256(seq.encode()).hexdigest()

        status = "PASS"
        reasons = []
        if name.startswith("GRCh38"):
            if n != GRCH38_EXACT:
                status = "FAIL"
                reasons.append(f"GRCh38 != {GRCH38_EXACT} bp")
            if n_pct > 90:
                status = "FAIL"
                reasons.append(f"GRCh38 overwhelmingly N ({n_pct:.1f}%)")
        else:
            if n > MAX_HAP:
                status = "FAIL"
                reasons.append(f"> {MAX_HAP} bp (full-contig risk)")
            if n < MIN_HAP:
                status = "FAIL"
                reasons.append(f"< {MIN_HAP} bp (partial coverage)")
            if n_pct > 50:
                status = "FAIL"
                reasons.append(f"high N ({n_pct:.1f}%)")

        rows.append({
            "sample": name.split("#")[0],
            "haplotype": name.split("#")[1] if "#" in name else "0",
            "fasta_header": name,
            "sequence_bp": n,
            "n_percent": round(n_pct, 3),
            "gc_percent": round(gc_pct, 3),
            "sha256": sha,
            "status": status,
            "reasons": reasons,
        })
        if status == "FAIL":
            failed.append((name, "; ".join(reasons)))

    if total > MAX_TOTAL:
        failed.append(("TOTAL", f"total {total} > {MAX_TOTAL}"))

    overall = "PASS" if not failed else "FAIL"
    return {"status": overall, "total_bp": total, "sequences": rows,
            "failures": [{"name": n, "reason": r} for n, r in failed]}


def main():
    fasta = sys.argv[1] if len(sys.argv) > 1 else \
        "work/final_input/chr21_1mb_multi.fa"
    if not os.path.exists(fasta):
        print(f"FATAL: {fasta} not found", file=sys.stderr)
        sys.exit(2)

    result = analyze(fasta)
    os.makedirs("results/validation", exist_ok=True)
    with open("results/validation/final_input_preflight.json", "w") as f:
        json.dump(result, f, indent=2)

    with open("results/validation/final_input_preflight.tsv", "w",
              newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=[
            "sample", "haplotype", "fasta_header", "sequence_bp",
            "n_percent", "gc_percent", "sha256", "status"])
        w.writeheader()
        for r in result["sequences"]:
            w.writerow({k: r[k] for k in w.fieldnames})

    print(f"preflight status={result['status']} total={result['total_bp']:,} bp")
    for r in result["sequences"]:
        print(f"  {r['fasta_header']:35s} {r['sequence_bp']:>10d} bp  "
              f"N%={r['n_percent']:.2f}  {r['status']}")
    if result["status"] == "FAIL":
        for f in result["failures"]:
            print(f"  FAIL: {f['name']}: {f['reason']}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
