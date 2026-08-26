"""
build_all_chunks.py — Reference-anchored chunk FASTA preparation.

Only the reference is sliced by coordinate. Haplotypes are included fully.
PGGB handles alignment within each chunk.
"""
import argparse
import csv
import os
import subprocess


def load_manifest(path="work/chunks/chunk_manifest.tsv"):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true",
                        help="Run PGGB on each chunk")
    parser.add_argument("--input",
                        default="results/preparation/chr21_multi.fa")
    args = parser.parse_args()

    chunks = load_manifest()
    if not chunks:
        print("No chunk manifest. Run make_chunks.py first.")
        return

    multi = args.input
    if not os.path.exists(multi):
        print(f"Not found: {multi}")
        return

    os.makedirs("work/chunks", exist_ok=True)

    # Load sequences
    seqs = {}
    cur_h, cur_s = None, []
    with open(multi) as f:
        for line in f:
            if line.startswith(">"):
                if cur_h:
                    seqs[cur_h] = "".join(cur_s)
                cur_h = line[1:].strip()
                cur_s = []
            else:
                cur_s.append(line.strip())
        if cur_h:
            seqs[cur_h] = "".join(cur_s)

    print(f"Loaded {len(seqs)} sequences, {len(chunks)} chunks")

    for c in chunks:
        cid = c["chunk_id"]
        chrom = c["chrom"]
        st = int(c["reference_start"])
        en = int(c["reference_end"])
        op = f"work/chunks/{cid}.fa"

        if os.path.exists(op) and os.path.getsize(op) > 0:
            print(f"  EXISTS {cid}")
            continue

        with open(op, "w") as fout:
            for hdr, seq in seqs.items():
                if hdr.startswith("GRCh38"):
                    # Reference: slice by known coordinate
                    cs = seq[st:en]
                    fout.write(f">{hdr}\n")
                    for i in range(0, len(cs), 80):
                        fout.write(cs[i:i + 80] + "\n")
                else:
                    # Haplotypes: write FULL sequence
                    # PGGB aligns them within the chunk context
                    fout.write(f">{hdr}\n")
                    for i in range(0, len(seq), 80):
                        fout.write(seq[i:i + 80] + "\n")

        print(f"  {cid}: ref sliced, haps full -> {op}")

    if args.execute:
        for c in chunks:
            cid = c["chunk_id"]
            gfa = f"work/chunks/{cid}.gfa"
            if os.path.exists(gfa):
                print(f"  SKIP {cid}: exists")
                continue
            subprocess.run(
                ["bash", "pipeline/parallel/build_chunk.sh", cid],
                check=False)

    print(f"Done. {len(chunks)} chunks.")


if __name__ == "__main__":
    main()