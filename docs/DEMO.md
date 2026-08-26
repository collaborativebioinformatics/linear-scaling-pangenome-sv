# Demo Walkthrough

## 30-Second Demo Script

1. Here are the same five genomic paths.
2. The left graph was built conventionally in one PGGB run.
3. The right graph was divided into overlapping regions and constructed in parallel.
4. Our algorithm then reassembled those regional graphs.
5. We measure whether that approach is faster.
6. More importantly, we test whether it preserves the biology.
7. Click a haplotype to follow the same sample through both graphs.
8. Green sites agree. Orange sites identify boundary/representation differences.

## Running Locally

```bash
make demo    # generates synthetic data, builds graphs, merges, exports JSON
make test    # runs 17+ tests
make web     # starts Next.js dev server at http://localhost:3000
```

## Running on DNAnexus

See [DNANEXUS.md](DNANEXUS.md) for cloud setup and execution.

## Expected Output

After `make demo`:

```
work/demo/
├── baseline.gfa              # Monolithic synthetic graph (502 nodes)
├── merged.gfa                # Merged chunk graph (599 nodes)
├── chunk_manifest.tsv        # 3 overlapping chunks
└── chunks/
    ├── chunk_0001.gfa
    ├── chunk_0002.gfa
    └── chunk_0003.gfa

web/public/data/latest.json   # Web visualizer data
```