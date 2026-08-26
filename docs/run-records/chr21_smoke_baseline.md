# Chromosome 21 smoke-test baseline

## Run

- Date: 2026-08-26
- Target: GRCh38 chr21:20,000,000-21,000,000
- Coordinate system: zero-based, half-open
- Method: monolithic PGGB
- PGGB version: 4225c6c
- Threads: 16
- Parameters: `-n 5 -p 90 -s 5000 -k 29 -j 0 -e 0`
- Runtime: 114 seconds

## Input paths

| Path | Length (bp) |
|---|---:|
| GRCh38#0#chr21 | 1,000,000 |
| HG00673#1#JAHBBZ020000061.1 | 1,001,738 |
| HG00673#2#JAHBBY020000010.1 | 1,001,667 |
| HG00733#1#JAHEPQ020000018.1 | 1,001,586 |
| HG00733#2#JAHEPP020000058.1 | 1,001,742 |

Each haplotype mapping covered the complete 1 Mb GRCh38 query with MAPQ 60 and forward orientation.

## Baseline graph

- Segments: 10,384
- Links: 14,178
- Paths: 5
- Walks: 0
- GFA size: approximately 1.5 MB
- SHA-256: `239bec903c275da4f56fd40a8410d5be3002bf1241bf56ef6170137a59f96561`
- DNAnexus location: `/graphs/baseline/chr21_20000000_21000000/`

The obsolete PGGB `-w 50000` option was removed because it is unsupported by PGGB version 4225c6c.
