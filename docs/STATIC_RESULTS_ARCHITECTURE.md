# Static Results Architecture

## Data Flow

```
              EXPENSIVE, RUN ONCE
                   |
       full 1 Mb DNAnexus experiment
                   |
  +----------------+----------------+
  |                |                |
timing        equivalence       Truvari
  |                |                |
  +----------------+----------------+
                   |
            frozen final JSON
                   |
           visualization selector
                   |
         25-50 kb compact graph
                   |
              Next.js site
```

## FULL 1 MB (benchmark region)

The following metrics are computed once, on the **full 1 Mb** target
(`GRCh38 chr21:20,000,000-21,000,000`):

- **timing** — baseline wall time, parallel wall time, sum worker seconds,
  orchestration overhead, speedup
- **graph comparison** — baseline vs stitched nodes/edges/components
- **path equivalence** — exact path SHA256 matches per haplotype
- **boundary comparison** — stitch boundary PASS/WARN/FAIL counts
- **Truvari** — baseline-vs-stitched variant preservation

These are frozen in `results/final_run/*.json` and copied to
`web/public/data/final/*.json`.

## SELECTED SMALL REGION (visualization only)

The Cytoscape canvas renders a **fixed representative subregion**
(default 50 kb) selected once by `pipeline/export/select_visualization_region.py`.

This region is **for visualization only**. It does NOT determine timing,
equivalence, or Truvari metrics.

## Commands

| Command | Purpose | Recomputes? |
|---------|---------|-------------|
| `make freeze-web-results` | Export existing results to static JSON | No (export only) |
| `make run-real-benchmark` | Run PGGB + stitch on DNAnexus | Yes (explicit) |
| `make run-variant-comparison` | Run vg + Truvari | Yes (explicit) |
| `npm run build` / `npm run dev` | Build/serve the site | **Never** |

## Web Data

The website loads ONLY `web/public/data/final/*.json`:

```
manifest.json
timing.json
graph_comparison.json
path_comparison.json
boundary_comparison.json
variant_comparisons.json
visualization_region.json
```

No DNAnexus calls. No Python backend. No genomics computation.
No raw GFA/FASTA/VCF/BAM/CRAM in `web/public`.

## Sample Selection Behavior

Selecting a sample/haplotype changes **path highlighting only**.
The graph topology is loaded once. No data refetch. No recomputation.
