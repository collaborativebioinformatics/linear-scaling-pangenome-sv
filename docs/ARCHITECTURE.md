# Architecture

## Overview

```
HPRC Release 2 FASTAs
         │
         ▼
   PREPARE chr21
         │
    ┌────┴────┐
    │         │
    ▼         ▼
 MONOLITHIC  PARALLEL CHUNKS
    │         │
    ▼     ┌───┼───┐
  PGGB   PGGB PGGB PGGB
    │     │   │   │
    ▼     └───┼───┘
 baseline     │
    │         ▼
    │       MERGE
    │         │
    └────┬────┘
         ▼
     VALIDATE
         │
    ┌────┼────┐
    │    │    │
   paths VCF topology
    │    │    │
    └────┼────┘
         ▼
     BENCHMARK
         │
         ▼
       JSON
         │
         ▼
  WEB VISUALIZER
```

## Key Modules

- `pipeline/merge/gfa.py` — GFA data model (no external deps)
- `pipeline/merge/merge_graphs.py` — Graph merging strategies
- `pipeline/parallel/make_chunks.py` — Chunk creation
- `pipeline/benchmark/` — Validation & metrics
- `pipeline/export/` — GFA to JSON for web
- `web/` — Next.js visualization

## Data Flow

GitHub → DNAnexus Cloud Workstation → DNAnexus Project Storage → JSON → Web