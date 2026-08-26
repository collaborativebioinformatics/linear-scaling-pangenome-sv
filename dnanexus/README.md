# DNAnexus Deployment

## Architecture

```
GitHub (code)
   │
   ▼ git clone
DNAnexus Cloud Workstation (Ubuntu 24.04 compute)
   │
   ├── python3 scripts/fetch_hprc_index.py  ──►  official HPRC index
   ├── python3 scripts/download_hprc.py --execute  ──►  S3 → local
   ├── dx upload work/downloads/*.fa --destination /data/hprc/  ──►  persistent
   │
   ├── python3 pipeline/parallel/make_chunks.py  ──►  chunk manifest
   ├── bash dnanexus/run_parallel_chunks.sh  ──►  scatter-gather dx jobs
   │
   ├── python3 pipeline/merge/merge_graphs.py  ──►  merge step
   ├── dx upload results/merge/merged.gfa --destination /graphs/merged/
   │
   ├── python3 pipeline/benchmark/build_report.py  ──►  metrics JSON
   └── dx upload results/benchmark/report.json --destination /benchmark/

DNAnexus Project Storage (persistent)
   /data/hprc/          HPRC assemblies
   /data/reference/     GRCh38
   /graphs/baseline/    monolithic PGGB
   /graphs/chunks/      parallel PGGB per chunk
   /graphs/merged/      stitched graph
   /benchmark/          metrics reports
   /web/                web JSON → web/public/data/latest.json
```

## Quick Start (inside DNAnexus Cloud Workstation)

```bash
# 1. Launch workstation
dx run --instance-type mem3_ssd1_v2_x16 --ssh app-cloud_workstation

# 2. Clone and setup
git clone <repo-url>
cd pangenome-parallel
bash dnanexus/setup_workstation.sh

# 3. Fetch HPRC index & download assemblies
python3 scripts/fetch_hprc_index.py
python3 scripts/download_hprc.py --execute

# 4. Upload to persistent storage
bash dnanexus/create_project_dirs.sh
bash dnanexus/upload_inputs.sh

# 5. Run pipeline
python3 pipeline/parallel/make_chunks.py
bash dnanexus/run_parallel_chunks.sh

# 6. Merge & benchmark
python3 pipeline/merge/merge_graphs.py
python3 pipeline/benchmark/graph_stats.py
python3 pipeline/benchmark/build_report.py

# 7. Upload results
dx upload results/merge/merged.gfa --destination /graphs/merged/
dx upload results/benchmark/report.json --destination /benchmark/
```

## Environment

`DX_PROJECT_CONTEXT_ID` is automatically set inside Cloud Workstations.
All `dnanexus/*.sh` scripts detect it automatically.

## Instance Types

Configure via `DX_INSTANCE_TYPE` environment variable:

```bash
export DX_INSTANCE_TYPE=mem3_ssd1_v2_x32  # for larger assemblies
```

## Data Flow

| Step | Source | Destination | Method |
|------|--------|-------------|--------|
| Fetch index | GitHub (HPRC official) | `work/manifests/` | `python3 scripts/fetch_hprc_index.py` |
| Download | HPRC public S3 | `work/downloads/` | `aws s3 cp --no-sign-request` |
| Upload | `work/downloads/` | `/data/hprc/` | `dx upload` |
| Prepare | `/data/hprc/` | `/data/prepared/` | `pipeline/prepare/` |
| Chunk PGGB | `/data/prepared/` | `/graphs/chunks/` | `dx run pggb` |
| Merge | `/graphs/chunks/` | `/graphs/merged/` | `pipeline/merge/` |
| Results | `/graphs/merged/` | local → web | `dx download` → `scripts/sync_web_results.py` |
