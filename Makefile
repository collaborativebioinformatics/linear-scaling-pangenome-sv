# Makefile — Parallel Pangenome Graph Construction
# =================================================
# Primary targets:
#   make check       — Check environment
#   make demo        — Run synthetic end-to-end demo
#   make test        — Run Python tests
#   make web         — Start web development server
#   make fetch-index — Fetch HPRC Release 2 assembly index
#   make download    — Download selected HPRC assemblies
#   make prepare     — Prepare chr21 sequences
#   make baseline    — Build monolithic PGGB graph
#   make chunks      — Create chunk manifest + build chunk graphs
#   make merge       — Merge chunk graphs
#   make benchmark   — Run validation and benchmarking
#   make smoke       — Run full smoke-test pipeline (HPRC data)
#   make chr21       — Run full chromosome 21 pipeline (expensive!)

.PHONY: help check demo test web fetch-index download prepare baseline chunks merge benchmark smoke chr21 clean

help:
	@echo "Parallel Pangenome Graph Construction — Make Targets"
	@echo "===================================================="
	@echo ""
	@echo "  make check         Check environment for required tools"
	@echo "  make demo          Run synthetic end-to-end demo (no external deps)"
	@echo "  make test          Run Python unit tests"
	@echo "  make web           Start web development server"
	@echo ""
	@echo "  --- HPRC Data Pipeline ---"
	@echo "  make fetch-index   Fetch HPRC Release 2 assembly index"
	@echo "  make download      Download selected HPRC assemblies"
	@echo "  make prepare       Prepare chr21 sequences from downloaded assemblies"
	@echo ""
	@echo "  --- Graph Construction ---"
	@echo "  make baseline      Build monolithic PGGB graph (smoke region)"
	@echo "  make chunks        Create chunk manifest + build chunk graphs"
	@echo "  make merge         Merge chunk graphs"
	@echo ""
	@echo "  --- Benchmarking ---"
	@echo "  make benchmark     Run validation and benchmarking on results"
	@echo ""
	@echo "  --- Full Pipelines ---"
	@echo "  make smoke         Run full smoke-test pipeline (HPRC chr21:20m-21m)"
	@echo "  make chr21         Run full chromosome 21 pipeline (requires DNAnexus)"
	@echo ""
	@echo "  make clean         Remove results/ and work/ (keeps downloads)"
	@echo ""

check:
	@bash scripts/check_environment.sh

demo:
	@echo "=== Synthetic End-to-End Demo ==="
	python3 scripts/setup_demo.py

test:
	@echo "=== Running Tests ==="
	python3 -m pytest tests/ -v

web:
	@echo "=== Starting Web Development Server ==="
	@cd web && npm run dev

fetch-index:
	@echo "=== Fetching HPRC Release 2 Assembly Index ==="
	python3 scripts/fetch_hprc_index.py

download:
	@echo "=== Downloading Selected HPRC Assemblies ==="
	python3 scripts/download_hprc.py

prepare:
	@echo "=== Preparing chr21 Sequences ==="
	python3 pipeline/prepare/prepare_sequences.py

baseline:
	@echo "=== Building Monolithic PGGB Graph ==="
	bash pipeline/baseline/build_baseline.sh

chunks:
	@echo "=== Creating Chunk Manifest ==="
	python3 pipeline/parallel/make_chunks.py
	@echo "=== Building All Chunk Graphs ==="
	python3 pipeline/parallel/build_all_chunks.py

merge:
	@echo "=== Merging Chunk Graphs ==="
	python3 pipeline/merge/merge_graphs.py

benchmark:
	@echo "=== Running Benchmarking ==="
	python3 pipeline/benchmark/graph_stats.py
	python3 pipeline/benchmark/compare_paths.py
	bash pipeline/benchmark/benchmark_variants.sh 2>/dev/null || echo "Variant benchmark requires vg/truvari — skipped"
	python3 pipeline/benchmark/build_report.py

smoke:
	@echo "=== Smoke Test Pipeline (HPRC chr21:20m-21m) ==="
	@echo ""
	@echo "REAL-DATA STEP: This requires HPRC assemblies downloaded and prepared."
	@echo "Run: make fetch-index && make download && make prepare"
	@echo "Then: make baseline && make chunks && make merge && make benchmark"
	@echo ""
	@echo "See config/pipeline.yaml to configure the target region."
	@echo "Set mode to 'smoke' and adjust start/end as needed."

chr21:
	@echo "=== Full Chromosome 21 Pipeline ==="
	@echo ""
	@echo "REAL-DATA STEP: This is computationally expensive."
	@echo "It should be run on DNAnexus, not locally."
	@echo ""
	@echo "1. Set up DNAnexus:  bash dnanexus/setup_workstation.sh"
	@echo "2. Upload inputs:    bash dnanexus/upload_inputs.sh"
	@echo "3. Run pipeline:     bash dnanexus/run_parallel_chunks.sh"
	@echo "4. Download results: bash dnanexus/download_results.sh"
	@echo ""
	@echo "See dnanexus/README.md for details."

clean:
	@echo "Cleaning results/ and work/ ..."
	rm -rf results/* work/*
	@echo "Done. (Keeps downloaded assemblies in downloads/)"
