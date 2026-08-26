.PHONY: help check demo test web setup fetch-index download prepare baseline chunks merge benchmark smoke chr21 clean deploy
TAB := $(shell printf '\t')

help:
	@echo "Targets: check demo test web deploy"
	@echo "  Data:   fetch-index download prepare"
	@echo "  Graph:  baseline chunks merge benchmark smoke chr21"
	@echo "  Utils:  setup clean"

check:; @bash scripts/check_environment.sh
demo:; @echo "=== Demo ===" && python3 scripts/setup_demo.py
test:; @python3 -m pytest tests/ -v
web:; @cd web && npm run dev
setup:; @pip install pyyaml pytest && cd web && npm install
fetch-index:; @python3 scripts/fetch_hprc_index.py
download:; @python3 scripts/download_hprc.py
prepare:; @python3 pipeline/prepare/prepare_sequences.py
baseline:; @bash pipeline/baseline/build_baseline.sh
chunks:; @python3 pipeline/parallel/make_chunks.py && python3 pipeline/parallel/build_all_chunks.py
merge:; @python3 pipeline/merge/merge_graphs.py
benchmark:; @python3 pipeline/benchmark/graph_stats.py && python3 pipeline/benchmark/compare_paths.py && bash pipeline/benchmark/benchmark_variants.sh 2>/dev/null || true && python3 pipeline/benchmark/build_report.py
smoke:; @bash scripts/run_smoke_test.sh
chr21:; @echo "REAL-DATA STEP: run on DNAnexus. See dnanexus/README.md"
clean:; @rm -rf results/* work/* && echo "Cleaned"
deploy:; @cd web && npx vercel --prod
