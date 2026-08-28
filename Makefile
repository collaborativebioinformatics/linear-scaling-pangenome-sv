.PHONY: help check demo test web setup fetch-index verify-index download prepare-ref prepare prepare-seq baseline chunks merge benchmark linear dipcall svim-asm variants clean deploy dnanexus-setup dnanexus-pipeline freeze-web-results run-real-benchmark run-variant-comparison

help:
	@echo "Targets:"
	@echo "  DNAnexus:     make dnanexus-setup && make dnanexus-pipeline"
	@echo "  Local test:   make demo && make test && make web"
	@echo "  HPRC data:    make fetch-index && make download && make prepare-ref && make prepare"
	@echo "  Graph build:  make baseline && make chunks && make merge && make benchmark"
	@echo "  Validation:   make linear && make variants"

check:;	@bash scripts/check_environment.sh
demo:;	@echo "=== Demo ===" && python3 scripts/setup_demo.py && python3 pipeline/export/build_web_dataset.py && python3 scripts/sync_web_results.py
test:;	@python3 -m pytest tests/ -v
web:;	@cd web && npm run dev
setup:;	@pip install pyyaml pytest && cd web && npm install

dnanexus-setup:;	@bash dnanexus/setup_workstation.sh
dnanexus-pipeline:;	@bash dnanexus/run_pipeline.sh --upload

fetch-index:;	@python3 scripts/fetch_hprc_index.py
verify-index:;	@python3 scripts/fetch_hprc_index.py --verify-urls
download:;	@python3 scripts/download_hprc.py --execute
prepare-ref:;	@bash scripts/prepare_reference.sh
prepare: prepare-seq
prepare-seq:;	@python3 pipeline/prepare/prepare_sequences.py
baseline:;	@bash pipeline/baseline/build_baseline.sh
chunks:;	@python3 pipeline/parallel/make_chunks.py && python3 pipeline/parallel/build_all_chunks.py --execute
merge:;	@python3 pipeline/merge/merge_graphs.py
linear: dipcall
dipcall:;	@bash pipeline/linear/run_dipcall.sh
svim-asm:;	@bash pipeline/linear/run_svim_asm.sh
variants:;	@bash pipeline/benchmark/benchmark_variants.sh
benchmark:;	@python3 pipeline/benchmark/graph_stats.py && python3 pipeline/benchmark/compare_paths.py && python3 pipeline/benchmark/build_report.py

# ---- Static results architecture ----
freeze-web-results:;	@python3 pipeline/export/validate_baseline_paths.py && python3 pipeline/export/freeze_web_results.py
run-real-benchmark:;	@echo "PGGB baseline + parallel chunks + stitch (see dnanexus/run_pipeline.sh)" && bash dnanexus/run_pipeline.sh --upload
run-variant-comparison:;	@bash pipeline/benchmark/benchmark_variants.sh

clean:;	@rm -rf results/* work/* web/public/data/baseline.json web/public/data/merged.json && echo "Cleaned"
deploy:;	@cd web && npx vercel --prod