# Analysis pipeline for EFSA TXG-MAP deliverable D3.1.
# Every stage reads and writes inside this repository; see README.md.

M1_PY = cd analysis/method1 && python3
M3_PY = cd analysis/method3 && python3

.PHONY: help method1 method1-fetch method1-analyse method3-snapshot method3-coverage \
        method3-genesets builder-genesets comparison figures-clean

help:
	@echo "make method1           - Method 1: fetch from the pinned AOP-Wiki RDF, then analyse"
	@echo "make method1-analyse   - Method 1 from the cached raw tables, without the endpoint"
	@echo "make method3-snapshot  - Method 3: snapshot the molAOP Builder corpus and KE frame"
	@echo "make method3-coverage  - Method 3: coverage of the four networks and FIG10"
	@echo "make method3-genesets  - Method 3: KE gene sets from the approved pathway mappings"
	@echo "make builder-genesets  - Methods 2 and 3 gene sets from the builder's GMT exports"
	@echo "make comparison        - the cross-method comparison and FIG11"
	@echo "make figures-clean     - the same figures without the provenance strip"
	@echo "Set OFFLINE=1 to skip live BridgeDb calls where a stage supports it."

method1-fetch:
	$(M1_PY) m1_01_fetch.py $(if $(ALLOW_DRIFT),--allow-drift,)

method1-analyse:
	$(M1_PY) m1_02_frames.py
	$(M1_PY) m1_03_harmonise.py $(if $(OFFLINE),--skip-bridgedb,)
	$(M1_PY) m1_04_stats.py
	$(M1_PY) m1_05_figures.py

method1: method1-fetch method1-analyse

method3-snapshot:
	$(M3_PY) m3_01_snapshot.py
	$(M3_PY) m3_02_frame.py

method3-coverage:
	$(M3_PY) m3_12_coverage.py
	$(M3_PY) m3_15_figures.py

method3-genesets:
	$(M3_PY) m3_13_genesets.py

builder-genesets:
	cd analysis && python3 builder_genesets.py $(if $(OFFLINE),--skip-bridgedb,)

comparison:
	cd analysis/comparison && python3 cmp_01_compare.py
	cd analysis/comparison && python3 cmp_02_source_drift.py
	cd analysis/comparison && python3 cmp_03_cascade.py
	cd analysis/comparison && python3 cmp_05_hub_kes.py
	cd analysis/comparison && python3 cmp_04_figures.py

figures-clean: export D31_FIG_CLEAN=1
figures-clean:
	$(M1_PY) m1_05_figures.py
	$(M3_PY) m3_15_figures.py
	cd analysis/comparison && python3 cmp_04_figures.py
	python3 analysis/figures/dia_01_schematics.py
