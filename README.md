# Structure_prediction

Code, structures, and results supporting the manuscript *A Controlled Benchmark of Ten Antibody Structure Predictors: Accuracy, Stereochemistry and Superposition Sensitivity*.

## Repository layout

- **`scripts/`** — the benchmarking pipeline: chain identification (`chain_detection.py`), the main scoring pipeline for both superposition protocols and all downstream metrics (`benchmark_master_pipeline.py`, `rescore_full.py`), CDR-H3/CDR-L3 IMGT-based scoring (`score_cdrh3_imgt.py`, `score_cdrl3_imgt.py`), interface-contact recovery (`interface_recovery_worker.py`), stereochemical validation (`rama_clash_worker.py`), FASTA/input preparation (`input_fasta.py`), and dataset validation (`Data validate.py`). See `scripts/README_metrics_setup.md` for the metrics-computation environment setup.

- **`structures/`** — the reference and predicted structures used for scoring: one subfolder per model (`ABodyBuilder3`, `AbFold`, `AlphaFold2`, `Boltz2`, `Chai1`, `ESMFold2`, `Ibex`, `IgFold`, `IntelliFold`, `OpenFold3`) plus `Reference`, each with one structure file per target (58 targets per model/reference, matched Fv-only input as described in Methods).

- **`data/`** — per-model per-target results tables, the full statistical outputs referenced in the manuscript (Wilcoxon signed-rank tests, bootstrap confidence intervals, stereochemistry, CDR-H3/L3 sensitivity analyses, the duplicate-target robustness check), the computational provenance record, and the manuscript's figures. See `data/README.md` for a full file-by-file description.

- **`requirements.txt`** — Python dependencies for the scripts in `scripts/`. PyMOL (Open-Source PyMOL 3.1.0) must be installed separately (see file for instructions).

## Citation

If you use this benchmark, please cite the manuscript (see the main text for full citation details once published).
