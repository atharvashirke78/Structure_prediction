# Supplementary Materials — Antibody Structure Prediction Benchmark

This folder contains the per-target and per-model datasets underlying the tables and figures in the manuscript.

- **S1_interface_fingerprint_per_target.csv** — Per-target, per-model heavy/light interface-contact fingerprint results (Section 3.6, Table 3a). Columns: PDB_ID, Model, ref_n_pairs (number of reference interface-contact pairs), pred_n_pairs, n_shared, jaccard, fraction_ref_recovered, side_chains_reconstructed (True for IgFold only; see Section 2.4).

- **S2_per_target_rmsd_all_metrics.csv** — Per-target, per-model Fv RMSD (conventional/masked and unfiltered/zero), CDR-H3 RMSD, and CDR-L3 RMSD, for all 10 models across the 58-structure reference set. Underlies Table 1, Table 1a, and Section 3.1-3.4.

- **S3_table1_summary.csv** — Per-model summary statistics (median conventional/unfiltered Fv RMSD, delta median/max, unfiltered RMSD max, CDR-H3 RMSD median) as reported in Table 1.

- **S4_wilcoxon_conventional.csv / S5_wilcoxon_unfiltered.csv** — Full pairwise Wilcoxon signed-rank test results (all 45 model pairs), including the matched-pairs rank-biserial effect size (`rank_biserial_r`) for every pair, unadjusted p-values, and Holm-Bonferroni-corrected p-values, under each protocol (Section 3.3).

- **S6_bootstrap_ci_fv.csv / S7_bootstrap_ci_cdr_h3.csv** — Bootstrap 95% confidence intervals (10,000 resamples) for each model's median Fv RMSD (both protocols) and median CDR-H3 RMSD (Table 1a).

- **S8_rejected_fraction_per_target.csv** — Per-target, per-model rejected-Cα-pair fraction under the conventional protocol's iterative outlier rejection, and the corresponding conventional-unfiltered RMSD delta (Section 4.1, Table 2a).

- **S9_stereochemistry_raw.csv** — Per-target, per-model Ramachandran-outlier percentage and MolProbity-style clashscore (Section 3.5, Table 1b, Table 2).


- **S10_figure5_confidence_vs_rmsd_data.csv** — Per-target confidence score (0-100 scale) and conventional Fv RMSD for each of the eight models with a per-target confidence score. This is the single authoritative dataset underlying both Table 7 (confidence-RMSD Pearson correlations) and Figure 5; the two are computed from and report identically against this file.

- **S11_per_model_metric_coverage.csv** — Per-model count of valid (non-missing) observations for each metric reported anywhere in the manuscript (predictions, both Fv RMSD protocols, CDR-H3 RMSD, CDR-L3 RMSD, TM-score, interface fingerprint, stereochemistry), out of the full 58-target reference set. `N_tm_score_independently_reverifiable` flags the three models (Chai-1, ESMFold2, OpenFold3) whose full 58-target TM-score raw output is still independently retrievable from storage; for the other seven models, only a later targeted correction for specific outlier targets was preserved (see Section 2.4 in the manuscript for the corresponding caveat).

- **S12_computational_provenance.csv** — Exact software versions and access dates/endpoints for every tool and third-party service used in this benchmark: PyMOL (with the precise `cmd.super()` call for each of the two superposition protocols), ANARCI, Arpeggio, PDBFixer, TM-align (tmtools), the CCTBX/MolProbity-equivalent stereochemistry tools, the ColabFold/MMseqs2 MSA server, the OpenFold3 inference API, and the statistical software and random seed used for the Wilcoxon tests and bootstrap confidence intervals. Model checkpoint versions and exact release dates for the other nine structure-prediction models are not uniformly recoverable across this project's multi-month timeline and are not reported individually; each model was run once via its respective inference pathway (local installation, hosted API, or third-party server), as described in Section 2.2 of the manuscript.

- **S13_heavy_chain_only_sensitivity_analysis.csv** — Sensitivity analysis for the four model-target cells in which a prediction is missing its light chain entirely (IntelliFold and OpenFold3 on 1MCO, OpenFold3 on 5DK3, IgFold on 6MTS). Compares the primary manuscript statistics (these four cells treated as missing) against the originally computed heavy-chain-only fallback statistics (these four cells scored as heavy-chain-only predictions against the full heavy+light reference) for each of the three affected models. The largest difference is OpenFold3's unfiltered median (0.032 Å); all other differences are at or below 0.011 Å. None of these differences changes any statistically significant pairwise comparison reported in Section 3.3 (see also Section 2.4 and Section 4.6 of the manuscript).

All files use PDB identifiers matching the SAbDab reference structures described in Section 2.1, and model-name codes matching those used in the corresponding pipeline scripts (e.g., AB_Fold_3 = AbFold, Alpha_Fold_2 = AlphaFold 2, Intellifold = IntelliFold, Chai1 = Chai-1).


## Missing values: what they mean and why

Several columns in these files contain missing (NaN/empty) values. Each case has a specific, documented cause — none reflect a computation failure on a target the model otherwise handled correctly.

**S1_interface_fingerprint_per_target.csv — 4 of 580 rows (0.7%) are NaN** (`pred_n_pairs`, `n_shared`, `jaccard`, `fraction_ref_recovered`): Intellifold and OpenFold3 on 1MCO, OpenFold3 on 5DK3, and IgFold on 6MTS. In every one of these four cases, that specific model's own predicted structure for that target contains only a heavy (H) chain, with no light (L) chain present at all — confirmed by direct inspection of the corrected PDB files. Since the interface-fingerprint metric measures heavy/light contacts, there is no light chain to select and no interface to score; Arpeggio's `/L//` selection genuinely fails with "entity not found" because that entity does not exist in the file, not because of a bug in the scoring pipeline. Both 1MCO and 5DK3 have reference structures with both chains present (verified by direct inspection), and in both cases most models correctly predicted both chains — 8 of 9 other models did so for 1MCO, and 9 of 9 other models did so for 5DK3. These are therefore genuine, model-specific prediction failures (Intellifold and OpenFold3 for 1MCO; OpenFold3 for 5DK3; IgFold for 6MTS), not properties of the reference structures. This corrects an earlier assumption in this project's working notes that 1MCO/5DK3 were reference-level heavy-chain-only edge cases; that assumption is false for all four of these model-target pairs (see also the manuscript's note in Section 4.6).

**S2_per_target_rmsd_all_metrics.csv — `CDR_L3_RMSD` is sparse by design, not by error.** AB_Fold_3, ABodyBuilder3, and Ibex return a usable CDR-L3 value for only 1 of 58 targets each; this is disclosed explicitly in the manuscript (Section 2.4) as a genuine model/pipeline limitation — these three models' light-chain CDR-L3 loop region could not be reliably isolated across the reference set with the residue-window approach used here, so CDR-L3 RMSD is reported for these three models only where usable data exists, and the manuscript's main text and Table 1 do not include a CDR-L3 comparison for this reason. IntelliFold and OpenFold3's extra CDR-L3 NaN values correspond to the same heavy-chain-only prediction gaps described above for S1 (a model can't have a CDR-L3 loop RMSD if it has no light chain in that particular prediction). IgFold has two CDR-L3 NaN values in S2 (6MTS, 4FQC); only 6MTS matches an S1 missing-light-chain case — 4FQC's prediction has both chains present, so that NaN reflects a separate case where the CDR-L3 residue window failed to isolate a usable loop region in that specific prediction, not a missing chain. `CDR_H3_RMSD` has 1-2 NaN values per model across the panel, corresponding to targets where the CDR-H3 loop was unresolved in either the reference or the prediction (for example, a genuinely disordered/missing loop in the crystal structure) rather than a scoring failure. `Masked_Fv_RMSD` and `Zero_Fv_RMSD` have zero missing values for any model across all 58 targets — every model produced a scoreable Fv (heavy-chain, or heavy+light where both are present) for every target.

## Note: Constant-domain coverage on the two full-length reference structures

Two of the 58 reference structures in this panel — 1MCO and 5DK3 — are full-length antibodies with resolved constant domains (heavy chain ~428-441 residues, light chain ~216-218 residues in the reference), rather than isolated Fv fragments. This note documents each model's predicted chain length and composition for these two targets, as a record of why a hinge-region (VH-CH1 / VL-CL orientation) comparison was not attempted in the main manuscript, rather than as a finding about hinge-region accuracy itself — no rotational angle or hinge-related error metric was computed for any model on any target in this project.

Chain length (Cα count) and chain identity for each model's prediction, both targets:

| Model | 1MCO: H / L (total) | 5DK3: H / L (total) | Full-length H+L on both targets? |
|---|---|---|---|
| AlphaFold 2 | 428 / 216 (644) | 444 / 218 (662) | Yes |
| Boltz-2 | 428 / 216 (644) | 444 / 218 (662) | Yes |
| Chai-1 | 428 / 216 (644) | 444 / 218 (662) | Yes |
| ESMFold2 | 428 / 216 (644) | 444 / 218 (662) | Yes |
| AbFold | 130 / 130 (260) | 130 / 130 (260) | No (Fv-length both targets) |
| ABodyBuilder3 | 117 / 110 (227) | 120 / 111 (231) | No (Fv-length both targets) |
| IgFold | 117 / 110 (227) | 120 / 111 (231) | No (Fv-length both targets) |
| Ibex | 117 / 216 (333) | 120 / 111 (231) | No (1MCO: Fv-length heavy chain, full-length light chain; 5DK3: Fv-length both) |
| IntelliFold | 216 / — (216, H only) | 444 / 218 (662) | No (H-only on 1MCO, full-length on 5DK3 — inconsistent between targets) |
| OpenFold3 | 216 / — (216, H only) | 218 / — (218, H only) | No (H-only on both targets) |

Only 4 of the 10 models (AlphaFold 2, Boltz-2, Chai-1, ESMFold2) produced a near-full-length prediction with both heavy and light chains present for both of these two targets. This observation is specific to these two targets and these ten model deployments; it is not a general characterization of which models are architecturally capable of full-length prediction (several of the models classified here as "Fv-length" — AbFold, ABodyBuilder3, IgFold, Ibex — are Fv-only architectures by design, as noted in Section 2.1 of the manuscript, and were never expected to produce a full-length structure regardless of input). IntelliFold's and OpenFold3's within-model inconsistency between the two targets (full-length vs. heavy-chain-only, or heavy-chain-only on both) is more notable, since neither is an Fv-only architecture by design; the cause was not investigated further here.

Given that only 4 of 10 models had a heavy+light full-length prediction on both of only 2 available targets, a hinge-region comparison across the panel would have covered at most a 4-model, n=2-target subset — far too small to support a general claim and inconsistent with the statistical treatment used for every other metric in this benchmark (Section 2.4, Section 3.3). This is why the main manuscript's scope is restricted to the Fv region (Section 2.1) and why no hinge-angle or constant-domain orientation metric appears anywhere in this project's results.


## Figures

The `figures/` subfolder contains the final, manuscript-numbered image files for all 9 figures used in the paper, matching the numbering in the current draft:

- **Figure1_workflow_diagram.png** — benchmarking workflow (Section 2)
- **Figure2_ranking_reversal.png** — model ranking reversal between protocols (Section 3.1)
- **Figure3_max_unfiltered_rmsd.png** — maximum unfiltered Fv RMSD by model (Section 3.2)
- **Figure4_bootstrap_ci.png** — bootstrap confidence intervals for Fv RMSD and CDR-H3 RMSD (Section 3.3)
- **Figure5_stereochemical_validity_paired.png** — Ramachandran-outlier rate and clashscore, median vs. worst-case (Section 3.5)
- **Figure6_abfold_clash.png** — AbFold's prediction for target 1R24, steric-clash hotspot highlighted (Section 3.5)
- **Figure7_alphafold2_clash.png** — AlphaFold 2's prediction for the same target 1R24, for direct comparison with Figure 6 (Section 3.5)
- **Figure8_abodybuilder3_clash.png** — ABodyBuilder3's prediction for the same target 1R24, for direct comparison with Figure 6 (Section 3.5)
- **Figure9_confidence_vs_rmsd.png** — confidence score vs. Fv RMSD, one panel per model (Section 4.3)

Figures 6-8 were originally a single three-panel composite image; they were split into three separate figures, one per model, at the user's request.

## Data revision note (2026-08-23)

Table 1, Table 1a (now the S3/S6/S7 datasets below), and all downstream statistics were recomputed treating the four heavy-chain-only cells listed under S13 as missing data rather than as heavy-chain-only fallback comparisons; this is a small correction (at most a 0.03 Å shift in any displayed median (OpenFold3's unfiltered median; all others ≤0.011 Å)) documented in the manuscript's Section 2.4, Section 4.6, and S13. The two standalone `wilcoxon_effect_sizes_conventional.csv` / `wilcoxon_effect_sizes_unfiltered.csv` files present in an earlier version of this folder have been removed as redundant — their `rank_biserial_r` column is already included in S4/S5 below.

- **S15_all_model_structures.tar.gz** — The corrected, chain-relabeled prediction structures (PDB format, H/L chain labeling, 1..N per-chain residue numbering) for all ten benchmarked models plus the reference structure, for all 58 targets (11 folders × 58 files = 638 structure files, ~37 MB compressed): `AB_Fold_3/` (AbFold), `ABodyBuilder3/`, `Alpha_Fold_2/` (AlphaFold 2), `Boltz_2/`, `Chai1/`, `ESMFold2/`, `Ibex/`, `IgFold/`, `Intellifold/` (IntelliFold), `OpenFold3/`, and `Reference/`. These are the exact structure files scored to produce every RMSD, CDR-H3/L3, TM-score, and stereochemistry result reported in the manuscript and in S1-S9 above. An earlier, now-superseded `ESM_2/` model directory (a prior single-sequence ESM-2-based predictor dropped from the current ten-model panel in favor of ESMFold2) is not included here.
