# M33-Corbelli Paper 7 Complete Package

**Package name:** `m33_corbelli_paper7_complete_package.zip`

**Lane label:** `M33-Corbelli source-native thick-disk reconstruction`

**Included:** final board, locked Hermes result, MOND comparator, baryonic-convention sensitivity diagnostic, outer-disk disequilibrium audit, origin-trail addendum, single-system age assumption note, and archive-lock notes.

## Current locked status

- Final reconstructed Format B board: included.
- Hermes full-board result: pass/elevated under both bracketed Amber-age lanes.
- MOND fixed-simple comparator: clean on the same frozen board.
- M33 interpretation: boundary-case roadmap; Hermes main disk is clean/near-clean, outer disk is the localized stress zone.

## Guardrails

- No board edits, reruns, or model reinterpretations are performed in this package.
- No nested older ZIP packages are included. Current files are included directly.
- Files named `prelock` or `age_anchor_audit` are retained only as provenance for how the age bracket was derived. The controlling age decision is `m33_age_bracketed_amber_lock_addendum.md` and `m33_age_bracketed_amber_lock.csv`.
- This board is **source-native**, not SPARC-equivalent: SPS/pixel-SED stellar mass profile, radially varying M/L provenance, helium factor 1.33, H2 prescription, flaring stellar disk, gas half-thickness 0.5 kpc, and `Vbul=0`.
- The `f≈1.21` baryonic-convention branch is a SPARC-like unit-M/L diagnostic analogue, **not** a literal SPARC rotmod/unit-M/L board.

## Folder map

- `00_final_board/` — Final reconstructed board, table provenance, kernel, mass closure, overlay benchmark, and age-bracket lock.
- `01_hermes_result_locked/` — Locked Hermes-only execution package contents.
- `02_mond_comparator/` — MOND comparator execution on the exact same frozen board.
- `03_baryonic_convention_sensitivity/` — M/L and stellar-mass convention sensitivity diagnostic; includes the SPARC-like-not-SPARC-equivalent clarification.
- `04_outer_disk_disequilibrium_audit/` — Outer-disk disequilibrium and radial-mask audit.
- `05_interpretive_notes/` — Single-system age assumption note, origin-trail addendum, and accepted archive lock.
- `99_archive_locks/` — Canon/archive-lock memo for Hermes-only result preservation.
- `file_manifest.csv` — machine-readable manifest with path, byte size, SHA256, and description for every included file.

## File inventory

### `00_final_board/`
- `00_final_board/README_final_board_package.md` — Board-freeze package README from the bracketed-age final board package, renamed to avoid conflicting with the master README.
- `00_final_board/board_freeze_manifest.json` — Board-freeze package manifest from the bracketed-age final board package, renamed to board-level manifest.
- `00_final_board/corbelli_fig12_overlay_reconstruction.png` — Published-figure overlay benchmark image comparing reconstructed components against Corbelli et al. Figure 12 top-panel source-native curves.
- `00_final_board/m33_age_bracketed_amber_lock.csv` — Machine-readable Amber bracketed t50 age lock: primary 6.0±1.2 Gyr and conservative 7.1±1.3 Gyr.
- `00_final_board/m33_age_bracketed_amber_lock_addendum.md` — Current controlling age-decision addendum: both age lanes are locked for sensitivity use and must be reported.
- `00_final_board/m33_age_conversion_memo_prelock.md` — Age-conversion analysis memo that produced the candidate values; retained as provenance only. The controlling age document is the bracketed lock addendum.
- `00_final_board/m33_corbelli2014_table1_import.csv` — Imported Corbelli 2014 Table 1 engineering table with R, Vobs, uncertainty, HI surface density, and stellar surface density.
- `00_final_board/m33_corbelli_age_anchor_audit.md` — Initial age-source audit retained as provenance only. The controlling age decision is the bracketed lock addendum.
- `00_final_board/m33_corbelli_board_summary.json` — Compact JSON summary of source-native board status, age bracket, and reconstruction metadata.
- `00_final_board/m33_corbelli_engineering_memo.md` — Format B reconstruction engineering memo documenting source-native board construction from public profiles.
- `00_final_board/m33_corbelli_final_board_freeze_memo.md` — Current final board-freeze memo with bracketed Amber age update and no-run status at board freeze.
- `00_final_board/m33_corbelli_kernel_provenance_review.md` — Kernel/provenance review explaining annular quadrature thick-disk kernel, Casertano-equivalence target, and SPARC non-equivalence.
- `00_final_board/m33_corbelli_kernel_specification.md` — Frozen provisional annular thick-disk kernel specification used to compute Vgas and Vdisk.
- `00_final_board/m33_corbelli_mass_closure_table.csv` — Mass closure checks for stellar mass, H2 mass, HI profile, total gas, and total baryonic mass.
- `00_final_board/m33_corbelli_overlay_benchmark_report.md` — Published-figure overlay benchmark report supporting provisional kernel freeze.
- `00_final_board/m33_corbelli_reconstructed_board.csv` — Final reconstructed board CSV: R, Vobs, errV, Vgas, Vdisk, Vbul, Vbar and provenance fields.
- `00_final_board/m33_corbelli_reconstruction_kernel.py` — Python implementation of the annular thick-disk reconstruction kernel.
- `00_final_board/m33_corbelli_remaining_limitations.md` — Remaining limitations and caveats for the M33-Corbelli source-native reconstruction.

### `01_hermes_result_locked/`
- `01_hermes_result_locked/corbelli_fig12_overlay_reconstruction.png` — Overlay benchmark image carried in the locked Hermes package for traceability.
- `01_hermes_result_locked/m33_age_bracketed_amber_lock.csv` — Age bracket table carried in the locked Hermes execution package.
- `01_hermes_result_locked/m33_baryonic_components_vbar_context.png` — Baryonic components / Vbar context plot from the Hermes execution.
- `01_hermes_result_locked/m33_corbelli_final_board_freeze_memo.md` — Board-freeze memo carried in the locked Hermes execution package.
- `01_hermes_result_locked/m33_corbelli_kernel_specification.md` — Kernel specification carried in the locked Hermes execution package.
- `01_hermes_result_locked/m33_corbelli_reconstructed_board_input.csv` — Exact reconstructed board input used for the locked Hermes execution.
- `01_hermes_result_locked/m33_first_external_hermes_execution_report.md` — Hermes-only execution report under both bracketed Amber age lanes.
- `01_hermes_result_locked/m33_hermes_age_bracket_results.csv` — Hermes result table for t50=6.0±1.2 and 7.1±1.3 Gyr age lanes.
- `01_hermes_result_locked/m33_hermes_age_uncertainty_sensitivity.csv` — Hermes age-uncertainty sensitivity table around the two locked age brackets.
- `01_hermes_result_locked/m33_hermes_chi2_contribution_vs_radius.png` — Hermes radial chi-squared contribution plot.
- `01_hermes_result_locked/m33_hermes_execution_script.py` — Script used for the locked Hermes-only M33 execution.
- `01_hermes_result_locked/m33_hermes_execution_summary.json` — Machine-readable summary of the locked Hermes execution.
- `01_hermes_result_locked/m33_hermes_gate_profile.png` — Hermes gate profile plot for the locked M33 board.
- `01_hermes_result_locked/m33_hermes_per_point_predictions.csv` — Per-radius Hermes predictions, residuals, phi, beta, and diagnostic values for both age lanes.
- `01_hermes_result_locked/m33_hermes_radial_failure_structure.csv` — Hermes radial mask metrics showing main-disk performance and outer-tail underprediction.
- `01_hermes_result_locked/m33_hermes_residuals_vs_radius.png` — Hermes residuals versus radius plot.
- `01_hermes_result_locked/m33_observed_vs_hermes_two_ages.png` — Observed rotation curve vs Hermes predictions for both age brackets.

### `02_mond_comparator/`
- `02_mond_comparator/m33_hermes_vs_mond_residual_comparison.png` — Hermes vs MOND residual comparison plot on the same frozen board.
- `02_mond_comparator/m33_mond_comparator_execution_report.md` — Fixed-simple MOND comparator execution report. No M/L tuning and no a0 fitting.
- `02_mond_comparator/m33_mond_comparator_per_point_predictions.csv` — Per-radius MOND predictions/residuals plus Hermes comparison values.
- `02_mond_comparator/m33_mond_comparator_result_summary.csv` — MOND result summary table with N, chi2, reduced chi2, residual metrics, and classification.
- `02_mond_comparator/m33_mond_comparator_run_config.json` — Machine-readable MOND run configuration: fixed simple interpolation, fixed a0, no M/L tuning.
- `02_mond_comparator/m33_mond_radial_failure_structure.csv` — MOND radial mask metrics and comparison to Hermes radial structure.
- `02_mond_comparator/m33_mond_residuals_vs_radius.png` — MOND residuals versus radius plot.
- `02_mond_comparator/m33_observed_baryons_hermes_mond_context.png` — Observed, baryons, Hermes, and MOND context plot.
- `02_mond_comparator/m33_observed_vs_mond_fixed_simple.png` — Observed rotation curve vs fixed-simple MOND prediction.
- `02_mond_comparator/m33_radial_chi2_contribution_comparison.png` — Radial chi-squared contribution comparison for Hermes and MOND.

### `03_baryonic_convention_sensitivity/`
- `03_baryonic_convention_sensitivity/m33_baryonic_convention_provenance_clarification_addendum.md` — Provenance clarification: f≈1.21 is a SPARC-like diagnostic analogue, not a true SPARC-equivalent unit-M/L board.
- `03_baryonic_convention_sensitivity/m33_baryonic_convention_sensitivity_report.md` — Final M33 baryonic-convention sensitivity diagnostic report with provenance clarification included.
- `03_baryonic_convention_sensitivity/m33_outer_disk_hermes_diagnostic_inversion_primary_age.csv` — Outer-disk diagnostic inversion for Hermes under the primary 6.0 Gyr age lane.
- `03_baryonic_convention_sensitivity/m33_vdisk_scaling_chi2_sensitivity.png` — Chi-squared sensitivity plot as Vdisk is scaled.
- `03_baryonic_convention_sensitivity/m33_vdisk_scaling_dense_sweep_primary_age.csv` — Dense Vdisk scaling sweep for the primary age lane.
- `03_baryonic_convention_sensitivity/m33_vdisk_scaling_outer_residuals.png` — Outer residual behavior versus Vdisk scaling plot.
- `03_baryonic_convention_sensitivity/m33_vdisk_scaling_series_both_ages.csv` — Vdisk scaling series for both locked age brackets, with recomputed g98/beta and Hermes/MOND metrics.

### `04_outer_disk_disequilibrium_audit/`
- `04_outer_disk_disequilibrium_audit/m33_chi2_disequilibrium_zone_overlay.png` — Overlay plot showing radial chi-squared contribution with disequilibrium/outer-disk zones marked.
- `04_outer_disk_disequilibrium_audit/m33_outer_disk_disequilibrium_audit_memo.md` — Outer-disk disequilibrium audit memo mapping Hermes stress zone to warp, break, disturbed HI, and accretion-sensitive features.
- `04_outer_disk_disequilibrium_audit/m33_outer_disk_literature_feature_map.csv` — Radius-feature map for M33 outer-disk complexity literature.
- `04_outer_disk_disequilibrium_audit/m33_outer_disk_mask_metrics.csv` — Radial mask metrics for Hermes and MOND under main-disk and outer-disk cuts.
- `04_outer_disk_disequilibrium_audit/m33_residuals_disequilibrium_zone_overlay.png` — Residual plot with disequilibrium/outer-disk zones marked.

### `05_interpretive_notes/`
- `05_interpretive_notes/m33_origin_trail_archive_lock.md` — Archive classification lock for origin-trail addendum: moderate support for multi-origin gravitational-biography boundary case.
- `05_interpretive_notes/m33_outer_disk_origin_trail_addendum.md` — Origin-trail addendum assessing accretion-sensitive, disturbed, or multi-origin evidence for M33 outer disk/HI.
- `05_interpretive_notes/m33_single_system_age_assumption_note.md` — Interpretive note on M33 and the single-system age assumption / dynamic-aging roadmap.

### `99_archive_locks/`
- `99_archive_locks/m33_corbelli_hermes_only_canon_archive_lock.md` — Canon lock memo preserving the original Hermes-only M33 execution status and checksum.

### Root files
- `README.md` — this file. Lists every included file and its purpose.
- `file_manifest.csv` — manifest with checksums and descriptions.

## Excluded as superseded or duplicate

- Older unbracketed board-freeze ZIPs and earlier reconstruction package ZIPs.
- The unupdated baryonic-convention sensitivity package; the updated provenance-corrected individual files are included instead.
- Nested Hermes/MOND/disequilibrium package ZIPs; their current individual files are included directly.
- The separate age-sweep thought-experiment package, because it was not requested for this consolidated Paper 7 package.
