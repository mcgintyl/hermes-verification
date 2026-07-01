# Paper 10 — Two Layers of Galaxy Aging in MaNGA DynPop

**Full title:** *Two Layers of Galaxy Aging in MaNGA DynPop: A Principle-Level Proof of Concept*

**Zenodo DOI:** [10.5281/zenodo.21088606](https://doi.org/10.5281/zenodo.21088606)

**Primary dataset:** MaNGA DynPop DR17 public catalogs
**Claim level:** proof of concept / signal documentation, not a detection claim

## What this package contains

This is the reproducibility and verification package for Paper 10. It tests a principle-level
prediction — that at fixed stellar mass, younger galaxies show larger mass-to-light discrepancies
between dynamical and stellar-population estimates — using three joined public MaNGA DynPop DR17
products (JAM dynamical catalogs, stellar-population/star-formation-history catalogs, and
circular-velocity-curve tables). After quality cuts, 5,952 galaxies are split into eight equal-count
stellar-mass bins, and the youngest and oldest quartiles are compared within each bin.

The result is a **two-layer age structure**:

- **Layer 1 (visible SPS layer).** The raw mass-to-light discrepancy (DML) is larger for young
  galaxies in 8/8 mass bins (observed SPS) and 6/8 bins (intrinsic SPS). A component decomposition
  shows this is mostly driven by the SPS denominator — older stellar populations have higher stellar
  mass-to-light ratios, exactly as standard stellar-population synthesis predicts. This is not a discovery.
- **Layer 2 (controlled dynamical residual).** Under a collinearity-free strict-control model that
  removes SPS mass-to-light, metallicity, and structural variables, a smaller dynamical residual
  persists in **7/8 mass bins** across all three SPS control sets. Parametric NFW/gNFW dark-matter
  fractions do not show the same signal, and the strongest residual tracks the dust-carrying SPS
  definition. The paper presents this as a proof of concept and an invitation to specialist
  replication, not a detection claim.

The main paper's primary controlled-residual result is the **strict-control specification** (7/8 bins).
The original kitchen-sink specification (8/8 bins) is retained as supplementary sensitivity only.

## Folder map

- `manga_two_layers_paper10_final.md` — the main paper.
- `manga_two_layers_paper10_supplementary_appendix.md` — detailed methodology, literature passes,
  control-variable specifications, limitation analyses, and sensitivity documentation.
- `00_strict_control_rerun/` — **primary result.** Strict-control rerun outputs: the collinearity-free
  reduced-control model, controlled dynamical-residual scoreboards and per-bin tables, model-term and
  age-coefficient tables, the reduced-spec sensitivity grid, the handoff note, and the rerun pipeline script.
- `01_clean_pipeline_v0_6/` — original clean-pipeline (v0.6) outputs: the DML sledgehammer report and
  scoreboards, the M/L split / controlled-residual hardening report (kitchen-sink specification),
  figures, the simulation/assembly-bias literature pass, and the working draft with revision notes.
- `file_manifest.csv` — machine-readable manifest with path, byte size, SHA256, and description for every file.

## `00_strict_control_rerun/` — primary controlled-residual result

The strict-control specification removes the v0.7 collinearity issue (including raw `Sigma_Re` and
`logSigma_Re` together, and duplicate metallicity and flattening terms). The shared controls are
mass-weighted metallicity (`sp_MW_Metal_Re`), MGE ellipticity (`Eps_MGE`), and `logSigma_Re`, plus
stellar mass, size, Sérsic index, and `Lambda_Re`. Young galaxies retain higher controlled dynamical
M/L residuals in **7/8 mass bins** for the intrinsic, observed, and both-SPS control sets; the
magnitude test is at the shuffle floor and the count-only p is ~0.03–0.04. One bin (bin 5,
median logM ≈ 10.663) is the lone non-decisive bin in all three sets.

- `manga_paper10_strict_control_primary_summary.csv` — primary young-vs-old scoreboard for the three SPS control sets.
- `manga_paper10_strict_control_residual_scoreboard.csv` — controlled-residual scoreboard with bootstrap intervals and shuffle-null statistics.
- `manga_paper10_strict_control_residual_perbin.csv` — per-bin controlled residuals across the eight mass bins.
- `manga_paper10_strict_control_model_terms.csv` — standardized T50 coefficient, robust SE/t/p, R², cross-validated ΔR² from age, and Spearman ρ.
- `manga_paper10_strict_control_age_coefficients_rawY.csv` — age (T50) coefficients on the raw dynamical-M/L outcome per 1 SD.
- `manga_paper10_strict_control_sensitivity_grid.csv` — reduced-spec sensitivity grid; all 12 variants give 7/8 bins.
- `manga_paper10_strict_control_claude_handoff.md` — handoff note documenting the specification choice and main-text recommendation.
- `run_clean_controls_fast.py` — the strict-control rerun pipeline script.

## `01_clean_pipeline_v0_6/` — original clean-pipeline outputs

- `manga_sledgehammer_final_report.md`, `manga_sledgehammer_final_scoreboard.csv`, `manga_sledgehammer_final_bin_results.csv` — primary DML sledgehammer split (8/8 observed-SPS, 6/8 intrinsic-SPS) and the secondary fDM checks.
- `manga_ml_split_hardening_final_report.md`, `manga_ml_split_hardening_scoreboard.csv`, `manga_ml_split_hardening_control_residual_scoreboard.csv`, `manga_ml_split_hardening_per_bin.csv` — M/L component split and the kitchen-sink controlled-residual result (8/8; supplementary sensitivity).
- `manga_simulation_assembly_bias_pass.md` — simulation/assembly-bias literature pass: no exact mock-MaNGA/JAM/SPS comparator located; the residual is unbenchmarked against cosmological mocks, not proven anomalous relative to ΛCDM.
- `manga_two_layers_clean_pipeline_draft_v0_6.md`, `manga_two_layers_v0_6_revision_notes.md`, `manga_two_layers_v0_5_to_v0_6.diff` — working draft, revision notes, and diff (retained as provenance; superseded by the final paper).
- `manga_sledgehammer_final_*.png`, `manga_ml_split_hardening_*.png` — figures (kitchen-table, mass-adjusted, bin-difference, decomposition, controlled-residual plots).
- `manifest_sha256.txt` — SHA256 checksums for the v0.6 bundle as originally released.

## Guardrails

- DML is a JAM/SPS mass-to-light construction, **not** a board-complete observed rotation curve. This
  is a principle-level test, not a strict test of the Hermes rotation-curve equation. Treating a
  JAM-inferred circular-velocity curve as a SPARC-style observed rotation curve would collapse the
  board-completeness discipline established in prior work.
- The raw DML signal is mostly SPS-driven and is not clean gravitational evidence. Only the controlled
  dynamical residual is the interesting layer, and it is modest, one-bin-sensitive, and could still be
  produced by IMF variation, assembly bias, dust conventions, cold gas, or JAM covariance.
- No cosmological-mock comparator has been run. The residual is documented for independent groups to
  test, refine, or falsify — it is not claimed to be new physics.

See the main paper and supplementary appendix for the full result, caveats, and hand-off.
