# Read this before using the Hermes values in this folder

This folder is the **historical record** of the MOND comparator run on the M33
board. Its **MOND** values are current: MOND does not use the Hermes gate, and the
board, a0 = 3702.813 (km/s)^2/kpc, the simple interpolating function and the
errV^2 + 386 scoring are the ones every version of Paper 7 reports (MOND full
board 0.12, R > 10 kpc 0.23).

Its **Hermes** values are not current. The Hermes comparison columns
(`Vhermes_*`, `residual_model_minus_obs_hermes_*` and `chi2_contribution_hermes_*`
in `m33_mond_comparator_per_point_predictions.csv`), the Hermes rows of
`m33_mond_comparator_result_summary.csv` and `m33_mond_radial_failure_structure.csv`,
the Hermes numbers in `m33_mond_comparator_execution_report.md` (for example a full
board of 3.126 and an outer share of 86.3%) and the Hermes curves in the plots were
copied from the first Hermes run, whose gate copy was superseded (see
`../01_hermes_result_locked/README_SUPERSEDED_GATE.md`).

**For all reporting and reproduction purposes, use:**

- `../01b_hermes_result_canonical_gate/` for the canonical-gate Hermes predictions
  (full board 3.22 and 3.39, not 3.13 and 3.30);
- `../expected/table2_masks.csv` for Table 2, both models, every mask and residual
  column, on the canonical gate;
- `python paper7/reproduce.py`, which rebuilds and checks both.
