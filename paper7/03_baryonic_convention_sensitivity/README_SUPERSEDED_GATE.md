# Read this before using the Hermes values in this folder

This folder is the **historical record** of the first Vdisk-scaling (stellar
mass-to-light convention) sweep. Every **Hermes** value in it was computed with
the superseded gate copy of the first M33 run (see
`../01_hermes_result_locked/README_SUPERSEDED_GATE.md`):

- the Hermes columns of `m33_vdisk_scaling_series_both_ages.csv`,
  `m33_vdisk_scaling_dense_sweep_primary_age.csv` and
  `m33_outer_disk_hermes_diagnostic_inversion_primary_age.csv`;
- the gate-output columns: `phi_last` and `r_knee` of
  `m33_vdisk_scaling_series_both_ages.csv`, and `phi` of the inversion file (the
  knee radius was one of the three gate deviations);
- the Hermes numbers in `m33_baryonic_convention_sensitivity_report.md`, including
  its minimum (chi2nu 2.3813 at Vdisk scale f = 1.331) and its baseline (3.126);
- the plots.

On the canonical gate the lowest Hermes chi2nu on the scaling is **2.41 at
Vdisk x 1.368** (t50 = 6.0 Gyr; 2.48 at x 1.415 for 7.1 Gyr), and the SPARC-like
unit-M/L analogue (x 1.210) gives **2.55** (2.67 at 7.1 Gyr). Paper 7 v3
recomputed Table 3 on the canonical gate but kept the scale of the superseded
sweep's minimum on its 0.025 grid (1.325), so it printed 2.45 as the lowest value,
and it kept the analogue's superseded value ("approximately 2.5"). Version 4
corrects both.

Still current: the **MOND** columns (MOND does not use the gate); the scale,
M/L factor and global M*/L columns (M*/L(3.6 um) = 0.683 of the source-native
disk); the recomputed g98 and beta; and the provenance clarification that the
x 1.21 branch is a SPARC-like diagnostic analogue, not a SPARC-equivalent board
(`m33_baryonic_convention_provenance_clarification_addendum.md`).

**For all reporting and reproduction purposes, use:**

- `../expected/table3_vdisk_sweep.csv`: the canonical-gate sweep, Vdisk scale 0.5
  to 3.0 in steps of 0.0025, refined to 0.0001 within 0.005 of each age's minimum;
- `../expected/table3_rows_printed.csv`: the rows Table 3 prints;
- `python paper7/reproduce.py`, which rebuilds and checks both.
