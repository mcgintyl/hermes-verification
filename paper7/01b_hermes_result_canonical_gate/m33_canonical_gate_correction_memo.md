# M33 canonical-gate correction

**Date:** 2026-07-13

## What this is

The Paper 7 M33 result in `../01_hermes_result_locked/` was produced by
`m33_hermes_execution_script.py`, which contains its own copy of the Hermes gate
function phi(R). That copy deviates in three numerical details from the canonical
gate used for the 133-galaxy SPARC sample in Paper 1 (`paper1/hermes_gate_phi.py`):

| Gate step | Canonical (Paper 1) | M33 script (superseded) |
|---|---|---|
| Knee radius | interpolated crossing of a_knee = 1585 | raw radius of first point below 1585 |
| Savitzky-Golay edge mode | `mirror` | `interp` (SciPy default) |
| Curvature epsilon | `1e-4 * median(|dV/dR|)` | flat `1e-6` |

Paper 7's stated methodology is a frozen chassis: the same Hermes implementation
used for SPARC, with nothing changed for external boards. The three deviations
above are inconsistent with that statement. They are not required by the Format B
data (they operate on the already-reconstructed R and g_bar arrays, not on the
board reconstruction). They are incidental differences in a separately-typed
execution script.

## The fix

This directory regenerates the M33 Hermes result using the **canonical** gate
(`paper1/hermes_gate_phi.py`), so a single frozen gate is used throughout the
project. The M33 baryonic board, ages, beta values (beta_primary = 1.83421,
beta_conservative = 1.69876), error model (errV^2 + sigma_int^2, sigma_int^2 =
386), and the MOND comparator are all unchanged. Only phi(R) changes.

## Effect on results

Point-by-point phi changes by up to 0.625, but the effect on every reported
quantity is small and crosses no interpretive threshold.

| Quantity | Locked (superseded gate) | Canonical gate |
|---|---:|---:|
| Full-board chi2nu, primary t50 = 6.0 | 3.13 | 3.22 |
| Full-board chi2nu, conservative t50 = 7.1 | 3.30 | 3.39 |
| Main disk chi2nu, R <= 10 (primary) | 0.80 | 0.93 |
| Main disk chi2nu, R <= 8 (primary) | 0.64 | 0.78 |
| Outer tail chi2nu, R > 10 (primary) | 5.79 | 5.86 |
| Outer share of total chi2 | 86.3% | 84.6% |
| Outer mean residual (primary) | -51.6 km/s | -51.9 km/s |
| MOND full-board chi2nu | 0.12 | 0.12 (unchanged) |

Conclusions are unchanged: Hermes ports non-catastrophically, the inner disk is
clean (chi2nu < 1), the residual concentrates in the outer disk (~85% of chi2),
MOND fits cleanly and wins the board, and convention scaling narrows but does not
eliminate the outer-disk residual. Full per-mask and convention-sweep values are
in `m33_hermes_mask_summary_canonical.csv` and `m33_hermes_per_point_predictions_canonical.csv`.

## Provenance

The original locked archive in `../01_hermes_result_locked/` is retained
unchanged as the historical record of what was run. This directory supersedes its
phi/chi2 values for all reporting purposes. The `phi_superseded_locked` column in
the per-point CSV carries the original values for direct comparison.
