# Changelog

All notable changes to this repository are recorded here: what changed, when,
and why. Entries are newest first.

## 2026-08-16

### Added

- `paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv` — new `psi` column holding the
  pre-computed wear variable, appended as the last column so existing column
  positions are unchanged. Computed as

      psi = t50_gyr * g98 / 46654

  where 46654 is c/(2π) expressed in Gyr·(km/s)²/kpc — the compound units that
  `t50_gyr × g98` already carries. Verified to reproduce the shipped `beta_eff`
  column via β = π·exp(−ψ) − 1/√(2π) to a maximum absolute difference of
  5.6 × 10⁻¹⁶ across all 133 galaxies. No existing cell was modified.

  **Why:** recomputing ψ from the symbolic form ψ = 2π·t₅₀·g₉₈/c invites
  substituting the SI speed of light (299,792,458 m/s) for c. That is ~1023×
  too large for these units, which drives ψ to ~0 and collapses β onto its
  ceiling of π − 1/√(2π) = 2.7427 for every galaxy. An external reviewer hit
  exactly this, despite the unit conversion being stated in the paper. Shipping
  ψ removes the need to recompute it.

  Note the divisor is 46654 rather than the rounded c = 293,136 used with a
  factor of 2π; the rounded form costs ~1.1 × 10⁻⁶ in β and would no longer
  reproduce `beta_eff` exactly.

### Changed

- `docs/stage1_calibration_set.md` — corrected the T1/T2/T3 terminology to match
  Paper 1 v6. The three classes are **age-state classes**, not data-quality or
  data-reliability classes; a simple cutoff on t₅₀ (≈4.05 and ≈8.25 Gyr)
  reproduces 131 of the 133 assignments. The `Tier` column name in
  `Table_S1_Age_Methods.csv` and the `T1`/`T2`/`T3` labels are unchanged, since
  they are the join key between the paper and the supplementary tables.
  (commit `717e97a`)

### Notes

- Paper 1 was revised to v6 and republished on Zenodo under the concept DOI
  `10.5281/zenodo.18809176` (age-state class terminology; removal of a broken
  Appendix A cross-reference; font and subtitle colour). No data, equations, or
  results changed. The repository hosts verification code and data only — paper
  PDFs are not mirrored here, so no repository change was required.
