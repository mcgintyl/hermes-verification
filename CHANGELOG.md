# Changelog

All notable changes to this repository are recorded here: what changed, when,
and why. Entries are newest first.

## 2026-08-16

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
