# Changelog

All notable changes to this repository are recorded here: what changed, when,
and why. Entries are newest first.

## 2026-09-12

### Added

Paper 9 is published (Zenodo concept DOI 10.5281/zenodo.22719297) and its
reproduction package now lives in the repository. **No existing code logic, no
expected value and no published number changed.**

- `paper9/` — **new.** The Paper 9 reproduction package, byte-identical to the
  archive published with the paper and verified against its own
  `MANIFEST_SHA256.txt` (47 entries, no mismatches, no unlisted files).
  `reproduce_clusters.py` rebuilds Table 2 from the CLASH lensing profiles and
  baryon models and checks every cluster number the paper quotes (24 checks);
  `reproduce_galaxies.py` rebuilds Table 1 from the SPARC rotation curves
  (10 checks). SPARC is fetched and checksum-verified by `fetch_sparc.py`, not
  redistributed.
- `README.md` — a `paper9/` row in the folder table and a `## Paper 9` section
  with the commands, the file inventory, and the two conventions that travel with
  the package.
- `.github/workflows/verify.yml` — `paper9` added to the compile check, plus
  two steps running both halves. The galaxy step reuses the SPARC copy the
  workflow already downloads (`fetch_sparc.py --from-dir "$SPARC_DIR"`) instead of
  fetching from astroweb on every run.
- `.gitignore` — excludes `paper9/data/sparc/`, `paper9/results/` and the
  downloaded archive, so a local run cannot accidentally commit redistributed
  SPARC data.

### Documented

- `paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv` is now annotated in the
  README: its `chi2nu_mond` column was computed with **a0 = 3700 exactly**, not
  3702.8. Verified to 3.6e-15 across all 133 galaxies, while both 3702.x
  conversions miss all 133 of them. Annotated rather than recomputed, because it
  is a published artifact.
- `README.md` and `verify_gates.py` — a scope note on the word "canonical".
  Paper 9 ships the **chain-rule** gate, `paper9/hermes_clusters/gate.py`, the
  convention that reproduces Paper 1's published *score* column, while
  `paper1/hermes_gate_phi.py` reproduces the published *phi* values. Adding a
  second gate to the repository made the earlier wording, "used for every result
  in the project", literally false.

  **Why:** the 1.312/1.323 question has already cost several rounds, and the
  published export carries one column from each convention. With two gate files
  in the repository the distinction has to be stated where someone will actually
  read it, not only in `docs/`.

## 2026-09-06

### Documented

Documentation only. **No code logic, no expected values, and no numerical
results changed.** The intrinsic variance floor `σ_int² = 386 (km/s)²` was
present in the code with only the bare comment "intrinsic variance floor" and no
statement of its origin anywhere in the repository. An external reviewer
reproducing the pipeline annotated it in their own code as "ACCORDING TO THE
AUTHOR," which is the correct response to an undocumented constant.

- `docs/chi2_error_floor.md` — **new.** Provenance, both distortions, the
  Appendix A sensitivity table, and how to reproduce every figure from the
  repository's own code.
- `paper1/verify_hermes.py` — provenance and consequences recorded above
  `SIGMA_INT_SQ`. The value is unchanged at `386.0`.
- `README.md` — new "A note on the χ² error floor" block beside the existing
  1.312/1.323 note, linking to the docs file.

  **Why:** 386 is not arbitrary — it is SPARC's own ~0.12 dex
  radial-acceleration-relation scatter carried into velocity space. At fixed R,
  `g ∝ V²/R`, so 0.12 dex in `g` is 0.06 dex in `V`, giving `dV/V = 0.148`;
  `19.65 / 0.148` implies a characteristic velocity of 132.6 km/s against SPARC's
  median Vobs of 134.0 km/s. The floor is the data source's own reported scatter
  at the sample's own median velocity, adopted in preference to inventing a
  Hermes-specific one.

  Two properties are now stated explicitly, both measured over 133 galaxies /
  3073 points using the repository's own functions:

  | property | measurement |
  |---|---|
  | floor exceeds `errV²` | 2975/3073 points (96.8%) |
  | median share of variance from the floor | 94.9% |
  | Spearman(χ² at floor 0, χ² at floor 386) | **0.71** — the floor reorders, not just rescales |
  | Spearman(χ² at floor 386, Vmax) | **0.705** — the score partly tracks galaxy mass |
  | floor as % of peak Vobs, dwarfs / giants | 34.0% / 8.7% |

  The second property is the one with consequences for readers: a fixed χ²
  cutoff used to grade fit quality will partly encode galaxy size. 97% of dwarfs
  but 23% of giants fall below χ² = 2, and 13 of the 14 galaxies above χ² = 8 are
  massive spirals. Neither property affects the Hermes-vs-MOND comparison, since
  the same floor is applied to both models.

## 2026-09-01

### Corrected

Documentation only. **No code logic, no expected values, and no numerical
results changed.** Three files carried an explanation of the 1.312 ↔ 1.323
baseline-median gap that has since been refuted.

- `paper1/verify_appendix_a.py` — the comment above `TOL_CHI2_MEDIAN` said the
  tolerance was widened "to absorb documented platform floating-point drift
  across NumPy/SciPy versions." Replaced with the derivative-lineage
  explanation. `TOL_CHI2_MEDIAN` stays at `0.15`; only the justification changed.
- `README.md` — replaced the "tolerances absorb platform floating-point drift"
  sentence with an explicit note on the two shear conventions.
- `docs/gate_version_history.md` — added the **third** gate lineage (the Paper 1
  scoring gate) and corrected a misleading bullet.

  **Why:** the old explanation was wrong. The 1.312 ↔ 1.323 gap is not version
  noise — it is two algebraically identical but numerically different forms of
  the gate's shear term on SPARC's non-uniform radial grid:

      published Paper 1 (1.3115) :  s = |(R/V) · dV/dR|     chain rule, linear grid
      current repo gate (1.3226) :  s = |d ln V / d ln R|   log grid

  Verified for this entry by independent re-implementation:

  | check | result |
  |---|---|
  | chain-rule shear vs published `chi2nu_configg` | **133/133**, max 3.6e-13, median 1.311534725 |
  | log-grid shear vs published `chi2nu_configg` | 0/133, max 1.48e-01, median 1.322593 |
  | log-grid shear vs published `phi_last` | **133/133**, max 4.2e-13 |
  | chain-rule shear vs published `phi_last` | 1/133, max 5.3e-02 |
  | scores, signed vs clipped `g_bar` | difference exactly **0.0** (score-degenerate) |

  Two consequences worth recording. First, the drift explanation was refuted by
  experiment, not argument: the same verifier under the original production
  environment (Python 3.13.2 / numpy 2.4.3 / scipy 1.17.1) is bit-identical to
  the pinned one, 0.0 across all 133 galaxies. Second, the perfect
  complementarity above means the published export joins **two gate lineages** —
  `phi_last` from the log-grid gate, `chi2nu_configg` from the chain-rule gate.
  A `verify_gates.py` pass confirms gate outputs (`phi_last`) only and must not
  be read as reproducing Paper 1's scores; that inference is what allowed the
  mismatch to persist. The negative-`g_bar` clipping convention cannot be
  recovered from the scores, since signed and clipped are numerically identical.

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
