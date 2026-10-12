# Changelog

All notable changes to this repository are recorded here: what changed, when,
and why. Entries are newest first.

## 2026-10-12

### Added

`paper7/` is now a reproduction package for Paper 7, the field audit (Zenodo concept DOI
10.5281/zenodo.20550091), built for version 4. Before this, the folder held the M33
board and the run records but nothing that checked a printed number. The canonical-gate
convention sweep behind Table 3 existed nowhere in the repository.

- `paper7/reproduce.py` rebuilds and checks every M33 number of version 4 (178 checks;
  exit 1 on any failure). It needs no download and Python 3.10 or later. It covers:
  - the board, rebuilt from Corbelli et al. (2014) Table 1 through the frozen kernel
    to floating-point precision (about 1e-14 km/s; checked to 1e-9), and the 840 kpc
    distance Appendix A states;
  - the gate against the canonical phi of v3, and the scoring statements;
  - Table 2 (24 rows), Table 2b (the floor sweep) and Table 3 (a Vdisk scaling sweep
    with its minimum at both ages);
  - Appendix A at 3 dp: it reproduces every chi2/N cell of Tables 2 and 3 and every
    prose number except the floor-free Hermes values of Table 2b (118.51 and 87.54
    against the printed 118.48 and 87.49), which need the CSV. A version 4 text must
    print its Appendix A at 3 dp;
  - every M33 number in the prose, including Section 6's Paper 8 refinement and, as
    support, that neither of its two parts alone closes the outer residual;
  - the chain-rule gate of Paper 1's published scores, which gives the same chi2/N to
    2 dp in all 38 Hermes cells of Tables 2 and 3;
  - the shipped hashes and tables. Operating-system and editor files (`.DS_Store`,
    `Thumbs.db`, `desktop.ini`, `*.pyc`, `*.swp`) are ignored.

  It runs in CI.
- `paper7/expected/` holds the tables version 4 prints, at full precision:
  - Table 2, Table 2b, the full Table 3 sweep and the 11 printed Table 3 rows (each
    at its printed scale; the minimum row at 1.368, the 3-dp rounding of the fine-grid
    minimum 1.3677);
  - Appendix A with Vgas and Vdisk to 3 dp;
  - every checked prose number.
- `paper7/MANIFEST_SHA256.txt` (Paper 8 style), checked by `reproduce.py`, and
  `paper7/make_manifest.py`, which refreshes it together with `file_manifest.csv`.
- `paper7/paper7_main.md`: the version 3 text as deposited (10.5281/zenodo.21347359).
  The version 4 text replaces it after its deposit.
- `README_SUPERSEDED_GATE.md` notes in `02_mond_comparator/`,
  `03_baryonic_convention_sensitivity/` and `04_outer_disk_disequilibrium_audit/`.
  Their Hermes values come from the first run's superseded gate copy; their MOND values
  are current.
- `paper7/.gitignore`: keeps `paper7/results/` out of git.

### Changed

- `paper7/README.md` is rewritten. It was the May 2026 package README and did not
  mention `01b_hermes_result_canonical_gate/`. It now gives:
  - the command and the checks;
  - the conventions;
  - a folder map that marks the canonical and superseded folders;
  - version notes for v1 to v4.
- `paper7/file_manifest.csv` now lists every file. Before, it omitted
  `01b_hermes_result_canonical_gate/` (the source of v3's numbers),
  `01_hermes_result_locked/README_SUPERSEDED_GATE.md` and `README.md`. Rows whose
  Hermes values are on the superseded gate say so.
- `.gitattributes`: `paper7/** -text`, so the hashed files keep their bytes.
- `.github/workflows/verify.yml`: compile-checks `paper7` and runs
  `paper7/reproduce.py`. The header comment now matches what CI runs.
- `README.md` (root): CI now runs the Paper 7 checks. The Paper 7 rows name
  `reproduce.py`, `01b_hermes_result_canonical_gate/`, `expected/`,
  `MANIFEST_SHA256.txt` and the Zenodo concept DOI, and mark the Hermes values of
  `01_` to `04_` as the superseded gate's, pointing to `paper7/README.md`. The
  "frozen data package, no reruns" sentence is replaced: `reproduce.py` reruns
  everything.

### Documented

Found while building the package. No committed data file or deposited text changes;
version 4 corrects the printed v3 numbers below.

- Version 3's Table 3 kept the superseded sweep's minimum scale (1.325) as "the lowest
  Hermes scaling", so Section 4 printed 2.45. On the canonical gate the minimum is 2.41 at
  Vdisk x 1.368, and the unit-M/L analogue is 2.55, not "approximately 2.5".
  Version 4 corrects both.
- The 01b per-point file was written with beta rounded to 5 dp (1.83421, 1.69876).
  Its V_Hermes values differ from the full-precision closed form by up to 7e-5 km/s.
  The file is unchanged: Paper 8 and `verify_gates.py` read its phi.

## 2026-10-11

### Changed

`paper8/` now ships the version 3 paper texts. Version 3 was deposited on Zenodo on
2026-10-10 (version DOI 10.5281/zenodo.23291143, concept DOI 10.5281/zenodo.20699990).
Each version's texts are committed after its Zenodo deposit.

- `paper8/paper8_main.md` and `paper8/paper8_supplementary_materials.md` are now the
  version 3 texts, the sources of the deposited PDFs.
  - Both pin commit d3355f4; this commit changes none of its code or data.
  - In the supplement, the six multiplication asterisks in one Supplement E sentence are
    escaped. Unescaped, Markdown pairs them into italics: the first PDF export printed
    `phiSWpiexp(-psi)` for `phi*S*W*pi*exp(-psi)`. The deposited PDF prints them correctly.
- `paper8/README.md`: the file list and version notes now describe the version 3 texts
  instead of saying version 3 is being prepared.
- `paper8/MANIFEST_SHA256.txt`: regenerated for the changed files.
- `README.md`: the Paper 8 file table now says version 3.

## 2026-10-07

### Added

`paper8/` is now a reproduction package for Paper 8 (Zenodo concept DOI
10.5281/zenodo.20699990). Before this, the folder held only the operator, and
`verify_eta.py` said it "reproduces the Paper 8 results" while its advertised
`--sparc_dir/--age_table` options were ignored and nothing was scored.

- `paper8/reproduce.py` rebuilds the paper's SPARC and M33 numbers from public inputs
  and checks them (176 checks; exit 1 on any failure). Candidates 1 and 2 are not
  rebuilt. It covers:
  - Table 5.1;
  - the chain-rule gate comparison;
  - stellar M/L of 0.5 for disks and 0.7 for bulges;
  - candidates 3 to 5 and other closed-form constants;
  - Supplements B and C;
  - M33 on the gate as corrected in Paper 7 v3, plus MOND on the same board;
  - the shipped CSV files.

  It runs in CI.
- `paper8/expected/` holds the paper's ten production CSV files, public for the first
  time (the Zenodo records of versions 1 and 2 hold only the PDFs):
  - the seven SPARC tables as produced in May 2026;
  - the three M33 tables, regenerated by `reproduce.py` on the corrected gate. The
    production copies were computed on Paper 7's superseded gate.
- `paper8/README.md`, `paper8/MANIFEST_SHA256.txt`, and the paper texts
  (`paper8_main.md`, `paper8_supplementary_materials.md`, as published in version 2).

### Fixed

- `paper8/verify_eta.py`:
  - The docstring no longer claims reproduction or a command-line interface.
  - Velocities use V = sqrt(max(g, 0) R), as in Paper 1's published per-galaxy
    export. The old signed form changed UGC 01281 only.
  - MOND uses a0 = 1.2e-10 m/s^2 with the exact parsec (3702.813, was 3703.0). No
    printed number changes.
  - "BTFR-safe" is dropped, unused constants are removed, and pi is listed as
    inherited.
- `README.md`:
  - The Paper 8 row and section give the concept DOI; they gave the version-1 DOI
    20699991, or none.
  - The "BTFR-safe" and "all constants are derived" wording, which the paper
    does not support, is replaced.
  - The section names the gate Paper 8 uses.
  - The opening's "Every script runs in CI" now lists what CI actually runs.
- `.github/workflows/verify.yml`: runs `paper8/reproduce.py`; the header comment now
  matches what CI runs.
- `.gitattributes`: `paper8/** -text`, so the hashed files keep their bytes.
  `.gitignore`: `paper8/results/`.

## 2026-10-01

### Fixed

Found in a cold-run audit of Paper 9 on 2026-09-30. **No computation, no expected
value and no published number changed:** `reproduce_clusters.py` 24/24 with Table 2
still byte-identical (`47ddd6b7…`), `reproduce_galaxies.py` 10/10.

- `requirements.txt` — added `pandas==2.3.1`. CI had failed on every commit since
  `paper9/` was added (71e896c, 07bc908): the workflow installs only this file, and
  the Paper 9 scripts import pandas.
- `README.md` — the 2026-09-12 entry below and the Paper 9 section said `paper9/` is
  byte-identical to a package archived or published with the paper. It is not: the
  Zenodo record (10.5281/zenodo.22719297) holds the paper PDFs only, and the
  reproduction code and data are published here, in `paper9/`. The section now says
  so. The opening "Nothing here is fitted" now names the one exception, Paper 9's
  per-cluster amplitude K. The Paper 10 DOI is now the concept DOI (21088605); the
  version DOI it replaced (21088606) pointed at the superseded v1.
- `paper9/` — package v2:
  - `README.md` states the SPARC input as Paper 1 does, the default normalisation
    Υ* = 1.0 used as given for both models, instead of "the published mass-to-light
    ratios"; adds a short version note.
  - `hermes_clusters/galaxies.py` — `reproduce_galaxies.py` run before
    `fetch_sparc.py` now prints its "Run: python fetch_sparc.py" message instead
    of a `FileNotFoundError` traceback.
  - `reproduce_clusters.py` — the K = 1 reference line is labelled chi2 per node
    (it divides by 51 nodes, not by the 44 degrees of freedom).
  - `MANIFEST_SHA256.txt` regenerated for the three changed files.

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
