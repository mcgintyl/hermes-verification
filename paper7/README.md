# Paper 7 reproduction package

*Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing*
(McGinty 2026). Zenodo, all versions: 10.5281/zenodo.20550091.

One command rebuilds the M33 board from Corbelli et al. (2014)'s Table 1 and checks every
M33 number of version 4 of the paper: Tables 2, 2b and 3, Appendix A and the prose.

```bash
python paper7/reproduce.py
```

Nothing is downloaded: every input is in this repository. The script needs Python 3.10 or
later (the frozen kernel's type annotations use the `X | Y` syntax), NumPy, SciPy and
pandas (`requirements.txt`; the frozen kernel reads its table with pandas) and runs in
about ten seconds. CI runs it on every push.

It prints one line per check with the place in version 4 of the paper that reports the
number, the value it got, the printed value and the tolerance:
- "Sec 4 para 6" is the sixth paragraph of Section 4.
- "(support)" marks a check behind a statement whose number is not printed.
- "[INFO]" marks a value printed for information only (C.1); it has no pass/fail and is
  not counted as a check.

It writes `paper7/results/` (git-ignored) and exits 1 if any check fails. A printed number
passes when the recomputed value is within half a unit of its last printed digit; ties
are allowed.

Version 4 is in preparation. Until it is deposited, `paper7_main.md` is the version 3
text (see Version notes).

## What it checks

| Group | What |
|---|---|
| H, shipped files | `MANIFEST_SHA256.txt` matches every file and lists them all (operating-system and editor files such as `.DS_Store`, `Thumbs.db`, `desktop.ini`, `*.pyc` and `*.swp` are ignored). `file_manifest.csv` agrees with it. SHA-256 of the Table 1 import, the board, the kernel, the age lock and the canonical per-point file (Paper 8's `reproduce.py` and `verify_gates.py` read it). |
| B, board | Vgas and Vdisk rebuilt from Table 1 through the frozen kernel, against the board CSV, to floating-point precision (about 1e-14 km/s; checked to 1e-9). The 58 radii from 0.24 to 22.72 kpc, Vbulge = 0, and the kernel constants and distance that Appendix A states (840 kpc, the 01b script's `D0_KPC`). Mass closure (support). |
| G, gate and scoring | The log-grid gate against the canonical phi of v3. g98 and beta. The canonical per-point predictions. The correction from 3.13 to 3.22 (Sec 2.1). a0 (Sec 2.5). The floor statements of Sec 2.7. |
| T2, Table 2 | All 24 rows: N, R range, chi2/N, chi2 fraction, mean, median and RMS residual. |
| T2b, Table 2b | All 24 cells of the floor sweep (sigma_int^2 = 0, 100, 386, 900). |
| P, prose | Abstract and Sections 4 and 6. This covers no-floor values, the inner overshoot and mid-disk deficit against errV alone, the age bracket, the unit-M/L analogue and the Paper 8 refinement (P.37, P.38). Section 6's "neither alone closes the residual" (support): with Paper 8's operator, the density release alone leaves the R > 10 chi2/N at 2.60, and a younger outer age alone (0 Gyr beyond 10 kpc) at 4.12 (P.39, P.40); both stay above 1. |
| T3, Table 3 | The minimum of the Vdisk scaling sweep at both ages. The 11 printed rows, the M/L columns and the caption's claim. |
| A, Appendix A | The 3-dp table reproduces 3.22, 3.39, 0.93, 5.86 and MOND 0.12 and 0.23 (A.2 to A.7), every chi2/N cell of Tables 2 and 3 (A.11) and every prose number that P.1 to P.38 check except the floor-free Hermes 118 (A.12), each at its printed precision. It does not reproduce the floor-free Hermes cells of Table 2b: it gives 118.51 and 87.54 where the table prints 118.48 and 87.49 (A.13, A.14). Those, and some residual and chi2-fraction cells, need the full-precision board CSV. Every cell of the Appendix A in `paper7_main.md` equals the board at its printed precision, and a version 4 text must print Vgas and Vdisk to 3 dp (A.10). |
| C, chain-rule gate | Paper 1's chain-rule gate (`paper9/hermes_clusters/gate.py`) gives the same chi2/N to 2 dp in all 38 Hermes chi2/N cells of Tables 2 and 3 (Sec 2.1). C.1 prints the largest difference between the two gates' phi on the board (0.0045) for information. |
| S, SPARC context | The SPARC medians and fractions that Sec 4 para 5 quotes, from the shipped Paper 1 table (Hermes) and Paper 8 table (MOND, a0 3702.813). |
| F, shipped tables | Each regenerated table in `results/` against `expected/`, cell by cell: numbers to 1e-9, relative or absolute, and text exactly. |

Not rebuilt:
- **Section 5's numbers** (the Repetto case and the SPARC-overlap THINGS and dwarf runs).
  Their inputs were extracted from published figures and are not in this repository;
  version 4 quotes them from the original runs.
- **Tables 1 and 1b**, a record of the source audit with no computed values.

## Conventions

Each of these changes a checked value if altered.

- **Board:** g_bar = (Vdisk^2 + Vbulge^2 + sign(Vgas) Vgas^2) / R. g98 is the 98th percentile of the board's g_bar: 2653.950 (km/s)^2/kpc.
- **Gate:** the log-grid gate in `paper1/hermes_gate_phi.py`, imported, not copied. The chain-rule gate of Paper 1's published scores is checked separately (group C). The two gates are described in `docs/gate_version_history.md`.
- **Hermes:** beta = pi exp(-t50 g98 / 46654) - 1/sqrt(2 pi), at t50 = 6.0 Gyr (primary) and 7.1 Gyr (conservative). V = sqrt(max(g_bar (1 + beta phi), 0) R).
- **MOND:** simple interpolating function. a0 = 1.2e-10 m/s^2 with the exact parsec, 3702.813 (km/s)^2/kpc.
- **chi2/N:** the mean over the N points of (Vobs - V)^2 / (errV^2 + 386). N is not reduced.
- **Residual:** V_model - Vobs.
- **Table 3:**
  - Vdisk is multiplied by the scale, and g98, beta and phi are recomputed at each scale.
  - The M/L factor is scale^2. The global M*/L at 3.6 um is 0.682647 scale^2, where 0.682647 = 4.9e9 / 10^(0.4 (24.60 + 3.24 - 3.2)): Corbelli's stellar mass, total [3.6] = 3.2, distance modulus 24.60, solar 3.24.
  - The sweep runs from 0.5 to 3.0 in steps of 0.0025. Within 0.005 of each age's coarse minimum it runs in steps of 0.0001.
  - Each printed row is evaluated at its printed scale. The minimum row is evaluated at 1.368, the 3-dp rounding of the fine-grid minimum at 1.3677. At 1.3677 the row's chi2/N cells are the same, but its global M*/L, Hermes outer mean and MOND outer mean would print 1.277, -42.58 and 5.12 instead of 1.278, -42.57 and 5.13.
  - The unit-M/L analogue is Vdisk x 1.210 (sqrt(1/0.683)).

## Files

| File | What it is |
|---|---|
| `reproduce.py` | Rebuilds and checks the paper's M33 numbers |
| `expected/` | The tables version 4 prints, at full precision (see below) |
| `MANIFEST_SHA256.txt` | SHA-256 of every file in this folder except itself, `results/`, `__pycache__/` and operating-system or editor files; `reproduce.py` checks it (H.1, H.2) |
| `file_manifest.csv` | Path, bytes, SHA-256 and a description of every file except itself and `MANIFEST_SHA256.txt`. Kept because the paper's Data Availability cites it for the board's SHA-256 (H.10); `reproduce.py` checks it agrees with `MANIFEST_SHA256.txt` (H.3) |
| `make_manifest.py` | Refreshes both manifests (descriptions are kept in `file_manifest.csv`; a new file needs a row first); `--check` compares without writing |
| `paper7_main.md` | The version 3 text, the source of the version 3 PDF on Zenodo (10.5281/zenodo.21347359) |
| `.gitignore` | Keeps `results/` out of git |

`expected/`:

- `table2_masks.csv`: Table 2, the 24 rows.
- `table2b_floor_sweep.csv`: Table 2b, the four floors.
- `table3_vdisk_sweep.csv`: the full canonical-gate sweep behind Table 3, 1,193 scales.
- `table3_rows_printed.csv`: the 11 rows Table 3 prints, each at its printed scale (the minimum row at 1.368, the 3-dp rounding of 1.3677).
- `appendix_a_board_3dp.csv`: Appendix A as version 4 prints it, with Vgas and Vdisk to 3 dp.
- `prose_numbers.csv`: every number checked outside the tables, with its location, full-precision value, printed value and tolerance.

## Folder map

| Folder | What it holds | Status |
|---|---|---|
| `00_final_board/` | The inputs (see below) and the board-freeze memos | Current |
| `01_hermes_result_locked/` | The first Hermes run, from May 2026 | Superseded gate; historical record (`README_SUPERSEDED_GATE.md`) |
| `01b_hermes_result_canonical_gate/` | The M33 Hermes run on the canonical gate, from Paper 7 v3 onwards | Current |
| `02_mond_comparator/` | The MOND run on the same board | MOND current; Hermes comparison columns superseded |
| `03_baryonic_convention_sensitivity/` | The first Vdisk-scaling sweep and its report | Hermes values superseded; MOND and M/L columns current |
| `04_outer_disk_disequilibrium_audit/` | Radial masks and the outer-disk literature map | Hermes values superseded; MOND and the literature map current |
| `05_interpretive_notes/`, `99_archive_locks/` | Interpretive notes and the archive lock of the first run | Historical; their Hermes numbers are the superseded gate's (for example 3.126) |

`00_final_board/` holds:
- Corbelli et al. (2014) Table 1 as imported (`m33_corbelli2014_table1_import.csv`);
- the frozen thick-disk kernel (`m33_corbelli_reconstruction_kernel.py`);
- the 58-radius board (`m33_corbelli_reconstructed_board.csv`);
- the age lock (`m33_age_bracketed_amber_lock.csv`: t50 = 6.0 +/- 1.2 and 7.1 +/- 1.3 Gyr);
- the mass-closure table.

Two things to keep in mind about `00_final_board/`:
- The controlling age decision is `m33_age_bracketed_amber_lock_addendum.md` with the age-lock CSV. Files named `*prelock*` or `*age_anchor_audit*` are provenance only: they record how the bracket was derived.
- The board is source-native, not SPARC-equivalent: Corbelli's pixel-SED stellar mass, helium factor 1.33, an H2 prescription, a flaring stellar disk and no bulge. The Vdisk x 1.210 case is a SPARC-like diagnostic analogue, not a SPARC board.

`01b_hermes_result_canonical_gate/` holds:
- `compute_m33_hermes_canonical.py`;
- the per-point predictions;
- the mask summary;
- the correction memo;
- the radius-uncertainty note.

The current values of every superseded column are in `expected/` and `01b_`. The first run's gate was a separately coded copy that differed from Paper 1's in three details:
- the knee radius;
- the Savitzky-Golay edge mode;
- the curvature regularisation.

Two things about `01b_` that `reproduce.py` makes explicit:
- The per-point file was written with beta rounded to 5 dp (1.83421 and 1.69876, as the correction memo quotes). Its V_Hermes values therefore differ from the full-precision closed form by up to 7e-5 km/s (G.6, G.7). No printed number changes.
- The correction memo's last "Effect on results" sentence says the convention-sweep values are in its two CSV files. They are not. The canonical-gate sweep is `expected/table3_vdisk_sweep.csv`.

The memo is left as written; its bytes are part of the v3 record. Do not change the per-point file: `paper8/reproduce.py` (M.0a) and `verify_gates.py` read it. Do not change the board either: it is the hashed input of every M33 number (H.5), and the paper's Data Availability gives its SHA-256 (H.10).

## Version notes

**Versions 1 and 2** printed the M33 numbers of the first run, on the superseded gate copy:
- full board 3.13 and 3.30;
- R <= 10 kpc 0.80;
- R > 10 kpc 5.79;
- outer share 86%;
- a best continuous scaling of 2.38 and a unit-M/L analogue of "approximately 2.5".

Version 1 (May 2026) was not deposited. Zenodo holds version 2 (10.5281/zenodo.20550092, 2026-06-05) and version 3.

**Version 3** (10.5281/zenodo.21347359, 2026-07-13) corrected the gate:
- Every M33 Hermes number was recomputed on Paper 1's gate (`01b_`): full board 3.22 and 3.39, R <= 10 kpc 0.93, R > 10 kpc 5.86, outer share 85%. MOND was unchanged (0.12).
- It added Appendix A, the board with Vgas and Vdisk printed to 2 dp.

The correction was not carried through Table 3 and its paragraph:
- v3 kept the 1.325 row, the minimum of the superseded sweep on its 0.025 grid, labelled as the lowest Hermes scaling, and printed its canonical value, 2.45, as "the best continuous stellar scaling".
- It kept "approximately 2.5" for the unit-M/L analogue.
- In that row, v3 printed the global M*/L as 1.199 (recomputed: 1.198) and the MOND outer mean as 3.43 (recomputed: 3.42498).

Separately, Appendix A at 2 dp gives 3.38, not 3.39, at 7.1 Gyr and 0.92, not 0.93, for R <= 10 kpc (A.8, A.9). Every other M33 number in v3 reproduces on the canonical gate, including all 168 numeric cells of Table 2.

**Version 4** (in preparation) is checked here. Table 2 and the cited numbers (3.22, 3.39, 0.93, 5.86, 85%, MOND 0.12 and 0.23) are unchanged. It adds or corrects:
- **Table 2b and Sec 2.7, the error floor.** sigma_int^2 = 386 (km/s)^2 exceeds errV^2 at all 58 radii. Without it, Hermes scores 118 and MOND 2.88.
- **Table 3.** The 1.325 row is replaced by the canonical minimum, 2.41 at Vdisk x 1.368 (t50 = 6.0 Gyr; 2.48 at x 1.415 for 7.1 Gyr). The row is evaluated at its printed scale, 1.368, the 3-dp rounding of the fine-grid minimum at 1.3677 (see Conventions). The unit-M/L analogue becomes 2.55, with MOND 0.19.
- **Sec 4.** New statements against errV alone: Hermes runs 11-17 km/s high at 0.46-1.45 kpc and 19-35 km/s low at 4.6-9.5 kpc, and 61% of its chi2 lies beyond 10 kpc. Also the +/-1.2 Gyr age range: full board 3.05-3.40, R > 10 kpc 5.51-6.20.
- **Sec 6.** A cross-check with the Paper 8 refinement: its density release combined with a younger effective age for the outer gas lowers the R > 10 kpc chi2/N from 5.86 to 1.30-1.53 (P.37, P.38), and neither alone closes the residual (P.39, P.40, support).
- **Appendix A** with Vgas and Vdisk to 3 dp, so the printed table reproduces 3.39 and 0.93. At its printed precision it reproduces every chi2/N cell of Tables 2 and 3 and every prose number except the floor-free 118 of Sec 4 (A.11, A.12). The floor-free Hermes values of Table 2b (118.48 and 87.49) and some residual and chi2-fraction cells need the CSV: from the 3-dp table the two floor-free values come out 118.51 and 87.54 (A.13, A.14).

Each version's paper text is committed here after its Zenodo deposit. `paper7_main.md` (SHA-256 aa02551b...) is the version 3 text, byte for byte the source of the deposited PDF. The version 4 text will replace it after the v4 deposit. When it does, first edit the `paper7_main.md` description in `file_manifest.csv` (it says "version 3"; `make_manifest.py` keeps descriptions as they are), then refresh the manifests with `python paper7/make_manifest.py`.

Check A.10 parses Appendix A from whichever text is shipped and compares every cell with the board at its printed precision:
- the v3 text gives 2 dp for Vgas and Vdisk, and passes;
- a v4 text (recognised by the words "Version 4" anywhere in it, as its version note has) must give 3 dp, or A.10 fails.

The places in the check lines refer to the version 4 text.
