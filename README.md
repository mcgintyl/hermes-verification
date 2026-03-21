# Hermes Verification Repository

Independent verification tools for four papers by McGinty (2026):

| Folder | Paper | What it verifies |
|--------|-------|-----------------|
| `paper1/` | *The Hermes Equation: Galaxy Rotation Curves from Age with No Free Constants* | Rotation curve fits + age derivations for 133 galaxies |
| `paper2/` | *Testing Gravity at Encounter Speed: Geometric Prediction of Earth Flyby Anomalies* | Geometric scores + sign predictions for 12 Earth flybys |
| `paper4/` | *The Emergent Plane at Cosmological Scale: Chirality Prediction and Exclusion Limits* | Internal consistency of 10 supplementary data files, 185 checks |
| `paper5/` | *Weak Lensing Pilot: Age-Dependent Shear Signal (KiDS x GAMA)* | Pipeline output consistency, 142 checks across 8 test groups |
| `docs/` | Supplementary materials | Age derivation audit trail (151 methods, 77 sources) |

## Requirements

- Python 3.8+
- numpy >= 1.20 (Paper 1 only)
- scipy >= 1.7 (Paper 1 only)
- Papers 2, 4, 5 use only the standard library (math, csv, re)

```
pip install -r requirements.txt
```

---

## Paper 1 — The Hermes Equation (`paper1/`)

### The Equation

The model predicts galaxy rotation curves from two inputs per galaxy: stellar age and peak baryonic acceleration.

**Age score:**

    beta = pi * exp(-2*pi * t50 * g98 / c) - 1 / sqrt(2*pi)

- `t50` — half-mass stellar age in Gyr
- `g98` — 98th-percentile baryonic acceleration in (km/s)^2/kpc
- `c` — speed of light (in matching units: c/(2*pi) = 46,654 Gyr*(km/s)^2/kpc)

**Model acceleration at each radius R:**

    g_model(R) = g_bar(R) * [1 + beta * phi(R)]

The gate function phi(R) opens in low-acceleration outer regions and closes in dense inner regions. All gate parameters are globally fixed — no per-galaxy tuning.

### Data

1. **SPARC rotation curve files** — the `Rotmod_LTG` folder from the [SPARC database](https://zenodo.org/records/16284118) ([![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16284118.svg)](https://doi.org/10.5281/zenodo.16284118))
2. **Galaxy age table** — included as `paper1/ages_133.csv`

### Rotation Curve Verification

```
python paper1/verify_hermes.py --sparc /path/to/Rotmod_LTG --ages paper1/ages_133.csv
python paper1/verify_hermes.py --sparc /path/to/Rotmod_LTG --ages paper1/ages_133.csv --csv results.csv
```

**Expected output:**

```
Hermes Equation Verification  (133 galaxies, 0 missing)
===============================================================================================
Galaxy             t50        g98     beta    N    chi2_H   chi2_M  winner
-----------------------------------------------------------------------------------------------
CamB               4.1      319.0   2.6558    9     0.607    0.997  Hermes
...
-----------------------------------------------------------------------------------------------
Median chi2_nu                  Hermes:    1.323   MOND:    1.143
Trimmed mean (5%)               Hermes:    2.498   MOND:    3.295
Hermes wins                     68 / 133 (51.1%)
Below chi2=5.0                  Hermes: 112 (84.2%)   MOND: 95 (71.4%)
```

### Age Derivation Verification

`paper1/verify_ages.py` independently re-derives the stellar half-mass age (t₅₀) for each galaxy using three documented conversion paths. Every intermediate step is printed for auditing.

**Path 1 — Cumulative SFH Interpolation** (Resolved CMD / SFH)

    t50 = t_young + (F_young - 0.5) / (F_young - F_old) * (t_old - t_young)

7 galaxies verified — **100% match**.

**Path 2 — sSFR / Birthrate Parameter** (Hα, UV, IR/radio SFR)

    b = sSFR × T_Hubble        (T_Hubble = 13.8 Gyr)

Then map b → t₅₀ using the audit's empirical calibration from 34 anchor galaxies. 30 galaxies verified — **87% match**.

**Path 3 — Broadband Color → SPS Lookup** (B-V, B-R, FUV-NUV)

Interpolate observed colors against calibration tables built from the audit's documented BC03 conversions. 29 galaxies verified — **90% match**.

**Overall: 91% match** (73/80 verified within ±0.5 Gyr). Remaining mismatches trace to documented multi-step protocol overrides.

```
python paper1/verify_ages.py                          # all 137 galaxies, full steps
python paper1/verify_ages.py --quiet                  # compact table only
python paper1/verify_ages.py --galaxy "NGC 55"        # single galaxy with steps
python paper1/verify_ages.py --path 1                 # only cumulative SFH path
python paper1/verify_ages.py --csv age_results.csv    # write results to CSV
```

### Paper 1 Files

| File | Description |
|---|---|
| `paper1/verify_hermes.py` | Rotation curve verification — reads SPARC data, computes Hermes + MOND fits |
| `paper1/verify_ages.py` | Age derivation verification — re-derives t₅₀ via three conversion paths |
| `paper1/hermes_gate_phi.py` | Standalone gate function phi(R), verified to floating-point precision |
| `paper1/ages_133.csv` | Age table for 133 galaxies (galaxy, t50 in Gyr, g98 in (km/s)^2/kpc) |

---

## Paper 2 — Flyby Anomaly Geometric Score (`paper2/`)

### The Score

A geometric scoring function classifies all 12 documented Earth flybys (1990–2013) using three fixed parameters and no mission-specific tuning:

    S = |Δcos δ| × exp(−h/H) × (V₀/V∞)^p

**Fixed parameters:**
- `H = 2500 km` — altitude scale height
- `V₀ = 10 km/s` — reference velocity
- `p = 1.0` — velocity exponent

**Inputs per flyby:**
- `δ_in, δ_out` — incoming/outgoing asymptotic declinations (degrees)
- `h` — perigee altitude above Earth's surface (km)
- `V∞` — hyperbolic excess speed (km/s)

**Sign prediction:** The sign of `cos(δ_in) − cos(δ_out)` predicts the direction of the anomaly (speed increase vs decrease).

### Results

- **11/12 flybys correctly classified** (anomaly vs null). Every confirmed null is predicted null; every confirmed anomaly is predicted anomaly.
- **5/5 sign predictions correct** for all detected anomalies (including low-confidence Cassini).
- **Juno paradox:** Juno scores higher than Galileo I but showed no anomaly — explained by the GRACE absorption hypothesis (see paper Section 5).
- **All 12 computed scores match published values** to floating-point precision.

### Usage

```
python paper2/verify_flyby.py                      # full verification, all 12 flybys
python paper2/verify_flyby.py --mission NEAR        # single mission with detailed steps
python paper2/verify_flyby.py --p 0.5               # sensitivity check with p=0.5
python paper2/verify_flyby.py --csv results.csv     # write results to CSV
```

### Paper 2 Files

| File | Description |
|---|---|
| `paper2/verify_flyby.py` | Flyby score verification — computes S for all 12 flybys, checks against published table |
| `paper2/flyby_data.csv` | Input data for 12 Earth flybys (Section 3.3 of the paper), every value source-traced |

---

## Paper 4 — Chirality at Cosmological Scale (`paper4/`)

### What It Tests

The chirality paper tested 273,055 spin-labeled galaxies against 26,111 galaxy clusters and 82,458 galaxy groups for chirality coherence in dense environments. The paper's conclusion: no signal at any scale. The verification tool checks that all 10 supplementary data files are internally consistent and that every number in the paper's tables matches the data.

### Checks Performed (185 total)

1. **Cross-match counts** — ALFALFA membership totals self-consistent (31,502 sources, 39 deep-void-interior)
2. **Histogram integrity** — bin counts sum to reported N for each void membership group
3. **Summary statistics** — N values consistent across files
4. **Cluster chirality** — CW + CCW = N_galaxies for every cut; all p-values > 0.01 (null result)
5. **Excess variance** — chi-squared degrees of freedom match cluster counts; chi2/df ratios near 1.0
6. **Per-void summary** — 117 deepest voids; source counts consistent with cross-match totals
7. **Dark candidates** — zero candidates found (header-only file, as documented)
8. **Paper Table 1** — all published galaxy/system counts match CSV data exactly
9. **Binomial p-values** — recomputed from CW/N using normal approximation; all consistent

### Usage

```
python paper4/verify_chirality.py              # run all 185 checks
python paper4/verify_chirality.py --verbose    # show every intermediate step
```

### Paper 4 Files

| File | Description |
|---|---|
| `paper4/verify_chirality.py` | Verification tool — 9 check groups, 185 individual checks |
| `paper4/paper4_step2_crossmatch_counts.csv` | ALFALFA x VAST membership counts |
| `paper4/paper4_hi_property_summary_stats.csv` | HI property summary statistics by void membership |
| `paper4/paper4_hist_logMHI.csv` | Histogram: log HI mass by void group |
| `paper4/paper4_hist_W50.csv` | Histogram: velocity width by void group |
| `paper4/paper4_hist_logMHI_over_Ms.csv` | Histogram: gas fraction by void group |
| `paper4/paper4_per_void_summary.csv` | Per-void summary for 117 deepest voids |
| `paper4/paper4_dark_candidates.csv` | Dark galaxy candidates (zero found) |
| `paper4/paper4_pure_cluster_member_cuts_summary.csv` | Chirality tests per cluster membership cut |
| `paper4/paper4_excess_variance_summary_by_richness.csv` | Excess variance by richness bin and radius |

---

## Paper 5 — Weak Lensing Pilot (`paper5/`)

### What It Tests

The lensing pilot measured galaxy-galaxy lensing signals around 9,004 GAMA lenses using KiDS DR4 source shapes, split by stellar age (Dn4000) to test whether younger stellar populations produce stronger shear at fixed mass. The verification tool checks that all pipeline output files are internally consistent.

### Checks Performed (142 total)

1. **Summary slopes (Table 1)** — significance = |slope|/error recomputed; 6 per-R radial bins present per mass bin
2. **Lens counts across files** — N per (mass_bin, Dn4000 slice) consistent across group summary, median table, and profile files (60 cross-checks)
3. **Combined beta** — inverse-variance weighted mean recomputed from 6 per-bin slopes
4. **Cross-shear null** — all 6 mass bins have cross-component slopes consistent with zero (< 3 sigma)
5. **Diagnostic summary** — 9,004 lenses, 21,162,005 sources, S/N = 7.2, null test chi2/dof values match
6. **Profile completeness** — every (mass_bin, slice) combination has exactly 12 radial bins
7. **Total lens count** — 9,004 in per-lens file; quintile slices balanced within each mass bin
8. **Stacked signal** — Young DS = +7.36, Old DS = +4.62, ratio = 1.59, difference = 0.7 sigma

### Usage

```
python paper5/verify_lensing.py              # run all 142 checks
python paper5/verify_lensing.py --verbose    # show every intermediate step
```

### Paper 5 Files

| File | Description |
|---|---|
| `paper5/verify_lensing.py` | Verification tool — 8 check groups, 142 individual checks |
| `paper5/FULL_summary_slopes_by_massbin.csv` | Table 1: per-mass-bin slope of DeltaSigma vs Dn4000 |
| `paper5/FULL_summary_slopesX_by_massbin.csv` | Cross-shear (45-deg) slopes — null test |
| `paper5/FULL_fits_slope_by_R_massbin.csv` | Per-radial-bin slopes for each mass bin |
| `paper5/FULL_fits_slopeX_by_R_massbin.csv` | Per-radial-bin cross-shear slopes |
| `paper5/FULL_group_summary_massbin_dn4000_q5.csv` | Group summary: N lenses per (mass, Dn4000 slice) |
| `paper5/FULL_dn4000_medians_by_massbin_slice.csv` | Dn4000 median/mean per (mass, slice) |
| `paper5/FULL_deltasigma_profiles_massbin_dn4000_q5.csv` | Full DeltaSigma profiles (12 R bins x 30 groups) |
| `paper5/FULL_lenses_with_massbin_dn4000slice_q5.csv` | Per-lens catalog (9,004 lenses with all properties) |
| `paper5/combined_stacked_summary.txt` | Stacked young vs old DeltaSigma result |
| `paper5/diagnostic_summary.txt` | Pipeline diagnostic: null tests, S/N, configuration |

---

## Supplementary Documentation (`docs/`)

| File | Description |
|---|---|
| `docs/supp_methods_age_derivation_structured.md` | Full age derivation procedures for all 133 galaxies (151 methods, 77 published sources) |
| `docs/galaxy_age_method_index.csv` | Per-galaxy index: method class, tier, confidence, and source citations |

---

## Citation

```
McGinty, L. A. (2026). The Hermes Equation: Galaxy Rotation Curves from Age
with No Free Constants.

McGinty, L. A. (2026). Testing Gravity at Encounter Speed: Geometric
Prediction of Earth Flyby Anomalies and the GRACE Absorption Hypothesis.

McGinty, L. A. (2026). The Emergent Plane at Cosmological Scale: A Geometric
Chirality Prediction, Empirical Exclusion Limits, and a Falsifiable Roadmap.

McGinty, L. A. (2026). Weak Lensing Pilot: Age-Dependent Shear Signal
in the Emergent Plane Framework (KiDS x GAMA).
```

The SPARC rotation curve data used in Paper 1 is archived on Zenodo:

```
Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models
for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.
doi:10.5281/zenodo.16284118
```

## License

MIT
