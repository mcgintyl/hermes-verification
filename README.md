# Hermes Verification Repository

Independent verification tools for two papers by McGinty (2026):

| Folder | Paper | What it verifies |
|--------|-------|-----------------|
| `paper1/` | *The Hermes Equation: Galaxy Rotation Curves from Age with No Free Constants* | Rotation curve fits + age derivations for 133 galaxies |
| `paper2/` | *Testing Gravity at Encounter Speed: Geometric Prediction of Earth Flyby Anomalies* | Geometric scores + sign predictions for 12 Earth flybys |
| `docs/` | Supplementary materials | Age derivation audit trail (151 methods, 77 sources) |

## Requirements

- Python 3.8+
- numpy >= 1.20 (Paper 1 only)
- scipy >= 1.7 (Paper 1 only)
- Paper 2 uses only the standard library (math, csv)

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
```

The SPARC rotation curve data used in Paper 1 is archived on Zenodo:

```
Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models
for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.
doi:10.5281/zenodo.16284118
```

## License

MIT
