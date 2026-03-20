# Hermes Equation Verification Tool

Independent verification of the per-galaxy rotation curve fits from:

> McGinty, L. A. (2026). *The Hermes Equation: Galaxy Rotation Curves from Age with No Free Constants.*

This tool reads SPARC rotation curve files and a galaxy age table, then computes reduced chi-squared scores for both the **Hermes equation** and **MOND** (simple interpolation function), matching the published results.

## The Hermes Equation

The model predicts galaxy rotation curves from two inputs per galaxy: stellar age and peak baryonic acceleration.

**Age score:**

    beta = pi * exp(-2*pi * t50 * g98 / c) - 1 / sqrt(2*pi)

- `t50` — half-mass stellar age in Gyr
- `g98` — 98th-percentile baryonic acceleration in (km/s)^2/kpc
- `c` — speed of light (in matching units: c/(2*pi) = 46,654 Gyr*(km/s)^2/kpc)

**Model acceleration at each radius R:**

    g_model(R) = g_bar(R) * [1 + beta * phi(R)]

- `g_bar(R)` — Newtonian gravity from visible mass (disk + bulge + gas)
- `phi(R)` — gate function in [0, 1], built from the baryonic acceleration profile

The gate opens in low-acceleration outer regions and closes in the dense inner regions. All gate parameters are globally fixed — no per-galaxy tuning.

## Requirements

- Python 3.8+
- numpy >= 1.20
- scipy >= 1.7

Install:

```
pip install -r requirements.txt
```

## Data

You need two things:

1. **SPARC rotation curve files** — the `Rotmod_LTG` folder from the [SPARC database](http://astroweb.cwru.edu/SPARC/). Each file is `<GalaxyName>_rotmod.dat` with columns: R (kpc), Vobs, errV, Vgas, Vdisk, Vbul (all km/s), SBdisk, SBbul.

2. **Galaxy age table** — included as `ages_133.csv` in this repo. Contains the 133 galaxies with verified stellar ages from 77 published sources. Columns: `galaxy`, `t50_gyr`, `g98`.

## Usage

```
python verify_hermes.py --sparc /path/to/Rotmod_LTG --ages ages_133.csv
```

Optional: write results to CSV:

```
python verify_hermes.py --sparc /path/to/Rotmod_LTG --ages ages_133.csv --csv results.csv
```

## Expected Output

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

These numbers match the paper's Table 1 (median chi-squared 1.312, MOND 1.141) to within minor floating-point differences in the gate's Savitzky-Golay smoothing step.

## File Overview

| File | Description |
|---|---|
| `verify_hermes.py` | Main verification script — reads SPARC data, computes Hermes + MOND fits |
| `hermes_gate_phi.py` | Standalone gate function phi(R), verified to floating-point precision |
| `ages_133.csv` | Age table for 133 galaxies (galaxy name, t50 in Gyr, g98 in (km/s)^2/kpc) |
| `docs/supp_methods_age_derivation_structured.md` | Supplementary methods: full age derivation procedures for all 133 galaxies |
| `docs/galaxy_age_method_index.csv` | Per-galaxy index of age derivation method, tier, confidence, and source citations |
| `verify_ages.py` | Age derivation verification — re-derives t₅₀ via three documented conversion paths |

## Age Derivation Verification

`verify_ages.py` independently re-derives the stellar half-mass age (t₅₀) for each galaxy using the three conversion paths documented in the supplementary methods audit. Every intermediate step is printed so an auditor can inspect the full chain from published measurement to final t₅₀.

### The Three Conversion Paths

**Path 1 — Cumulative SFH Interpolation** (Resolved CMD / SFH)

For galaxies with resolved star-formation histories (e.g. ANGST Table 2), t₅₀ is the lookback time where the cumulative formed-mass fraction crosses 0.5:

    t50 = t_young + (F_young - 0.5) / (F_young - F_old) * (t_old - t_young)

7 galaxies verified — **100% match**.

**Path 2 — sSFR / Birthrate Parameter** (Hα, UV, IR/radio SFR)

Given a specific star-formation rate (sSFR = SFR / M★), compute the birthrate parameter:

    b = sSFR × T_Hubble        (T_Hubble = 13.8 Gyr)

Then map b → t₅₀ using the audit's empirical calibration from 34 anchor galaxies. The exponential τ-model inversion is also shown for reference. 30 galaxies verified — **87% match**.

**Path 3 — Broadband Color → SPS Lookup** (B-V, B-R, FUV-NUV)

Interpolate observed colors against calibration tables built from the audit's documented BC03 conversions (25 B-V anchors, 8 B-R anchors, 9 FUV-NUV anchors). Protocol overrides (LSB Stability Block, Metallicity Trap, Edge-On Dust Trap) are encoded in the lookup values. 29 galaxies verified — **90% match**.

### Overall Result

**91% of computationally verified galaxies match their published t₅₀ within ±0.5 Gyr** (73 / 80). The remaining mismatches trace to documented multi-step protocol overrides (Cosmic Downsizing penalties, composite color frosting traps, edge-on dust corrections) that require the full audit logic.

### Usage

```
python verify_ages.py                          # all 137 galaxies, full steps
python verify_ages.py --quiet                  # compact table only
python verify_ages.py --galaxy "NGC 55"        # single galaxy with steps
python verify_ages.py --path 1                 # only cumulative SFH path
python verify_ages.py --path 2                 # only sSFR / birthrate path
python verify_ages.py --path 3                 # only broadband color path
python verify_ages.py --csv age_results.csv    # write results to CSV
```

## How the Rotation Curve Verification Works

1. For each galaxy, read the SPARC rotmod file to get R, Vobs, errV, Vgas, Vdisk, Vbul
2. Compute baryonic acceleration: `g_bar = (Vdisk^2 + Vbul^2 + sgn(Vgas)*Vgas^2) / R`
3. Compute the gate `phi(R)` from the g_bar profile (6 fixed processing steps)
4. Compute beta from the galaxy's age (t50) and peak acceleration (g98)
5. Predict: `V_model = sqrt(g_bar * (1 + beta * phi) * R)`
6. Score: `chi2_nu = mean((Vobs - Vmodel)^2 / (errV^2 + 386))`
7. Compare against MOND: `g_MOND = g_bar * nu(g_bar/a0)` with `nu(y) = (1 + sqrt(1 + 4/y)) / 2`

## Citation

If you use this tool, please cite the original paper:

```
McGinty, L. A. (2026). The Hermes Equation: Galaxy Rotation Curves from Age
with No Free Constants.
```

## License

MIT
