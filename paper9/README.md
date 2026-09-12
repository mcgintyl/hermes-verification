# Paper 9 reproduction package

This package reproduces the results of *The Hermes Equation: Cluster Lensing from a Frozen Galaxy Gate* (L. A. McGinty, Paper 9 of the Hermes series) from the published input data: both the galaxy baseline of Table 1 and the cluster results of Table 2. Appendix A of the paper describes the cluster method; this package is the code and data behind both halves.

## Run it

```
pip install -r requirements.txt
python reproduce_clusters.py
```

For the galaxy half, fetch the SPARC rotation curves first (they are not redistributed here; see *Data and licences* below):

```
python fetch_sparc.py
python reproduce_galaxies.py
```

Each script takes a few seconds, prints one line per check, and exits with status 0 only if every check passes. Between them they write five files to `results/`:

- `table2_seven.csv`: Table 2 of the paper. It should be byte-identical to `expected/table2_seven.csv` (SHA-256 `47ddd6b7867018f8b42b55a422d768205381384471c0240f918eb285e16bb96c`).
- `boards.csv`: per-cluster details (nodes, radial range, baryonic accelerations, nodes beyond the X-ray data, scores).
- `checks.txt`: every cluster check, with the computed values.
- `table1_galaxies.csv`: the per-galaxy scores behind Table 1, for all 133 galaxies.
- `checks_galaxies.txt`: every galaxy check, with the computed values.

## What it reproduces

| where in the paper | what |
|---|---|
| Table 1 | the galaxy baseline on 133 SPARC galaxies and 3,073 points: median χ²/N (Hermes 1 1.311535, MOND 1.142686), the 5% trimmed means (2.495, 3.295), the median \|ΔV\| pooled over points (20.6, 25.8 km/s) and per galaxy (18.9, 20.2 km/s), the catastrophes with χ²/N > 10 (12, 15), and the head-to-head record (68 to 65) |
| §3 | that the two are statistically indistinguishable on those galaxies: paired Wilcoxon on log χ² (p = 0.9365) and an exact sign test on the head-to-head record (p = 0.8624) |
| Table 2 | χ², χ²/dof and the amplitude K for Hermes 1 and MOND on each of the seven relaxed clusters; the totals (Hermes 1 71.132/44 = 1.617, MOND 101.787/44 = 2.313); median and mean K; Hermes 1 winning 5 of 7 |
| §4 | the scores of the two excluded clusters quoted in the text (MACSJ0744, RXJ1347); the composition of the baryonic carrier |
| §5 | the A611 worked example; the K ranges; the 1/r⁴ gas alternative (Hermes 1 1.96, MOND 1.45); one K for all clusters; each cluster predicted with K from the others; K calibrated on three clusters and applied blind to the other four; the test of whether K varies between clusters; the baryon-closure estimates (3–8%, 74%, and the outermost-node baryon fractions of 0.44–0.76) |
| §6 | identical χ² at t₅₀ = 5, 10 and 13 Gyr |
| Appendix A | node admission (70 lensing nodes, 51 scored) and X-ray coverage (33 of the 51 nodes beyond R_max,X) |

`reproduce_galaxies.py` also checks all 133 per-galaxy Hermes 1 scores against the published Paper 1 values in `expected/paper1_per_galaxy.csv`; they agree to 4 × 10⁻¹³.

## Layout

```
reproduce_clusters.py   the cluster half: runs everything and checks it against the paper
fetch_sparc.py          downloads the SPARC rotation curves and verifies their checksums
reproduce_galaxies.py   the galaxy half: runs everything and checks it against the paper
hermes_clusters/
  gate.py               the Hermes 1 gate phi(R), Paper 1's chain-rule construction
  boards.py             baryonic carrier and lensing-node admission
  model.py              Hermes 1 and MOND responses at cluster scale; the GLS amplitude fit
  galaxies.py           rotation-curve reader, Hermes 1 and MOND velocities, the chi^2 score
data/
  famaey/               baryonic inputs (Famaey, Pizzuti & Saltas)
  mistele/              lensing profiles and correlation matrices (Mistele et al.)
  rmax_x.csv            radial extent of the X-ray data (Mistele et al. Table 1)
  ages_133.csv          t50 and g98 for the 133 galaxies
  sparc_manifest.csv    the file name, SHA-256 and size of each rotation curve used
  sparc/                the rotation curves themselves, once fetch_sparc.py has run
  SOURCES.md            citations, licences and checksums
expected/               the published Table 2, the per-galaxy Paper 1 scores, and
                        full-precision reference values
```

## Method in brief

At galaxy scale, for each of the 133 galaxies:

1. The baryonic acceleration is g_bar = [V_disk² + V_bul² + sign(V_gas) V_gas²] / R, from the SPARC rotation-curve file, at the published mass-to-light ratios.
2. The gate φ(R) is computed from g_bar by the same chain-rule construction used at cluster scale.
3. Hermes 1 predicts V = √(g_bar [1 + β φ] R), with β = π e^(−ψ) − 1/√(2π) and ψ = t₅₀ g₉₈ / 46654.
4. MOND uses the simple interpolation function, ν(y) = ½ [1 + √(1 + 4/y)] with y = g_bar / a₀.
5. Both are scored as χ²/N = mean of (V_obs − V_model)² / (errV² + σ_int²), with σ_int² = 386 (km/s)². This is a mean over points, not a reduced χ² over degrees of freedom: neither model has a free parameter per galaxy here.

At cluster scale, for each cluster:

1. The baryonic carrier, M_bar = M_BCG+companions + M_gas (1 + f_gal), is computed on Famaey et al.'s 50-point radial grid (0.05 Mpc to 1.2 r₂₀₀).
2. Mistele et al.'s lensing nodes that lie within that grid are admitted. For RXJ2129, the innermost three nodes are also left out, as Mistele et al. warn they may not be reliable.
3. At each node, M_bar is interpolated log-log from the grid.
4. The gate φ is computed on the grid from g_bar = G M_bar / r², then interpolated to the nodes.
5. The Hermes 1 prediction is M_pred = M_bar (1 + K q), with q = φ (π e^(−ψ) − 1/√(2π)) and ψ = t₅₀ g₉₈ / 46654, where g₉₈ is the 98th percentile of g_bar over the nodes and t₅₀ = 10 Gyr.
6. MOND uses the simple interpolation function, ν = ½ [1 + √(1 + 4 a₀ / g_bar)], so M_pred = M_bar (1 + K (ν − 1)).
7. Each model gets one amplitude K ≥ 0 per cluster, fitted by generalised least squares against the covariance built from Mistele et al.'s correlation matrix and statistical errors. This costs one degree of freedom per cluster.

`boards.py` is copied verbatim from the board-build script used for the paper. `gate.py` is the same gate function that produces Paper 1's galaxy results, and both halves of this package call it.

## Assumptions to be aware of

- **t₅₀ = 10 Gyr for every cluster.** It has no effect on χ², because K absorbs ψ exactly (cluster check 3).
- **The 3–8% baryon-closure estimate assumes that stars make up 10–15% of a cluster's baryons.** This is a benchmark taken from the literature, not a measurement from these data. The carrier's own stellar share is about 7.5%.
- **The gas beyond R_max,X is the best-fit profile continued outward, as in Famaey et al.'s code.** Mistele et al.'s alternative, a 1/r⁴ tail from R_max,X, reverses the Hermes 1 vs MOND ranking (cluster check 6).
- **The two halves of the paper convert a₀ = 1.2 × 10⁻¹⁰ m s⁻² to (km/s)²/kpc slightly differently.** Table 1 was computed with a₀ = 3702.789, from the rounded unit constant 1 (km/s)²/kpc = 3.2408 × 10⁻¹⁴ m s⁻²; Table 2 with a₀ = 3702.813, from the exact parsec conversion. The two differ by 6.4 × 10⁻⁶ relative. Each half of this package uses the constant its own published table used, so both tables reproduce exactly. The difference moves the galaxy MOND median χ²/N by 1 × 10⁻⁵ (1.142686 against 1.142696) and changes nothing the paper quotes to three decimals; `reproduce_galaxies.py` prints the median under both conversions so you can see the size of it.

## Data and licences

The cluster inputs are unmodified copies of files from two Zenodo records, both released under CC BY 4.0: Famaey, Pizzuti & Saltas (https://doi.org/10.5281/zenodo.15299349) and Mistele et al. (https://doi.org/10.5281/zenodo.15476959).

The galaxy rotation curves are from SPARC: Lelli, McGaugh & Schombert (2016), AJ 152, 157, at http://astroweb.case.edu/SPARC/. They are **not** redistributed in this package. `fetch_sparc.py` downloads the public archive `Rotmod_LTG.zip`, keeps only the 133 files the paper scores, and checks every one of them against `data/sparc_manifest.csv`, so you can be sure you are running on the same bytes the published result used. SPARC's terms ask that users of the data cite that paper.

`data/SOURCES.md` has the full citations and a checksum for every file.

The code in this package (`reproduce_clusters.py`, `reproduce_galaxies.py`, `fetch_sparc.py` and `hermes_clusters/`) is released under the MIT licence; see `LICENSE`.

## Tested with

Python 3.10.14, numpy 1.26.4, scipy 1.15.3 and pandas 2.3.1, on Windows 10.
