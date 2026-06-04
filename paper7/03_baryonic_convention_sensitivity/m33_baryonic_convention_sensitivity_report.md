# M33 Baryonic Convention Sensitivity Diagnostic

**Locked result status:** the archived Hermes and MOND results are unchanged. This is a diagnostic branch only.

**Board:** M33-Corbelli source-native thick-disk reconstruction.

## Effective baseline M/L estimate

Using the Corbelli stellar mass baseline $M_* = 4.9 \times 10^9 M_\odot$, total M33 $[3.6] = 3.2$, distance modulus 24.60, and solar 3.4/3.6 micron magnitude 3.24 gives:

- $L_{3.6} \approx 7.18 \times 10^9 L_\odot$
- $M_*/L_{3.6} \approx 0.683$

Thus a SPARC-like unit $M/L=1$ convention would be roughly a velocity scaling of $f=\sqrt{1/0.683} \approx 1.210$, not 1.5x.


### Provenance clarification: this is not a literal SPARC-equivalent unit-M/L board

The $f\approx1.21$ branch is a **SPARC-like unit-M/L diagnostic analogue**, not a true SPARC-equivalent rotmod convention lane.

Reason: Corbelli's public Table 1 supplies a stellar **mass surface-density profile** derived from pixel-SED / stellar-population modeling. It does **not** supply the raw 3.6 μm / W1 light profile needed to rebuild a literal SPARC-style $\Upsilon_\star=1$ stellar component from luminosity. Therefore the $M/L=1$ comparison in this diagnostic is based on an estimated effective near-IR $M/L\approx0.683$ and a corresponding velocity rescaling $f=\sqrt{1/0.683}\approx1.21$.

Archive label for this branch:

> **SPARC-like unit-M/L diagnostic analogue / stellar-mass convention sensitivity branch.**

Do **not** describe it as:

> SPARC-equivalent rotmod reconstruction, literal $\Upsilon_\star=1$ board, or replacement for the locked source-native M33 board.

The locked board remains **M33-Corbelli source-native thick-disk reconstruction**: source-native SPS/pixel-SED stellar mass profile, radially varying M/L provenance, helium factor 1.33, H$_2$ prescription, flaring stellar disk, gas half-thickness 0.5 kpc, and $V_{\rm bul}=0$.

## Primary age scaling series

|   scale_Vdisk |   implied_ML_factor |   implied_global_ML_W1 |      g98 |   beta |   Hermes_chi2nu |   Hermes_outer_chi2nu |   Hermes_outer_mean_residual |   MOND_chi2nu |   MOND_outer_chi2nu |   MOND_outer_mean_residual |
|--------------:|--------------------:|-----------------------:|---------:|-------:|----------------:|----------------------:|-----------------------------:|--------------:|--------------------:|---------------------------:|
|           1   |                1    |                  0.683 |  2653.95 |  1.834 |           3.126 |                 5.794 |                      -51.622 |         0.119 |               0.228 |                     -9.394 |
|           1.5 |                2.25 |                  1.536 |  5887.88 |  1.074 |           2.469 |                 3.531 |                      -39.927 |         0.994 |               0.322 |                     10.365 |
|           2   |                4    |                  2.731 | 10415.4  |  0.424 |           3.009 |                 2.373 |                      -32.38  |         4.269 |               2.139 |                     29.947 |
|           2.5 |                6.25 |                  4.267 | 16236.5  | -0.01  |           3.612 |                 1.428 |                      -23.589 |         9.916 |               5.568 |                     49.014 |
|           3   |                9    |                  6.144 | 23351.1  | -0.243 |           4.671 |                 0.724 |                      -11.86  |        17.926 |              10.493 |                     67.598 |

## Older bracket scaling series

|   scale_Vdisk |   implied_ML_factor |   implied_global_ML_W1 |      g98 |   beta |   Hermes_chi2nu |   Hermes_outer_chi2nu |   Hermes_outer_mean_residual |   MOND_chi2nu |   MOND_outer_chi2nu |   MOND_outer_mean_residual |
|--------------:|--------------------:|-----------------------:|---------:|-------:|----------------:|----------------------:|-----------------------------:|--------------:|--------------------:|---------------------------:|
|           1   |                1    |                  0.683 |  2653.95 |  1.699 |           3.296 |                 6.109 |                      -53.049 |         0.119 |               0.228 |                     -9.394 |
|           1.5 |                2.25 |                  1.536 |  5887.88 |  0.883 |           2.48  |                 4.042 |                      -42.903 |         0.994 |               0.322 |                     10.365 |
|           2   |                4    |                  2.731 | 10415.4  |  0.245 |           2.687 |                 2.891 |                      -35.999 |         4.269 |               2.139 |                     29.947 |
|           2.5 |                6.25 |                  4.267 | 16236.5  | -0.133 |           3.02  |                 1.737 |                      -26.553 |         9.916 |               5.568 |                     49.014 |
|           3   |                9    |                  6.144 | 23351.1  | -0.309 |           4.014 |                 0.813 |                      -13.744 |        17.926 |              10.493 |                     67.598 |

## Continuous primary-age crossover/optimum

|   best_f |   best_chi2nu |   best_outer_chi2nu |   best_beta |   best_g98 |   best_ML_factor |   best_global_ML_W1 |   chi2nu_cross_f | notes                                                             |
|---------:|--------------:|--------------------:|------------:|-----------:|-----------------:|--------------------:|-----------------:|:------------------------------------------------------------------|
|    1.331 |        2.3813 |              4.0985 |      1.3286 |    4650.09 |           1.7716 |              1.2094 |              nan | grid/fine global minimum over f=0.5-5; no chi2nu≈1 crossing found |

## Outer-disk inversion summary

|   outer_N |   median_required_Vbar |   median_measured_Vbar |   median_required_to_measured_Vbar_ratio |   median_required_Vdisk_scale |   median_required_stellar_ML_factor |   p90_required_stellar_ML_factor |   max_required_stellar_ML_factor |   baseline_Mstar |   median_equivalent_Mstar_if_uniform |   median_extra_Mstar_if_uniform |   scale_to_best_hermes |   ML_factor_to_best_hermes |   implied_Mstar_to_best_hermes |   baseline_global_ML_W1 |   best_global_ML_W1 |
|----------:|-----------------------:|-----------------------:|-----------------------------------------:|------------------------------:|------------------------------------:|---------------------------------:|---------------------------------:|-----------------:|-------------------------------------:|--------------------------------:|-----------------------:|---------------------------:|-------------------------------:|------------------------:|--------------------:|
|        27 |                  81.95 |                 46.565 |                                    1.839 |                         2.189 |                               4.793 |                            6.518 |                            7.246 |          4.9e+09 |                          2.34856e+10 |                     1.85856e+10 |                  1.331 |                      1.772 |                    8.68065e+09 |                   0.683 |               1.209 |

## Interpretation

Increasing the stellar disk helps Hermes at first, but only partially. The primary-age global minimum occurs at $f \approx 1.331$, equivalent to a stellar mass/M/L factor of 1.772 over the source-native baseline. That corresponds to $M_* \approx 8.68 \times 10^9 M_\odot$ and $M/L_{3.6} \approx 1.21$.

Hermes does not reach $\chi^2_\nu \approx 1$ anywhere in the tested continuous range. A heavier disk reduces the outer underprediction, but increasing $g_{98}$ also lowers $\beta$, and the global fit develops inner/mid-disk overshoot. MOND is very sensitive to over-heavy disks in the normal way: it fits the baseline board cleanly and becomes worse as the stellar disk is scaled upward.
