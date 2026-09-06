# The chi-squared error floor (σ_int² = 386)

Hermes scores rotation curves with

    χ²_ν  =  (1/N) Σ (Vobs − Vmodel)² / (errV² + σ_int²),     σ_int² = 386 (km/s)²

`σ_int = √386 = 19.65 km/s`. This note records where that number comes from and
two properties of it that a reader needs before interpreting any χ² in this
repository. Neither property is an error, and both are visible in Paper 1
Appendix A — but until now neither was written down here, and an external
reviewer reasonably annotated the constant as "according to the author."

## Where 386 comes from

It is **SPARC's own scatter budget in velocity space**, not a value fitted to
Hermes. The SPARC team report a total scatter of ~0.12 dex on the radial
acceleration relation. At fixed radius `g ∝ V²/R`, so a scatter of 0.12 dex in
`g` is 0.06 dex in `V`:

    dV/V  =  10^0.06 − 1  =  0.148

Setting `19.65 = 0.148 · V` implies a characteristic velocity of **132.6 km/s**.
SPARC's median observed velocity across the 133-galaxy sample is **134.0 km/s** —
agreement to about 1%.

| SPARC scatter | implied characteristic V |
|---|---|
| 0.11 dex | 145.5 km/s |
| **0.12 dex** | **132.6 km/s** |
| 0.13 dex | 121.7 km/s |

The reasoning behind the choice: the published `errV` column captures measurement
noise only, not the systematics that dominate rotation-curve work — distance,
inclination, mass-to-light ratio, asymmetric drift. Some floor is needed. Rather
than invent a Hermes-specific one, the floor adopts the uncertainty the data
source already reports for itself. It is not a field standard; other groups
absorb the same systematics into nuisance parameters instead.

## Property 1: the floor dominates the error budget

Measured over all 133 galaxies / 3073 points:

| | |
|---|---|
| σ_int | 19.65 km/s |
| median `errV` | 4.57 km/s (IQR 2.89 – 7.70) |
| points where `386 > errV²` | 2975 / 3073 (**96.8%**) |
| median share of variance from the floor | **94.9%** |
| points where the floor supplies >90% of variance | 68.6% |

For most of the sample the score is therefore close to a rescaled sum of squared
residuals rather than a conventional χ². **Absolute χ² values in this repository
are not directly comparable to published fits that use `errV` alone.** Compare
either at matched floor values (Appendix A tabulates 0 / 100 / 386 / 900) or with
a floor-free statistic such as the median absolute velocity residual.

## Property 2: the floor is not scale-free

`σ_int` is one global constant applied to galaxies spanning a factor of ~20 in
rotation velocity, so it is a much larger correction for a dwarf than for a
giant. This does not merely compress the scores — it **reorders** them:

    Spearman( χ² at floor 0 , χ² at floor 386 )  =  0.71

A pure rescaling would give ~1.0.

| sample | floor as % of peak Vobs | median χ² at 386 | median χ² at 0 |
|---|---|---|---|
| dwarfs, Vmax < 80 km/s (n=31) | 34.0% | 0.38 | 29.4 |
| mid, 80–150 km/s (n=50) | 18.4% | 0.77 | 21.8 |
| giants, Vmax > 150 km/s (n=52) | 8.7% | 3.78 | 93.4 |

Consequently `Spearman(χ² at 386, Vmax) = 0.705`: at this floor the score
correlates strongly with galaxy mass. Any fixed χ² cutoff used to classify fit
quality will therefore encode galaxy size as well as model performance — 97% of
dwarfs but only 23% of giants fall below χ² = 2, and 13 of the 14 galaxies above
χ² = 8 are massive spirals. Classification schemes built on these scores should
be mass-normalised, or built on a floor-free residual, or stated as
size-conditional.

## What the floor does not affect

The same floor is applied to Hermes and to the MOND comparator, so the
model-vs-model comparison carries no advantage from the choice. Paper 1
Appendix A Table A1 reports both models at four floor values, and the ordering of
the two models on each statistic is stable across all of them:

| floor | Hermes median | MOND median | Hermes trim-10% | MOND trim-10% | Hermes win rate |
|---|---|---|---|---|---|
| 0 | 37.79 | 45.44 | 76.55 | 83.73 | 0.519 |
| 100 | 4.07 | 3.79 | 6.57 | 8.21 | 0.511 |
| 386 | 1.32 | 1.14 | 2.11 | 2.71 | 0.511 |
| 900 | 0.61 | 0.53 | 0.96 | 1.26 | 0.511 |

Note that the two statistics disagree, at every floor above zero: MOND has the
better **median**, Hermes the better **trimmed mean** and the better win rate.
That is a statement about the shape of the two χ² distributions — Hermes has
fewer extreme failures, MOND a better typical galaxy — and it is not caused by
the floor.

The strongest floor-independent comparison is the median absolute velocity
residual: **Hermes 20.67 km/s, MOND 25.82 km/s**.

## Reproducing the numbers in this note

Every figure above comes from the repository's own code — `paper1/verify_hermes.py`
provides `read_rotmod`, `hermes_model`, `mond_model` and `chi2_nu`, and `chi2_nu`
takes the floor as its `sigma_int_sq` argument. The Appendix A table is checked in
CI by `paper1/verify_appendix_a.py`.

Reported medians here are from the current repository gate (log-grid shear), so
the floor-386 Hermes median reads 1.323 rather than Paper 1's published 1.312.
That difference is the derivative-convention lineage documented in
[`gate_version_history.md`](gate_version_history.md), unrelated to the floor.
