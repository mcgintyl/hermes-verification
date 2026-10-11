# Supplementary Materials

**A Wear-Activated Density-Release Refinement of the Hermes Gravity Equation: A Paper 1 Addendum**

**Louis Albert McGinty**

mcgintyl@grinnell.edu

Independent Researcher, Team Hermes

AI-assisted research; full disclosure in main text

---

# Supplement A: Candidate Rejection Path and Refinement Constraints


## 1. Motivation

See the main text, Section 1.


## 2. Refinement Constraints

Any modification to an equation with no per-galaxy parameters risks converting it into a fitting exercise. The Hermes equation's credibility rests on its rigidity: once the constants were locked, they were never adjusted per galaxy. Paper 8 inherits that discipline. The refinement is only meaningful if it survives under the same constraints that made Paper 1 credible. These constraints are not procedural decoration. They are the evidentiary standard.

The surviving operator meets all of the following requirements; several were adopted during the search (Section 3).

First, no per-galaxy parameters. The refinement may not introduce any quantity that is adjusted, optimized, or selected on a galaxy-by-galaxy basis. Every galaxy receives the same equation with the same constants applied identically.

Second, no refitting of Paper 1 constants. The baseline constants, the coherence ceiling pi, the light-speed normalization in psi, the floor 1/sqrt(2pi), and the gate threshold a_k = 1585, are inherited and frozen. The refinement may introduce new operators, but it may not modify existing constants.

Third, use only frozen Hermes constants where possible. Any new constant introduced by the refinement should be derived from the existing constant set rather than independently optimized against the sample. This limits the effective degrees of freedom to discrete choices of form, not continuous fitted scalars.

Fourth, apply the refinement to the positive retained-coherence term, not to the net age score. The age score beta has a physical floor at -1/sqrt(2pi) representing fully eroded coherence. Any density-release mechanism must act on the surviving positive coherence (pi times exp(-psi)) rather than on the net score, to avoid amplifying the floor into unphysical regimes.

Fifth, preserve or improve SPARC aggregate behavior without introducing casualties. In this study, a material casualty is defined as any galaxy degraded by more than 0.5 chi-squared-nu. The refinement must also avoid creating new catastrophic cases (chi-squared-nu > 5 or chi-squared-nu > 10). Improvement in the median is desirable but secondary to safety.

Sixth, avoid a direct Baryonic Tully-Fisher Relation collision. The density-release numerator must be a universal constant, not the galaxy-specific peak acceleration g_98. This constraint covers the numerator only: the size of the release still depends on g_98 through psi, in pi*exp(-psi) and in W(psi), as the Paper 1 age score does. Any density-release term that anchors the outer-disk behavior to the peak baryonic acceleration g_98 risks injecting surface-brightness scatter into the BTFR, which is empirically one of the tightest scaling relations in extragalactic astronomy.

Seventh, fit no constant to M33. M33 motivated this search and every candidate was scored on it, so it is a motivating case here, not an independent test; its result was cited against candidate 4 only after that candidate had failed on SPARC.

Eighth, report all failed candidates. Each failure identifies a structural requirement that the final form must satisfy.


## 3. Candidate Rejection Path

Five candidates were tested against the requirements of Section 2 before the surviving operator. Each failed at least one, and each failure shaped the operator that survived; other forms that also pass are noted in the main text, Section 4.

The first candidate replaced the raw peak baryonic acceleration g_98 in the wear exponent psi with a logarithmically compressed version, intended to smooth a population-level kink in psi across g_98 quartiles (not shown). The compression narrowed the spread of psi across the sample, but psi still sorted by density quartile, since the compression is monotone in g_98. It also substantially degraded SPARC performance, producing 53 to 58 galaxies with chi-squared-nu above 5.0 depending on the compression variant. The high-wear branch of the equation, which governs the oldest and densest galaxies near the coherence floor, is load-bearing. Smoothing it away destroyed fits that were already working.

The second candidate left psi intact and instead introduced a radial density-release factor eta(R) that allowed retained coherence to express more strongly in low-density outer regions. A sign guard set eta = 1 whenever the net age score beta was negative, preventing the release from amplifying the coherence floor. This candidate produced the strongest M33 result of any version tested: the outer-disk chi-squared-nu dropped from 5.79 to 0.72 (superseded Paper 7 v2 gate; not rerun). However, it created 8 material SPARC casualties and 3 new catastrophic cases. The failures traced to two structural flaws. The density-release anchored its outer asymptote to the galaxy-specific peak acceleration g_98, creating a BTFR-collision risk (Section 2, sixth). And the hard sign guard at beta = 0 created a topological discontinuity: a galaxy at beta = +0.001 received the full density release while a galaxy at beta = -0.001 received none, with no physical basis for the boundary. The lesson: the density seam is real and exploitable, but coupling it to galaxy-specific core density is structurally wrong, and a binary sign guard is not a physical transition mechanism.

The third candidate addressed both flaws. It decoupled the outer asymptote from g_98 by replacing it with a universal numerator a_u = a_k/pi, derived from the existing gate threshold. It applied the density release to the positive decay term pi times exp(-psi) rather than to the net age score, eliminating the sign guard entirely: as psi grows large, the positive term decays toward zero and the release has progressively less to amplify, so the floor is reached continuously without a discontinuity. This version reduced SPARC casualties from 8 to 1 and created no new catastrophic cases. The single remaining casualty, UGC 05750, is a very low-wear, high-retention system near the coherence ceiling. The refinement added outer support to a galaxy that did not need it because it had barely begun the decay process. The lesson: the BTFR decoupling and the positive-term application were kept, and a mechanism to suppress activation in low-wear systems was added next (W, Section 4).

The fourth candidate proposed that the physical core and the low-density branch represent orthogonal domains whose coherence contributions should combine via Pythagorean sum rather than linear addition. The resulting expression, P(R) = exp(-psi) times sqrt(pi^2 + W(psi)^2 times pi times a_k / Gamma(g_bar)), is mathematically elegant: it eliminates the logistic transition function, the max floor, and the nested eta structure in favor of a single geometric expression. It produced 15 material SPARC casualties and 4 new catastrophic cases. The Pythagorean sum did not shut off fast enough in the dense core, so it leaked a parasitic boost there. In addition, the pi cap on total P(R) strangled the low-density asymptote: in M33's outer disk, P(R) saturated at pi regardless of the component age, erasing the entire component-age diagnostic axis. The lesson: this smooth blend across the full density range failed; the surviving operator keeps the two regimes distinct and selects between them rather than averaging them.

The fifth candidate reframed the problem as a thermodynamic phase transition between a bound physical state and an unbound low-density state. The probability of being in the unbound state was modeled as the product of wear activation W(psi) and a Fermi-Dirac-style spatial transition S(g). The low-density multiplier used Pythagorean geometry only within the unbound state, avoiding the continuous-superposition error. On M33 it gave an outer chi-squared-nu of about 1.10 at the 2 Gyr component-age bracket (logistic width 1/pi; the surviving form gives 1.53 there), but it created 14 material SPARC casualties and 2 new catastrophic cases. The reported values use the theoretically motivated width 1/pi; an isolation test with the reused Paper 1 width 0.30 gave the same casualty count, so the structural form itself was too permissive. The casualties clustered in the acceleration range between g_bar of about 594 and a_k = 1585 (km/s)^2/kpc, a transition zone where the phase switch had activated but the low-density branch was not yet dominant. The lesson: the phase-transition framing captures part of the density-seam behavior, but smooth amplitude scaling in the transition zone adds support where the data do not require it. The surviving operator withholds the low-density contribution until the low-density branch is actually stronger than the local three-dimensional field.

The surviving refinement is therefore not the most elegant candidate tested, nor the best on median: candidate 3 (median chi-squared-nu 1.115) and candidate 5 (1.053) both score lower than the surviving form (1.189). Of the five candidates and the surviving form, it is the only one that meets all eight constraints, so the choice rests on casualties, not on the median. As archived, it is candidate 3 with the wear gate W added. W suppresses the release in low-wear systems, which, because psi is proportional to t50 g_98 and g_98 varies far more than t50 across the sample, are mostly low-g_98 galaxies. The universal numerator removes g_98 from the density contrast, although g_98 still sets the release amplitude through psi. The positive-term application avoids the sign-guard discontinuity. The max-selected density contrast, already present in candidate 3, withholds the release until Gamma(g_bar) < a_u; candidates 4 and 5 replaced it with smooth alternatives and failed. These features answer specific empirical failures, but they are not the only forms that pass (main text, Section 4). The logarithmic compression, the logistic S and the pi cap were also present in candidate 3; no failure forced them, and removing any one of them still leaves zero casualties.


## 4. The Locked Refinement Form

The baseline Hermes equation (McGinty 2026a) predicts the gravitational acceleration at each measured radius R as:

g_model(R) = g_bar(R) * [1 + beta * phi(R)]

where g_bar is the baryonic acceleration from visible mass, beta is the age score, and phi(R) is the spatial gate. The age score is computed from the wear variable psi as:

beta = pi * exp(-psi) - 1/sqrt(2*pi)

where psi = 2*pi * t_50 * g_98 / c, with t_50 the galaxy's stellar age, g_98 the 98th-percentile peak baryonic acceleration, and c expressed in the same unit convention used in Paper 1. The gate phi(R) is constructed from the baryonic acceleration profile through six fixed processing steps described in Paper 1. All of these components are inherited without modification.

The refinement introduces a density-release operator eta_WA(R) that modifies how retained coherence couples to low-density outer regions. It is applied to the positive decay term, not to the net age score:

g_model(R) = g_bar(R) * [1 + phi(R) * (pi * exp(-psi(R)) * eta_WA(R) - 1/sqrt(2*pi))]

When eta_WA = 1 everywhere, this reduces exactly to the baseline equation. The refinement adds structure only where eta_WA deviates from unity.

The operator eta_WA(R) is built from the following fixed components.

**Logarithmic compression.** The baryonic acceleration at each radius is compressed logarithmically:

Gamma(g) = a_k * ln(1 + max(g, 0) / a_k)

where a_k = 1585 (km/s)^2/kpc is the gate threshold inherited from Paper 1. For high accelerations Gamma grows logarithmically, but this does not change the operator: above g_bar ≈ 594 the max(1, ·) floor already sets eta_U = 1. Where the release is active, Gamma is up to 15% below g, so the compression slightly strengthens the release and moves its onset from g_bar ≈ 504.5 to about 594. For low accelerations, Gamma approximately equals g, preserving the linear scaling in the outer disk.

**Wear activation.** The refinement is suppressed in galaxies with little wear (small psi):

W(psi) = 1 - exp(-2*pi^2 * psi)

For low-wear systems (psi close to zero), W is small and the release is suppressed; across the 133 galaxies W runs from 0.22 to 1, and because psi is proportional to t50 g_98, low wear mostly means low g_98 rather than young age. For systems with significant wear, W saturates rapidly toward unity. The coefficient 2*pi^2 is a closed-form expression of the inherited constant pi; interpretation is deferred to Supplement E. This component was introduced after the third candidate (Section 3) revealed that the refinement over-activated in low-wear systems near the coherence ceiling.

**Spatial transition.** A logistic factor scales the density release along the radial profile (the release itself switches on only below g_bar ≈ 594 (km/s)^2/kpc, through the max function below, where S is already at least 0.889):

S(g) = [1 + exp((g - a_k) / (0.30 * a_k))]^(-1)

S falls toward 0 as g rises past a_k. The center a_k is the Paper 1 knee threshold; the width ratio 0.30 reuses the width of the Paper 1 gate logistic, which there acts on x = R/r_knee (centered at 1.40), not on g. No fitted parameter is introduced.

**Universal numerator.** The density contrast that drives the release is anchored to a universal constant rather than to the galaxy-specific peak acceleration:

a_u = a_k / pi (approximately 504.5 (km/s)^2/kpc)

This decoupling was required by the second candidate's failure (Section 3): anchoring the outer asymptote to g_98 created a BTFR-collision risk. The value a_k/pi is derived from existing constants. When Gamma(g_bar) falls below a_u, the contrast term enters the square-root branch sqrt(a_u/Gamma). The interpretation of this branch is discussed in Supplement E.

**The local density-release operator.** At each radius, the operator evaluates whether the low-density branch exceeds the local three-dimensional field:

eta_U(R) = min[pi, 1 + S(g_bar(R)) * (max(1, sqrt(a_u / Gamma(g_bar(R)))) - 1)]

The max function selects the larger of two values: unity, corresponding to no density-release beyond baseline, or sqrt(a_u/Gamma), corresponding to the low-density branch. When Gamma exceeds a_u, the max evaluates to unity and the release contributes nothing beyond baseline. When Gamma falls below a_u, the low-density branch dominates and the release activates. The logistic S(g) modulates the transition. The min function caps the total operator at pi, preventing the local multiplier from growing without bound in the deep outer disk. Unlike the rejected Pythagorean candidate (Section 3), the cap is applied to the multiplier, not to the total positive coherence term.

**The wear-gated operator.** The local operator is modulated by the global wear activation:

eta_WA(R) = 1 + W(psi_sys) * [eta_U(R) - 1]

When W is zero (psi = 0), eta_WA equals 1 regardless of local conditions. When W is unity (well-worn systems), eta_WA equals eta_U. The wear gate uses the systemic wear variable psi_sys computed from the galaxy's global stellar age, not the local psi(R). This keeps the release off in galaxies with little wear.

**Component-age extension.** For galaxies with evidence of distinct dynamical components, the positive decay term can use a spatially varying effective age:

psi(R) = 2*pi * t_eff(R) * g_98 / c

For the standard SPARC execution, psi(R) = psi_sys at all radii. For the M33 component-age diagnostic, the positive decay term uses psi(R) derived from t_eff(R), with t_eff equal to the global stellar age for the inner disk (R <= 10 kpc) and a younger effective age for the outer disk (R > 10 kpc). The wear gate W remains defined by the systemic bracket. The component-age extension does not modify any constant or operator; it modifies only the age input for regions where independent observations, read after the M33 residual was seen (Paper 7), motivate a multi-component interpretation.

**Summary of constants.**

| Constant | Value | Role | Status |
|----------|-------|------|--------|
| pi | 3.14159... | Coherence ceiling, operator cap | Inherited |
| e | 2.71828... | Exponential decay | Inherited |
| c | Paper 1 units | Light-speed normalization in psi | Inherited |
| 1/sqrt(2pi) | 0.3989... | Coherence floor | Inherited |
| a_k | 1585 (km/s)^2/kpc | Gate threshold, transition center | Inherited |
| 0.30 | -- | Logistic width ratio (Paper 1 gate-logistic width, reused in g) | Reused |
| a_u = a_k/pi | 504.5 (km/s)^2/kpc | Universal density-release numerator | Derived |
| 2pi^2 | 19.74 | Wear activation rate | Derived |

The inherited constants are taken unchanged from Paper 1, and no continuous parameter was fitted to SPARC or M33. The two derived constants are closed-form expressions of inherited constants, chosen by trial on the 133-galaxy board; other closed forms also pass (main text, Section 4). No per-galaxy parameters are introduced on SPARC; the M33 diagnostic's outer age is a free input (Supplement D).

---

# Supplement B: Signed Residual Audit

The wear-activated eta refinement was audited for directional bias in the outer disk. Zero-casualty chi-squared results can mask systematic sign flips if a correction shifts residuals from under-prediction to over-prediction without changing the variance.

For the outermost 5 measured points per galaxy, the median (across galaxies) of the mean signed residual V_model - V_obs shifted from -20.51 km/s (baseline) to -10.80 km/s (corrected); negative means under-predicted. The median galaxy remains under-predicted after correction. Of 133 galaxies, 16 experienced sign flips in the outermost 5 points, all from under-prediction to over-prediction. Five of the 16 flipped galaxies now overshoot by more than +5 km/s (DDO 161, UGC 06983, UGC 12732, NGC 2366, UGC 05005). Across all 133, 29 galaxies have an outer-5 mean residual above +5 km/s after correction, up from 22 at baseline. No galaxy flipped from over-prediction to under-prediction; this holds by construction, since eta_WA >= 1 can only raise the model velocity.

The correction reduces outer-disk starvation, and the median galaxy stays under-predicted, but the shift is one-sided by construction, so the flipped galaxies are kept as a watchlist (paper8_sparc_outer5_sign_flip_watchlist.csv). Per-galaxy signed residuals are in paper8_sparc_outer_signed_residual_audit_last3_last5.csv.

---

# Supplement C: Improved-Galaxy Characterization

Thirteen galaxies improved by more than 0.5 chi-squared-nu under the wear-activated eta. Most share an extended, gas-rich, low-density outer disk, which is where the operator is built to act.

Compared to the remaining 120 galaxies, the improved class has: more negative baseline outer residuals (median -41.83 vs -16.76 km/s), greater radial extent beyond the knee (median Rlast/rknee 6.95 vs 2.44, with 12/13 exceeding 2.0), higher outer gas fraction (0.518 vs 0.204), lower disk surface brightness (median over each galaxy's radii where SPARC's 3.6 micron SBdisk is nonzero: 1.115 vs 14.92 L_sun/pc^2), lower outer-to-peak ratio rho = g_bar(R_last)/g_98 (0.036 vs 0.148), and so more Chimeras (Paper 1's term for rho < 0.02; 6/13 vs 21/120).

Two galaxies (NGC 6674 and NGC 0289) improve numerically but remain poor fits and should not be treated as clean rescues.

Full characterization data is provided in paper8_sparc_13_improved_galaxies_characterization.csv.

---

# Supplement D: M33 Component-Age Sweep and Literature Notes

## Sweep Results

The wear-activated eta was combined with a component-resolved age model on M33, using the inner disk (R <= 10 kpc) at Paper 7's global stellar age and the outer disk (R > 10 kpc) at younger effective ages. Neither eta alone nor component age alone closes the outer-disk residual. Together they cut the outer chi2_nu by about three quarters, but every outer point stays under-predicted, and the outer chi2_nu keeps falling as the outer age is lowered, so the sweep does not select an age. For reference, MOND gives outer (R > 10 kpc) chi2_nu = 0.23 (full board 0.12) on the same board with no age input (Paper 7). All values use the canonical Paper 1 gate, as corrected for M33 in Paper 7 v3.

Primary bracket (inner t50 = 6.0 Gyr):

|     Outer age (Gyr)     | Full chi2_nu | Inner (R <= 10) chi2_nu | Outer (R > 10) chi2_nu | Far outer (R >= 15) chi2_nu |
| :---------------------: | :----------: | :---------------------: | :--------------------: | :-------------------------: |
| 6.0 (uniform, eta only) |    1.682     |          0.885          |         2.597          |            3.150            |
|           3.0           |    1.303     |          0.885          |         1.782          |            2.290            |
|           2.0           |    1.186     |          0.885          |         1.532          |            2.018            |
|           1.5           |    1.131     |          0.885          |         1.413          |            1.886            |
|           1.0           |    1.077     |          0.885          |         1.297          |            1.756            |

Conservative bracket (inner t50 = 7.1 Gyr):

|     Outer age (Gyr)     | Full chi2_nu | Inner (R <= 10) chi2_nu | Outer (R > 10) chi2_nu | Far outer (R >= 15) chi2_nu |
| :---------------------: | :----------: | :---------------------: | :--------------------: | :-------------------------: |
| 7.1 (uniform, eta only) |    1.844     |          0.914          |         2.912          |            3.474            |
|           3.0           |    1.317     |          0.914          |         1.780          |            2.288            |
|           2.0           |    1.201     |          0.914          |         1.531          |            2.016            |
|           1.5           |    1.145     |          0.914          |         1.411          |            1.884            |
|           1.0           |    1.091     |          0.914          |         1.295          |            1.754            |

Full results in paper8_m33_wear_eta_component_age_summary.csv.

Baseline Hermes without eta: component age alone (even at newborn outer limit, psi = 0 for R > 10 kpc) reaches only chi2_nu = 4.12 in the outer disk.

## Literature Notes

Published work on M33's outer gas allows, but does not establish, an outer-component age of 1 to 2 Gyr, read as the effective age of the outer gas layer rather than the stellar half-mass age of the whole outer disk. As in Paper 7, this literature was read after the residual pattern was seen.

Putman et al. (2009) found M33's HI mass extends well beyond the star-forming disk with extended warps, an arc from the northern warp to the disk, and a southern cloud with a filament back to the galaxy. Their orbit analysis favored tidal disruption by M31 at 1 to 3 Gyr ago.

Semczuk et al. (2018) simulated a 37 kpc M31-M33 passage, a scenario proposed to explain star formation peaks in both galaxies about 2 Gyr ago; their best match to M33 comes 0.75 Gyr after pericenter. Tides alone reproduced only the inner gaseous warp; adding ram pressure from M31's hot gas halo gave an S-shaped warp qualitatively similar to the observed one.

Corbelli and Burkert (2024), using orbits compatible with Gaia-EDR3 proper motions, find M33 is probably on first infall (probability above 70%); any earlier M31 passage was at about 100 kpc or more, over 3 Gyr ago, and very unlikely to have made the warp. They suggest that more recent gas accretion, possibly from a cosmic filament, misaligned the outer disk.

Barker et al. (2011) found an age-gradient reversal in M33's outer disk: the field at 9.1 kpc has a mean age of approximately 3 +/- 1 Gyr, while the farther field at 11.6 kpc is older at approximately 7 +/- 2 Gyr. Only the farther field lies in this paper's outer region (R > 10 kpc), beyond the disk break; its stars are old, so the younger component age refers to the outer gas layer, not to these stars.

The literature does not give a single clean answer. The 1 to 3 Gyr timing comes from older work favoring a recent M31 encounter. Newer proper-motion work favors first infall, and suggests the outer-disk misalignment comes from gas accretion, possibly from a cosmic filament, which sets no single event age. The 1 to 2 Gyr bracket is therefore an input motivated by this literature, not an age the sweep measures.

---

# Supplement E: Dimensional-Transition Interpretation (Provisional)

This supplement describes a provisional interpretive reading of the wear-activated eta operator. These notes are not required to evaluate the empirical results in the main text. They are provided for readers interested in the structural implications of the operator's form.

In the candidates tested (Supplement A), smooth superposition of the core and low-density regimes failed, while sharp branch-selection survived. The surviving operator contains a max function that selects between two scaling behaviors:

When Gamma(g_bar) exceeds a_u: the multiplier is unity and gravity follows the baseline three-dimensional field.

When Gamma(g_bar) falls below a_u: the multiplier enters the square-root branch, and for g_bar well below a_u the refinement adds about A\*sqrt(a_u\*g_bar) to the acceleration, with A = phi\*S\*W\*pi\*exp(-psi), the form of MOND's deep-regime term. Where g_bar falls as 1/R^2 this added part falls as 1/R; the rest of the field still scales with g_bar, and below g_bar of about 49 (km/s)^2/kpc (17 of the 3,073 SPARC points, in 8 galaxies; the two points with g_bar <= 0 get no release) the pi cap fixes the multiplier, so the scaling returns to 1/R^2.

A 1/R force law is the scaling expected from flux spreading across a two-dimensional surface rather than a three-dimensional volume. The SPARC and M33 audits are consistent with a dimensional-transition interpretation in which the equation selects between three-dimensional and two-dimensional propagation geometries depending on local density.

The region between a_u (corresponding to raw g_bar approximately 594) and a_k (1585) represents a transition zone where the spatial switch S(g) has activated but the low-density branch is not yet dominant. Candidates that boosted this zone created SPARC casualties. The surviving operator withholds support until the low-density branch is actually stronger than the local field.

The wear activation coefficient 2*pi^2 equals the surface area of a unit 3-sphere. This identity is noted as a structural observation; whether it reflects deeper physical significance or coincidence remains an open question.

This interpretation is structurally supported by the candidate rejection pattern but is not independently established. It should be treated as provisional. We do not claim direct observation of a two-dimensional physical substrate. Rather, we note a correspondence: the max selection floor and the low-density square-root branch of the surviving operator have the form of a dimensional collapse. Other forms also pass the SPARC constraints (main text, Section 4), so the data do not single out this signature.

---

# Supplement F: Production Data Package and CSV Manifest

The CSV files below are deposited with this record on Zenodo. paper8/reproduce.py in github.com/mcgintyl/hermes-verification (commit d3355f4) rebuilds them from public data, together with every SPARC and M33 result of this version except the candidate 1 and 2 entries; paper8/verify_eta.py implements the operator.

| Filename | Contents | Rows |
|----------|----------|------|
| paper8_sparc_133_model_comparison.csv | Per-galaxy results: baseline, MOND, wear-activated eta chi2_nu, deltas, flags | 133 |
| paper8_sparc_per_radius_baseline_mond_wear_eta.csv | Per-radius predictions for all three models | 3,073 |
| paper8_sparc_aggregate_metrics.csv | Summary statistics for all models | 3 |
| paper8_iteration_comparison_metrics.csv | Aggregate metrics for baseline Hermes, MOND, candidates 3 to 5 (candidate 5 at logistic width 1/pi) and the final operator; candidates 1 and 2 not included | 6 |
| paper8_m33_wear_eta_component_age_summary.csv | M33 chi2_nu by inner and outer age, four masks (canonical gate) | 10 |
| paper8_m33_wear_eta_component_age_summary_long.csv | M33 results by mask and outer age, long format, with residual statistics | 40 |
| paper8_m33_wear_eta_component_age_radial_profile.csv | M33 per-radius residual profiles | 580 |
| paper8_sparc_outer_signed_residual_audit_last3_last5.csv | Signed outer residuals, outermost 3 and 5 points | 266 (two rows per galaxy: tail_N=3 and tail_N=5) |
| paper8_sparc_13_improved_galaxies_characterization.csv | Physical characterization of 13 improved galaxies | 13 |
| paper8_sparc_outer5_sign_flip_watchlist.csv | Galaxies with outer sign flips | 16 |

The code runs in the repository's pinned environment (Python 3.10, NumPy 1.26.4, SciPy 1.15.3). The recomputed baseline median (1.323) differs from Paper 1 (1.312) because the two papers use different derivative conventions in the gate's shear stage; see the recomputation note in Section 5 of the main text.
