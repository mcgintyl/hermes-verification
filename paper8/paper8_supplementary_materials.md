# Supplementary Materials

**The Diagnostic Horizon: A Wear-Activated Density-Release Refinement of the Hermes Gravity Relation**

**Louis Albert McGinty**

mcgintyl@grinnell.edu

Independent Researcher, Team Hermes

AI-assisted research; full disclosure in main text

---

# Supplement A: Candidate Rejection Path and Refinement Constraints


## 1. Motivation

The Hermes equation (McGinty 2026a) predicts galaxy rotation curves from the observed baryonic acceleration profile, with its retained-coherence term fixed by stellar age and peak baryonic acceleration. Applied to 133 galaxies from the SPARC database with verified stellar ages, the equation achieves a median reduced chi-squared of 1.312 with zero adjustable parameters per galaxy. It is competitive with MOND in typical performance and produces substantially fewer catastrophic failures: 84% of galaxies score below 5.0, compared to 71% for MOND.

That baseline performance was established on a locked framework. The equation's constants were derived from staged tests on galaxy subsets and then frozen. No constant was adjusted per galaxy, and the full 133-galaxy sample was tested without modification. This rigidity is what makes the equation's successes meaningful and its failures informative.

Paper 7 (McGinty 2026b) tested the equation on M33, one of the nearest large spiral galaxies and a system not included in the original SPARC sample. The inner disk, within 10 kpc of center, fit cleanly: chi-squared-nu = 0.80. The outer disk, beyond 10 kpc, did not. The outer region produced chi-squared-nu = 5.79, carrying approximately 86% of the total chi-squared budget. The equation was not failing everywhere. It was failing in a specific regime: the extended, low-density outer disk where baryonic acceleration transitions from moderate to very low values.

A separate investigation examined whether the Hermes wear exponent psi approached a smooth population distribution capable of supporting power-law behavior in inverse coherence. Instead, psi is structured by peak baryonic acceleration quartile: galaxies in different density classes produce distinct wear distributions rather than a single continuous population. The kink in the population distribution occurs in the same density regime where M33's per-galaxy residual concentrates.

These are unlikely to be two unrelated problems. The population-level density kink and the per-galaxy radial failure appear to reflect the same structural limitation viewed at two scales. The baseline equation handles much of the high-density and low-density SPARC population competitively, but these diagnostics suggest that the transition regime, where a galaxy moves from one density domain to the other, requires a more selective coupling mechanism for the retained coherence to express correctly.

This paper presents and tests a zero-parameter density-transition refinement. The goal is not to replace the Paper 1 equation. It is to preserve the locked chassis while adding a constrained selector for the density-seam regime exposed by Paper 7 and the population-level diagnostics. The refinement must satisfy the same rigidity that made Paper 1 credible: no per-galaxy parameters, no refitting of existing constants, no material degradation of the 133-galaxy SPARC board. These constraints are specified in Section 2 and enforced throughout.


## 2. Refinement Constraints

Any modification to a zero-parameter equation risks converting it into a fitting exercise. The Hermes equation's credibility rests on its rigidity: once the constants were locked, they were never adjusted per galaxy. Paper 8 inherits that discipline. The refinement is only meaningful if it survives under the same constraints that made Paper 1 credible. These constraints are not procedural decoration. They are the evidentiary standard.

The following requirements governed every candidate tested in this study.

First, no per-galaxy parameters. The refinement may not introduce any quantity that is adjusted, optimized, or selected on a galaxy-by-galaxy basis. Every galaxy receives the same equation with the same constants applied identically.

Second, no refitting of Paper 1 constants. The baseline constants, the coherence ceiling pi, the light-speed normalization in psi, the floor 1/sqrt(2pi), and the gate threshold a_k = 1585, are inherited and frozen. The refinement may introduce new operators, but it may not modify existing constants.

Third, use only frozen Hermes constants where possible. Any new constant introduced by the refinement should be derived from the existing constant set rather than independently optimized against the sample. This limits the effective degrees of freedom to the architectural choices, not to fitted scalars.

Fourth, apply the refinement to the positive retained-coherence term, not to the net age score. The age score beta has a physical floor at -1/sqrt(2pi) representing fully eroded coherence. Any density-release mechanism must act on the surviving positive coherence (pi times exp(-psi)) rather than on the net score, to avoid amplifying the floor into unphysical regimes.

Fifth, preserve or improve SPARC aggregate behavior without introducing casualties. In this study, a material casualty is defined as any galaxy degraded by more than 0.5 chi-squared-nu. The refinement must also avoid creating new catastrophic cases (chi-squared-nu > 5 or chi-squared-nu > 10). Improvement in the median is desirable but secondary to safety.

Sixth, avoid Baryonic Tully-Fisher Relation collision. The asymptotic flat rotation velocity must not depend on galaxy-specific core density. Any density-release term that anchors the outer-disk behavior to the peak baryonic acceleration g_98 risks injecting surface-brightness scatter into the BTFR, which is empirically one of the tightest scaling relations in extragalactic astronomy.

Seventh, preserve M33 as an external stress test, not as a tuning target. M33 is used to evaluate the refinement's behavior on a system with known outer-disk complications, not to optimize the refinement's form. Results on M33 are reported but do not influence the equation's constants or structure.

Eighth, report all failed candidates. The rejection path constrains the surviving operator at least as strongly as the surviving operator constrains the data. Each failure identifies a specific structural requirement that the final form must satisfy. Omitting them would obscure the evidentiary basis for the surviving equation.


## 3. Candidate Rejection Path

The final refinement was not selected first and justified after the fact. It was discovered by constrained elimination. Five candidates were tested against the requirements of Section 2. Each one failed, and each failure taught the surviving equation something about its own necessary structure.

The first candidate replaced the raw peak baryonic acceleration g_98 in the wear exponent psi with a logarithmically compressed version, intended to smooth the population-level density kink identified in the motivation. The compression succeeded at the population level: the wear distribution became substantially less structured by density quartile. But it substantially degraded SPARC performance, producing 53 to 58 galaxies with chi-squared-nu above 5.0 depending on the compression variant. The high-wear branch of the equation, which governs the oldest and densest galaxies near the coherence floor, is load-bearing. Smoothing it away destroyed fits that were already working. The lesson: the density kink in psi is not a scaling artifact that can be removed by transformation. It reflects genuine structural heterogeneity in the galaxy population.

The second candidate left psi intact and instead introduced a radial density-release factor eta(R) that allowed retained coherence to express more strongly in low-density outer regions. A sign guard set eta = 1 whenever the net age score beta was negative, preventing the release from amplifying the coherence floor. This candidate produced the strongest M33 result of any version tested: the outer-disk chi-squared-nu dropped from 5.79 to 0.72. However, it created 8 material SPARC casualties and 3 new catastrophic cases. The failures traced to two structural flaws. The density-release anchored its outer asymptote to the galaxy-specific peak acceleration g_98, creating a BTFR-collision risk (Section 2, sixth). And the hard sign guard at beta = 0 created a topological discontinuity: a galaxy at beta = +0.001 received the full density release while a galaxy at beta = -0.001 received none, with no physical basis for the boundary. The lesson: the density seam is real and exploitable, but coupling it to galaxy-specific core density is structurally wrong, and a binary sign guard is not a physical transition mechanism.

The third candidate addressed both flaws. It decoupled the outer asymptote from g_98 by replacing it with a universal numerator a_u = a_k/pi, derived from the existing gate threshold. It applied the density release to the positive decay term pi times exp(-psi) rather than to the net age score, eliminating the sign guard entirely: as psi grows large, the positive term decays toward zero and the release has progressively less to amplify, so the floor is reached continuously without a discontinuity. This version reduced SPARC casualties from 8 to 1 and created no new catastrophic cases. The single remaining casualty, UGC 05750, is a very low-wear, high-retention system near the coherence ceiling. The refinement added outer support to a galaxy that did not need it because it had barely begun the decay process. The lesson: the BTFR decoupling and positive-term application are both necessary components of any viable refinement, but the correction also requires a mechanism to suppress activation in near-newborn systems.

The fourth candidate proposed that the physical core and the low-density branch represent orthogonal domains whose coherence contributions should combine via Pythagorean sum rather than linear addition. The resulting expression, P(R) = exp(-psi) times sqrt(pi^2 + W(psi)^2 times pi times a_k / Gamma(g_bar)), is mathematically elegant: it eliminates the logistic transition function, the max floor, and the nested eta structure in favor of a single geometric expression. It produced 15 material SPARC casualties and 4 new catastrophic cases. The Pythagorean sum did not shut off fast enough in the dense core, leaking a parasitic boost into galaxies that were already well-fit. More critically, the pi cap on total P(R) strangled the low-density asymptote: in M33's outer disk, P(R) saturated at pi regardless of the component age, erasing the entire component-age diagnostic axis. The lesson: the core and the low-density branch are not continuously superposed geometries. Smooth blending across the full density range fails empirically. The data prefer an operator that keeps the two regimes distinct and selects between them rather than averaging them.

The fifth candidate reframed the problem as a thermodynamic phase transition between a bound physical state and an unbound low-density state. The probability of being in the unbound state was modeled as the product of temporal wear activation W(psi) and a Fermi-Dirac-style spatial transition S(g). The low-density multiplier used Pythagorean geometry only within the unbound state, avoiding the continuous-superposition error. This version produced the strongest M33 result among the phase-transition candidates (outer chi-squared-nu approximately 1.02 at the 2 Gyr component-age bracket), but it created 14 material SPARC casualties and 2 new catastrophic cases. An isolation test confirmed that replacing the inherited logistic width (0.30) with the theoretically motivated value 1/pi made no difference to the casualty count; the structural form itself was too permissive. The casualties clustered in the acceleration range between the kinematic collapse threshold and the decoherence boundary, a transition zone where the phase switch had activated but the low-density branch was not yet dominant. The lesson: the phase-transition framing captures part of the density-seam behavior, but smooth amplitude scaling in the transition zone adds support where the data do not require it. The surviving operator must withhold the low-density contribution until the low-density branch is actually stronger than the local three-dimensional field.

The surviving refinement is therefore not the most elegant candidate tested. It is the only candidate that satisfies all eight constraints simultaneously. Its wear activation suppresses the release in near-newborn systems. Its universal numerator decouples the asymptote from core density. Its positive-term application avoids the sign-guard discontinuity. Its max-selected density contrast withholds the low-density branch until it dominates the local field. Each of these features was forced into the equation by a specific empirical failure, not by theoretical preference. The equation's structure was shaped by 133 galaxies of constraint pressure.


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

**Logarithmic compression.** The baryonic acceleration at each radius is compressed to regularize extreme density contrasts between galaxy cores and outer disks:

Gamma(g) = a_k * ln(1 + max(g, 0) / a_k)

where a_k = 1585 (km/s)^2/kpc is the gate threshold inherited from Paper 1. For high accelerations, Gamma grows logarithmically rather than linearly, preventing the density contrast from dominating the operator in bulge-heavy systems. For low accelerations, Gamma approximately equals g, preserving the linear scaling in the outer disk.

**Wear activation.** The refinement is suppressed in galaxies that have accumulated little temporal wear:

W(psi) = 1 - exp(-2*pi^2 * psi)

For near-newborn systems (psi close to zero), W is approximately zero and the refinement has no effect. For systems with significant wear, W saturates rapidly toward unity. The coefficient 2*pi^2 is a closed-form expression of the inherited constant pi; interpretation is deferred to Supplement E. This component was introduced after the third candidate (Section 3) revealed that the refinement over-activated in low-wear systems near the coherence ceiling.

**Spatial transition.** A logistic transition function governs where the density release activates along the radial profile:

S(g) = [1 + exp((g - a_k) / (0.30 * a_k))]^(-1)

When g is well above a_k, S is near zero and the release is suppressed. When g is well below a_k, S is near unity and the release is fully active. The transition center and width (a_k and 0.30 * a_k) are inherited from the Paper 1 gate construction. No new parameters are introduced.

**Universal numerator.** The density contrast that drives the release is anchored to a universal constant rather than to the galaxy-specific peak acceleration:

a_u = a_k / pi (approximately 504.5 (km/s)^2/kpc)

This decoupling was required by the second candidate's failure (Section 3): anchoring the outer asymptote to g_98 created a BTFR-collision risk. The value a_k/pi is derived from existing constants. When Gamma(g_bar) falls below a_u, the contrast term enters the square-root branch sqrt(a_u/Gamma). The interpretation of this branch is discussed in Supplement E.

**The local density-release operator.** At each radius, the operator evaluates whether the low-density branch exceeds the local three-dimensional field:

eta_U(R) = min[pi, 1 + S(g_bar(R)) * (max(1, sqrt(a_u / Gamma(g_bar(R)))) - 1)]

The max function selects the larger of two values: unity, corresponding to no density-release beyond baseline, or sqrt(a_u/Gamma), corresponding to the low-density branch. When Gamma exceeds a_u, the max evaluates to unity and the release contributes nothing beyond baseline. When Gamma falls below a_u, the low-density branch dominates and the release activates. The logistic S(g) modulates the transition. The min function caps the total operator at pi, preventing the local multiplier from growing without bound in the deep outer disk. Unlike the rejected Pythagorean candidate (Section 3), the cap is applied to the multiplier, not to the total positive coherence term.

**The wear-gated operator.** The local operator is modulated by the global wear activation:

eta_WA(R) = 1 + W(psi_sys) * [eta_U(R) - 1]

When W is zero (near-newborn systems), eta_WA equals 1 regardless of local conditions. When W is unity (well-worn systems), eta_WA equals eta_U. The wear gate uses the systemic wear variable psi_sys computed from the galaxy's global stellar age, not the local psi(R). This prevents the density release from acting on galaxies that have not yet accumulated significant temporal stress.

**Component-age extension.** For galaxies with evidence of distinct dynamical components, the positive decay term can use a spatially varying effective age:

psi(R) = 2*pi * t_eff(R) * g_98 / c

For the standard SPARC execution, psi(R) = psi_sys at all radii. For the M33 component-age diagnostic, the positive decay term uses psi(R) derived from t_eff(R), with t_eff equal to the global stellar age for the inner disk (R <= 10 kpc) and a younger effective age for the outer disk (R > 10 kpc). The wear gate W remains defined by the systemic bracket. The component-age extension does not modify any constant or operator; it modifies only the age input for regions where independent observations motivate a multi-component interpretation.

**Summary of constants.**

| Constant | Value | Role | Status |
|----------|-------|------|--------|
| pi | 3.14159... | Coherence ceiling, operator cap | Inherited |
| e | 2.71828... | Exponential decay | Inherited |
| c | Paper 1 units | Light-speed normalization in psi | Inherited |
| 1/sqrt(2pi) | 0.3989... | Coherence floor | Inherited |
| a_k | 1585 (km/s)^2/kpc | Gate threshold, transition center | Inherited |
| 0.30 | -- | Logistic width ratio | Inherited |
| a_u = a_k/pi | 504.5 (km/s)^2/kpc | Universal density-release numerator | Derived |
| 2pi^2 | 19.74 | Wear activation rate | Derived |

No constant in the refinement was newly fitted to the 133-galaxy SPARC board or to M33. The inherited constants are taken unchanged from Paper 1. The two derived constants are closed-form expressions of inherited constants. No per-galaxy parameters are introduced at any stage.

---

# Supplement B: Signed Residual Audit

The wear-activated eta refinement was audited for directional bias in the outer disk. Zero-casualty chi-squared results can mask systematic sign flips if a correction shifts residuals from under-prediction to over-prediction without changing the variance.

For the outermost 5 measured points per galaxy, the median signed residual shifted from -20.51 km/s (baseline) to -10.80 km/s (corrected). The median galaxy remains under-predicted after correction. Of 133 galaxies, 16 experienced sign flips in the outermost 5 points, all from under-prediction to over-prediction. Five galaxies exceeded +5 km/s material overshoot (DDO 161, UGC 06983, UGC 12732, NGC 2366, UGC 05005). No galaxies flipped from over-prediction to under-prediction.

The correction reduces outer-disk starvation without systematically overshooting the sample. The one-sided drift is documented as a watchlist item. Full per-galaxy signed residual data is provided in paper8_sparc_outer_signed_residual_audit_last3_last5.csv.

---

# Supplement C: Improved-Galaxy Characterization

Thirteen galaxies improved by more than 0.5 chi-squared-nu under the wear-activated eta. These galaxies share a physically coherent profile.

Compared to the remaining 120 galaxies, the improved class has: more negative baseline outer residuals (median -41.83 vs -16.76 km/s), greater radial extent beyond the knee (median Rlast/rknee 6.95 vs 2.44, with 12/13 exceeding 2.0), higher outer gas fraction (0.518 vs 0.204), lower median disk surface brightness (1.12 vs 14.92 L/pc^2), lower outer-density ratio (0.036 vs 0.148), and enrichment in Chimera-flagged systems (6/13 vs 21/120).

The improvement targets extended, gas-rich, low-density outer disks where baseline Hermes was structurally under-supporting the tail. Two galaxies (NGC 6674 and NGC 0289) improve numerically but remain poor fits and should not be treated as clean rescues.

Full characterization data is provided in paper8_sparc_13_improved_galaxies_characterization.csv.

---

# Supplement D: M33 Component-Age Sweep and Literature Notes

## Sweep Results

The wear-activated eta was combined with a component-resolved age model on M33, using the inner disk (R <= 10 kpc) at the global stellar age and the outer disk (R > 10 kpc) at younger effective ages. Neither eta alone nor component age alone closes the outer-disk residual. Together they largely resolve it.

Primary bracket (inner t50 = 6.0 Gyr):

|     Outer age (Gyr)     | Full chi2_nu | Inner (R <= 10) chi2_nu | Outer (R > 10) chi2_nu | Far outer (R >= 15) chi2_nu |
| :---------------------: | :----------: | :---------------------: | :--------------------: | :-------------------------: |
| 6.0 (uniform, eta only) |    1.571     |          0.758          |         2.505          |            2.999            |
|           3.0           |    1.192     |          0.758          |         1.691          |            2.141            |
|           2.0           |    1.076     |          0.758          |         1.442          |            1.870            |
|           1.5           |    1.021     |          0.758          |         1.323          |            1.739            |
|           1.0           |    0.967     |          0.758          |         1.208          |            1.611            |

Conservative bracket (inner t50 = 7.1 Gyr):

|     Outer age (Gyr)     | Full chi2_nu | Inner (R <= 10) chi2_nu | Outer (R > 10) chi2_nu | Far outer (R >= 15) chi2_nu |
| :---------------------: | :----------: | :---------------------: | :--------------------: | :-------------------------: |
| 7.1 (uniform, eta only) |    1.740     |          0.800          |         2.820          |            3.324            |
|           3.0           |    1.214     |          0.800          |         1.689          |            2.138            |
|           2.0           |    1.098     |          0.800          |         1.441          |            1.868            |
|           1.5           |    1.043     |          0.800          |         1.322          |            1.737            |
|           1.0           |    0.989     |          0.800          |         1.207          |            1.609            |

Full results in paper8_m33_wear_eta_component_age_summary.csv.

Baseline Hermes without eta: component age alone (even at newborn outer limit) reaches only chi2_nu = 4.05 in the outer disk. The baseline ceiling is structural.

## Literature Notes

The required outer-component age of 1 to 2 Gyr is consistent with published gas dynamics work, when framed as the effective coherence or dynamical age of the outer gas layer rather than the stellar half-mass age of the whole outer disk.

Putman et al. (2009) found M33's HI mass extends well beyond the star-forming disk with warps, arcs, and a southern filament. Their orbit analysis favored tidal disruption by M31 at 1 to 3 Gyr ago.

Semczuk et al. (2018) modeled a recent M31-M33 interaction at approximately 2 Gyr, reproducing aspects of the gaseous warp when tidal forcing was combined with ram pressure from M31's hot gas halo.

Corbelli and Burkert (2024) used Gaia-EDR3 proper motions to argue against a recent close M31 passage, favoring first infall. They attribute M33's outer-disk misalignment to recent cold gas accretion from a cosmic filament, with an estimated inflow scale of order 1 solar mass per year.

Barker et al. (2011) found an age-gradient reversal in M33's outer disk: the field at 9.1 kpc has a mean age of approximately 3 +/- 1 Gyr, while the farther field at 11.6 kpc is older at approximately 7 +/- 2 Gyr.

The literature does not give a single clean answer. Older work supports 1 to 3 Gyr tidal disruption. Newer proper-motion work pushes toward first infall plus filament accretion. Both support a recently disturbed or refreshed outer gas component. The M33 component-age result is defensible as a diagnostic bracket, not as a precise age assignment.

---

# Supplement E: Dimensional-Transition Interpretation (Provisional)

This supplement describes a provisional interpretive reading of the wear-activated eta operator. These notes are not required to evaluate the empirical results in the main text. They are provided for readers interested in the structural implications of the operator's form.

The rejected candidates (Supplement A) show that smooth superposition of the core and low-density regimes fails empirically, while sharp branch-selection survives. The surviving operator contains a max function that selects between two scaling behaviors:

When Gamma(g_bar) exceeds a_u: the multiplier is unity and gravity follows the baseline three-dimensional field.

When Gamma(g_bar) falls below a_u: the multiplier enters the square-root branch sqrt(a_u/Gamma), which produces an effective force scaling proportional to 1/R rather than the Newtonian 1/R^2.

A 1/R force law is the scaling expected from flux spreading across a two-dimensional surface rather than a three-dimensional volume. The SPARC and M33 audits are consistent with a dimensional-transition interpretation in which the equation selects between three-dimensional and two-dimensional propagation geometries depending on local density.

The region between a_u (corresponding to raw g_bar approximately 594) and a_k (1585) represents a transition zone where the spatial switch S(g) has activated but the low-density branch is not yet dominant. Candidates that boosted this zone created SPARC casualties. The surviving operator withholds support until the low-density branch is actually stronger than the local field.

The remaining components of the operator also correspond to specific geometric quantities. The wear activation coefficient 2*pi^2 equals the surface area of a unit 3-sphere. The operator cap pi equals the ratio of a circle's circumference to its diameter. These identities are noted as structural observations. Whether they reflect deeper physical significance or coincidence remains an open question.

This interpretation is structurally supported by the candidate rejection pattern but is not independently established. It should be treated as provisional. We do not claim direct observation of a two-dimensional physical substrate. Rather, we report a mathematical isomorphism: the empirical operators required to preserve the SPARC sample, specifically the max selection floor and the low-density square-root branch, map onto the geometry of a dimensional collapse. The data selected for this mathematical signature; the dimensional-transition framework provides its conceptual interpretation.

---

# Supplement F: Production Data Package and CSV Manifest

All tabular outputs and verification materials required to audit the reported comparisons are provided in the accompanying data package.

| Filename | Contents | Rows |
|----------|----------|------|
| paper8_sparc_133_model_comparison.csv | Per-galaxy results: baseline, MOND, wear-activated eta chi2_nu, deltas, flags | 133 |
| paper8_sparc_per_radius_baseline_mond_wear_eta.csv | Per-radius predictions for all three models | ~3073 |
| paper8_sparc_aggregate_metrics.csv | Summary statistics for all models | 3 |
| paper8_iteration_comparison_metrics.csv | Aggregate metrics for all six candidate iterations | 6 |
| paper8_m33_wear_eta_component_age_summary.csv | M33 results by mask and outer age | Multiple |
| paper8_m33_wear_eta_component_age_radial_profile.csv | M33 per-radius residual profiles | Multiple |
| paper8_sparc_outer_signed_residual_audit_last3_last5.csv | Signed outer residuals, outermost 3 and 5 points | 266 (two rows per galaxy: tail_N=3 and tail_N=5) |
| paper8_sparc_13_improved_galaxies_characterization.csv | Physical characterization of 13 improved galaxies | 13 |
| paper8_sparc_outer5_sign_flip_watchlist.csv | Galaxies with outer sign flips | 16 |

Computational environment: Python 3.x, NumPy, SciPy. The recomputed baseline median (1.323) differs from Paper 1 (1.312) because the two papers evaluate the gate's shear stage with different derivative conventions: Paper 1 used a chain-rule derivative of $V_{\rm sm}$, while this paper uses a log-grid finite difference of $\ln V_{\rm sm}$. Each convention reproduces its corresponding per-galaxy score vector to better than $3 \times 10^{-13}$ absolute (see Paper 1 v6, Methods).

The Hermes verification repository is available at github.com/mcgintyl/hermes-verification.
