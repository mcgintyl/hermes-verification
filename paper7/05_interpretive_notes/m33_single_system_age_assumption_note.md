# M33 and the Single-System Age Assumption

**Archive status:** Interpretive/provenance note only  
**Lane:** M33-Corbelli source-native thick-disk reconstruction  
**Action:** No Hermes rerun. No MOND rerun. No board edits. No age edits. No gravity-comparison reinterpretation.

## Purpose

This note records the interpretive boundary exposed by the M33-Corbelli external board. The locked model results are not changed: Hermes passes globally but remains elevated, with the penalty dominated by outer-disk underprediction; fixed simple MOND fits the same frozen full board cleanly. The purpose here is to frame what the radial localization of the Hermes residual may imply for future Hermes development.

## Locked empirical facts

The M33 board uses the frozen source-native Corbelli reconstruction: pixel-SED/SPS stellar mass, radially varying stellar-mass provenance, helium factor 1.33, H2 prescription, flaring stellar disk, gas half-thickness 0.5 kpc, and Vbul = 0. No post-result reconstruction changes are allowed.

The locked Hermes results are:

| Region / lane | Result |
|---|---:|
| Full board, primary t50 = 6.0 Gyr | chi2nu = 3.126 |
| Full board, older bracket t50 = 7.1 Gyr | chi2nu = 3.296 |
| Main disk, R <= 8 kpc | chi2nu approx 0.64 |
| Main disk, R <= 10 kpc | chi2nu approx 0.80 |
| Outer disk, R > 10 kpc | chi2nu approx 5.79 |
| Outer mean residual | approx -52 to -53 km/s |

The locked MOND comparator result is clean on the same frozen board:

| Comparator | Result |
|---|---:|
| Fixed simple MOND, no M/L tuning | chi2nu = 0.119 |

The baryonic-convention diagnostic does not rescue Hermes. Raising the stellar disk mass improves Hermes but does not produce a clean fit; the outer disk would require a stellar-mass surplus far larger than a plausible SPARC-like convention shift.

The age diagnostic also does not rescue Hermes. Even if the age is allowed to float, Hermes improves only by pushing toward the physical beta ceiling. This means the M33 stress is not simply an age-anchor error.

## Interpretation

Hermes currently applies one board-level age, t50, to the entire galaxy. The local gate varies with radius, but beta is a single global scalar. This is a resource-limited first implementation. It treats the galaxy as one coherent gravitational biography.

M33 may not be such a clean single-system object. The inner/native disk and the outer warped, break-dominated, accretion-sensitive disk may not share the same formation history or dynamical state. The key observation is that Hermes does not fail evenly across M33. It performs cleanly on the main disk and fails primarily in the outer disk.

This makes M33 a plausible boundary case for the single-global-age assumption. The residual structure is spatially associated with an independently documented outer-disk regime: H I warp, kinematic position-angle twist, perturbed outer gas, low-brightness anomalous H I, wide and multiple-peak H I profiles, asymmetric motions, outer stellar disk break, far-outer resolved stellar-population fields, extended disturbed H I structures, and accretion or interaction scenarios.

The cautious causal chain is:

1. Hermes uses one global t50 and one derived beta for the whole board.
2. M33's inner disk is well described by that global treatment.
3. The outer disk is not merely farther out; it is structurally and kinematically distinct in the literature.
4. The Hermes residual appears where that outer complexity begins to dominate.
5. Therefore, the residual may reflect a breakdown of the single-system age assumption and/or the single circular-equilibrium disk assumption, rather than a random full-galaxy failure.

This is a plausible interpretation, not a proof. The strongest supported claim is spatial association plus model-localized stress. The evidence does not justify saying that Hermes “wins by failing,” nor does it justify dismissing MOND’s clean fit.

## MOND comparator status

The MOND result must remain intact. Fixed simple MOND fits the same frozen full board cleanly, with no M/L tuning. That fact is part of the M33 result and should not be softened.

The correct interpretation is not that MOND failed to notice disequilibrium. The present audit cannot establish that. A clean MOND fit may mean MOND is correctly capturing the effective acceleration law on this board, or it may mean MOND is flexible enough to absorb outer-disk kinematic structure that Hermes cannot absorb. Distinguishing those possibilities requires additional tests beyond the present paper, such as independent two-dimensional velocity-field modeling, exclusion masks for non-circular regions, or component-resolved dynamical treatment.

## Dynamic/component-age roadmap

A future Hermes implementation could replace the single global beta with a component-resolved or radial-zone treatment. Such a model would not simply choose a better age after seeing the result. It would require independent data assigning age/formation state by component or radial zone, for example:

- inner disk age / t50,
- outer warped disk age or accretion epoch,
- stellar-break population age,
- gas-accretion or interaction timescale,
- equilibrium quality flags from two-dimensional kinematics.

That extension is outside the executable scope of the present work. It would require new data standards and new pre-declared rules. It should not be retroactively applied to the frozen M33 board.

One important caution: the previous age sweep showed that making the entire board younger does not solve the M33 outer residual. Therefore, “component age” alone may not be sufficient. The future extension likely needs both component-resolved age and component-resolved dynamical status, especially for warped or accreting outer gas.

## Paper-safe language

> The M33 result suggests a limitation of the current single-global-age implementation. Hermes fits the main disk cleanly but underpredicts the outer disk, and the residual zone overlaps independently documented warp, disk-break, disturbed-H I, and accretion-sensitive structures. This does not overturn the MOND comparator result, which remains clean on the same frozen board. Instead, M33 identifies a likely boundary condition: galaxies with multi-origin or dynamically disturbed outer components may require a component-resolved age and dynamical-status treatment rather than a single board-level t50. In this sense, the M33 stress result motivates a dynamic-aging extension of Hermes, while remaining outside the executable scope of the present work.

## Recommended archive classification

**M33-Corbelli source-native thick-disk reconstruction:** valid external Format B execution; Hermes pass/elevated; MOND clean; residual localization indicates outer-disk boundary case.

**Interpretive label:** Single-system age-assumption boundary case.

**Confidence assessment:**

| Claim | Confidence | Reason |
|---|---|---|
| Hermes residual is localized in the outer disk | High | Direct from locked radial chi2 masks |
| The outer residual overlaps documented M33 outer-disk complexity | High | Warp, break, disturbed H I, and far-outer fields are independently published |
| Single-global-age treatment may be insufficient for M33 | Moderate | Plausible from localized stress, but not directly proven |
| Dynamic/component-age treatment would solve M33 | Low to moderate | Possible roadmap, but age-only diagnostics show beta ceiling remains limiting |
| MOND’s clean result is invalid because of disequilibrium | Not established | MOND clean fit is locked; disequilibrium interpretation requires further kinematic tests |

## Final doctrine

M33 should be presented as a boundary case and research-roadmap case, not as a solved Hermes success. The result is valuable because it separates three facts that can coexist:

1. Hermes can port non-catastrophically to a source-native external board.
2. Hermes can still show a real, localized stress where the board enters a physically distinct outer regime.
3. MOND can fit the same frozen board cleanly.

The paper-level lesson is that a single board-level t50 is a strong simplifying assumption. It may be adequate for coherent main disks, but galaxies with warped, accretion-sensitive, or multi-origin outer components may require component-resolved age and dynamical-status treatment before a rigid historical gravity model can be interpreted at full radial range.
