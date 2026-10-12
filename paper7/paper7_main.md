# Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing

*A portability audit of the Hermes gravity framework*

**Louis McGinty**

mcgintyl@grinnell.edu

Independent Researcher, Team Hermes

AI-assisted research; full disclosure in AI Disclosure section

## Abstract

The Hermes gravity relation predicts galaxy rotation curves from observed baryonic structure and a galaxy-age term with no per-galaxy free parameters. Previous work evaluated Hermes on a 133-galaxy age-verified subset of the SPARC database. This paper documents the attempt to test Hermes on independent, non-SPARC data.

The attempt exposed a data standardization problem. Fair cross-model rotation-curve testing requires board-complete inputs: observed velocities, uncertainties, and independent per-radius gas, stellar disk, and bulge components with documented mass-to-light ratio conventions. We audited public databases, survey products, author files, and doctoral theses. Most candidate paths failed as model boards before any model could be run, each for a different reason: calibration-set overlap, incompatible velocity targets, missing baryonic components, privately held data, and convention mismatches. The BIG-SPARC project independently identifies this limitation but is not yet publicly available.

One external board was reconstructed: M33, from the published baryonic mass model of Corbelli et al. (2014). On this Format B board, Hermes produced a non-catastrophic but elevated fit (χ²ν ≈ 3.2), while Modified Newtonian Dynamics fit cleanly (χ²ν ≈ 0.1). The Hermes residual concentrated in the outer disk, a region independently documented as dynamically complex. Baryonic-convention and age diagnostics narrowed but did not eliminate the gap.

The result is one external board, not a general verdict. Its value lies in three findings: the public data ecosystem is less board-complete than commonly assumed, Hermes failures can localize to physically meaningful boundaries, and rigid models can expose structured input problems that flexible models may absorb. External tests of any baryonic gravity model will remain as much a test of the data infrastructure as of the model itself until standardized board-complete datasets become available.


# 1. Introduction

A model calibrated on a single dataset must eventually be tested on independent data. The Hermes gravity relation (McGinty 2026) predicts galaxy rotation curves from observed baryonic structure and a galaxy-age term, with no free parameters fitted per galaxy. Its initial test base was the SPARC database (Lelli, McGaugh & Schombert 2016), a 175-galaxy sample with near-infrared photometry, HI/Hα rotation curves, and published baryonic mass models. This paper documents the attempt to test Hermes outside that calibration set.

The obvious next step, applying Hermes to galaxies from other surveys, proved harder than anticipated. Fair cross-model testing requires what we term board-complete data: observed rotation curves with independent per-radius gas, stellar disk, and bulge components, uncertainties, and documented mass-to-light ratio conventions. We audited public databases, survey products, author-supplied files, doctoral theses, and small-sample studies. Most candidate sources failed to meet these requirements before any model could be run.

Hermes is a model of how a galaxy's gravity changes as it ages. To test it, or to test any model that builds gravity from visible matter, you have to see that matter in its separate parts: the gas, the stellar disk, and the bulge. Each is a different kind of matter in a different place, measured a different way. Gas comes from radio observations of neutral hydrogen and sits mostly in the outskirts, where the rotation-curve discrepancy is largest. The disk and bulge come from near-infrared light and describe how the stars are arranged, which is itself one of the things that shifts as a galaxy ages. None of the three can be recovered from the others, because they are not versions of one measurement. They are independent measurements of independent things.

That is why a dataset missing a component is not a lower-quality input an analyst can patch. It is a different test. The model has to reproduce the observed motion from the galaxy's actual visible mass. Remove a component and you are no longer comparing the model to the galaxy, you are comparing it to a guess about the part you could not measure. A model with adjustable parameters can hide some of that guess inside its fit. Hermes has none, so the guess goes straight into the prediction as error. Board-complete, in this paper, means all three components present and independently measured. Anything short of that is not a weaker version of the same test. It is a different test under the same name.

One external board was successfully reconstructed: M33, built from the published baryonic mass model of Corbelli et al. (2014). M33 was tested under the same frozen methodology used for SPARC, with Modified Newtonian Dynamics (MOND; Milgrom 1983) as a comparator on the same inputs. Hermes produced a non-catastrophic but elevated result. MOND fit cleanly. The Hermes residual was not spatially random; it concentrated in M33's outer disk, a region independently documented as dynamically complex.

This paper is not a broad validation of Hermes outside SPARC. It is a documentation of what happened when an independent team attempted external validation: the infrastructure barriers encountered, the single external test completed, and the diagnostic information that emerged from both successes and failures. The public data ecosystem is less board-complete than commonly assumed, M33 provides one hard-won external case, and the structure of the Hermes result, where and why it fails, becomes part of the finding.


# 2. Methods

Hermes takes a galaxy's observed baryonic mass structure and a galaxy-age term as inputs, then predicts the rotation curve without fitting per-galaxy parameters. In practice, this requires per-radius baryonic components at each sampled radius: observed rotation velocity, observational uncertainty, gas contribution, stellar disk contribution, and bulge contribution where present. If any component is missing or convention-incompatible, the test becomes a reconstruction problem rather than a direct model evaluation. The equation, pipeline, and age-source hierarchy used here follow McGinty (2026). External boards required new provenance checks, but no model terms were added and no per-galaxy parameters were fitted.

Hermes applies an age-dependent correction to the baryonic (Newtonian) acceleration at each radius. The full derivation is in McGinty (2026); the equations below show only what is being tested here.

$$g_{\text{model}}(R) = g_{\text{bar}}(R) \times [1 + \beta \cdot \varphi(R)] \tag{1}$$

$$\beta = \pi \cdot \exp(-\psi) - 1/\sqrt{2\pi} \tag{2}$$

$$\psi = 2\pi \cdot t_{50} \cdot g_{98} / c \tag{3}$$

Here $g_{\text{bar}}(R)$ is the acceleration from the combined visible mass (gas, stellar disk, and bulge); $\varphi(R)$ is the gate, a profile running from 0 to 1 built from the shape of $g_{\text{bar}}(R)$ by the fixed procedure in McGinty (2026); $\beta$ is the global age score; $\psi$ is the dimensionless wear variable; $t_{50}$ is the age at which the galaxy had formed half its present-day stars; $g_{98}$ is the 98th-percentile baryonic acceleration; and $c$ is the speed of light. The constants are global and fixed; none is fitted per galaxy.

The following subsections define the frozen test chassis, board-completeness classes, SPARC-overlap screen, mass-to-light convention, MOND comparator, and execution sequence used throughout the audit.

**2.1 Frozen chassis.** All external tests use the same Hermes implementation applied to SPARC. No element of the equation, pipeline, age assignment, or diagnostic protocol was changed between SPARC evaluation and external testing. No per-galaxy parameter tuning is permitted. No post-result modifications to the board or model are permitted. If the frozen implementation fails on an external board, the failure is reported as-is. As part of post-publication reproducibility verification, the M33 external board was found to have been executed with a separately-coded copy of the pipeline whose gate function contained three incidental numerical deviations from the SPARC implementation (knee-radius definition, Savitzky-Golay edge mode, and curvature regularization). These deviations are not required by the Format B reconstruction; they have been corrected, and all M33 results reported here use the single canonical gate. The correction changes M33's full-board χ²ν from 3.13 to 3.22 and alters no conclusion of this paper. Both the original and corrected M33 results are archived (github.com/mcgintyl/hermes-verification).

**2.2 Board completeness classification.** External datasets vary in what they publish. We classify candidate boards by schema completeness. Format A: the source provides machine-readable per-radius arrays comparable to SPARC's rotmod format and can be tested directly. Format B: the source publishes sufficient numeric ingredients to reconstruct the required components without figure digitization, with documented provenance. Format C: the source is missing at least one critical numeric component and cannot be tested without author-supplied data, figure digitization, or external reconstruction.

**2.3 SPARC overlap discipline.** Any galaxy claimed as a non-SPARC external test must be verified against the full 175-galaxy SPARC manifest by filename, known aliases, and coordinates. SPARC-overlap galaxies are excluded from external validation claims. Where SPARC-overlap galaxies were tested with foreign baryonic decompositions, those results are reported separately as pipeline robustness comparisons.

**2.4 Mass-to-light ratio convention.** Hermes was calibrated on SPARC using the database's provided baryonic component velocities, where stellar disk and bulge components are supplied at unit mass-to-light normalization. For external boards where the raw near-infrared light profile is available, the SPARC-calibration convention is applied as the primary test. Where the source publishes only a stellar mass profile from population synthesis modeling, the source-native mass-to-light ratio is used and the convention difference is documented. In no case is the mass-to-light ratio optimized per galaxy.

**2.5 MOND comparator.** Modified Newtonian Dynamics (Milgrom 1983) is run on every frozen board as a comparator, using the same baryonic inputs, radial grid, and uncertainties as Hermes. We use the simple interpolation function, applied to the same $g_{\text{bar}}(R)$ Hermes receives:

$$g_{\text{MOND}}(R) = g_{\text{bar}}(R) \times \nu(g_{\text{bar}}/a_0) \tag{4}$$

$$\nu(y) = \tfrac{1}{2} + \sqrt{\tfrac{1}{4} + 1/y}, \qquad a_0 = 1.2 \times 10^{-10}\ \text{m s}^{-2} \tag{5}$$

where $y = g_{\text{bar}}/a_0$, and equation (5) is the boost form of the standard simple interpolation function $\mu(x) = x/(1+x)$. The acceleration scale $a_0$ is held fixed and no mass-to-light ratio adjustment is applied. If both Hermes and MOND produce elevated residuals on the same board, the result is treated first as a board-closure or input-compatibility problem before being interpreted as model-specific failure.

**2.6 Execution sequence.** A model test begins only after the board is frozen. The sequence is: source audit, SPARC-overlap screening, schema classification, board reconstruction where needed, mass closure verification, board freeze with provenance lock, Hermes execution, MOND execution on the same frozen board, and diagnostic triage. The Hermes result is archived before the MOND comparator is run.

**2.7 Scoring convention.** Fits are scored with reduced chi-squared (χ²ν), using the same uncertainty convention as our first paper (McGinty 2026), including its fixed error-bar floor, applied here without change. No per-galaxy parameters are fitted, so no fitted-parameter penalty is applied. Because χ²ν depends on that convention, we also report the plain velocity residuals in km/s (Table 2), so the result does not rest on the chi-squared alone.


# 3. The Data Standardization Problem

Testing any gravity model on rotation curves requires a board-complete set of inputs at each sampled radius: the observed rotation velocity, its uncertainty, and independent contributions from gas, stellar disk, and bulge. These components are not interchangeable with summary statistics, integrated masses, or model-inferred circular velocities. This requirement is not specific to Hermes. MOND, dark matter halo models, and any baryonic-input gravity framework face the same constraint: without independent per-radius components, the test is incomplete.

SPARC (Lelli, McGaugh & Schombert 2016) provides exactly this schema for 175 galaxies. We attempted to locate comparable data outside SPARC. The following summarizes what we found.

**McGaugh pre-SPARC files.** Eight low-surface-brightness galaxies from earlier work by the SPARC collaboration were available in SPARC-compatible format. All eight, however, are already present in the SPARC rotmod manifest. The lane succeeded as a schema test but failed as an independence test: the data are clean, but they are not independent of the calibration set.

**MaNGA DynPop.** The MaNGA DynPop survey provides stellar population histories and circular velocity curves for thousands of galaxies (Zhu et al. 2025; Lu et al. 2023). Its age products are strong. However, the velocity curves are derived from Jeans Anisotropic Modelling of stellar kinematics, not from observed HI or Hα rotation. The survey does not provide independent per-radius gas, disk, or bulge components. A JAM-inferred circular velocity curve is a model-derived gravitational potential estimate, not an observed rotation measurement. This distinction blocks MOND as well as Hermes.

**Noordermeer early-type disks.** Noordermeer et al. (2007) present observed HI rotation curves for 19 early-type disk galaxies, a subset of which are not in SPARC. However, the per-radius baryonic decomposition is not published in the journal paper, its supplement, or the associated VizieR catalogue. Mass decompositions appear only as plotted figures in the lead author's doctoral thesis (Noordermeer 2006). While figure digitization is accepted practice in the field, it introduces uncertainties that models with free parameters can absorb through fitting; a zero-parameter model propagates those uncertainties directly into the prediction, making digitized inputs more hazardous than they would be for a flexible model. A separate MOND analysis (Sanders & Noordermeer 2007) fits the mass-to-light ratio jointly with the rotation curve, introducing circularity that disqualifies those components as independent inputs.

**Namumba Local Group dwarfs.** Namumba et al. (2017, 2018, 2019) publish HI rotation curves and gas profiles for four non-SPARC dwarf galaxies: NGC 6822, Sextans A, Sextans B, and IC 10. The gas side is strong. The blocker is the stellar component: the WISE 3.4 μm radial profiles used in the mass models are not tabulated. For Sextans A and Sextans B, these profiles were provided by T. Jarrett via private communication and were never publicly released. Jarrett has since passed away, and recovery now depends on whether archived copies remain with the original authors or institutions.

**Repetto et al. (2015).** Two non-SPARC galaxies were available with published baryonic components. The schema was nominally complete, but preflight indicated a severe gas-closure mismatch: the gas curve for UGC 4284 implied only about 0.3% of the independently reported HI mass, equivalent to a roughly 300-fold discrepancy. A stellar mass-to-light ratio convention mismatch compounded the problem. The lane was reclassified as a methodology case study; its diagnostic implications are discussed in Section 5.

**GHASP.** The GHASP survey (Epinat et al. 2008) provides Hα rotation curves for over 170 galaxies. Published tables include kinematic parameters but not per-radius baryonic decomposition. A data request for component tables is pending.

These six paths represent distinct failure modes: calibration-set overlap, wrong target type, missing components, unrecoverable or privately held data, convention incompatibility, and unpublished decompositions. The diversity of failure modes is itself a finding. External testing is not blocked by a single missing dataset; it is blocked by the absence of a shared standard for publishing board-complete baryonic decompositions alongside observed kinematics.

This problem is recognized within the field. The BIG-SPARC project (Haubner et al. 2024), currently in development as SPARC's successor, will provide homogeneously derived HI rotation curves, surface brightness profiles, and mass models for approximately 4,000 galaxies. Its motivating rationale explicitly cites the heterogeneous nature of existing rotation curve data and SPARC's limited sample size. BIG-SPARC is not yet publicly available. In correspondence, its project lead did not identify an obvious alternative public database with complete per-radius baryonic decompositions (K. Haubner, personal communication, 2026). This does not criticize BIG-SPARC; it underscores why such a project is needed.

The consequence is that external validation currently proceeds galaxy by galaxy, from small-sample studies whose authors happened to publish sufficient numeric ingredients. That one external board was reconstructed in this paper reflects this scarcity, not a lack of effort.




**Table 1. Candidate-lane failure modes.**

| Lane | Source | Attempted | Board-ready | Failure mode | Status |
|---|---|---:|---:|---|---|
| McGaugh pre-SPARC lane | Pre-SPARC collaboration files (Lelli, McGaugh & Schombert) | 8 | 0 | All eight candidates already present in SPARC rotmod manifest; schema succeeded but independence failed. | Rejected as outside-SPARC lane |
| MaNGA DynPop VII + II | Zhu et al. (2025); Lu et al. (2023) | ∼6,065 CVC / ∼10,296 SP/SFH observations | 0 | JAM-inferred circular velocity is a model-derived potential target, not observed H I/Hα rotation; no independent per-radius Vgas, Vdisk, or Vbulge. | Methodology lane only: JAM potential target |
| Noordermeer ETD | Noordermeer et al. (2007); Noordermeer (2006, 2008); Sanders & Noordermeer (2007) | 19 | 0 | Observed H I/optical target exists, but no machine-readable per-radius Vgas, Vdisk, or Vbulge; thesis mass models are figure/scalar only; MOND-fitted M/L creates circularity risk. | Quarantined / insufficient baryonic closure |
| Namumba Local Group dwarfs | Namumba et al. (2017, 2018, 2019) | 4 | 0 | Good H I rotation curves and gas profiles, but missing WISE 3.4 μm radial stellar profiles or ROTMOD Vdisk arrays; IC 10 has fitted-M/L caution. | Format C; narrow author data recovery only |
| Repetto micro-board | Repetto et al. (2015); Di Teodoro & Fraternali (2014) | 2 | 0 strict; methodology only | Components syntactically exist, but gas mass closure fails; M/L convention mismatch makes native input unsafe. | Methodology case study / invalid as strict baryonic-decomposition test |
| GHASP | Epinat et al. (2008) | >170 | 0 | Hα rotation curves and kinematic parameters are published, but per-radius baryonic decompositions are not. | Data request pending / unpublished decompositions |

# 4. M33 as the First External Board

M33 (NGC 598) was the first candidate that survived the standardization audit as a reconstructable external board. Corbelli et al. (2014) publish a baryonic mass model for M33 including observed HI rotation velocities, uncertainties, gas and stellar surface density profiles, and adopted disk geometry. From these ingredients, a Format B board was reconstructed with documented provenance: mass closure was verified against published integrated masses, the reconstructed rotation components were benchmarked against the source paper's figures, and the stellar disk contribution was computed using a frozen thick-disk kernel. The galaxy-age input was derived from resolved stellar population studies (Williams et al. 2009; Barker et al. 2011) and bracketed at two values: a primary estimate of 6.0 ± 1.2 Gyr and a conservative older estimate of 7.1 ± 1.3 Gyr. M33 was screened against the SPARC manifest and treated as a non-SPARC external board. Unlike SPARC galaxies, M33 did not arrive as a direct rotmod-style board. It became testable only after documented Format B reconstruction, so the result carries both scientific value and provenance limits that a native SPARC board would not.

The board was frozen before model execution. Hermes was run first, and the result was archived with a cryptographic hash before MOND was run on the same frozen inputs.

MOND fits M33 cleanly. On 58 radial points spanning 0.24 to 22.72 kpc, MOND produces χ²ν = 0.12 with a median absolute velocity residual of 3.7 km/s. This is a clean fit under the scoring convention used here. We treat it as a very low-residual result, not a formal likelihood statement; values far below unity may reflect conservative or correlated uncertainties. Because a₀, the interpolation function, the baryonic board, and the mass-to-light convention are all fixed before execution, the low value is not overfitting in the per-galaxy parameter-tuning sense, since there are no per-galaxy parameters to tune.

Hermes produces a non-catastrophic but elevated result: χ²ν = 3.22 at the primary age, χ²ν = 3.39 at the conservative age. The median absolute residual is approximately 32 km/s. Hermes does not catastrophically fail M33, but it does not match MOND's performance on the full board.

This result should be read in context. On the original SPARC calibration set, Hermes and MOND perform comparably overall, with Hermes achieving a median χ²ν of 1.312 and fewer catastrophic failures (84% of galaxies below χ²ν = 5.0, compared to 71% for MOND), while MOND outperforms Hermes on some individual galaxies within the SPARC sample as well (McGinty 2026). M33 is not the first case in which MOND produces a better fit than Hermes on a specific galaxy. It should therefore be read as a portability and boundary-condition case, not as a unique reversal.

The radial structure of the Hermes residual, however, is not uniform. Within the main disk (R ≤ 10 kpc), Hermes produces χ²ν = 0.93 with a mean residual of −7.2 km/s. Beyond 10 kpc, Hermes produces χ²ν = 5.86 with a mean residual of −51.9 km/s. The outer disk contributes approximately 85% of the total Hermes χ². MOND does not reproduce the same large-amplitude outer-disk residual.

Convention sensitivity was tested by scaling the stellar disk contribution across a range of mass-to-light ratio factors. The source-native mass-to-light ratio used by Corbelli et al. (2014) is derived from pixel-level spectral energy distribution modeling, not from a raw near-infrared light profile, so a literal SPARC unit-M/L reconstruction is not possible for M33. A SPARC-like unit-M/L diagnostic analogue reduces the Hermes full-board χ²ν to approximately 2.5, and the best continuous stellar scaling produces χ²ν = 2.45. Convention adjustment narrows the gap but does not eliminate it: the outer-disk residual remains the dominant contributor even at the most favorable scaling, and higher scalings begin to degrade the MOND fit.

The outer-disk stress is spatially correlated with documented features in the M33 literature. The HI disk warps and the position angle twists beyond approximately 8 kpc (Corbelli et al. 2014; Koch et al. 2018). A far-outer disk break with age-gradient reversal occurs near 9 to 12 kpc (Barker et al. 2011). Extended warps, diffuse gas, and a southern cloud or filament structure reach to 22 kpc (Putman et al. 2009). Anomalous low-brightness HI, broad or multiple velocity peaks, and asymmetric motions are documented throughout the outer disk (Kam et al. 2017). A recent analysis interprets the outer-disk misalignment as evidence of gas accretion from a cosmic filament rather than a past M31 passage (Corbelli & Burkert 2024). M33 also sits within the broader M31-M33 interaction environment, with published tidal-interaction models placing possible encounters on Gyr timescales comparable to the age ranges relevant to the M33 board (Bekki 2008; Tepper-García et al. 2019).

This context matters because Hermes is a coherence-decay model. Its current implementation assigns one global age to the entire galaxy, producing one coherence-decay state. That assumes the visible mass distribution belongs to a single historical system. Where M33's inner disk behaves most like a coherent main disk, Hermes fits cleanly. Where the outer disk contains material whose origin, dynamical status, or coherence history may differ from the inner disk, the single-age assumption is the thing under stress. This interpretation is post-hoc: the outer-disk literature was examined after the residual pattern was observed, not before. Confirming whether the correlation between Hermes residuals and documented outer-disk complexity is systematic rather than coincidental requires testing on additional galaxies with varying degrees of outer-disk disturbance.

M33 is therefore not a simple success or a generic failure. It is a boundary case. Hermes ports non-catastrophically and fits the main disk, but loses cleanly to MOND on the full board. The localized outer-disk residual is consistent with the expected behavior of a single-age coherence-decay model applied to a galaxy whose outer regions are dynamically and historically complex. MOND's clean fit remains a real and reported result, but the two models are not equivalently constrained: MOND's interpolation function has no structural ceiling on its low-acceleration boost, while Hermes caps its coherence-decay amplification, making certain outer-disk velocities physically unreachable. That difference affects how their respective failures should be read. The strongest defensible claim is that M33 identifies a boundary of the current Hermes implementation: galaxies with multi-origin or dynamically disturbed outer components may require component-resolved age and dynamical-status treatment before a coherence-decay model can be fairly evaluated across the full radial range.

Because this is one reconstructed external board, it cannot establish broad external success or broad external failure. If the result had gone the other way, the same caution would apply: one galaxy does not validate a framework any more than one galaxy invalidates it. Its value lies in what it shows about the testing process and the structure of one boundary case. The depth of the M33 analysis in this section should not be mistaken for statistical breadth. It reflects the scarcity documented in Section 3: when only one external board can be reconstructed cleanly, that board necessarily carries more interpretive weight than any single galaxy should be asked to bear.




**Table 2. M33 Hermes/MOND result summary for the full board and radial masks.**

| Mask | R range (kpc) | N | Model | χ²ν | χ² fraction (%) | Mean residual (km/s) | Median abs residual (km/s) | RMS residual (km/s) |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| Full board | 0.24–22.72 | 58 | Hermes primary t50=6.0 | 3.22 | 100.00 | -28.02 | 32.00 | 38.99 |
| Full board | 0.24–22.72 | 58 | Hermes conservative t50=7.1 | 3.39 | 100.00 | -29.49 | 33.35 | 39.91 |
| Full board | 0.24–22.72 | 58 | MOND fixed simple | 0.12 | 100.00 | -4.12 | 3.69 | 7.67 |
| Main disk, R ≤ 8 kpc | 0.24–7.57 | 27 | Hermes primary t50=6.0 | 0.78 | 11.24 | -4.18 | 10.55 | 17.46 |
| Main disk, R ≤ 8 kpc | 0.24–7.57 | 27 | Hermes conservative t50=7.1 | 0.78 | 10.75 | -5.71 | 10.37 | 17.49 |
| Main disk, R ≤ 8 kpc | 0.24–7.57 | 27 | MOND fixed simple | 0.02 | 9.31 | 0.15 | 2.58 | 3.08 |
| Warp zone, R > 8 kpc | 8.06–22.72 | 31 | Hermes primary t50=6.0 | 5.35 | 88.76 | -48.78 | 49.61 | 50.79 |
| Warp zone, R > 8 kpc | 8.06–22.72 | 31 | Hermes conservative t50=7.1 | 5.65 | 89.25 | -50.20 | 51.03 | 52.10 |
| Warp zone, R > 8 kpc | 8.06–22.72 | 31 | MOND fixed simple | 0.20 | 90.69 | -7.84 | 7.81 | 10.09 |
| Main disk, R ≤ 10 kpc | 0.24–9.53 | 31 | Hermes primary t50=6.0 | 0.93 | 15.41 | -7.21 | 14.61 | 19.09 |
| Main disk, R ≤ 10 kpc | 0.24–9.53 | 31 | Hermes conservative t50=7.1 | 0.96 | 15.14 | -8.73 | 13.03 | 19.39 |
| Main disk, R ≤ 10 kpc | 0.24–9.53 | 31 | MOND fixed simple | 0.02 | 10.66 | 0.47 | 2.58 | 3.07 |
| Outer tail, R > 10 kpc | 10.02–22.72 | 27 | Hermes primary t50=6.0 | 5.86 | 84.59 | -51.92 | 57.00 | 53.37 |
| Outer tail, R > 10 kpc | 10.02–22.72 | 27 | Hermes conservative t50=7.1 | 6.17 | 84.86 | -53.33 | 58.26 | 54.68 |
| Outer tail, R > 10 kpc | 10.02–22.72 | 27 | MOND fixed simple | 0.23 | 89.34 | -9.39 | 8.55 | 10.75 |
| Barker break outer, R ≥ 9.1 kpc | 9.53–22.72 | 28 | Hermes primary t50=6.0 | 5.71 | 85.50 | -50.99 | 54.71 | 52.63 |
| Barker break outer, R ≥ 9.1 kpc | 9.53–22.72 | 28 | Hermes conservative t50=7.1 | 6.02 | 85.84 | -52.41 | 56.04 | 53.95 |
| Barker break outer, R ≥ 9.1 kpc | 9.53–22.72 | 28 | MOND fixed simple | 0.22 | 89.34 | -9.04 | 8.25 | 10.55 |
| Barker S2 far outer, R ≥ 11.6 kpc | 11.97–22.72 | 23 | Hermes primary t50=6.0 | 6.39 | 78.61 | -55.21 | 57.67 | 56.15 |
| Barker S2 far outer, R ≥ 11.6 kpc | 11.97–22.72 | 23 | Hermes conservative t50=7.1 | 6.70 | 78.51 | -56.54 | 59.05 | 57.43 |
| Barker S2 far outer, R ≥ 11.6 kpc | 11.97–22.72 | 23 | MOND fixed simple | 0.25 | 84.93 | -10.11 | 8.90 | 11.42 |
| Putman tidal-radius exposure, R ≥ 15 kpc | 15.39–22.72 | 16 | Hermes primary t50=6.0 | 7.16 | 61.31 | -60.83 | 59.48 | 61.13 |
| Putman tidal-radius exposure, R ≥ 15 kpc | 15.39–22.72 | 16 | Hermes conservative t50=7.1 | 7.46 | 60.80 | -62.05 | 60.78 | 62.35 |
| Putman tidal-radius exposure, R ≥ 15 kpc | 15.39–22.72 | 16 | MOND fixed simple | 0.32 | 74.83 | -11.96 | 13.22 | 13.05 |

*Note: χ² fraction reports each mask's share of a given model's own total χ². For MOND, high outer fractions occur despite very small absolute residuals because MOND's total χ² is itself small.*


**Table 3. Convention sensitivity diagnostic sweep.** This table shows how both models respond to varying stellar mass-to-light ratio assumptions. It is not a parameter optimization: the source-native result (first row) is the primary reported result. No "best" scaling is adopted. The sweep demonstrates that convention choice affects both models and that no scaling eliminates the Hermes outer-disk residual without degrading the MOND fit.

| Vdisk scale | M/L factor | Global M/L_W1 | t50 (Gyr) | Hermes χ²ν | Hermes outer χ²ν | Hermes outer mean residual (km/s) | MOND χ²ν | MOND outer χ²ν | MOND outer mean residual (km/s) | Note |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1.000 | 1.00 | 0.683 | 6.0 | 3.22 | 5.86 | -51.92 | 0.12 | 0.23 | -9.39 | Baseline/source-native |
| 1.000 | 1.00 | 0.683 | 7.1 | 3.39 | 6.17 | -53.33 | 0.12 | 0.23 | -9.39 | Conservative age bracket |
| 1.325 | 1.76 | 1.199 | 6.0 | 2.45 | 4.17 | -43.48 | 0.41 | 0.09 | 3.43 | Lowest-Hermes diagnostic scaling (not adopted) |
| 1.500 | 2.25 | 1.536 | 6.0 | 2.53 | 3.56 | -40.05 | 0.99 | 0.32 | 10.36 | Discrete sweep |
| 1.500 | 2.25 | 1.536 | 7.1 | 2.52 | 4.06 | -43.01 | 0.99 | 0.32 | 10.36 | Discrete sweep |
| 2.000 | 4.00 | 2.731 | 6.0 | 3.04 | 2.38 | -32.39 | 4.27 | 2.14 | 29.95 | Discrete sweep |
| 2.000 | 4.00 | 2.731 | 7.1 | 2.70 | 2.89 | -36.01 | 4.27 | 2.14 | 29.95 | Discrete sweep |
| 2.500 | 6.25 | 4.267 | 6.0 | 3.61 | 1.43 | -23.60 | 9.92 | 5.57 | 49.01 | Discrete sweep |
| 2.500 | 6.25 | 4.267 | 7.1 | 2.97 | 1.75 | -26.69 | 9.92 | 5.57 | 49.01 | Discrete sweep |
| 3.000 | 9.00 | 6.144 | 6.0 | 4.68 | 0.73 | -11.95 | 17.93 | 10.49 | 67.60 | Discrete sweep |
| 3.000 | 9.00 | 6.144 | 7.1 | 4.03 | 0.82 | -13.86 | 17.93 | 10.49 | 67.60 | Discrete sweep |

# 5. Hermes as a Diagnostic Instrument

A model with no free parameters cannot absorb mismatches or errors in its inputs. When a dark matter halo model encounters a baryonic decomposition with a missing gas component or an incompatible mass-to-light ratio convention, its optimizer adjusts the halo parameters to compensate. The fit may look acceptable, and the input problem remains hidden. A zero-parameter model has no such buffer. If the inputs are wrong, the prediction is wrong, and the residual is visible. This section argues that the resulting failures carry diagnostic information.

**The Repetto case.** Two non-SPARC galaxies (UGC 4284 and UGC 11861) from Repetto et al. (2015) were tested with their published baryonic components. Both Hermes and MOND produced elevated residuals on the native inputs: Hermes χ²ν = 9.57 and 12.73, MOND χ²ν = 3.07 and 5.34. The fact that MOND also showed stress triggered the dual-failure triage described in Section 2.5: when both models fail on the same board, the inputs are suspect before either model is interpreted.

Investigation revealed two compounding problems. The published gas curve for UGC 4284 implied approximately 0.3% of the galaxy's independently reported HI mass (Di Teodoro & Fraternali 2014), a roughly 300-fold discrepancy. The stellar mass-to-light ratio convention used by Repetto et al. was 3 to 5 times lower than the SPARC calibration convention. After applying a declared gas-closure and mass-to-light convention diagnostic adjustment, the board shifted from a MOND sweep to a mixed result: Hermes χ²ν = 2.92 vs MOND 2.99 on UGC 4284, MOND χ²ν = 1.27 vs Hermes 4.32 on UGC 11861. The board was reclassified as a methodology case study rather than a valid model test.

Both problems would be easier to hide in a dark matter halo fit. An optimizer with halo mass, concentration, and core radius as free parameters could compensate for the missing gas and the lighter stellar convention, producing an acceptable fit on unsafe inputs. Hermes could not compensate. Its residual forced the investigation that found the actual problems.

**Pipeline robustness on SPARC-overlap galaxies.** To test whether Hermes failures follow recognizable patterns, we ran the frozen Hermes chassis on SPARC-overlap galaxies using baryonic decompositions from external surveys (THINGS, Walter et al. 2008 and de Blok et al. 2008; LITTLE THINGS, Oh et al. 2015; Iorio et al. 2017; Lelli et al. 2014). These are not external validation tests, because the galaxies are in SPARC. They are pipeline robustness tests: same galaxies, different baryonic inputs.

On THINGS spirals (13 galaxies), Hermes performed better on 7 and MOND on 6, with three catastrophic failures for each model. On dwarf galaxies across three separate external reduction pipelines (12 runs), neither model produced a catastrophic failure. The contrast is informative: spiral galaxies, with their larger disks, bulges, and greater sensitivity to decomposition choices, generated structured Hermes failures that fell into identifiable subclasses. Dwarf galaxies, with simpler baryonic structure, did not. The Hermes chassis remained stable when decomposition pathologies were absent. This pattern is difficult to reconcile with the hypothesis that Hermes is overfit to the SPARC data reduction pipeline: if the equation had memorized SPARC-specific formatting rather than responding to physical inputs, foreign decompositions should produce failures regardless of galaxy type.

The spiral failures sorted into three repeatable subclasses. In one, the foreign decomposition flattened the outer baryonic slope enough to close the Hermes gate prematurely, but the equation could still reach the observed curve if the gate was forced open (NGC 3198, χ²ν = 3.62 on THINGS baryons, 2.84 on SPARC baryons). In another, both the gate and the outer tail were starved simultaneously (NGC 2403). In a third, the outer curve was unreachable regardless of gate state (NGC 2841, χ²ν = 34.3). These are not random scatter. They are classifiable decomposition-dependent failure modes.

**M33 as an observation.** Section 4 showed that the Hermes residual on M33 concentrates beyond 10 kpc, where the literature documents warp, accretion-sensitive gas, and multi-origin structure. Whether this spatial correlation is systematic or coincidental cannot be determined from one galaxy. It is noted here as an observation consistent with the diagnostic framework, not as independent confirmation of it.

**The general point.** A zero-parameter model that fails informatively is not the same as a model that fails randomly. When failures are localizable, classifiable, and correlated with known features of the input data or the target galaxy, the model is functioning as a diagnostic instrument: it flags where assumptions break. This does not make every Hermes failure meaningful, and it does not substitute for a model that fits well. But it suggests that including a rigid model in the testing toolkit provides information that flexible models structurally cannot.


# 6. Discussion and Future Work

This paper set out to test Hermes outside its SPARC calibration set. Three outcomes emerged.

First, the public data ecosystem does not readily support board-complete external testing. Most candidate paths failed before any model could be run, each for a distinct reason: calibration-set overlap, incompatible target definitions, missing baryonic components, unrecoverable or privately held data, and convention mismatches. The BIG-SPARC project is being developed to address this gap, but is not yet available. In the interim, external validation proceeds galaxy by galaxy from scattered small-sample studies.

Second, one external board was reconstructed from published data on M33. On this Format B board, Hermes produced a non-catastrophic but elevated fit, while MOND fit cleanly. The Hermes residual concentrated in M33's outer disk, a region independently documented as dynamically complex and potentially shaped by Local Group interaction history. This is a boundary case for the current single-global-age implementation, not a general external verdict.

Third, the rigidity of a zero-parameter model proved diagnostically useful. In the Repetto methodology case study, Hermes residuals exposed a gas-closure discrepancy and a mass-to-light ratio convention mismatch that would have been easier to hide in a flexible halo fit. On M33, the localized residual structure correlated with documented astrophysical features. These examples suggest that rigid models can flag input and boundary problems that fitted models may absorb.

The primary limitation of this work is sample size. One reconstructed external board cannot establish broad success or broad failure. M33 is a Format B reconstruction, not a native SPARC-style board, and carries provenance limits accordingly.

Two developments are needed for external testing to mature. On the data side, board-complete datasets with standardized per-radius baryonic decompositions, observed kinematics, uncertainties, and documented mass-to-light ratio conventions must become publicly available at scale. BIG-SPARC represents the most promising path toward that goal. On the model side, M33 suggests that galaxies with multi-origin or dynamically disturbed components may require component-resolved age and dynamical-status treatment. The current single-global-age implementation is the first version of the Hermes framework, not the final form. Extending it to handle radially resolved formation histories and equilibrium quality flags is a natural next step, but requires datasets and resources beyond the current project's scope.

This work was carried out by a small independent team without institutional infrastructure. That context does not excuse limitations, but it frames what was achievable: a careful audit of the data landscape, one honestly reported external test, and a diagnostic framework for interpreting rigid-model failures.

The result is therefore not a verdict but a map. Hermes did not receive broad external validation here, and M33 does not overturn MOND. What this work shows is that independent testing requires board-complete data, that such data are scarce, and that a rigid model can expose structured failure modes that flexible models may absorb. The next stage requires larger standardized datasets and a component-resolved treatment of galaxy history. Until those exist, external tests of any baryonic gravity model will remain as much a test of the data infrastructure as of the model itself.


---

# Acknowledgments

We thank Konstantin Haubner for correspondence regarding the BIG-SPARC project and for confirming the current status of publicly available board-complete rotation-curve databases. We thank Brenda Namumba for responding to our data request and for efforts to locate archived WISE stellar profiles. This work was supported by private family funding.

# AI Disclosure

This work was conducted by Team Hermes, a human-AI collaborative research group. AI systems (Claude/Anthropic, ChatGPT/OpenAI, Gemini/Google) were used as research collaborators for data analysis, literature review, methodology development, drafting, and editorial review. The human author retains responsibility for all scientific claims, interpretations, and final approval. The collaborative methodology is documented in McGinty (2026).

# Data Availability

The SPARC database is publicly available (Lelli et al. 2016). The M33 reconstruction board, Hermes execution outputs, MOND comparator outputs, and all diagnostic products are archived with SHA256 provenance locks in the Hermes verification repository (github.com/mcgintyl/hermes-verification) and will be deposited in a Zenodo supplement upon publication. The pipeline robustness test outputs for SPARC-overlap galaxies are archived in the same repository. The full reconstructed M33 board (R, Vobs, errV, Vgas, Vdisk, Vbulge) is printed in Appendix A and provided in machine-readable form as m33_corbelli_reconstructed_board.csv in the same repository.

# Appendix A: M33 Reconstructed Board

The table below is the full per-radius board used for every M33 result in this paper. It was reconstructed from the baryonic mass model of Corbelli et al. (2014): observed HI rotation velocities and uncertainties, the helium-corrected gas component, and the stellar disk component computed with a uniform-thickness annular-quadrature kernel from the source geometry. The stellar mass-to-light convention is the source-native Corbelli (2014) BVIgi pixel-SED stellar surface density, not the SPARC unit mass-to-light normalization (see Section 2.4). M33 was modeled with no bulge, so Vbulge is zero at all radii. Radii are in kpc and all velocities in km/s. The same board is provided in machine-readable form as m33_corbelli_reconstructed_board.csv in the verification repository, with SHA256 provenance locks.

| R (kpc) | Vobs (km/s) | errV (km/s) | Vgas (km/s) | Vdisk (km/s) | Vbulge (km/s) |
|---:|---:|---:|---:|---:|---:|
| 0.24 | 37.3 | 6.2 | 3.04 | 20.97 | 0.0 |
| 0.28 | 37.9 | 5.5 | 3.51 | 24.31 | 0.0 |
| 0.46 | 47.1 | 3.9 | 5.43 | 34.55 | 0.0 |
| 0.64 | 53.5 | 2.8 | 7.06 | 40.96 | 0.0 |
| 0.73 | 55.1 | 3.3 | 7.78 | 43.04 | 0.0 |
| 0.82 | 58.5 | 2.5 | 8.46 | 44.64 | 0.0 |
| 1.08 | 66.2 | 3.3 | 10.12 | 49.13 | 0.0 |
| 1.22 | 69.4 | 2.8 | 10.77 | 51.21 | 0.0 |
| 1.45 | 74.6 | 5.0 | 11.67 | 54.03 | 0.0 |
| 1.71 | 77.9 | 2.9 | 12.80 | 56.16 | 0.0 |
| 1.87 | 81.7 | 2.5 | 13.62 | 56.80 | 0.0 |
| 2.20 | 86.8 | 3.0 | 15.68 | 58.08 | 0.0 |
| 2.28 | 90.1 | 2.1 | 16.25 | 58.40 | 0.0 |
| 2.69 | 94.4 | 3.9 | 18.30 | 59.99 | 0.0 |
| 2.70 | 95.4 | 1.8 | 18.33 | 60.02 | 0.0 |
| 3.12 | 99.2 | 2.6 | 19.54 | 61.15 | 0.0 |
| 3.18 | 98.7 | 2.0 | 19.70 | 61.24 | 0.0 |
| 3.53 | 101.3 | 2.4 | 21.02 | 61.69 | 0.0 |
| 3.66 | 101.5 | 1.3 | 21.69 | 61.75 | 0.0 |
| 4.15 | 106.3 | 2.8 | 24.03 | 61.80 | 0.0 |
| 4.64 | 109.4 | 1.4 | 26.10 | 61.37 | 0.0 |
| 5.13 | 108.8 | 1.6 | 27.89 | 60.60 | 0.0 |
| 5.62 | 107.3 | 2.7 | 29.75 | 59.57 | 0.0 |
| 6.11 | 108.2 | 2.2 | 31.91 | 58.36 | 0.0 |
| 6.60 | 109.8 | 1.6 | 34.20 | 57.04 | 0.0 |
| 7.09 | 110.1 | 1.7 | 36.52 | 55.64 | 0.0 |
| 7.57 | 111.1 | 1.6 | 38.51 | 54.23 | 0.0 |
| 8.06 | 113.0 | 2.6 | 39.51 | 52.77 | 0.0 |
| 8.55 | 113.9 | 2.4 | 39.48 | 51.32 | 0.0 |
| 9.04 | 115.1 | 3.2 | 38.98 | 49.87 | 0.0 |
| 9.53 | 116.3 | 3.1 | 38.02 | 48.42 | 0.0 |
| 10.02 | 119.1 | 2.7 | 36.83 | 46.94 | 0.0 |
| 10.51 | 121.0 | 2.7 | 35.71 | 45.58 | 0.0 |
| 10.99 | 121.5 | 2.1 | 34.75 | 44.46 | 0.0 |
| 11.48 | 118.6 | 2.0 | 33.94 | 43.46 | 0.0 |
| 11.97 | 118.7 | 1.9 | 33.27 | 42.56 | 0.0 |
| 12.46 | 117.2 | 1.9 | 32.59 | 41.75 | 0.0 |
| 12.95 | 116.2 | 2.0 | 31.90 | 41.00 | 0.0 |
| 13.44 | 118.3 | 2.4 | 31.21 | 40.31 | 0.0 |
| 13.93 | 119.0 | 2.4 | 30.55 | 39.66 | 0.0 |
| 14.41 | 121.3 | 2.9 | 30.03 | 39.06 | 0.0 |
| 14.90 | 121.4 | 3.7 | 29.69 | 38.47 | 0.0 |
| 15.39 | 120.3 | 3.4 | 29.33 | 37.93 | 0.0 |
| 15.88 | 121.9 | 4.2 | 28.83 | 37.41 | 0.0 |
| 16.37 | 126.3 | 7.2 | 28.39 | 36.91 | 0.0 |
| 16.86 | 126.3 | 4.5 | 28.08 | 36.43 | 0.0 |
| 17.35 | 127.2 | 5.5 | 27.67 | 35.98 | 0.0 |
| 17.84 | 126.2 | 4.7 | 27.19 | 35.55 | 0.0 |
| 18.32 | 124.2 | 12.2 | 26.67 | 35.13 | 0.0 |
| 18.81 | 127.2 | 13.3 | 26.17 | 34.73 | 0.0 |
| 19.30 | 120.2 | 14.2 | 25.92 | 34.34 | 0.0 |
| 19.79 | 121.8 | 13.2 | 25.90 | 33.98 | 0.0 |
| 20.28 | 136.0 | 19.1 | 25.85 | 33.65 | 0.0 |
| 20.77 | 128.3 | 14.3 | 25.55 | 33.32 | 0.0 |
| 21.26 | 127.4 | 15.4 | 25.25 | 33.01 | 0.0 |
| 21.74 | 120.1 | 13.1 | 25.04 | 32.75 | 0.0 |
| 22.23 | 112.2 | 14.9 | 24.91 | 32.53 | 0.0 |
| 22.72 | 119.6 | 13.4 | 24.85 | 32.36 | 0.0 |

# References

Barker, M. K., et al. 2011, *MNRAS*, 410, 504.

Bekki, K. 2008, *MNRAS*, 390, L24.

Corbelli, E., et al. 2014, *A&A*, 572, A23.

Corbelli, E., & Burkert, A. 2024, *A&A*, 685, A38.

de Blok, W. J. G., et al. 2008, *AJ*, 136, 2648.

Di Teodoro, E. M., & Fraternali, F. 2014, *A&A*, 567, A68.

Epinat, B., et al. 2008, *MNRAS*, 390, 466.

Haubner, K., et al. 2024, "BIG-SPARC: The new SPARC database," IAU Proceedings, arXiv:2411.13329.

Haubner, K. 2026, personal communication.

Iorio, G., et al. 2017, *MNRAS*, 466, 4159.

Kam, S. Z., et al. 2017, *AJ*, 154, 41.

Koch, E. W., et al. 2018, *MNRAS*, 479, 2505.

Lelli, F., McGaugh, S. S., & Schombert, J. M. 2016, *AJ*, 152, 157.

Lelli, F., Verheijen, M. A. W., & Fraternali, F. 2014, *A&A*, 566, A71.

Lu, S., et al. 2023, *MNRAS*, 526, 1022.

McGinty, L. 2026, "The Hermes Gravity Equation," Team Hermes, Zenodo. DOI: 10.5281/zenodo.19904551.

Milgrom, M. 1983, *ApJ*, 270, 365.

Namumba, B., et al. 2017, *MNRAS*, 472, 3761.

Namumba, B., et al. 2018, *MNRAS*, 478, 487.

Namumba, B., et al. 2019, *MNRAS*, 490, 3365.

Noordermeer, E. 2006, PhD thesis, University of Groningen.

Noordermeer, E. 2008, *MNRAS*, 385, 1359.

Noordermeer, E., et al. 2007, *MNRAS*, 376, 1513.

Oh, S.-H., et al. 2015, *AJ*, 149, 180.

Putman, M. E., et al. 2009, *ApJ*, 703, 1486.

Repetto, P., et al. 2015, *MNRAS*, 451, 353.

Sanders, R. H., & Noordermeer, E. 2007, *MNRAS*, 379, 702.

Tepper-García, T., Bland-Hawthorn, J., & Li, D. 2019, "The M31/M33 tidal interaction: a hydrodynamic simulation," arXiv:1911.11911.

Walter, F., et al. 2008, *AJ*, 136, 2563.

Williams, B. F., et al. 2009, *ApJ*, 695, L15.

Zhu, K., et al. 2025, *ApJS*, 280, 55.
