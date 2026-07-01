# Two Layers of Galaxy Aging in MaNGA DynPop: Supplementary Appendix

**Role:** This document provides detailed confound analyses, literature passes, control-variable specifications, and limitation discussions supporting the main paper. The main paper uses a collinearity-free strict-control specification as its primary result (7/8 bins across all three SPS control sets). Some supplementary tables preserve earlier kitchen-sink sensitivity outputs and should be read as sensitivity documentation, not as the primary scoreboard.

**Primary dataset:** MaNGA DynPop DR17 public catalogs  
**Primary analysis date:** 2026-06-26  
**Claim level:** proof of concept / signal documentation, not detection claim

---

## Abstract

See the main paper for the authoritative abstract. This supplementary appendix provides the detailed methodology, literature passes, control-variable specifications, limitation analyses, and sensitivity documentation supporting the main paper's result.

---

## 1. What this note is and is not

This note is a proof-of-concept report. It is not a final paper model, not a hierarchical analysis, and not a claim that MaNGA proves the Hermes framework.

The purpose is narrower. We asked whether MaNGA DynPop contains an age-organized mass-to-light pattern when the data are arranged in the simplest defensible way: fixed stellar mass, youngest quartile versus oldest quartile, medians, bootstrap intervals, and a within-bin shuffle null.

That blunt test passed. Then we asked the harder question: does the result live only in the stellar-population model, or does a smaller dynamical residual remain after the stellar-population layer is stripped away?

The answer is mixed and useful.

The raw DML signal mostly lives in the stellar-population denominator. That is the expected layer. Old stellar populations have higher stellar mass-to-light ratios. No new physics is needed to explain that local mechanism. The standard account is allowed to own the visible layer.

But after explicitly controlling for SPS mass-to-light, metallicity, and standard structural variables, a smaller dynamical residual remains in the same predicted age direction. That is the interesting layer. It is not decisive. It is not large enough to carry the framework by itself. But it is clean enough to document and hand off.

The note therefore makes one empirical claim and one interpretive proposal.

The empirical claim is:

> MaNGA DynPop shows a two-layer age structure in mass-to-light space: a dominant stellar-population layer and a smaller controlled dynamical residual layer. The first is standard and expected. The second remains after the first is explicitly stripped away in this pipeline and motivates independent testing.

The interpretive proposal is:

> The conventional split between “stellar-population effect” and “dynamical residual” is methodologically useful, but it may not be ontological. The two layers may be separate bookkeeping categories, or they may be two age-organized expressions of a broader galactic aging process. This note does not prove the one-box interpretation. It only shows why it is worth testing.

It does not make these claims:

- MaNGA confirms the Hermes equation.
- MaNGA proves age-dependent gravity.
- Standard galaxy formation has no possible explanation.
- The DML result is a clean gravitational observable.
- The raw SPS layer is anomalous.
- The controlled residual survives every possible SPS, IMF, JAM, or simulation test.

The result is smoke, not fire. The cleaner fire test remains weak lensing, because lensing measures the bending of background light rather than a mass-to-light quantity built from stellar dynamics and stellar population modeling.

---

## 2. Why MaNGA is a principle-level test

The locked Hermes rotation-curve equation requires board-complete rotation-curve inputs: observed velocities, observational uncertainties, and independent per-radius gas, stellar disk, and bulge components under known mass-to-light conventions. MaNGA DynPop does not provide that kind of board. Its circular velocities are inferred from Jeans Anisotropic Modelling of stellar kinematics. Those curves are valuable, but they are not observed H I/Hα rotation curves with independent baryonic component arrays.

That distinction matters. A strict equation-level Hermes test cannot be run on MaNGA without changing the target. Treating a JAM-inferred circular-velocity curve as a SPARC-style observed rotation curve would collapse the board-completeness discipline established in the external-data audit.

The MaNGA test in this note is therefore different. It does not ask whether the Hermes equation fits individual MaNGA rotation curves. It asks whether a principle-level prediction appears in a large independent stellar-dynamical population.

The principle is simple:

> If the mass-to-light discrepancy is age-organized, then age should carry information about dynamical mass-to-light structure after stellar mass and standard structural variables are controlled.

In the broader Hermes framework, the corresponding physical interpretation is called the **gravity gap**: the part of a galaxy's dynamical behavior that exceeds what its visible matter would naively predict. In this MaNGA note we mostly avoid that phrase in public-facing summary sections, because the observable here is not a board-complete rotation-curve discrepancy. It is a JAM/SPS mass-to-light construction.

This is closer in logic to the weak-lensing pilot than to the original rotation-curve work. The weak-lensing pilot did not insert the Hermes rotation-curve equation into a lensing pipeline. It extracted the age principle and tested whether the gravitational signal changed with an age proxy at fixed mass. MaNGA is treated the same way here.

---

## 3. Related work and novelty check

This note does not enter empty territory. MaNGA DynPop was built precisely to connect stellar dynamics and stellar populations in a large public sample, and the core ingredients used here are already established products of that program. The novelty is not that age, stellar mass-to-light ratio, dynamical mass-to-light ratio, dark-matter fraction, and circular-velocity-curve shape are related. The novelty, if it survives a full ADS/full-text audit, is the specific decision test used here: a fixed-stellar-mass young-versus-old DML split, followed by an explicit numerator/denominator separation and controlled dynamical-residual hardening.

### 3.1 Closest prior DynPop work

**DynPop I** provides the dynamical side of this note. It presents the quality-assessed JAM catalogue for more than 10,000 MaNGA DR17 galaxies, including masses, mass-to-light ratios, density-profile quantities, and dark-matter components under multiple JAM assumptions. This is the source family for the JAM dynamical numerator used here.

**DynPop II** provides the stellar-population side. It analyzes stellar populations, gradients, and star-formation histories for roughly 10,000 MaNGA galaxies and explicitly links these quantities to dynamical properties. Its headline results already show the expected covariance: stellar age, metallicity, and stellar-population mass-to-light ratio vary systematically with rotation and velocity dispersion. That is the standard layer this note must not pretend to discover.

**DynPop III** studies dynamical scaling relations for about 6,000 galaxies with reliable dynamical models. It reports that total mass-to-light ratio within one effective radius increases with velocity dispersion and with older stellar populations, that total density profiles steepen with galaxy age at fixed velocity dispersion, and that dark-matter fraction within one effective radius is more tightly related to velocity dispersion than to stellar mass. This is close in spirit to the present note because it already connects age, dynamical mass-to-light, density slope, and dark-matter fraction.

**DynPop V** is the closest methodological neighbor. It estimates dark-matter fraction using both JAM decomposition and a comparison between total dynamical mass-to-light ratios and SPS stellar mass-to-light ratios. It also compares JAM and SPS stellar mass-to-light ratios to explore IMF variation, finding that the stellar mass excess factor increases with velocity dispersion, with weak positive correlations with age and no metallicity correlation in the best early-type sample. This paper is the strongest warning against overclaiming our DML result: JAM-SPS offsets are already used in the literature as IMF and dark-matter diagnostics.

**DynPop VII** connects circular-velocity-curve morphology to structure and stellar populations. It reports that declining CVCs dominate among massive, early-type, bulge-dominated, old, metal-rich, early-quenched galaxies, while rising CVCs prevail among younger disk-dominated systems with ongoing star formation. It also identifies bulge-to-total mass ratio, dark-matter fraction within one effective radius, and bulge Sérsic index as the three governing parameters for the unified bulge-disk-halo CVC model. This means age-dependent CVC structure is already in the DynPop literature, but under a conventional bulge-disk-halo interpretation.

### 3.2 The standard visible layer: SPS M/L should age

The dominant raw layer in this note is expected under standard stellar-population modeling. Stellar-population synthesis is built around age-dependent stellar evolution. A simple stellar population is described as a coeval population at a single metallicity, and its spectrum is explicitly a function of age and metallicity. Composite galaxy spectra then fold those simple populations through star-formation history, chemical evolution, dust, and IMF assumptions. In plain language: age is not an afterthought in SPS. It is one of the main knobs that determines how much stellar mass is assigned per unit light.

This matters because the DML quantity subtracts SPS stellar \(M_*/L\):

\[
DML = \log_{10}(M/L)_{dyn} - \log_{10}(M_*/L)_{SPS}.
\]

If old galaxies have higher SPS \(M_*/L\), then old galaxies will automatically be pushed toward smaller DML at fixed dynamical \(M/L\). That is not an anomaly. It is the expected visible layer. Bell and de Jong showed that stellar \(M/L\) varies substantially within and among galaxies and correlates strongly with integrated stellar-population color. Bruzual and Charlot, Vazdekis et al., and Conroy review the same underlying machinery: SPS models use stellar-evolution tracks, stellar libraries, metallicity, IMF, and star-formation history to infer stellar masses and mass-to-light ratios from galaxy light.

The MaNGA-specific literature says the same thing in the language of DynPop. DynPop II reports that stellar age, metallicity, and stellar \(M_*/L\) all decrease with increasing galaxy rotation, while higher-\(\sigma_e\), quenched systems have higher \(M_*/L\) and earlier star-formation histories. Therefore the raw sledgehammer result should be read first as a successful recovery of a known stellar-population covariance: old systems contain older stellar populations, and older stellar populations are dimmer per unit stellar mass.

This is why the hardening test is necessary. The question is not whether the visible SPS layer exists. It does. The question is whether any age-organized dynamical residual remains after that layer, and its obvious covariates, are granted every advantage.

This distinction matters for the framing. A local explanation is not the same thing as a complete partition of the phenomenon. Stellar evolution can explain why old stellar populations have higher SPS \(M_*/L\). That does not by itself decide whether the age structure stops at the stellar-population layer or continues into the dynamical layer. The present note grants the standard SPS explanation first, then asks whether the filing system is complete.

### 3.3 JAM/SPS offsets are already an IMF and dark-matter diagnostic

The second literature point is equally important. Offsets between JAM dynamical \(M/L\) and SPS stellar \(M_*/L\) are not new, and they are not automatically gravitational anomalies. They are already used to study dark-matter fractions, IMF variation, and modeling systematics.

DynPop V is the closest anchor. It estimates dark-matter fraction within \(R_e\) using both direct JAM decomposition and comparison of \((M/L)_{JAM}\) with \((M_*/L)_{SPS}\). In its best early-type subsample, it compares JAM and SPS stellar mass-to-light ratios to infer a stellar mass-excess factor, \(\alpha_{IMF}\). It finds that \(\alpha_{IMF}\) increases with velocity dispersion, detects weak positive correlations with age, and finds no metallicity correlation. It also reports only marginal support from IMF-sensitive spectral features, explicitly warning that one or both methods may be incomplete.

That warning is not isolated. ATLAS\(^{3D}\)-based work found that dynamically inferred IMF normalization is mainly tied to velocity dispersion and only weakly to stellar-population age. Other comparisons between dynamical and spectroscopic IMF estimates show that the two approaches can agree in their broad trend with velocity dispersion while disagreeing galaxy-by-galaxy. The lesson for this note is straightforward: a controlled dynamical residual can be real and still be caused by IMF systematics, dark-matter decomposition, abundance effects, JAM assumptions, aperture differences, or covariance among galaxy structure variables.

Therefore the residual layer should be framed cautiously:

> After the standard SPS age layer is removed, a smaller age-organized dynamical residual remains in this pipeline. Existing JAM/SPS and IMF literature provides several conventional mechanisms that could contribute to such a residual. The result is worth documenting because it survives the first obvious stripping step, not because that survival uniquely identifies new gravity.

The framing difference is not a rejection of that literature. It is a question about whether the literature's categories are complete. The conventional categories separate the stellar-population layer from the dynamical residual because they are measured and modeled by different tools. That separation is useful for analysis. The one-box interpretation asks whether nature itself respects the same separation.

### 3.4 Simulation and assembly-bias comparator status

The strongest standard-model counterargument is assembly history. In \(\Lambda\)CDM, age is not expected to be a meaningless label. Halo formation time, quenching history, morphology, stellar velocity dispersion, dark-matter fraction, and merger history are all coupled. A smaller age-organized dynamical residual could therefore arise without new gravity if old and young galaxies at fixed stellar mass occupy different halo histories, different structural states, or different IMF/systematic regimes.

The targeted literature pass found two important facts.

First, MaNGA-like simulation infrastructure exists. The iMaNGA project generates mock SDSS-IV/MaNGA integral-field spectroscopic observations from IllustrisTNG/TNG50, including MaNGA-like instrumental effects and spectral fitting recovery. iMaNGA is therefore the closest public path toward the correct comparison. However, the iMaNGA papers located in this pass focus on recovered kinematics, stellar ages, metallicities, star-formation histories, and stellar-population gradients, not on rerunning the MaNGA DynPop JAM-plus-SPS DML residual pipeline or reproducing the fixed-mass young/old quartile test used here.

Second, MaNGA and nearby-galaxy dynamics have already been compared to cosmological simulations in adjacent ways. DynPop VI compares MaNGA density slopes to TNG50/TNG100 and reports that the simulations do not reproduce the observed galaxy mass distribution, attributing the mismatch to overestimated dark-matter fraction, constant IMF assumptions, and/or excessive adiabatic contraction. Earlier MaNGA density-slope work compared MaNGA to EAGLE, Illustris, and IllustrisTNG and found that all simulations predicted shallower slopes for massive high-\(\sigma\) galaxies. SEAGLE strong-lensing comparisons likewise show that simulation predictions for central dark-matter fractions depend strongly on the simulation: EAGLE can agree with SLACS, while Illustris and IllustrisTNG are lower than observed, with feedback prescriptions likely driving the differences.

This means the present result is not protected from a standard-model explanation. Simulations already show that inner density slopes, dark-matter fractions, and stellar/dynamical structure are sensitive to feedback, assembly history, and IMF choices. But the pass did not locate the specific comparator that would close the question: a mock MaNGA/DynPop-like sample from EAGLE, TNG, or Illustris, processed through comparable JAM/SPS measurements, split into the same fixed-stellar-mass young/old quartiles, and tested for the same controlled dynamical \(M/L\) residual.

The correct public wording is therefore:

> We did not perform a cosmological mock comparison. Existing MaNGA-like and simulation-comparison literature shows that \(\Lambda\)CDM assembly history, feedback, IMF assumptions, and structural covariance can plausibly generate age-linked dynamical residuals. However, we did not locate a published mock-MaNGA analysis that applies the exact DML component split and controlled young/old residual test reported here. This remains the strongest standard-model hardening task.

This should be treated as a limitation, not as support for the framework. The residual may be a genuine new layer, or it may be an assembly-bias/feedback/IMF/JAM covariance effect that requires a specialist mock pipeline to quantify. This note documents the pattern and states the required comparator; it does not claim that comparator has been passed.

### 3.5 What appears new in this note

The novelty check did not locate a prior paper that reports the exact analysis performed here:

1. define

\[
DML = \log_{10}(M/L)_{dyn,MFL} - \log_{10}(M_*/L)_{SPS}
\]

using the MaNGA DynPop JAM and SP/SFH products;

2. split the `Qual >= 1` DynPop sample into eight equal-count stellar-mass bins;

3. within each bin, compare the youngest 25% and oldest 25% by `T50`, discarding the middle 50%;

4. report a bin-by-bin young-versus-old median scoreboard with bootstrap confidence intervals and within-bin age-shuffle p-values;

5. repeat the split separately on the JAM dynamical numerator and the SPS denominator;

6. test whether a smaller JAM dynamical age residual remains after controlling for SPS mass-to-light, metallicity, stellar mass, size, Sérsic structure, flattening, \(\lambda_R\), and velocity dispersion.

That is the specific contribution. It is an organization of known DynPop quantities into a blunt age-decision test, followed by a hardening step that grants the obvious SPS explanation before asking whether any dynamical age residual remains.

### 3.6 Correct novelty claim

The public claim should be restrained:

> To our knowledge, prior MaNGA DynPop papers have not reported this exact young-versus-old fixed-stellar-mass DML split with a subsequent numerator/denominator hardening test and controlled dynamical-residual split.

The public claim should not be:

> No one has studied age, mass-to-light ratios, dark-matter fraction, or CVC shape in MaNGA before.

They have. That is the point. The two-layer note is useful only if it stands on top of that prior work rather than pretending the prior work does not exist.

### 3.7 Search-status caveat

This novelty check was performed at the level of accessible paper titles, abstracts, searchable full-text snippets, and targeted phrase searches around MaNGA DynPop, JAM/SPS mass-to-light comparison, age, DML, IMF, and young/old splits. It is sufficient to support a cautious "to our knowledge" statement. It is not a formal claim that every ADS citation, conference proceeding, thesis, or unpublished notebook has been exhaustively ruled out.


## 4. Data and pipeline

### 4.1 Files ingested

Three public MaNGA DynPop products were used:

1. `SDSSDR17_MaNGA_JAM_v2.fits`  
   This is the JAM dynamical catalog. It supplies stellar kinematic and dynamical quantities, including \(\lambda_{R_e}\), JAM dynamical mass-to-light measures, total mass estimates, stellar mass estimates, and dark-matter fractions from NFW/gNFW decompositions.

2. `DynPop2_SP_SFH_v2.hdf5.zip`  
   This is the stellar-population and star-formation-history catalog. It supplies age clocks and SPS mass-to-light estimates, including `T50`, `T90`, luminosity-weighted age, mass-weighted age, intrinsic SPS \(M_*/L\), observed SPS \(M_*/L\), metallicities, and quality information.

3. `SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt`  
   This is the circular-velocity-curve table associated with the MaNGA DynPop VII release. It supplies circular velocities at characteristic radii and a `Qual` flag. The table defines `Qual` as JAM model quality and states that only `Qual >= 1` has reliable CVC measurements.

### 4.2 Join validation

The three products joined cleanly.

| Check | Result |
|---|---:|
| JAM rows | 10,296 |
| SP/SFH rows | 10,296 |
| CVC rows | 10,296 |
| JAM ↔ SP/SFH row-order identity match by `plateIFU` and `mangaid` | 1.000 |
| JAM ↔ CVC row-order identity match by `PlateIFU` and `Mangaid` | 1.000 |
| `Qual >= 1` rows | 6,065 |
| Primary finite decision sample | 5,952 |

The primary decision sample required `Qual >= 1`, finite `T50`, finite stellar mass, finite `DML_mfl_int`, and finite `DML_mfl_obs`.

### 4.3 Primary variables

The main age variable is `T50`. Larger `T50` means older assembly. In seven of eight mass bins, the old-quartile median `T50` is pinned at approximately 12.85 Gyr, which is the DynPop star-formation-history grid ceiling. The practical comparison is therefore young galaxies versus age-saturated galaxies rather than young versus moderately old.

The primary mass-to-light discrepancy variables are:

\[
DML_{int} = \log_{10}(M/L)_{dyn,MFL} - \log_{10}(M_*/L)_{SPS,int}
\]

and

\[
DML_{obs} = \log_{10}(M/L)_{dyn,MFL} - \log_{10}(M_*/L)_{SPS,obs}.
\]

In the working table these are called `DML_mfl_int` and `DML_mfl_obs`.

The M/L split test separates those terms into:

- `mfl_cyl_log_ML_dyn`: JAM dynamical \(\log M/L\)
- `sp_ML_int_Re`: intrinsic SPS \(\log M_*/L\)
- `sp_ML_obs_Re`: observed SPS \(\log M_*/L\)

Secondary checks use:

- `gnfw_cyl_fdm_Re`
- `nfw_cyl_fdm_Re`

These dark-matter-fraction variables are not treated as headline outcomes because the earlier age tests showed little independent age signal after standard controls.

---

## 5. Sledgehammer design

The first test deliberately avoids complex modeling.

The sample is split into eight equal-count stellar-mass bins. With 5,952 galaxies, each mass bin contains 744 galaxies.

Inside each stellar-mass bin:

1. sort galaxies by `T50`;
2. define the youngest 25% as the young group;
3. define the oldest 25% as the old group;
4. discard the middle 50%;
5. compare young and old medians.

Each mass bin therefore contains 186 young galaxies and 186 old galaxies.

The EP-favorable direction for DML is:

\[
DML_{young} > DML_{old}.
\]

This follows the framework prediction that younger galaxies should retain a larger mass-to-light discrepancy, while older galaxies should have a smaller one. In framework language, that discrepancy is interpreted as the gravity gap; in this MaNGA note we treat it operationally as DML, not as a direct gravity measurement.

For each outcome and mass bin, the pipeline reports:

- young N;
- old N;
- young median;
- old median;
- young-old median difference;
- bootstrap 95% interval;
- whether the sign is EP-favorable.

The null test shuffles `T50` labels within each stellar-mass bin, rebuilds the young/old quartile split, and recomputes both the number of EP-favorable bins and the summed/mean median difference. The test uses 5,000 shuffles. With the standard +1 correction, the minimum reportable p-value is approximately 0.0002.

---

## 6. Primary sledgehammer result: DML passes

The blunt DML split is numerically strong for the observed-SPS definition and positive overall for intrinsic SPS.

| Outcome | Young > old bins | Mean bin young-old median difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| `DML_mfl_obs` | 8 / 8 | +0.0742 | +0.0669 | [+0.0498, +0.0828] | 0.004 | 0.0002 |
| `DML_mfl_int` | 6 / 8 | +0.0606 | +0.0772 | [+0.0562, +0.0959] | 0.14 | 0.0002 |

The bin-count p reports how often shuffled data produce as many or more EP-favorable bins; the magnitude p reports how often shuffled data produce a summed young-old median difference as large or larger. The observed-SPS result is significant by both metrics. The intrinsic-SPS result is significant by magnitude but not by bin count alone, because two low-mass bins have young below old.

![Kitchen-table DML plot](manga_sledgehammer_final_kitchen_table_DML.png)

![Overall mass-adjusted DML plot](manga_sledgehammer_final_overall_mass_adjusted_DML.png)

The nearest-mass paired check agrees with the main split. For `DML_mfl_int`, young galaxies win 881 of 1,488 pairs, or 59.2%. For `DML_mfl_obs`, young galaxies win 913 of 1,488 pairs, or 61.4%. Both paired shuffle p-values are at the 0.0002 floor.

The secondary fDM checks do not show the same signal.

| Outcome | EP-favorable bins | Mean young-old difference | Joint shuffle p |
|---|---:|---:|---:|
| `gnfw_cyl_fdm_Re` | 5 / 8 | +0.0011 | 0.3007 |
| `nfw_cyl_fdm_Re` | 4 / 8 | -0.0035 | 0.5135 |

This matters. The age signal is strong in the DML mass-to-light discrepancy, not in the parametric NFW/gNFW dark-matter fraction within \(R_e\).

---

## 7. Hardening test: separate the two sides of DML

The DML construction has two sides:

\[
DML = \log(M/L)_{dyn} - \log(M_*/L)_{SPS}.
\]

A positive young-old DML difference can happen because young galaxies have higher dynamical mass-to-light. That would be a gravity-side signal.

It can also happen because old galaxies have higher SPS stellar mass-to-light. Since the SPS term is subtracted, a strong old>SPS effect automatically makes DML higher for young galaxies.

The hardening test therefore repeats the same sledgehammer design on each component separately.

### 7.1 Raw component split

| Outcome | Direction that would enlarge young DML | Bins in that direction | Mean young-old difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Directional shuffle p |
|---|---|---:|---:|---:|---:|---:|
| JAM dynamical \(\log M/L\) | young > old | 0 / 8 | -0.1045 | -0.1075 | [-0.1200, -0.0920] | 1.0000 |
| SPS intrinsic \(\log M_*/L\) | young < old | 8 / 8 | -0.1519 | -0.1693 | [-0.1790, -0.1587] | 0.0002 |
| SPS observed \(\log M_*/L\) | young < old | 8 / 8 | -0.1579 | -0.1569 | [-0.1679, -0.1478] | 0.0002 |

The raw dynamical result goes opposite the simple gravity-side expectation. Old galaxies have higher raw JAM dynamical \(M/L\) than young galaxies in all eight mass bins.

The SPS result is cleaner and larger. Old galaxies have higher stellar-population \(M_*/L\) in all eight bins for both intrinsic and observed SPS definitions. Because that term is subtracted in DML, it drives the raw DML young>old result.

![M/L split kitchen-table plot](manga_ml_split_hardening_kitchen_table.png)

![DML component decomposition plot](manga_ml_split_hardening_DML_decomposition.png)

The first hardening conclusion is therefore unavoidable:

> The raw DML sledgehammer result is mostly carried by the stellar-population denominator. It cannot be treated as clean gravity evidence.

That does not end the test. It defines the next question.

---

## 8. Controlled dynamical residual

The second hardening question is whether dynamical \(M/L\) still carries an age residual after the obvious stellar-population explanation is removed.

The controlled models predict raw JAM dynamical \(M/L\) using standard variables and then ask whether `T50` still contributes.

Controls include:

- SPS \(M_*/L\) (observed, intrinsic, or both, depending on the model);
- mass-weighted and light-weighted metallicity (`sp_MW_Metal_Re`, `sp_LW_Metal_Re`);
- stellar mass (`nsa_sersic_mass`, already log10);
- size (`logRe`, already log10);
- Sérsic index (`nsa_sersic_n`);
- two flattening / axis-ratio terms (`nsa_sersic_ba` and `Eps_MGE`);
- \(\lambda_R\) (`Lambda_Re`);
- velocity dispersion in both raw and log form (`Sigma_Re` and `logSigma_Re`).

No independent stellar surface-density proxy (such as \(\log M - 2\log R_e\)) was included. Earlier drafts listed "surface density" among the controls; that was an error in the description and has been corrected. The double-inclusion of related measures (both metallicity weightings, both flattening terms, both dispersion forms) makes the specification conservative but also introduces collinearity: including a variable and its logarithm in the same linear model inflates variance and can destabilize individual coefficients. This contributes to the specification sensitivity documented in Section 8.3, where alternative reasonable control sets shift a borderline bin.

Because larger `T50` means older, a negative standardized `T50` coefficient means older galaxies have lower dynamical \(M/L\) after controls. Equivalently, young galaxies retain higher controlled dynamical residuals.

### 8.1 Age terms in controlled dynamical M/L models

| Control model | N | Standardized `T50` coefficient | Robust p-value | Model \(R^2\) | Read |
|---|---:|---:|---:|---:|---|
| Controls + intrinsic SPS \(M/L\) | 5,952 | -0.0215 | \(8.1 \times 10^{-12}\) | 0.4799 | older lower after controls |
| Controls + observed SPS \(M/L\) | 5,952 | -0.0346 | \(4.9 \times 10^{-24}\) | 0.5057 | older lower after controls |
| Controls + both SPS \(M/L\) | 5,952 | -0.0357 | \(1.4 \times 10^{-24}\) | 0.5058 | older lower after controls |

The coefficients are small. Their p-values are tiny because the sample is large. The practical effect size is better judged by out-of-sample gain and residual split.

### 8.2 Cross-validated gain from age

| Control model | Controls-only \(R^2\) | Controls + age \(R^2\) | \(\Delta R^2\) from age |
|---|---:|---:|---:|
| Controls + intrinsic SPS \(M/L\) | 0.4701 | 0.4753 | +0.0052 |
| Controls + observed SPS \(M/L\) | 0.4862 | 0.5009 | +0.0146 |
| Controls + both SPS \(M/L\) | 0.4871 | 0.5008 | +0.0137 |

The predictive improvement is modest. It is not a large effect. It is also not zero.

### 8.3 Residual sledgehammer

The same young/old quartile split was then applied to the controlled dynamical residuals.

| Controlled residual outcome | Young > old bins | Mean young-old residual | Mass-adjusted young-old residual | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| After observed SPS controls | 7-8 / 8 | +0.0359 | +0.0358 | [+0.0270, +0.0475] | 0.004-0.035 | 0.0002 |
| After both SPS controls | 7-8 / 8 | +0.0342 | +0.0311 | [+0.0216, +0.0436] | 0.004-0.035 | 0.0002 |
| After intrinsic SPS controls | 6 / 8 | +0.0210 | +0.0188 | [+0.0092, +0.0312] | 0.14 | 0.0002 |

The bin count for the observed-SPS and both-SPS controlled residuals is specification-sensitive: the original control set (which includes doubled metallicity, flattening, and velocity-dispersion terms) produces 8/8; a simplified but reasonable alternative produces 7/8, with the flipped bin sitting near zero. The direction and magnitude of the residual are stable across specifications. The intrinsic-SPS control set is positive in 6/8 bins with a significant magnitude but a nonsignificant bin count. The residual is therefore strongest for the SPS definition that includes dust attenuation (see Section 10.4).

![Controlled dynamical residual plot](manga_ml_split_hardening_control_residuals.png)

This is the central hardening result:

> The raw DML signal is mostly SPS-driven, but a smaller controlled dynamical M/L age residual remains after SPS, metallicity, and standard structural controls.

The residual is smaller than the raw DML signal. It is not a standalone proof. But it is the part of the MaNGA result that remains scientifically interesting after the obvious objection is granted.

As a distribution-free robustness check, Spearman rank correlation between `T50` and the controlled dynamical residual (observed-SPS controls) gives \(\rho = -0.119\), \(p = 2 \times 10^{-20}\). The negative sign confirms the quartile-split direction: older galaxies have lower controlled dynamical residuals. This continuous test does not depend on the quartile split, the binning, or the median statistic.

---

## 9. The two-layer interpretation

*Note: The main paper (Section 4) contains the condensed version of this interpretation. The expanded discussion below provides additional context.*

The cleanest interpretation is not that MaNGA proves age-dependent gravity.

It is that MaNGA shows a two-layer age structure and gives us a reason to test whether the usual two-box filing system is complete.

In conventional analysis, the two effects are naturally filed separately. The stellar-population \(M_*/L\) trend belongs to SPS modeling and stellar evolution. The JAM residual belongs to dynamical modeling, IMF variation, dark-matter fraction, anisotropy, covariance, and structural systematics. That separation is methodologically useful. It is not necessarily ontological.

### Layer 1: the visible stellar-population layer

Old galaxies have old stars. Old stellar populations are dimmer per unit stellar mass. Stellar-population synthesis models therefore assign higher \(M_*/L\) to older systems. The literature pass confirms that this is not merely plausible; it is the standard expectation of SPS modeling.

That layer is standard. This analysis does not dispute it. In the raw sledgehammer split, this layer carries most of the DML result. The analysis should lean into this rather than hide it: the visible layer is the part current astrophysics already explains.

This is where the wrinkle analogy helps, provided it is used carefully. A wrinkle can be explained by collagen loss, UV damage, elastin breakdown, hydration, and cellular repair. Those explanations are real. But explaining a wrinkle does not make it unrelated to aging. It identifies one pathway through which aging becomes visible.

The SPS layer should be treated the same way. Stellar evolution explains why older stellar populations have higher \(M_*/L\). That does not automatically prove that the age structure ends at the SPS layer. It proves only that the visible layer has a known stellar-population mechanism.

### Layer 2: the controlled dynamical residual layer

After the stellar-population layer is explicitly controlled, a smaller dynamical residual remains in the same age direction. Young galaxies retain higher controlled dynamical residuals than old galaxies in nearly every mass bin, depending on the SPS control set.

This layer is the reason the result is worth documenting.

But it has to be stated with the main guardrail intact: raw JAM dynamical \(M/L\) alone did not show the simple young-greater-than-old pattern. In the raw component split, old galaxies had higher JAM dynamical \(M/L\) in all eight mass bins. The raw DML split passed because the SPS \(M_*/L\) age effect was stronger. Only after SPS, metallicity, and structural covariance were controlled did a smaller dynamical residual reappear in the predicted direction.

That is less dramatic than saying “both raw layers point the same way,” but it is more accurate and harder to attack.

### The one-box interpretation

We call the broader reading the one-box interpretation:

> The visible stellar-population layer and the controlled dynamical residual may not be two unrelated age correlations. They may be two age-organized layers of galactic aging measured through different instruments.

This is an interpretation, not a result. The result is the measured two-layer pattern. The interpretation asks whether the standard partition is complete.

The standard explanation may be locally correct but globally incomplete. This note does not fight the standard explanation for Layer 1. It grants it. Then it asks whether Layer 2 should be filed as unrelated noise, IMF/SPS covariance, assembly history, or as a second age-organized layer that deserves joint interpretation.

---

## 10. Limitations

### 10.1 DML is not a clean gravitational observable

`DML` is built from JAM dynamical mass-to-light minus SPS stellar mass-to-light. Both sides carry model assumptions. Neither side is equivalent to direct weak lensing or a board-complete observed rotation curve.

The raw result should not be described as a gravity measurement. It is a mass-to-light discrepancy.

### 10.2 SPS age structure is a major confound

The largest layer of the signal is exactly what standard stellar-population modeling predicts: old stellar populations have higher stellar mass-to-light ratios. That is not a problem for the analysis, but it does limit the claim level.

The literature grounding makes this unavoidable. SPS models are designed to infer stellar mass, age, metallicity, star-formation history, dust effects, and IMF assumptions from galaxy light. A raw DML split that subtracts SPS \(M_*/L\) will inherit that age structure. Any public version must state clearly that the raw DML split is mostly SPS-driven.

This is a concession, not a retreat. The SPS layer is the visible layer. The proof-of-concept begins only after that layer is granted and removed from the dynamical question.

### 10.3 IMF variation: open limitation, unfavorable direction

If the stellar initial mass function varies systematically with galaxy age, velocity dispersion, or morphology, part of the residual could reflect stellar-population modeling mismatch rather than dynamical physics.

This is not hypothetical. DynPop V and prior ATLAS\(^{3D}\) work use JAM/SPS \(M/L\) offsets to study IMF variation, and they find a strong velocity-dispersion connection with weaker age dependence. The present pipeline controls for velocity dispersion and SPS \(M/L\), but it does not settle IMF shape, abundance-ratio, aperture, or spectral-feature systematics.

However, the directional prediction from IMF variation works against the observed residual, not with it. The standard consensus (Conroy & van Dokkum, 2012; ATLAS\(^{3D}\) XX; DynPop V) is that older, high-dispersion galaxies possess bottom-heavy IMFs with excess low-mass stars. A bottom-heavy IMF means that fixed-IMF SPS models underestimate the true stellar mass for old galaxies, which inflates the dynamical-to-stellar mass ratio and pushes DML higher for old systems. Standard unmodeled IMF variation therefore predicts an old-greater-than-young DML residual. The controlled residual measured here shows the opposite. The age signal in this pipeline is surviving against the expected IMF gradient, not riding it.

### 10.4 Dust conventions matter

Intrinsic and observed SPS mass-to-light definitions behave differently. The observed-SPS DML split is visually cleaner, but observed SPS also carries dust and attenuation conventions. The intrinsic-SPS version is less visually perfect but still positive overall.

This correlation deserves direct attention: the strongest controlled residual tracks the SPS definition that carries dust, and young galaxies tend to be dustier. If observed SPS \(M_*/L\) systematically underestimates stellar mass for dusty young galaxies, the controlled dynamical residual would be inflated for the young quartile by dust mismatch rather than dynamical age structure. The intrinsic-SPS residual, which removes the dust pathway, is smaller and achieves only 6/8 bins. The paper retains both definitions so the reader can assess how much of the residual tracks dust rather than dynamics.

### 10.5 JAM covariance remains open

The same stellar kinematic maps contribute to \(\lambda_R\), velocity dispersion, and JAM dynamical quantities. Shared measurement inputs can manufacture or tighten correlations. A final paper would need mock kinematic maps and full covariance treatment.

### 10.6 Structural covariance is strong

Age, mass, size, surface density, morphology, Sérsic structure, metallicity, velocity dispersion, and orbital organization are all entangled in galaxy populations. The controlled residual models reduce this problem but do not eliminate it.

### 10.7 fDM mildly cuts against a gravitational reading

The NFW and gNFW dark-matter fractions within \(R_e\) do not show the same age signal. For a framework that interprets the mass-to-light discrepancy as gravitational, the absence of a corresponding signal in the parametric dark-matter fraction is not neutral. It mildly cuts against the gravitational interpretation, because a gravitational age effect large enough to appear in the controlled DML residual might be expected to leave some trace in the halo decomposition as well. This does not refute the DML residual, but it blocks any simple claim that MaNGA directly shows age-dependent dark-matter fraction, and it should be weighed honestly against the residual when assessing the result.

### 10.8 No cosmological-mock replication yet

The strongest standard-model comparison would use mock MaNGA observations from hydrodynamic simulations, processed through a comparable analysis. MaNGA-like simulation tools and MaNGA-simulation comparisons already exist, including iMaNGA/TNG50 and DynPop density-slope comparisons to TNG50/TNG100, EAGLE, Illustris, and IllustrisTNG. Those works show that simulations can be compared to MaNGA-style stellar-population and dynamical quantities, and they also show that inner mass distributions and dark-matter fractions remain sensitive to feedback, IMF assumptions, and simulation choice.

What has not been done in this note, and what we did not locate in the literature pass, is the exact comparator needed here: a mock MaNGA/DynPop-like sample processed through comparable JAM/SPS observables, then split by age at fixed stellar mass and tested for the same controlled dynamical \(M/L\) residual.

Therefore the residual could still be produced by:

- ordinary \(\Lambda\)CDM assembly bias;
- halo concentration and quenching history;
- feedback-driven structure changes;
- IMF variation tied to velocity dispersion;
- morphology and merger history;
- JAM covariance;
- unmodeled cold gas mass covarying with age;
- mock/observational pipeline differences.

This is the main hardening step that remains outside the present pipeline. Until this comparison is performed, the residual should be described as unexplained by this analysis, not unexplained by \(\Lambda\)CDM.

### 10.9 This is not DESI weak lensing

Weak lensing is the cleaner external test because it measures light deflection. MaNGA is useful because it is public, rich, and already contains ages, dynamics, and stellar-population outputs. But it is entangled with SPS and JAM assumptions in a way lensing is not.


### 10.10 The novelty is narrow

The relevant DynPop literature already shows that stellar population age, stellar-population mass-to-light ratio, dynamical mass-to-light ratio, total density slope, circular-velocity-curve shape, velocity dispersion, dark-matter fraction, and IMF-sensitive JAM-SPS offsets are mutually entangled. The present note should not imply otherwise. Its contribution is the deliberately blunt fixed-mass young/old DML split and the two-layer hardening sequence, not the discovery that MaNGA contains age-dependent dynamical scaling relations.

### 10.11 The one-box interpretation is not a measured fact

The analysis measures a two-layer pattern. It does not measure a common mechanism behind the two layers. The one-box interpretation is a way to organize the question, not an answer to it. A standard-model account could still explain the residual through IMF variation, assembly history, structural covariance, SPS mismatch, dust treatment, cold gas fraction, JAM assumptions, or simulation-predicted age correlations. The paper should invite those tests rather than preempt them.

### 10.12 Unmodeled cold gas mass

The SPS denominator in DML accounts for stellar mass, but the JAM dynamical numerator measures total mass driving the kinematics, including dark matter and cold gas (\(H I\) and \(H_2\)). Young, star-forming galaxies have systematically higher cold gas fractions (e.g., Saintonge et al., 2017; Catinella et al., 2018), so unmodeled gas mass will naturally inflate DML for the young quartile relative to the old. While cold gas mass is generally subdominant to stellar mass within the inner effective radius where these DynPop quantities are evaluated, any systematic covariance between gas fraction and age remains a standard-model contributor to the controlled residual. Future mock comparators should include gas-fraction scaling.

## 11. Hand-off statement

*Note: The main paper (Section 5) contains the authoritative hand-off. The expanded version below provides additional context.*

This analysis reaches the practical limit of what can be responsibly claimed from the present independent pipeline.

We can show that the public MaNGA DynPop catalogs contain a reproducible two-layer age structure. We can show that the raw DML result is mostly carried by SPS mass-to-light. We can show that a smaller controlled dynamical residual remains after that layer is stripped away. We can show that the result points in the direction predicted by the broader age-dependent gravity framework.

We cannot, from this pipeline alone, decide whether the residual is new physics, SPS/IMF mismatch, JAM covariance, structural assembly bias, or some combination of these.

That question requires domain infrastructure this team does not currently have: stellar-population synthesis expertise, IMF-systematics modeling, JAM covariance treatment, cosmological simulation mocks, and independent gravitational measurements.

We therefore present this note as a proof of concept and an invitation. Independent groups with the relevant tools can test, refine, or falsify the result. The goal is not to own a MaNGA anomaly. The goal is to examine whether galaxy aging is being incompletely described when the stellar-population layer and dynamical residual layer are treated as unrelated by default. The standard partition may be correct. It may also be a useful filing system rather than a complete description of the phenomenon.

The framework that motivated this test extends beyond galactic dynamics. This note marks the practical boundary of the present team's contribution to the gravitational evidence, and the broader research program continues in the domains the framework was built to examine.

---

## 12. Data and reproducibility package

The analysis package currently includes:

| File | Role |
|---|---|
| `manga_dynpop_merged_thin_firstpass.csv` | joined working table |
| `manga_firstpass_pipeline_smoke_test.md` | ingestion and join validation |
| `manga_sledgehammer_final_report.md` | primary DML sledgehammer report |
| `manga_sledgehammer_final_scoreboard.csv` | DML scoreboard |
| `manga_sledgehammer_final_bin_results.csv` | per-bin DML table |
| `manga_sledgehammer_final_outputs.zip` | DML sledgehammer outputs |
| `manga_ml_split_hardening_final_report.md` | M/L split and controlled residual report |
| `manga_ml_split_hardening_scoreboard.csv` | raw M/L split scoreboard |
| `manga_ml_split_hardening_control_residual_scoreboard.csv` | controlled residual scoreboard |
| `manga_ml_split_hardening_outputs.zip` | M/L split outputs |

Analysis code and verification outputs are available at the companion GitHub repository (github.com/mcgintyl/hermes-verification). A full reproducibility package including the joined working table, scoreboards, and a manifest distinguishing raw public inputs from derived outputs will be added to the repository with the Zenodo release.

---

## 13. Figure list

**Figure 1. Primary DML sledgehammer kitchen-table plot.**  
Young and old median DML values across eight stellar-mass bins, shown separately for intrinsic and observed SPS definitions.

**Figure 2. Overall mass-adjusted DML young/old comparison.**  
Two-bar summary of the mass-adjusted young and old DML medians with bootstrap intervals.

**Figure 3. M/L component split.**  
The DML numerator and denominator separated into raw JAM dynamical \(M/L\), intrinsic SPS \(M_*/L\), and observed SPS \(M_*/L\).

**Figure 4. DML component decomposition.**  
Shows why the raw DML signal is dominated by the SPS denominator.

**Figure 5. Controlled dynamical M/L residuals.**  
Young-versus-old residual split after SPS, metallicity, and standard structural controls.

---

## 14. Detailed archived tables and sensitivity scoreboards

### Table 1. Primary sledgehammer scoreboard

| Outcome | Young > old bins | Mean young-old median difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| `DML_mfl_obs` | 8 / 8 | +0.0742 | +0.0669 | [+0.0498, +0.0828] | 0.004 | 0.0002 |
| `DML_mfl_int` | 6 / 8 | +0.0606 | +0.0772 | [+0.0562, +0.0959] | 0.14 | 0.0002 |

### Table 2. Raw M/L split scoreboard

| Outcome | Young > old bins | Old > young bins | Direction relevant to DML | Directional bins | Mass-adjusted young-old difference | Bootstrap 95% CI |
|---|---:|---:|---|---:|---:|---:|
| JAM dynamical \(\log M/L\) | 0 / 8 | 8 / 8 | young > old | 0 / 8 | -0.1075 | [-0.1200, -0.0920] |
| SPS intrinsic \(\log M_*/L\) | 0 / 8 | 8 / 8 | old > young | 8 / 8 | -0.1693 | [-0.1790, -0.1587] |
| SPS observed \(\log M_*/L\) | 0 / 8 | 8 / 8 | old > young | 8 / 8 | -0.1569 | [-0.1679, -0.1478] |

### Table 3. Controlled dynamical residual scoreboard

| Controlled residual outcome | Young > old bins | Mean young-old residual | Mass-adjusted young-old residual | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| After observed SPS controls | 7-8 / 8 | +0.0359 | +0.0358 | [+0.0270, +0.0475] | 0.004-0.035 | 0.0002 |
| After both SPS controls | 7-8 / 8 | +0.0342 | +0.0311 | [+0.0216, +0.0436] | 0.004-0.035 | 0.0002 |
| After intrinsic SPS controls | 6 / 8 | +0.0210 | +0.0188 | [+0.0092, +0.0312] | 0.14 | 0.0002 |

Bin counts for the observed-SPS and both-SPS controlled residuals are specification-sensitive. The original control set produces 8/8; a simplified alternative produces 7/8, with one borderline bin near zero. Residual direction and magnitude are stable across specifications.

---

## 15. Conclusion

*Note: The main paper's hand-off (Section 5) is the authoritative conclusion. The expanded version below preserves additional detail from the supplementary analysis.*

MaNGA does not provide a strict board for the Hermes rotation-curve equation. It does provide a large public population in which age, stellar population, internal dynamics, and dynamical mass-to-light can be tested together.

When the data are split bluntly by age at fixed stellar mass, young galaxies show larger DML than old galaxies. The observed-SPS result is the headline: 8/8 mass bins in the predicted direction, with bin-count p of 0.004 and magnitude p at the 5,000-trial floor. The intrinsic-SPS result is positive overall at 6/8 bins but does not reach significance by bin count alone.

The hardening test shows why that result must be handled cautiously. Most of the raw DML age signal is produced by the SPS stellar mass-to-light denominator. That is expected standard astrophysics, not a discovery: older stellar populations have higher stellar \(M_*/L\), and SPS models are built to encode that. Raw JAM dynamical \(M/L\) alone does not show the simple young-greater-than-old pattern.

But the controlled signal does not vanish when that textbook layer is stripped away. After controlling for SPS mass-to-light, metallicity, and standard structural variables, a smaller dynamical mass-to-light residual remains in the same age direction. Existing JAM/SPS, IMF, assembly-history, and simulation literature gives conventional ways this residual could arise. A targeted pass found adjacent MaNGA-like simulation infrastructure, but not the exact mock-MaNGA/JAM/SPS residual replication required here. The result therefore remains a proof-of-concept target for specialist testing rather than a detection claim.

The interpretive contribution is the filing-system question. Standard analysis treats the visible SPS layer and the smaller dynamical residual as separate boxes: stellar evolution over here, dynamical/IMF/systematic residual over there. That separation is useful. It may even be correct. But it is not automatically guaranteed by nature. This note asks whether the two layers should be tested jointly as age-organized expressions of galaxy aging, rather than dismissed separately because one layer has a known local mechanism and the other is small.

That is the finding. Not proof. Not a claim of solved gravity. A two-layer age pattern in a public dynamical catalog, with the obvious stellar-population layer acknowledged and a smaller residual layer documented for independent groups to test.

Whether the filing system that separates those two layers reflects the structure of the phenomenon, or merely the structure of the instruments used to measure it, remains open.

---

## 16. Remaining hardening and future work

Two hardening steps would strengthen the result before a final public version. The standard-model mock comparator remains a future-work item rather than a prerelease requirement.

### 16.1 Quantitative SPS expectation and dust convention check

The medium-depth literature pass establishes the qualitative point: old stellar populations should have higher SPS \(M_*/L\), and JAM/SPS offsets are already used as IMF/dark-matter diagnostics. One quantitative addition would strengthen the note:

- estimate the expected age-driven SPS \(M_*/L\) difference across the actual young-old quartile split using the DynPop SP/SFH values themselves;
- report how much of the raw DML young-old difference is numerically explained by the SPS denominator in each mass bin;
- clarify the dust-convention difference between intrinsic and observed SPS \(M_*/L\);
- preserve the DynPop V warning that JAM/SPS comparison can mix IMF, dark matter, and modeling systematics.

### 16.2 Simulation and assembly-bias comparator status

The literature pass found adjacent but not exact standard-model comparators.

Existing infrastructure:

- iMaNGA provides mock MaNGA-like IFU data cubes from IllustrisTNG/TNG50 and demonstrates recovery of stellar kinematics, ages, metallicities, and star-formation histories.
- DynPop VI compares MaNGA density slopes to TNG50/TNG100 and finds simulation mismatches in the observed mass distribution and dark-matter fraction.
- Earlier MaNGA density-slope work compares MaNGA to EAGLE, Illustris, and IllustrisTNG and finds simulations too shallow for massive high-\(\sigma\) galaxies.
- SEAGLE strong-lensing comparisons show that inner dark-matter fractions differ across EAGLE, Illustris, and IllustrisTNG, with feedback physics likely controlling the mismatch.

What remains missing:

- a mock MaNGA/DynPop-like catalog processed through comparable JAM/SPS mass-to-light measurements;
- the same fixed-stellar-mass young/old quartile split;
- the same DML numerator/denominator separation;
- the same controlled dynamical residual test;
- a reported mock distribution of residual amplitudes against which the observed +0.0188 to +0.0358 controlled young-old residual can be compared.

The defensible summary is therefore:

> A cosmological mock comparison has not been performed. Existing \(\Lambda\)CDM and hydrodynamic simulation literature provides plausible standard routes for an age-linked dynamical residual through assembly history, feedback, IMF variation, morphology, and JAM/SPS covariance. The present note therefore does not claim the residual is anomalous relative to \(\Lambda\)CDM; it identifies the exact mock-MaNGA comparator required to decide that question.

### 16.3 Novelty summary

Based on the completed literature passes:

> We did not locate a prior MaNGA DynPop publication reporting this exact fixed-mass young/old DML split, DML component separation, and controlled dynamical-residual age split. Closely related DynPop work has already established age, mass-to-light, velocity-dispersion, CVC-shape, dark-matter-fraction, and IMF-related covariances. The present contribution is therefore a proof-of-concept organization of existing public quantities and a one-box interpretive question: are the visible SPS layer and the controlled dynamical residual merely separate bookkeeping categories, or two age-organized layers that deserve to be tested together?

## 17. Supplementary references

Zhu, K., Lu, S., Cappellari, M., Li, R., Mao, S., & Gao, L. (2023). *MaNGA DynPop -- I. Quality-assessed stellar dynamical modelling from integral-field spectroscopy of 10K nearby galaxies: a catalogue of masses, mass-to-light ratios, density profiles and dark matter.* arXiv:2304.11711.

Lu, S., Zhu, K., Cappellari, M., Li, R., Mao, S., & Xu, D. (2023). *MaNGA DynPop -- II. Global stellar population, gradients, and star-formation histories from integral-field spectroscopy of 10K galaxies: link with galaxy rotation, shape, and total-density gradients.* arXiv:2304.11712.

Zhu, K., Lu, S., Cappellari, M., Li, R., Mao, S., Gao, L., & Ge, J. (2023). *MaNGA DynPop -- III. Stellar dynamics versus stellar population relations in 6000 early-type and spiral galaxies: Fundamental Plane, mass-to-light ratios, total density slopes, and dark matter fractions.* arXiv:2304.11714.

Lu, S., Zhu, K., Cappellari, M., Li, R., Mao, S., & Xu, D. (2024). *MaNGA DynPop -- V. The dark-matter fraction versus stellar velocity dispersion relation and stellar initial mass function variations in galaxies: dynamical models and full spectrum fitting of integral-field spectroscopy.* arXiv:2309.12395.

Zhu, K., Cappellari, M., Mao, S., Lu, S., Li, R., Shi, Y., Simon, D. A., Fu, Y., & Wang, X. (2025). *MaNGA DynPop -- VII. A Unified Bulge-Disk-Halo Model for Explaining Diversity in Circular Velocity Curves of 6000 Spiral and Early-Type Galaxies.* arXiv:2503.06968.

Bruzual, G., & Charlot, S. (2003). *Stellar population synthesis at the resolution of 2003.* Monthly Notices of the Royal Astronomical Society, 344, 1000-1028. arXiv:astro-ph/0309134.

Bell, E. F., & de Jong, R. S. (2001). *Stellar mass-to-light ratios and the Tully-Fisher relation.* The Astrophysical Journal, 550, 212-229. arXiv:astro-ph/0011493.

Vazdekis, A., Sánchez-Blázquez, P., Falcón-Barroso, J., Cenarro, A. J., Beasley, M. A., Cardiel, N., Gorgas, J., & Peletier, R. F. (2010). *Evolutionary stellar population synthesis with MILES. I. The base models and a new line index system.* Monthly Notices of the Royal Astronomical Society, 404, 1639-1671. arXiv:1004.4439.

Conroy, C. (2013). *Modeling the panchromatic spectral energy distributions of galaxies.* Annual Review of Astronomy and Astrophysics, 51, 393-455. arXiv:1301.7095.

Cappellari, M., et al. (2012). *Systematic variation of the stellar initial mass function in early-type galaxies.* Nature, 484, 485-488.

Cappellari, M., et al. (2013). *The ATLAS3D project — XX. Mass-size and mass-sigma distributions of early-type galaxies: bulge fraction drives kinematics, mass-to-light ratio, molecular gas fraction and stellar initial mass function.* Monthly Notices of the Royal Astronomical Society, 432, 1862-1893.

McDermid, R. M., et al. (2014). *The ATLAS3D project — XXX. Connection between dynamically derived IMF normalisation and stellar population parameters.* The Astrophysical Journal Letters, 792, L37. arXiv:1408.3189.

Smith, R. J. (2014). *Variations in the initial mass function in early-type galaxies: a critical comparison between dynamical and spectroscopic results.* Monthly Notices of the Royal Astronomical Society, 443, L69-L73. arXiv:1403.6114.

Nanni, L., Thomas, D., Trayford, J., Maraston, C., Neumann, J., Law, D. R., Hill, L., Pillepich, A., Yan, R., Chen, Y., & Lazarz, D. (2022). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. I. Construction and analysis of the mock data cubes.* arXiv:2203.11575.

Nanni, L., Thomas, D., Trayford, J., Maraston, C., Neumann, J., Law, D. R., Hill, L., Pillepich, A., Yan, R., Chen, Y., & Lazarz, D. (2022). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. II. The catalogue.* arXiv:2211.13146.

Nanni, L., Neumann, J., Thomas, D., Maraston, C., Trayford, J., Lovell, C. C., Law, D. R., Yan, R., & Chen, Y. (2023). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. III. Stellar metallicity drivers in MaNGA and TNG50.* arXiv:2309.14257.

Li, R., Li, H., Shao, S., Lu, S., Zhu, K., Wang, C., Gao, L., Mao, S., Dutton, A. A., Ge, J., Wang, Y., Leauthaud, A., Zheng, Z., & Bundy, K. (2019). *SDSS-IV MaNGA: The Inner Density Slopes of nearby galaxies.* arXiv:1903.09282.

Li, S., Li, R., Zhu, K., Lu, S., Cappellari, M., Mao, S., Wang, C., & Gao, L. (2023). *MaNGA DynPop -- VI. Matter density slopes from dynamical models of 6000 galaxies versus cosmological simulations: the interplay between baryonic and dark matter.* arXiv:2310.13278.

Mukherjee, S., Koopmans, L. V. E., Tortora, C., Schaller, M., Metcalf, R. B., Schaye, J., & Vernardos, G. (2021). *SEAGLE-III: Towards resolving the mismatch in the dark-matter fraction in early-type galaxies between simulations and observations.* arXiv:2110.07615.

Tinker, J., Wetzel, A., Conroy, C., & Mao, Y.-Y. (2016). *Halo Histories vs. Galaxy Properties at z=0, I: The Quenching of Star Formation.* arXiv:1609.03388.

Lovell, M. R., Pillepich, A., Genel, S., Nelson, D., Springel, V., Pakmor, R., Marinacci, F., Weinberger, R., Torrey, P., Vogelsberger, M., et al. (2018). *The fraction of dark matter within galaxies from the IllustrisTNG simulations.* arXiv:1801.10170.

de Graaff, A., Franx, M., Bell, E. F., Bezanson, R., Schaller, M., Schaye, J., & van der Wel, A. (2022). *A common origin for the Fundamental Plane of quiescent and star-forming galaxies in the EAGLE simulations.* arXiv:2207.13491.

Montero-Dorta, A. D., Chaves-Montero, J., Artale, M. C., & Favole, G. (2021). *On the influence of halo mass accretion history on galaxy properties and assembly bias.* arXiv:2105.05274.

Xu, X., & Zheng, Z. (2018). *Galaxy assembly bias of central galaxies in the Illustris simulation.* arXiv:1812.11210.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.* The Astronomical Journal, 152, 157. arXiv:1606.09251.

Conroy, C., & van Dokkum, P. G. (2012). *The Stellar Initial Mass Function in Early-Type Galaxies from Absorption Line Spectroscopy. II. Results.* The Astrophysical Journal, 760, 71. arXiv:1205.6473.

Saintonge, A., Catinella, B., Tacconi, L. J., Kauffmann, G., Genzel, R., Cortese, L., Davé, R., Fletcher, T. J., Graciá-Carpio, J., Kramer, C., Hainline, L. J., Janowiecki, S., Lutz, D., Rosario, D., Schiminovich, D., Schuster, K., Staguhn, J., Sturm, E., & Wuyts, S. (2017). *xCOLD GASS: The Complete IRAM 30 m Legacy Survey of Molecular Gas for Galaxy Evolution.* The Astrophysical Journal Supplement Series, 233, 22. arXiv:1710.02157.

Catinella, B., Saintonge, A., Janowiecki, S., Cortese, L., Davé, R., Lemonias, J. J., Cooper, A. P., Schiminovich, D., Hainline, L. J., Fabello, S., Geréb, K., Kilborn, V., & Wang, J. (2018). *xGASS: total cold gas scaling relations and molecular-to-atomic gas ratios of galaxies in the local Universe.* Monthly Notices of the Royal Astronomical Society, 476, 875-895. arXiv:1802.02373.
