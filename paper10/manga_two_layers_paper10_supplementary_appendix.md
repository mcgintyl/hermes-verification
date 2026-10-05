# Two Layers of Galaxy Aging in MaNGA DynPop: Supplementary Appendix

**Role:** This document provides detailed confound analyses, literature passes, control-variable specifications, and limitation discussions supporting the main paper. The main paper uses a reduced strict-control specification with galaxy size in kpc as its primary result (Section 8.4); Sections 8.1 to 8.3 document the same specification with angular size, the primary result of earlier versions (7/8 bins across all three SPS control sets). Section 14, Table 3 preserves the earlier kitchen-sink specification as sensitivity documentation, not as the primary scoreboard.

**Primary dataset:** MaNGA DynPop DR17 public catalogs  
**Primary analysis date:** 2026-06-29 (strict-control run; size and distance checks October 2026)  
**Claim level:** proof of concept / signal documentation, not detection claim

---

## 1. What this note is and is not

This note is a proof-of-concept report. It is not a hierarchical analysis and not a claim that MaNGA proves the Hermes framework.

The purpose is narrower. We asked whether MaNGA DynPop contains an age-organized mass-to-light pattern when the data are arranged in the simplest defensible way: fixed stellar mass, youngest quartile versus oldest quartile, medians, bootstrap intervals, and a within-bin shuffle null.

That blunt test passed. Then we asked the harder question: does the result live only in the stellar-population model, or does a smaller dynamical residual remain after the stellar-population layer is stripped away?

The answer is mixed and useful.

The raw DML signal mostly lives in the stellar-population denominator. That is the expected layer. Old stellar populations have higher stellar mass-to-light ratios. No new physics is needed to explain that local mechanism. The standard account is allowed to own the visible layer.

But after explicitly controlling for SPS mass-to-light, metallicity, and standard structural variables, a smaller dynamical residual remains: a linear age term keeps the predicted sign in every full-sample fit, while whether the young-old split shows it, and in how many mass bins, depends on how galaxy size and distance are controlled and on the MaNGA sample (Section 8.4). That is the interesting layer. It is not decisive. It is not large enough to carry the framework by itself. But it is clean enough to document and hand off.

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

The Hermes rotation-curve equation (McGinty, 2026a) requires board-complete rotation-curve inputs: observed velocities, observational uncertainties, and independent per-radius gas, stellar disk, and bulge components under known mass-to-light conventions. MaNGA DynPop does not provide that kind of board. Its circular velocities are inferred from Jeans Anisotropic Modelling of stellar kinematics. Those curves are valuable, but they are not observed H I/Hα rotation curves with independent baryonic component arrays.

That distinction matters. A strict equation-level Hermes test cannot be run on MaNGA without changing the target. Treating a JAM-inferred circular-velocity curve as a SPARC-style observed rotation curve (Lelli, McGaugh, & Schombert, 2016) would collapse the board-completeness discipline established in our Field Audit (McGinty, 2026b).

The MaNGA test in this note is therefore different. It does not ask whether the Hermes equation fits individual MaNGA rotation curves. It asks whether a principle-level prediction appears in a large independent stellar-dynamical population.

The principle is simple:

> If the mass-to-light discrepancy is age-organized, then age should carry information about dynamical mass-to-light structure after stellar mass and standard structural variables are controlled.

In the broader Hermes framework, the corresponding physical interpretation is called the **gravity gap**: the part of a galaxy's dynamical behavior that exceeds what its visible matter would naively predict. The main paper avoids that phrase, because the observable here is not a board-complete rotation-curve discrepancy. It is a JAM/SPS mass-to-light construction.

This is closer in logic to the weak-lensing pilot (McGinty, 2026c) than to the original rotation-curve work (McGinty, 2026a). The weak-lensing pilot did not insert the Hermes rotation-curve equation into a lensing pipeline. It extracted the age principle and tested whether the gravitational signal changed with an age proxy at fixed mass. MaNGA is treated the same way here.

---

## 3. Related work and novelty check

This note does not enter empty territory. MaNGA DynPop was built precisely to connect stellar dynamics and stellar populations in a large public sample, and the core ingredients used here are already established products of that program. The novelty is not that age, stellar mass-to-light ratio, dynamical mass-to-light ratio, dark-matter fraction, and circular-velocity-curve shape are related. The novelty, if it survives a full ADS/full-text audit, is the specific decision test used here: a fixed-stellar-mass young-versus-old DML split, followed by an explicit numerator/denominator separation and controlled dynamical-residual hardening.

### 3.1 Closest prior DynPop work

**DynPop I** provides the dynamical side of this note. It presents the quality-assessed JAM catalog for more than 10,000 MaNGA DR17 galaxies, including masses, mass-to-light ratios, density-profile quantities, and dark-matter components under multiple JAM assumptions. This is the source family for the JAM dynamical numerator used here.

**DynPop II** provides the stellar-population side. It analyzes stellar populations, gradients, and star-formation histories for roughly 10,000 MaNGA galaxies and explicitly links these quantities to dynamical properties. Its headline results already show the expected covariance: stellar age, metallicity, and stellar-population mass-to-light ratio vary systematically with rotation and velocity dispersion. That is the standard layer this note must not pretend to discover.

**DynPop III** studies dynamical scaling relations for about 6,000 galaxies with reliable dynamical models. It reports that total mass-to-light ratio within one effective radius increases with velocity dispersion and with older stellar populations, that total density profiles steepen with galaxy age at fixed velocity dispersion, and that dark-matter fraction within one effective radius is more tightly related to velocity dispersion than to stellar mass. This is close in spirit to the present note because it already connects age, dynamical mass-to-light, density slope, and dark-matter fraction.

**DynPop V** is the closest methodological neighbor. It estimates dark-matter fraction using both JAM decomposition and a comparison between total dynamical mass-to-light ratios and SPS stellar mass-to-light ratios. It also compares JAM and SPS stellar mass-to-light ratios to explore IMF variation, finding that the stellar mass excess factor increases with velocity dispersion, with weak positive correlations with age and no metallicity correlation in the best early-type sample. At fixed \(\sigma_e\) it finds that neither the JAM \(M/L\) nor the JAM dark-matter fraction shows a clear age correlation, while its SPS-based dark-matter fraction, \(1-(M_*/L)_{SPS}/(M/L)_{JAM}\), still varies with age, possibly through gas in the JAM mass. With DynPop V's dust-corrected JAM \(M/L\), that fraction equals \(1-10^{-DML}\) for our observed-SPS DML, so our raw layer restates a published DynPop V result at fixed stellar mass instead of fixed \(\sigma_e\). This paper is the strongest warning against overclaiming our DML result: JAM-SPS offsets are already used in the literature as IMF and dark-matter diagnostics.

**DynPop VII** connects circular-velocity-curve morphology to structure and stellar populations. It reports that declining CVCs dominate among massive, early-type, bulge-dominated, old, metal-rich, early-quenched galaxies, while rising CVCs prevail among younger disk-dominated systems with ongoing star formation. It also identifies bulge-to-total mass ratio, dark-matter fraction within one effective radius, and bulge Sérsic index as the three governing parameters for the unified bulge-disk-halo CVC model. This means age-dependent CVC structure is already in the DynPop literature, but under a conventional bulge-disk-halo interpretation.

### 3.2 The standard visible layer: SPS M/L should age

The dominant raw layer in this note is expected under standard stellar-population modeling. Stellar-population synthesis is built around age-dependent stellar evolution. A simple stellar population is described as a coeval population at a single metallicity, and its spectrum is explicitly a function of age and metallicity. Composite galaxy spectra then fold those simple populations through star-formation history, chemical evolution, dust, and IMF assumptions. In plain language: age is not an afterthought in SPS. It is one of the main knobs that determines how much stellar mass is assigned per unit light.

This matters because the DML quantity subtracts SPS stellar \(M_*/L\):

\[
DML = \log_{10}(M/L)_{dyn} - \log_{10}(M_*/L)_{SPS}.
\]

If old galaxies have higher SPS \(M_*/L\), then old galaxies will automatically be pushed toward smaller DML at fixed dynamical \(M/L\). That is not an anomaly. It is the expected visible layer. Bell and de Jong showed that stellar \(M/L\) varies substantially within and among galaxies and correlates strongly with integrated stellar-population color. Bruzual and Charlot and Vazdekis et al. build, and Conroy reviews, the same underlying machinery: SPS models use stellar-evolution tracks, stellar libraries, metallicity, IMF, and star-formation history to infer stellar masses and mass-to-light ratios from galaxy light.

The MaNGA-specific literature says the same thing in the language of DynPop. DynPop II reports that stellar age, metallicity, and stellar \(M_*/L\) all decrease with increasing galaxy rotation, while higher-\(\sigma_e\), quenched systems have higher \(M_*/L\) and earlier star-formation histories. Therefore the raw fixed-mass young-versus-old DML result (the "sledgehammer" split of Sections 5 and 6) should be read first as a successful recovery of a known stellar-population covariance: old systems contain older stellar populations, and older stellar populations are dimmer per unit stellar mass.

This is why the hardening test is necessary. The question is not whether the visible SPS layer exists. It does. The question is whether any age-organized dynamical residual remains after that layer, and its obvious covariates, are granted every advantage.

### 3.3 JAM/SPS offsets are already an IMF and dark-matter diagnostic

The second literature point is equally important. Offsets between JAM dynamical \(M/L\) and SPS stellar \(M_*/L\) are not new, and they are not automatically gravitational anomalies. They are already used to study dark-matter fractions, IMF variation, and modeling systematics.

DynPop V is the closest anchor. It estimates dark-matter fraction within \(R_e\) using both direct JAM decomposition and comparison of \((M/L)_{JAM}\) with \((M_*/L)_{SPS}\). In its best early-type subsample, it compares JAM and SPS stellar mass-to-light ratios to infer a stellar mass-excess factor, \(\alpha_{IMF}\). It finds that \(\alpha_{IMF}\) increases with velocity dispersion, detects weak positive correlations with age, and finds no metallicity correlation. It also reports only marginal support from IMF-sensitive spectral features, noting that this casts doubt on the validity of one or both methods to measure the IMF.

That warning is not isolated. ATLAS\(^{3D}\) work found that dynamically inferred IMF normalization varies systematically among early-type galaxies (Cappellari et al., 2012), mainly with velocity dispersion (Cappellari et al., 2013), and only weakly with stellar-population age (McDermid et al., 2014). A direct comparison of dynamical and spectroscopic IMF estimates for the same 34 galaxies found that both favor an IMF heavier than the Milky Way's on average and both show a raw trend with velocity dispersion, yet the two estimates do not correlate galaxy-by-galaxy, and once velocity dispersion and [Mg/Fe] are fitted together the dynamical estimate tracks velocity dispersion while the spectroscopic one tracks [Mg/Fe] (Smith, 2014). The lesson for this note is straightforward: a controlled dynamical residual can be real and still be caused by IMF systematics, dark-matter decomposition, abundance effects, JAM assumptions, aperture differences, or covariance among galaxy structure variables.

Therefore the residual layer should be framed cautiously:

> After the standard SPS age layer is removed, a smaller age-organized dynamical residual remains in this pipeline. Existing JAM/SPS and IMF literature provides several conventional mechanisms that could contribute to such a residual. The result is worth documenting because it survives the first obvious stripping step, not because that survival uniquely identifies new gravity.

### 3.4 Simulation and assembly-bias comparator status

The strongest standard-model counterargument is assembly history. In \(\Lambda\)CDM, age is not expected to be a meaningless label. Halo formation time, quenching history, morphology, stellar velocity dispersion, dark-matter fraction, and merger history are all coupled. In hydrodynamic simulations, star-formation history and quenching track halo formation time or mass-accretion history at fixed halo mass (Xu & Zheng, 2020; Montero-Dorta et al., 2021). In SDSS central galaxies, the measured link between halo formation history and quenching is small but statistically significant above a stellar mass of about \(10^{10}\,h^{-2}\,M_\odot\) and absent below it (Tinker et al., 2017). A smaller age-organized dynamical residual could therefore arise without new gravity if old and young galaxies at fixed stellar mass occupy different halo histories, different structural states, or different IMF/systematic regimes.

The targeted literature pass found two important facts.

First, MaNGA-like simulation infrastructure exists. The iMaNGA project (Nanni et al., 2022, 2023, 2024) generates mock SDSS-IV/MaNGA integral-field spectroscopic observations from IllustrisTNG/TNG50, including MaNGA-like instrumental effects and spectral fitting recovery. iMaNGA is therefore the closest public path toward the correct comparison. However, the iMaNGA papers located in this pass focus on recovered kinematics, stellar ages, metallicities, star-formation histories, and stellar-population gradients, not on rerunning the MaNGA DynPop JAM-plus-SPS DML residual pipeline or reproducing the fixed-mass young/old quartile test used here.

Second, MaNGA and nearby-galaxy dynamics have already been compared to cosmological simulations in adjacent ways. DynPop VI compares MaNGA density slopes to TNG50/TNG100 and reports that the simulations do not reproduce the observed galaxy mass distribution, attributing the mismatch to an overestimated dark-matter fraction in the simulations, possibly due to a constant IMF and excessive adiabatic contraction. Earlier MaNGA density-slope work (Li et al., 2019) compared MaNGA to EAGLE, Illustris, and IllustrisTNG and found that all simulations predicted shallower slopes for massive high-\(\sigma\) galaxies. SEAGLE strong-lensing comparisons likewise show that simulation predictions for central dark-matter fractions depend strongly on the simulation: EAGLE can agree with SLACS, while Illustris and IllustrisTNG are lower than observed, with feedback prescriptions likely driving the differences.

This means the present result is not protected from a standard-model explanation. Simulations already show that inner density slopes, dark-matter fractions, and stellar/dynamical structure are sensitive to feedback, assembly history, and IMF choices (Lovell et al., 2018; Montero-Dorta et al., 2021; Mukherjee et al., 2022; de Graaff et al., 2023). But the pass did not locate the specific comparator that would close the question: a mock MaNGA/DynPop-like sample from EAGLE, TNG, or Illustris, processed through comparable JAM/SPS measurements, split into the same fixed-stellar-mass young/old quartiles, and tested for the same controlled dynamical \(M/L\) residual.

We therefore state the comparator status as follows:

> We did not perform a cosmological mock comparison. Existing MaNGA-like and simulation-comparison literature shows that \(\Lambda\)CDM assembly history, feedback, IMF assumptions, and structural covariance can plausibly generate age-linked dynamical residuals. However, we did not locate a published mock-MaNGA analysis that applies the exact DML component split and controlled young/old residual test reported here. This remains the strongest standard-model hardening task.

This should be treated as a limitation, not as support for the framework. The residual may be a genuine new layer, or it may be an assembly-bias/feedback/IMF/JAM covariance effect that requires a specialist mock pipeline to quantify. This note documents the pattern and states the required comparator; it does not claim that comparator has been passed.

### 3.5 What appears new in this note

The novelty check did not locate a prior paper that reports the exact analysis performed here:

1. define

\[
DML = \log_{10}(M/L)_{dyn,MFL} - \log_{10}(M_*/L)_{SPS}
\]

using the MaNGA DynPop JAM and SP/SFH products, where MFL denotes the DynPop I mass-follows-light JAM model (a constant total \(r\)-band \(M/L\)) fitted with a cylindrically aligned velocity ellipsoid (JAM\(_{cyl}\));

2. split the `Qual >= 1` DynPop sample into eight equal-count stellar-mass bins;

3. within each bin, compare the youngest 25% and oldest 25% by `T50`, discarding the middle 50%;

4. report a bin-by-bin young-versus-old median scoreboard with bootstrap confidence intervals and within-bin age-shuffle p-values;

5. repeat the split separately on the JAM dynamical numerator and the SPS denominator;

6. test whether a smaller JAM dynamical age residual remains after controlling for SPS mass-to-light, metallicity, stellar mass, galaxy size (in kpc for the main result; angular in Sections 8.1 to 8.3), Sérsic structure, flattening, \(\lambda_R\), and velocity dispersion.

That is the specific contribution. It is an organization of known DynPop quantities into a blunt age-decision test, followed by a hardening step that grants the obvious SPS explanation before asking whether any dynamical age residual remains.

### 3.6 Correct novelty claim

Our novelty claim is deliberately restrained:

> To our knowledge, prior MaNGA DynPop papers have not reported this exact young-versus-old fixed-stellar-mass DML split with a subsequent numerator/denominator hardening test and controlled dynamical-residual split.

We do not claim that age, mass-to-light ratios, dark-matter fraction, or CVC shape are unstudied in MaNGA. They have been studied (Sections 3.1 to 3.3), and the two-layer result is useful only if it stands on top of that prior work.

### 3.7 Search-status caveat

This novelty check was performed at the level of accessible paper titles, abstracts, searchable full-text snippets, and targeted phrase searches around MaNGA DynPop, JAM/SPS mass-to-light comparison, age, DML, IMF, and young/old splits. It is sufficient to support a cautious "to our knowledge" statement. It is not a formal claim that every ADS citation, conference proceeding, thesis, or unpublished notebook has been exhaustively ruled out.


## 4. Data and pipeline

### 4.1 Files ingested

Three public MaNGA DynPop products were used:

1. `SDSSDR17_MaNGA_JAM_v2.fits`  
   This is the JAM dynamical catalog. It supplies stellar kinematic and dynamical quantities, including \(\lambda_{R_e}\), JAM dynamical mass-to-light measures, total mass estimates, stellar mass estimates, and dark-matter fractions from NFW/gNFW decompositions.

2. `DynPop2_SP_SFH_v2.hdf5.zip`  
   This is the stellar-population and star-formation-history catalog. It supplies age clocks and SPS mass-to-light estimates, including `T50`, `T90`, luminosity-weighted age, mass-weighted age, intrinsic SPS \(M_*/L\), observed SPS \(M_*/L\), metallicities, and fit diagnostics (`SNR_Re`, `chi2_Re`); the quality flag `Qual` comes from the JAM and CVC tables, which agree on every row.

3. `SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt`  
   This is the circular-velocity-curve table associated with the MaNGA DynPop VII release. It supplies circular velocities at characteristic radii and a `Qual` flag (identical, galaxy by galaxy, to the JAM catalog's). No circular-velocity quantity enters the results reported here. The table defines `Qual` as JAM model quality and states that only `Qual >= 1` has reliable CVC measurements.

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

The primary decision sample required `Qual >= 1`, finite `T50`, finite stellar mass, finite `DML_mfl_int`, and finite `DML_mfl_obs`. The 10,296 rows are PlateIFU entries of 10,160 distinct MaNGA IDs, and the 5,952 entries of the decision sample belong to 5,875 distinct MaNGA IDs, since some galaxies were observed more than once. Keeping one entry per MaNGA ID leaves every bin count in main-paper Table 1 unchanged.

### 4.3 Primary variables

The main age variable is `T50`. Larger `T50` means older assembly. It is the lookback time at which a galaxy reached half of its present-day stellar mass (DynPop II; Lu et al., 2023, Table B1), reported in steps of about 0.16 Gyr. In seven of eight mass bins the old-quartile median is 12.85 Gyr, the most common value in the sample (278 of 5,952 galaxies); in the other it is 13.01 Gyr. That is the mode, not a grid limit: 605 galaxies lie above it, up to 14.13 Gyr, and the DynPop II template ages reach 15.85 Gyr. The old quartiles are nonetheless compressed: every old-quartile value is at least 12.0 Gyr, with interquartile ranges of 0.3 to 0.6 Gyr against 1.8 to 3.2 Gyr for the young quartiles, and DynPop II notes that age differences are difficult to measure above about 8 Gyr. The practical comparison is therefore young galaxies versus a narrow band of the oldest galaxies rather than young versus moderately old.

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

- `mfl_cyl_log_ML_dyn`: JAM dynamical \(\log M/L\) from the mass-follows-light model (column prefixes `mfl_cyl_`, `nfw_cyl_` and `gnfw_cyl_` name the DynPop I mass model, each fitted with the cylindrically aligned JAM\(_{cyl}\) velocity ellipsoid)
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

The predicted direction for DML is:

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
- whether the sign is in the predicted direction.

The null test shuffles `T50` labels within each stellar-mass bin, rebuilds the young/old quartile split, and recomputes both the number of bins in the predicted direction and the summed/mean median difference. The test uses 5,000 shuffles. With the standard +1 correction, the minimum reportable p-value is approximately 0.0002. Because `T50` is discrete, 15 of the 16 quartile cuts fall inside a group of tied ages; ties are broken by stellar mass (lower mass sorts as younger).

The mass-adjusted difference subtracts from every young-quartile and old-quartile galaxy the median of all 744 galaxies in its mass bin, pools the eight bins, and takes the young median minus the old median; its 95% interval comes from bootstrap resampling of the pooled young and old groups.

---

## 6. Raw DML sledgehammer result

The blunt DML split is numerically strong for the observed-SPS definition and positive overall for intrinsic SPS.

| Outcome | Young > old bins | Mean bin young-old median difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| `DML_mfl_obs` | 8 / 8 | +0.0742 | +0.0669 | [+0.0498, +0.0828] | 0.006 | 0.0002 |
| `DML_mfl_int` | 6 / 8 | +0.0606 | +0.0772 | [+0.0562, +0.0959] | 0.15 | 0.0002 |

The bin-count p reports how often shuffled data produce as many or more young > old bins; the magnitude p reports how often shuffled data produce a summed young-old median difference as large or larger. The observed-SPS result is significant by both metrics. The intrinsic-SPS result is significant by magnitude but not by bin count alone, because two low-mass bins have young below old.

The nearest-mass paired check agrees with the main split. For `DML_mfl_int`, young galaxies win 881 of 1,488 pairs, or 59.2%. For `DML_mfl_obs`, young galaxies win 913 of 1,488 pairs, or 61.4%. Both paired shuffle p-values are at the 0.0002 floor.

The secondary fDM checks do not show the same signal.

| Outcome | Young > old bins | Mean young-old difference | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|
| `gnfw_cyl_fdm_Re` | 5 / 8 | +0.0011 | 0.35 | 0.48 |
| `nfw_cyl_fdm_Re` | 4 / 8 | -0.0035 | 0.62 | 0.65 |

This matters. The age signal is strong in the DML mass-to-light discrepancy, not in the parametric NFW/gNFW dark-matter fraction within \(R_e\).

---

## 7. Hardening test: separate the two sides of DML

The DML construction has two sides:

\[
DML = \log(M/L)_{dyn} - \log(M_*/L)_{SPS}.
\]

A positive young-old DML difference can happen because young galaxies have higher dynamical mass-to-light. That would be a gravity-side signal.

It can also happen because old galaxies have higher SPS stellar mass-to-light. Since the SPS term is subtracted, a strong old > young SPS effect automatically makes DML higher for young galaxies.

The hardening test therefore repeats the same sledgehammer design on each component separately.

### 7.1 Raw component split

| Outcome | Direction that would enlarge young DML | Bins in that direction | Mean young-old difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Directional magnitude p |
|---|---|---:|---:|---:|---:|---:|
| JAM dynamical \(\log M/L\) | young > old | 0 / 8 | -0.1045 | -0.1075 | [-0.1200, -0.0920] | 1.0000 |
| SPS intrinsic \(\log M_*/L\) | young < old | 8 / 8 | -0.1519 | -0.1693 | [-0.1790, -0.1587] | 0.0002 |
| SPS observed \(\log M_*/L\) | young < old | 8 / 8 | -0.1579 | -0.1569 | [-0.1679, -0.1478] | 0.0002 |

The raw dynamical result goes opposite the simple gravity-side expectation. Old galaxies have higher raw JAM dynamical \(M/L\) than young galaxies in all eight mass bins.

The SPS result is cleaner and larger. Old galaxies have higher stellar-population \(M_*/L\) in all eight bins for both intrinsic and observed SPS definitions. Because that term is subtracted in DML, it drives the raw DML young>old result.

The first hardening conclusion is therefore unavoidable:

> The raw DML sledgehammer result is mostly carried by the stellar-population denominator. It cannot be treated as clean gravity evidence.

That does not end the test. It defines the next question.

---

## 8. Controlled dynamical residual

The second hardening question is whether dynamical \(M/L\) still carries an age residual after the obvious stellar-population explanation is removed.

The controlled models predict raw JAM dynamical \(M/L\) using standard variables and then ask whether `T50` still contributes. Sections 8.1 to 8.3 report the strict specification with angular size, the primary result of earlier versions of this paper; Section 8.4 repeats it with galaxy size in kpc, which the main paper now uses as its primary result, and adds distance and MaNGA-subsample checks. Section 14, Table 3 keeps the earlier kitchen-sink specification as sensitivity documentation.

The strict specification has these controls:

- SPS \(M_*/L\) (observed, intrinsic, or both, depending on the model);
- mass-weighted metallicity (`sp_MW_Metal_Re`);
- stellar mass (`nsa_sersic_mass`, already log10);
- angular size (`logRe`, log10 of the MGE half-light radius `Re_arcsec_MGE` in arcsec, not converted to kpc);
- Sérsic index (`nsa_sersic_n`);
- MGE ellipticity (`Eps_MGE`);
- \(\lambda_R\) (`Lambda_Re`);
- log velocity dispersion (`logSigma_Re`).

The earlier kitchen-sink specification also included light-weighted metallicity (`sp_LW_Metal_Re`), a second flattening term (`nsa_sersic_ba`), and raw velocity dispersion (`Sigma_Re`) alongside its logarithm. Including a variable and its logarithm, or two measures of one property, inflates the variance of individual coefficients. The strict specification drops those duplicates but keeps the ordinary correlations among galaxy properties, and its both-SPS model deliberately keeps the intrinsic and observed SPS terms together. No separate stellar surface-density term was included; with \(R_e\) in arcsec, \(\log M - 2\log R_e\) is a combination of the mass and size controls, but a physical surface density would also need \(\log D_A\) (Section 8.4).

Because larger `T50` means older, a negative `T50` coefficient means older galaxies have lower dynamical \(M/L\) after controls. Equivalently, young galaxies retain higher controlled dynamical residuals.

### 8.1 Age terms in controlled dynamical M/L models

| Control model | N | `T50` coefficient (dex per 1 SD of `T50`) | Robust (HC3) p-value | Model \(R^2\) | Read |
|---|---:|---:|---:|---:|---|
| Controls + intrinsic SPS \(M/L\) | 5,952 | -0.0267 | \(3.6 \times 10^{-16}\) | 0.4258 | older lower after controls |
| Controls + observed SPS \(M/L\) | 5,952 | -0.0404 | \(6.0 \times 10^{-32}\) | 0.4560 | older lower after controls |
| Controls + both SPS \(M/L\) | 5,952 | -0.0416 | \(4.7 \times 10^{-31}\) | 0.4561 | older lower after controls |

The coefficients are small. Their p-values are tiny because the sample is large. The practical effect size is better judged by out-of-sample gain and residual split.

### 8.2 Cross-validated gain from age

| Control model | Controls-only \(R^2\) | Controls + age \(R^2\) | \(\Delta R^2\) from age |
|---|---:|---:|---:|
| Controls + intrinsic SPS \(M/L\) | 0.4126 | 0.4212 | +0.0086 |
| Controls + observed SPS \(M/L\) | 0.4303 | 0.4513 | +0.0209 |
| Controls + both SPS \(M/L\) | 0.4320 | 0.4512 | +0.0192 |

The predictive improvement is modest. It is not a large effect. It is also not zero.

### 8.3 Residual sledgehammer

The same young/old quartile split was then applied to the controlled dynamical residuals.

| Controlled residual outcome | Young > old bins | Mean young-old residual | Mass-adjusted young-old residual | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| After observed SPS controls | 7 / 8 | +0.0393 | +0.0387 | [+0.0292, +0.0494] | 0.031 | 0.0002 |
| After both SPS controls | 7 / 8 | +0.0347 | +0.0344 | [+0.0231, +0.0442] | 0.034 | 0.0002 |
| After intrinsic SPS controls | 7 / 8 | +0.0238 | +0.0225 | [+0.0114, +0.0342] | 0.037 | 0.0002 |

All three SPS control sets give 7/8 bins. The one negative bin is the same in all three (the sixth of the eight mass bins, median NSA \(\log M \approx 10.66\)), and its bootstrap interval crosses zero in each. All 12 variants of the strict set (mass- or light-weighted metallicity, MGE ellipticity or NSA Sérsic axis ratio, three SPS sets) also give 7/8; they do not vary the size term (Section 8.4). The kitchen-sink set (Section 14, Table 3) gives 8/8, 8/8, and 6/8 for the observed, both, and intrinsic sets, so between these two specifications the bin count is specification-sensitive while the direction and the magnitude test are not; size and distance controls change those as well (Section 8.4). The residual is largest for the observed-SPS definition, which, like the JAM \(M/L\), is defined per unit of observed (dust-attenuated) light (Section 10.4).

With angular size, the hardening result is:

> The raw DML signal is mostly SPS-driven, but a smaller controlled dynamical M/L age residual remains after SPS, metallicity, and standard structural controls.

The residual is smaller than the raw DML signal. It is not a standalone proof. But it is the part of the MaNGA result that remains scientifically interesting after the obvious objection is granted.

As a distribution-free robustness check, Spearman rank correlation between `T50` and the controlled dynamical residual gives \(\rho = -0.117\) (\(p = 1.7 \times 10^{-19}\)) with observed-SPS controls, \(-0.103\) (\(p = 1.4 \times 10^{-15}\)) with both, and \(-0.065\) (\(p = 5.4 \times 10^{-7}\)) with intrinsic SPS; the kitchen-sink set gives \(-0.119\), \(-0.110\), and \(-0.067\). The negative sign confirms the quartile-split direction: older galaxies have lower controlled dynamical residuals. This continuous test does not depend on the quartile split, the binning, or the median statistic.

### 8.4 Physical size (main-paper specification), distance, and MaNGA subsample

The main paper's primary specification is the strict specification above with galaxy size in kpc. The size control of Sections 8.1 to 8.3 is angular, and the 12-variant grid of Section 8.3 does not vary it. At fixed stellar mass the young quartile lies farther away than the old quartile in every mass bin (median angular-diameter distance `DA` 1.2 to 1.7 times larger) and is 2 to 3.5 times more often a MaNGA Secondary-sample target (`target` = 1; the sample has 2,776 Primary, 2,337 Secondary, and 839 color-enhanced entries), so an angular size carries part of a distance difference into the controls. The first table repeats the strict specification with size in kpc (converted with `DA`), with \(\log D_A\) added, and within the MaNGA Primary and Secondary samples; the second gives the bin counts of the full 12-variant grid under each size treatment.

| Specification | N | Observed SPS | Both SPS | Intrinsic SPS |
|---|---:|---:|---:|---:|
| Strict, size in arcsec (Sections 8.1 to 8.3; earlier versions) | 5,952 | 7/8, +0.0387, t = -11.8 | 7/8, +0.0344, t = -11.6 | 7/8, +0.0225, t = -8.2 |
| Strict, size in kpc (main paper) | 5,952 | 7/8, +0.0188, t = -10.2 | 6/8, +0.0162, t = -10.2 | 3/8, -0.0008, t = -5.5 |
| Strict plus \(\log D_A\) (either size unit) | 5,952 | 6/8, +0.0088, t = -8.9 | 3/8, +0.0057, t = -9.0 | 3/8, -0.0084, t = -4.2 |
| MaNGA Primary only, size in arcsec | 2,776 | 6/8, +0.0209, t = -5.9 | 6/8, +0.0217, t = -7.0 | 6/8, +0.0125, t = -5.5 |
| Primary only, size in kpc | 2,776 | 4/8, +0.0110, t = -4.7 | 4/8, +0.0138, t = -6.0 | 4/8, +0.0029, t = -3.7 |
| MaNGA Secondary only, size in arcsec | 2,337 | 4/8, -0.0032, t = -5.3 | 3/8, -0.0088, t = -4.0 | 2/8, -0.0192, t = -0.8 |
| Secondary only, size in kpc | 2,337 | 4/8, -0.0079, t = -4.5 | 3/8, -0.0133, t = -3.3 | 2/8, -0.0228, t = -0.2 |

Cells give young > old bins, the mass-adjusted residual, and the HC3 \(t\) of a linear `T50` term added to the same controls; subsample rows refit within the subsample. With \(\log R_e\) (arcsec) and \(\log D_A\) both in the model their coefficients are nearly equal (0.54 and 0.61, intrinsic SPS), as a dependence on physical size would produce.

| Metallicity, flattening term | Size in arcsec | Size in kpc | Either size plus \(\log D_A\) |
|---|---:|---:|---:|
| Mass-weighted, MGE ellipticity (strict) | 7, 7, 7 | 7, 6, 3 | 6, 3, 3 |
| Mass-weighted, NSA Sérsic axis ratio | 7, 7, 7 | 6, 6, 3 | 6, 6, 3 |
| Light-weighted, MGE ellipticity | 7, 7, 7 | 7, 7, 7 | 6, 6, 5 |
| Light-weighted, NSA Sérsic axis ratio | 7, 7, 7 | 7, 7, 7 | 6, 6, 6 |

Cells give the young > old bins (of 8) for the observed, both, and intrinsic SPS sets of the 12-variant grid, ties broken by stellar mass; the first row is the strict specification of the table above. The HC3 \(t\) of the linear `T50` term lies between -4.2 and -12.7 in all 48 fits.

With physical size, the main-paper specification, the observed-SPS residual keeps 7/8 bins and the both-SPS residual gives 6/8, each at about half its angular-size amplitude, and the intrinsic-SPS residual vanishes; with distance also controlled only the observed-SPS residual keeps a majority (6/8, about a quarter of the angular-size amplitude) and the both-SPS residual falls to 3/8. In the MaNGA Primary sample the kpc specification gives 4/8 in all three sets (6/8 with angular size, at 54 to 63% of the full-sample angular-size amplitude). With size in kpc, Spearman correlations of the unbinned residual with `T50` give \(\rho = -0.095\) (\(p \approx 2 \times 10^{-13}\)), \(-0.083\) (\(p \approx 1 \times 10^{-10}\)), and \(-0.027\) (\(p = 0.04\)) for the observed, both, and intrinsic sets. The linear `T50` term stays significant in every full-sample variant, including the 12-variant grid under all four size treatments (HC3 \(t\) from -4.2 to -12.7), and in the Primary sample; for intrinsic SPS with distance controlled it is the only evidence left (rank correlation with `T50` -0.003, \(p = 0.81\)). In the Secondary sample the bin test finds no young > old residual with either size unit: with angular size the intrinsic-SPS difference is reversed (-0.0192, bootstrap 95% CI -0.033 to -0.008), with size in kpc the intrinsic- and both-SPS differences are reversed (-0.0228, CI -0.042 to -0.014; -0.0133, CI -0.027 to -0.004), and the linear term is not significant for intrinsic SPS (\(t\) = -0.8 and -0.2). Under kpc size the 12-variant grid gives 3/8 to 7/8: the intrinsic 3/8 needs mass-weighted metallicity, and with light-weighted metallicity all three sets give 7/8. With distance controlled the grid gives 3/8 to 6/8: the three 3/8 cells all use mass-weighted metallicity, and 9 of the 12 variants give 5/8 or 6/8; with light-weighted metallicity the observed- and both-SPS residuals give 6/8, while the intrinsic-SPS residual falls to a tenth to a quarter of its angular-size amplitude. No distance-controlled variant reaches the 7/8 that the bin-count test needs (a 6/8 count gives count p 0.13 to 0.14 in the sensitivity runs). The magnitude test separates the strict sets: with distance controlled the observed-SPS residual's magnitude p is 0.010 (bootstrap 95% CI +0.001 to +0.017), the both-SPS residual's 0.055 (its interval crosses zero), and the intrinsic-SPS residual's 0.89; with size in kpc the observed and both sets stay at the 0.0002 floor and the intrinsic set gives 0.52. The distance-controlled counts also depend on how tied `T50` values at the quartile cuts are assigned; every table here orders them by stellar mass. With the ties broken at random instead (200 draws per variant), the strict specification with distance controlled gives 4/8 to 6/8 (median 5/8) for observed SPS and 3/8 to 6/8 (median 4/8) for both SPS, and none reaches 7/8 under either rule, while the angular-size counts stay at 7/8 in 598 of 600 draws (8/8 in the other two) and the kpc counts move by at most one bin. The mass-adjusted residuals of the strict distance-controlled rows move by at most 0.004 under random ties, and the linear `T50` term does not use the quartiles. Because distance and `T50` are correlated (Spearman \(\rho = -0.18\)), the \(\log D_A\) control also absorbs part of any real age signal, so its bin counts are partly conservative; with `T50` in the model, \(\log D_A\) still has an effect beyond physical size (\(t\) = 3.8 to 4.5). Whether the distance link reflects MaNGA's target selection or a distance-dependent systematic is not established; a replication should control distance or analyze the MaNGA samples separately.

---

## 9. The two-layer interpretation

*Note: The main paper (Section 4) contains the condensed version of this interpretation. The expanded discussion below provides additional context.*

### Layer 1: the visible stellar-population layer

Old galaxies have old stars. Old stellar populations are dimmer per unit stellar mass. Stellar-population synthesis models therefore assign higher \(M_*/L\) to older systems. The literature pass confirms that this is not merely plausible; it is the standard expectation of SPS modeling.

That layer is standard. This analysis does not dispute it. In the raw sledgehammer split, this layer carries most of the DML result. We lean into this rather than hide it: the visible layer is the part current astrophysics already explains.

This is where the wrinkle analogy helps, provided it is used carefully. A wrinkle can be explained by collagen loss, UV damage, elastin breakdown, hydration, and cellular repair. Those explanations are real. But explaining a wrinkle does not make it unrelated to aging. It identifies one pathway through which aging becomes visible.

The SPS layer should be treated the same way. Stellar evolution explains why older stellar populations have higher \(M_*/L\). That does not automatically prove that the age structure ends at the SPS layer. It proves only that the visible layer has a known stellar-population mechanism.

### Layer 2: the controlled dynamical residual layer

After the stellar-population layer is explicitly controlled, a smaller dynamical residual remains in a linear age term, which keeps the predicted sign in every full-sample fit. With galaxy size in kpc, young galaxies retain higher controlled dynamical residuals than old galaxies in 7 of 8 mass bins for the observed-SPS set and 6 of 8 for the both-SPS set, but not for the intrinsic-SPS set (Section 8.4).

This layer is the reason the result is worth documenting.

But it has to be stated with the main guardrail intact: raw JAM dynamical \(M/L\) alone did not show the simple young-greater-than-old pattern. In the raw component split, old galaxies had higher JAM dynamical \(M/L\) in all eight mass bins. The raw DML split passed because the SPS \(M_*/L\) age effect was stronger. Only after SPS, metallicity, and structural covariance were controlled did a smaller dynamical residual reappear in the predicted direction.

That is less dramatic than saying “both raw layers point the same way,” but it is more accurate and harder to attack.

### The one-box interpretation

We call the broader reading the one-box interpretation:

> The visible stellar-population layer and the controlled dynamical residual may not be two unrelated age correlations. They may be two age-organized layers of galactic aging measured through different instruments.

This is an interpretation, not a result. The result is the measured two-layer pattern. The interpretation asks whether the standard filing system is complete.

The standard explanation may be locally correct but globally incomplete. This note does not fight the standard explanation for Layer 1. It grants it. Then it asks whether Layer 2 should be filed as unrelated noise, IMF/SPS covariance, assembly history, or as a second age-organized layer that deserves joint interpretation.

---

## 10. Limitations

### 10.1 DML is not a clean gravitational observable

`DML` is built from JAM dynamical mass-to-light minus SPS stellar mass-to-light. Both sides carry model assumptions. Neither side is equivalent to direct weak lensing or a board-complete observed rotation curve.

The raw result should not be described as a gravity measurement. It is a mass-to-light discrepancy.

### 10.2 SPS age structure is a major confound

The largest layer of the signal is exactly what standard stellar-population modeling predicts: old stellar populations have higher stellar mass-to-light ratios. That is not a problem for the analysis, but it does limit the claim level.

The literature grounding makes this unavoidable. SPS models are designed to infer stellar mass, age, metallicity, star-formation history, dust effects, and IMF assumptions from galaxy light. A raw DML split that subtracts SPS \(M_*/L\) will inherit that age structure. We state plainly that the raw DML split is mostly SPS-driven (Section 7.1).

This is a concession, not a retreat. The SPS layer is the visible layer. The proof-of-concept begins only after that layer is granted and removed from the dynamical question.

### 10.3 IMF variation: open limitation, direction not established

If the stellar initial mass function varies systematically with galaxy age, velocity dispersion, or morphology, part of the residual could reflect stellar-population modeling mismatch rather than dynamical physics.

This is not hypothetical. DynPop V and prior ATLAS\(^{3D}\) work use JAM/SPS \(M/L\) offsets to study IMF variation, and they find a strong velocity-dispersion connection with weaker age dependence. The present pipeline controls for velocity dispersion and SPS \(M/L\), but it does not settle IMF shape, abundance-ratio, aperture, or spectral-feature systematics.

The direction of any IMF bias therefore depends on how the IMF varies at fixed controls. The controlled residual is the age dependence of the JAM/SPS offset at fixed controls, so an IMF heavier in older galaxies at fixed controls would push it toward old greater than young, the opposite of what we measure; the literature does not establish that premise. Dynamical studies tie the IMF normalization mainly to velocity dispersion (ATLAS\(^{3D}\) XX; DynPop V), which the controls absorb, and the weak positive age trend of McDermid et al. (2014) loses significance for \(\sigma > 130\) km s\(^{-1}\). DynPop V finds a weak positive correlation of its IMF mass-excess factor with age in early-type galaxies, without holding \(\sigma_e\) fixed; in its full sample the JAM \(M/L\) shows no age correlation at fixed \(\sigma_e\) while the SPS \(M_*/L\) depends strongly on age, which it suggests may reflect gas in late-type galaxies. Spectroscopic studies tie bottom-heavy IMFs to [Mg/Fe] as well as velocity dispersion (Conroy & van Dokkum, 2012); [Mg/Fe] is not controlled here, and the two approaches disagree for the same galaxies (Smith, 2014). The sign of any IMF contribution at fixed controls is therefore not established, and IMF variation remains an open explanation for the residual.

### 10.4 Dust conventions matter

Intrinsic and observed SPS mass-to-light definitions share the same SPS stellar mass and differ only in luminosity: dust-free for the intrinsic definition, attenuated (Milky Way foreground plus internal dust) for the observed one (DynPop II). The JAM dynamical \(M/L\) is normalized to SDSS \(r\)-band light not corrected for Galactic or internal extinction (DynPop I), which is why DynPop V corrects it before comparing it with intrinsic SPS. In its dust treatment the observed-SPS DML is therefore the like-for-like comparison, and \(DML_{int}\) equals \(DML_{obs}\) plus each galaxy's attenuation term. The raw intrinsic split is weaker by bin count (6/8 against 8/8) but still positive overall.

This reverses the reading in which observed SPS is the dust-contaminated definition. With intrinsic-SPS controls the attenuation term stays uncontrolled in the dynamical \(M/L\); the both-SPS set controls it, since its two SPS terms span the intrinsic value plus that term. That raises the residual from +0.0225 to +0.0344 with size in arcsec (observed SPS: +0.0387) and from -0.0008 to +0.0162 with size in kpc (observed SPS: +0.0188). Young galaxies are not uniformly dustier here: in the low-mass bins the old quartile carries the larger attenuation term. The attenuation fitted by the SPS models therefore does not create the young-greater-than-old residual. Dust those fits do not model, such as dust lanes in the MGE photometry, is not tested here.

### 10.5 JAM covariance remains open

The same stellar kinematic maps contribute to \(\lambda_R\), velocity dispersion, and JAM dynamical quantities. Shared measurement inputs can manufacture or tighten correlations. Separating a physical signal from this shared-input covariance would require mock kinematic maps and a full covariance treatment, which this note does not attempt.

### 10.6 Structural covariance is strong

Age, mass, size, surface density, morphology, Sérsic structure, metallicity, velocity dispersion, and orbital organization are all entangled in galaxy populations. The controlled residual models reduce this problem but do not eliminate it.

### 10.7 fDM mildly cuts against a gravitational reading

The NFW and gNFW dark-matter fractions within \(R_e\) do not show the same age signal. For a framework that interprets the mass-to-light discrepancy as gravitational, the absence of a corresponding signal in the parametric dark-matter fraction is not neutral. It mildly cuts against the gravitational interpretation, because a gravitational age effect large enough to appear in the controlled dynamical residual might be expected to leave some trace in the halo decomposition as well. This does not refute the DML residual, but it blocks any simple claim that MaNGA directly shows age-dependent dark-matter fraction, and it should be weighed honestly against the residual when assessing the result.

### 10.8 No cosmological-mock replication yet

The strongest standard-model comparison would use mock MaNGA observations from hydrodynamic simulations, processed through a comparable analysis. Section 3.4 reviews the MaNGA-like simulation tools and comparisons that already exist.

What has not been done in this note, and what we did not locate in the literature pass, is the exact comparator needed here: a mock MaNGA/DynPop-like sample processed through comparable JAM/SPS observables, then split by age at fixed stellar mass and tested for the same controlled dynamical \(M/L\) residual.

Therefore the residual could still be produced by:

- ordinary \(\Lambda\)CDM assembly bias;
- halo concentration and quenching history;
- feedback-driven structure changes;
- IMF variation tied to velocity dispersion;
- morphology and merger history;
- JAM covariance;
- unmodeled cold gas mass covarying with age;
- distance-dependent selection or measurement systematics (Section 8.4);
- mock/observational pipeline differences.

This is the main hardening step that remains outside the present pipeline. Until this comparison is performed, the residual should be described as unexplained by this analysis, not unexplained by \(\Lambda\)CDM.

### 10.9 This is not DESI weak lensing

Weak lensing is the cleaner external test because it measures light deflection. MaNGA is useful because it is public, rich, and already contains ages, dynamics, and stellar-population outputs. But it is entangled with SPS and JAM assumptions in a way lensing is not.


### 10.10 The novelty is narrow

The relevant DynPop literature already shows that stellar population age, stellar-population mass-to-light ratio, dynamical mass-to-light ratio, total density slope, circular-velocity-curve shape, velocity dispersion, dark-matter fraction, and IMF-sensitive JAM-SPS offsets are mutually entangled. The present note should not imply otherwise. Its contribution is the deliberately blunt fixed-mass young/old DML split and the two-layer hardening sequence, not the discovery that MaNGA contains age-dependent dynamical scaling relations.

### 10.11 The one-box interpretation is not a measured fact

The analysis measures a two-layer pattern. It does not measure a common mechanism behind the two layers. The one-box interpretation is a way to organize the question, not an answer to it. A standard-model account could still explain the residual through IMF variation, assembly history, structural covariance, SPS mismatch, dust treatment, cold gas fraction, JAM assumptions, or simulation-predicted age correlations. This note invites those tests rather than preempting them.

### 10.12 Unmodeled cold gas mass

The SPS denominator in DML accounts for stellar mass, but the JAM dynamical numerator measures total mass driving the kinematics, including dark matter and cold gas (\(H I\) and \(H_2\)). Young, star-forming galaxies have systematically higher cold gas fractions (e.g., Saintonge et al., 2017; Catinella et al., 2018), so unmodeled gas mass will naturally inflate DML for the young quartile relative to the old. While cold gas mass is generally subdominant to stellar mass within the inner effective radius where these DynPop quantities are evaluated, any systematic covariance between gas fraction and age remains a standard-model contributor to the controlled residual. Future mock comparators should include gas-fraction scaling.

## 11. Hand-off statement

See the main paper, Section 5. Deciding what the residual is requires domain infrastructure this team does not currently have: stellar-population synthesis expertise, IMF-systematics modeling, JAM covariance treatment, cosmological simulation mocks, and independent gravitational measurements. The goal is not to own a MaNGA anomaly.

---

## 12. Data and reproducibility package

The package is the `paper10/` folder of github.com/mcgintyl/hermes-verification (commit 1094bdd). That commit contains the joined working table (10,296 rows), `build_merged_table.py`, which rebuilds that table from the three public files of Section 4.1 after checking their checksums, the analysis scripts, the archived outputs, and a manifest with the SHA-256 hash of every file. `python reproduce.py` reruns every analysis and checks each analysis number quoted in the main paper and this appendix against the regenerated values. The code behind the outputs of the earlier v0.6 pipeline (Sections 6 and 7 and Table 3; archived in `01_clean_pipeline_v0_6/`) was not preserved; `run_sledgehammer.py` and the kitchen-sink variant of `run_size_distance_sensitivity.py` re-implement it and reproduce every v0.6 point estimate exactly, while shuffle p-values, bootstrap intervals, and cross-validation scores agree only to within their random-draw precision. The archived v0.6 values remain the record.

---

## 13. Figures

This version carries no figures. The plots from the v0.6 pipeline remain in the package (`01_clean_pipeline_v0_6/`).

---

## 14. Detailed archived tables and sensitivity scoreboards

### Table 1. Raw DML sledgehammer scoreboard

| Outcome | Young > old bins | Mean young-old median difference | Mass-adjusted young-old difference | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| `DML_mfl_obs` | 8 / 8 | +0.0742 | +0.0669 | [+0.0498, +0.0828] | 0.006 | 0.0002 |
| `DML_mfl_int` | 6 / 8 | +0.0606 | +0.0772 | [+0.0562, +0.0959] | 0.15 | 0.0002 |

### Table 2. Raw M/L split scoreboard

| Outcome | Young > old bins | Old > young bins | Direction relevant to DML | Directional bins | Mass-adjusted young-old difference | Bootstrap 95% CI |
|---|---:|---:|---|---:|---:|---:|
| JAM dynamical \(\log M/L\) | 0 / 8 | 8 / 8 | young > old | 0 / 8 | -0.1075 | [-0.1200, -0.0920] |
| SPS intrinsic \(\log M_*/L\) | 0 / 8 | 8 / 8 | old > young | 8 / 8 | -0.1693 | [-0.1790, -0.1587] |
| SPS observed \(\log M_*/L\) | 0 / 8 | 8 / 8 | old > young | 8 / 8 | -0.1569 | [-0.1679, -0.1478] |

### Table 3. Controlled dynamical residual scoreboard, kitchen-sink specification (sensitivity only)

| Controlled residual outcome | Young > old bins | Mean young-old residual | Mass-adjusted young-old residual | Bootstrap 95% CI | Bin-count p | Magnitude p |
|---|---:|---:|---:|---:|---:|---:|
| After observed SPS controls | 8 / 8 | +0.0359 | +0.0358 | [+0.0270, +0.0475] | 0.004 | 0.0002 |
| After both SPS controls | 8 / 8 | +0.0342 | +0.0311 | [+0.0216, +0.0436] | 0.004 | 0.0002 |
| After intrinsic SPS controls | 6 / 8 | +0.0209 | +0.0188 | [+0.0092, +0.0312] | 0.14 | 0.0002 |

These rows use the kitchen-sink control set, which adds light-weighted metallicity (`sp_LW_Metal_Re`), a second flattening term (`nsa_sersic_ba`), and raw `Sigma_Re` to the strict controls of Section 8, with size in arcsec. Compared with the strict specification (7/8 bins in all three sets, Section 8.3), it gains one bin in the observed-SPS and both-SPS sets and loses one in the intrinsic set; the direction of the residual is unchanged. Its `T50` coefficients for the observed, both, and intrinsic sets are -0.0346, -0.0357, and -0.0215 dex per 1 SD of `T50` (model \(R^2\) 0.4799 to 0.5058), with cross-validated gains from age of +0.0146, +0.0137, and +0.0052.

---

## 15. Conclusion

*The main paper's hand-off (Section 5) is the authoritative conclusion.*

The raw DML split is mostly the expected stellar-population layer. After that layer is controlled, a smaller dynamical residual remains in a linear age term; whether the young-old split shows it, and in how many mass bins, depends on how galaxy size and distance are controlled and on the MaNGA sample (Section 8.4).

That is the finding. Not proof. Not a claim of solved gravity. A two-layer age pattern in a public dynamical catalog, with the obvious stellar-population layer acknowledged and a smaller residual layer documented for independent groups to test.

Whether the filing system that separates those two layers reflects the structure of the phenomenon, or merely the structure of the instruments used to measure it, remains open.

---

## 16. Remaining hardening and future work

Two hardening steps would strengthen the result. The larger of them, the standard-model mock comparator, needs specialist infrastructure and remains future work.

### 16.1 Quantitative SPS expectation

The measured SPS \(M_*/L\) difference across the young-old split is in Section 7.1, the per-bin split of the raw DML difference into its dynamical and SPS parts is in the package (v0.6 hardening report, "DML component decomposition by bin"), and Section 10.4 sets out the dust convention. A model-based estimate of the expected age-driven SPS \(M_*/L\) difference remains open; it should keep the DynPop V warning that JAM/SPS comparison mixes IMF, dark matter, and modeling systematics.

### 16.2 Simulation and assembly-bias comparator status

The literature pass found adjacent but not exact standard-model comparators (Section 3.4). What remains missing:

- a mock MaNGA/DynPop-like catalog processed through comparable JAM/SPS mass-to-light measurements;
- the same fixed-stellar-mass young/old quartile split;
- the same DML numerator/denominator separation;
- the same controlled dynamical residual test;
- a reported mock distribution of residual amplitudes against which the observed controlled young-old residuals of Sections 8.3 and 8.4 can be compared.

## 17. Supplementary references

Zhu, K., Lu, S., Cappellari, M., Li, R., Mao, S., & Gao, L. (2023). *MaNGA DynPop -- I. Quality-assessed stellar dynamical modelling from integral-field spectroscopy of 10K nearby galaxies: a catalogue of masses, mass-to-light ratios, density profiles and dark matter.* Monthly Notices of the Royal Astronomical Society, 522, 6326-6353. arXiv:2304.11711.

Lu, S., Zhu, K., Cappellari, M., Li, R., Mao, S., & Xu, D. (2023). *MaNGA DynPop -- II. Global stellar population, gradients, and star-formation histories from integral-field spectroscopy of 10K galaxies: link with galaxy rotation, shape, and total-density gradients.* Monthly Notices of the Royal Astronomical Society, 526, 1022-1045. arXiv:2304.11712.

Zhu, K., Lu, S., Cappellari, M., Li, R., Mao, S., Gao, L., & Ge, J. (2024). *MaNGA DynPop -- III. Stellar dynamics versus stellar population relations in 6000 early-type and spiral galaxies: Fundamental Plane, mass-to-light ratios, total density slopes, and dark matter fractions.* Monthly Notices of the Royal Astronomical Society, 527, 706-730. arXiv:2304.11714.

Lu, S., Zhu, K., Cappellari, M., Li, R., Mao, S., & Xu, D. (2024). *MaNGA DynPop -- V. The dark-matter fraction versus stellar velocity dispersion relation and stellar initial mass function variations in galaxies: dynamical models and full spectrum fitting of integral-field spectroscopy.* Monthly Notices of the Royal Astronomical Society, 530, 4474-4492. arXiv:2309.12395.

Zhu, K., Cappellari, M., Mao, S., Lu, S., Li, R., Shi, Y., Simon, D. A., Fu, Y., & Wang, X. (2025). *MaNGA DynPop -- VII. A Unified Bulge-Disk-Halo Model for Explaining Diversity in Circular Velocity Curves of 6000 Spiral and Early-Type Galaxies.* The Astrophysical Journal Supplement Series, 280, 55. arXiv:2503.06968.

Bruzual, G., & Charlot, S. (2003). *Stellar population synthesis at the resolution of 2003.* Monthly Notices of the Royal Astronomical Society, 344, 1000-1028. arXiv:astro-ph/0309134.

Bell, E. F., & de Jong, R. S. (2001). *Stellar mass-to-light ratios and the Tully-Fisher relation.* The Astrophysical Journal, 550, 212-229. arXiv:astro-ph/0011493.

Vazdekis, A., Sánchez-Blázquez, P., Falcón-Barroso, J., Cenarro, A. J., Beasley, M. A., Cardiel, N., Gorgas, J., & Peletier, R. F. (2010). *Evolutionary stellar population synthesis with MILES. I. The base models and a new line index system.* Monthly Notices of the Royal Astronomical Society, 404, 1639-1671. arXiv:1004.4439.

Conroy, C. (2013). *Modeling the panchromatic spectral energy distributions of galaxies.* Annual Review of Astronomy and Astrophysics, 51, 393-455. arXiv:1301.7095.

Cappellari, M., et al. (2012). *Systematic variation of the stellar initial mass function in early-type galaxies.* Nature, 484, 485-488.

Cappellari, M., et al. (2013). *The ATLAS3D project — XX. Mass-size and mass-sigma distributions of early-type galaxies: bulge fraction drives kinematics, mass-to-light ratio, molecular gas fraction and stellar initial mass function.* Monthly Notices of the Royal Astronomical Society, 432, 1862-1893.

McDermid, R. M., et al. (2014). *Connection between dynamically derived initial mass function normalization and stellar population parameters.* The Astrophysical Journal Letters, 792, L37. arXiv:1408.3189.

Smith, R. J. (2014). *Variations in the initial mass function in early-type galaxies: a critical comparison between dynamical and spectroscopic results.* Monthly Notices of the Royal Astronomical Society, 443, L69-L73. arXiv:1403.6114.

Nanni, L., Thomas, D., Trayford, J., Maraston, C., Neumann, J., Law, D. R., Hill, L., Pillepich, A., Yan, R., Chen, Y., & Lazarz, D. (2022). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. I. Construction and analysis of the mock data cubes.* Monthly Notices of the Royal Astronomical Society, 515, 320-338. arXiv:2203.11575.

Nanni, L., Thomas, D., Trayford, J., Maraston, C., Neumann, J., Law, D. R., Hill, L., Pillepich, A., Yan, R., Chen, Y., & Lazarz, D. (2023). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. II. The catalogue.* Monthly Notices of the Royal Astronomical Society, 522, 5479-5499. arXiv:2211.13146.

Nanni, L., Neumann, J., Thomas, D., Maraston, C., Trayford, J., Lovell, C. C., Law, D. R., Yan, R., & Chen, Y. (2024). *iMaNGA: mock MaNGA galaxies based on IllustrisTNG and MaStar SSPs. III. Stellar metallicity drivers in MaNGA and TNG50.* Monthly Notices of the Royal Astronomical Society, 527, 6419-6438. arXiv:2309.14257.

Li, R., Li, H., Shao, S., Lu, S., Zhu, K., Wang, C., Gao, L., Mao, S., Dutton, A. A., Ge, J., Wang, Y., Leauthaud, A., Zheng, Z., Bundy, K., & Brownstein, J. R. (2019). *SDSS-IV MaNGA: The Inner Density Slopes of nearby galaxies.* Monthly Notices of the Royal Astronomical Society, 490, 2124-2138. arXiv:1903.09282.

Li, S., Li, R., Zhu, K., Lu, S., Cappellari, M., Mao, S., Wang, C., & Gao, L. (2024). *MaNGA DynPop -- VI. Matter density slopes from dynamical models of 6000 galaxies versus cosmological simulations: the interplay between baryonic and dark matter.* Monthly Notices of the Royal Astronomical Society, 529, 4633-4649. arXiv:2310.13278.

Mukherjee, S., Koopmans, L. V. E., Tortora, C., Schaller, M., Metcalf, R. B., Schaye, J., & Vernardos, G. (2022). *SEAGLE-III: Towards resolving the mismatch in the dark-matter fraction in early-type galaxies between simulations and observations.* Monthly Notices of the Royal Astronomical Society, 509, 1245-1251. arXiv:2110.07615.

Tinker, J. L., Wetzel, A. R., Conroy, C., & Mao, Y.-Y. (2017). *Halo histories versus galaxy properties at z = 0. I. The quenching of star formation.* Monthly Notices of the Royal Astronomical Society, 472, 2504-2516. arXiv:1609.03388.

Lovell, M. R., Pillepich, A., Genel, S., Nelson, D., Springel, V., Pakmor, R., Marinacci, F., Weinberger, R., Torrey, P., Vogelsberger, M., et al. (2018). *The fraction of dark matter within galaxies from the IllustrisTNG simulations.* Monthly Notices of the Royal Astronomical Society, 481, 1950-1975. arXiv:1801.10170.

de Graaff, A., Franx, M., Bell, E. F., Bezanson, R., Schaller, M., Schaye, J., & van der Wel, A. (2023). *A common origin for the Fundamental Plane of quiescent and star-forming galaxies in the EAGLE simulations.* Monthly Notices of the Royal Astronomical Society, 518, 5376-5402. arXiv:2207.13491.

Montero-Dorta, A. D., Chaves-Montero, J., Artale, M. C., & Favole, G. (2021). *On the influence of halo mass accretion history on galaxy properties and assembly bias.* Monthly Notices of the Royal Astronomical Society, 508, 940-949. arXiv:2105.05274.

Xu, X., & Zheng, Z. (2020). *Galaxy assembly bias of central galaxies in the Illustris simulation.* Monthly Notices of the Royal Astronomical Society, 492, 2739-2754. arXiv:1812.11210.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.* The Astronomical Journal, 152, 157. arXiv:1606.09251.

McGinty, L. A. (2026a). *The Hermes Equation: Galaxy Rotation Curves from Age and Density with No Free Constants.* Zenodo. DOI (concept): 10.5281/zenodo.18809176.

McGinty, L. A. (2026b). *Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing.* Zenodo. DOI (concept): 10.5281/zenodo.20550091.

McGinty, L. A. (2026c). *Age-Dependent Weak Lensing in Disk Galaxies: A Pilot Measurement and DESI Roadmap.* Zenodo. DOI (concept): 10.5281/zenodo.19154823.

Conroy, C., & van Dokkum, P. G. (2012). *The Stellar Initial Mass Function in Early-Type Galaxies from Absorption Line Spectroscopy. II. Results.* The Astrophysical Journal, 760, 71. arXiv:1205.6473.

Saintonge, A., Catinella, B., Tacconi, L. J., Kauffmann, G., Genzel, R., Cortese, L., Davé, R., Fletcher, T. J., Graciá-Carpio, J., Kramer, C., Heckman, T. M., Janowiecki, S., Lutz, K., Rosario, D., Schiminovich, D., Schuster, K., Wang, J., Wuyts, S., Borthakur, S., Lamperti, I., & Roberts-Borsani, G. W. (2017). *xCOLD GASS: The Complete IRAM 30 m Legacy Survey of Molecular Gas for Galaxy Evolution Studies.* The Astrophysical Journal Supplement Series, 233, 22. arXiv:1710.02157.

Catinella, B., Saintonge, A., Janowiecki, S., Cortese, L., Davé, R., Lemonias, J. J., Cooper, A. P., Schiminovich, D., Hummels, C. B., Fabello, S., Geréb, K., Kilborn, V., & Wang, J. (2018). *xGASS: total cold gas scaling relations and molecular-to-atomic gas ratios of galaxies in the local Universe.* Monthly Notices of the Royal Astronomical Society, 476, 875-895. arXiv:1802.02373.
