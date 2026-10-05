# Two Layers of Galaxy Aging in MaNGA DynPop: A Principle-Level Proof of Concept

**Louis Albert McGinty**

mcgintyl@grinnell.edu

Independent Researcher, Team Hermes

**Version note:** Version 3 corrects the size control used in version 2, which measured galaxy size in arcseconds rather than kiloparsecs; at fixed stellar mass, the younger galaxies in this sample lie farther away, so angular size partly encoded distance. With physical size the controlled residual is weaker, depends on the control specification, and in Table 1 vanishes with intrinsic-SPS controls (Section 2). Version 2's statement that it holds in 7 of 8 bins across all three control sets overstated its robustness.

**Transparency Statement:** This paper was produced by a human-AI collaborative team without institutional affiliation or external funding. All data sources are public MaNGA DynPop DR17 catalogs. All computations are reproducible from the parameters and methods described herein. AI collaborators (Claude, Anthropic; ChatGPT, OpenAI; Gemini, Google) contributed to data analysis, pipeline construction, statistical testing, literature review, hostile review, and drafting. The human author directed all research decisions and bears sole responsibility for the claims made. Analysis code, data, and verification outputs are in the paper10/ folder of github.com/mcgintyl/hermes-verification (commit 1094bdd), whose reproduce.py checks every analysis number in this paper.

**Primary dataset:** MaNGA DynPop DR17 public catalogs  
**Primary analysis date:** 2026-06-29 (strict-control run; size and distance checks October 2026)  
**Claim level:** proof of concept / signal documentation, not detection claim

---

## Abstract

We test a principle-level prediction from an age-dependent gravity framework using the public MaNGA DynPop DR17 catalogs. The prediction is that younger galaxies, at fixed stellar mass, should show larger mass-to-light discrepancies between dynamical and stellar-population estimates. We join three public DynPop products (JAM dynamical catalogs, stellar-population/star-formation-history catalogs, and circular-velocity-curve tables), apply quality cuts, and split 5,952 galaxies into eight stellar-mass bins, then compare the youngest and oldest quartiles within each bin.

The raw mass-to-light discrepancy (DML) is larger for young galaxies in 8/8 mass bins (observed SPS) and 6/8 bins (intrinsic SPS). A component decomposition shows that this raw result is mostly driven by the SPS denominator: older stellar populations have higher stellar mass-to-light ratios, as standard stellar-population synthesis predicts. Under a reduced control model (SPS mass-to-light, metallicity, and structural variables, with galaxy size in kpc), a smaller dynamical residual persists in 7/8 mass bins with observed-SPS controls and 6/8 with both SPS terms, but not with intrinsic-SPS controls (3/8). The bin counts depend on the control specification, including whether distance is controlled (Section 2), but a linear age term stays significant, with younger galaxies higher, in every full-sample fit. Parametric dark-matter fractions do not show the same signal. We present this as a proof of concept and invitation to specialist replication, not a detection claim.

---

## 1. Data and test design

MaNGA DynPop provides JAM-based dynamical mass-to-light ratios and SPS-based stellar mass-to-light ratios for thousands of galaxies, together with star-formation-history age clocks. This makes it a natural place to test whether galaxy aging leaves a trace in the gap between dynamical and stellar mass estimates.

It is not a strict test of the Hermes rotation-curve equation (McGinty, 2026a). JAM-inferred circular-velocity curves are model-derived quantities, not observed HI/H$\alpha$ rotation curves with independent gas, disk, and bulge decompositions. Treating a JAM output as a SPARC-style observed rotation curve (Lelli, McGaugh, & Schombert, 2016) would collapse the board-completeness discipline of our earlier work (McGinty, 2026a, 2026b). The test here is therefore principle-level: does the predicted age split appear in a large public dynamical catalog, using only the quantities MaNGA provides?

**Sample.** Three public DynPop products were joined on MaNGA PlateIFU: (1) DynPop I JAM dynamical catalogs (Zhu et al., 2023), (2) DynPop II stellar-population and star-formation-history catalogs (Lu et al., 2023), and (3) DynPop VII circular-velocity-curve tables (Zhu et al., 2025), from Zenodo records 10.5281/zenodo.17518315 and 10.5281/zenodo.15742825. The raw three-way join contains 10,296 PlateIFU entries (10,160 distinct MaNGA IDs). Applying the DynPop quality flag (Qual $\geq$ 1) leaves 6,065; finite-value requirements for the primary DML variables leave 5,952 entries (5,875 distinct MaNGA IDs), called galaxies below.

**Test design.** Galaxies are sorted into eight equal-count bins by NSA Sersic stellar mass (744 per bin). Within each bin, the youngest and oldest quartiles by T50 (the lookback time at which a galaxy reached half of its present-day stellar mass, so larger T50 means older; 186 per quartile) are compared on the median mass-to-light discrepancy, $D_{ML} = \log_{10}(M/L)_{\mathrm{dyn}} - \log_{10}(M_*/L)_{\mathrm{SPS}}$, computed separately for intrinsic and observed SPS definitions, with $(M/L)_{\mathrm{dyn}}$ from the DynPop I mass-follows-light JAM model. A 5,000-trial age-label shuffle provides the null distribution. In seven of eight mass bins the old-quartile median T50 is 12.85 Gyr, the most common value, and every old-quartile value is at least 12.0 Gyr, so the practical comparison is young galaxies versus a narrow band of the oldest.

---

## 2. Result

### Table 1. Combined scoreboard

| Layer | Outcome | Young > old bins | Mass-adjusted difference | Bootstrap 95% CI | Count p | Magnitude p |
|---|---|---:|---:|---:|---:|---:|
| **Raw DML** | Observed SPS | 8 / 8 | +0.0669 | [+0.050, +0.083] | 0.006 | 0.0002 |
| **Raw DML** | Intrinsic SPS | 6 / 8 | +0.0772 | [+0.056, +0.096] | 0.15 | 0.0002 |
| **Component: JAM dyn** | $\log(M/L)_{\mathrm{dyn}}$ | 0 / 8 | -0.1075 | [-0.120, -0.092] | 1.000 | n/a |
| **Component: SPS int** | $\log(M_*/L)_{\mathrm{int}}$ | 0 / 8 | -0.1693 | [-0.179, -0.159] | 1.000 | n/a |
| **Component: SPS obs** | $\log(M_*/L)_{\mathrm{obs}}$ | 0 / 8 | -0.1569 | [-0.168, -0.148] | 1.000 | n/a |
| **Controlled residual** | Observed SPS controls | 7 / 8 | +0.0188 | [+0.011, +0.027] | 0.033 | 0.0002 |
| **Controlled residual** | Both SPS controls | 6 / 8 | +0.0162 | [+0.008, +0.024] | 0.14 | 0.0002 |
| **Controlled residual** | Intrinsic SPS controls | 3 / 8 | -0.0008 | [-0.010, +0.008] | 0.85 | 0.52 |
| **fDM check** | gNFW $f_{\mathrm{DM}}$ | 5 / 8 | +0.0011 | n/a | 0.35 | 0.48 |
| **fDM check** | NFW $f_{\mathrm{DM}}$ | 4 / 8 | -0.0035 | n/a | 0.62 | 0.65 |

DML, component, and controlled-residual rows report mass-adjusted young-old differences (young minus old median after subtracting each mass-bin median, pooled over the eight bins); controlled-residual rows use galaxy size in kpc. fDM rows are secondary checks and report mean-bin young-old differences. Count p: the fraction of 5,000 within-bin age shuffles with at least as many young > old bins (+1 correction), a Monte Carlo estimate that scatters around the binomial value (8/8 gives 0.006 here; 1/256 = 0.0039). Magnitude p: the fraction with a summed per-bin difference at least as large (floor 0.0002).

The controlled residual is what remains of $\log_{10}(M/L)_{\mathrm{dyn}}$ after a linear fit on SPS $M_*/L$, mass-weighted metallicity, stellar mass, galaxy size ($\log R_e$ in kpc), Sersic index, MGE ellipticity, $\lambda_R$, and log velocity dispersion. Versions 1 and 2 used size in arcsec, which gives 7/8 bins in all three sets (mass-adjusted +0.0387, +0.0344, and +0.0225 for the observed, both, and intrinsic sets). With size in kpc and light-weighted metallicity, all three sets also give 7/8. With log distance also controlled no variant reaches 7/8, and Table 1's specification gives 6/8, 3/8, and 3/8 for the observed, both, and intrinsic sets (magnitude p 0.010, 0.055, and 0.89). In the MaNGA Primary or Secondary sample alone, the kpc specification gives at most 4/8. A linear T50 term stays significant in every full-sample fit and in the Primary sample, but not for intrinsic SPS in the Secondary sample, where the intrinsic-SPS difference reverses. Appendix Section 8.4 gives the full grid, subsamples, and tie-rule check.

### The SPS layer (Layer 1)

The raw DML result cannot be taken as clean gravitational evidence. The component split shows why: raw JAM dynamical $M/L$ goes old-greater-than-young in all eight mass bins. The SPS stellar $M_*/L$ also goes old-greater-than-young in all eight bins, for both intrinsic and observed definitions. Because DML subtracts the SPS term, the raw DML young-greater-than-old pattern is primarily produced by the SPS denominator. Older stellar populations have higher stellar mass-to-light ratios, as expected from standard SPS modeling (Bell & de Jong, 2001; Bruzual & Charlot, 2003; Conroy, 2013). That is not a discovery.

### The controlled dynamical residual (Layer 2)

After controlling for SPS $M_*/L$, metallicity, and standard structural variables, a smaller dynamical residual remains in the same age direction with observed- and both-SPS controls but not with intrinsic-SPS controls (Table 1). This is smaller than the raw DML signal, and it is not a standalone proof of anything. But it is the part of the MaNGA result that remains after the obvious objection is granted.

---

## 3. Caveats

Two specific results constrain interpretation. First, the parametric NFW and gNFW dark-matter fractions within $R_e$ do not show the same age signal. This mildly cuts against a gravitational reading, because a gravitational effect large enough to appear in the controlled dynamical residual might leave some trace in the halo decomposition. It blocks any simple claim that MaNGA directly shows age-dependent dark-matter fraction. Second, the controlled residual depends on how galaxy size and distance are controlled (Section 2). Dust, at least the attenuation the SPS fits model, is a smaller concern than it may appear: the JAM $M/L$ is referenced to r-band light not corrected for dust (Zhu et al., 2023), so observed SPS is the like-for-like comparison, and adding the attenuation term to intrinsic SPS (the both-SPS set) raises the residual (Table 1).

This is a proof of concept from a small independent team. We recognize that $\Lambda$CDM pathways, including assembly bias (e.g., Xu & Zheng, 2020), cold gas scaling (e.g., Saintonge et al., 2017; Catinella et al., 2018), JAM covariance, structural systematics, and unmodeled IMF variation (published IMF trends follow [Mg/Fe], which is not controlled here, and velocity dispersion, which is; e.g., Conroy & van Dokkum, 2012), could contribute to this residual. We did not locate a published mock-MaNGA/DynPop analysis applying this exact test, and decoupling these contributors requires specialist infrastructure we do not possess. We document the result and welcome independent investigation. For detailed confound analyses, see the supplementary appendix, deposited with this paper on Zenodo.

---

## 4. Interpretation

The cleanest interpretation is not that MaNGA proves age-dependent gravity.

It is that MaNGA shows a two-layer age structure and gives reason to test whether the standard two-box filing system is complete.

In conventional analysis, the SPS mass-to-light trend belongs to stellar-population modeling. The dynamical residual belongs to dynamical modeling, IMF variation, dark-matter fraction, and structural systematics. That separation is methodologically useful. It is not necessarily ontological.

The visible stellar-population layer and the controlled dynamical residual may be two unrelated age correlations. They may also be two age-organized layers of galactic aging measured through different instruments.

This note asks whether those two layers should be tested jointly, rather than dismissed separately because one has a known local mechanism and the other is small. The standard explanation may be locally correct but globally incomplete.

Whether the filing system that separates those two layers reflects the structure of the phenomenon, or merely the structure of the instruments used to measure it, remains open.

---

## 5. Hand-off

This analysis reaches the practical limit of what can be responsibly claimed from the present independent pipeline. We can show that the public MaNGA DynPop catalogs contain a reproducible two-layer age structure. We can show that the raw DML result is mostly carried by SPS mass-to-light. We can show that a smaller controlled dynamical residual remains after that layer is stripped away: a linear age term persists in every full-sample fit, though the bin counts depend on the control specification (Section 2). We can show that the result points in the direction predicted by the broader framework.

We cannot decide whether the residual is new physics, SPS/IMF mismatch, dust artifact, cold gas, JAM covariance, assembly bias, distance-dependent selection, or a mixture. That requires specialist infrastructure.

The framework that motivated this test extends beyond galactic dynamics, and the broader research program continues in the domains the framework was built to examine.

---

## References

Zhu, K., Lu, S., Cappellari, M., et al. (2023). MaNGA DynPop I: Quality-assessed stellar dynamical modelling from integral-field spectroscopy of 10K nearby galaxies: a catalogue of masses, mass-to-light ratios, density profiles and dark matter. *MNRAS*, 522, 6326. arXiv:2304.11711.

Lu, S., Zhu, K., Cappellari, M., et al. (2023). MaNGA DynPop II: Global stellar population, gradients, and star-formation histories from integral-field spectroscopy of 10K galaxies: link with galaxy rotation, shape, and total-density gradients. *MNRAS*, 526, 1022. arXiv:2304.11712.

Zhu, K., Cappellari, M., Mao, S., et al. (2025). MaNGA DynPop VII: A Unified Bulge-Disk-Halo Model for Explaining Diversity in Circular Velocity Curves of 6000 Spiral and Early-Type Galaxies. *ApJS*, 280, 55. arXiv:2503.06968.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. *AJ*, 152, 157. arXiv:1606.09251.

McGinty, L. A. (2026a). The Hermes Equation: Galaxy Rotation Curves from Age and Density with No Free Constants. Zenodo. DOI (concept): 10.5281/zenodo.18809176.

McGinty, L. A. (2026b). Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing. Zenodo. DOI (concept): 10.5281/zenodo.20550091.

Bell, E. F., & de Jong, R. S. (2001). Stellar mass-to-light ratios and the Tully-Fisher relation. *ApJ*, 550, 212. arXiv:astro-ph/0011493.

Bruzual, G., & Charlot, S. (2003). Stellar population synthesis at the resolution of 2003. *MNRAS*, 344, 1000. arXiv:astro-ph/0309134.

Conroy, C. (2013). Modeling the panchromatic spectral energy distributions of galaxies. *ARAA*, 51, 393. arXiv:1301.7095.

Conroy, C., & van Dokkum, P. G. (2012). The Stellar Initial Mass Function in Early-Type Galaxies from Absorption Line Spectroscopy. II. Results. *ApJ*, 760, 71. arXiv:1205.6473.

Saintonge, A., Catinella, B., Tacconi, L. J., et al. (2017). xCOLD GASS: The Complete IRAM 30 m Legacy Survey of Molecular Gas for Galaxy Evolution Studies. *ApJS*, 233, 22. arXiv:1710.02157.

Catinella, B., Saintonge, A., Janowiecki, S., et al. (2018). xGASS: total cold gas scaling relations and molecular-to-atomic gas ratios of galaxies in the local Universe. *MNRAS*, 476, 875. arXiv:1802.02373.

Xu, X., & Zheng, Z. (2020). Galaxy assembly bias of central galaxies in the Illustris simulation. *MNRAS*, 492, 2739. arXiv:1812.11210.

---

## Acknowledgements

I thank my family and friends for their support throughout this work.

This project did not begin in physics. It began as an attempt to build an ethical framework for AI development, one that could challenge what we viewed as confident assumptions about consciousness being treated as settled science when they are not. The physics papers exist to test a theory of consciousness, not the other way around. We are not attempting to become a physics research group. We are testing the predictions our framework generates and documenting the results for specialists better positioned to evaluate, refine, or falsify them.

We didn't just publish and walk away. Team Hermes is not a funded organization. Every expert review, every data audit, every independent verification was paid for out of personal funds because we wanted to get it right before asking anyone else to look. Our low profile is not modesty. It is a deliberate decision. We are developing and refining this framework quietly until it has enough independent verification to withstand the scrutiny that visibility would bring. We understand the difficulty and professional risk that may come with engaging with an independent, AI-assisted research program. We hope that our work provides a useful contribution as we look to move beyond physics to the other areas the framework was built to examine. The author can be reached at [mcgintyl@grinnell.edu](mailto:mcgintyl@grinnell.edu).
