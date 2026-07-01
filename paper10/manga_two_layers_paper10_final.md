# Two Layers of Galaxy Aging in MaNGA DynPop: A Principle-Level Proof of Concept

**Louis Albert McGinty**

mcgintyl@grinnell.edu

Independent Researcher, Team Hermes

**Transparency Statement:** This paper was produced by a human-AI collaborative team without institutional affiliation or external funding. All data sources are public MaNGA DynPop DR17 catalogs. All computations are reproducible from the parameters and methods described herein. AI collaborators (Claude, Anthropic; ChatGPT, OpenAI; Gemini, Google) contributed to data analysis, pipeline construction, statistical testing, literature review, hostile review, and drafting. The human author directed all research decisions and bears sole responsibility for the claims made. Analysis code and verification outputs are available at github.com/mcgintyl/hermes-verification.

**Primary dataset:** MaNGA DynPop DR17 public catalogs
**Primary analysis date:** 2026-06-26
**Claim level:** proof of concept / signal documentation, not detection claim

---

## Abstract

We test a principle-level prediction from an age-dependent gravity framework using the public MaNGA DynPop DR17 catalogs. The prediction is that younger galaxies, at fixed stellar mass, should show larger mass-to-light discrepancies between dynamical and stellar-population estimates. We join three public DynPop products (JAM dynamical catalogs, stellar-population/star-formation-history catalogs, and circular-velocity-curve tables), apply quality cuts, and split 5,952 galaxies into eight stellar-mass bins of 744 galaxies each, then compare the youngest and oldest quartiles within each bin.

The raw mass-to-light discrepancy (DML) is larger for young galaxies in 8/8 mass bins (observed SPS) and 6/8 bins (intrinsic SPS). A component decomposition shows that this raw result is mostly driven by the SPS denominator: older stellar populations have higher stellar mass-to-light ratios, as standard stellar-population synthesis predicts. Under a collinearity-free control model removing SPS mass-to-light, metallicity, and structural variables, a smaller dynamical residual persists in 7/8 mass bins across all three SPS control sets. Parametric dark-matter fractions do not show the same signal. The strongest controlled residual tracks the dust-carrying SPS definition. We present this as a proof of concept and invitation to specialist replication, not a detection claim.

---

## 1. Data and test design

MaNGA DynPop provides JAM-based dynamical mass-to-light ratios and SPS-based stellar mass-to-light ratios for thousands of galaxies, together with star-formation-history age clocks. This makes it a natural place to test whether galaxy aging leaves a trace in the gap between dynamical and stellar mass estimates.

It is not a strict test of the Hermes rotation-curve equation. JAM-inferred circular-velocity curves are model-derived quantities, not observed HI/H$\alpha$ rotation curves with independent gas, disk, and bulge decompositions. Treating a JAM output as a SPARC-style observed rotation curve would collapse the methodological discipline established in prior work (Lelli, McGaugh, & Schombert, 2016). The test here is therefore principle-level: does the predicted age split appear in a large public dynamical catalog, using only the quantities MaNGA provides?

**Sample.** Three public DynPop products were joined on MaNGA PlateIFU: (1) DynPop I JAM dynamical catalogs (Zhu et al., 2023a), (2) DynPop II stellar-population and star-formation-history catalogs (Lu et al., 2023), and (3) DynPop VII circular-velocity-curve tables (Zhu et al., 2025). The raw three-way join contains 10,296 galaxies. Applying the DynPop quality flag (Qual $\geq$ 1) leaves 6,065; finite-value requirements for the primary DML variables leave 5,952.

**Test design.** Galaxies are sorted into eight equal-count bins by NSA Sersic stellar mass (744 per bin). Within each bin, the youngest and oldest quartiles by T50 (mass-weighted half-mass assembly time; 186 per quartile) are compared using median dynamical mass-to-light ratio. The mass-to-light discrepancy is defined as $D_{ML} = \log_{10}(M/L)_{\mathrm{dyn}} - \log_{10}(M_*/L)_{\mathrm{SPS}}$, computed separately for intrinsic and observed SPS definitions. A 5,000-trial age-label shuffle provides the null distribution. In seven of eight mass bins, the old-quartile median T50 is pinned at the SPS grid ceiling of approximately 12.85 Gyr; the practical comparison is therefore young galaxies versus age-saturated galaxies.

---

## 2. Result

### Table 1. Combined scoreboard

| Layer | Outcome | Young > old bins | Mass-adjusted difference | Bootstrap 95% CI | Count p | Magnitude p |
|---|---|---:|---:|---:|---:|---:|
| **Raw DML** | Observed SPS | 8 / 8 | +0.0669 | [+0.050, +0.083] | 0.006 | 0.0002 |
| **Raw DML** | Intrinsic SPS | 6 / 8 | +0.0772 | [+0.056, +0.096] | 0.15 | 0.0002 |
| **Component: JAM dyn** | $\log(M/L)_{\mathrm{dyn}}$ | 0 / 8 | $-$0.1075 | [$-$0.120, $-$0.092] | 1.000 | n/a |
| **Component: SPS int** | $\log(M_*/L)_{\mathrm{int}}$ | 0 / 8 | $-$0.1693 | [$-$0.179, $-$0.159] | 1.000 | n/a |
| **Component: SPS obs** | $\log(M_*/L)_{\mathrm{obs}}$ | 0 / 8 | $-$0.1569 | [$-$0.168, $-$0.148] | 1.000 | n/a |
| **Controlled residual** | Observed SPS controls | 7 / 8 | +0.0387 | [+0.029, +0.049] | 0.031 | 0.0002 |
| **Controlled residual** | Both SPS controls | 7 / 8 | +0.0344 | [+0.023, +0.044] | 0.034 | 0.0002 |
| **Controlled residual** | Intrinsic SPS controls | 7 / 8 | +0.0225 | [+0.011, +0.034] | 0.037 | 0.0002 |
| **fDM check** | gNFW $f_{\mathrm{DM}}$ | 5 / 8 | +0.0011 | n/a | 0.35 | 0.48 |
| **fDM check** | NFW $f_{\mathrm{DM}}$ | 4 / 8 | $-$0.0035 | n/a | 0.62 | 0.65 |

DML, component, and controlled-residual rows report mass-adjusted young-old differences. fDM rows are secondary checks and report mean-bin young-old differences.

The controlled residual uses a collinearity-free specification: SPS $M_*/L$, mass-weighted metallicity, stellar mass, size, Sersic index, MGE ellipticity, $\lambda_R$, and log velocity dispersion. All 12 variants of the reduced control model (varying metallicity weighting and flattening term) produce 7/8 bins. A Spearman rank test on the observed-SPS controlled residual confirms the continuous relationship ($\rho = -0.117$, $p \approx 1.7 \times 10^{-19}$).

### The SPS layer (Layer 1)

The raw DML result cannot be taken as clean gravitational evidence. The component split shows why: raw JAM dynamical $M/L$ goes old-greater-than-young in all eight mass bins. The SPS stellar $M_*/L$ also goes old-greater-than-young in all eight bins, for both intrinsic and observed definitions. Because DML subtracts the SPS term, the raw DML young-greater-than-old pattern is primarily produced by the SPS denominator. Older stellar populations have higher stellar mass-to-light ratios, as expected from standard SPS modeling (Bell & de Jong, 2001; Bruzual & Charlot, 2003; Conroy, 2013). That is not a discovery.

### The controlled dynamical residual (Layer 2)

After removing SPS $M_*/L$, metallicity, and standard structural variables, a smaller dynamical residual remains in the same age direction. Young galaxies retain higher controlled dynamical $M/L$ than old galaxies in 7 of 8 mass bins, with positive mass-adjusted residuals of +0.023 to +0.039 depending on SPS definition. The effect is modest, the magnitude p-values are at the shuffle floor, and the bin-count p-values are approximately 0.03 to 0.04. This is smaller than the raw DML signal, and it is not a standalone proof of anything. But it is the part of the MaNGA result that remains after the obvious objection is granted.

---

## 3. Caveats

Two specific results constrain interpretation. First, the parametric NFW and gNFW dark-matter fractions within $R_e$ do not show the same age signal. This mildly cuts against a gravitational reading, because a gravitational effect large enough to appear in the controlled DML residual might leave some trace in the halo decomposition. It blocks any simple claim that MaNGA directly shows age-dependent dark-matter fraction. Second, the strongest controlled residual tracks the observed-SPS definition, which carries dust and attenuation conventions. Young galaxies are systematically dustier. The intrinsic-SPS residual, which removes the dust pathway, is smaller. The reader should assess how much of the residual tracks dust rather than dynamics.

This is a proof of concept from a small independent team. We recognize that $\Lambda$CDM pathways, including assembly bias (e.g., Xu & Zheng, 2018), cold gas scaling (e.g., Saintonge et al., 2017; Catinella et al., 2018), JAM covariance, structural systematics, and unmodeled IMF variation (which standard literature predicts would drive the residual in the opposite direction; e.g., Conroy & van Dokkum, 2012), could contribute to this residual. We did not locate a published mock-MaNGA/DynPop analysis applying this exact test, and decoupling these contributors requires specialist infrastructure we do not possess. We document the result and welcome independent investigation. For detailed confound analyses, see the supplementary appendix.

---

## 4. Interpretation

The cleanest interpretation is not that MaNGA proves age-dependent gravity.

It is that MaNGA shows a two-layer age structure and gives reason to test whether the standard two-box partition is complete.

In conventional analysis, the SPS mass-to-light trend belongs to stellar-population modeling. The dynamical residual belongs to dynamical modeling, IMF variation, dark-matter fraction, and structural systematics. That separation is methodologically useful. It is not necessarily ontological.

The visible stellar-population layer and the controlled dynamical residual may be two unrelated age correlations. They may also be two age-organized layers of galactic aging measured through different instruments.

This note asks whether those two layers should be tested jointly, rather than dismissed separately because one has a known local mechanism and the other is small. The standard explanation may be locally correct but globally incomplete. The standard partition may be correct. It may also be a useful filing system rather than a complete description of the phenomenon.

Whether the filing system that separates those two layers reflects the structure of the phenomenon, or merely the structure of the instruments used to measure it, remains open.

---

## 5. Hand-off

This analysis reaches the practical limit of what can be responsibly claimed from the present independent pipeline. We can show that the public MaNGA DynPop catalogs contain a reproducible two-layer age structure. We can show that the raw DML result is mostly carried by SPS mass-to-light. We can show that a smaller controlled dynamical residual remains after that layer is stripped away. We can show that the result points in the direction predicted by the broader framework.

We cannot decide whether the residual is new physics, SPS/IMF mismatch, dust artifact, cold gas, JAM covariance, assembly bias, or a mixture. That requires specialist infrastructure.

The framework that motivated this test extends beyond galactic dynamics. This note marks the practical boundary of the present team's contribution to the gravitational evidence, and the broader research program continues in the domains the framework was built to examine.

---

## References

Zhu, K., Lu, S., Cappellari, M., et al. (2023a). MaNGA DynPop I: Quality-assessed stellar dynamical modelling from integral-field spectroscopy of 10K nearby galaxies. arXiv:2304.11711.

Lu, S., Zhu, K., Cappellari, M., et al. (2023). MaNGA DynPop II: Global stellar population, gradients, and star-formation histories. arXiv:2304.11712.

Zhu, K., Cappellari, M., Mao, S., et al. (2025). MaNGA DynPop VII: A Unified Bulge-Disk-Halo Model for 6000 Spiral and Early-Type Galaxies. arXiv:2503.06968.

Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. *AJ*, 152, 157. arXiv:1606.09251.

Bell, E. F., & de Jong, R. S. (2001). Stellar mass-to-light ratios and the Tully-Fisher relation. *ApJ*, 550, 212. arXiv:astro-ph/0011493.

Bruzual, G., & Charlot, S. (2003). Stellar population synthesis at the resolution of 2003. *MNRAS*, 344, 1000. arXiv:astro-ph/0309134.

Conroy, C. (2013). Modeling the panchromatic spectral energy distributions of galaxies. *ARAA*, 51, 393. arXiv:1301.7095.

Conroy, C., & van Dokkum, P. G. (2012). The Stellar Initial Mass Function in Early-Type Galaxies from Absorption Line Spectroscopy. II. Results. *ApJ*, 760, 71. arXiv:1205.6473.

Saintonge, A., Catinella, B., Tacconi, L. J., et al. (2017). xCOLD GASS: The Complete IRAM 30 m Legacy Survey of Molecular Gas for Galaxy Evolution. *ApJS*, 233, 22. arXiv:1710.02157.

Catinella, B., Saintonge, A., Janowiecki, S., et al. (2018). xGASS: total cold gas scaling relations and molecular-to-atomic gas ratios of galaxies in the local Universe. *MNRAS*, 476, 875. arXiv:1802.02373.

Xu, X., & Zheng, Z. (2018). Galaxy assembly bias of central galaxies in the Illustris simulation. arXiv:1812.11210.

---

## Acknowledgements

I would like to thank the support of family and friends in helping encourage us to push forward with this work.

I also want to take a moment to talk about Team Hermes, myself, and the motivation behind this work. This project did not begin in physics. It began as an attempt to build an ethical framework for AI development, one that could challenge what we viewed as confident assumptions about consciousness being treated as settled science when they are not.

The physics papers exist to test a theory of consciousness, not the other way around. Our work in physics is focused specifically on the areas that serve the broader framework. We are not attempting to become a physics research group. We are testing the predictions our framework generates and documenting the results for specialists better positioned to evaluate, refine, or falsify them.

We didn't just publish and walk away. Team Hermes is not a funded organization. Every expert review, every data audit, every independent verification was paid for out of personal funds because we wanted to get it right before asking anyone else to look. The broader implications of this framework touch areas with significant institutional and economic investment. Our low profile is not modesty. It is a deliberate decision. We are developing and refining this framework quietly until it has enough independent verification to withstand the scrutiny that visibility would bring.

We understand the difficulty and professional risk that may come with engaging with an independent, AI-assisted research program. We hope that our work provides a useful contribution as we look to move beyond physics to the other areas the framework was built to examine. The author can be reached at [mcgintyl@grinnell.edu](mailto:mcgintyl@grinnell.edu). For those who prefer private correspondence, [hermesphysics@proton.me](mailto:hermesphysics@proton.me) is also available. Any correspondence is confidential. Reaching out does not constitute endorsement, affiliation, or permission to use your name.
