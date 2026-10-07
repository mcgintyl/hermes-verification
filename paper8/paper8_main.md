# A Wear-Activated Density-Release Refinement of the Hermes Gravity Equation: A Paper 1 Addendum

**Louis Albert McGinty**, mcgintyl@grinnell.edu, Independent Researcher, Team Hermes

AI-assisted research; full disclosure in Section 6.

---

## Abstract

The Hermes equation (McGinty 2026a) predicts galaxy rotation curves from stellar age and peak baryonic acceleration with zero adjustable parameters per galaxy. Two independent diagnostics, an outer-disk residual concentration in M33 (Paper 7) and a density-structured population kink in the wear exponent, point to a coupling limitation in the transition regime between high-density and low-density baryonic environments. This addendum presents a wear-activated density-release refinement that preserves the locked Paper 1 chassis. Five candidate operators were tested against eight constraints including zero per-galaxy parameters, no refitting of inherited constants, and no material SPARC casualties. Each candidate failed for a specific structural reason; the surviving operator was forced into its form by constrained elimination, not theoretical preference. The refinement introduces no newly fitted constants: its two derived quantities ($a_u = a_k/\pi$ and $2\pi^2$) are closed-form expressions of inherited Paper 1 constants. Applied to 133 SPARC galaxies, the refinement improves the median reduced $\chi^2_\nu$ from 1.323 to 1.189 with zero material casualties, reduces the MOND loss rate from 48.9% to 42.1%, and preserves the baseline catastrophic count at 21 (versus 38 for MOND). Combined with a literature-supported component-age model on M33, the outer-disk $\chi^2_\nu$ improves from 5.79 to 1.21. The full candidate rejection path, signed residual audit, M33 component-age sweep, and production data package are provided for independent verification.

---

## 1. Motivation

The Hermes equation (McGinty 2026a) predicts galaxy rotation curves from the observed baryonic acceleration profile, with its retained-coherence term fixed by stellar age and peak baryonic acceleration. Applied to 133 SPARC galaxies with verified stellar ages, the equation achieves a median $\chi^2_\nu$ of 1.312 with zero adjustable parameters per galaxy.

Two independent diagnostics point to a specific limitation. Paper 7 (McGinty 2026b) tested the equation on M33 and found the inner disk fits cleanly ($\chi^2_\nu = 0.80$) while the outer disk beyond 10 kpc carries 86% of the total $\chi^2$ budget ($\chi^2_\nu = 5.79$). Separately, a population-level investigation found that the wear exponent $\psi$ is structured by peak baryonic acceleration quartile rather than smoothly distributed. The kink in the population distribution occurs in the same density regime where M33's per-galaxy residual concentrates.

The baseline equation remains competitive across much of the sample, but these two diagnostics point to a specific transition-regime limitation. This addendum presents and tests a zero-parameter density-transition refinement. The goal is not to replace the Paper 1 equation but to extend its reach into regimes where extended low-density outer disks expose a coupling limitation.


## 2. Refinement Constraints

Any refinement had to satisfy five practical constraints: no per-galaxy tuning, no refitting of Paper 1 constants, no material SPARC casualties (defined as any galaxy degraded by more than 0.5 $\chi^2_\nu$), no BTFR-collision structure (the outer asymptote must not depend on galaxy-specific core density), and M33 retained as a stress test rather than a tuning target.

Additionally, the refinement was required to act on the positive retained-coherence term ($\pi \cdot e^{-\psi}$) rather than on the net age score, to avoid amplifying the coherence floor into unphysical regimes. New constants were restricted to closed-form expressions of existing Hermes constants. All failed candidates were documented because the rejection path constrains the surviving operator at least as strongly as the data constrain the equation.

Full constraint specification is provided in the Supplementary Materials.


## 3. Candidates Tested

Five candidates were tested. Each failed for a specific structural reason. The surviving refinement was not selected first and justified after the fact. It was discovered by constrained elimination.

| # | Candidate | SPARC casualties | M33 outer $\chi^2_\nu$ | Failed constraint | Structural lesson |
|---|-----------|:---:|:---:|---|---|
| 1 | Scalar $\psi$ smoothing | 53-58 above $\chi^2_\nu > 5$ | ~4.14-5.08 | SPARC safety | High-wear branch is load-bearing |
| 2 | Sign-protected density-seam | 8 + 3 new catastrophic | 0.72 | BTFR collision risk, discontinuity | Seam is real; $g_{98}$ coupling is wrong |
| 3 | BTFR-safe positive-term $\eta$ | 1 | 2.50 | SPARC safety | Near-newborn over-activation; needs wear suppression |
| 4 | Pythagorean geometric projection | 15 + 4 new catastrophic | Component-age axis erased | Smooth superposition fails | Regimes must stay distinct |
| 5 | Expected-value phase transition | 14 + 2 new catastrophic | 1.02 (at 2 Gyr outer age) | SPARC safety | Transition-zone permissiveness; must withhold until low-density branch dominates |

Each candidate taught the surviving equation something about its necessary structure. The wear activation (from candidate 3's failure), the universal numerator (from candidate 2's failure), the positive-term application (from candidate 2's discontinuity), and the max-selected density contrast (from candidates 4 and 5's smooth-superposition failures) were all forced into the equation by specific empirical rejections.

Full candidate-by-candidate analysis is provided in the Supplementary Materials.


## 4. The Surviving Refinement

The baseline equation is modified by a density-release operator $\eta_\text{WA}(R)$ applied to the positive decay term:

$$g_\text{model}(R) = g_\text{bar}(R) \left[1 + \phi(R) \left(\pi \, e^{-\psi(R)} \cdot \eta_\text{WA}(R) - \frac{1}{\sqrt{2\pi}}\right)\right]$$

When $\eta_\text{WA} = 1$ everywhere, this reduces to the Paper 1 equation.

The operator is built from:

**Logarithmic compression:**

$$\Gamma(g) = a_k \ln\!\left(1 + \frac{\max(g,\, 0)}{a_k}\right)$$

**Wear activation:**

$$W(\psi) = 1 - e^{-2\pi^2 \psi}$$

**Spatial transition:**

$$S(g) = \left[1 + \exp\!\left(\frac{g - a_k}{0.30 \, a_k}\right)\right]^{-1}$$

**Universal numerator:**

$$a_u = \frac{a_k}{\pi} \approx 504.5 \;\text{(km/s)}^2\text{/kpc}$$

**Local operator:**

$$\eta_U(R) = \min\!\left[\pi,\; 1 + S\bigl(g_\text{bar}(R)\bigr) \cdot \left(\max\!\left(1,\; \sqrt{\frac{a_u}{\Gamma(g_\text{bar}(R))}}\right) - 1\right)\right]$$

**Wear-gated operator:**

$$\eta_\text{WA}(R) = 1 + W(\psi_\text{sys}) \cdot \bigl[\eta_U(R) - 1\bigr]$$

The wear gate $W$ uses the systemic $\psi_\text{sys}$ from the galaxy's global stellar age. The positive decay term uses $\psi(R)$, which equals $\psi_\text{sys}$ for standard SPARC execution. For the M33 component-age diagnostic, $\psi(R)$ uses a spatially varying effective age (inner disk at global $t_{50}$, outer disk at a younger effective age), while $W$ remains systemic.

**Constants:**

| Constant | Value | Status |
|----------|-------|--------|
| $\pi$, $e$, $c$, $1/\sqrt{2\pi}$, $a_k$, 0.30 | Per Paper 1 | Inherited |
| $a_u = a_k/\pi$ | 504.5 (km/s)$^2$/kpc | Derived |
| $2\pi^2$ | 19.74 | Derived |

No constant in the refinement was newly fitted to the 133-galaxy SPARC board or to M33. The two derived constants are closed-form expressions of inherited constants. While no continuous parameters were fitted, the specific algebraic forms ($a_k/\pi$ and $2\pi^2$) were selected because they uniquely satisfied the zero-casualty constraint across all five candidate iterations (Section 3). No per-galaxy parameters are introduced.


## 5. Results

**Recomputation note.** The baseline median reported here (1.323) differs from Paper 1's reported value (1.312) because the two papers evaluate the gate's shear stage with different derivative conventions: Paper 1 used a chain-rule derivative of $V_{\rm sm}$, while this paper uses a log-grid finite difference of $\ln V_{\rm sm}$. Each convention reproduces its corresponding per-galaxy score vector to better than $3 \times 10^{-13}$ absolute (see Paper 1 v6, Methods). All comparisons in this paper use baseline, MOND, and wear-activated $\eta$ recomputed within the same computational environment under the log-grid convention.

### 5.1 SPARC All-133

| Model | Median $\chi^2_\nu$ | $\chi^2_\nu > 5$ | $\chi^2_\nu > 10$ | Worsened > 0.5 vs baseline | Improvements > 0.5 | Lower $\chi^2_\nu$ than MOND |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|
| Baseline Hermes | 1.323 | 21 | 12 | 0 | 0 | 68/133 |
| MOND (recomputed) | 1.143 | 38 | 15 | 44 | 44 | n/a |
| Wear-activated $\eta$ | 1.189 | 21 | 12 | 0 | 13 | 77/133 |

The wear-activated $\eta$ does not optimize the median as aggressively as MOND, but it preserves the baseline board: zero material casualties and no new catastrophic cases. The catastrophic count ($\chi^2_\nu > 5$) remains at 21, compared to MOND's 38. The win rate against MOND improves from 51.1% to 57.9%. The median absolute velocity residual is 17.84 km/s, compared to 19.05 for baseline and 20.23 for MOND.

The 13 improved galaxies share a physically coherent profile: extended, gas-rich, low-density outer disks where baseline Hermes was under-supporting the tail. A signed residual audit confirms the correction reduces outer-disk starvation without systematically overshooting the sample. Full characterization and audit details are provided in the Supplementary Materials.

### 5.2 M33 Component-Age Diagnostic

Neither the density-release $\eta$ alone nor a component-resolved age alone is sufficient to close M33's outer-disk residual. Together they largely resolve it.

| Configuration | Outer ($R > 10$ kpc) $\chi^2_\nu$ |
|--------------|:---:|
| Baseline (no $\eta$, uniform age 6.0 Gyr) | 5.79 |
| $\eta$ alone (uniform age 6.0 Gyr) | 2.50 |
| Component age alone (no $\eta$, outer 0 Gyr) | 4.05 |
| $\eta$ + outer age 2.0 Gyr | 1.44 |
| $\eta$ + outer age 1.0 Gyr | 1.21 |

The required outer-component age of 1 to 2 Gyr is consistent with published gas dynamics work on M33. Putman et al. (2009) place M33's tidal disruption timescale at 1 to 3 Gyr. Semczuk et al. (2018) model a recent M31 interaction at approximately 2 Gyr. Corbelli and Burkert (2024) favor cosmic filament accretion over M31 tidal passage but support a recently refreshed outer gas layer. The component age should be read as the effective coherence or dynamical age of the outer gas layer, not as the stellar half-mass age of the whole outer disk, which Barker et al. (2011) find to be old (mean stellar age approximately 7 Gyr at 11.6 kpc) with an age-gradient reversal beyond the disc break.

The equation's failure under a single-age assumption is consistent with a multi-component interpretation. When tested with a component-resolved effective outer-gas age consistent with published accretion timescales (1 to 2 Gyr), the residual structure largely resolves. The equation does not fit this age; it constrains the diagnostic bracket within which the outer-disk residual closes.

Full M33 radial profiles and literature review are provided in the Supplementary Materials.


## 6. Discussion and Limitations

This addendum does not replace Paper 1. It reports a constrained refinement candidate and its audit results.

The refinement was forced into its surviving form by empirical pressure across 133 galaxies. Each component of the operator traces to a specific candidate failure. The operator structure admits interpretive readings discussed in the Supplementary Materials, but these interpretations remain provisional and are not required to evaluate the empirical results.

Several limitations apply. The refinement has not been tested on independent non-SPARC rotation-curve boards beyond the reconstructed M33 case. The signed residual audit shows a small one-sided positive drift (16/133 sign flips, all under-to-over). BTFR scatter under the refinement has not been tested with a dedicated analysis. The inherited logistic width (0.30) has not been re-derived. The M33 component-age result is one external case, not a validation of multi-component modeling in general.

The full production data package, including per-galaxy results for all three models, M33 component-age sweeps, and the signed residual audit, is provided in the Supplementary Materials for independent verification.

**AI disclosure.** AI systems (Claude by Anthropic, ChatGPT by OpenAI, Gemini by Google) were used for candidate generation, adversarial review, code/report generation and statistical-audit support, and drafting support. The author directed all research questions, constraints, interpretation, and final editorial decisions. Numerical outputs are supplied in reproducible files for independent audit.


## References

Barker, M. K. et al. (2011). The star formation history in the far outer disc of M33. Monthly Notices of the Royal Astronomical Society, 410, 504. arXiv:1008.0760.

Corbelli, E. and Burkert, A. (2024). The outskirts of M33: Tidally induced distortions versus signatures of gas accretion. arXiv:2402.16957.

McGinty, L. A. (2026a). The Hermes Equation: Galaxy Rotation Curves from Age with No Free Constants. Zenodo. DOI: 10.5281/zenodo.19904551.

McGinty, L. A. (2026b). Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing. Zenodo. DOI: 10.5281/zenodo.20550092.

Putman, M. E. et al. (2009). The Disruption and Fueling of M33. arXiv:0812.3093.

Semczuk, M. et al. (2018). Tidally induced morphology of M33 in hydrodynamical simulations of its recent interaction with M31. arXiv:1804.04536.


## Supplementary Materials

Supplementary Materials are provided as an accompanying document containing the candidate rejection path and refinement constraints, signed residual audit, improved-galaxy characterization, M33 component-age sweep and literature notes, dimensional-transition interpretation, and the production data package and CSV manifest.
