# A Wear-Activated Density-Release Refinement of the Hermes Gravity Equation: A Paper 1 Addendum

**Louis Albert McGinty**, mcgintyl@grinnell.edu, Independent Researcher, Team Hermes

AI-assisted research; full disclosure in Section 6.

**Version note:** Version 3 recomputes the M33 results with the gate as corrected in Paper 7 v3; earlier versions used the superseded gate, so the outer-disk $\chi^2_\nu$ now falls from 5.86 to 1.30–1.53 rather than from 5.79 to 1.21. It also withdraws the statements that the operator was forced by constrained elimination and that its derived constants were uniquely selected, treats M33 as a motivating case rather than an independent test, adds a comparison at the more usual stellar mass-to-light ratios (where MOND fits better), corrects the summary of the M33 literature, and deposits the production data, with a pinned commit of the code that rebuilds them. The SPARC results at unit stellar mass-to-light ratio are unchanged.

---

## Abstract

The Hermes equation (McGinty 2026a) predicts galaxy rotation curves from the baryonic acceleration profile, stellar age and peak baryonic acceleration with zero adjustable parameters per galaxy. In M33 its residual concentrates in the extended, low-density outer disk (McGinty 2026b). This addendum presents a wear-activated density-release refinement that acts only where the baryonic acceleration is low and keeps every Paper 1 constant. Five candidate operators were tested against eight constraints, including no per-galaxy parameters and no material SPARC casualties. Each failed for a specific structural reason and shaped the surviving operator, but the constraints do not single it out: its two new constants ($a_u = a_k/\pi$ and $2\pi^2$) are closed-form expressions of Paper 1 constants, chosen by trial on the same 133 galaxies. Applied to 133 SPARC galaxies at unit stellar mass-to-light ratio, the refinement improves the median reduced $\chi^2_\nu$ from 1.323 to 1.189 (MOND 1.143) with zero material casualties, reduces the share of galaxies MOND fits better from 48.9% to 42.1%, and keeps the number above $\chi^2_\nu = 5$ at 21 (versus 38 for MOND); at the more usual ratios of 0.5 for disks and 0.7 for bulges, MOND fits better. Combined with a literature-motivated younger effective age for M33's outer disk (1 to 2 Gyr), the outer-disk $\chi^2_\nu$ improves from 5.86 to 1.30–1.53; MOND reaches 0.23 on the same board with no age input. The production CSV files are deposited with this record, and the code that rebuilds them is public (Supplement F).

---

## 1. Motivation

The Hermes equation (Paper 1; McGinty 2026a) predicts galaxy rotation curves from the observed baryonic acceleration profile, with its retained-coherence term fixed by the wear exponent $\psi \propto t_{50}\,g_{98}$ (stellar age times peak baryonic acceleration). Applied to 133 SPARC galaxies with verified stellar ages, the equation achieves a median $\chi^2_\nu$ of 1.312 with zero adjustable parameters per galaxy.

Paper 7 (McGinty 2026b) tested the equation on M33 and found the inner disk fits cleanly ($\chi^2_\nu = 0.93$) while the outer disk beyond 10 kpc carries about 85% of the total $\chi^2$ budget ($\chi^2_\nu = 5.86$).

The goal of this addendum is not to replace the Paper 1 equation but to improve its fits to extended, low-density outer disks.


## 2. Refinement Constraints

Any refinement had to satisfy eight constraints, several adopted during the search (Supplement A, Section 2). Five are practical: no per-galaxy tuning, no refitting of Paper 1 constants, no material SPARC casualties (defined as any galaxy degraded by more than 0.5 $\chi^2_\nu$), no baryonic Tully-Fisher relation (BTFR) collision structure (the release's density contrast must not be anchored to the galaxy's peak acceleration $g_{98}$, although $g_{98}$ still enters through $\psi$), and no constant fitted to M33. M33 motivated this search and every candidate was scored on it, so it is a motivating case here, not an independent test.

The other three govern form and reporting. The refinement must act on the positive term $\pi e^{-\psi}$ of the age score rather than on the net score, so that it cannot amplify the floor $-1/\sqrt{2\pi}$; this requirement was added after candidate 2 (Section 3). New constants were restricted to closed-form expressions of existing Hermes constants. All failed candidates are reported (Supplement A).


## 3. Candidates Tested

Five candidates were tested before the operator of Section 4. Each failed at least one constraint.

| # | Candidate | SPARC casualties | SPARC median $\chi^2_\nu$ | M33 outer $\chi^2_\nu$ | Failed constraint |
|---|-----------|:---:|:---:|:---:|---|
| 1 | Scalar $\psi$ smoothing | 53–58 with $\chi^2_\nu > 5$ (baseline 21) | n/a | ~4.14–5.08* | SPARC safety |
| 2 | Sign-protected density-seam | 8 + 3 new catastrophic | n/a | 0.72* | SPARC safety, BTFR collision risk, discontinuity |
| 3 | BTFR-safe positive-term $\eta$ | 1 | 1.115 | 2.59 (uniform 6.0 Gyr) | SPARC safety |
| 4 | Pythagorean geometric projection | 15 + 4 new catastrophic | 1.361 | 4.12 at every outer age (axis erased) | SPARC safety |
| 5 | Expected-value phase transition | 14 + 2 new catastrophic | 1.053 | 1.10 (2 Gyr outer age) | SPARC safety |
| Final | Wear-activated $\eta$ (Section 4) | 0 | 1.189 | 1.53 (2 Gyr outer age) | none |

*Superseded Paper 7 v2 gate; not rerun. All other M33 values use the Paper 7 v3 gate. New catastrophic: newly above $\chi^2_\nu = 5$.

Supplement A gives the candidate-by-candidate analysis and says which features of the surviving operator answer a failure and which do not.


## 4. The Surviving Refinement

The baseline equation is modified by a density-release operator $\eta_\text{WA}(R)$ applied to the positive decay term:

$$
g_\text{model}(R) = g_\text{bar}(R) \left[1 + \phi(R) \left(\pi \, e^{-\psi(R)} \cdot \eta_\text{WA}(R) - \frac{1}{\sqrt{2\pi}}\right)\right]
$$

When $\eta_\text{WA} = 1$ everywhere, this reduces to the Paper 1 equation.

The operator is built from:

**Logarithmic compression:**

$$
\Gamma(g) = a_k \ln\!\left(1 + \frac{\max(g,\, 0)}{a_k}\right)
$$

**Wear activation:**

$$
W(\psi) = 1 - e^{-2\pi^2 \psi}
$$

**Spatial transition:**

$$
S(g) = \left[1 + \exp\!\left(\frac{g - a_k}{0.30 \, a_k}\right)\right]^{-1}
$$

**Universal numerator:**

$$
a_u = \frac{a_k}{\pi} \approx 504.5 \;\text{(km/s)}^2\text{/kpc}
$$

**Local operator:**

$$
\eta_U(R) = \min\!\left[\pi,\; 1 + S\bigl(g_\text{bar}(R)\bigr) \cdot \left(\max\!\left(1,\; \sqrt{\frac{a_u}{\Gamma(g_\text{bar}(R))}}\right) - 1\right)\right]
$$

**Wear-gated operator:**

$$
\eta_\text{WA}(R) = 1 + W(\psi_\text{sys}) \cdot \bigl[\eta_U(R) - 1\bigr]
$$

The wear gate $W$ uses the systemic $\psi_\text{sys}$ from the galaxy's global stellar age. The positive decay term uses $\psi(R)$, which equals $\psi_\text{sys}$ for standard SPARC execution. For the M33 component-age diagnostic, $\psi(R)$ uses a spatially varying effective age (inner disk at global $t_{50}$, outer disk at a younger effective age), while $W$ remains systemic.

**Constants:**

| Constant | Value | Status |
|----------|-------|--------|
| $\pi$, $e$, $c$, $1/\sqrt{2\pi}$, $a_k = 1585$ (km/s)$^2$/kpc | Per Paper 1 | Inherited |
| 0.30 | Paper 1's gate-logistic width (in $R/r_\text{knee}$), reused as a fractional width in $g$ | Reused |
| $a_u = a_k/\pi$ | 504.5 (km/s)$^2$/kpc | Derived |
| $2\pi^2$ | 19.74 | Derived |

The two derived constants are closed-form expressions of inherited constants, and no continuous parameter was fitted to SPARC or M33. The forms were chosen by trial on the 133-galaxy board ($a_u$ with candidate 3, $W$ after it) and are not unique: other closed forms also give zero casualties, for example $a_k/e$ with $2\pi^2$ (median 1.079) or $a_k/\pi$ with $8\pi$ (median 1.154). One galaxy, UGC 05750, sets the zero-casualty limit on both constants. No per-galaxy parameters are introduced on SPARC.


## 5. Results

**Recomputation note.** The baseline median here (1.323) differs from Paper 1's 1.312 because the gate's shear stage uses a log-grid finite difference of $\ln V_{\rm sm}$ here and a chain-rule derivative of $V_{\rm sm}$ in Paper 1; both conventions are described in Paper 1 v6, Section 2, and each reproduces its per-galaxy scores to within $10^{-12}$. On Paper 1's chain-rule convention the same comparison gives a baseline median of 1.312 and a refinement median of 1.150, with zero casualties, 13 improvements and 76 wins over MOND. All three models are recomputed on the same 133 SPARC galaxies (Lelli et al. 2016) at SPARC's unit stellar mass-to-light ratio, as in Paper 1, with the log-grid gate for both Hermes models; $\chi^2_\nu = \chi^2/N$, with $\sigma_\text{int}^2 = 386$ (km/s)$^2$ added to each squared velocity error. MOND (Milgrom 1983) uses the simple interpolating function with $a_0 = 1.2 \times 10^{-10}$ m s$^{-2}$; Paper 1 rounded $a_0$ to 3700 (km/s)$^2$/kpc, hence its MOND median of 1.141.

### 5.1 SPARC All-133

| Model | Median $\chi^2_\nu$ | $\chi^2_\nu > 5$ | $\chi^2_\nu > 10$ | Worsened > 0.5 vs baseline | Improvements > 0.5 | Lower $\chi^2_\nu$ than MOND |
|-------|:---:|:---:|:---:|:---:|:---:|:---:|
| Baseline Hermes | 1.323 | 21 | 12 | 0 | 0 | 68/133 |
| MOND (recomputed) | 1.143 | 38 | 15 | 44 | 44 | n/a |
| Wear-activated $\eta$ | 1.189 | 21 | 12 | 0 | 13 | 77/133 |

Under the refinement the catastrophic count ($\chi^2_\nu > 5$) stays at 21, with no galaxy newly above 5 or 10, against MOND's 38. Taking each galaxy's median absolute velocity residual, the median across galaxies is 17.84 km/s for the refinement, 19.05 for baseline and 20.23 for MOND; pooled over all 3,073 points, as in Paper 1 Table A3, the values are 18.83, 20.67 and 25.82 km/s. All of these values are at SPARC's unit stellar mass-to-light ratio. At the ratios more usually adopted at 3.6 $\mu$m (0.5 for disks, 0.7 for bulges, with $g_{98}$ recomputed), the refinement still improves on the baseline (median 1.90 to 1.20, zero casualties), but MOND fits better (median 0.48, as in Section 3 of Paper 9, McGinty 2026c), and the counts above $\chi^2_\nu = 5$ reverse: 35 for the refinement against 6 for MOND.

Most of the 13 improved galaxies have extended, gas-rich, low-density outer disks where baseline Hermes under-predicted the tail, which is where the operator is built to act; two (NGC 6674, NGC 0289) improve but remain poor fits (Supplement C). In the outermost five points the median signed residual moves from $-20.51$ to $-10.80$ km/s, so the median galaxy is still under-predicted; 16 galaxies flip from under- to over-prediction (none the other way, by construction), and five of them overshoot by more than 5 km/s (Supplement B).

### 5.2 M33 Component-Age Diagnostic

Neither the density-release $\eta$ alone nor a component-resolved age alone closes M33's outer-disk residual. Together they bring the outer $\chi^2_\nu$ from 5.86 to 1.30–1.53.

| Configuration | Outer ($R > 10$ kpc) $\chi^2_\nu$ |
|--------------|:---:|
| Baseline (no $\eta$, uniform age 6.0 Gyr) | 5.86 |
| $\eta$ alone (uniform age 6.0 Gyr) | 2.60 |
| Component age alone (no $\eta$, outer 0 Gyr) | 4.12 |
| $\eta$ + outer age 2.0 Gyr | 1.53 |
| $\eta$ + outer age 1.0 Gyr | 1.30 |

M33 values use the same Paper 1 gate as Section 5.1, as corrected for M33 in Paper 7 v3.

For this one galaxy the outer age is a free input, and the age step at 10 kpc is the same radius used to score the outer disk. The data do not select the age: the outer $\chi^2_\nu$ keeps falling as the outer age approaches zero (1.08 at 0 Gyr), and all 27 outer points stay under-predicted at every age tested. The 1 to 2 Gyr bracket is motivated by published work on M33's outer gas, which allows but does not establish it (Putman et al. 2009; Semczuk et al. 2018; Corbelli and Burkert 2024; Supplement D); as in Paper 7, that literature was read after the residual was seen. It is an effective age for the outer gas layer, not a stellar age: M33's outer stars are old (Barker et al. 2011).

On the same board, MOND reaches an outer $\chi^2_\nu$ of 0.23 (full board 0.12) with no age input (Paper 7); Supplement D gives the full sweep.


## 6. Discussion and Limitations

This addendum reports a constrained refinement candidate and its audit results. We do not adopt the refinement elsewhere in the program: the Paper 1 equation is simpler, and the refinement's gains rest on choices the data do not fix.

Supplement E offers a provisional reading of the operator's form; it is not needed to evaluate the results.

Several limitations apply. The operator was chosen on the same 133 galaxies it is scored on, so the SPARC gains are in-sample, and the constraints do not single it out: candidate 3 has a lower median (1.115) at the cost of one casualty (UGC 05750), and the operator without $S$ has median 1.162 with none. The refinement has not been tested on an independent rotation-curve board: M33 motivated it and was scored throughout (Section 2). Its behavior at cluster scale is not reported: Paper 9 (McGinty 2026c) tests cluster lensing with the unrefined Paper 1 equation. BTFR scatter under the refinement has not been tested with a dedicated analysis. The logistic width ratio (0.30, reused from the Paper 1 gate) has not been re-derived. The M33 component-age result is one galaxy, not a validation of multi-component modeling in general.

The production CSV files are deposited with this record and listed in Supplement F; paper8/reproduce.py in github.com/mcgintyl/hermes-verification (commit d3355f4) rebuilds them from public data, together with every SPARC and M33 result of this version except the candidate 1 and 2 entries.

**AI disclosure.** AI systems (Claude by Anthropic, ChatGPT by OpenAI, Gemini by Google) were used for candidate generation, adversarial review, code/report generation and statistical-audit support, and drafting support. The author directed all research questions, constraints, interpretation, and final editorial decisions.


## References

Barker, M. K. et al. (2011). The star formation history in the far outer disc of M33. Monthly Notices of the Royal Astronomical Society, 410, 504. arXiv:1008.0760.

Corbelli, E. and Burkert, A. (2024). The outskirts of M33: Tidally induced distortions versus signatures of gas accretion. Astronomy & Astrophysics, 685, A38. arXiv:2402.16957.

Lelli, F., McGaugh, S. S. and Schombert, J. M. (2016). SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves. The Astronomical Journal, 152, 157. arXiv:1606.09251.

McGinty, L. A. (2026a). The Hermes Equation: Galaxy Rotation Curves from Age and Density with No Free Constants. Zenodo. Version of 2026-09-01, DOI: 10.5281/zenodo.22237082; all versions (concept): 10.5281/zenodo.18809176.

McGinty, L. A. (2026b). Field Audit: Board-Complete Data and the Limits of External Gravity-Model Testing. Zenodo. Version of 2026-07-13 (v3), DOI: 10.5281/zenodo.21347359; all versions (concept): 10.5281/zenodo.20550091.

McGinty, L. A. (2026c). The Hermes Equation: Cluster Lensing from a Frozen Galaxy Gate. Zenodo. Version of 2026-10-03, DOI: 10.5281/zenodo.23112608; all versions (concept): 10.5281/zenodo.22719297.

Milgrom, M. (1983). A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis. The Astrophysical Journal, 270, 365.

Putman, M. E. et al. (2009). The Disruption and Fueling of M33. The Astrophysical Journal, 703, 1486. arXiv:0812.3093.

Semczuk, M. et al. (2018). Tidally induced morphology of M33 in hydrodynamical simulations of its recent interaction with M31. The Astrophysical Journal, 864, 34. arXiv:1804.04536.


## Supplementary Materials

Supplements A to F accompany this record.
