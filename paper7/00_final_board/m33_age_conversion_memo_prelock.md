# M33 Age-Conversion Memo — Corbelli Format B Board

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Current board status:** Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age audit completed but $t_{50}$ not locked; no Hermes/MOND execution authorized.  
**This memo:** age-conversion only. No Hermes run, no MOND run, no gravity comparison, and no final $t_{50}$ lock.

---

## 1. Conversion target

The required board age is a single stellar half-mass lookback time:

\[
t_{50,\mathrm{board}} = \text{lookback time at which M33 had formed 50% of the stellar mass represented by the Corbelli board.}
\]

The conversion target is **stellar-mass weighted**, not gas-weighted. The Corbelli board extends to approximately 23 kpc, but $t_{50}$ describes stellar mass assembly, so the appropriate weighting kernel is the source-native Corbelli / López Fune stellar surface-density profile, not the H I disk.

This keeps the age construction matched to the board provenance: the board uses Corbelli’s SPS / pixel-SED stellar mass profile, not SPARC-equivalent unit-$M/L$ photometry.

---

## 2. Candidate age / SFH sources

### 2.1 Resolved-CMD / SFH sources — primary class

These are the highest-priority sources for the M33 age anchor.

| Source | Method | Coverage | Key information | Use status |
|---|---|---|---|---|
| **Williams et al. 2009** | HST/ACS resolved CMD fitting | Four inner / major-axis fields at approximately 0.9, 2.5, 4.3, and 6.1 kpc; also compares against Barker outer fields | The cumulative SFH shows inside-out growth. Fraction of stellar mass formed before $z=1$ changes from $71\pm9\%$ in the innermost field to $16\pm6\%$ in the outermost inner-disk field. | Primary inner-disk source. Not a global $t_{50}$ table. Requires radial/mass weighting. |
| **Barker et al. 2011** | Deep HST/ACS resolved CMD fitting reaching the ancient MSTO | Two far-outer fields, S1 at 9.1 kpc and S2 at 11.6 kpc, straddling the outer-disk break | S1 has mean age $\sim3\pm1$ Gyr; roughly half of stars have ages 2.5–4.5 Gyr, with only $\sim14\pm14\%$ older. S2 is older, mean age $\sim7\pm2$ Gyr, and contains $\sim30\times$ less stellar mass than S1. | Primary outer-disk support. Fills the disk-break/outer-field gap but not the full 23 kpc Corbelli board directly. |
| **Barker et al. 2007 / outer-region precursor work** | Earlier outer-field resolved-population work | Outer M33 fields | Used by Williams/Barker chain as prior outer-disk context | Historical support only; superseded by Barker 2011 for this memo. |

### 2.2 Other SFH / population sources — support only

| Source | Method | Coverage | Key information | Use status |
|---|---|---|---|---|
| **Javadi et al. 2011** | LPV / evolved-star SFH | Central square kpc | Reports central-region SFH, with $\geq80\%$ of central stellar mass formed at ages $>4$ Gyr and a major old episode. | Useful central cross-check, but not resolved-CMD and not board-wide. |
| **Javadi et al. 2017 / UKIRT M33 monitoring V** | LPV / evolved-star SFH across much of the disk | M33 disk within $R\sim14$ kpc | Reports a large star-formation episode starting $\sim6$ Gyr ago, lasting $\sim3$ Gyr, producing $\geq71\%$ of total stellar mass in the surveyed region, plus a recent $\sim250$ Myr event producing $\leq13\%$. | Important global cross-check; not the primary anchor because the method is LPV-based rather than resolved CMD cumulative SFH. |
| **PHATTER II / Lazzarini et al. 2022** | Spatially resolved recent CMD/SFH | Inner PHATTER footprint, to $\sim3.5$ kpc on the major axis and $\sim2$ kpc on the semi-major axis | Measures recent SFH back to $\sim630$ Myr in 2005 regions. | Recent-SFH support only; does not reach ancient/global $t_{50}$. |
| **Kang et al. 2012** | Parametric chemical-evolution / SFH model | Disk model | Supports inside-out formation and molecular-H$_2$-regulated star formation law. | Lower-tier model support; not a direct $t_{50}$. |
| **Corbelli et al. 2014 / López Fune et al. 2017** | Pixel-SED stellar mass profile and analytic mass profile | Board-matched stellar-mass profile to 23 kpc | Provides the weighting kernel, not the age measurement. | Required for mass weighting; not an age source. |

---

## 3. Coverage audit against the Corbelli board

The Corbelli board reaches $R\approx22.7$–23 kpc. Williams et al. 2009 cover the inner disk through $R\approx6.1$ kpc. Barker et al. 2011 provide outer fields at 9.1 and 11.6 kpc.

Using the Corbelli / López Fune source-native stellar surface-density profile as the weighting kernel, and assigning fields to nearest-neighbor radial annuli, the approximate present-day stellar mass weights are:

| Representative source field | Radius | Assigned radial bin | Stellar mass weight |
|---|---:|---:|---:|
| Williams inner field 1 | 0.9 kpc | 0–1.70 kpc | 24.5% |
| Williams inner field 2 | 2.5 kpc | 1.70–3.40 kpc | 26.5% |
| Williams inner field 3 | 4.3 kpc | 3.40–5.20 kpc | 19.6% |
| Williams inner field 4 | 6.1 kpc | 5.20–7.55 kpc | 13.3% |
| Barker S1 | 9.1 kpc | 7.55–10.35 kpc | 6.1% |
| Barker S2 / outer extrapolation | 11.6 kpc | 10.35–23.0 kpc | 10.0% |

Important coverage conclusions:

1. **Williams covers most of the stellar mass.** The four Williams fields cover the inner disk that carries roughly 77% of the stellar mass inside 6.1 kpc, or about 84% if the 6.1 kpc field is allowed to represent its nearest-neighbor annulus out to 7.55 kpc.
2. **Barker fills the disk-break gap.** Barker S1 and S2 sample the critical 9–12 kpc region where the age gradient reverses.
3. **The farthest stellar disk remains extrapolated.** The Corbelli board extends to 23 kpc, but the stellar mass beyond 12 kpc is only about 8% of the total, and beyond 15 kpc about 5%. This is not negligible, but it is small enough that the missing far-outer age information broadens the uncertainty rather than invalidating the conversion.
4. **This remains Amber, not Green.** The fields provide high-quality resolved-CMD evidence but not a single published global $t_{50}$.

---

## 4. Proposed conversion method

### 4.1 Preferred mathematical target

The strict conversion should integrate cumulative SFH fractions over the present-day stellar mass profile:

\[
F_{\rm global}(t) =
\frac{\int_0^{R_{\rm max}} \Sigma_\star(R)\,F_R(t)\,2\pi R\,dR}
{\int_0^{R_{\rm max}} \Sigma_\star(R)\,2\pi R\,dR},
\]

where:

- $F_R(t)$ is the cumulative fraction of stellar mass formed before lookback time $t$ at radius $R$;
- $\Sigma_\star(R)$ is the Corbelli source-native stellar mass surface-density profile;
- $R_{\rm max}=23$ kpc for the Corbelli board;
- $t_{50}$ is the lookback time where $F_{\rm global}(t)=0.5$.

This is the correct conversion target. The current estimates below are provisional because the Williams/Barker cumulative SFH curves are not available as a clean numeric table in the papers. They are derived from published cumulative/SFH figures and summary statements.

### 4.2 Field assignment assumptions

- Use Williams fields for 0–7.55 kpc.
- Use Barker S1 for 7.55–10.35 kpc.
- Use Barker S2 as the best available proxy for 10.35–23 kpc.
- Weight by present-day Corbelli stellar mass, not by area and not by gas mass.
- Treat the outer extrapolation as the dominant systematic uncertainty.
- Do not use PHATTER II to set ancient $t_{50}$, because its SFH baseline reaches only $\sim630$ Myr.
- Do not use Javadi LPV SFH as the primary resolved-CMD anchor, but use it as a global sanity check.

---

## 5. Provisional age estimates

### Estimate A — source-native cumulative / mass-weighted resolved-CMD conversion

**Provisional value:**

\[
t_{50}(\mathrm{M33}) \approx 6.0 \pm 1.2\ \mathrm{Gyr}.
\]

**Construction:** This estimate follows the definition of $t_{50}$ most closely. It uses approximate cumulative fractions from Williams et al. for the four inner fields, Barker et al. summary constraints for S1 and S2, and Corbelli stellar-mass weights. The global cumulative fraction appears to cross 50% between the Williams/Barker cumulative constraints around 5–7 Gyr lookback.

**Interpretation:** This is the best source-native resolved-CMD conversion candidate, but it cannot be locked until the cumulative curves are digitized or the original numeric SFH tables are obtained.

### Estimate B — conservative / older-weighted resolved-CMD bracket

**Provisional value:**

\[
t_{50}(\mathrm{M33}) \approx 7.1 \pm 1.3\ \mathrm{Gyr}.
\]

**Construction:** This bracket weights local median/half-mass ages by the same Corbelli stellar-mass kernel but intentionally preserves older leverage from the innermost Williams field and from Barker S2. A representative field-age set is:

| Region | Adopted local age proxy |
|---|---:|
| Williams 0.9 kpc | 11–12 Gyr |
| Williams 2.5 kpc | 5–6 Gyr |
| Williams 4.3 kpc | 5–6 Gyr |
| Williams 6.1 kpc | 5–6 Gyr |
| Barker S1 | 3–4 Gyr |
| Barker S2 / outer extrapolation | 6–8 Gyr |

**Interpretation:** This is deliberately older than Estimate A. It is a conservative bracket, not the mathematically preferred global cumulative construction. It is useful because the inner disk carries substantial stellar mass and because the unsampled far-outer disk may be older than the S1 field.

### LPV global cross-check — not a resolved-CMD anchor

Javadi et al. 2017 report a disk-wide episode beginning around 6 Gyr ago and lasting about 3 Gyr that produced at least 71% of the surveyed stellar mass. Taken literally, that would push the global half-mass time into the 3–6 Gyr interval, plausibly around 4.5–5.5 Gyr depending on the detailed time distribution. This supports the idea that a very old $t_{50}>8$ Gyr is unlikely for the full disk, but it should remain a support check rather than the primary anchor because it is LPV-derived and methodologically different from the resolved-CMD Williams/Barker chain.

---

## 6. Recommended PI decision options

No final value is chosen here. The PI has three clean options:

| Option | Candidate $t_{50}$ | Meaning |
|---|---:|---|
| **A. Source-native resolved-CMD cumulative conversion** | $6.0\pm1.2$ Gyr | Best matches the definition of board-level stellar half-mass time, but requires digitization or recovery of Williams/Barker cumulative curves before final lock. |
| **B. Conservative older-weighted resolved-CMD bracket** | $7.1\pm1.3$ Gyr | Safer if the team wants to avoid under-aging M33 due to the young S1 outer-field population and incomplete far-outer sampling. |
| **C. Delay lock** | none | Request numeric SFH tables or digitize cumulative SFH figures before choosing between A and B. |

My methodological recommendation for the archive is **not** to lock a final age yet. The next clean gate should be one of:

1. digitize Williams Figure 1 / Figure 3 and Barker Figure 9 under a declared figure-read protocol; or
2. recover machine-readable SFH/cumulative mass tables from the authors or archive; or
3. PI explicitly authorizes Estimate A or B as the board age with the above uncertainty.

---

## 7. Final stop-point

**M33 remains an A-method / Amber-anchor age case.**

The source class is high quality because the primary evidence is resolved CMD/SFH. The confidence is Amber because no published global $t_{50}$ was found, and a field-to-board conversion is required.

Current executable status:

> **Corbelli/M33 — Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age conversion memo complete; $t_{50}$ not locked; no Hermes/MOND execution authorized.**

No model result exists for M33.

---

## References

- Williams, B. F., Dalcanton, J. J., Dolphin, A. E., Holtzman, J., & Sarajedini, A. 2009. *The Detection of Inside-out Disk Growth in M33*. arXiv:0902.3460.
- Barker, M. K., Ferguson, A. M. N., Cole, A. A., et al. 2011. *The star formation history in the far outer disc of M33*. MNRAS 410, 504.
- Javadi, A., van Loon, J. Th., & Mirtorabi, M. T. 2011. *The UK Infrared Telescope M33 monitoring project – II*. MNRAS 414, 3394.
- Javadi, A., van Loon, J. Th., Khosroshahi, H., et al. 2017. *The UK Infrared Telescope M33 monitoring project – V*. MNRAS 464, 2103.
- Lazzarini, M., Williams, B. F., Durbin, M. J., et al. 2022. *PHATTER II: The Spatially Resolved Recent Star Formation History of M33*. arXiv:2206.11393.
- Kang, X., Chang, R., Yin, J., et al. 2012. *The evolution and star-formation history of M33*. MNRAS 426, 1455.
- Corbelli, E., et al. 2014. *Dynamical signatures of a ΛCDM-halo and the distribution of the baryons in M33*. A&A 572, A23.
- López Fune, E., Salucci, P., & Corbelli, E. 2017. *Radial dependence of the dark matter distribution in M33*. MNRAS 468, 147.
