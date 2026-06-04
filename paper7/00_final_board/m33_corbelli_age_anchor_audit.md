# M33 Age-Anchor Audit — Source Review Only

**Galaxy:** M33 / NGC 598 / UGC 1117  
**Lane:** M33-Corbelli source-native thick-disk reconstruction  
**Action:** age-anchor audit only  
**Stop condition:** no Hermes run, no MOND run, no gravity comparison

## Executive finding

M33 has excellent high-tier resolved-CMD star-formation-history sources, but I did **not** find a directly published, galaxy-integrated `t50` value ready to drop into the Hermes chassis.

The highest-quality age route is therefore:

> **Resolved CMD / SFH source-reconstruction:** combine Williams et al. (2009) inner-disk cumulative SFHs with Barker et al. (2011) outer-disk SFHs, weight them by the Corbelli stellar surface-density profile adopted in the board, and solve for the lookback time at which the cumulative formed stellar mass reaches 50%.

Current audit status:

> **Age source quality: Green. Age value lock: not yet.**

A provisional expectation from the literature is that M33’s global `t50` is likely in the **intermediate-age range, roughly 5–7 Gyr**, because the disk grew strongly between 10 and 5 Gyr ago and much of the disk outside the innermost field formed after `z=1`. This bracket is **not locked**. It should be converted into a numeric `t50` only by a declared mass-weighted reconstruction.

## Source hierarchy result

| Source | Method | What it provides | `t50` status | Tier |
|---|---|---|---|---|
| Williams et al. 2009, *The Detection of Inside-Out Disk Growth in M33* | HST/ACS resolved CMD modeling with MATCH | Four inner-disk cumulative SFHs; radial inside-out growth; percentage formed before `z=1`; scale-length evolution with time | Regional cumulative SFH, not a direct global `t50` table | Green source / pending reconstruction |
| Barker et al. 2011, *The star formation history in the far outer disc of M33* | Deep HST/ACS resolved CMD modeling | Two outer-disk fields at 9.1 and 11.6 kpc; mean ages and cumulative SFH behavior across the break | Outer-field constraints; not direct global `t50` | Green source / pending reconstruction |
| Lazzarini et al. 2022, PHATTER II | Spatially resolved recent CMD-based SFH | 2005 regions in the PHATTER footprint; recent SFH back to ~630 Myr | Too recent for lifetime `t50`; useful recent-SFR context only | Supplementary, not `t50` anchor |
| Kang et al. 2012 | Parametric chemical/SFH model | Inside-out formation model constrained by M33 observables | Model-dependent; lower priority than resolved CMD | Amber fallback only |
| Corbelli et al. 2014 | SPS / pixel-SED stellar mass map | Stellar mass surface-density profile for board weighting | Mass-weighting support, not an age source | Support only |

## Highest-tier source: Williams et al. 2009

Williams et al. present HST/ACS resolved photometry of four fields along M33’s major axis and use CMD modeling to derive SFHs. The paper reports that the fraction of stellar mass formed before `z=1` changes from `71 ± 9%` in the innermost field to `16 ± 6%` in the outermost field. It also reports that the disk scale length grew from `1.0 ± 0.1 kpc` about 10 Gyr ago to `1.8 ± 0.1 kpc` at times more recent than 5 Gyr ago.

Age relevance:

- This is the strongest source because it uses resolved CMD/SFH methods.
- It directly shows that M33’s stellar disk grew inside-out.
- It does not provide a single galaxy-integrated `t50` value.
- It provides the inner-disk cumulative SFH curves needed for a mass-weighted `t50` reconstruction.

Source-specific caution:

- The authors note that the high time resolution at 2–12 Gyr is more robust in the outer three major-axis fields than in the shallower innermost field.
- They tested binning, distance, and reddening assumptions and report that the inside-out trend remains.

## Outer-disk source: Barker et al. 2011

Barker et al. extend the resolved-CMD SFH picture to two far-outer fields at 9.1 and 11.6 kpc. The paper finds that the majority of stars in the two fields combined formed at `z < 1`. The inner outer-disk field S1 has mean age approximately `3 ± 1 Gyr`, with about half its stars in the 2.5–4.5 Gyr age range and only about `14 ± 14%` older than 4.5 Gyr. The outer field S2 is older, with mean age approximately `7 ± 2 Gyr`, but contains roughly 30 times less stellar mass than S1.

Age relevance:

- This is also high-tier resolved CMD/SFH evidence.
- It is essential for not over-weighting the inner disk when reconstructing a global `t50`.
- Because S2 is much lower mass, it should affect the global mass-weighted `t50` weakly despite being older.

## Recent-SFH source: PHATTER II / Lazzarini et al. 2022

PHATTER II measures M33’s spatially resolved recent SFH over the HST PHATTER footprint using CMD fitting in 2005 regions, reaching back to about 630 Myr. It reports a mean SFR over the last 100 Myr and a current SFR over the last 10 Myr.

Age relevance:

- Strong recent-SFH source.
- Not suitable for lifetime `t50`, because the time baseline stops at ~630 Myr.
- Useful as a check on present/recent SFR, not as the Hermes age anchor.

## Recommended age-construction recipe

No `t50` should be locked until this reconstruction is performed:

1. Use Williams et al. 2009 cumulative SFH curves for the inner radial fields.
2. Use Barker et al. 2011 cumulative SFH curves for the outer fields.
3. Interpolate cumulative SFH as a function of radius.
4. Weight each radial bin by the same Corbelli `Sigma_star(R)` profile used in the board.
5. Integrate cumulative formed stellar mass over radius.
6. Solve for the lookback time where cumulative formed mass reaches 50%.
7. Propagate a bracket from the Williams/Barker uncertainties and from the coarse old-age binning.

Preferred label after reconstruction:

> **Resolved CMD / SFH mass-weighted reconstruction — Green source, Amber numerical lock if figure digitization is required; Green numerical lock only if source tables or author data are obtained.**

## Provisional bracket, not locked

The literature supports an intermediate-age global half-mass time, but not a precise number yet.

Evidence:

- Williams et al. find most disk stars outside the innermost field formed after `z=1`.
- The disk scale length grew substantially from 10 Gyr ago to the recent past, with most growth occurring between 10 and 5 Gyr ago.
- Barker S1, which is closer to the high-surface-density outer edge than S2, has mean age ~3 Gyr and about half its stars in the 2.5–4.5 Gyr bin.
- Barker S2 is older but much less massive.

Audit bracket:

> **Provisional `t50` expectation: ~5–7 Gyr. Not locked.**

This bracket should not be used in Hermes until PI approval of a formal age-construction pass.

## Age-audit verdict

| Gate | Status |
|---|---|
| Direct published global `t50` found | No |
| Resolved CMD/SFH sources found | Yes |
| Source quality | Green |
| Numerical age lock | Pending reconstruction |
| Recommended fallback | Do not use lower-tier fallback unless resolved-CMD reconstruction fails |
| Current model authorization | None |

Final age status:

> **M33 age anchoring is promising and high-tier, but not locked. The next age step is a mass-weighted resolved-SFH reconstruction from Williams et al. 2009 + Barker et al. 2011 using the Corbelli stellar-mass profile.**

## References

- Williams, B. F., Dalcanton, J. J., Dolphin, A. E., Holtzman, J., & Sarajedini, A. 2009, *The Detection of Inside-Out Disk Growth in M33*, ApJ Letters, 695, L15–L19.
- Barker, M. K., Ferguson, A. M. N., Cole, A. A., et al. 2011, *The star formation history in the far outer disc of M33*, MNRAS, 410, 504–516.
- Lazzarini, M., Williams, B. F., Durbin, M. J., et al. 2022, *PHATTER II: The Spatially Resolved Recent Star Formation History of M33*, ApJ, 934, 76.
- Kang, X., Chang, R., Yin, J., et al. 2012, *The evolution and star-formation history of M33*, MNRAS, 426, 1455–1464.
- Corbelli, E., Thilker, D., Zibetti, S., Giovanardi, C., & Salucci, P. 2014, A&A, 572, A23.
