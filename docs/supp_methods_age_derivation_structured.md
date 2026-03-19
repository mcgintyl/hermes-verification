# Supplementary Methods — Age Derivation (Hermes / Paper 1)

**Built from:** `galaxy_method_audit (4).md` (151 detailed method entries) and the paper’s 9-class method taxonomy.  
**Purpose:** Auditor-traceable documentation of how each published measurement is converted into **t₅₀** (stellar half-mass time) for the Paper 1 / Hermes dataset.

## Key definitions

- **t₅₀ (Gyr):** lookback time when the cumulative formed stellar mass reached 50% of today’s formed mass (mass-weighted median assembly time).
- **Tier labels:** As recorded in the audit entries (T1/T2/T3) and data-gap/poison classifications where applicable.

## Galaxy index (quick lookup)

A machine-readable index mapping **Galaxy → Method # → Method Class → t₅₀/Tier** was generated from the audit’s per-method tables:

- `galaxy_age_method_index.csv`

(Use this index to jump from a galaxy name to the exact Method entry below.)

---
## Table of contents
- [Resolved CMD / SFH](#resolved-cmd--sfh)
- [IFU / Spectroscopic](#ifu--spectroscopic)
- [Parametric SED fitting](#parametric-sed-fitting)
- [Broadband color → SPS](#broadband-color--sps)
- [UV survey photometry](#uv-survey-photometry)
- [Hα / SFR-based](#ha--sfr-based)
- [sSFR / mass-based](#ssfr--mass-based)
- [IR / radio SFR](#ir--radio-sfr)
- [Compiled / SPARC-derived](#compiled--sparc-derived)

---
## Resolved CMD / SFH

**Conversion:** A resolved SFH provides cumulative formed mass fraction \(F(t)\) versus lookback time \(t\).  
We define \(t_{50}\) as the lookback time where \(F(t_{50})=0.5\). If the source gives binned cumulative fractions, we linearly interpolate between the two bins that bracket 0.5:
\[
t_{50}=t_i+\frac{0.5-F_i}{F_{i+1}-F_i}\,(t_{i+1}-t_i).
\]
If the source directly reports \(\tau_{50}\) or “median age” (50% cumulative mass), we take it directly (unit conversion only).

### Worked example

**Worked example (Method 151 — tabulated τ₅₀):**  
        This method uses a published tabulated value `log(τ₅₀) = 9.71`, which implies:
        \[
        	au_50 = 10^9.71\,\mathrm{yr} = 5.13	imes 10^9\,\mathrm{yr} = 5.13\,\mathrm{Gyr}.
        \]
        The audit records the resulting \(t_50pprox 5.1\,\mathrm{Gyr}\) for the galaxy in this entry.

        ### Method 151: Albers+2019 MNRAS 490, Table 2 — Resolved CMD SFH (Tier A Gold Standard) → tabulated log(τ₅₀) = 9.71 → t₅₀ = 5.1 Gyr [WLM — Isolated Local Group dIrr]
- Data source(s): Albers et al. (2019) MNRAS 490, 5538, Table 2: ACS field log(τ₅₀) = 9.71 → t₅₀ = 5.1 Gyr, log(τ₉₀) = 9.02 → t₉₀ = 1.05 Gyr. UVIS field (outer halo): log(τ₅₀) = 10.01 → t₅₀ = 10.2 Gyr. Deep HST resolved stellar populations reaching below oldest MSTO. Cross-ID: UGCA 444 = WLM = DDO 221. Isolated dIrr at edge of Local Group zero-velocity surface.
- Method: **Tier A Gold Standard.** Tabulated τ₅₀ directly from resolved CMD SFH — no model inversion required. ACS field used as global proxy (authors state it forms ~100× more stellar mass than UVIS field). Mass-weighted correction: t₅₀,global ≈ (0.99 × 5.1) + (0.01 × 10.2) ≈ 5.15 Gyr (ancient halo structurally real but statistically irrelevant). t₉₀ = 1.05 Gyr rules out T1 (not young enough) and T3 (not passive — ~10% formed in last Gyr). Inside-out growth: ancient halo + younger star-forming disk. "Delayed Assembly" dwarf irregular.
- Method validated: Single use (first Albers+2019 application)
- Known issues: Single-field global proxy (justified by 100:1 mass dominance). Radial gradient exists but mass-irrelevant.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGCA 444 (WLM / DDO 221) | 5.1 Gyr | T2 (Active / Intermediate) | HIGH | No | — |

---

## Screen Checks

### Screen 1 — Protocol 18 (The Outshining Veto — τ_main Bracket Test)
**Protocol 18 (revised):** When Smith et al. (2022) reports burst parameters but does NOT tabulate τ_main or t₅₀, you CANNOT reverse-engineer τ_main — not via t₉₀ math, not via morphological priors. Must bracket: compute t₅₀ at τ_main = 0.2 and 10 Gyr. If entire bracket stays in one tier → lock. If it crosses → Data Gap (unless independent Tier C saves it). Exception: f_burst ≥ 0.50 → T1 safe regardless.

**Three Smith et al. galaxies audited (all with 0.10 ≤ f_burst < 0.50):**

| Galaxy | f_burst | Bracket (τ_main 0.2→10) | Tiers crossed | Independent Tier C? | Final status |
|--------|---------|------------------------|---------------|-------------------|--------------|
| NGC 5229 | 0.20 | 13.5→7.4 Gyr | T2/T3 | **YES** — B-V=0.48, Hα, gas 500% → T2 | **T2/~5.0 Gyr (Tier C fallback)** |
| UGC 07866 | 0.43 | 13.2→3.0 Gyr | T1/T2/T3 | **Yes but burst-compromised** — Hα→4.9 Gyr killed (fading burst snapshot) | **DATA GAP** |
| UGCA 442 | 0.41 | 13.3→3.7 Gyr | T1/T2/T3 | **No** — no independent age chain | **DATA GAP** |

**Retired methods:** (1) t₉₀ budget analysis (treats Bayesian summaries as algebraically closed; dead-vs-active distinction is ~1% of mass, within SED fitting uncertainty). (2) Morphological prior for τ_main (morphology is transient in dwarfs; S0s prove disks/bars in ancient systems). Both independently confirmed retired by Pro and Deep Think.

### Screen 2 — Frosting Fallacy
**No frosting calculation used as backbone age evidence. All clear.**

Frosting-adjacent corrections caught by pipeline: IC 2574 (7.0→11.7), NGC 55 (7.2→10.8), NGC 300 (8.0→11.1), NGC 2366 (7.0→11.0), NGC 2976 (8.0→10.0), NGC 3109 (7.0→11.5), NGC 3741 (6.0→11.1), NGC 4068 (6.8→9.1), NGC 4214 (9.0→11.2), DDO 168 (3.3→7.7). NGC 2915 BCD decomposition (5.5→8.5). NGC 3198 T3→T2 (9.5→7.0). NGC 3893 T3→T2 (8.0→6.5, phantom data). NGC 3992 T3→T2 (8.5→7.2, Soft Valley + UV). NGC 4010 raw→corrected (7.0→6.0, edge-on dust). NGC 4088 T3→T2 (9.0→6.8, fourth Tully two-line trap). NGC 4389 T3→T2 (8.4→7.2, Infrared Supremacy). NGC 5055 T3→T2 (8.5→7.6, NUV-r Active boundary). NGC 5907 T3→T2 (10.0→7.0, Edge-On Dust Trap). Opposite direction (age underestimated by frosting-like mechanisms): NGC 7217 composite frosting trap (8.6→10.5, SF ring dragged global color blue). T1→T2 correction: UGC 00634 (DDO 7) T1/3.3→T2/6.0 ("Blue Dwarf Trap" — naive blue color mapping contradicted by source author's ≥10 Gyr constant SFH). AGN correction: NGC 2273 raw b = 0.54 → corrected b ≈ 0.25 (Warm Dust Diagnostic: S100/S60 = 1.49 = AGN torus heating inflated 60μm). IC 356 naive b ≈ 0.11 → true b < 0.1 (Cold Dust Diagnostic: S100/S60 = 4.18 = cirrus-dominated). UGC 04278 reclassified 5.0→6.0 Gyr (Edge-On Bifurcation / LSB Stability Block — Transparent/Primitive class). UGC 04305 (Holmberg II) starburst frosting on ancient cake: looks young (giant HII regions) but f₁₀ = 0.81 → 81% mass >10 Gyr. NGC 2552 (UGC 04325) "optical frosting" from post-burst B/A stars → (B-V)₀ = 0.25 (extremely blue) but backbone established. UGC 4483 Deep Think initial 10.0 Gyr T3 overturned → 5.0 Gyr T2 (Onset vs Dominance — RR Lyrae prove ancient onset, not ancient mass dominance). UGC 05005 T1/3.3 Gyr rejected → T2/6.0 Gyr (LSB Stability Block / Metallicity Trap at Blue Limit — B-V = 0.35 = bluest in sample but reflects Z ~ 0.1 Z☉, not youth).

### Screen 3 — Host Color Contradiction
**No unresolved contradictions.**

Resolved: F563-V1 (2014 value rejected, 3 sources agree); F583-1 (Bell rejected, de Blok adopted); F567-2 (column misread corrected); F579-V1 (subset average → per-object). NGC 2683 B-V uncited — bypassed via proxy. NGC 2955 SFRS column misread (ratios → fluxes). NGC 2841 UV-only upgraded to full proxy stack. NGC 3198 T3/9.5 rejected (νM contradicts). NGC 3769 disk/μ₀ colors → corrected to global (Tully+1996 top line). NGC 3877 citation corrected (Bianchi → Bouquin). NGC 3893 phantom data replaced. NGC 3953 disk/μ₀ → global (third Tully+1996 two-line trap). NGC 3992 T3/8.5 rejected — Soft Valley Protocol + UV confirms T2. NGC 4010 raw Table 4 → corrected Table 5 (edge-on dust, B-K' 4.14→3.37). NGC 4068 uncitable figure read replaced with McQuinn+2010b tabulated data. NGC 4088 disk/μ₀ → global (fourth Tully two-line trap, T3/9.0→T2/6.8). NGC 4138 α/Fe assignment reversed (counter-rotating disk is MORE α-enhanced). NGC 4214 upgraded 9.0→11.2 via Weisz+2011 Table 2 (Astra too conservative). NGC 4217 UV catalog conflict (Bouquin vs Bianchi) resolved via Backbone Tiebreaker. NGC 4389 T3/8.4 rejected — B-R only, overridden by B-K'/R-K/UV (Infrared Supremacy). NGC 4559 B-V corrected 0.49→0.45 (table column mixup with IRAS 12μm). NGC 5055 T3/8.5 rejected — NUV-r = 3.60 < 4.0 = Blue Cloud, not Green Valley; FUV-NUV = 0.54 confirms Active. NGC 5229 aggregator B-V = 0.66 rejected (aperture artifact); Makarova (1999) B-V = 0.48 preferred; Tier C superseded by Tier B (Smith+2022). NGC 5585 original source (VizieR J/AZh/88/342) rejected — region-by-region photometry, not integrated; replaced with Bouquin+2018 GALEX + RC3. NGC 5907 T3/10 rejected — uncited B-R = 1.41 is uncorrected edge-on dust artifact; GALEX FUV-NUV = 0.43 + RC3 (B-V)° = 0.52 confirm active. NGC 6503 distance inconsistency corrected (Strickland D = 5.2 vs SPARC Cepheid D = 6.26 Mpc; SFR rescaled, b 0.43→0.62, age 7.0→6.5). NGC 6674 aggregator B-V ≈ 0.84 rejected; McGaugh & Schombert (2014) B-V = 0.57 preferred (peer-reviewed). NGC 6789 Pro's peak-rate extrapolation (b = 3.8 → T2/7.7) corrected via Deep Think duty cycle challenge → McQuinn burst fraction 3.9% → T3/9.0. NGC 6946 UV discrepancy (Gil de Paz FUV-NUV = 0.31 vs SINGS 0.684) resolved as intrinsic vs observed; both Active; conservative floor adopted. NGC 7217 (UGC 11914) Astra's 8.6 Gyr (from global u-r = 2.32) corrected to 10.5 Gyr — Color Frosting Trap from young SF ring; component decomposition reveals >85% mass in 10–13 Gyr spheroid; "Tier 1" claim downgraded to Tier C. NGC 7331 red optical B-V = 0.87 would suggest T3, but GALEX FUV-NUV = 0.70 confirms Active (Blue Survival Principle); dust + massive bulge at i ≈ 70° explains red optical. UGC 00634 (DDO 7) T1/3.3 Gyr rejected — "Blue Dwarf Trap": (B-V)₀ = 0.39 is blue but source author (van Zee 2001 §3.3.1) explicitly states dIs are "not young systems" with ≥10 Gyr constant SF; b ≈ 0.51 confirms not newborn; reclassified T2/6.0. UGC 00731 (DDO 9) 1981 de Vaucouleurs B-V ≈ 0.76: rejected — physics conflict (red color contradicts active Hα in metal-poor dwarf); SPARC + Kaisin SFR used instead. UGC 01230 B-V = 0.42 → 4.9 Gyr reclassified to 6.0 Gyr — "Metallicity Trap" + "LSB Stability Block": naive solar-Z BC03 mapping underestimates age of metal-poor LSBs; created "Redder is Younger" inversion with UGC 00634 (B-V = 0.39 → 6.0). UGC 01281 B-V = 1.10 **REJECTED** (Source Veto) — author (Barazza & Binggeli 2002) explicitly states "clearly caused by internal dust absorption"; GALEX UV rescue (FUV-NUV = 0.34) confirms Active. UGC 02023 mass discrepancy (Gavilán log 8.62 vs SPARC log 8.82) resolved by adopting SPARC for consistency; b recalculated 0.24→0.15. IC 356 (UGC 02953) optical colors completely rejected — "IC 342 zone" heavy MW extinction; bypassed via IRAS far-IR + Cold Dust Diagnostic (S100/S60 = 4.18). NGC 2273 (UGC 03546) HI gas 5.3% would suggest T3, but SB(r)a bar funnels gas to center (HI→H₂ conversion invisible to 21cm); IRAS confirms active dust; raw b = 0.54 corrected to 0.25 via AGN Warm Dust Diagnostic. Contrast IC 356 (same HI%, no bar, T3). UGC 03580 Sa morphology would suggest quenched early-type, but log M★ = 9.82 too low for morphological quenching — Scale Paradox resolved: bulge visible but too shallow to suppress 66% gas disk. UGC 04278 HyperLEDA/aggregator B-V = 0.35 superseded by GALEX Atlas Table 4 B-V = 0.44; reclassified 5.0→6.0 Gyr via Edge-On Bifurcation + LSB Stability Block. NGC 2552 (UGC 04325) HyperLEDA B-V ≈ 0.40 superseded by Hunter+2009 (B-V)₀ = 0.25; post-burst diagnostic resolves extreme blueness as optical frosting. UGC 4483 Astra's 8.8 Gyr (T3) and Deep Think's initial 10.0 Gyr (T3) both overturned by Pro's challenge — Sacchi+2021 explicit regional t₅₀ ≈ 5 Gyr; "Onset vs Dominance" protocol: RR Lyrae prove t_start, not t₅₀. ANGST f₁₀ = 0.00 is depth artifact. UGC 04499 gas 142% > T1 threshold but t_dbl = 9.5 >> 4 Gyr; Local vs Global τ paradox (outer HI dynamically decoupled). UGC 05005 Astra's T1/3.3 Gyr rejected — B-V = 0.35 (bluest in sample) reflects metallicity trap, not youth; Pro's 4.6 Gyr revised to 6.0 Gyr (LSB Stability Block). UGC 05716 Astra's 5.0 Gyr (color only) corrected to 5.3 Gyr via metallicity correction — B-V = 0.42 at 25% solar Z corresponds to ~5.3, not 5.0. UGC 05721 (NGC 3274) aggregator optical colors (HyperLEDA B-V = 0.97 vs OpenNGC B-V = 0.37) **BOTH REJECTED** — "Aggregator Veto": >0.2 mag internal conflict disqualifies optical axis; SFRS GALEX UV (FUV-NUV = 0.208) adopted as sole classifier. NGC 3104 (UGC 05414) distance mismatch (van Zee 8.23 vs SPARC 9.40 Mpc) — minor, doesn't affect b calculation significantly. UGC 05750 Pro's 6.0 Gyr overridden to 5.7 Gyr by Twin Test with NGC 2552 (identical gas/b/UV). UGC 05829 GHASP Table B1 Hα flux was lower limit only (≥0.015 M☉/yr) — superseded by Gavilán+2013 true SFR (0.066 M☉/yr). "Broken Twin" with UGC 05716: identical mass/gas but opposite dynamical states. NGC 3432 (UGC 05986) i = 90° edge-on: absolute Hα flux severely attenuated (2–3×) → raw b is lower limit only; EW(Hα) = 64 Å robust (line+continuum both attenuated). Pro's 6.2 Gyr overridden → 5.8 Gyr via Edge-On Twin Test with NGC 3104.

### Screen 4 — Method Tier Violation
**No violations found.**

Notable non-violations: NGC 0891 (Pro's T3 overruled by activity — "Activity Trumps Age"); KK 251 ("Quality Floor" — Bronze data insufficient for T1). NGC 7217 (UGC 11914): Astra claimed "Tier 1" — downgraded to Tier C by Pro/Deep Think (component SSP spectroscopy ≠ tabulated cumulative SFH). Method rank corrected, not a violation of pipeline integrity. NGC 2273 (UGC 03546): raw IRAS b = 0.54 corrected to 0.25 via AGN Warm Dust Diagnostic — methodological correction (AGN contamination removed), not tier violation. UGC 4483: Deep Think initially ruled 10.0 Gyr (T3) — overturned by Pro's challenge citing Sacchi+2021 explicit regional t₅₀. Demonstrates pipeline self-correction: "Onset vs Dominance" distinction. Not a tier violation — a physics correction.

### Screen 5 — BCD / Starburst Contamination
**NGC 1705 is a BCD (Data Gap). NGC 2915 is a BCD but M★ from Werk et al. (not SPARC proxy) — NOT TRIGGERED. NGC 6789 is a BCD — M★ from McQuinn+2010 resolved CMD (not SPARC proxy); burst mass fraction 3.9% confirmed as frosting, not structural — NOT TRIGGERED. PGC 51017 (SBS 1415+437) is a BCD/WR — M★ from SPARC L[3.6] (Υ★ = 0.5 default), but age from Aloisi+2005 resolved CMD (Tier 1); gas fraction 259% independently confirms T1 — NOT TRIGGERED (age not derived from SPARC M★ alone). UGC 4483 is a starburst dwarf — age from Sacchi+2021 Tier 1 resolved stellar populations (HB + RR Lyrae + regional t₅₀), not SPARC M★ alone; "Late Bloomer" pathway established — NOT TRIGGERED. No verified galaxy uses BCD-contaminated SPARC M★ as sole age proxy. All clear.**

---

## Flagged Galaxies — Detailed Report

(None.)

---

## Updated Tally

| Tier | Count | Flagged/Demoted | Net |
|------|-------|-----------------|-----|
| T1 | 5 (D631-7, ESO 444-G084, NGC 3521, PGC 51017, UGC 00731) | 0 | 5 |
| T2 | 129 (incl. 2 PENDING: F568-V1, F571-V1) | 0 | 129 |
| T3 | 33 (IC 356, IC 2574, NGC 55, NGC 300, NGC 801, NGC 1167, NGC 2366, NGC 2403, NGC 2683, NGC 2841, NGC 2915, NGC 2976, NGC 3109, NGC 3741, NGC 3898, NGC 3900, NGC 3953, NGC 4013, NGC 4068, NGC 4157, NGC 4214, NGC 4217, NGC 5005, NGC 5289, NGC 5533, NGC 5985, NGC 6789, NGC 7217, NGC 7814, UGC 02885, UGC 04305, UGC 06614, UGC 11455) | 0 | 33 |
| Data Gap | 11 (DDO 170, ESO116-G012, ESO563-G021, F583-4, NGC 1705, UGC 02259, UGC 06818, UGC 07125, UGC 07866, UGCA 281, UGCA 442) | 0 | 11 |
| Poison Data | 1 (UGC 06930 — Aperture Fallacy, T1/4.0 Gyr purged) | N/A | 0 |

---

## Notes

- **Protocol 18 (τ_main Bracket Test) — FULLY APPLIED AND RESOLVED.** All three Smith et al. (0.10 ≤ f_burst < 0.50) galaxies bracket-tested. NGC 5229 saved by Tier C fallback (T2); UGCA 442 → Data Gap (no Tier C); UGC 07866 → Data Gap (Tier C burst-compromised per Age Checker Claude ruling). t₉₀ budget analysis and morphological prior for τ_main both RETIRED.
- Methods validated (≥2 uses): M2 (2×), M10 (4×), M12 (5×), M13 (2×), M17 (6×), M48 (2×), M51 (4×), Wyder+2009 LSB (2×: M102+M107), HaGS James+2004 (4×: M109+M111+M117+M119), Di Teodoro & Fraternali 2014 (5×), Wiegert radio (2×: M120+M121), Liu+2024 Table D2 (3×: M123+M127+…), LSB Stability Block (8+ applications), Noordermeer+2007 Table A3 (2×: M137+M140), Lelli+2013/2014 Table B.4 (4×: M144+M147+M148+…), Hallenbeck HIghMass (2×: M139+M146), Smith+2022 CIGALE sfh2exp (3× tested, 0 valid standalone — all require bracket test or Tier C fallback)
- Nine "Hidden T3" discoveries: IC 2574, NGC 55, NGC 2366, NGC 3109, NGC 3741 (ANGST Table 2). NGC 300, NGC 2976 (ANGST figure reads). NGC 4068 (McQuinn+2010b). NGC 4214 (Weisz+2011, "Starburst Illusion").
- Seven T3→T2 corrections: NGC 3198 (νM contradicts), NGC 3893 (phantom data), NGC 3992 (Soft Valley + UV), NGC 4088 (fourth Tully two-line trap), NGC 4389 (Infrared Supremacy), NGC 5055 (NUV-r Active boundary), NGC 5907 (Edge-On Dust Trap, 10→7.0)
- One T1→T2 correction: UGC 00634 (Blue Dwarf Trap, 3.3→6.0)
- NGC 3521: first massive spiral T1 — "Rejuvenation" via merger ~3 Gyr ago
- NGC 2915 BCD Bi-Modal Trap: burst/backbone decomposition (5.5→8.5 Gyr)
- NGC 4010/4013 Differential Correction pair validates Tully+1996 edge-on method
- NGC 4051 Seyfert 1: AGN contamination ruled out (cold dust diagnostic S_100/S_60 = 3.35)
- NGC 4138: counter-rotating disks — spectroscopic component integration. α/Fe rules out primordial burst
- NGC 4157: borderline T3 (B-K' = 3.81, +0.01 above threshold) confirmed by R-K = 2.56
- NGC 4183: Blue Anchor (B-K' = 2.43, R-K = 1.72 — bluest spiral in pipeline). Edge-on validation suite spans 2.43–3.99
- NGC 4217: Backbone Tiebreaker Protocol — UV catalog conflict resolved by R-K = 2.50 (Structure > UV for edge-on)
- NGC 4559: Face-On Blue Anchor (FUV-NUV = 0.20, bluest UV in pipeline). Gil de Paz+2007 source protocol approved
- NGC 5005: "Anemic Spiral" Anchor — spiral morphology persists past passive threshold (b ≈ 0.1, gas 1.3%). Completes T1→T2→T3 spiral evolutionary sequence
- NGC 5055 (M63): NUV-r < 4.0 = Active boundary formalized. Derivative (FUV-NUV) vs Integral (NUV-r) principle
- NGC 5229: **Smith+2022 bracket crosses T2/T3 — saved by independent Tier C (B-V = 0.48, Hα, gas 500%).** Protocol 18 bracket test: τ_main 0.2→10 gives t₅₀ 13.5→7.4 Gyr. Smith-derived 5.8 Gyr INVALID. Tier C fallback: Makarova B-V = 0.48 + Kaisin Hα + gas 500% all converge on T2/~5.0 Gyr. Confidence downgraded HIGH→MEDIUM
- NGC 5371: "Super-Spiral" — most massive star-forming spiral (log M★ = 11.6). Mass Penalty (Protocol 20) for Downsizing. Cirrus contamination noted for IR SFR
- NGC 5585: Chain-of-custody correction — region-by-region photometry (H II knots) rejected, replaced with integrated GALEX + RC3. "Aperture/Granularity Error." Cosmic Downsizing: bulgeless Sd younger than bulged Scd at same color
- NGC 5907: "Edge-On Dust Trap" formalized — for i > 80°, raw optical measures dust column, not age. UV binary check: FUV-NUV < 0.9 = active regardless of optical. Uncited B-R = 1.41 rejected
- NGC 5985: "Suppressed Super-Spiral" — second T3 pathway (Suppression vs Exhaustion). Gas 8% but engine broken. Seyfert AGN contamination of IRAS SFR. T3 Quenching Modes established
- NGC 6015: Conservative mass-weighted anchor. Observed B-V = 0.56 → 6.0 Gyr; corrected (B-V)_c = 0.41 is much bluer. "Skin vs Skeleton" — retain conservative estimate for Scd morphology
- NGC 6195: "Fading Super-Spiral" — b ≈ 0.25 despite 10.7% gas, morphological quenching by massive bulge. Cold dust diagnostic (S_100/S_60 = 4.8) distinguishes "Honest Red" (old stars) from "Dusty Disguise" (young stars + dust)
- NGC 6503: "The Lonely Galaxy" — Local Void isolation → "retarded evolution" with 27% gas retention. Distance consistency correction required (Strickland D = 5.2 vs SPARC Cepheid D = 6.26 Mpc). Environmental contrast: void (6.5 Gyr) vs cluster (7.0 Gyr)
- UGC 05999: "Fuel-Retention Anchor" — at constant b (~1.5), gas fraction governs age. 120% gas → 5.5 Gyr. Perfect interpolation between UGC 04499 (142%, 5.4) and UGC 05750 (66%, 5.7). ~0.1 Gyr per ~20% gas at constant activity. Wyder+2009 GALEX LSB Survey second application
- UGC 06399: **"Pro Override Anchor."** Pro classified T3/10.2 Gyr using exponential decay toy model — **REJECTED.** Toy model baked return fraction into b, used Salpeter without Kroupa correction. Standard b = 0.30, log sSFR = −10.66 — fails ALL T3 criteria (needs < −11.0, gas < 10%, b < 0.1). Fading sequence: DDO 87 (0.40, 7.3) → UGC 06399 (0.30, 7.6) → UGC 03205 (0.23, 7.9)
- UGC 06446: **"IMF Correction Anchor."** Pro flagged T1 (t_dbl = 3.5 Gyr) — **CORRECTED.** Two errors: mixed Salpeter SFR with Kroupa mass, mixed distance scales (though t_dbl is distance-invariant since M★∝D² and SFR∝D²). Corrected t_dbl = 5.5 Gyr → T2. b ≈ 2.5 with 279% gas. Second-youngest T2 at 4.6 Gyr
- UGC 06614: **"Giant LSB Bifurcation."** Distinct from dwarf LSBs: massive luminous bulge (M_B = −20.7) + diffuse disk. Dwarf LSBs (no bulge, 6.0 Gyr via Stability Block) vs Giant LSBs (old bulge dominates mass → T3). Inferred log sSFR ≈ −10.8 (Transition Zone). Astra's 8.5 adjusted to 8.0 Gyr (boundary). MEDIUM confidence — no direct M★, no UV
- UGC 06628: **"Iso-Age Diagonal Anchor."** Face-on (i=20°) provides pristine Hα — b is accurate, not lower limit. Opposing vectors cancel: higher b (1.3 vs NGC 3104's 1.21) pulls younger, lower gas (80% vs 102%) pulls older → same 5.8 Gyr. Pro's 5.6 Gyr (return fraction) overridden
- UGC 06667: **"Mass-to-Light Override."** No direct SFR found. M★/L_B = 0.5 proves vigorous state (young O/B stars dominate light). (B-V) = 0.65 is edge-on dust artifact (i=89°). Pro's 7.0 Gyr (LOW) overridden. Fuel Retention Gradient: UGC 05999 (120%, 5.5) → UGC 06667 (116%, 5.7) → NGC 3432 (113%, 5.8)
- NGC 3900 (UGC 06786): **"Morphological Quenching Anchor."** Fossil S0. Lb/Ld = 3.91 (bulge 4× disk). Central Hα in ABSORPTION (dead core >1 Gyr). Passes all T3: log sSFR = −11.2, b = 0.09, t_dbl = 153 Gyr. 14% gas present but dynamically stabilized (Toomre Q >> 1). Fifth Di Teodoro & Fraternali application
- NGC 3898 (UGC 06787): **"Relative Terminology Trap + Optical Lock."** Astra classified T2/7.0 Gyr based on "blue for Sa" — REJECTED. Noordermeer+2007 Platinum Source: B-R = 1.40 → absolutely RED on Hermes Grid (T3 > 1.3). "Blue for Sa" = still red in absolute terms. Optical Lock: massive early-type spirals with B-R > 1.3 → T3 null hypothesis. Matches NGC 2841. New protocols: Quantitative Supremacy, B-R Absolute Scale
- UGC 06818: **DATA GAP — "Ghost Data Protocol."** Excel entry shows sSFR = 5.3×10⁻¹¹ (log −10.28) → 7.0 Gyr T2, but M★ and SFR source values NOT FOUND. Likely HECATE/GSWLC-2 but catalog row unlocated. Deep Think: "Methodology ≠ Data." Pre-approved for 7.0 Gyr T2 MEDIUM once source retrieved and log sSFR confirmed ∈ [−10.0, −10.6]
- UGC 06917: "LSB Stability Block." B-V = 0.53 from dual sources (Moffat & Rahvar 2013 + Zonoozi & Haghi 2010, both UMa lineage). Astra's 6.2 Gyr adjusted to 6.0 per LSB Stability Block — color variations in LSBs are metallicity-driven, not age-driven. "Stereoscopic Verification" eliminates measurement error
- UGC 06923: **"Late-Stage Starburst Anchor."** High b (2.2) on depleted backbone (56% gas) = rejuvenation, not primordial burst. Pro flagged T1/3.7 Gyr (uncorrected Salpeter) → OVERRIDDEN. (B-V) = 0.42 confirms old+young mix. Iso-Age Diagonal: matches UGC 04499 (142%, b=1.46) at 5.4 Gyr — opposing vectors cancel
- UGC 06930: **POISON DATA — "Aperture Fallacy."** Excel T1/4.0 Gyr derived from Watson+2012 Σ_SFR (21″ aperture) scaled to full D₂₅ disk area. Central aperture density ≠ disk-average density. 2.19 dex swing proves value is noise. Deep Think: automatic REJECT. New protocol: aperture-scaled sSFR = Poison Data. Requires global photometry to re-admit
- UGC 06973 (IC 750): **"Quenching Hierarchy Anchor."** Activity > Gas Fraction. Gas at IC 356-level (6.5%) but b = 0.37 → fails BOTH T3 criteria. Pro proposed T3/9.4 Gyr → overridden. Terminal Circumnuclear Burn: last gas funneled to nucleus, outer disk dead, core still burning. Oldest active spiral at 8.2 Gyr. (B-V) = 0.64 is dust, not old stars
- UGC 06983: **"Radio IMF Scaling Anchor."** Radio SFR (M>5 M☉) requires ×3.44 net multiplier (×5.5 Salpeter total ÷1.6 Kroupa). Bar-driven nuclear burst in LSB disk — optical (B-R=1.05) sees quiet disk, radio sees buried engine. High-b Fuel Retention: at b≈2.5, gas fraction governs age → 5.1 Gyr
- UGC 07089: **"Bivariate Interpolation Anchor."** 2D (b, Gas%) calibration. Pro flagged T1/3.8 Gyr (uncorrected Salpeter ×5.5) → OVERRIDDEN by Radio IMF ×3.44. Bounded by NGC 2552 (67%, b=1.35, 5.7) and NGC 1003 (~50%, b=1.00, 6.0) → 5.9 Gyr. (B-R) = 0.70 remarkably blue for i=80°
- UGC 07125: **DATA GAP — "Twin Degeneracy."** 341% gas, i=90° (perfectly edge-on). T2 confirmed via Structural Veto (diffuse Sm, sub-critical). Age uncalibrated: could be 4.5 Gyr (bursting like UGC 05829) or 5.3 Gyr (simmering like UGC 05716). No SFR in any catalog. Edge-on extinction wall blocks all optical diagnostics. Needs radio or far-IR
- UGC 07151 (NGC 4144): **"Phoenix Protocol — Identity Contamination Recovery."** Original Excel used PGC 38218 = NGC 4140 (UGC 7063, Dec +01°) — WRONG OBJECT. Durbala+2020 footprint (Dec ≤ 36.4°) cannot contain UGC 7151 (Dec +46°). De Novo Provenance via Liu+2024 Table D2: log M★ = 9.33, log sSFR = −10.11. Excel 6.0 Gyr was "accidentally correct" (Stopped Clock). New protocol established: Identity Contamination → Phoenix recovery
- UGC 07232 (NGC 4190): **"b-Parameter Activity Discriminator."** Deep Think initially classified T3/9.0 Gyr ("8% burst → 92% ancient") — **Pro challenged and won.** b ≈ 3 proves substantial historical assembly (not fossil+burst where b would be 20–50). ~64% mass formed < 6 Gyr → T2. Aligned with study-peer NGC 5204. McQuinn+2015 Tier 1 CMD data (HST/ACS)
- UGC 07261 (NGC 4204): **"Multiwavelength Harmony Anchor."** Face-on (i=30°) eliminates dust — HaGS uncorrected Hα (0.103) matches SFRS UV+IR (0.107) perfectly. HaGS Overcorrection Trap: fixed A(Hα)=1.1 inflates SFR 2.75× for dust-poor systems. Bar-driven starburst in SBdm disk. b ≈ 2.8, 158% gas → 4.8 Gyr
- UGC 07323 (NGC 4242): **"Stochastic Flicker Anchor."** UV SFR (0.17) > Hα SFR (0.10) in patchy Sdm disk — galaxy in temporary 10 Myr Hα trough. UV (100 Myr) is correct tracer for t₅₀. Pro's Return-Fraction Trap: b = SFR×(1-R)×T/M★ gives artificially low b=0.80; correct b = sSFR×T = 1.14. Pro's 7.6 Gyr overridden → 6.4 Gyr. At b~1.0, 35% gas → older than 50% gas (NGC 1003, 6.0)
- UGC 07399 (NGC 4288): **"Clean Verification — Liu+2024 Platinum Source."** All protocols passed (Identity, Aperture, Ghost, Footprint). b ≈ 2.3 from log sSFR = −9.79. Youthful Pull below 6.0 equilibrium → 5.4 Gyr. Confirms Liu+2024 as verified stream for M106/UMa sector
- UGC 07524 (NGC 4395): **"Structural Veto + Adversarial Review Anchor."** Highest confirmed b (≈3.6) in T2. t_dbl = 3.8 Gyr technically passes T1 — **Structural Veto:** mature 10⁹ M☉ disk with SMBH and spiral arms cannot be "embryonic." Deep Think claimed AGN contamination → **Pro rebutted** (nuclear UV = 0.08% of galaxy FUV, O'Neill+2006). Deep Think conceded. Stochastic Flicker explains UV/Hα split
- UGC 07559 (DDO 126): **"Double Error Detection Anchor."** Error 1: Astra misread R_D = 0.87 kpc (disk scale length, Elmegreen & Hunter 2015) as (B-V) = 0.87 → false T3. Error 2: Pro used M_HI instead of M★ in sSFR denominator → b ≈ 0.9 (wrong). Corrected: true (B-V)₀ = 0.29 (blue!), b = 2.03. b/t_dbl interlock: b = 13.8/t_dbl must hold as consistency check
- UGC 07577 (DDO 125): **"Twin Anchor Standardization."** B-V = 0.59 from dual sources (FIGGS + Makarova 1999). Astra's 6.7 Gyr standardized to 6.0 Gyr — accepting 6.7 for B-V=0.59 while UGC 00128 (B-V=0.60) sits at 6.0 creates "Bluer is Older" artifact. CVn I cloud dIrr in Stability Block
- UGC 07603 (NGC 4455): **"Multi-Wavelength Extinction Trap Anchor."** At i=78°, Hα heavily attenuated but M★ (3.6μm) unaffected → false low b=0.73. EW(Hα) = 35 Å is extinction-proof (line/continuum cancel) → true b ≈ 1.0 (steady-state). 11HUGS Table 3. Fuel Retention: 137% gas at b~1.0 → 5.7 Gyr
- UGC 07608: **"Grid Cap / Ceiling Veto."** Grasha+2013 SED: constant-SF duration = 5 Gyr — but this is model grid MAXIMUM, paper explicitly warns "galaxies may have longer duration." Value is t ≥ 5 Gyr (censored datum), not point estimate. Astra's T1/2.5 Gyr rejected — if t ≥ 5 Gyr, T1 is statistically excluded. Standardized to LSB Stability Block (6.0 Gyr)
- UGC 07690: **"Interpolation Principle."** B-V = 0.38±0.05 from Dunn (2015) Table 6 UBVR photometry. Falls between UGC 05005 (B-V=0.35, 6.0) and UGC 00634 (B-V=0.39, 6.0) — both at 6.0 Gyr. Assigning 5.9 creates "Redder is Younger" inversion vs UGC 05005. Error bar spans entire Blue LSB Block → false precision. No Second Axis → Stability Block
- UGC 07866 (IC 3687 / DDO 141): **⚠️ DATA GAP — Protocol 18 + Tier C burst-compromised.** Smith T1/3.2 killed (bracket spans T1/T2/T3). Tier C Hα → 4.9 Gyr also killed: Age Checker Claude ruling — Smith burst parameters still trustworthy, Fading Burst Trap still holds, Hα measures declining burst snapshot not secular rate. Both age paths compromised. Data Gap until Tier A resolved CMD.
- UGC 08286 (NGC 5023): **"Scale-Height Extinction Differential Anchor."** At i≈90°, Hα (O-stars in thin midplane) catastrophically dust-attenuated; FUV (B-stars in thicker scale height) partially bypasses dust. FUV/Hα = 2.2× is diagnostic. Pro's HALOGAS-based 7.0 Gyr overridden — FUV Kroupa = 0.069 M☉/yr gives b = 1.51 (vigorous, not declining). Bivariate: same gas as NGC 3104 (102%) but higher b → 5.6 Gyr
- UGC 08490 (NGC 5204): **"Excel Error Correction + Mass-Sequence Alignment."** Excel B-V = 0.59 was contamination (likely from UGC 07577). Pro verified B-V = 0.41 (RC3 via Larsen & Richtler 2000). Deep Think rejected Stability Block (6.0) for this low-mass SA(s)m — bluer than NGC 5585 (0.46, 5.0), cannot be older. Anchored to late-type Sm sequence at 5.0 Gyr
- UGC 08699 (NGC 5289): **"Dust Trap — Morphology-Anchored T3 (NOT anchor)."** Noordermeer+2007 Table A3 B-R = 1.31±0.26. i = 77° creates extreme dust reddening — blue edge (1.05) permits T1 ages. No internal extinction correction. T3 maintained via morphology prior (S0/a) only. 8.8 Gyr is upper limit. MEDIUM-LOW confidence. Needs independent proxy for hardening
- UGC 08837 (Holmberg IV): **"False Precision Ruling."** de los Reyes & Kennicutt (2019) B-V = 0.43±0.04. Astra's 4.5 Gyr creates "Redder is Younger" inversion vs NGC 5204 (0.41, 5.0). Δ(B-V) = 0.02 < σ = 0.04 → statistically indistinguishable. Bounded by NGC 5204 (0.41) and NGC 5585 (0.46), both at 5.0 → Interpolation Principle → 5.0 Gyr. IBm dwarf
- UGC 09037: **"HIghMass — Data Correction + Birthrate Inversion Fix."** Hallenbeck+2016 Table 1: log sSFR = −9.53, b ≈ 4.1 (Vigorous). Astra's sSFR was 2× too low → wrong 6.2 Gyr. Pro's 5.5 Gyr overridden — Deep Think: b = 4.1 galaxy cannot be OLDER than b = 2.3 galaxy (UGC 07399, 5.4). b-Parameter Grid: b > 3 → Active T2 Floor (5.0 Gyr)
- UGC 09133 (NGC 5533): **"Forensic Distinction — Robust T3 Anchor."** Noordermeer+2007 Table A3 B-R = 1.48±0.13. KEY: i = 53° (moderate) vs UGC 08699 at 77° — seeing intrinsic stellar population, not dust "sunset." Blue edge (1.35) stays firmly T3 (~7 Gyr min). Independent confirmation: "hot disk" velocity dispersions (Noordermeer 2008) require Gyr-scale gravitational scattering. First confident T3 anchor from Noordermeer dataset. MEDIUM-HIGH confidence
- UGC 09992: **"Local Maximum Avoidance."** Gavilán+2013 Table 3: B-V = 0.37±0.04, sSFR = 8.6×10⁻¹¹ (log −10.07), b ≈ 1.19. Astra's 6.2 Gyr creates "spike" between UGC 05005 (0.35, 6.0) and UGC 07690 (0.38, 6.0) — physically invalid local maximum. Also b > 1 opposes older age. Standardized to 6.0 Gyr (Stability Block)
- UGC 10310 (DDO 204 / Arp 2): **"Second Axis Override."** van Zee (2001) Table 1: (B-V)₀ = 0.39±0.05, (U-B)₀ = −0.27±0.08. B-V matches UGC 00634 (6.0 Gyr) but U-B = −0.27 is extremely blue "Youth Flag" — proves young hot stars absent in 00634. Twin Anchor only applies when ALL indices match. SFR/⟨SFR⟩_past ≈ 1.9. Standardized to 5.0 Gyr (Active Cohort). Astra's 4.6 Gyr = False Precision vs NGC 5204
- UGC 11455: **"Cosmic Downsizing Anchor."** MANGROVE/GLADE + SPARC: log M★ = 11.29, gas ~7%, b = 0.23. First log M★ > 11 galaxy in recent queue. τ-model robust at high mass (no burst traps). Gas fraction 7% independently vetoes T2 — 93% baryons already stars. Even 4× dust correction (edge-on) stays T3. 10.3 Gyr. Catalog-level SFR, MEDIUM confidence
- UGC 11557: **"T2/T3 Boundary Anchor."** Lelli+2013 Table B.4: SFR = 0.486 (Salpeter) → 0.330 (Kroupa ×0.68). b = 0.75 (gently declining, not quenched), gas 43% (healthy), EW = 35 Å (active), τ_local = 4.0 Gyr (sustainable). All four physical state indicators favor T2 despite t₅₀ = 7.8 Gyr being near 8.0 boundary. Face-on (i ≈ 30°) = clean Hα. Boundary Decision Framework established
- UGC 11820: **"Twin Anchor + Blue Inversion."** van Zee & Haynes (2006) Table 6: (B-V)₀ = 0.37±0.02. Perfect Twin Match with UGC 09992 — identical B-V, both LSB, high gas, subsolar Z. Blue Inversion Physics: LSB dwarfs appear bluer than Active T2 (5.0 Gyr) despite being older (6.0 Gyr) because dust-poor/metal-poor → "naked" stellar population. 6.0 Gyr Stability Block
- UGC 12506: **"High Spin Delayed Assembly."** Hallenbeck+2014 HIghMass Table 1: log M★ = 10.46, log SFR = 0.40. **Column Mixup Corrected:** Astra read GF = 1.41 as SFR. Corrected sSFR = 8.7×10⁻¹¹ (log −10.06), b ≈ 1.20. False Precision Twin with UGC 09992 (Δ log sSFR = 0.01). Despite massive disk, high spin λ ≈ 0.15 delayed assembly → evolutionarily matches LSB dwarfs. Pro's 6.6 Gyr overridden (b > 1 cannot be older than equilibrium). 6.0 Gyr
- UGC 12632: **"Steady-State Archetype."** Lelli+2013 Table B.4 + SPARC (distance-harmonized). SFR = 0.046 (Kroupa, D²-corrected), b = 0.97 (essentially unity = constant SFH). For constant SFR: t₅₀ = half Hubble time ≈ 7.0 Gyr. Gas 268% validates flat SFH as sustainable (not transient). "Enormous fuel tank, narrow fuel line." Initial premature Data Gap corrected. HIGH confidence — analytically exact mapping
- UGC 12732: **"Metal-Poor EW Calibration Anchor — Galaxy #100."** Lelli+2014 Table B.4 + SPARC. b = 0.96, gas 439% (most extreme in sample). EW(Hα) = 88 Å is NOT a starburst — at low Z, old continuum is faint (metal-poor RGB) and massive stars produce more ionizing UV → Starburst99 predicts 80–100 Å for constant SFH. Same b as UGC 12632 but double EW = metallicity effect. Galaxy optically thin (i = 39°, low Z). 7.0 Gyr T2 HIGH
- UGCA 281 (I Zw 36 / Mrk 209 / Haro 29): **⚠️ DATA GAP — Prior T3/9.0 RESCINDED.** Reconciliation Exemplar. Feb 5 T3 built on (1) B-R ≈ 1.0 narrative citation — contradicted by Gil de Paz Table 3 (B-R ≈ 0.54); (2) Frosting Fallacy — f_burst(200 Myr) ≈ 8.5% constrains only present, not 0.2–4 Gyr assembly window. K-band mass (4×10⁷ M☉) resolves 3.6μm contamination but age remains unconstrained. Every constraint degenerate across T1/T2/T3 without Tier A/B data. BCD = pipeline's most dangerous galaxy type
- UGCA 442 (HIPASS J2343-31): **⚠️ DATA GAP — Both τ_main inference methods retired by Protocol 18.** Smith+2022: f_burst = 0.41. Bracket test: τ_main 0.2→10 gives 13.3→3.7 Gyr (spans T1/T2/T3). t₉₀ budget analysis RETIRED (treats Bayesian summaries as closed). Morphological prior RETIRED (morphology transient in dwarfs). No independent Tier C age chain. Feb 5 T1/2.9 and Feb 15 T3/9.0 both formally rescinded
- UGCA 444 (WLM / DDO 221): **"Tier A Gold Standard."** Albers+2019 Table 2: resolved CMD SFH reaching below oldest MSTO. Tabulated log(τ₅₀) = 9.71 → 5.1 Gyr. ACS field as global proxy (100:1 mass dominance over UVIS halo). t₉₀ = 1.05 Gyr rules out T1 (not young enough) and T3 (not passive). Isolated Local Group dIrr — no environmental effects. Inside-out growth. 5.1 Gyr T2 HIGH
- NGC 6674: McGaugh & Schombert (2014) elevated to Tier 1 for SPARC galaxies. Aggregator B-V ≈ 0.84 rejected (aperture bias or correction mismatch)
- NGC 6789: **BCD Bi-Modal Trap #2.** "Duty Cycle Error" — Pro extrapolated Lelli peak b = 3.8 over secular timescales; Deep Think challenged with McQuinn burst fraction = 3.9%. BCDs: do NOT extrapolate peak SFR. Second BCD case after NGC 2915
- NGC 6946: "Fireworks Galaxy" — UV discrepancy (0.31 vs 0.684) resolved as intrinsic vs observed. **Discrepant Sources Same Tier protocol:** if both values agree on tier, use observed/redder for conservative floor. 10+ supernovae confirm SF
- NGC 7217 (UGC 11914): "Bulge-Dominated Fossil Anchor." **Color Frosting Trap (Composites)** — global u-r = 2.32 dragged blue by SF ring (<5% mass). Component reconstruction: spheroid >85% mass at 10–13 Gyr. SSP ages ≤ mass-weighted ages for mixed populations. Astra's "Tier 1" downgraded to Tier C. b << 1 confirms profoundly quenched
- NGC 7331: **"Blue Survival Principle"** — if UV survives dust screen (FUV-NUV < 0.9), T2 confirmed. Dust always reddens; intrinsic UV bluer than observed 0.70. Red B-V = 0.87 from dust + massive Sb bulge at i ≈ 70°. "Milky Way's twin"
- NGC 7793: **Bluest UV in entire sample** (FUV-NUV = 0.17). "Hyper-Active Blue Limit." HST/GHOSTS resolved RGB confirms ancient skeleton under active skin. **Hyper-Active Triad:** NGC 4183 (IR Blue Limit), NGC 4559 (UV Control), NGC 7793 (Absolute UV Limit) = Zero Point of Active Spiral Sequence. PLATINUM confidence
- NGC 7814: **"Fossil Spiral"** — endpoint of disk galaxy evolution. b ≈ 0.03 (dead). Gas 1.6% (exhausted). **Quenching Triad completed:** Anemic (NGC 5005, dying) → Suppressed (NGC 5985, choked) → Fossil (NGC 7814, dead). WISE 22μm dust-robust for edge-on
- PGC 51017 (SBS 1415+437): **7th T1 galaxy.** BCD/WR "Late Bloomer." HST resolved CMD (Tier 1 data). RGB age 2.2 Gyr via Standard Model. Gas fraction 259% vetoes T3. **Methodological lock:** SPARC distance = Aloisi model distance → must accept associated age
- UGC 00128: **First LSB protocol.** "Transparency Principle" — LSB low dust makes optical colors reliable. **Luminosity-weighting priority** over area-weighting for mass tracking. Dual-source cross-check (de Blok+1995 = van der Hulst+1993). "Slow Burn" evolution
- UGC 00191: "Active Steady Builder." b = 1.66 (forming faster than lifetime average) but t_dbl = 8.3 Gyr → secular, not burst → T2. **T1/T2 discriminator: t_dbl < 3 Gyr or gas > 200% for T1.** Gas 129% insufficient alone
- UGC 00634 (DDO 7): **"Blue Dwarf Trap" — T1 REJECTED.** (B-V)₀ = 0.39 mapped to 3.3 Gyr by Astra. Source author (van Zee 2001) explicitly states dIs are "not young systems" with ≥10 Gyr constant SF. b ≈ 0.51 confirms not newborn. **Source Author Priority ("Context Rule"):** author SFH interpretation supersedes naive color→age mapping
- UGC 00731 (DDO 9): **T1 #5.** Gas fraction M_HI/M★ ≈ 11 (>90% baryons in gas). **Protocol 19: Gas Fraction Override** — M_HI/M★ > 5 strengthens T1 when t_dbl near boundary. b ≈ 3.76. 1981 optical color rejected (physics conflict). "Embryonic" assembly
- UGC 00891: **Key T1/T2 test case.** 229% gas but NOT T1 — t_dbl = 17.8 Gyr (slow assembly). "Gas-Rich Steady State." Same B-V = 0.46 as UGC 00191 but b = 0.78 vs 1.66 → 6.5 vs 5.0 Gyr. **Birthrate is finer clock when colors identical**
- UGC 01230: **"Metallicity Trap" + "LSB Stability Block."** B-V = 0.42 → 4.9 Gyr reclassified to 6.0 Gyr. Metal-poor LSBs appear bluer at same age. Standardized with UGC 00634 and UGC 00128 at 6.0 Gyr equilibrium
- UGC 01281: **Edge-On Dust Trap + UV Rescue.** B-V = 1.10 rejected (Source Veto — author explicitly flags dust). GALEX FUV-NUV = 0.34 (bluest edge-on in sample). Aligned with NGC 7793 (6.5 Gyr)
- UGC 02023 (DDO 25): **"Metallicity Trap" + Gas Veto.** B-V = 0.44 looks blue but 12+log(O/H) = 8.02 (20% solar) → metal-poor mirage. b = 0.15 (inefficient). Gas 73% vetoes T3. Aligned with F571-8 inefficient anchor
- UGC 02259: **Data Gap.** T2 confirmed via Grand Design morphology + Hα. No tabulated SFR or color. "Grand Design Dwarf" requires Gyr of stable rotation → rules out T1
- NGC 1156 (UGC 02455): **"Starburst Anchor."** b ≈ 2.8 (highest in session!) but NOT T1. **"Fuel Limit Discriminator":** M_gas/M★ = 0.44 → cannot double mass → evolved + rejuvenated, not primitive. B-V = 0.58 confirms old backbone. "Middle-aged galaxy having a mid-life crisis"
- NGC 1167 (UGC 02487): **"Ancient S0 Anchor" — OLDEST IN SAMPLE (10.8 Gyr).** First S0 lenticular. CALIFA IFU Tier 1 data: 83% mass > 6 Gyr. Downsizing at Cosmic Noon. **Morphological Quenching:** 7.3% gas exists but dynamically stable (bulge + shear → high Toomre Q). HIGHEST confidence
- UGC 02885 ("Rubin's Galaxy"): **sSFR Diagnostic formally established.** log sSFR = −11.47, τ ≈ 300 Gyr. **sSFR thresholds:** T1 > −9.5, T2 −9.5 to −11, T3 < −11. PLATINUM confidence. One of largest known spirals (~5× MW mass)
- UGC 02916: **"Twin Convergence"** with NGC 5033 — both 7.4 Gyr via opposite mechanisms (NGC 5033 starved, UGC 02916 stabilized). **Morphological Stabilization:** Sab bulge suppresses 37% gas → b = 0.55 (expected ~0.8). Pre-S0 phase
- IC 356 (UGC 02953, Arp 213): **"Dust-Bypassed Anchor."** In "IC 342 zone" — heavy MW extinction → optical colors unreliable. IRAS far-IR bypasses dust. **Cold Dust Diagnostic:** S100/S60 = 4.18 → cirrus-dominated (old-star heating, not SF). True b < 0.1. Gas 5.9%. T3/8.8 Gyr. Zone of Avoidance Protocol established
- UGC 03205: **"Fading Sab."** Completes **Sab Evolutionary Track:** Stabilized (UGC 02916, 37% gas, 7.4 Gyr) → Fading (UGC 03205, 17% gas, 7.9 Gyr) → Quenched (IC 356, 5.9% gas, 8.8 Gyr). **Mass-Morphology Cancellation** with NGC 6195 — identical b ≈ 0.23–0.25 → identical 7.9 Gyr despite different mass/type
- NGC 2273 (UGC 03546): **"Nuclear Ring Anchor."** SB(r)a Seyfert 2. **Warm Dust Diagnostic:** S100/S60 = 1.49 → AGN torus inflates IR → raw b = 0.54, corrected b ≈ 0.25. **Hidden H₂ via bar:** only 5.3% HI but bar funnels gas to center → HI→H₂ → nuclear ring SF. **Low HI ≠ Quenched when bars present.** Contrast IC 356 (same HI%, no bar, T3). **Warm/Cold Dust Diagnostic codified:** < 2.0 = AGN, 2.0–3.5 = SF, > 4.0 = cirrus
- UGC 03580: **"Blue Early-Type Anchor."** Sa morphology but 66% gas, b ≈ 0.88 — late-type metrics. **Scale Paradox:** at low mass (log 9.82), bulge is prominent for Sa classification but too shallow to suppress disk. **"Morphological Quenching Requires Mass"** — at log M★ < 10, even prominent bulges can't stabilize gas. EW(Hα) = 38.4 Å confirms genuine SF. Two-component system: older compact bulge + massive star-forming disk
- UGC 04278: **"Transparent Edge-On."** B-V = 0.44 despite edge-on = remarkably blue = low dust/metallicity. **Edge-On Bifurcation protocol:** Dusty/Evolved (B-V ≳ 0.8, 6.5–7.0 Gyr) vs Transparent/Primitive (B-V ≈ 0.4–0.5, 6.0 Gyr). Reclassified 5.0→6.0 Gyr (LSB Stability Block). Aggregator B-V = 0.35 rejected → GALEX Atlas B-V = 0.44 adopted
- UGC 04305 (Holmberg II, DDO 50): **ANGST Tier 1.** f₁₀ = 0.81 → 11.5 Gyr (T3). Famous "Swiss Cheese" HI + giant HII regions — looks young but 81% of mass >10 Gyr. "Starburst frosting on ancient cake." Bi-Modal Trap extends to dIrr (Im), not just BCDs
- NGC 2552 (UGC 04325): **"Post-Burst Anchor."** **Post-Burst Diagnostic established:** SFR(FUV) = 0.115 > SFR(Hα) = 0.079 → declining from burst peak. (B-V)₀ = 0.25 (extremely blue) is "optical frosting" from B/A stars, not extreme youth. **Codified:** FUV > Hα = declining, FUV < Hα = rising, FUV ≈ Hα = steady. HyperLEDA superseded by Hunter+2009
- UGC 4483: **"Late Bloomer" — NOT "Zombie."** Sacchi+2021 Tier 1: RR Lyrae prove >10 Gyr onset but regional t₅₀ ≈ 5 Gyr (outer). **Onset vs Dominance protocol:** t_start ≠ t₅₀. Deep Think initially ruled 10.0 Gyr T3 — **overturned by Pro's challenge.** ANGST f₁₀ = 0.00 is depth artifact. Two BCD pathways: Zombie (ancient mass dominates, T3) vs Late Bloomer (ancient onset + delayed assembly, T2)
- UGC 04499: **"Simmering Dwarf."** Gas 142% > T1 threshold BUT t_dbl = 9.5 >> 4 Gyr. **Local vs Global Paradox:** τ_global = 18.6 vs τ_local = 6.6 Gyr — outer HI dynamically stable, decoupled from SF. "Enormous fuel tank, narrow fuel line." Multi-parameter triangulation → 5.4 Gyr
- UGC 05005: **T1/3.3 REJECTED → LSB Stability Block 6.0 Gyr.** B-V = 0.35 = BLUEST in sample. **"Metallicity Trap at Blue Limit":** extreme blueness measures low opacity of metal-poor stars, not extreme youth. LSB Stability Block universality affirmed: B-V 0.35–0.60 all → 6.0 Gyr. To break block: need FUV-NUV < 0.25 or b >> 1
- NGC 2985 (UGC 05253): **"Terminal T2 Anchor" at 8.0 Gyr** — absolute edge of Active Sequence. log sSFR ≈ −10.93 (borderline!) but **Hα Tiebreaker:** EW = 14.4 Å + gas 19% → engine still sparking. Sab sequence complete: UGC 02916 (7.4) → UGC 03205 (7.9) → NGC 2985 (8.0) → IC 356 (8.8 T3). Three galaxies at 8.0 Gyr boundary via different mechanisms (NGC 1090, NGC 2985, F571-8)
- NGC 3104 (UGC 05414): **"Gold Standard Anchor."** ONLY galaxy with direct literature b measurement. SPARC b ≈ 1.13 vs van Zee internal b ≈ 1.21 — excellent convergence from independent methods. Validates entire b-parameter framework. Gas 102% but t_dbl = 12 Gyr → T2. **Color-Age anchor:** (B-V)₀ = 0.37 → 5.8 Gyr
- UGC 05716: **"Extreme Simmering Anchor."** 372% gas — highest gas fraction in T2! But t_dbl = 12.5 >> 4 Gyr. **Sub-Critical Envelope:** shallow potential + sub-critical density + low Z (25% solar, poor cooling) = gas trapped. More gas ≠ T1 when gas can't collapse. Metallicity correction: B-V = 0.42 at 25% Z → 5.3 Gyr (not 5.0)
- UGC 05721 (NGC 3274): **Hyper-Active Mass Gradient.** FUV-NUV = 0.208 (Hyper-Active). Conflicting aggregator optical (0.97 vs 0.37) **both REJECTED** ("Aggregator Veto" — >0.2 mag conflict → disqualified). UV-only classification. Positioned in mass gradient: massive Sd (6.5) → intermediate (6.0) → dwarf (5.0). "Bluer is Older" inversion = Cosmic Downsizing
- UGC 05750: **"LSB Twin Anchor."** **Twin Test with NGC 2552:** identical gas (66% vs 67%), b (1.4 vs 1.35), FUV-NUV (0.15 both) → must match at 5.7 Gyr. Pro's 6.0 overridden. LSB morphology = diffuse gas, sub-critical density, can't collapse rapidly. Clean UV tracer (LSB = dust-poor)
- UGC 05764 (DDO 83): **"Fading Dwarf" — Standard verification.** B-V = 0.53 + b ≈ 0.38 + log sSFR = −10.56. No traps. Redness is stellar (evolved population), not dust. Forming at ~40% of past average → past peak, winding down. 7.6 Gyr distinguishes from "Ancient/Fossil" T3 systems. No Astra baseline existed. MEDIUM confidence (Tier 2 data only)
- UGC 05829 (Spider Galaxy): **YOUNGEST T2 IN SAMPLE at 4.5 Gyr.** b = 2.33 (highest in T2!), (B-V) = 0.21, t_dbl = 4.3 (barely above T1 threshold). **"Broken Twin Paradox" with UGC 05716:** identical mass/gas (363% vs 372%) but opposite states — simmering vs BURSTING. "Same fuel tank, different ignition conditions." Gas fraction is NOT destiny
- DDO 87 (UGC 05918): **"Declining Dwarf Anchor."** 254% gas but b = 0.40, t_dbl = 34 Gyr. **Sub-Critical Density:** log Σ_SFR = −3.16 (far below Kennicutt-Schmidt threshold). "Starved by lack of density, not lack of gas." **Gas-Rich Dwarf Sequence completed:** Bursting(2.33/4.5) → Primitive(1.66/5.0) → Simmering(1.1/5.3) → Steady(0.78/6.5) → Declining(0.40/7.3)
- NGC 3432 (UGC 05986): **"Edge-On Twin Anchor."** i = 90°. EW(Hα) = 64 Å (robust — line+continuum both attenuated). Absolute flux severely attenuated → raw b is lower limit. **Physical Twin to NGC 3104:** same gas (113% vs 102%), same true b (~1.2). "Tilt NGC 3104 by 90° = NGC 3432." Pro's 6.2 Gyr overridden
- Anchors: NGC 0891 "Oldest Active" (7.8), NGC 1003 "Steady State" (6.0), NGC 1090 "Massive Inefficient" (8.0), **NGC 1167 "Ancient S0" (10.8)**, NGC 2841 "Red Limit/Fossil" (9.0), **NGC 2985 "Terminal T2" (8.0)**, NGC 2998 "Massive Steady State" (6.8), NGC 3726 "Blue Spiral" (5.5), NGC 4013 "Red End" (9.0), NGC 4183 "Blue Limit" (5.2), NGC 5005 "Anemic Spiral" (9.3), NGC 5033 "Fading" (7.4), NGC 5371 "Super-Spiral" (7.5), NGC 5985 "Suppressed" (10.0), NGC 6195 "Fading Super-Spiral" (7.9), NGC 7217 "Bulge-Dominated Fossil" (10.5), NGC 7793 "Hyper-Active Blue Limit" (6.5), NGC 7814 "Fossil Spiral" (10.5), **IC 356 "Dust-Bypassed" (8.8)**, **NGC 3104 "Gold Standard" (5.8)**, **UGC 02885 "Rubin's Galaxy" (9.0)**, **UGC 03580 "Blue Early-Type" (7.0)**, **UGC 05716 "Extreme Simmering" (5.3)**, **UGC 05829 "Youngest T2" (4.5)**, **DDO 87 "Declining Dwarf" (7.3)**
- Two PENDING Deep Think: F568-V1, F571-V1 (both T2)
- Error patterns: column misreads (NGC 0100, NGC 2955, NGC 4559), phantom data (NGC 3893), disk/global confusion (NGC 3769, NGC 3953, NGC 4088 — four instances), citation errors (NGC 3877), raw/corrected table confusion (NGC 4010), uncitable figure reads (NGC 4068), α/Fe reversal (NGC 4138), B-R-only classification (NGC 4389), NUV-r misinterpretation (NGC 5055), aggregator photometry artifacts (NGC 5229, NGC 6674), region-vs-integrated source error (NGC 5585), uncited uncorrected edge-on optical (NGC 5907), distance inconsistency (NGC 6503), BCD duty cycle extrapolation (NGC 6789), UV intrinsic/observed confusion (NGC 6946), composite color frosting trap + tier claim inflation (NGC 7217), Blue Dwarf Trap (UGC 00634), physics conflict with optical color (UGC 00731), **Metallicity Trap / LSB age inversion (UGC 01230)**, **Source Veto edge-on dust (UGC 01281)**, **mass source discrepancy (UGC 02023)**

---

*Audit in progress — awaiting remaining galaxy profiles.*
*Method Checker Claude, Team Hermes*

### Method entries (verbatim from the audit)

### Method 8: sSFR (Wang+2017 LVHIS) → constant SFH mass-return → t₅₀ [Direct metallicity + TRGB constraint]
- Data source(s): Wang et al. (2017) LVHIS Table 1: SFR = 0.0031 ± 0.0005 M☉/yr, log M★ = 6.9
- SPS model: BC03, Chabrier IMF, constant SFH
- Metallicity: Direct — Lee et al. (2003): 12+log(O/H) = 7.40; López-Sánchez et al. (2012): 7.37 ± 0.20 (~5% solar)
- TRGB detected (Karachentsev+2002): proves ≥1–2 Gyr old stars exist
- Mass-return fraction: f_return = 0.41
- Method validated: Variant of Method 6 (same math, different source/recycling)
- Known issues: SFR tracer variance ~4×; M★ systematic ~0.20 dex
- Strengths: All SFR tracers yield t₅₀ < 4 Gyr — T1 robust across entire scatter

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| ESO 444-G084 | 2.2 (+1.8/−1.0) Gyr | T1 (borderline) | MEDIUM | No | — |

---

### Method 17: ANGST resolved-star CMD → cumulative SFH fractions (Weisz+2011 Table 2) → t₅₀
- Data source(s): Weisz et al. (2011) Table 2 — per-galaxy cumulative mass fractions
- SPS model: ANGST CMD-fitting pipeline (resolved stars). **Tier 1 data.**
- Method validated: Yes (6 applications)
- Known issues: IC 2574 limited field coverage; NGC 55 limited A₂₅ coverage

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| IC 2574 | 11.7 Gyr (f₁₀ = 0.86) | T3 | HIGH | No | — |
| NGC 55 | 10.8 Gyr (f₁₀ = 0.63) | T3 | HIGH | No | — |
| NGC 2366 | 11.0 Gyr (f₁₀ = 0.67) | T3 | HIGH | No | — |
| NGC 3109 | 11.5 Gyr (f₁₀ = 0.79) | T3 | HIGH | No | — |
| NGC 3741 | 11.1 Gyr (f₁₀ = 0.68) | T3 | HIGH | No | — |
| NGC 4214 | 11.2 Gyr (f₁₀ = 0.71) | T3 | HIGH | No | — |

All: "Hidden T3" — Astra estimates superseded by Table 2. NGC 4214 = famous starburst dwarf, "Starburst Illusion" — current activity masks ancient backbone (f₆ = 0.95 → 95% of mass >6 Gyr). Williams+2011 cross-check (~75% >8 Gyr) consistent.

---

### Method 24: Resolved-star CMD → cumulative SFH figure read (Gogarten+2010 Figure 8) → t₅₀
- Data source(s): Gogarten et al. (2010) Figure 8: cumulative SFH to 5.4 kpc
- SPS model: CMD-fitting pipeline (resolved stars). **Tier 1 data.**
- Method validated: Single use (distinct from Method 17 — different paper, figure read vs table)
- Known issues: Figure-read (±0.8 Gyr); radial coverage to 5.4 kpc only (~14 kpc galaxy)
- Correction: Astra's 8.0 Gyr was a misread. Pro: ~0.35 at 12 Gyr, ~0.53 at 11 Gyr → t₅₀ ≈ 11.1. "Hidden T3."

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 300 | 11.1 Gyr (range 10.3–11.9) | T3 | HIGH | No | — |

---

### Method 29: Resolved-star CMD → spatially-resolved median age profile (Williams+2013) → t₅₀
- Data source(s): Williams et al. (2013), ApJ 765, 120 (ANGST XI). Spatially-resolved CMD SFHs to ~11 disk scale lengths. "Median age" defined as 50% cumulative mass = t₅₀.
- SPS model: CMD-fitting pipeline (resolved stars). **Tier 1 data.**
- Method validated: Single use (distinct from Methods 17/24 — different paper, direct median age readout vs cumulative fraction interpolation)
- Known issues: ±2.5 Gyr systematic uncertainty (affects value, not tier — even at 10.8 − 2.5 = 8.3, still T3)
- Strengths: First T3 where Astra's initial classification stood without revision. "Remarkably undisturbed disk" — minimal radial age gradient confirms uniformly old population.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2403 | 10.8 Gyr | T3 | HIGH (tier) / MEDIUM (value) | No | — |

---

### Method 35: Resolved-star CMD → cumulative SFH figure read (Williams+2010 ANGST IV Figure 10) → t₅₀ [Outside-In Quenching proxy]
- Data source(s): Williams et al. (2010), ApJ 709, 135, Figure 10: Outer Field cumulative SFH. 0.5 crossing at ~10 Gyr lookback.
- SPS model: CMD-fitting pipeline (resolved stars). **Tier 1 data** (for Outer Field).
- Method validated: Single use (distinct from Methods 17/24/29 — different paper, outer-field proxy with quenching physics)
- Known issues: Outer Field = only ~3% of stellar mass. Inner-1 (86% mass) depth-limited beyond 3 Gyr. Tier justified by authors' "similar ancient populations at all radii" + inside-out assembly physics. Astra misread figure as 8.0 Gyr; Pro corrected to ~10.0 Gyr.
- M81 Group member undergoing gas stripping → outside-in quenching exposes ancient backbone in outskirts.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2976 | 10.0 Gyr | T3 | MEDIUM | No | — |

---

### Method 50: McQuinn+2010b tabulated ⟨SFR⟩ + M★ → mass fraction → t₅₀ [Resolved CMD-derived]
- Data source(s): McQuinn et al. (2010b) ApJ 724, 49. Table 1: ⟨SFR⟩(0–6 Gyr) = 0.010 M☉/yr. Table 2: Total M★ = 3.2×10⁸ M☉. Derived: mass formed <6 Gyr = 6×10⁷ M☉ = ~19%. Mass formed >6 Gyr = ~81%.
- Method: 81% of mass older than 6 Gyr → t₅₀ must exceed 6 Gyr by definition. Uniform assumption in ancient bin → t₅₀ ≈ 9.1 Gyr (conservative; declining SFH would push older). Resolved CMD underlies tabulated values (Platinum quality).
- Method validated: Single use (distinct from Method 17 — different paper, SFR-averaged bins vs cumulative fractions)
- Known issues: Astra's figure read (6.8 Gyr) **REJECTED** — cumulative curve for NGC 4068 not shown in cited paper (Figure 14 is ESO 154-023 and NGC 784). McQuinn classifies NGC 4068 as "starburst dwarf" — refers to current rate, not mass-weighted age. "We classify the Cake, not the Frosting."

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4068 | 9.1 Gyr (range 8.5–10.5) | T3 (Hidden T3 #8) | MEDIUM | No | — |

---

### Method 75: Aloisi+2005 HST/ACS resolved CMD (RGB/AGB/TRGB) + SPARC distance consistency → t₅₀ [Late Bloomer BCD]
- Data source(s): Aloisi et al. (2005) ApJL 631, L45. HST/ACS resolved CMD: RGB + AGB + TRGB detected. I_TRGB = 26.67, (V-I)_TRGB = 1.28. Mass fraction ≥1.3 Gyr: ≥80%. Standard Model RGB age: 2.2 Gyr at D = 13.6 Mpc (Z = 0.001, A_V = 0). SPARC: L[3.6] = 0.155×10⁹ L☉, M_HI = 0.201×10⁹ M☉, D = 13.60 Mpc. M★ ≈ 7.75×10⁷ M☉. Gas fraction 259%.
- Method: **Tier 1 data (resolved CMD).** SPARC uses D = 13.6 Mpc = Aloisi Standard Model distance → must accept associated age (can't accept distance but reject age). ≥80% mass in RGB at ≈2.2 Gyr. Gas fraction 259% vetoes T3 — definitionally young ("Late Bloomer"). BCD/Wolf-Rayet galaxy with "young light over old mass" but even "old" mass is only ~2 Gyr. **Methodological lock:** distance and age are coupled outputs from same model.
- Method validated: Single use (seventh T1 galaxy)
- Known issues: Paper allows RGB ages 1.3 Gyr to Hubble time — resolution via distance/model consistency, not direct measurement. MEDIUM confidence due to model dependence.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| PGC 51017 (SBS 1415+437) | 2.2 Gyr | T1 (Late Bloomer BCD) | MEDIUM | No | — |

---

### Method 78: van Zee (2001) Table 1 (B-V)₀ + Table 2 birthrate + §3.3.1 author SFH interpretation → t₅₀ [Blue Dwarf Trap — T1 REJECTED]
- Data source(s): van Zee (2001) AJ 121, 2003. Table 1: (B-V)₀ = 0.39±0.03, (U-B)₀ = −0.30±0.06. Table 2: SFR = 0.073 M☉/yr, ⟨SFR⟩_past = 0.144 M☉/yr → b ≈ 0.51. §3.3.1: "global colors most consistent with nearly constant star formation for at least ~10 Gyr"; "dI galaxies are not young systems."
- Method: Astra mapped (B-V)₀ = 0.39 → BC03 → 3.3 Gyr (T1). **REJECTED.** b ≈ 0.51 (forming at half historical average — not a newborn galaxy signature). Source author explicitly states dI galaxies are not young, have ≥10 Gyr constant SF. **"Blue Dwarf Trap":** blue colors in metal-poor, gas-rich dwarfs mimic youth but reflect low metallicity and gas richness, not young mass-weighted age. Constant SFH over 12 Gyr → t₅₀ = 12/2 = 6.0 Gyr. **Source Author Priority ("Context Rule"):** author's astrophysical interpretation supersedes naive color→age lookup.
- Method validated: Single use
- Known issues: **T1→T2 RECLASSIFICATION.** Astra's 3.3 Gyr rejected. Blue colors ≠ young mass. Same "Skin vs Skeleton" / "Derivative vs Integral" principle. Eighth correction of this type.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 00634 (DDO 7) | 6.0 Gyr (range 5.0–7.0) | T2 (Blue Dwarf Trap resolved) | HIGH | No | — |

---

### Method 93: Weisz+2011 ANGST Table 2 cumulative mass fractions → interpolated t₅₀ [Starburst Frosting on Ancient Cake]
- Data source(s): Weisz et al. (2011) ApJ 739, 5 — ANGST Table 2 (Ho II row): f₁₀ = 0.81, f₆ = 0.81, f₃ = 0.82, f₂ = 0.85, f₁ = 0.89. SPARC: UGC 04305 = Holmberg II = DDO 50. D = 6.92 Mpc (from SPARC not stated in profile but known). Multi-field HST pointings merged after deduplication.
- Method: **Tier 1 data (resolved-star CMD, ANGST — oMSTO).** f₁₀ = 0.81 > 0.50 → t₅₀ > 10 Gyr → **T3 locked (HIGH confidence).** Interpolation: t₅₀ = 14 − (2/f₁₀) = 14 − 2.47 = **11.5 Gyr.** Grid consistent: NGC 2366 (f₁₀ = 0.67, 11.0 Gyr), UGC 04305 (f₁₀ = 0.81, 11.5 Gyr), IC 2574 (f₁₀ = 0.86, 11.7 Gyr). Famous for giant HII regions and "Swiss Cheese" HI morphology — looks young but 81% of mass >10 Gyr old. **"Starburst frosting on ancient cake"** — same pattern as NGC 2915. Bi-Modal Trap applies to dIrr (Im), not just BCDs. Deep Think validated interpolation formula as physically grounded (constant SFR assumption from Cosmic Dawn to 10 Gyr).
- Method validated: ANGST interpolation (fourth application after NGC 2366, NGC 55, IC 2574, NGC 300, etc.)
- Known issues: Numeric t₅₀ MEDIUM confidence (interpolation within broad bin); tier HIGH confidence (directly from f₁₀ > 0.5).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 04305 (Holmberg II) | 11.5 Gyr | T3 (Ancient dIrr) | HIGH (tier) / MEDIUM (numeric) | No | — |

---

### Method 95: Sacchi+2021 ApJ 911, 62 resolved stellar populations (HB + RR Lyrae + regional t₅₀) + ANGST Table 2 cross-check → Onset vs Dominance → t₅₀ [Late Bloomer — NOT Zombie]
- Data source(s): Sacchi et al. (2021) ApJ 911, 62: HB stars detected, 6 candidate RR Lyrae (>10 Gyr age floor), 87% mass >1 Gyr old. Regional t₅₀: Regions 3&4 (outer) ~5 Gyr, Regions 1&2 (inner) ~2–3 Gyr. Weisz+2011 ANGST Table 2: f₁₀ = 0.00 (depth artifact), f₆ = 0.91.
- Method: **Tier 1 data (resolved stellar populations).** **Critical distinction: Onset vs Dominance.** RR Lyrae prove ancient stars *exist* (t_start > 10 Gyr) but do NOT prove they *dominate the mass*. Sacchi+2021 explicitly states regional t₅₀ ≈ 5 Gyr (outer). **Deep Think initially ruled 10.0 Gyr (T3) — OVERTURNED by Pro's challenge:** "RR Lyrae prove clock started >10 Gyr ago, not that ancient epoch dominates mass." Deep Think conceded: "I committed a fundamental error by conflating Age of Onset with Median Mass Age." **"Late Bloomer" pathway established** — distinct from "Zombie" BCDs (NGC 2915): ancient onset + delayed accumulation → T2, not T3. Outer regions used for backbone characterization (avoid "Starburst Bias" from active cores). 87% mass >1 Gyr constrains against T1. ANGST f₁₀ = 0.00 is depth artifact (HST couldn't resolve oMSTO at this distance).
- Method validated: Single use
- Known issues: Astra's original 8.8 Gyr (T3) and Deep Think's initial 10.0 Gyr (T3) both overturned. Pro's challenge sustained.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 4483 | 5.0 Gyr | T2 (Late Bloomer) | HIGH | No | — |

---

### Method 124: McQuinn+2015 MNRAS 450, 3886, Table 2 resolved CMD burst fraction + b-parameter back-calculation → t₅₀ [Tier 1 data — HST/ACS CMD]
- Data source(s): McQuinn et al. (2015) Table 2: burst mass fraction = 8% (upper limit), burst duration ≥ 380 Myr. HST/ACS resolved CMD showing RGB + young stars. D = 2.83±0.17 Mpc (TRGB). Cross-ID: UGC 07232 = NGC 4190.
- Derived: b_recent ≈ 3.0, peak SFR ≈ 0.014 M☉/yr. Back-calculated: SFR_{t<6Gyr} ≈ 4.7×10⁻³ M☉/yr, mass in last 6 Gyr ≈ 2.8×10⁷ M☉ out of 4.4×10⁷ total → ~64% formed < 6 Gyr.
- Method: **Deep Think initially classified T3/9.0 Gyr** ("8% burst → 92% ancient") — **Pro challenged and won.** b ≈ 3 proves substantial historical mass assembly, not ancient backbone with frosting. If b were truly fossil+burst, b would be 20–50. b ≈ 3 means past average was ~30% of burst rate. >50% mass formed < 6 Gyr → t₅₀ must be < 6.0 Gyr → T2. Aligned with study-peer NGC 5204 (also 5.0 Gyr T2).
- Method validated: Single use (Tier 1 CMD data — highest data quality)
- Known issues: MEDIUM confidence. Deep Think's T3 overturned by Pro.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07232 (NGC 4190) | 5.0 Gyr | T2 (Active) | MEDIUM | No | — |

---

### Method 134: Smith+2022 MNRAS 515, 3270, Table 4 SED parameters (f_burst, τ_burst, t_burst) → **DATA GAP** — Protocol 18 + Tier C compromised by burst contamination
- Data source(s): Smith et al. (2022) Table 4: log M★ = 7.75, τ_burst = 399 Myr, t_burst = 929 Myr, **f_burst = 0.43**, t₉₀ = 439 Myr. CIGALE sfh2exp framework: t_main fixed at 10 Gyr, τ_main explored on grid (0.2, 0.5, 1, 5, 10 Gyr). Cross-ID: UGC 07866 = IC 3687 = DDO 141. D = 4.57 Mpc (TRGB).
- Derived: τ_main ≈ 4.4 Gyr (inferred from grid, NOT directly tabulated). t₅₀ ≈ 3.2 Gyr reconstructed from parameters — **INVALID per Protocol 18.**
- Method: **DATA GAP — Both age paths compromised.**
  - **Smith-derived T1/3.2 Gyr: KILLED by Protocol 18.** Bracket test: τ_main = 0.2 → t₅₀ ≈ 13.2 Gyr (T3); τ_main = 5 → 5.1 Gyr (T2); τ_main = 10 → 3.0 Gyr (T1). Spans all three tiers. f_burst = 0.43 < 0.50 → T1 exception does not apply.
  - **Tier C fallback T2/4.9 Gyr: COMPROMISED by burst contamination.** Ruling by Age Checker Claude: Smith burst parameters (f_burst, τ_burst, t_burst) remain trustworthy — Protocol 18 only kills our ability to read the old component. The Fading Burst Trap argument still holds: Hα is measuring a transient burst state that's already declining, not the galaxy's secular assembly rate. The Tier C T2/4.9 Gyr measures a snapshot, not a reliable t₅₀.
  - **Data Gap** until Tier A resolved CMD data becomes available.
- **Formal Method Hierarchy ratified:** A (tabulated cumulative SFH) > B (parametric SFH with mass fractions) > C (SFR proxies → t_dbl → b) > D (twin inference / structural)
- Method validated: N/A (both age paths compromised)
- Known issues: ⚠️ DATA GAP. Smith T1 killed (Protocol 18). Tier C T2 killed (burst contamination). Prior T2/4.9 Gyr deprecated. Prior T1/3.2 Gyr deprecated.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07866 (IC 3687 / DDO 141) | [NULL] | **DATA GAP** | NONE | N/A | Protocol 18 bracket spans T1/T2/T3; Tier C Hα measures fading burst snapshot, not secular rate. Both age paths compromised. |

---

### Method 149: UGCA 281 (I Zw 36 / Mrk 209 / Haro 29) — DATA GAP [Reconciliation Exemplar — Prior T3/9.0 RESCINDED]
- Data source(s): Karachentsev et al. (2021) Table 1: log LK = 7.82, M★(K-band) ≈ 4.0×10⁷ M☉, log SFR(Hα) = −1.29, log SFR(FUV) = −1.78, SFR(Hα)/SFR(FUV) ≈ 3.1 ("rising burst"). Schulte-Ladbeck et al. (2001): HST/NICMOS confirms RGB/AGB (old stars exist). SPARC: Quality 3 (Low), L[3.6] contaminated by nebular + hot dust. Gil de Paz & Madore (2004) Table 3: B-R ≈ 0.54–0.59±0.13–0.15 (tabulated CCD photometry). Noeske et al. (2004): B-R ≈ 1.0 (narrative citation of Papaderos 1996a). Lelli et al. (2014b) Table 1: 12+log(O/H) = 7.77 (Z ≈ 1/15 Z☉). Gas fraction (K-band corrected): ~155%. Blue Compact Dwarf / Wolf-Rayet galaxy.
- Method: **DATA GAP — Prior Lock Rescinded.**
  - **Feb 5 classification (T3/9.0 Gyr) RESCINDED:** Built on (1) B-R ≈ 1.0 from narrative citation, now contradicted by Gil de Paz Table 3 (B-R ≈ 0.54); (2) Frosting calc (f_burst ≈ 8.5% over 200 Myr) — but this is **Frosting Fallacy**: constrains only present burst, NOT the 0.2–4 Gyr assembly window that determines t₅₀.
  - **Feb 15 T2/~6 Gyr REJECTED:** SPARC M★ contaminated (3.6μm BCD trap), no concrete SFR chain.
  - **K-band mass resolves denominator** (M★ ≈ 4×10⁷, 2.4× lower than SPARC) but age of that mass remains unconstrained.
  - Every constraint is degenerate across T1/T2/T3. Without Tier A (resolved cumulative SFH) or Tier B (f_burst from SED), cannot classify.
- **Protocols established:**
  - **Frosting Fallacy:** f_burst over last N Myr constrains only present, not 0.2–4 Gyr assembly window. Cannot leap from "current burst is small" to "backbone is ancient."
  - **Tabulated > Narrative:** Gil de Paz Table 3 (with uncertainties) overrides Noeske→P96a narrative citation.
  - **Duplicate Detection:** Feb 15 re-investigation should have found Feb 5 lock before greenfield work.
- Method validated: Single use (BCD Data Gap)
- Known issues: ⚠️ RECLASSIFICATION from T3 to Data Gap. Pipeline self-correction.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGCA 281 (I Zw 36 / Mrk 209) | [NULL] | DATA GAP | NONE | N/A | BCD — all constraints degenerate; prior T3 rescinded |

---

### Method 150: Smith+2022 Table 4 f_burst=0.41, τ_burst=100, t_burst=504, t₉₀=369 → **DATA GAP** [Both inference methods retired — Protocol 18]
- Data source(s): Smith et al. (2022) MNRAS 515, 3270, Table 4: log M★ = 7.26, f_burst = 0.41, t_burst = 504 Myr, τ_burst = 100 Myr, t₉₀ = 369 Myr. SPARC: L[3.6] = 0.140×10⁹ L☉, M_HI = 0.263×10⁹ M☉, D = 4.35 Mpc (TRGB-class), i = 64°, Vflat = 56.4 km/s, Q = 1. Gas fraction (SPARC) ~376%. Cross-ID: UGCA 442 = HIPASS J2343-31. Sm dwarf.
- Method: **DATA GAP — Both τ_main inference methods retired by Protocol 18.**
  - **Feb 5 T1/2.9 Gyr (morphological prior): RETIRED.** Morphology is transient in dwarfs and doesn't uniquely constrain mass assembly history.
  - **Feb 15 T3/9.0 Gyr (t₉₀ budget analysis): RETIRED.** Treats Bayesian summary statistics as algebraically closed. The dead-vs-active main component distinction is ~1% of total mass — well within SED fitting uncertainty.
  - **Bracket test (Protocol 18):** τ_main = 0.2 → t₅₀ ≈ 13.3 Gyr (T3); τ_main = 5 → 5.8 Gyr (T2); τ_main = 10 → 3.7 Gyr (T1). **Bracket spans T1 through T3 → DATA GAP.**
  - f_burst = 0.41 < 0.50 → burst does not hold mass majority → T1 exception does not apply.
  - No independent Tier C age chain available.
- **Protocol established: Protocol 18 ratified.** For Smith et al. (0.10 ≤ f_burst < 0.50), you CANNOT reverse-engineer τ_main — not via t₉₀ math, not via morphological priors. Must bracket and test tier stability.
- Method validated: N/A (methods retired)
- Known issues: ⚠️ RECLASSIFIED T3→DATA GAP. Prior T1 (Feb 5) and T3 (Feb 15) both formally rescinded.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGCA 442 | [NULL] | **DATA GAP** | NONE | N/A | Protocol 18 — bracket spans T1/T2/T3; both τ_main inference methods retired |

---



## IFU / Spectroscopic

**Conversion:** Spectroscopy typically provides mass fractions in coarse age bins (e.g., “young”, “intermediate”, “old”) or a cumulative mass fraction by age.  
We compute cumulative mass fraction from the youngest bin upward and identify the bin where the cumulative crosses 0.5. If age bins are wide, we treat the bin as uniform and interpolate within that bin to obtain a scalar \(t_{50}\).

### Worked example

**Worked example (Method 38 — IFU mass fractions → median age):**  
        The paper reports total stellar mass fractions by age bin (old/intermediate/young).  
        Since the “old” bin is below 50% while “old + intermediate” exceeds 50%, \(t_50\) lies in the intermediate bin and is taken as the median age within that bin (as recorded in the audit).

        ### Method 38: IFU spectroscopic decomposition (Coccato+2018) → mass fractions → median age → t₅₀
- Data source(s): Coccato et al. (2018) MNRAS 477, 1958, Figure 7: Total mass fractions — Old (≥7 Gyr) = 39.3 ± 5.3%, Intermediate (~3 Gyr) = 58.8 ± 5.2%, Young (≤1 Gyr) = 1.9 ± 0.1%. VIRUS-W IFU (105×55″).
- Method: Old < 50% → median (t₅₀) must cross into Intermediate bin → t₅₀ ≈ 3.2 Gyr. Mean age (5.7 Gyr) differs from median — pipeline uses median.
- Method validated: Single use
- Known issues: FOV limited to inner galaxy (bulge + inner disc). Deep Think ruled this establishes **upper limit** on t₅₀ via inside-out formation (outer disc younger). Bulge is 85% intermediate — pseudo-bulge/merger remnant, consistent with "bubble halo" tidal features.
- **First massive spiral in T1.** "Rejuvenation" pathway: merger ~3 Gyr ago rebuilt bulge and triggered widespread SF.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3521 | 3.2 Gyr | T1 (Rejuvenated Spiral) | MEDIUM | No | — |

---

### Method entries (verbatim from the audit)

### Method 30: GC spectroscopy (Proctor+2008) + boxy bulge morphology (Kuzio de Naray+2009) → proxy age → t₅₀
- Data source(s): Proctor et al. (2008) MNRAS 385, 1709 Table 4: 14/15 GCs aged 10–12 Gyr, 1 GC at 3.3 Gyr (central minor starburst). Kuzio de Naray et al. (2009) AJ 138, 1082: boxy bulge confirmed at all wavelengths.
- SPS model: Not BC03 — direct spectroscopic ages of individual GCs (Platinum quality)
- Method validated: Single use
- Known issues: Proxy-based (GC ages ≠ field star t₅₀ directly); original B-V = 0.85 was uncited and unreliable (edge-on dust); age assigned from GC system + structural proxy, not integrated light
- Key physics: Edge-on (i ≈ 84°) makes optical colors fail. GC spectroscopy and bulge morphology bypass dust. Boxy/peanut bulge requires bar + secular buckling over multiple Gyr.
- Deep Think "UFO Rescue" — recovered from data gap via proxy evidence.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2683 | ~10.0 Gyr | T3 | HIGH | No | — |

---

### Method 53: Spectroscopic component decomposition (Pizzella+2014) → mass-weighted t₅₀ [Counter-rotating disks]
- Data source(s): Pizzella et al. (2014) A&A 570, A79. VIRUS-W IFU spectroscopy. Lick indices + α/Fe-sensitive SSP grids. Main (co-rotating) disk: ~6.6 Gyr, [α/Fe] = +0.08. Counter-rotating disk: ~1.1 Gyr, [α/Fe] = +0.24, concentrated in star-forming ring.
- Method: Tier-1 spectroscopic input → Tier-2 conversion. Main disk dominates most radii. Global t₅₀ ≈ f_main × t_main + f_counter × t_counter ≈ 6.2 Gyr. Low α/Fe of main disk rules out primordial burst (T3) scenario.
- Method validated: Single use
- Known issues: Astra α/Fe assignment **REVERSED** (corrected: counter-rotating disk is MORE α-enhanced, not less). Exact mass partition not tabulated — light dominance used as proxy. Large uncertainty on main disk (±3.6 Gyr) but chemical abundances support central estimate.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4138 | 6.2 Gyr (range 5.0–7.5) | T2 | MEDIUM | No | — |

---

### Method 71: Fabricius+2014 structural decomposition + Sil'chenko+2011 component SSP ages + De Looze+2011 SFR → mass-weighted component reconstruction → t₅₀ [Bulge-Dominated Fossil]
- Data source(s): Fabricius et al. (2014) MNRAS 441, 2212: u-r = 2.32 (global, 9′ aperture, extinction-corrected), disc light fraction ~20%, spheroid ~80%, disc M/L < spheroid M/L. Sil'chenko et al. (2011) MNRAS 414, 3645: bulge SSP 10–13 Gyr, inner disc SSP ~5 Gyr, outer disc ~2–3 Gyr, outer SF ring ~1 Gyr. De Looze et al. (2011) MNRAS 416, 2712 Table 3: SFR(24μm) = 0.30±0.11, SFR(FUV+24μm) = 0.49±0.10 M☉/yr.
- Method: **Color Frosting Trap (Composites)** — Astra mapped global u-r = 2.32 → BC03 → 8.6 Gyr. Too low: young outer SF ring contributes disproportionate blue light despite <5% of mass. Spheroid contains >85% of mass at 10–13 Gyr SSP. Mass-weighted reconstruction: if >85% of mass is 10–13 Gyr, median mass-assembly time is mathematically bound to ancient spheroid. SSP ages underestimate mass-weighted ages for mixed populations (SSP ≤ mass-weighted, always). SFR cross-check: b << 1 (0.07–0.23) confirms profoundly quenched. **Astra's "Tier 1" claim downgraded to Tier C** — component spectroscopy ≠ tabulated cumulative SFH.
- Method validated: Single use
- Known issues: Tier C method rank. No tabulated cumulative SFH. M★ not directly tabulated (estimated from morphological class). SSP-to-mass-weighted conversion involves assumptions. Ring SF is real but localized (<5% mass).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 11914 (NGC 7217) | 10.5 Gyr | T3 (Bulge-Dominated Fossil) | HIGH | No | — |

---

### Method 85: CALIFA IFU mass fractions by age bin (MNRAS 484, 4298 Table 3) + B/T decomposition + Downsizing → t₅₀ [Ancient S0 — OLDEST IN SAMPLE]
- Data source(s): CALIFA IFU analysis (MNRAS 484, 4298) Table 3: Bulge mass fractions: Young (<1 Gyr) = 0.00±0.00, Intermediate (1–6 Gyr) = 0.07±0.03, Old (>6 Gyr) = 0.93±0.03. Disc mass fractions: Young = 0.01±0.01, Intermediate = 0.21±0.04, Old = 0.78±0.04. Table 2: B/T(r-band) = 0.34. Emission-line: low EW(Hα), BPT Seyfert → AGN, not SF. SPARC: L[3.6] = 489.955×10⁹ L☉, M_HI = 17.963×10⁹ M☉, D = 69.10 Mpc, Vflat = 332 km/s. M★ ≈ 2.45×10¹¹ M☉ (log 11.39). Gas fraction ~7.3%.
- Method: **Tier 1 data (IFU stellar populations).** Combined f_old = 0.34×0.93 + 0.66×0.78 = **0.83** (83% of mass > 6 Gyr). t₅₀ falls 40% into "old" bin; linear estimate = 9.2 Gyr, but **Downsizing** for super-massive galaxies (log M★ ≈ 11.4, formation at z ~ 2–3 "Cosmic Noon") skews distribution toward old edge → t₅₀ = 10.8 Gyr. **Morphological Quenching** explains gas paradox: 7.3% gas exists as extended HI ring but massive bulge + rotational shear → high Toomre Q → gas is dynamically stable, cannot collapse. "Engine broken." First S0 lenticular in sample. **HIGHEST confidence.**
- Method validated: Single use
- Known issues: Downsizing correction is physics-based inference within old bin. Gas presence doesn't affect T3 when morphologically quenched.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 1167 (UGC 02487) | 10.8 Gyr | T3 (Ancient S0 — OLDEST) | HIGHEST | No | — |

---



## Parametric SED fitting

**Conversion:** Parametric SED fitting returns an assumed SFH form \( \mathrm{SFR}(t;\theta)\) (e.g., constant, exponential \(\tau\), or two-component).  
Define cumulative formed mass:
\[
M(t)=\int_0^t \mathrm{SFR}(t')\,dt'.
\]
Solve for \(t_\*\) such that \(M(t_\*)=0.5\,M(T)\), where \(T\) is the model galaxy age; then
\[
t_{50,\mathrm{lookback}}=T-t_\*.
\]
When SED parameters are degenerate (outshining / “Protocol 18”), we do not invert \(\tau_{\rm main}\) from incomplete summaries; we apply the documented bracket test and/or Tier‑C fallback exactly as recorded in the audit.

### Worked example

**Worked example (Method 60 — outshining / bracket test):**  
        This entry shows why incomplete SED summaries cannot be inverted to a unique \(t_50\) when the main component is degenerate.  
        The audit documents the bracket outcomes and the Tier‑C fallback used to lock the final \(t_50\) (with the raw SED parameters shown in the entry).

        ### Method 60: Smith+2022 Table 4 CIGALE sfh2exp SED parameters → f_burst = 0.20 → **Smith bracket crosses T2/T3 — TIER C FALLBACK** [First Smith et al. galaxy]
- Data source(s): Smith et al. (2022) MNRAS 515, 3270. Table 4: log M★ = 7.48, f_burst = 0.20, t_burst = 1627 Myr, τ_burst = 464 Myr, t₉₀ = 1122 Myr. Supporting: Makarova (1999) B-V = 0.48, Kaisin & Karachentsev (2007) Hα SFR = 0.014 M☉/yr. Gas fraction ~500%.
- Method: **Protocol 18 bracket test:** τ_main = 0.2 → t₅₀ ≈ 13.5 Gyr (T3); τ_main = 5 → 9.3 Gyr (T3); τ_main = 10 → 7.4 Gyr (T2). **Bracket crosses T2/T3 boundary.** Smith-derived t₅₀ = 5.8 Gyr is INVALID (depended on inferred τ_main).
- **Independent Tier C saves classification:** Makarova (1999) B-V = 0.48 slots into Active T2 sequence (cf. NGC 5585 B-V = 0.46, 5.0 Gyr). Kaisin Hα confirms active star formation. Gas fraction ~500% with active Hα = not quenched. SBd morphology with extended disk. All Tier C indicators converge on T2 independently of Smith. **Tier C age: ~5.0 Gyr T2.**
- **Note:** Original t₅₀ = 5.8 Gyr (Smith) replaced by Tier C estimate ~5.0 Gyr. Tier unchanged (T2). Confidence downgraded from HIGH to MEDIUM (Tier C, not Tier B).
- Method validated: Single use (Protocol 18 first bracket test)
- Known issues: Aggregator B-V = 0.66 **REJECTED** (aperture artifact for edge-on); peer-reviewed Makarova B-V = 0.48 preferred. Edge-on SBd. Smith bracket crosses T2/T3 but Tier C independently confirms T2.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5229 (UGC 8550) | ~5.0 Gyr | T2 (Steady Edge-On Builder) | MEDIUM | No | Smith bracket crosses T2/T3; saved by independent Tier C (B-V = 0.48, Hα, gas 500%) |

---

### Method entries (verbatim from the audit)

### Method 58: Noll+2009 CIGALE SED Table 2 (log M★, log SFR) → sSFR → b → anchor interpolation → t₅₀
- Data source(s): Noll et al. (2009), CIGALE SED fits for SINGS sample. Table 2: log M★ = 10.77 ± 0.11, log SFR = 0.31 ± 0.20, A_FUV = 1.88 mag.
- Method: log sSFR = −10.46 ± 0.23 → b ≈ 0.48. Even worst-case (−0.23 dex) gives b ≈ 0.28, still 3× above T3 threshold. Interpolation between ESO 079-G014 (b ≈ 0.6, 7.2 Gyr) and fading line (b = 0.4, 7.6 Gyr). Deep Think adjusted 7.3→7.4 Gyr.
- Method validated: Single use
- Known issues: ±0.23 dex uncertainty does not threaten tier. Significant dust (A_FUV = 1.88) confirms gas-rich ISM still present (contrast NGC 5005's 1.3%).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5033 | 7.4 Gyr | T2 (Winding Down) | MEDIUM | No | — |

---

### Method 132: Grasha+2013 Table 2 SED fit (constant-SF duration) → Grid Cap / Ceiling Veto → LSB Stability Block → t₅₀ [Censored Datum]
- Data source(s): Grasha et al. (2013) Table 2: UGC 7608 = PGC 41063. Best-fit SFH = C (constant), Age = 5 Gyr (duration), EW(Hα) = 49.5 Å, IRX = −0.501. Alternate: SFH = B (burst), Age = 0.1 Gyr. D = 7.76 Mpc.
- Method: **Grid Cap / Ceiling Veto.** Grasha's model grid caps constant-SF duration at 5 Gyr — paper explicitly warns "galaxies may have longer duration" and pile-up at 5 Gyr is "artificial." Value is t ≥ 5 Gyr (lower limit), NOT t = 5 Gyr (point estimate). Astra proposed T1/2.5 Gyr via t₅₀ = T/2 — **REJECTED.** If t ≥ 5 Gyr, then t₅₀ ≥ 2.5 Gyr, which statistically excludes T1 (<4 Gyr) since true duration likely >> 5 Gyr. LSB dwarf → Stability Block (6.0 Gyr).
- **NEW PROTOCOL: Grid Cap / Ceiling Veto** — When derived age parameter pins to boundary of model's parameter space, treat as censored datum (lower limit for age). Reject any tier claim conflicting with limit direction. For dwarfs/LSBs, revert to Class Prior (6.0 Gyr).
- Method validated: Single use (first Grasha application)
- Known issues: Astra's T1/2.5 Gyr rejected. Dual solution ambiguity (constant vs burst) unresolved but burst alternative (0.1 Gyr) also problematic.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07608 | 6.0 Gyr | T2 (LSB Stability Block / Grid Cap) | HIGH | No | — |

---



## Broadband color → SPS

**Conversion:** For a set of observed integrated colors \(\mathbf{C}_\mathrm{obs}\) (e.g., \(B\!-\!V\), \(B\!-\!R\), \(V\!-\!I\), \(B\!-\!K\)), we compare against a precomputed SPS grid (BC03) spanning:
- metallicity \(Z\),
- SFH family (SSP, constant SF, exponential \(\tau\)),
- age steps.

We choose the best-matching model by nearest-neighbor or \(\chi^2\) in color-space:
\[
\chi^2=\sum_k \frac{(C_k-C_{k,\mathrm{obs}})^2}{\sigma_k^2},
\]
using whatever subset of colors exists for that galaxy. The inferred model’s \(t_{50}\) is taken from the same SFH by integrating the model SFH (as above), or from the grid’s tabulated half-mass time when available.  
Where the audit specifies protocol overrides (e.g., LSB metallicity corrections, Valley/Soft-Valley tie-breakers), those are part of the documented conversion.

### Worked example

**Worked example (Method 23 — BVR colors → BC03):**  
        From the published magnitudes:
        \[
        (B-V) = B - V,\quad (B-R)=B - R.
        \]
        These colors are then matched to the BC03 grid under the method’s stated metallicity and SFH priors, yielding the recorded \(t_50\).

        ### Method 23: BVR integrated color (Kassin+2006 Table 3) → BC03 τ-model (near-solar Z) → t₅₀
- Data source(s): Kassin et al. (2006) Table 3: B = 11.05 ± 0.03, V = 10.44 ± 0.03, R = 10.01 ± 0.03
- SPS model: BC03, Chabrier IMF, declining τ-model
- Metallicity: Assumed near-solar (bright SBbc spiral prior)
- Method validated: Single use
- Known issues: Metallicity assumed; age-Z-dust degeneracy
- Two-color consistency: B-V = 0.61, B-R = 1.04 agree

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 289 | 6.7 Gyr (range 5.5–7.5) | T2 | MEDIUM | No | — |

---

### Method entries (verbatim from the audit)

### Method 1: B-V integrated color (Begum+2003) → BC03 τ-model → t₅₀
- Data source(s): Begum, Chengalur & Hopp (2003) — B, V magnitudes within Holmberg isophote
- SPS model: BC03, Chabrier IMF, τ = 1 Gyr baseline
- Metallicity: Berg et al. (2012) L-Z relation (estimated, not measured)
- Dust: Foreground (MW) only
- Method validated: Untested (single use)
- Known issues: No photometric uncertainties; metallicity estimated; single-color degeneracy

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| CamB | 4.1 (+2.2/−1.7) Gyr | T2 (borderline T1/T2) | MEDIUM-LOW | No | — |

---

### Method 3: B-V integrated color (Begum+2008) → BC03 τ-model → t₅₀ [Direct spectroscopic metallicity]
- Data source(s): Begum, Chengalur & Karachentsev (2008) Table 1 — B-V color
- SPS model: BC03, Chabrier IMF, τ = 2.0 Gyr fiducial (bracket 1.5–3.0)
- Metallicity: Direct spectroscopic — Vaduvescu et al. (2007) Table 4: 12+log(O/H) = 7.81 ± 0.20
- Dust: Foreground; high galactic latitude makes it negligible
- Method validated: Single use
- Known issues: No photometric uncertainty published; extinction correction status unverified

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| D631-7 | 3.6 (+1.6/−0.9) Gyr | T1 | MEDIUM | No | — |

---

### Method 4: B-R integrated color (Schulte-Ladbeck+1998) → BC03 extended SFH → t₅₀ [Direct T_e metallicity]
- Data source(s): Schulte-Ladbeck & Hopp (1998) Table 4: B-R(fit) = 0.89
- SPS model: BC03, Chabrier IMF, extended/quasi-continuous SFH
- Metallicity: Direct T_e — Berg et al. (2012) Table 4: 12+log(O/H) = 7.87 ± 0.02
- Dust: b ≈ +51° → E(B-V) ≲ 0.02 (negligible)
- Method validated: Single use
- Known issues: No explicit BC03 grid inversion shown; age inferred from physics direction (low Z + red → older)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| DDO 64 | 6.2 (±2) Gyr | T2 | MEDIUM | No | — |

---

### Method 5: Dual-color (FUV-NUV + B-V) → BC03 continuous SFH → t₅₀ [Spectroscopic metallicity]
- Data source(s): Gil de Paz et al. (2007) Table 3 (GALEX FUV/NUV) + Hunter & Elmegreen (2006) Table 3 via VizieR (B-V₀)
- SPS model: BC03, Chabrier IMF, extended/continuous SFH
- Metallicity: Spectroscopic — Kennicutt & Skillman (2001): 12+log(O/H) = 7.67 ± 0.05
- Dust: Foreground corrected; FUV-NUV ≈ 0 confirms low internal dust
- Method validated: Single use
- Known issues: SFH shape degeneracy (cannot distinguish T = 8 vs 12 Gyr); bursty SFHs can mimic continuous colors

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| DDO 154 | 5.0 (±1.5) Gyr | T2 | MEDIUM-LOW | No | — |

---

### Method 10: Multi-color "SED Lock" (de Blok+1995, 4 colors) → BC03 extended SFH → t₅₀ [Per-object strong-line metallicity]
- Data source(s): de Blok, van der Hulst & Bothun (1995) Table 4: U-B, B-V, B-R, V-I (luminosity-weighted and area-weighted)
- SPS model: BC03, Chabrier IMF, extended/continuous SFH
- Metallicity: Per-object strong-line — Kuzio de Naray et al. (2004) Table 6 (F561-1); van den Hoek+2000 Table 1 (F568-V1)
- Dust: MW foreground corrected; internal not corrected (LSB low-dust prior)
- Method validated: Yes (4 applications)
- Known issues: BC03 grid interpolation not directly reproduced; strong-line metallicity, not direct T_e
- Notes: Area-weighted primary for F561-1/F565-V2/F568-V1; luminosity-weighted primary for F574-2

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F561-1 | 5.2 Gyr (range 4.0–7.0) | T2 | MEDIUM | No | — |
| F565-V2 | ~5 Gyr (range 4–7) | T2 | MEDIUM-LOW | No | — |
| F568-V1 ⏳ | 5.7 (+1.1/−0.9) Gyr | T2 | MEDIUM | No | — |
| F574-2 | 6.0 (+1.0/−1.0) Gyr | T2 | MEDIUM | No | — |

⏳ F568-V1 PENDING Deep Think review.

---

### Method 11: V-I color (Pildis+1997) → explicit BC03 τ-model interpolation → t₅₀ [Direct T_e metallicity]
- Data source(s): Pildis+1997 Table 2: V-I = 0.82 for D563-4 (= F563-1)
- SPS model: BC03 Padova 1994, Chabrier IMF, τ = 5 Gyr fiducial (bracket 3–5)
- Metallicity: Direct T_e — de Blok & van der Hulst (1998): 12+log(O/H) = 8.19 direct, 8.01 ± 0.06 strong-line. Z ≈ 0.004.
- Explicit BC03 math shown and verified by Deep Think ("Platinum Standard for Traceability")
- Sensitivity grid: Z = {0.004, 0.008} × τ = {3, 5} → t₅₀ range 4.09–5.67 Gyr
- Method validated: Distinct from Method 2 (different Z source, different τ, explicit grid reproduction)
- Known issues: Single-color degeneracy; gas-phase vs stellar metallicity

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F563-1 | 5.6 Gyr (range 4.0–7.5) | T2 | MEDIUM | No | — |

---

### Method 12: B-V color (LSB photometry catalogs) → BC03 mapping with LSB metallicity correction → t₅₀
- Data source(s): Various — de Blok & McGaugh (1997) Table 1; McGaugh & Bothun (1994) Table 2; McGaugh & de Blok (1997) Table 1; van den Hoek+2000 Table 4; de Blok+1996 Appendix A2.4
- SPS model: BC03 color mapping with metallicity correction (metal-poor LSBs bluer at same age)
- Metallicity: Assumed ~0.1 Z☉ except F571-V1 (per-object van den Hoek+2000)
- Method validated: Yes (5 applications)
- Known issues: No explicit BC03 grid interpolation; metallicity mostly assumed; relies partly on session "LSB Color Ladder"
- Corrections: F567-2 column misread (0.54 → 0.67); F583-1 Bell value rejected (0.33 → 0.50 from B-R)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F563-V1 | 6.5 Gyr | T2 | MEDIUM | No | — |
| F563-V2 | 6.0 Gyr | T2 | MEDIUM | No | — |
| F567-2 | 7.2 Gyr (range 5.8–8.8) | T2 (borderline T3) | MEDIUM-LOW | No | — |
| F571-V1 🔄 | 5.5 (+1.3/−1.0) Gyr | T2 | MEDIUM-LOW | No | — |
| F583-1 | 5.0 Gyr (range 4.5–5.5) | T2 | MEDIUM | No | — |

🔄 F571-V1 PENDING Deep Think review.

---

### Method 15: Schombert+2011 color + birthrate parameter + gas fraction → t₅₀
- Data source(s): Schombert, Maciel & McGaugh (2011) Table 1 (B-V = 0.48), Table 2 (b = 0.14, f_g = 0.80)
- Schombert's "b" is log-scaled sSFR (b = 0.14 → dimensionless b ≈ 1.9)
- Method validated: Single use
- Known issues: No per-object metallicity; inclination may redden B-V
- Key finding: Gas fraction (80%) does NOT force T1. t_dbl > 4 Gyr = "simmering" (T2). Protocol 19.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F574-1 | 5.5 Gyr | T2 (Gas-Rich LSB) | HIGH (tier) / MEDIUM (age) | No | — |

---

### Method 16: BVRK per-object photometry (Bell+2000) → multi-color BC03 with B-K constraint → t₅₀
- Data source(s): Bell et al. (2000) Table 7: B, V, R, K₀ with uncertainties. B-K = 3.61 ± 0.13.
- SPS model: BC03, Chabrier IMF, extended SFH, sub-solar Z prior
- Method validated: Single use
- Known issues: Monotonicity enforcement — raw 6.8 Gyr revised to 7.5 for sequence consistency
- Correction: Astra used subset-average B-V ≈ 0.56; corrected to per-object 0.70

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F579-V1 | 7.5 Gyr (range 6.5–8.5) | T2 (borderline T3) | MEDIUM-LOW | No | — |

---

### Data Gap: F583-4
- No measured multi-band color, no UV, no SFR. B-V = "*" (missing). Assumed B-R unstable (0.78 vs 0.9). Coordinate ambiguity (~30 arcmin).

---

### Method 19: Single B-V color (FIGGS/Begum+2008) → BC03 with assumed Z → t₅₀ [Quality Floor]
- Data source(s): Begum et al. (2008) FIGGS Table 1: B-V = 0.37
- SPS model: BC03, Chabrier IMF, extended SFH, sub-solar Z prior
- Method validated: Single use
- Known issues: Bronze quality — no error bars, no metallicity, single color
- Key ruling: Central 3.6 Gyr is nominally T1 but "Quality Floor" prevents T1 lock on Bronze data

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| KK 251 | 3.6 Gyr (range 2.5–5.0) | T2 (borderline T1) | LOW | No | — |

---

### Method 32: GALEX UV (FUV-NUV) + optical B-V (Gil de Paz+2007) → BC03 color mapping → t₅₀
- Data source(s): Gil de Paz et al. (2007) Table 3: FUV-NUV = 0.43–0.46. Table 4: B-V = 0.67 ± 0.14.
- SPS model: BC03 color mapping (UV + optical)
- Method validated: Single use
- Known issues: B-V uncertainty ±0.14 is large; dust reddening biases older; proxy-based age
- Deep Think: "Hotspot Spiral" — bar instability rules out T1, UV slope rules out T3.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2903 | 6.2 Gyr | T2 (Main Sequence) | MEDIUM | No | — |

---

### Method 36: B-V color (McGaugh & Schombert 2014) + IRAS SFR (Seigar 2005) → b + Downsizing penalty → t₅₀
- Data source(s): SPARC Table 1. McGaugh & Schombert (2014): B-V = 0.45. Seigar (2005): S_60μm = 1.58 Jy, S_100μm = 4.63 Jy, SFR = 4.39 ± 0.89 M☉/yr (Salpeter).
- Derived: log M★ ≈ 10.88, b ≈ 0.50–0.80 (Kroupa/Salpeter), gas fraction 31%
- Method validated: Single use
- Known issues: FIR-only misses warm dust (b ≈ 0.50 is lower limit; true b likely 0.7–0.8). Downsizing penalty +1.3 Gyr applied for mass (Δlog M★ ≈ 0.4 vs NGC 3726).
- "Massive Steady State" anchor — blue (0.45) but older than lower-mass analogs due to earlier core formation.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2998 | 6.8 Gyr | T2 (Massive Steady State) | HIGH | No | — |

---

### Method 39: B-V color (Conselice+2000 Table 2, RC3-based) → BC03 mapping → t₅₀
- Data source(s): Conselice, Bershady & Jangren (2000) Table 2: B-V = 0.48. SAB(r)c at 17.0 Mpc.
- SPS model: BC03 color mapping; b ≈ 1.2–1.5 inferred from blue color
- Method validated: Single use
- Known issues: Single-color degeneracy; no sSFR cross-check retrieved. Blue color argues against significant dust.
- "Blue Spiral Anchor" — contrasts with NGC 2903 "Dusty Anchor" (same activity, different dust/morphology).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3726 | 5.5 Gyr | T2 (Vigorous Spiral) | MEDIUM | No | — |

---

### Method 40: Tully+1996 Table 5 global multi-color (B-R, B-I, B-K') + Pak+2014 GALEX UV → BC03 → t₅₀
- Data source(s): Tully et al. (1996) AJ 112, 2471 Table 5 (top line = global): B-R = 1.01, B-I = 1.49, B-K' = 3.10. Pak et al. (2014): NUV-r ≈ 2.85, FUV-NUV ≈ 0.42.
- Method validated: Single use
- Known issues: Interacting system (Arp 280 with NGC 3769A) — interactions may bias colors younger. No resolved CMD. Astra used bottom line (disk/μ₀ colors) instead of top line (global) — corrected.
- B-K' = 3.10 is in T2 regime (T3 requires >4.0 per "Infrared Chasm").

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3769 | 6.0 Gyr (range 4.5–7.5) | T2 | MEDIUM | No | — |

---

### Method 41: Bouquin+2018 GALEX UV + RC3 B-V + 2MASS B-K → BC03 → t₅₀ [Edge-on dust caveat]
- Data source(s): Bouquin et al. (2018) ApJS 234, 18: FUV-NUV = 0.72 ± 0.014. RC3: B-V = 0.63. 2MASS LGA (Jarrett+2003): B-K ≈ 3.19.
- Method validated: Single use
- Known issues: **Edge-on** (Sc) — dust preferentially reddens UV, making FUV-NUV unreliable as age indicator. UV/optical discordance (red UV, intermediate optical) = classic dust reddening signature. B-K = 3.19 confirms T2 regime. 7.3 Gyr is an **upper bound** (dust creates "false olds").
- Citation corrected: Astra attributed to "Bianchi et al. 2018" — correct author is Bouquin.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3877 | 7.3 Gyr (range 5.5–8.0) | T2 (upper bound, dust) | MEDIUM | No | — |

---

### Method 42: Tully+1996 Table 5 (B-R, B-K') + Hernández-Toledo+2001 Table 3 (B-V) + Bouquin+2018 UV → t₅₀ [Phantom data correction]
- Data source(s): Tully+1996 Table 5: B-R = 0.94, B-K' = 3.21. Hernández-Toledo & Puerari (2001) A&A 379, 546 Table 3: (B-V)_oT = 0.51. Bouquin+2018: FUV-NUV = 0.50.
- Method validated: Single use
- Known issues: **Original T3/8.0 REJECTED** — spreadsheet claimed B-V = 0.64, B-R = 1.13, neither reproducible from any citable source ("Phantom Data" protocol). Verified colors are significantly bluer. Interacting with NGC 3896 (KPG 302). B-K' = 3.21 confirms T2 regime.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3893 | 6.5 Gyr (range 5.0–7.5) | T2 | MEDIUM | No | — |

---

### Method 44: Kassin+2006 Table 3 integrated B, R, K photometry → B-R, B-K → BC03 → t₅₀
- Data source(s): Kassin, de Jong & Pogge (2006) ApJS 162, 80 Table 3: B = 11.65 ± 0.05, R = 10.77 ± 0.05, K = 8.47. B-R = 0.88 ± 0.07, B-K = 3.18.
- Method validated: Single use
- Known issues: Table 3 footnote says uncorrected for Galactic extinction (high latitude b ≈ +66° → minimal). Correction would make colors bluer → even younger. SAURON central spectroscopy exists but doesn't provide galaxy-wide t₅₀.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3949 | 5.4 Gyr (range 4.5–6.5) | T2 | MEDIUM | No | — |

---

### Method 45: Tully+1996 Table 5 global colors (B-R, B-K') → Valley Rule → t₅₀
- Data source(s): Tully+1996 Table 5 (global/top line): B-R = 1.25, B-I = 1.85, B-K' = 3.73.
- Method: **Valley Rule** (new protocol): 3.3 < B-K < 4.0 AND B-R > 1.15 → T3. B-K' = 3.73 exceeds LSB Control (3.61 → 7.5 Gyr). Massive Spiral Prior → 8.8 Gyr.
- Method validated: Single use
- Known issues: Astra used disk/μ₀ colors (B-R=1.35, B-K'=3.96) instead of global — corrected. No UV cross-check. First galaxy in the "Valley" (B-K 3.3–3.9).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3953 | 8.8 Gyr (range 8.0–9.5) | T3 (Valley Rule) | MEDIUM-HIGH | No | — |

---

### Method 47: Kassin+2006 + Tully+1996 multi-color + NED UV → Soft Valley Protocol → t₅₀
- Data source(s): Kassin+2006 Table 3: B = 10.66, R = 9.45 → B-R = 1.21. Tully+1996 Table 5: B-R = 1.21, B-K' = 3.42 (cross-source agreement). NED SED: UV_1650 = 12.126, UV_2400 = 11.533 → UV color = 0.59.
- Method: **Soft Valley Protocol** (new): B-K 3.4–3.8 = Valley Zone → second axis mandatory. UV color 0.59 < 0.9 → Active → T2. Red B-R = 1.21 attributed to bar/bulge dominance, not ancient mass.
- Method validated: Single use
- Known issues: Original T3/8.5 **REJECTED**. Soft Valley Protocol established via Deep Think + Pro collaboration after methodology dispute. Red optical color anomalous for T2 but explained by bar structure.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3992 (M109) | 7.2 Gyr (range 6.5–7.8) | T2 (Soft Valley, UV-confirmed) | MEDIUM-HIGH | No | — |

---

### Method 48: Tully+1996 Table 5 edge-on corrected colors (B-R, B-K') → Soft Valley Protocol → t₅₀ [Differential Correction pair]
- Data source(s): Tully+1996 Table 5 (inclination/extinction-corrected magnitudes). NGC 4010: B-R = 0.96, B-K' = 3.37 (raw Table 4: B-K' = 4.14, Δ = −0.77 mag). NGC 4013: B-R = 1.39, B-K' = 3.99. Both edge-on (~88–90°).
- Method: Table 5 corrections applied for edge-on dust. NGC 4010 B-K' = 3.37 < 3.4 → Strong T2 (no second axis). NGC 4013 B-K' = 3.99 > 3.8 → Strong T3 (no second axis). R-K = 2.60 for NGC 4013 confirms evolved backbone.
- Method validated: Yes (2 applications — Differential Correction pair)
- Known issues: NGC 4010 Astra used raw Table 4 (would have been T3 at B-K' = 4.14). Paper explicitly warns of dust patches. NGC 4013 correction validates Tully+1996 method — doesn't force all edge-on galaxies to same bin. NGC 4013 has reddest B-R (1.39) of any verified spiral.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4010 | 6.0 Gyr (range 5.0–7.0) | T2 | MEDIUM | No | — |
| NGC 4013 | 9.0 Gyr (range 8.0–10.0) | T3 | MEDIUM | No | — |

---

### Method 51: Tully+1996 Table 5 global colors (B-R, B-I, B-K') → Soft Valley Protocol → t₅₀ [Standard Ursa Major spirals]
- Data source(s): Tully+1996 Table 5 (inclination/extinction-corrected). All use global (top line) colors.
- Method: Apply Soft Valley Protocol. B-K' < 3.4 → Strong T2. B-K' > 3.8 → Strong T3. R-K backbone diagnostic as confirmation.
- Method validated: Yes (4 applications)
- Known issues: NGC 4085 interacting with NGC 4088 (caveat only). NGC 4100 = procedural milestone (Astra correctly avoided Line 2 Trap). NGC 4157 borderline T3 (B-K' = 3.81, just +0.01 above threshold) but R-K = 2.56 confirms evolved backbone. NGC 4183 = Blue Anchor (B-K' = 2.43, bluest spiral in pipeline; R-K = 1.72). Edge-on validation suite spans B-K' 2.43–3.81.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4085 | 6.6 Gyr (range 5.5–7.5) | T2 | MEDIUM | No | — |
| NGC 4100 | 7.0 Gyr (range 6.0–8.0) | T2 | MEDIUM | No | — |
| NGC 4157 | 8.4 Gyr (range 7.5–9.5) | T3 (borderline, R-K confirmed) | MEDIUM | No | — |
| NGC 4183 | 5.2 Gyr (range 4.0–6.5) | T2 (Blue Anchor) | MEDIUM | No | — |

---

### Method 52: Tully+1996 Table 5 global colors + Bouquin+2018 GALEX UV → Soft Valley Protocol → t₅₀ [disk/μ₀ correction]
- Data source(s): Tully+1996 Table 5 global: B-R = 1.04, B-K' = 3.32. Bouquin+2018 GALEX: FUV = 13.82, NUV = 13.21, FUV-NUV = 0.61.
- Method: Astra used μ₀ disk-central colors (B-K' = 3.58 → Valley Zone). Corrected to global (B-K' = 3.32 → Strong T2). UV confirms active (FUV-NUV = 0.61 < 0.9). Fourth Tully two-line trap instance.
- Method validated: Single use
- Known issues: Original T3/9.0 **REJECTED**. Interacting pair with NGC 4085 — both now T2, physically coherent ("Wet Interaction"). Deep Think issued SYSTEMIC ALERT: audit all Tully+1996 entries for Line 2 Trap.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4088 (Arp 18) | 6.8 Gyr (range 5.5–7.5) | T2 | MEDIUM | No | — |

---

### Method 54: Tully+1996 Table 5 Valley Zone + Backbone Tiebreaker Protocol → t₅₀ [UV catalog conflict]
- Data source(s): Tully+1996 Table 5 global: B-R = 1.27, B-K' = 3.77 (Valley Zone). Bouquin+2018: FUV-NUV = 0.77 (Active). Bianchi+2015: FUV-NUV = 0.95 (Passive). R-K = 2.50.
- Method: Valley Zone (3.4–3.8) triggers second axis. UV conflicted — Bouquin says T2, Bianchi says T3. **Backbone Tiebreaker Protocol** (new): In edge-on systems (i > 75°), when UV conflicts with backbone structure, Structure Wins. R-K = 2.50 matches Fossil T3 spirals (NGC 3953: 2.48, NGC 4013: 2.60). A fortiori: NGC 4217 redder than confirmed T3 NGC 3953 in both optical and infrared.
- Method validated: Single use
- Known issues: Borderline Valley Zone (3.77, near 3.8 threshold). UV catalog discrepancy of 0.18 mag straddles the 0.9 threshold. Edge-on (i = 80°) drives UV uncertainty. New protocol: "When UV conflicts with Structure, Structure Wins."

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4217 | 9.0 Gyr (range 8.0–10.0) | T3 (Backbone Tiebreaker) | HIGH | No | — |

---

### Method 56: Gil de Paz+2007 GALEX Atlas UV (Table 3) + optical B-V (Table 4) → BC03 SPS mapping → t₅₀
- Data source(s): Gil de Paz et al. (2007) ApJS 173, 185. Table 3: FUV = 12.07±0.01, NUV = 11.87±0.01 → FUV-NUV = 0.20±0.01. Table 4: B = 10.46±0.11, V = 10.01±0.12 → B-V = 0.45±0.16.
- Method: Map {FUV-NUV ≈ 0.2, B-V ≈ 0.45} onto BC03 tracks → mid-sequence "Active Spiral Equilibrium" → t₅₀ ≈ 6.3 Gyr. Hyper-Active UV signal (0.20 << 0.9 threshold). Table 3 (asymptotic) prioritized over Table 4 (aperture) per Deep Think source protocol.
- Method validated: Single use
- Known issues: Astra's B-V ≈ 0.49 corrected to 0.45 (table column mixup — 0.49 was IRAS 12μm flux). No direct CMD t₅₀. "Face-On Blue Anchor" for Phase II.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4559 | 6.3 Gyr (range 5.0–7.5) | T2 (Hyper-Active Control) | HIGH | No | — |

---

### Method 62: Bouquin+2018 GALEX asymptotic UV + RC3 integrated B-V → BC03 + Downsizing → t₅₀ [Chain-of-custody correction]
- Data source(s): Bouquin+2018: FUV = 13.01±0.01, NUV = 12.78±0.01 → FUV-NUV = 0.23. RC3: B(T) = 11.20±0.14, V(T) = 10.74±0.14 → B-V = 0.46.
- Method: FUV-NUV = 0.23 << 0.9 → Very Active → Strong T2. Cosmic Downsizing: bulgeless Sd morphology → younger mass-weighted age than bulged Scd at same color (NGC 4559 at 6.3 Gyr has bulge; NGC 5585 bulgeless → 5.0 Gyr). **Chain-of-custody correction:** Astra's original source (VizieR J/AZh/88/342) contained region-by-region photometry of individual H II knots, not integrated galaxy colors — "Aperture/Granularity Error" (mistaking part for whole).
- Method validated: Single use
- Known issues: Original data source **REJECTED** (region-by-region, not integrated). Correct tier but wrong provenance. Astra's B-V ≈ 0.49 replaced by RC3 B-V = 0.46.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5585 | 5.0 Gyr (range 4.0–6.0) | T2 (Hyper-Active) | MEDIUM | No | — |

---

### Method 63: Bouquin+2018 GALEX asymptotic UV + RC3 extinction-corrected optical → Edge-On Dust Trap resolution → t₅₀
- Data source(s): Bouquin+2018: FUV = 14.06±0.01, NUV = 13.63±0.01 → FUV-NUV = 0.43. RC3: B(T)° = 9.70, V(T)° = 9.18 → (B-V)° = 0.52.
- Method: Astra's T3/10 Gyr based on uncited, uncorrected B-R = 1.41 **REJECTED** — classic "Edge-On Dust Trap" (i ≈ 87°). UV pierces dust: FUV-NUV = 0.43 < 0.9 → Active. Corrected (B-V)° = 0.52 confirms blue stellar population. **Edge-On Dust Trap principle formalized:** For i > 80°, raw optical colors measure dust column density, not stellar age. UV binary check: if FUV-NUV < 0.9, galaxy is active regardless of optical appearance. UV "Signal" (OB stars) overpowers "Noise" (dust attenuation).
- Method validated: Single use
- Known issues: Original T3/10 Gyr **REJECTED** — uncited, uncorrected optical. Seventh T3→T2 correction. B-R = 1.41 has no chain-of-custody.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5907 | 7.0 Gyr (range 6.0–8.0) | T2 (Edge-On Dust Trap resolved) | MEDIUM | No | — |

---

### Method 65: Hernández-Toledo+2007 Table 3 (observed BVRI) + Table 4 (corrected colors) → BC03 → t₅₀ [Conservative mass-weighted anchor]
- Data source(s): Hernández-Toledo et al. (2007) AJ 134, 2286. Table 3 (observed): B-V = 0.56, B-R = 1.08, B-I = 1.61. Table 4 (extinction-corrected): (B-V)_c = 0.41, (B-R)_c = 0.94, (B-I)_c = 1.33. Identity: CIG 710 = NGC 6015 (Table 2).
- Method: Observed B-V = 0.56 → BC03 → t₅₀ ≈ 6.0 Gyr. Corrected (B-V)_c = 0.41 is significantly bluer (among bluest in sample), would suggest younger age or even T1 borders. Conservative estimate retained at 6.0 Gyr to account for "Skin vs Skeleton" — blue disk (skin) vs older bulge (skeleton) in Scd morphology. **Principle: Always prefer extinction-corrected colors; use morphology to interpret mass-weighted age.**
- Method validated: Single use
- Known issues: Corrected colors suggest galaxy could be younger than 6.0 Gyr. Conservative approach prevents "age deflation" into T1 from blue disk. Isolated galaxy (CIG catalog).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6015 | 6.0 Gyr (range 5.0–7.0) | T2 (Conservative Anchor) | MEDIUM | No | — |

---

### Method 68: McGaugh & Schombert (2014) Table 1 B-V + J-Ks → BC03 → t₅₀ [SPARC source elevated]
- Data source(s): McGaugh & Schombert (2014) AJ 148, 77. Table 1: B-V = 0.57, J-Ks = 0.86.
- Method: B-V = 0.57 → BC03 → t₅₀ ≈ 6.2 Gyr. Comfortably T2 (T3 would need B-V ≳ 0.75). J-Ks = 0.86 consistent with moderately evolved spiral. **Aggregator B-V ≈ 0.84 REJECTED** — likely aperture bias (red bulge dominated) or correction-state mismatch. **McGaugh & Schombert (2014) elevated to Tier 1 (High Reliability)** for all SPARC galaxies — meticulous surface photometry for stellar mass derivation.
- Method validated: Single use
- Known issues: Aggregator discrepancy Δ ≈ 0.27 mag would shift galaxy from Blue Cloud to Red Sequence — demonstrates critical importance of source selection. SBb bar/bulge accounts for age anchoring at 6.2 rather than younger.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6674 | 6.2 Gyr (range 5.0–7.0) | T2 (Control) | HIGH | No | — |

---

### Method 72: Gil de Paz+2007 Table 3 GALEX UV asymptotic → FUV-NUV + Blue Survival Principle → t₅₀
- Data source(s): Gil de Paz+2007 ApJS 173, 185. Table 3: FUV = 13.42±0.01, NUV = 12.72±0.01 → FUV-NUV = 0.70. D25 check: FUV-NUV = 0.69 (consistent). Table 4: B = 10.35, V = 9.48 → B-V = 0.87.
- Method: FUV-NUV = 0.70 < 0.9 → Active → T2. **Blue Survival Principle:** dust extinction is wavelength-dependent (A_FUV > A_NUV > A_V), so dust always reddens. Since reddened UV still clears Active threshold, intrinsic SF must be even more vigorous. T3 is physically impossible — a fossil population would show FUV-NUV > 1.5, not 0.70. Red optical B-V = 0.87 attributed to dust in inclined system (~70°) + massive Sb bulge. Structural analog to NGC 5907 (Edge-On Dust Trap) but at moderate inclination.
- Method validated: Single use (second Gil de Paz Table 3 application after NGC 4559/M56)
- Known issues: B-V = 0.87 would suggest T3 if used alone — UV overrides. Galactic extinction corrected but not internal dust. "Milky Way's twin."

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 7331 | 7.0 Gyr (range 6.0–8.0) | T2 (Late T2 / Blue Survival) | HIGH | No | — |

---

### Method 73: Gil de Paz+2007 Table 3 GALEX UV asymptotic + resolved RGB context (GHOSTS) → t₅₀ [Hyper-Active Blue Limit]
- Data source(s): Gil de Paz+2007 Table 3: FUV = 11.17±0.01, NUV = 11.00±0.01 → FUV-NUV = 0.17. D25 check: FUV-NUV = 0.18 (consistent). Table 4: B = 9.63. Resolved-star context: Radburn-Smith+2012 (GHOSTS), Vlajić+2011 — well-populated RGB/old populations detected in outer disk.
- Method: FUV-NUV = 0.17 = **bluest UV in entire sample** (displaces NGC 4559 at 0.20). Extended/near-continuous SFH → mass-weighted t₅₀ naturally lands in mid-Gyr range. **"Skin vs Skeleton" principle:** UV = 0.17 measures active <100 Myr "skin," but HST-resolved RGB populations prove ancient "skeleton" exists. For galaxy with verified ancient substrate, mass-weighted age cannot drop below constant-SFH equilibrium (~6.5 Gyr). **Hyper-Active Triad established:** NGC 4183 (B-K' = 2.43, IR Blue Limit), NGC 4559 (FUV-NUV = 0.20, UV Control), NGC 7793 (FUV-NUV = 0.17, Absolute UV Limit) — defines "Zero Point" of Active Spiral Sequence.
- Method validated: Single use
- Known issues: Bluer than NGC 5585 (0.23) yet older (6.5 vs 5.0) — resolved by Cosmic Downsizing (NGC 5585 is bulgeless dwarf; NGC 7793 has ancient halo/thick disk). PLATINUM confidence from Deep Think.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 7793 | 6.5 Gyr (range 5.5–7.5) | T2 (Hyper-Active Blue Limit) | PLATINUM | No | — |

---

### Method 76: de Blok+1995 Table 4 luminosity-weighted B-V + van der Hulst+1993 cross-check → BC03 → t₅₀ [LSB Transparency Principle]
- Data source(s): de Blok, van der Hulst & Bothun (1995) MNRAS 274, 235. Table 4: B-V(luminosity-weighted) = 0.60, B-V(area-weighted) = 0.49. Cross-check: van der Hulst et al. (1993) AJ 106, 548: B-V(area) = 0.49±0.05 — exact match.
- Method: B-V(lum) = 0.60 → BC03 → t₅₀ ≈ 6.0 Gyr. **LSB Transparency Principle:** LSB galaxies have low H I surface density, low metallicity (~0.2–0.4 Z☉), low dust → "optically transparent" → optical colors are reliable age proxies without heavy extinction corrections. **Luminosity-weighting priority:** stellar mass follows luminosity more closely than geometric area in exponential disks; area-weighting oversamples blue "ethereal skin." Both weightings give T2.  "Slow Burn" evolution — LSBs evolve inefficiently, naturally intermediate age. Dual-source verification (two independent papers agree).
- Method validated: Single use (first LSB protocol application)
- Known issues: LSB low-dust assumption is working approximation. Both B-V values (0.60, 0.49) support T2 — tier is robust.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 00128 | 6.0 Gyr (range 5.0–7.0) | T2 (LSB Control) | HIGH | No | — |

---

### Method 80: Gavilán+2013 Table 3 B-V + SFR + SPARC mass → b/t_dbl T1/T2 discriminator → t₅₀ [Gas-Rich Steady State]
- Data source(s): Gavilán et al. (2013) Table 3: B-V = 0.46±0.03, SFR = 0.0100 M☉/yr, log M★ = 8.25. Cross-check: van Zee et al. B-V = 0.46±0.03 (exact match). SPARC: L[3.6] = 0.374×10⁹ L☉, M_HI = 0.428×10⁹ M☉, D = 10.20 Mpc. M★ ≈ 1.87×10⁸ M☉. Gas fraction ~229%.
- Method: b = 0.78 (near steady-state), t_dbl = 17.8 Gyr (>> T_universe). **KEY TEST CASE:** 229% gas but NOT T1 because t_dbl = 17.8 >> 4 Gyr. T1 requires BOTH high gas AND rapid assembly (short t_dbl). High gas + slow consumption = "Gas-Rich Steady State" → T2. Aligned with NGC 0024 anchor (b ≈ 0.70, 6.5 Gyr). Same B-V = 0.46 as UGC 00191 but different b (0.78 vs 1.66) → different ages (6.5 vs 5.0 Gyr). **Birthrate is finer clock when colors identical.**
- Method validated: Single use (second Gavilán+2013 application after UGC 00191)
- Known issues: Pro initially proposed 7.7 Gyr, Deep Think corrected to 6.5 based on b ≈ 0.78 near steady-state.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 00891 | 6.5 Gyr | T2 (Gas-Rich Steady State) | HIGH | No | — |

---

### Method 81: van der Hulst+1993 B-V + LSB Stability Block + Metallicity Trap correction → t₅₀
- Data source(s): van der Hulst et al. (1993) AJ 106, 548: B-V = 0.42±0.03 (mean disk), U-B = −0.18±0.03.
- Method: Astra mapped B-V = 0.42 → BC03 → 4.9 Gyr. **RECLASSIFIED to 6.0 Gyr.** Created "Redder is Younger" inversion with UGC 00634 (B-V = 0.39 → 6.0 Gyr) — physically inconsistent. **"Metallicity Trap":** LSB galaxies are metal-poor (Z ~ 0.1–0.3 Z☉); metal-poor stars burn hotter/bluer → 6 Gyr population at low Z exhibits B-V ≈ 0.40–0.45. Naive solar-metallicity BC03 mapping underestimates age by ~1.0–1.5 Gyr. **"LSB Stability Block":** Standardize blue LSBs to 6.0 Gyr = "Slow Burn" equilibrium (constant SFH over Hubble time → t₅₀ ≈ 12/2). **T1 Anomaly Filter for LSB/dwarfs:** to claim T1 or early T2, need explicit b >> 1 or FUV-NUV < 0.2.
- Method validated: Single use (third LSB application after UGC 00128, related to de Blok+1995)
- Known issues: 4.9→6.0 Gyr reclassification. Tier unchanged (T2 throughout).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 01230 | 6.0 Gyr (range 5.5–7.0) | T2 (LSB Stability Block) | HIGH | No | — |

---

### Method 84: Gil de Paz+2007 GALEX UV + GALEX Table 4 optical + Kaisin+2011 Hα SFR + SPARC mass → Fuel Limit discriminator → t₅₀ [Starburst Anchor]
- Data source(s): Gil de Paz+2007 Table 3: FUV = 12.69±0.01, NUV = 12.48±0.01 → FUV-NUV = 0.21. Table 4: B-V = 0.58±0.18. Kaisin et al. (2011) Table 3: SFR(Hα) = 0.43 M☉/yr, SFR(FUV) = 0.33 M☉/yr. SPARC: L[3.6] = 3.649×10⁹ L☉, M_HI = 0.803×10⁹ M☉, D = 6.92 Mpc. M★ ≈ 1.82×10⁹ M☉. Gas fraction ~44%.
- Method: b ≈ 2.8 (highest in session!), t_dbl ≈ 4.2 Gyr. **NOT T1 despite extreme activity.** **"Fuel Limit Discriminator":** M_gas/M★ = 0.44 — galaxy cannot double stellar mass with current reserves. If M_gas < M★, galaxy is evolved + experiencing rejuvenation, not primitive/embryonic. B-V = 0.58 (intermediate) confirms old stellar backbone. Contrast: PGC 51017 (M_gas > 2.5×M★, can double → T1) vs NGC 1156 (M_gas < M★, cannot double → T2). "Middle-aged galaxy having a mid-life crisis." Age pulled younger than steady-state (6.0) by high b, but backbone prevents reaching primitive regime (5.0).
- Method validated: Single use
- Known issues: B-V = 0.58 has large uncertainty (±0.18). t_dbl = 4.2 borderline near 4.0 threshold but Fuel Limit resolves.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 1156 (UGC 02455) | 5.5 Gyr | T2 (Starburst Anchor) | HIGH | No | — |

---

### Method 92: Gil de Paz+2007 GALEX Table 3 UV + Table 4 optical + Edge-On Bifurcation + LSB Stability Block → t₅₀ [Transparent Edge-On]
- Data source(s): Gil de Paz+2007 Table 3: FUV = 14.87±0.01, NUV = 14.50±0.01 → FUV-NUV = 0.38 (asymptotic). D25: FUV-NUV = 0.35. Table 4: B = 13.07±0.15, V = 12.63±0.16 → B-V = 0.44. HyperLEDA B-V = 0.35 **SUPERSEDED** by GALEX Atlas.
- Method: FUV-NUV = 0.38 < 0.9 → Active → T2. B-V = 0.44 despite edge-on → remarkably blue → **"Transparent Edge-On"** (low dust/metallicity). **"Edge-On Bifurcation" protocol established:** (1) Dusty/Evolved "Red" Edge-Ons (B-V ≳ 0.8 or rejected; 6.5–7.0 Gyr; UGC 01281, NGC 7331, NGC 5907) vs (2) **Transparent/Primitive "Blue" Edge-Ons (B-V ≈ 0.4–0.5; 6.0 Gyr; UGC 04278).** Blue edge-on → low metallicity → LSB/dwarf physics → **LSB Stability Block at 6.0 Gyr.** Resolves "Redder is Younger" inversion with UGC 00634 (B-V = 0.39 → 6.0). Pro proposed 5.0 Gyr, reclassified to 6.0 by Deep Think.
- Method validated: Single use (fourth Gil de Paz+2007 application)
- Known issues: Aggregator B-V = 0.35 superseded. 5.0→6.0 Gyr reclassification for consistency.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 04278 | 6.0 Gyr (range 5.5–7.0) | T2 (Transparent Edge-On) | HIGH | No | — |

---

### Method 97: McGaugh & Schombert (2014) Table 1 B-V → LSB Stability Block + Metallicity Trap → t₅₀ [T1 REJECTED — Bluest Galaxy]
- Data source(s): McGaugh & Schombert (2014) AJ 148, 77 Table 1: B-V = 0.35 (extinction-corrected), D = 57.1 Mpc, M_V = −18.7, M_[3.6] = −21.18.
- Method: **B-V = 0.35 — BLUEST in entire sample.** Astra claimed T1/3.3 Gyr → **REJECTED.** Pro proposed 4.6 Gyr → **revised to 6.0 Gyr.** **"Metallicity Trap at Blue Limit":** at very low metallicity (Z ~ 0.1 Z☉), MS turnoff is hotter/bluer for given age. A 6 Gyr constant-SFH population at low Z naturally exhibits B-V ≈ 0.35. **"Extreme blueness is a thermometer measuring low opacity, not a clock measuring extreme youth."** Reductio: assigning 4.6 Gyr would make UGC 05005 youngest galaxy in dataset without independent burst signature. **LSB Stability Block universality affirmed:** B-V range 0.35–0.60 all map to 6.0 Gyr for LSB class. **To break block and claim <5 Gyr, must show:** Hyper-Active UV (FUV-NUV < 0.25) OR Starburst b >> 1. UGC 05005 lacks both.
- Method validated: Second McGaugh & Schombert 2014 application (after NGC 6674)
- Known issues: T1/3.3 Gyr rejected. Pro's 4.6 Gyr revised upward. No second-axis confirmation available.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05005 | 6.0 Gyr | T2 (LSB Stability Block) | HIGH | No | — |

---

### Method 110: Sprayberry+1995 B-V + Rahman+2007 SFR(TIR) → inferred sSFR → Giant LSB Bifurcation → t₅₀ [GLSB distinct from dwarf LSB]
- Data source(s): Sprayberry et al. (1995) AJ 109, 558: B-V = 0.72 (per-object, citing McGaugh & Bothun 1994). Rahman et al. (2007) arXiv:0704.1483: SFR(TIR) = 0.88 M☉/yr, M_B = −20.70. Inferred: L_B ≈ 3×10¹⁰ L☉, Υ_B ≈ 2.0, M★ ≈ 6×10¹⁰ M☉, log sSFR ≈ −10.8, τ_dbl ≈ 68 Gyr.
- Method: **Giant LSB Bifurcation established.** Dwarf LSBs (no bulge, B-V 0.35–0.60, 6.0 Gyr via Stability Block) are physically distinct from Giant LSBs (massive luminous bulge + diffuse disk, B-V ~0.72, T3 physics). Mass-weighted age dominated by ancient bulge, not the LSB disk. Inferred sSFR = Transition Zone (−10.5 to −11.0). Astra's 8.5 Gyr adjusted to 8.0 Gyr for boundary status — bulge-dominated but disk trickle persists.
- Method validated: Single use
- Known issues: MEDIUM confidence — no direct M★, no UV axis. sSFR inferred via "napkin check." Bulge/disk decomposition assumed from GLSB class membership.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06614 | 8.0 Gyr (range 7.5–9.0) | T3 (Boundary / Giant LSB) | MEDIUM | No | — |

---

### Method 114: Noordermeer & van der Hulst (2007) MNRAS 376, 1480 surface photometry B-R → Optical Lock → t₅₀ [Relative Terminology Trap — Astra T2 REJECTED]
- Data source(s): Noordermeer & van der Hulst (2007): B = 11.57±0.11, R = 10.17±0.26 (dedicated surface photometry, Platinum Source). Derived B-R = 1.40 (m_lim), 1.31 (m_25). Extinction-corrected: (B-R)_lim = 1.37, (B-R)_25 = 1.28. All methods give B-R ≈ 1.3–1.4. van Driel et al. (1994): qualitative "blue for Sa" — NOT quantitative. SA(s)ab morphology.
- Method: **Astra classified T2/7.0 Gyr based on "blue for Sa" — REJECTED as Relative Terminology Trap.** "Blue for Sa" means bluer than typical Sa (B-R ~1.5–1.6), but B-R = 1.40 is still absolutely RED on the Hermes Grid (T3 threshold: B-R > 1.3). Physics does not care about Hubble Type labels. **Optical Lock protocol:** for massive early-type spirals with B-R > 1.3, null hypothesis is T3; burden of proof on active diagnostics (UV required to claim T2). Matches T3 anchor NGC 2841.
- Method validated: Single use (first Noordermeer+2007 photometry application)
- Known issues: Astra's T2/7.0 Gyr rejected. Aggregator colors rejected in favor of Platinum Source surface photometry.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3898 (UGC 06787) | 9.0 Gyr (range 8.5–10.0) | T3 (Early-Type Spiral) | HIGH | No | — |

---

### Method 116: Moffat & Rahvar (2013) Table 2 + Zonoozi & Haghi (2010) Table 1 B-V → LSB Stability Block → t₅₀ [Dual-source UMa photometry]
- Data source(s): Moffat & Rahvar (2013) MNRAS 436, 1439, Table 2: B-V = 0.53 (extinction-corrected), A_B = 0.098. Zonoozi & Haghi (2010) A&A 524, A53, Table 1: B-V = 0.53 (independent confirmation). Both trace to Verheijen/Sanders-McGaugh UMa photometry lineage. SBd morphology, LSB (Tully et al. description). D ≈ 15.5–18.6 Mpc.
- Method: BC03 τ-model: B-V = 0.53 → t₅₀ ≈ 5–7 Gyr, modal 6.2 Gyr. **Adjusted to 6.0 Gyr per LSB Stability Block** — color variations (0.35–0.60) in LSBs are metallicity-driven, not age-driven. Astra's 6.2 Gyr would create "Bluer is Older" artifact vs UGC 00128 (B-V=0.60→6.0). Dual-source "Stereoscopic Verification" eliminates measurement error. LSB transparency (A_B ≈ 0.1) makes optical color sufficient without UV.
- Method validated: Single use (first Moffat & Rahvar application)
- Known issues: No UV/sSFR cross-check. Both sources share UMa lineage (not fully independent). No explicit photometric uncertainties reported.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06917 | 6.0 Gyr | T2 (LSB Stability Block) | HIGH | No | — |

---

### Method 118: REJECTED — UGC 06930 (Poison Data / Aperture Fallacy)
- Data source(s): Watson et al. (2012) ApJ 751, 123: log Σ_SFR = −2.63 M☉/yr/kpc² (21″ aperture), log M★ = 9.1, D = 17.0 Mpc. Excel sheet: sSFR ≈ 7.8×10⁻¹⁰ yr⁻¹ (log −9.11) → t₅₀ = 4.0 Gyr (T1).
- Method: **REJECTED as Poison Data.** Pro traced Excel sSFR to Watson+2012 Σ_SFR scaled by D₂₅ disk area — NOT a direct catalog value. Watson's Σ_SFR measured within 21″ central aperture (CO beam strategy), NOT integrated over galaxy. Scaling central aperture density to full disk (~138″ radius) overestimates global SFR by orders of magnitude (2.19 dex swing from aperture to disk scaling). Deep Think ruling: "Aperture Fallacy" — equivalent to measuring city center population density and multiplying by state area.
- **NEW PROTOCOL: Aperture Fallacy / Poison Data** — Any sSFR derived by scaling aperture Σ_SFR to D₂₅ is automatically rejected. Poison Data (verifiably wrong methodology) distinct from Ghost Data (unverifiable source).
- Resolution: Requires global photometry (GALEX Atlas or RC3 integrated colors). Expected reversion to T2 (~6.0 Gyr).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06930 | — (REJECTED) | POISON DATA | — | N/A | Aperture Fallacy — 4.0 Gyr T1 purged |

---

### Method 130: FIGGS (2008MNRAS.383..809B) + Makarova (1999) dual-source B-V → LSB Stability Block → t₅₀ [Twin Anchor Standardization]
- Data source(s): FIGGS compilation table: B-V = 0.59 (tabulated, not figure-read). Makarova (1999A&AS..139..491M): B-V = 0.59 (B=12.95, V=12.36). Cross-ID: UGC 07577 = DDO 125 = PGC 40904. CVn I cloud member, HI-rich dIrr.
- Method: Dual-source color agreement (FIGGS + Makarova) eliminates measurement ambiguity. Astra proposed 6.7 Gyr. **Deep Think standardized to 6.0 Gyr** — "Twin Anchor" check: accepting 6.7 for B-V=0.59 while UGC 00128 (B-V=0.60) sits at 6.0 creates "Bluer is Older" artifact. Δ=0.01 mag is within photometric noise and metallicity variation. LSB Stability Block absorbs B-V = 0.35–0.60 range.
- Method validated: Multiple LSB Stability Block applications
- Known issues: Astra's 6.7 Gyr overridden. No sSFR cross-check (not required for Stability Block).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07577 (DDO 125) | 6.0 Gyr | T2 (LSB Stability Block) | HIGH | No | — |

---

### Method 136: Larsen & Richtler (2000) Table 2 RC3 B-V = 0.41 → Mass-Sequence Alignment (Sm) → t₅₀ [Excel Error Corrected]
- Data source(s): Larsen & Richtler (2000) Table 2 "Integrated properties… mostly taken from RC3": B-V = 0.41, U-B = −0.33, T = 9.0 (RC3 late-type). Cross-ID: UGC 08490 = NGC 5204 = PGC 47368. SA(s)m Magellanic Spiral.
- Method: **Excel Error Corrected.** Excel claimed B-V = 0.59 → 6.7 Gyr — likely spreadsheet contamination from UGC 07577 (also B-V = 0.59). Pro verified B-V = 0.41 from RC3. Deep Think rejected Pro's 6.0 Gyr (Stability Block) — NGC 5204 is a low-mass Sm system, not an LSB dwarf. B-V = 0.41 (bluer) cannot be older than NGC 5585 (B-V = 0.46, 5.0 Gyr) without violating Cosmic Downsizing. **Mass-Sequence Alignment:** anchor to late-type sequence at 5.0 Gyr rather than LSB block.
- Method validated: Single use (first Larsen & Richtler application)
- Known issues: Excel value wrong. Pro's 6.0 Gyr overridden by Deep Think.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 08490 (NGC 5204) | 5.0 Gyr | T2 (Mass-Sequence) | MEDIUM | No | — |

---

### Method 137: Noordermeer & van der Hulst (2007) MNRAS 376, 1480, Table A3 B-R = 1.31±0.26 + S4G morphology (S0/a) → Morphology-Anchored T3 → t₅₀ [Dust Trap — Upper Limit]
- Data source(s): Noordermeer & van der Hulst (2007) Table A3: B(m₂₅) = 13.45±0.11, R(m₂₅) = 12.14±0.24. B-R = 1.31±0.26 (total within D₂₅ isophote). S4G classification: S0/a with ring/spiral qualifiers. Inclination i = 77° (Table A2). Authors flag "internal absorption by dust" — Galactic foreground corrected but NO internal dust correction. Cross-ID: UGC 08699 = NGC 5289 = PGC 48749.
- Method: **Morphology-Anchored classification.** Color B-R = 1.31 maps to ~8.8 Gyr BUT: (1) i = 77° creates Dust Trap — disk dust preferentially reddens B more than R, biasing age older; (2) ±0.26 uncertainty permits T1-level ages at blue edge (B-R = 1.05 → ~0.9 Gyr per Lanyon-Foster+2007); (3) no internal extinction correction applied. T3 maintained via morphology prior (S0/a is characteristically old), NOT color certainty. **Reject "Backbone Anchor" designation.** 8.8 Gyr tagged as upper limit.
- Method validated: First Noordermeer Table A3 application
- Known issues: MEDIUM-LOW confidence. Dust-contaminated. Not suitable as anchor. Needs independent proxy (GSWLC sSFR or GALEX NUV-r) for hardening.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 08699 (NGC 5289) | 8.8 Gyr | T3 (Morphology-Anchored, upper limit) | MEDIUM-LOW | No | — |

---

### Method 138: de los Reyes & Kennicutt (2019) ApJ 872, Table 1 RC3 B-V = 0.43±0.04 → False Precision + Interpolation Principle → t₅₀ [Sm/IBm Cohort Standardization]
- Data source(s): de los Reyes & Kennicutt (2019) machine-readable Table 1: B-V = 0.43±0.04, B = 13.70, E(B-V) = 0.01, D = 8.30 Mpc. Type flag: d (dwarf). Cross-ID: UGC 08837 = DDO 185 = Holmberg IV = PGC 49448. IBm (Dwarf Irregular).
- Method: **False Precision Ruling + Interpolation Principle.** Astra proposed 4.5 Gyr. Pro flagged "Redder is Younger" inversion: B-V = 0.43 (redder) would be younger than NGC 5204 (B-V = 0.41, 5.0 Gyr). Δ(B-V) = 0.02 mag is HALF the error bar (±0.04) → galaxies are statistically indistinguishable in color space. Interpolation: B-V = 0.43 bounded by NGC 5204 (0.41, 5.0) and NGC 5585 (0.46, 5.0) → snap to 5.0 Gyr. BC03 Z=0.004 constant SFH: B-V ≈ 0.44 at t = 5.0 Gyr confirms calibration.
- **Protocol reinforced: False Precision** — if Δ(color) < σ(color), galaxies are statistically indistinguishable; standardize to common bin.
- Method validated: Second de los Reyes & Kennicutt application
- Known issues: Astra's 4.5 Gyr overridden. Hidden-T3 risk flagged (not backbone anchor).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 08837 (Holmberg IV) | 5.0 Gyr | T2 (Sm/IBm Cohort) | MEDIUM | No | — |

---

### Method 140: Noordermeer & van der Hulst (2007) Table A3 B-R = 1.48±0.13 + moderate inclination (53°) + hot disk dynamics → T3 Anchor → t₅₀ [Forensic Distinction from UGC 08699]
- Data source(s): Noordermeer & van der Hulst (2007) Table A3: B(m₂₅) = 12.66±0.09, R(m₂₅) = 11.13±0.09. B-R = 1.48±0.13 (representative of m₂₅ and m_lim values: 1.53 and 1.43 respectively). Morphology: Sab, bulge-dominated. Inclination i = 53° (moderate). Noordermeer (2008): "hot disk" with high velocity dispersions. Cross-ID: UGC 09133 = NGC 5533.
- Method: **Forensic Distinction from UGC 08699.** Key discriminator: inclination. At i = 53° (sec(i) = 1.7), dust attenuation is moderate — seeing intrinsic stellar population, not "sunset effect." Blue edge: 1.48 − 0.13 = 1.35 → still firmly T3 (~7 Gyr minimum), unlike UGC 08699 (blue edge 1.05 → T1 possible). Independent dynamical confirmation: hot disk velocity dispersions require billions of years of gravitational scattering. First robust T3 anchor from Noordermeer dataset.
- Method validated: Second Noordermeer Table A3 application
- Known issues: Astra's B-V ≈ 0.87 (RC3/HyperLEDA) not directly confirmed — Noordermeer gives B-R not B-V. Biblio correction: paper is MNRAS 376, 1480 (not 1513).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 09133 (NGC 5533) | 9.3 Gyr | T3 (Anchor) | MEDIUM-HIGH | No | — |

---



## UV survey photometry

**Conversion:** UV photometry enters via GALEX colors (primarily \( \mathrm{FUV}-\mathrm{NUV}\), sometimes \( \mathrm{NUV}-r\)).  
Depending on the source and method type, UV is used either:
1) as an SPS color axis combined with optical/IR colors (BC03 mapping), or
2) as an “activity classifier” (active vs passive thresholds) combined with other evidence (structure, dust geometry) per the audit protocols.

When a UV+optical BC03 mapping is used, the conversion follows the broadband SPS grid procedure, with UV colors included in \(\chi^2\).

### Worked example

**Worked example (Method 56 — GALEX UV + optical → BC03):**  
        Compute:
        \[
        \mathrm{FUV-NUV}=\mathrm{FUV}-\mathrm{NUV},\quad (B-V)=B-V,
        \]
        and match \(\{\mathrm{FUV-NUV}, (B-V)}\) to BC03 tracks under the entry’s assumptions to obtain \(t_50\) (as recorded).

        ### Method 56: Gil de Paz+2007 GALEX Atlas UV (Table 3) + optical B-V (Table 4) → BC03 SPS mapping → t₅₀
- Data source(s): Gil de Paz et al. (2007) ApJS 173, 185. Table 3: FUV = 12.07±0.01, NUV = 11.87±0.01 → FUV-NUV = 0.20±0.01. Table 4: B = 10.46±0.11, V = 10.01±0.12 → B-V = 0.45±0.16.
- Method: Map {FUV-NUV ≈ 0.2, B-V ≈ 0.45} onto BC03 tracks → mid-sequence "Active Spiral Equilibrium" → t₅₀ ≈ 6.3 Gyr. Hyper-Active UV signal (0.20 << 0.9 threshold). Table 3 (asymptotic) prioritized over Table 4 (aperture) per Deep Think source protocol.
- Method validated: Single use
- Known issues: Astra's B-V ≈ 0.49 corrected to 0.45 (table column mixup — 0.49 was IRAS 12μm flux). No direct CMD t₅₀. "Face-On Blue Anchor" for Phase II.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4559 | 6.3 Gyr (range 5.0–7.5) | T2 (Hyper-Active Control) | HIGH | No | — |

---

### Method entries (verbatim from the audit)

### Method 13: sSFR (Boissier+2008 GALEX LSB survey) → extended SFH → t₅₀
- Data source(s): Boissier et al. (2008) Table 2 (NUV), Table 3 (SFR_NUV), Table 5 (log M★)
- SPS model: Extended/constant SFH prior
- Metallicity: F568-1: none. F568-3: strong-line ~solar (de Blok & van der Hulst 1998)
- Method validated: Yes (2 applications)
- Known issues: Both had ~10× Astra arithmetic errors corrected via "Boissier Directive" (mandatory Table 3 SFR)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F568-1 | 7.1 Gyr (range 6.2–8.6) | T2 (borderline T3) | MEDIUM-LOW | No | — |
| F568-3 | 7.6 ± 1.3 Gyr (range 6.3–8.9) | T2 (borderline T3) | MEDIUM-LOW | No | — |

---

### Method 59: Muñoz-Mateos+2009 SINGS Table 6 UV/optical asymptotic photometry → FUV-NUV + NUV-r → t₅₀ [NUV-r < 4.0 = Active]
- Data source(s): Muñoz-Mateos et al. (2009) ApJ 703, 1569. Table 6 (asymptotic AB): FUV = 12.353±0.003, NUV = 11.815±0.002, r = 8.213±0.001. Derived: FUV-NUV = 0.54, NUV-r = 3.60, g-r = 0.66.
- Method: FUV-NUV = 0.54 < 0.9 → Active. NUV-r = 3.60 < 4.0 → Blue Cloud (not Green Valley). Both metrics Active → T2. Astra classified T3/8.5 from NUV-r alone. **New protocol formalized: NUV-r < 4.0 = Active boundary (never T3).** Derivative (FUV-NUV, current activity) vs Integral (NUV-r, accumulated history) — if derivative is Active, galaxy is physically active.
- Method validated: Single use
- Known issues: Original T3/8.5 **REJECTED**. Sixth T3→T2 correction. "Late T2" at 7.6 Gyr accounts for massive evolved bulge.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5055 (M63) | 7.6 Gyr (range 7.0–8.0) | T2 (Late T2, UV Active) | MEDIUM-HIGH | No | — |

---

### Method 70: Muñoz-Mateos+2009 SINGS Table 5 UV → FUV-NUV + discrepant source resolution → t₅₀ ["Fireworks Galaxy"]
- Data source(s): Muñoz-Mateos+2009 ApJ 703, 1569. Table 5 (asymptotic AB): FUV = 10.425±0.003, NUV = 9.741±0.002 → FUV-NUV = 0.684. Astra cited Gil de Paz+2007: FUV-NUV = 0.31 (unverified, likely intrinsic/de-reddened). Δ = 0.37 mag.
- Method: Both UV values (0.31 and 0.684) are < 0.9 → both Active → T2 tier robust regardless. **Discrepancy resolved:** 0.684 = observed (foreground-corrected only); 0.31 = likely intrinsic (internally de-reddened). Internal dust (A_FUV > A_NUV) causes significant reddening. **New protocol: Discrepant Sources Same Tier** — if both values agree on tier, tier confidence is HIGH; use observed/redder value for conservative age floor. 10+ supernovae since 1917 physically confirms vigorous SF. Conservative 6.0 Gyr (up from Astra's 5.4) for massive dusty Scd spiral.
- Method validated: Single use
- Known issues: Astra's Gil de Paz optical B-V = 0.501 not independently verified. SINGS BVRI columns blank for NGC 6946. UV discrepancy explained but not definitively resolved. Conservative age adopted pending optical verification.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6946 | 6.0 Gyr (range 5.0–7.0) | T2 (Fireworks Galaxy) | MEDIUM | No | — |

---

### Method 82: Lee+2011 Table 2 GALEX UV asymptotic → FUV-NUV + Edge-On Optical Exclusion (Source Veto) → t₅₀ [UV Rescue]
- Data source(s): Lee et al. (2011) ApJS 192, 6. Table 2: FUV = 15.08±0.09, NUV = 14.74±0.07 → FUV-NUV = 0.34±0.11. MW-corrected: 0.31±0.12. Barazza & Binggeli (2002) A&A 371, 389: B-V = 1.10 with explicit author statement: "Red colour is clearly caused by internal dust absorption." **B-V = 1.10 REJECTED (Source Veto).**
- Method: FUV-NUV = 0.34 < 0.9 → Active → T2. **"Optical Exclusion" protocol:** when source author explicitly flags optical color as dust artifact, data point is disqualified (zero weight). **"UV Rescue":** GALEX UV provides clean second axis for edge-on systems. FUV-NUV = 0.34 is bluest among edge-on sample (bluer than NGC 5907's 0.43 and NGC 7331's 0.70). Dwarf ISM has patchy dust / supernova "chimneys" allowing UV escape even at i ≈ 90°. Aligned with NGC 7793 (6.5 Gyr) — both hyper-active with evolved backbones. Dust lane implies chemical evolution, distinguishing from "Transparent/Primitive" LSBs at 6.0 Gyr.
- Method validated: Single use (second Lee+2011 application; third Edge-On Dust Trap after NGC 5907, NGC 7331)
- Known issues: MEDIUM confidence — data rescue from geometric constraint. UV alone doesn't uniquely determine mass-weighted t₅₀.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 01281 | 6.5 Gyr (range 5.5–7.5) | T2 (Edge-On UV Rescue) | MEDIUM | No | — |

---

### Method 101: Ashby+2013 SFRS Table 1 GALEX UV → FUV-NUV + Hyper-Active Mass Gradient → t₅₀ [Aggregator Veto]
- Data source(s): Ashby et al. (2013) SFRS Table 1 (NGC 3274 row, SFRS ID 92): FUV = 14.687±0.006, NUV = 14.479±0.003 → FUV-NUV = 0.208±0.007. E(B-V) = 0.024, D = 10.0 Mpc. Cross-ID: UGC 05721 = NGC 3274 = PGC 31122.
- Method: FUV-NUV = 0.208 < 0.25 → **Hyper-Active** → T2. Conflicting aggregator optical colors (HyperLEDA B-V = 0.97 vs OpenNGC B-V = 0.37, spread > 0.2 mag) → **both REJECTED ("Aggregator Veto":** when aggregator optical data exhibits >0.2 mag internal conflict, disqualified). SFRS is Tier 1 peer-reviewed; single high-quality UV axis superior to averaging noisy optical. **Hyper-Active Mass Gradient:** NGC 7793 (0.17, 6.5 Gyr, massive Sd) → NGC 4559 (0.20, 6.3 Gyr) → **UGC 05721 (0.21, 6.0 Gyr, intermediate)** → NGC 5585 (0.23, 5.0 Gyr, dwarf). "Bluer is Older" inversion at massive end reflects Cosmic Downsizing. 6.0 Gyr = equilibrium point in mass hierarchy. Extended constant SFH over Hubble time → t₅₀ ≈ 12/2 = 6.0 Gyr.
- Method validated: Single use (SFRS first application)
- Known issues: Optical colors completely rejected (aggregator conflict). UV-only classification.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05721 (NGC 3274) | 6.0 Gyr (range 5.5–6.5) | T2 (Hyper-Active) | HIGH | No | — |

---

### Method 102: Wyder+2009 GALEX LSB Survey Tables 4-6 UV + SFR + SPARC mass → Twin Test with NGC 2552 → t₅₀ [LSB Twin Anchor]
- Data source(s): Wyder et al. (2009) Table 4: FUV = 17.44±0.06, NUV = 17.29±0.06. Table 5: r = 15.66±0.13. Table 6: log SFR = −0.57 (Salpeter). SPARC: L[3.6] = 3.336×10⁹ L☉, M_HI = 1.099×10⁹ M☉, D = 58.70 Mpc. M★ ≈ 1.67×10⁹ M☉. Gas fraction ~66%. Derived: FUV-NUV = 0.15, NUV-r = 1.63, SFR(Kroupa) = 0.17 M☉/yr.
- Method: b ≈ 1.4, log sSFR = −10.00, t_dbl ≈ 9.9 Gyr. **"Twin Test" with NGC 2552:** identical gas fraction (66% vs 67%), birthrate (1.4 vs 1.35), UV color (FUV-NUV = 0.15 both) → must have identical age. Pro proposed 6.0 Gyr → **rejected** (inconsistent with twin); Deep Think refined to 5.7 Gyr matching NGC 2552. LSB morphology explains inefficiency: diffuse gas, sub-critical surface density, can't collapse rapidly. FUV-NUV = 0.15 (very blue) = clean UV tracer — LSB galaxies produce almost no dust. **Twin Test principle:** when two galaxies share same gas fraction (±5%), same b (±0.1), same UV signature → must have same age.
- Method validated: Single use (Wyder+2009 first application)
- Known issues: Pro's 6.0 Gyr overridden by Twin Test with NGC 2552.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05750 | 5.7 Gyr | T2 (Vigorous LSB Dwarf) | HIGH | No | — |

---

### Method 108: Magaña-Serrano+2020 Table 2 SFR(FUV) + Pak+2014 UV + SPARC mass → Kroupa-corrected b + Fading sequence → t₅₀ [Pro Override — Toy Model Trap]
- Data source(s): Magaña-Serrano et al. (2020) Table 2: SFR(FUV) = 0.04±0.01 M☉/yr (Salpeter), SFR(Hα) = 0.015±0.004, A_FUV = 0.10 mag, μB = 23.1 mag/arcsec² (LSB). Pak et al. (2014) Table 2: FUV-NUV = 0.28, NUV-r = 2.46. SPARC: L[3.6] = 2.296×10⁹ L☉, M_HI = 0.674×10⁹ M☉, D = 18.00 Mpc. M★ ≈ 1.15×10⁹ M☉ (log 9.06). Gas fraction ~59%.
- Derived: SFR(Kroupa) = 0.025 M☉/yr, log sSFR = −10.66, b ≈ 0.30, t_dbl ≈ 46 Gyr.
- Method: **Pro initially classified T3/10.2 Gyr using exponential decay toy model — REJECTED.** Toy model baked return fraction into b, used Salpeter SFR without Kroupa correction, assumed Big Bang start for SFH. Standard b = sSFR × T gives b ≈ 0.30. **Fails ALL T3 criteria:** log sSFR = −10.66 (> −11.0), gas 59% (> 10%), b = 0.30 (> 0.1). Fading sequence interpolation: DDO 87 (b=0.40, 7.3) → UGC 06399 (b=0.30, 7.6) → UGC 03205 (b=0.23, 7.9).
- Method validated: Single use
- Known issues: Pro's T3/10.2 Gyr rejected. LSB (μB = 23.1) sub-critical density explains fading despite 59% gas.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06399 | 7.6 Gyr | T2 (Fading LSB Dwarf) | HIGH | No | — |

---



## Hα / SFR-based

**Conversion:** When a global star formation rate is available from Hα (or EW(Hα) is used as a proxy), we compute:
- stellar mass \(M_\*\) (often from SPARC 3.6μm backbone),
- SFR (IMF-harmonized when needed),
- specific SFR \( \mathrm{sSFR}=\mathrm{SFR}/M_\*\),
- and/or doubling time \(t_{\rm dbl}=M_\*/\mathrm{SFR}\).

A simple constant-SF timescale mapping is:
\[
T \approx \frac{1}{(1-R)\,\mathrm{sSFR}},\quad t_{50}\approx \frac{T}{2},
\]
where \(R\) is the mass-return fraction (documented when used).  
When the audit uses birthrate \(b\), we use:
\[
b = \frac{\mathrm{SFR}_{\rm now}}{\langle \mathrm{SFR}\rangle_{\rm past}},
\]
and (for exponential SFHs) invert \(b\) to a \(\tau\) and compute \(t_{50}\) analytically or numerically; where the audit uses protocol tie-breakers (gas fraction, EW thresholds, dust geometry), those steps are explicitly recorded per-method.

### Worked example

**Worked example (Method 79 — Hα SFR → \(t_{dbl}\) → \(t_50\)):**  
        Using the published SFR and mass backbone, compute:
        \[
        t_{dbl}=rac{M_\*}{\mathrm{SFR}}.
        \]
        This entry documents \(t_{dbl}\) and uses the pipeline’s T1/T2 discriminator rules (gas fraction + doubling time) to lock the final \(t_50\).

        ### Method 79: SPARC L[3.6] → M★ + Kaisin & Karachentsev (2011) Hα SFR → assembly timescale + gas fraction override → t₅₀ [Delayed Assembly]
- Data source(s): SPARC: L[3.6] = 0.323×10⁹ L☉, M_HI = 1.807×10⁹ M☉, D = 12.50 Mpc. Υ★ = 0.5 → M★ = 1.615×10⁸ M☉. Kaisin & Karachentsev (2011): log SFR(Hα) = −1.36 → SFR = 0.044 M☉/yr. 1981 de Vaucouleurs B-V ≈ 0.76: (colon flag = uncertain) **REJECTED** — physics conflict with Hα activity (active metal-poor dwarf should be blue).
- Method: M_HI/M★ ≈ 11 (>90% baryons in gas = "embryonic"). t_dbl = 3.67 Gyr (near T1/T2 boundary). b ≈ 3.76 (forming at ~4× lifetime average). **Protocol 19: Gas Fraction Override** — M_HI/M★ > 5 strengthens T1 case when t_dbl is near boundary. Galaxy has barely touched baryon reservoir → Delayed Assembly, not rejuvenated fossil (fossils don't retain 10× stellar mass in gas). Stochasticity ruling: FUV unnecessary, gas reservoir ensures secular "Rising" trend. 1981 optical color rejected (physics conflict: red color + active Hα in metal-poor dwarf = unreliable photometry or extinction).
- Method validated: Single use
- Known issues: Proxy-based t₅₀ (assembly timescale), not direct CMD. t_dbl = 3.67 near 4.0 Gyr threshold. Single SFR timescale (Hα ~10 Myr). Optical color rejected.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 00731 (DDO 9) | 3.7 Gyr | T1 (Delayed Assembly) | MEDIUM | No | — |

---

### Method entries (verbatim from the audit)

### Method 2: V-I integrated color (Pildis+1997) → BC03 τ-model → t₅₀
- Data source(s): Pildis, Schombert & Eder (1997) Table 2 — average V-I color
- SPS model: BC03, Chabrier IMF, τ = 1.6–2.0 Gyr
- Metallicity: Assumed from LSB/dwarf priors (Z = 0.004–0.008)
- Dust: Foreground (MW) only
- Method validated: Yes (two applications; D564-8 has independent UV cross-check)
- Known issues: Single-color degeneracy; metallicity assumed; no internal dust correction

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| D512-2 | 4.3 (+1.6/−1.4) Gyr | T2 (borderline T1/T2) | MEDIUM-LOW | No | — |
| D564-8 | 5.2 (±1.5) Gyr | T2 | MEDIUM | No | — |

---

### Method 25: Regional outer-disk color (Hunter+2013) + global sSFR → BC03 declining τ-model → t₅₀
- Data source(s): Hunter+2013 Figure 8 caption: (B-V)₀ = 0.70 ± 0.01 (outer disk). Table 1: log M★ = 10.48, log SFR(Hα) = 0.50.
- SPS model: BC03, declining τ-model
- Metallicity: Assumed near-solar
- Method validated: Single use
- Known issues: Outer-disk only; SFH-prior sensitive at T2/T3 boundary; borderline (8.3 Gyr)
- Key physics: "Outer Disk Floor" — outer disk = youngest component in inside-out formation. SFR(Hα)/SFR(V) ≈ 0.09 = fossilization.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 801 | 8.3 Gyr (range 8.0–9.5) | T3 (borderline) | MEDIUM | No | — |

---

### Method 33: BCD burst/backbone decomposition — disk color (Meurer+1994) + mass fractions (Tang+2022) → t₅₀
- Data source(s): Meurer et al. (1994) AJ 107, 2021: (B-R)₀ = 1.65 for red diffuse population. Tang et al. (2022) arXiv:2210.06485: ~10% mass in last 1 Gyr, ~3% in recent burst (~50 Myr). Werk et al. (2010): M★ = 3.2×10⁸ M☉, SFR = 0.09 M☉/yr.
- Method: Decompose into burst (~10% mass, ~0.5 Gyr) + backbone (~90% mass, ~8+ Gyr). t₅₀ = median → must penetrate backbone → snaps to backbone age.
- Method validated: Single use
- Known issues: **BCD Bi-Modal Trap** — standard b-parameter enforcement does NOT apply. b ≈ 3.9 reflects luminosity of current burst, not mass distribution. Deep Think's original 7.25 Gyr was mean age (wrong statistic); Pro corrected to median = 8.5 Gyr.
- Screen 5 note: NGC 2915 is a BCD but M★ from Werk et al. (not SPARC 0.5×L[3.6]). Screen 5 NOT triggered.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2915 | 8.5 Gyr | T3 (Ancient / Rejuvenated) | MEDIUM | No | — |

---

### Method 69: McQuinn+2010 Paper II Table 2 burst mass fraction + Lelli+2014 b-parameter (corrected for duty cycle) → t₅₀ [BCD Bi-Modal Trap #2]
- Data source(s): McQuinn et al. (2010) ApJ 724, 49 (Paper II). Table 2: total M★ = 100 ± 30 × 10⁶ M☉, burst mass fraction = 3.9 ± 1.4%, burst duration = 480 ± 70 Myr. Table 1: ⟨SFR⟩(0-6 Gyr) = 3.9×10⁻³ M☉/yr. Lelli et al. (2014): b = 3.8 ± 1.3 (uses SFR_p = peak SFR over past 1 Gyr). Drozdovsky+2001: RGB detection, "not a galaxy in formation." García-Benito & Pérez-Montero (2012): SFH episodes at ~6–10 Gyr.
- Method: **BCD Duty Cycle Trap** — Pro initially extrapolated Lelli's peak b = 3.8 over 6 Gyr → 34% recent mass → t₅₀ ≈ 7.7 Gyr (T2). **Deep Think challenge:** BCDs operate on duty cycles (100–500 Myr bursts separated by quiescence). SFR_p is peak rate, not sustained. McQuinn Table 2: burst = only 3.9% of mass. Corrected: mean SFR over 0–6 Gyr gives ~23% formed in last 6 Gyr → ~77% predates 6 Gyr → t₅₀ ≈ 9.0 Gyr (T3). Second BCD case after NGC 2915. **Protocol reinforced:** Do NOT extrapolate peak SFR over secular timescales for BCDs.
- Method validated: Single use (second BCD application after NGC 2915)
- Known issues: Pro initially called T2/7.7, corrected to T3/9.0 after Deep Think challenge. "Duty Cycle Error" = extrapolating burst peak as sustained rate. Failed T1 candidate.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6789 (UGC 11425) | 9.0 Gyr | T3 (BCD, duty cycle corrected) | MEDIUM-HIGH | No | — |

---

### Method 77: Gavilán+2013 Table 3 B-V + SFR → b parameter + t_dbl → t₅₀ [Active Builder]
- Data source(s): Gavilán et al. (2013) Table 3: B-V = 0.46±0.03, SFR = 0.1320 M☉/yr, log M★ = 9.04, 12+log(O/H) = 8.10±0.05. SPARC: L[3.6] = 2.004×10⁹ L☉, M_HI = 1.343×10⁹ M☉, D = 17.10 Mpc. Gas fraction ~129%.
- Method: b = 1.66 (forming faster than lifetime average), t_dbl = 8.3 Gyr. **T1/T2 discriminator: t_dbl.** t_dbl = 8.3 >> 3 Gyr → secular growth, not burst → T2 despite high gas fraction (129%). Gas 129% is high but not extreme (T1 needs >200% or t_dbl < 3 Gyr). b > 1 requires younger t₅₀ than constant-SFR baseline (~7 Gyr) → adjusted to 5.0 Gyr. Pro initially proposed 6.0, Deep Think corrected to 5.0.
- Method validated: Single use
- Known issues: t₅₀ inferred from activity metrics, not directly measured. Gas fraction alone insufficient for T1 — need t_dbl discriminator.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 00191 | 5.0 Gyr | T2 (Active Steady Builder) | MEDIUM | No | — |

---

### Method 94: Hunter+2009 Table 2 (B-V)₀ + Table 3 (FUV-NUV)₀ + Table 6 dual SFR (FUV + Hα) → Post-Burst Diagnostic → t₅₀ [Post-Burst Anchor]
- Data source(s): Hunter et al. (2009) Table 2: (B-V)₀ = 0.25±0.01, (U-B)₀ = −0.22±0.01, log L(Hα) = 40.055. Table 3: (FUV-NUV)₀ = 0.15. Table 6: log SFR(FUV) = −0.94 → 0.115 M☉/yr, log SFR(Hα) = −1.10 → 0.079 M☉/yr. SPARC: L[3.6] = 2.026×10⁹ L☉, M_HI = 0.678×10⁹ M☉, D = 9.60 Mpc. M★ ≈ 1.01×10⁹ M☉. Gas fraction ~67%. Cross-ID: UGC 04325 = NGC 2552.
- Method: b ≈ 1.3 (FUV-based), t_dbl ≈ 8.8 Gyr. **"Post-Burst Diagnostic" established:** SFR(FUV) > SFR(Hα) → galaxy is **declining** from burst peak (~30% drop in last 10 Myr). FUV traces ~100 Myr average; Hα traces ~10 Myr instantaneous. Extremely blue (B-V)₀ = 0.25 is "optical frosting" from bright B/A stars left by recent peak, not indicator of extreme youth. **Post-Burst Diagnostic codified:** FUV > Hα = declining, FUV < Hα = rising, FUV ≈ Hα = steady state. Why bluer than UGC 00191 (b = 1.66, B-V = 0.46)? UGC 00191 is dust-reddened by ongoing burst; NGC 2552 is seen "naked" after peak. T1 fails: gas 67% < 100%, t_dbl 8.8 >> 4 Gyr. Age between primitive (5.0) and steady-state (6.0), backbone established but burst pulls younger.
- Method validated: Single use (Hunter+2009 first application)
- Known issues: HyperLEDA B-V ≈ 0.40 superseded by Hunter+2009 (B-V)₀ = 0.25.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2552 (UGC 04325) | 5.7 Gyr | T2 (Post-Burst Dwarf) | HIGH | No | — |

---

### Method 96: Lelli+Table B.4 SFR + EW(Hα) + τ_global/τ_local + SPARC mass → Local vs Global paradox → t₅₀ [Simmering Dwarf]
- Data source(s): Lelli et al. Table B.4: SFR = 0.082 M☉/yr, EW(Hα) = 49±5 Å, M_HI = 1.15×10⁹ M☉, τ_global = 18.6 Gyr, τ_local = 6.6 Gyr. SPARC: L[3.6] = 1.552×10⁹ L☉, M_HI = 1.100×10⁹ M☉, D = 12.50 Mpc. M★ ≈ 7.76×10⁸ M☉. Gas fraction ~142%.
- Method: b ≈ 1.46, log sSFR = −9.98, t_dbl ≈ 9.5 Gyr. Gas 142% > 100% (passes T1 gas threshold) BUT t_dbl = 9.5 >> 4 Gyr (fails rapid assembly). **"Local vs Global Paradox":** τ_global = 18.6 Gyr vs τ_local = 6.6 Gyr — 3× discrepancy. Outer HI envelope is dynamically stable (sub-critical Toomre Q), decoupled from star formation. Inner disk forms at healthy pace. **"Enormous fuel tank, narrow fuel line."** EW(Hα) = 49 Å confirms genuine SF. Multi-parameter triangulation: b = 1.46 → between 5.0 and 6.0 Gyr; gas 142% (primordial) → younger end; τ_global = 18.6 → simmering → 5.4 Gyr.
- Method validated: Single use (Lelli Table B.4 first application)
- Known issues: Gas fraction above T1 threshold but assembly timescale clearly T2.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 04499 | 5.4 Gyr | T2 (Simmering Dwarf) | HIGH | No | — |

---

### Method 98: Hameed & Devereux (2005) Table 3 Hα EW + L_FIR + SPARC mass → Hα Tiebreaker at T2/T3 boundary → t₅₀ [Terminal T2 Anchor]
- Data source(s): Hameed & Devereux (2005) Table 3: F(Hα+[NII]) = 2.8×10⁻¹² erg/s/cm², L(Hα+[NII]) = 16.9×10⁴⁰ erg/s, EW(Hα+[NII]) = 14.4 Å, L_FIR = 0.69×10¹⁰ L☉, D = 22.4 Mpc. SPARC: L[3.6] = 171.582×10⁹ L☉, M_HI = 16.396×10⁹ M☉, D = 22.90 Mpc. M★ ≈ 8.58×10¹⁰ M☉ (log 10.93). Gas fraction ~19%. Cross-ID: UGC 05253 = NGC 2985.
- Method: Multi-tracer SFR ≈ 1.0±0.3 M☉/yr (FIR + Hα consistent). b ≈ 0.16, log sSFR ≈ −10.93 (borderline!). **Hα Tiebreaker at T2/T3 boundary:** EW(Hα) = 14.4 Å → direct detection of OB stars → engine still sparking. Gas = 19% → not depleted. All T3 tests fail: HÎ± detected, gas > 10%, b > 0.1. **"Terminal T2 Anchor" at 8.0 Gyr** — absolute edge of Active Sequence. Three galaxies define 8.0 Gyr limit via different mechanisms: NGC 1090 (bar inefficiency), **NGC 2985 (bulge stabilization)**, F571-8 (LSB inefficiency). All have b ≈ 0.1–0.16, gas 19–35%, Hα detected. **Sab sequence complete:** UGC 02916 (7.4, Stabilized) → UGC 03205 (7.9, Fading) → NGC 2985 (8.0, Terminal) → IC 356 (8.8, Quenched T3).
- Method validated: Single use (Hameed & Devereux 2005 first application)
- Known issues: SFR derived from FIR + Hα (not directly tabulated as SFR). log sSFR at −10.93 very close to −11.0 threshold.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2985 (UGC 05253) | 8.0 Gyr | T2 (Terminal Sab) | HIGH | No | — |

---

### Method 100: Gavilán+2013 Table 3 B-V + SFR + metallicity + SPARC mass → Sub-Critical Envelope + metallicity correction → t₅₀ [Extreme Simmering Anchor]
- Data source(s): Gavilán et al. (2013) Table 3: B-V = 0.42±0.01, SFR = 0.02347 M☉/yr, 12+log(O/H) = 8.1±0.1 (~25% solar). SPARC: L[3.6] = 0.588×10⁹ L☉, M_HI = 1.094×10⁹ M☉, D = 21.30 Mpc. M★ ≈ 2.94×10⁸ M☉ (log 8.47). Gas fraction ~372%.
- Method: b ≈ 1.1, log sSFR = −10.10, t_dbl = 12.5 Gyr. **372% gas fraction — one of highest in sample!** BUT t_dbl = 12.5 >> 4 Gyr → fails T1 by factor of 3. **"Sub-Critical Envelope":** shallow potential (low M★) + sub-critical gas density (Toomre Q > 1) + low metallicity (poor cooling, no dust, UV unshielded → no H₂) = gas trapped in diffuse stable state. More gas ≠ more SF when gas can't collapse. Contrast PGC 51017 (259% gas, deep potential → T1). **Highest gas fraction in T2.** Astra's 5.0 Gyr (color only) revised to 5.3 Gyr via metallicity correction: at 25% solar, B-V = 0.42 corresponds to ~5.3 Gyr (not 5.0). **Metallicity color correction table:** Z = 100% → 5.0, Z = 25% → 5.3, Z = 10% → 5.5 for B-V = 0.42.
- Method validated: Third Gavilán+2013 application (after UGC 00191, UGC 00891)
- Known issues: Astra's 5.0 corrected to 5.3 via metallicity.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05716 | 5.3 Gyr | T2 (Extreme Simmering Dwarf) | HIGH | No | — |

---

### Method 103: Gavilán+2013 Table 3 B-V + log M★ + SFR (Kennicutt+2008) → BC03 τ-model → t₅₀ [Fading Dwarf — Standard Verification]
- Data source(s): Gavilán et al. (2013) MNRAS 434, 2491 Table 3: B-V = 0.53 (from RC3), log M★ = 8.28, SFR = 0.0053 M☉/yr (from Kennicutt+2008). Cross-ID: UGC 05764 = DDO 83 = Arp 267. Corroboration: Grasha+2013 Table 2: SFH type C (constant), duration = 5 Gyr (model cap = ≥5 Gyr), EW(Hα) = 42.7 Å, 12+log(O/H) = 7.95.
- Method: M★ = 1.905×10⁸ M☉, sSFR = 2.78×10⁻¹¹ yr⁻¹, log sSFR = −10.56, b ≈ 0.38. **Standard "Fading Dwarf" verification — no traps.** B-V = 0.53 redder than Stability Block (0.37) → redness is stellar (evolved population), not dust. b ≈ 0.38 = forming at ~40% of past average → past peak, slowly winding down. Deep Think ruled: maintain 7.6 Gyr, do NOT round to T3 boundary — 7.6 distinguishes "Fading" from "Ancient/Fossil" (T3 typically b < 0.1, t₅₀ > 9 Gyr). No Astra baseline existed for this galaxy — Pro conducted open search.
- Method validated: Third Gavilán+2013 Table 3 application (after UGC 05716, UGC 05829)
- Known issues: MEDIUM confidence — Tier 2 data only (no UV or resolved populations).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05764 (DDO 83) | 7.6 Gyr | T2 (Fading Dwarf) | MEDIUM | No | — |

---

### Method 104: Gavilán+2013 Table 3 SFR + (B-V) + SPARC mass → Broken Twin paradox → t₅₀ [YOUNGEST T2 — Burst Peak]
- Data source(s): Gavilán et al. (2013) Table 3: SFR = 0.066 M☉/yr, (B-V) = 0.21. SPARC: L[3.6] = 0.564×10⁹ L☉, M_HI = 1.023×10⁹ M☉, D = 8.64 Mpc. M★ ≈ 2.82×10⁸ M☉ (log 8.45). Gas fraction ~363%. GHASP Table B1 Hα flux ≥ 3.3×10⁻¹⁶ W/m² (lower limit — superseded by Gavilán).
- Method: b ≈ 2.33 (highest in T2 sample!), log sSFR = −9.77, t_dbl ≈ 4.3 Gyr. (B-V) = 0.21 (extremely blue — burst peak confirmed, massive O/B production). **"Broken Twin Paradox" with UGC 05716:** nearly identical mass (log 8.45 vs 8.47) and gas fraction (363% vs 372%) but completely opposite dynamical states — UGC 05716 simmering (b = 1.1), UGC 05829 bursting (b = 2.33). **"Same fuel tank, different ignition conditions."** Gas fraction is NOT destiny; dynamical state of gas determines activity. T1 fails: t_dbl = 4.3 > 4.0 Gyr (barely). Despite extreme b = 2.33, stellar mass already exists → rejuvenation, not primordial assembly. **Youngest T2 in entire sample** at 4.5 Gyr.
- Method validated: Fourth Gavilán+2013 Table 3 application
- Known issues: GHASP lower limit bypassed by Gavilán recovery. t_dbl = 4.3 just above T1 threshold (4.0).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05829 (Spider Galaxy) | 4.5 Gyr | T2 (Youngest T2 — Burst Peak) | HIGH | No | — |

---

### Method 105: Hunter & Elmegreen (2004) Table 3 SFR + log Σ_SFR + SPARC mass → Sub-Critical Density → t₅₀ [Declining Dwarf Anchor]
- Data source(s): Hunter & Elmegreen (2004) Table 3: log L(Hα) = 38.75±0.02, SFR = 0.0034 M☉/yr, log Σ_SFR = −3.16 (extremely low!), log τ = 10.84. SPARC: L[3.6] = 0.233×10⁹ L☉, M_HI = 0.297×10⁹ M☉, D = 7.66 Mpc. M★ ≈ 1.17×10⁸ M☉ (log 8.07). Gas fraction ~254%. Cross-ID: UGC 05918 = DDO 87 = PGC 32405.
- Method: b ≈ 0.40, log sSFR = −10.54, t_dbl ≈ 34 Gyr. **Sub-Critical Density physics:** 254% gas but log Σ_SFR = −3.16 (far below Kennicutt-Schmidt threshold). Gas spread over huge diffuse halo → surface density below Toomre critical threshold → can't self-gravitate → HI can't convert to H₂ → star formation suppressed despite abundant fuel. **"Starved by lack of density, not lack of gas."** Same gas regime as PGC 51017 (259%, T1) but completely opposite outcome — PGC 51017 has super-critical density. **Gas-Rich Dwarf Sequence completed:** b = 2.33/Bursting/4.5 Gyr → 1.66/Primitive/5.0 → 1.1/Simmering/5.3 → 0.78/Steady State/6.5 → **0.40/Declining/7.3 Gyr.** T1 fails by factor of 8.5 (t_dbl = 34 >> 4 Gyr).
- Method validated: Single use (Hunter & Elmegreen 2004 first application)
- Known issues: None significant.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| DDO 87 (UGC 05918) | 7.3 Gyr | T2 (Declining Dwarf) | HIGH | No | — |

---

### Method 106: Kennicutt+2008 11HUGS Table 3 EW(Hα) + SPARC mass/gas → Edge-On Twin Test with NGC 3104 → t₅₀ [Edge-On Twin Anchor]
- Data source(s): Kennicutt et al. (2008) 11HUGS Table 3: log f(Hα) = −11.34±0.03, EW(Hα) = 64±5 Å (strong), log L(Hα) = 40.47. SPARC: L[3.6] = 4.695×10⁹ L☉, M_HI = 2.667×10⁹ M☉, D = 8.63 Mpc, i = 90° (edge-on). M★ ≈ 2.35×10⁹ M☉ (log 9.37). Gas fraction ~113%. Cross-ID: UGC 05986 = NGC 3432 = PGC 32643.
- Method: Raw b ≈ 1.0–1.4 (lower limit due to edge-on extinction). SFR(raw) = 0.16–0.23 M☉/yr. **Edge-On Extinction Effect:** at i = 90°, absolute Hα flux attenuated 2–3× by dust lane, but **EW(Hα) is robust** (line + continuum both attenuated → ratio preserved). EW = 64 Å confirms vigorous SF. True b > 1.2 after correction. **Physical Twin Test with NGC 3104:** same gas fraction (113% vs 102%), same true activity (b ≈ 1.2) → same age. "If you tilt NGC 3104 by 90°, you get NGC 3432." Pro proposed 6.2 Gyr → **rejected** (violates sequence — raw b lower limit already ≈ 1.0, true b > 1.0, must be younger than 6.0 Gyr baseline). T1 fails even with 3× extinction correction (t_dbl would be ~4–5 Gyr, still not < 4).
- Method validated: Single use (Kennicutt+2008 11HUGS first application)
- Known issues: Pro's 6.2 Gyr overridden by Twin Test. Edge-on extinction makes absolute metrics lower limits only.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3432 (UGC 05986) | 5.8 Gyr | T2 (Vigorous Edge-On) | HIGH | No | — |

---

### Method 107: Wyder+2009 GALEX LSB Survey Tables 4-6 UV + SFR + SPARC mass → Fuel-Retention Gradient interpolation → t₅₀ [Simmering LSB]
- Data source(s): Wyder et al. (2009) Table 4: FUV = 17.03±0.08, NUV = 16.86±0.06. Table 5: r = 15.15±0.06. Table 6: log SFR = −0.54 (Salpeter). SPARC: L[3.6] = 3.384×10⁹ L☉, M_HI = 2.022×10⁹ M☉, D = 47.70 Mpc. M★ ≈ 1.69×10⁹ M☉ (log 9.23). Gas fraction ~120%.
- Derived: FUV-NUV = 0.17, NUV-r = 1.71, SFR(Kroupa) = 0.18 M☉/yr, b ≈ 1.5, log sSFR = −9.97, t_dbl ≈ 9.4 Gyr.
- Method: Same Wyder+2009 data pipeline as Method 102 (UGC 05750) but calibrated via **Fuel-Retention Gradient** instead of Twin Test. At constant b (~1.4–1.5), gas fraction governs age: more retained gas → suppressed early history → younger backbone. Interpolation: UGC 04499 (142%, 5.4 Gyr) → UGC 05999 (120%, 5.5 Gyr) → UGC 05750 (66%, 5.7 Gyr). Pro proposed 5.8 Gyr → overridden by Deep Think gradient calibration.
- Method validated: Second Wyder+2009 application (after UGC 05750); distinct inference chain (Fuel-Retention Gradient vs Twin Test)
- Known issues: Pro's 5.8 Gyr overridden. LSB sub-critical density prevents T1 despite 120% gas.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 05999 | 5.5 Gyr | T2 (Simmering LSB Dwarf) | HIGH | No | — |

---

### Method 109: James+2004 HaGS Table 3 SFR + Pak+2014 photometry + SPARC mass → IMF-corrected b → t₅₀ [IMF Trap + Distance Paradox]
- Data source(s): James et al. (2004) HaGS Table 3: SFR(Salpeter) = 0.1095±0.047 M☉/yr (at D = 10.5 Mpc), EW = 8.0 nm. Pak et al. (2014): g-r = 0.14, FUV-NUV = 0.21, (B-V) ≈ 0.36 (transformed). SPARC: L[3.6] = 0.988×10⁹ L☉, M_HI = 1.379×10⁹ M☉, D = 12.00 Mpc. M★ ≈ 4.94×10⁸ M☉ (log 8.69). Gas fraction ~279%.
- Derived: SFR scaled to D=12 Mpc, Kroupa-corrected → SFR ≈ 0.09 M☉/yr, b ≈ 2.5, log sSFR = −9.74, t_dbl ≈ 5.5 Gyr.
- Method: **Pro initially flagged T1 (t_dbl = 3.5 Gyr) — CORRECTED.** Two errors: (1) mixed D=10.5 and D=12 Mpc (t_dbl is distance-invariant since M★ ∝ D² and SFR ∝ D²), (2) divided Kroupa mass by Salpeter SFR without harmonizing IMF. Corrected t_dbl = 5.5 Gyr → safely T2. Interpolation: UGC 05829 (b=2.33, 4.5) → UGC 06446 (b=2.5, 4.6) → UGC 00191 (b=1.66, 5.0).
- Method validated: Single use (second HaGS application)
- Known issues: Pro's T1 alarm corrected. IMF harmonization critical.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06446 | 4.6 Gyr | T2 (Vigorous Starburst Dwarf) | HIGH | No | — |

---

### Method 111: James+2004 HaGS Table 3 SFR + SPARC mass → IMF-corrected b + Iso-Age Diagonal → t₅₀ [Face-On Precision]
- Data source(s): James et al. (2004) HaGS Table 3: SFR(Salpeter) = 0.2811±0.0837 M☉/yr, EW(Hα+[NII]) = 2.4 nm, D = 15.2 Mpc. SPARC: L[3.6] = 3.739×10⁹ L☉, M_HI = 1.500×10⁹ M☉, D = 15.10 Mpc, i = 20° (nearly face-on). M★ ≈ 1.87×10⁹ M☉ (log 9.27). Gas fraction ~80%.
- Derived: SFR(Kroupa) = 0.176 M☉/yr, b ≈ 1.3, log sSFR = −10.03, t_dbl ≈ 10.6 Gyr.
- Method: Face-on (i=20°) provides pristine Hα — no extinction correction needed, b is accurate (not lower limit). **Iso-Age Diagonal:** opposing vectors cancel — higher b (1.3 vs 1.21) pulls younger, lower gas (80% vs 102%) pulls older → same age as NGC 3104 (5.8 Gyr). Pro proposed 5.6 Gyr (return fraction method) → overridden.
- Method validated: Second HaGS application (after UGC 06446)
- Known issues: Pro's 5.6 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06628 | 5.8 Gyr | T2 (Vigorous Face-On Dwarf) | HIGH | No | — |

---

### Method 119: James+2004 HaGS Table 3 SFR + Sanders (1998) B-V + SPARC mass → distance + IMF corrected b → Quenching Hierarchy → t₅₀ [Terminal Circumnuclear Burn]
- Data source(s): James et al. (2004) HaGS Table 3: SFR(Salpeter) = 0.334±0.156 M☉/yr at D = 9.7 Mpc, EW = 2.6 nm. Sanders (1998): (B-V) = 0.64 (dust-reddened). SPARC: L[3.6] = 53.870×10⁹ L☉, M_HI = 1.753×10⁹ M☉, D = 18.00 Mpc, i = 71°. M★ ≈ 2.69×10¹⁰ M☉ (log 10.43). Gas fraction ~6.5%. Cross-ID: IC 750.
- Derived: SFR scaled D=9.7→18 Mpc (×3.44) then Kroupa (÷1.6) = 0.72 M☉/yr. b ≈ 0.37, log sSFR = −10.57, t_dbl ≈ 37 Gyr.
- Method: **Quenching Hierarchy: Activity > Gas Fraction.** Gas at IC 356-level (6.5%) but b = 0.37 (T2-level). Fails BOTH T3 criteria: log sSFR = −10.57 (> −11.0), b = 0.37 (> 0.1). Pro proposed T3/9.4 Gyr → overridden. **Terminal Circumnuclear Burn:** last gas funneled to nucleus via bar/interaction, outer disk dead, core still burning. (B-V) = 0.64 is dust reddening, not old population. Oldest active spiral: NGC 2985 (8.0) → UGC 06973 (8.2) → NGC 3900 (8.6, T3).
- Method validated: Fourth HaGS application
- Known issues: Pro's T3/9.4 Gyr overridden. Dusty inner region → HÎ±-derived SFR is strict lower limit.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06973 (IC 750) | 8.2 Gyr | T2 (Terminal Circumnuclear Burn) | HIGH | No | — |

---

### Method 126: Botticella+Table 2 UV SFR + HaGS HÎ± + GHASP B-V + SPARC mass → Stochastic Flicker (UV > Hα) + Fuel Retention → t₅₀ [Steady-State Disk]
- Data source(s): Botticella et al. Table 2: SFR(UV) = 0.17 M☉/yr, SFR(Hα) = 0.10, M★ = 2×10⁹ M☉ (matches SPARC). James+2004 HaGS Table 3: SFR(Hα) = 0.1088 M☉/yr at D=6.8 Mpc. GHASP: (B-V) = 0.46. SPARC: L[3.6] = 4.109×10⁹ L☉, M_HI = 0.722×10⁹ M☉, D = 8.00 Mpc, i = 47°. M★ ≈ 2.05×10⁹ M☉ (log 9.31). Gas fraction ~35%. Cross-ID: NGC 4242, SABdm.
- Derived: UV-based b ≈ 1.14, log sSFR = −10.08, t_dbl ≈ 12 Gyr.
- Method: **Stochastic Flicker Principle:** UV (100 Myr timescale) > Hα (10 Myr) in patchy Sdm disks — galaxy in temporary Hα trough. UV is correct tracer for t₅₀ calibration. **Pro's Return-Fraction Trap:** Pro used b = SFR×(1-R)×T/M★ → b=0.80 (artificially low). Correct: b = sSFR × T = 1.14. Return fraction already baked into M★. Pro's 7.6 Gyr overridden. Fuel Retention at constant b~1.0: UGC 07089 (68%, 5.9) → NGC 1003 (50%, 6.0) → UGC 07323 (35%, 6.4). (B-V) = 0.46 confirms steady-state mix.
- Method validated: Single use (first Botticella application)
- Known issues: Pro's 7.6 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07323 (NGC 4242) | 6.4 Gyr | T2 (Steady-State Disk) | HIGH | No | — |

---

### Method 128: Nandi+2023 Table 7 multi-tracer SFR + GALEX Atlas FUV/NUV + O'Neill+2006 nuclear UV + SPARC mass → Stochastic Flicker + Structural Veto → t₅₀ [Extreme Boundary — AGN Hypothesis Retracted]
- Data source(s): Nandi et al. (2023) Table 7: SFR(FUV) = 0.47, SFR(NUV) = 0.55, SFR(Hα) = 0.32, SFR(1.4GHz) = 0.52, SFR(SED) = 0.21 M☉/yr (all Salpeter). GALEX Atlas: FUV = 11.67±0.01, NUV = 11.54±0.01, FUV-NUV = 0.13, (B-V) = 0.46. O'Neill et al. (2006): nuclear UV AB ≈ 19.4 → 7.73 mag gap → AGN = 0.08% of galaxy FUV. SPARC: L[3.6] = 2.436×10⁹ L☉, M_HI = 1.779×10⁹ M☉, D = 4.74 Mpc, i = 46°. M★ ≈ 1.22×10⁹ M☉ (log 9.09). Gas fraction ~146%. Cross-ID: NGC 4395, Seyfert 1 nucleus.
- Derived: UV/Radio Kroupa average ≈ 0.32 M☉/yr, b ≈ 3.6, log sSFR = −9.58, t_dbl ≈ 3.8 Gyr.
- Method: **Adversarial Review:** Deep Think claimed AGN contamination inflated UV/Radio — **Pro rebutted** (nuclear UV = 0.08% of galaxy FUV, C-config radio avoids AGN-only). Deep Think conceded. UV/Hα divergence = Stochastic Flicker (100 Myr vs 10 Myr). **Structural Veto prevents T1:** despite t_dbl = 3.8 Gyr (technically T1), mature ~10⁹ M☉ disk with SMBH and spiral arms cannot be "embryonic." Late-Stage Rejuvenation, not primordial collapse. Highest confirmed b in T2 sequence. Bivariate: UGC 06446 (279%, 4.6) → NGC 4395 (146%, 4.7) → UGC 07261 (158%, 4.8).
- Method validated: Single use
- Known issues: Deep Think AGN claim retracted. Borderline t_dbl — Structural Veto is sole T2 discriminant.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07524 (NGC 4395) | 4.7 Gyr | T2 (Extreme Boundary Starburst) | HIGH | No | — |

---

### Method 131: Kennicutt+2008 11HUGS Table 3 EW(Hα) → extinction-proof b + SPARC mass/gas → Fuel Retention Gradient → t₅₀ [Multi-Wavelength Extinction Trap]
- Data source(s): Kennicutt et al. (2008) 11HUGS Table 3: EW(Hα+[NII]) = 35±4 Å, [NII]/Hα = 0.14, log L(Hα) = 39.75. SPARC: L[3.6] = 0.376×10⁹ L☉, M_HI = 0.258×10⁹ M☉, D = 4.70 Mpc, i = 78°. M★ ≈ 1.88×10⁸ M☉ (log 8.27). Gas fraction ~137%. Cross-ID: NGC 4455.
- Derived: Raw Kroupa SFR/M★ gives false b ≈ 0.73 (extinction-trapped). **EW(Hα) = 35 Å is extinction-proof** (line and continuum attenuated equally → ratio cancels). EW mapping: <20 Å = fading (b<0.5), 30–50 Å = steady-state (b~1.0), >100 Å = starburst (b>2.5). True b ≈ 1.0, t_dbl ≈ 13.8 Gyr.
- Method: **Multi-Wavelength Extinction Trap:** at i=78°, Hα is heavily attenuated but M★ (3.6μm) is unaffected → sSFR = SFR/M★ systematically suppressed → false "declining" signal. EW solution bypasses this. Fuel Retention: UGC 04499 (142%, b=1.46, 5.4) → NGC 4455 (137%, b~1.0, 5.7) → UGC 07089 (68%, b=1.06, 5.9).
- Method validated: First 11HUGS application
- Known issues: NOT NGC 7603 (Seyfert at UGC 12493).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07603 (NGC 4455) | 5.7 Gyr | T2 (Obscured Steady-State) | HIGH | No | — |

---

### Method 144: Lelli+2013 Table B.4 Hα SFR (Salpeter→Kroupa ×0.68) + SPARC M★ → b-parameter + τ-model → T2/T3 Boundary Anchor → t₅₀
- Data source(s): Lelli et al. (2013) Table B.4: SFR = 0.486 M☉/yr (Salpeter), EW(Hα) = 35±3 Å, M_HI = 2.5×10⁹ M☉, τ_global = 6.8 Gyr, τ_local = 4.0 Gyr. SPARC: log M★ = 9.78 (~6.0×10⁹ M☉), D = 23.7±4.4 Mpc (TF), i ≈ 30° (face-on). Cross-ID: UGC 11557 = PGC 64616 = MCG +10-29-005 = CGCG 304-005.
- Derived: SFR(Kroupa) = 0.330 M☉/yr (×0.68), t_dbl = 18.2 Gyr, b = 0.75, gas fraction ~43%.
- Method: **T2/T3 Boundary Anchor.** t₅₀ = 7.8 Gyr is within ~0.2 Gyr of 8.0 boundary, but physical state indicators all favor T2: b = 0.75 (gently declining, not quenched), gas fraction 43% (healthy, not exhausted <10%), EW = 35 Å (active), τ_local = 4.0 Gyr (efficient, sustainable). τ_global/τ_local = 1.7× is normal inside-out disk (NOT UGC 04499-style bottleneck). Face-on geometry (i ≈ 30°) ensures clean Hα. IMF: Lelli explicitly states Salpeter; ×0.68 = Kroupa (distinct from FUV ÷1.59 = Chabrier).
- **Boundary Decision Framework established:** When t₅₀ within ~1 Gyr of tier boundary, use physical state indicators (b, gas fraction, EW, τ_local) as tie-breakers.
- Method validated: Second Lelli+2013 Table B.4 application
- Known issues: TF distance ~19% uncertainty. Single SFR tracer (Hα only).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 11557 | 7.8 Gyr | T2 (Mature Steady Builder) | MEDIUM | No | — |

---

### Method 145: van Zee & Haynes (2006) ApJ 636, Table 6 (B-V)₀ = 0.37±0.02 → Twin Anchor (UGC 09992) → Stability Block → t₅₀ [Blue Inversion Physics]
- Data source(s): van Zee & Haynes (2006) Table 6: (B-V)₀ = 0.37±0.02, M_HI/L_B = 2.1 (high gas fraction), μ₀,B = 23.7 mag/arcsec² (LSB), 12+log(O/H) = 8.00±0.20 (~0.2–0.3 Z☉), D = 17.9 Mpc. Cross-ID: UGC 11820.
- Method: **Twin Anchor with UGC 09992** — identical B-V = 0.37, both LSB, high gas fraction, subsolar metallicity, no Second Axis signals. Straight Stability Block. **Blue Inversion Physics confirmed:** LSB dwarfs appear BLUER than Active T2 cohort (B-V = 0.41) despite being OLDER (6.0 vs 5.0 Gyr) because they are dust-poor and metal-poor — we see "naked" stellar population without reddening.
- Method validated: First van Zee & Haynes (2006) application
- Known issues: Colors corrected for Galactic extinction only, NOT internal (acceptable for LSB dwarf with low metallicity).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 11820 | 6.0 Gyr | T2 (Stability Block) | MEDIUM-HIGH | No | — |

---

### Method 147: Lelli+2013 Table B.4 Hα SFR (Salpeter→Kroupa ×0.68, distance-harmonized) + SPARC M★/gas → b ≈ 0.97 → Constant-SFH benchmark → t₅₀ [Steady-State Archetype]
- Data source(s): Lelli et al. (2013) Table B.4: SFR = 0.060 M☉/yr (Salpeter), EW(Hα) = 40±6 Å, D = 9.2±1.7 Mpc (TF). SPARC: L[3.6] = 1.301×10⁹ L☉, M_HI = 1.744×10⁹ M☉, D = 9.77 Mpc, i = 46°. Cross-ID: UGC 12632 = PGC 71596. Sm dwarf.
- Derived: M★ ≈ 6.5×10⁸ M☉ (log 8.81). SFR scaled to SPARC distance: 0.060 × (9.77/9.2)² = 0.068 → Kroupa ×0.68 = 0.046 M☉/yr. t_dbl = 14.1 Gyr, b = 0.97, gas fraction ~268%.
- Method: **Steady-State Archetype.** b ≈ 0.97 (essentially unity) = constant SFH over Hubble time. For constant SFR: t₅₀ = half Hubble time ≈ 6.85 Gyr; with b = 0.97 (mild decline), marginally older → 7.0 Gyr. Gas fraction 268% validates flat SFH is sustainable (not transient coincidence). "Enormous fuel tank, narrow fuel line" (same class as UGC 04499). Distance harmonization: SFR scaled from Lelli D=9.2 to SPARC D=9.77 via D² correction. T1 vetoed (t_dbl = 14.1 >> 4 Gyr despite 268% gas). T3 physically vetoed (b ≈ 1, gas 268%, EW 40 Å).
- Method validated: Third Lelli+2013 Table B.4 application
- Known issues: Initial premature Data Gap — sent Pro back to check SPARC/Lelli tables. TF distance uncertainty ~15–18%.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 12632 | 7.0 Gyr | T2 (Steady Builder) | HIGH | No | — |

---

### Method 148: Lelli+2014 Table B.4 Hα SFR (Salpeter→Kroupa ×0.68) + SPARC M★/gas → b ≈ 0.96 → Constant-SFH → t₅₀ [Metal-Poor EW Calibration Anchor — Galaxy #100]
- Data source(s): Lelli et al. (2014) Table B.4 / James et al. (2004): SFR = 0.086 M☉/yr (Salpeter), EW(Hα) = 88±9 Å, τ_global = 49.9 Gyr, τ_local = 10.6 Gyr. SPARC: L[3.6] = 1.667×10⁹ L☉, M_HI = 3.660×10⁹ M☉, D = 13.20 Mpc, i = 39° (face-on), Q = 1. Cross-ID: UGC 12732 = PGC 72087 = MCG+04-55-042 = CGCG 476-106. Sm dwarf.
- Derived: M★ ≈ 8.34×10⁸ M☉ (log 8.92), SFR(Kroupa) = 0.058 M☉/yr, t_dbl = 14.3 Gyr, b = 0.96, gas fraction ~439% (most extreme in sample).
- Method: **Steady-State Archetype.** b ≈ 0.96 = constant SFH → t₅₀ ≈ half Hubble time = 7.0 Gyr. Gas 439% validates flat SFH as sustainable. **Metal-Poor EW Resolution:** EW(Hα) = 88 Å appears high but is NOT a starburst indicator. At low metallicity, old continuum is intrinsically faint (metal-poor RGB burns hotter/bluer) and massive stars produce more ionizing UV → Starburst99 predicts EW ~80–100 Å for constant SFH at this Z. Galaxy is optically thin (i = 39°, low Z = no dust). T1 vetoed (t_dbl = 14.3 >> 4 Gyr). T3 physically vetoed (b ≈ 1, gas 439%, EW 88 Å).
- **Pipeline rule established:** Do not interpret high EW as burst evidence without checking metallicity context.
- Method validated: Fourth Lelli+2013/2014 Table B.4 application
- Known issues: Single SFR tracer (Hα only). SFR uncorrected for internal extinction (negligible at this i and Z).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 12732 | 7.0 Gyr | T2 (Steady Builder) | HIGH | No | — |

---



## sSFR / mass-based

**Conversion:** For methods where the primary data is \((M_\*, \mathrm{SFR})\) or sSFR directly, compute:
\[
\mathrm{sSFR}=\frac{\mathrm{SFR}}{M_\*},\quad t_{\rm dbl}=\frac{M_\*}{\mathrm{SFR}},\quad b \approx \mathrm{sSFR}\,T_{\rm U},
\]
with \(T_{\rm U}\approx 13.8\,\mathrm{Gyr}\) used as the lifetime scale in the audit’s birthrate logic.  
A constant-SF mapping gives \(t_{50}\approx T_{\rm U}/2\) at \(b\approx 1\), with younger/older deviations calibrated by the audit’s monotonic b→age sequence and tie-breaker rules (gas fraction, morphology/downsizing penalties, quenching diagnostics).

### Worked example

**Worked example (Method 57 — direct \(M_\*\) and SFR → sSFR threshold):**  
        Compute:
        \[
        \mathrm{sSFR}=rac{\mathrm{SFR}}{M_\*},\quad \log(\mathrm{sSFR}).
        \]
        This entry uses the audit’s quiescent/active thresholds (e.g., \(\log(\mathrm{sSFR})<-11\Rightarrow\) quiescent) plus corroborating evidence to assign \(t_50\).

        ### Method 57: Richards+2015 Table 1 (M★, SFR, gas) → sSFR → b → passive threshold → t₅₀ [Anemic Spiral]
- Data source(s): Richards et al. (2015) MNRAS, "Baryonic distributions in the dark matter halo of NGC 5005." Table 1: M★ = 9.16×10¹⁰ M☉, SFR = 0.67 M☉/yr, (B-R)₀ = 1.43, M_HI = 1.22×10⁹ M☉.
- Method: log sSFR = −11.14 → b ≈ 0.10 → crosses passive threshold (b < 0.1, log sSFR < −11.0). Gas fraction 1.3% confirms reservoir exhaustion. "Anemic Spiral" — spiral arms persist as density waves in old stellar disk but no longer form stars actively.
- Method validated: Single use
- Known issues: Morphology paradox (SABbc classified T3) resolved by van den Bergh "Anemic Spiral" classification. Completes T1→T2→T3 spiral evolutionary sequence.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5005 | 9.3 Gyr (range 8.5–10.0) | T3 (Anemic Spiral) | HIGH | No | — |

---

### Method entries (verbatim from the audit)

### Method 6: sSFR (Karachentsev+2017) → BC03 mass-return bookkeeping → t₅₀ [Spectroscopic metallicity + resolved stellar constraint]
- Data source(s): Karachentsev et al. (2017) Table 2: FUV = 14.6, log M★ = 8.75, log sSFR = −9.80
- SPS model: BC03, Chabrier IMF, constant SFH (primary)
- Metallicity: Direct spectroscopic — Lee et al. (2003): 12+log(O/H) = 8.06; confirmed by Alabi et al. (2025)
- Mass-return fraction: R ≈ 0.48
- Method validated: Single use
- Known issues: sSFR not a unique clock; dust attenuation on FUV unknown; no optical color cross-check
- Strengths: SFH degeneracy broken by Alabi+2025 finding old stars (>10 Gyr) in DDO 161

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| DDO 161 | 5.0 (±1.5) Gyr | T2 | MEDIUM-LOW | No | — |

---

### Method 26: 3.6μm mass decomposition (Fraternali+2011) + tabulated SFR (Yoon+2021) → sSFR → Downsizing-adjusted anchor → t₅₀
- Data source(s): Fraternali et al. (2011) Tables 4–5 (3.6μm disk+bulge) + Yoon et al. (2021) Table 1: SFR = 5.0 M☉/yr
- Derived: M★ ≈ 1.1–1.2 × 10¹¹ M☉, log sSFR ≈ −10.35, b ≈ 0.63
- Method validated: Single use
- Known issues: Edge-on SFR systematic; +0.6 Gyr Downsizing penalty is a calibration choice
- Correction: Pro called T3/8.1; Deep Think overruled — b ≈ 0.63 cannot be "Passive"
- "Oldest Active" anchor

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 0891 | 7.8 Gyr | T2 (Oldest Active) | MEDIUM | No | — |

---

### Method 31: HERACLES sSFR (Leroy+2013) + backbone colors (B-V, B-K) → sSFR-to-age mapping → t₅₀
- Data source(s): Leroy et al. (2013) HERACLES: log M★ = 10.7, r₂₅ = 14.2 kpc, ⟨ΣSFR⟩ = 1.4×10⁻³ M☉/yr/kpc². Gil de Paz+2007 / 2MASS: B-V ≈ 0.87, B-K ≈ 4.03.
- Derived: SFR ≈ 0.50 M☉/yr, sSFR ≈ 1.0×10⁻¹¹ yr⁻¹, t_dbl ≈ 100 Gyr
- Method validated: Single use
- Known issues: No resolved CMD (14 Mpc — too distant); sSFR aperture-dependent (0.75 r₂₅ convention); SFH shape assumed declining. Original UV-only (Astra) upgraded to full proxy stack.
- "Standard Massive Spiral" archetype — 10× more quiescent than NGC 801 (T3/8.3).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2841 | 9.0 Gyr (range 8–11) | T3 | MEDIUM | No | — |

---

### Method 37: Pezzulli+2015 νM (specific mass growth rate) + GALEX UV + optical B-V → sSFR limit check → t₅₀
- Data source(s): Pezzulli et al. (2015) MNRAS 451, 2324 Table 2: M★ = (32.0 ± 1.8)×10⁹ M☉, νM = (8.05 ± 0.81)×10⁻² Gyr⁻¹. Gil de Paz+2007 Table 3: FUV-NUV = 0.19 ± 0.02. Table 4: B-V ≈ 0.54.
- Method: νM > 1/T_universe (0.072 Gyr⁻¹) → SFH is constant or rising → declining/front-loaded history physically ruled out → max t₅₀ ≈ 7.0 Gyr.
- Method validated: Single use
- Known issues: Original spreadsheet claimed T3/9.5 Gyr — **REJECTED** as undocumented and physically inconsistent. No resolved CMD (14 Mpc). Could be Hidden T3 if future CMD data shows ancient backbone, but current data does not support.
- Precedent: NGC 2403 (blue Sc) is T3 with Platinum CMD data; NGC 3198 (blue Sc) defaults to T2 without it.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3198 | 7.0 Gyr (range 5.5–7.5) | T2 (Main Sequence) | MEDIUM | No | — |

---

### Method 61: SPARC/WISE M★ + Misiriotis+2004 SFR → sSFR → Mass Penalty (Downsizing) → t₅₀ [Super-Spiral]
- Data source(s): SPARC/WISE II Table 1: log M★ = 11.607. Misiriotis et al. (2004) Table 1: log SFR = 1.36 (model-derived from IR/submm). Derived: log sSFR = −10.25, b ≈ 0.8.
- Method: Raw sSFR mapping gives t₅₀ ≈ 6.8 Gyr, but this underestimates age for super-spirals (log M★ > 11). **Mass Penalty (Protocol 20)** applied: massive ancient bulge formed at z > 2 dominates mass budget; current disk activity is small fraction. +0.7 Gyr correction → 7.5 Gyr. Cirrus contamination (old-star-heated dust inflating IR SFR) also noted. Even worst-case SFR (halved): b ≈ 0.4 → still T2.
- Method validated: Single use
- Known issues: SFR model-derived (±0.30 dex). Cirrus contamination possible. Mass penalty physics-motivated but not empirically calibrated. Most massive star-forming spiral in sample.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5371 | 7.5 Gyr | T2 (Super-Spiral) | MEDIUM | No | — |

---

### Method 64: Di Teodoro & Fraternali (2014) Table 2 (M_bary, M_HI, SFR) → sSFR → Seyfert correction → Downsizing → t₅₀ [Suppressed Super-Spiral]
- Data source(s): Di Teodoro & Fraternali (2014) WHISP sample Table 2: M_HI = 10.76×10⁹ M☉, M_bary = 144.58×10⁹ M☉, SFR = 1.19 M☉/yr (IRAS). Derived: M★ = M_bary − 1.4×M_HI = 1.295×10¹¹ M☉ (log M★ ≈ 11.1). log sSFR = −11.04 (upper limit due to Seyfert AGN contamination). SPARC/WISE cross-check: log M★ ≈ 11.286.
- Method: b ≈ 0.1 at passive threshold. **Seyfert 1 correction:** IRAS SFR contaminated by AGN-heated dust → true stellar SFR strictly lower → true log sSFR likely < −11.3, deeper in passive regime. **Second T3 pathway established:** "Suppression" (fuel exists but engine broken, gas fraction ~8%) vs NGC 5005's "Exhaustion" (fuel depleted, gas fraction ~1%). Mechanisms: AGN feedback, virial shock heating, morphological quenching. Downsizing penalty (+0.7 Gyr over NGC 5005 base) for higher mass.
- Method validated: Single use
- Known issues: AGN contamination not precisely quantified. Suppression mechanism not determined. Downsizing penalty approximate.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 5985 (UGC 09969) | 10.0 Gyr | T3 (Suppressed Super-Spiral) | HIGH | No | — |

---

### Method 67: SPARC M★ + Strickland+2003 SFR_IR → distance-corrected sSFR → t₅₀ [Void isolation]
- Data source(s): SPARC (Lelli+2016): L[3.6] = 12.845×10⁹ L☉, M_HI = 1.744×10⁹ M☉, D = 6.26 Mpc (Cepheid). M★ ≈ 6.42×10⁹ M☉ (log M★ ≈ 9.81). Gas fraction 27%. Strickland et al. (astro-ph/0306598) Table 2: SFR_IR = 0.20 M☉/yr at D = 5.2 Mpc.
- Method: **Distance consistency correction** required — Strickland used D = 5.2 Mpc, SPARC uses D = 6.26 Mpc (Cepheid). SFR scales as D²: corrected SFR = 0.29 M☉/yr. Corrected: log sSFR = −10.35, b ≈ 0.62 (up from uncorrected b ≈ 0.43). **"Retarded evolution"** in Local Void: no environmental harassment → high gas retention (27%), sustained activity, delayed aging vs cluster galaxies (NGC 3917: 15% gas, b ≈ 0.4, 7.0 Gyr).
- Method validated: Single use
- Known issues: Distance correction changes b from 0.43→0.62 and age from 7.0→6.5 Gyr. Cepheid distance is gold standard. "The Lonely Galaxy" — void isolation as evolutionary modifier.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6503 | 6.5 Gyr | T2 (Field Spiral / Lonely Galaxy) | MEDIUM-HIGH | No | — |

---

### Method 74: Wang+2016 Table 1 WISE 22μm SFR + 2MASS K-band M★ → sSFR → t₅₀ [Fossil Spiral]
- Data source(s): Wang et al. (2016) MNRAS 457, 1385. Table 1: M★ = 7.02×10¹⁰ M☉ (2MASS K-band + color-dependent M/L), SFR = 0.14 M☉/yr (WISE 22μm W4, dust-robust). Derived: log sSFR = −11.70, b ≈ 0.027. Gas fraction ~1.6%.
- Method: b ≈ 0.03 → deeply quenched, forming at <3% of lifetime average. t_dbl ≈ 500 Gyr → galaxy is effectively dead. **"Fossil Spiral"** — endpoint of disk galaxy evolution via exhaustion pathway. Represents what Anemic spirals (NGC 5005) become after final gas consumption. WISE 22μm SFR is dust-robust for edge-on system (i ≈ 90°). **Quenching Triad completed:** Anemic (NGC 5005, b~0.1, dying) → Suppressed (NGC 5985, b~0.1, choked) → Fossil (NGC 7814, b~0.03, dead).
- Method validated: Single use
- Known issues: Edge-on (i ≈ 90°) carries some systematic uncertainty despite dust-robust methods. t₅₀ inferred from sSFR→SFH mapping, not directly measured. Tier separation unambiguous (0.7 dex below threshold).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 7814 | 10.5 Gyr | T3 (Fossil Spiral) | HIGH | No | — |

---

### Method 86: Carvalho+2024 direct M★ + SFR → sSFR diagnostic → t₅₀ [sSFR Protocol Established]
- Data source(s): Carvalho et al. (2024) A&A (arXiv:2410.16467): M★ = (4.83±1.52)×10¹¹ M☉, SFR = 1.63±0.72 M☉/yr. Both per-object, peer-reviewed.
- Method: sSFR = 3.37×10⁻¹² yr⁻¹, log sSFR = −11.47. Mass-doubling time τ ≈ 300 Gyr (25× age of universe). Mass growth ~0.34%/Gyr. **sSFR Diagnostic formally established:** log sSFR < −11 → Quiescent (T3). **sSFR thresholds:** T1 Starburst (> −9.5, τ < 3 Gyr), T2 Active (−9.5 to −11, 3 < τ < 100 Gyr), T3 Quiescent (< −11, τ > 100 Gyr). **"Massive Galaxy Insulation":** for giants (M★ > 10¹¹), sSFR supersedes UV because UV measures "frosting" temperature while sSFR measures growth efficiency relative to total mass. "Rubin's Galaxy" — one of largest known spirals (~5× Milky Way mass). PLATINUM confidence.
- Method validated: Single use
- Known issues: t₅₀ = 9.0 Gyr is inferred from quiescent status, not directly measured. Aligns with NGC 2841 (9.0 Gyr) as massive quiescent spiral.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 02885 ("Rubin's Galaxy") | 9.0 Gyr (range 8.0–11.0) | T3 (Quiescent Massive Spiral) | PLATINUM | No | — |

---

### Method 87: Di Teodoro & Fraternali (2014) Table A.1 IRAS SFR + SPARC mass → b + morphological stabilization → t₅₀ [Twin Convergence]
- Data source(s): Di Teodoro & Fraternali (2014) Table A.1: SFR = 2.45 M☉/yr (IRAS), M_HI = 23.12×10⁹ M☉, M_bar = 94.12×10⁹ M☉, M★(derived) = 6.18×10¹⁰ M☉, D = 68.02 Mpc. SPARC: L[3.6] = 124.153×10⁹ L☉, M_HI = 23.273×10⁹ M☉, D = 65.40 Mpc. M★ ≈ 6.21×10¹⁰ M☉ (log 10.79). Gas fraction ~37%.
- Method: b ≈ 0.55, log sSFR = −10.40, t_dbl ≈ 25 Gyr. **"Efficiency Paradox":** 37% gas should yield b ~ 0.8–1.0, but only 0.55. **Morphological Stabilization:** Sab massive bulge → deep potential well + rotational shear → high Toomre Q → suppressed fragmentation. Galaxy retains gas not because young, but because bulge prevents rapid consumption. **"Twin Convergence"** with NGC 5033: both land at 7.4 Gyr via opposite mechanisms — NGC 5033 (Sc, 14% gas, starved by environment) vs UGC 02916 (Sab, 37% gas, self-regulating via bulge). Same age, different paths. **Pre-S0 phase:** evolutionary stage just before complete quenching (NGC 1167).
- Method validated: Single use (second Di Teodoro & Fraternali 2014 application after NGC 5985)
- Known issues: IRAS SFR may include AGN contamination (Sab with massive bulge). Distance mismatch (65.4 vs 68.0 Mpc) minor.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 02916 | 7.4 Gyr | T2 (Stabilized Spiral) | HIGH | No | — |

---

### Method 89: Di Teodoro & Fraternali (2014) Table A.1 IRAS SFR + SPARC mass → b + Sab evolutionary track calibration → t₅₀ [Fading Sab]
- Data source(s): Di Teodoro & Fraternali (2014) Table A.1: SFR = 0.95 M☉/yr (IRAS), M_HI = 9.21×10⁹ M☉, M_bar = 65.30×10⁹ M☉, D = 47.62 Mpc. SPARC: L[3.6] = 113.642×10⁹ L☉, M_HI = 9.677×10⁹ M☉, D = 50.00 Mpc. M★ ≈ 5.68×10¹⁰ M☉ (log 10.75). Gas fraction ~17%.
- Method: b ≈ 0.23, log sSFR = −10.78, t_dbl ≈ 60 Gyr. Clearly T2 (above all T3 thresholds). **Sab Evolutionary Track completed:** Stabilized (UGC 02916, 37% gas, b = 0.55, 7.4 Gyr) → **Fading (UGC 03205, 17% gas, b = 0.23, 7.9 Gyr)** → Quenched (IC 356, 5.9% gas, b < 0.1, 8.8 Gyr). Gas depletion drives active→quenched transition; each ~halving of gas ≈ 0.5 Gyr aging. **Mass-Morphology Cancellation** with NGC 6195: different masses/types but identical b ≈ 0.23–0.25 → identical 7.9 Gyr.
- Method validated: Single use (third Di Teodoro & Fraternali 2014 application)
- Known issues: IRAS SFR may include minor AGN contamination for early-type. Distance mismatch minor (47.6 vs 50.0 Mpc).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 03205 | 7.9 Gyr | T2 (Fading Sab) | HIGH | No | — |

---

### Method 91: Di Teodoro & Fraternali (2014) Table A.1 SFR + Noordermeer B-R + Hameed Hα EW + SPARC → Scale Paradox → t₅₀ [Blue Early-Type]
- Data source(s): Di Teodoro & Fraternali (2014) Table A.1: SFR = 0.48 M☉/yr, M_HI = 3.81×10⁹ M☉, M_bar = 12.87×10⁹ M☉, M★ = 7.54×10⁹ M☉, D = 25.9 Mpc. Noordermeer & van der Hulst (2007) Table A3: B-R = 1.19±0.05. Hameed et al. Hα Table 3: EW(Hα+[NII]) = 38.4 Å (strong — >20 Å = robust SF). SPARC: L[3.6] = 13.266×10⁹ L☉, M_HI = 4.370×10⁹ M☉, D = 20.70 Mpc. M★ ≈ 6.63×10⁹ M☉ (log 9.82). Gas fraction ~66%.
- Method: b ≈ 0.88, log sSFR = −10.20. **"Scale Paradox":** Sa morphology but gas/activity metrics typical of late-type dwarfs. **Resolution:** at low mass (log M★ ≈ 9.8), bulge is prominent enough for visual Sa classification but too shallow to suppress extended disk. **"Blue Early-Type Galaxy"** class: older compact bulge + massive star-forming disk. **Morphological Quenching requires mass** — low-mass early-types can remain active. Two-component synthesis: disk (b ≈ 0.88, ~6.5 Gyr pull) + Sa bulge (~8+ Gyr pull) → weighted 7.0 Gyr. EW(Hα) = 38.4 Å confirms genuinely active (not AGN mirage).
- Method validated: Single use (fourth Di Teodoro & Fraternali 2014 application)
- Known issues: Distance mismatch (20.7 vs 25.9 Mpc). Two-component age weighting is approximate.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 03580 | 7.0 Gyr | T2 (Blue Early-Type) | HIGH | No | — |

---

### Method 112: Rahvar+2014 Table II M_B + SPARC M★ → M★/L_B override (no direct SFR) + Fuel Retention Gradient → t₅₀ [Mass-to-Light Bypass]
- Data source(s): Rahvar et al. (2014) Table II: M_B = −17.83, (B-V) = 0.65, D = 19.8 Mpc. SPARC: L[3.6] = 1.397×10⁹ L☉, M_HI = 0.809×10⁹ M☉, D = 18.00 Mpc, i = 89° (edge-on). M★ ≈ 6.99×10⁸ M☉ (log 8.84). Gas fraction ~116%. No SFR found in HaGS.
- Derived: L_B ≈ 1.4×10⁹ L☉, M★/L_B ≈ 0.5 (vigorous — young O/B stars dominate light).
- Method: **No direct SFR available.** M★/L_B = 0.5 proves b > 1.0 (fading populations have M★/L_B > 1.0). (B-V) = 0.65 is edge-on dust artifact (i=89°) — true color much bluer. Pro proposed 7.0 Gyr (LOW confidence) → **overridden.** Fuel Retention Gradient: UGC 05999 (120%, 5.5) → UGC 06667 (116%, 5.7) → NGC 3432 (113%, 5.8).
- Method validated: Single use
- Known issues: Pro's 7.0 Gyr overridden. No direct SFR — M/L proxy only. Edge-on dust renders optical color useless.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06667 | 5.7 Gyr | T2 (Vigorous Edge-On) | HIGH | No | — |

---

### Method 113: Di Teodoro & Fraternali (2014) Table A.1 SFR + Noordermeer spectroscopy (Hα absorption) + bulge decomposition (Lb/Ld) → Morphological Quenching → t₅₀ [Fossil S0]
- Data source(s): Di Teodoro & Fraternali (2014) Table A.1: SFR(FIR) = 0.24 M☉/yr. Noordermeer et al. (2007): central Hα in ABSORPTION (dead core, no ionizing radiation). Noordermeer & van der Hulst bulge decomposition: Lb/Ld = 3.91 (bulge 4× disk mass), Sérsic n = 5.5. OpenNGC/HyperLEDA: B-V = 0.85. SPARC: L[3.6] = 73.407×10⁹ L☉, M_HI = 5.030×10⁹ M☉, D = 29.30 Mpc. M★ ≈ 3.67×10¹⁰ M☉ (log 10.56). Gas fraction ~14%.
- Derived: log sSFR = −11.2, b ≈ 0.09, t_dbl ≈ 153 Gyr.
- Method: **Passes ALL T3 criteria:** log sSFR < −11.0, b < 0.1, Hα absorption (dead core). **Morphological Quenching:** Lb/Ld = 3.91 → deep central potential → extreme shear → Toomre Q >> 1 → gas dynamically stable despite 14% retention. Same mechanism as NGC 1167 (S0) but less massive. Interpolation: NGC 2985 (19%, b=0.16, 8.0) → NGC 3900 (14%, b=0.09, 8.6) → IC 356 (5.9%, b<0.1, 8.8).
- Method validated: Fifth Di Teodoro & Fraternali (2014) application
- Known issues: Distance mismatch (22.5 vs 29.3 Mpc) between DT14 and SPARC.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3900 (UGC 06786) | 8.6 Gyr | T3 (Fossil S0) | HIGH | No | — |

---

### Method 123: Liu+2024 Table D2 log M★ + log SFR(Hα/FUV) → sSFR → BC03 → t₅₀ [Phoenix Protocol — Identity Contamination Recovery]
- Data source(s): Liu et al. (2024) MNRAS, "Deep HI mapping of M106 group with FAST," Table D2: log M★ = 9.33, log SFR(Hα) = −0.84, log SFR(FUV) = −0.73. D = 7.2 Mpc (TRGB). Cross-ID: UGC 07151 = NGC 4144 (SABcd, Dec ~+46°).
- Derived: sSFR(Hα) = 6.76×10⁻¹¹, sSFR(FUV) = 8.71×10⁻¹¹, adopted midpoint log sSFR = −10.11, τ_dbl ≈ 12.9 Gyr.
- Method: **Phoenix Protocol applied.** Original Excel provenance was Identity-Contaminated: PGC 38218 = NGC 4140 (UGC 7063, Dec ~+01°), NOT UGC 7151 (Dec ~+46°). Durbala+2020 catalog doesn't cover Dec > 36.4° — footprint violation. Original data belonged to wrong galaxy. **De Novo Provenance** established via Liu+2024 (correct coordinates). Excel value (6.0 Gyr) was "accidentally correct" — Stopped Clock coincidence. New protocol: Identity Contamination is distinct failure mode (right data, wrong object).
- Method validated: Second Liu+2024 application (after UGC 07399)
- Known issues: MEDIUM confidence (sSFR proxy + edge-on caveat). Original provenance expunged.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07151 (NGC 4144) | 6.0 Gyr | T2 | MEDIUM | No | — |

---

### Method 127: Liu+2024 Table D2 log M★ + log SFR(Hα/FUV) → sSFR + b-parameter → t₅₀ [Clean Verification — M106 Group]
- Data source(s): Liu et al. (2024) Table D2: log M★ = 8.56 (M★ = 3.63×10⁸ M☉), log SFR(Hα) = −1.30 (0.0501 M☉/yr), log SFR(FUV) = −1.17 (0.0676 M☉/yr). Parodi+2018 Table B.1 confirms UGC 7399 = NGC 4288. Coordinates: RA 12:20:38, Dec +46:17:30. SFRs from Karachentsev & Kaisina (2013) — integrated flux, NOT surface density.
- Derived: Adopted midpoint log sSFR = −9.79, b ≈ 2.3, τ_dbl ≈ 6.2 Gyr.
- Method: Clean verification — all protocols passed (Identity, Aperture Fallacy, Ghost Data, Footprint). b > 2 exerts "Youthful Pull" below 6.0 Gyr equilibrium. 5.4 Gyr sits in Active Transition Zone. Liu+2024 confirmed as Platinum Source for M106/UMa sector.
- Method validated: Third Liu+2024 application
- Known issues: None. Deep Think upgraded confidence to HIGH.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07399 (NGC 4288) | 5.4 Gyr | T2 (Active) | HIGH | No | — |

---

### Method 135: LVGDB FUV SFR (Karachentsev & Kaisina 2013) + SPARC mass → Kroupa correction (÷1.59) + Scale-Height Extinction Differential → Bivariate Interpolation → t₅₀ [Edge-On FUV Baseline]
- Data source(s): Karachentsev & Kaisina (2013) LVGDB: SFR(FUV) = 0.11 M☉/yr (K98/Salpeter), SFR(Hα) = 0.049 M☉/yr. HALOGAS: SFR = 0.039 M☉/yr. SPARC: L[3.6] = 1.255×10⁹ L☉, M_HI = 0.642×10⁹ M☉, D = 6.50 Mpc (TRGB), i ≈ 90°. M★ ≈ 6.28×10⁸ M☉ (log 8.80). Gas fraction ~102%. Cross-ID: NGC 5023 = PGC 45849. GHOSTS HST program confirms mixed-age populations (Tier A partial, no global SFH tabulated).
- Derived: FUV Kroupa = 0.11 ÷ 1.59 = 0.069 M☉/yr. t_dbl = 9.1 Gyr, b = 1.51. FUV/Hα ratio = 2.2× (edge-on dust signature).
- Method: **Scale-Height Extinction Differential:** at i≈90°, Hα (from O-stars in thin midplane) catastrophically dust-attenuated; FUV (from B-stars migrated to thicker scale height) partially bypasses midplane dust. FUV/Hα > 2× is diagnostic. Pro initially used HALOGAS SFR → false b < 1.0 ("declining") and t₅₀ = 7.0 Gyr → **corrected** to FUV baseline. Bivariate: NGC 3104 (102%, b=1.16, 5.8) → UGC 08286 (102%, b=1.51, 5.6) → UGC 04499 (142%, b=1.46, 5.4). Same gas as NGC 3104 but higher b → younger; same b as UGC 04499 but less gas → older.
- Method validated: Single use (first LVGDB FUV application)
- Known issues: Pro's 7.0 Gyr (HÎ±-based) overridden. FUV itself partially attenuated at i≈90° → b=1.51 is lower limit. NGC 3104 anchor not formally checksummed.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 08286 (NGC 5023) | 5.6 Gyr | T2 (Vigorous Edge-On Disk) | HIGH | No | — |

---

### Method 139: Hallenbeck+2016 Table 1 log M★ + log SFR → sSFR → b-parameter gradient → t₅₀ [HIghMass — Data Correction + Birthrate Inversion Fix]
- Data source(s): Hallenbeck et al. (2016) Table 1: log M★ = 10.09, log SFR = 0.56. Cross-checked: Hallenbeck et al. (2014) Table 1: log M★ = 10.09, log SFR = 0.63, log sSFR = −9.46 (text: "sSFR = 3.5×10⁻¹⁰ yr⁻¹ for UGC 9037"). Cross-ID: UGC 09037.
- Derived: log sSFR = −9.53, sSFR ≈ 2.95×10⁻¹⁰ yr⁻¹, b ≈ 4.1 (Vigorous).
- Method: **Data Correction:** Astra's sSFR was ~2× too low (~1.4×10⁻¹⁰ vs verified ~3×10⁻¹⁰), yielding wrong 6.2 Gyr. **Birthrate Inversion Fix:** Pro proposed 5.5 Gyr — Deep Think overruled because b = 4.1 galaxy cannot be OLDER than b = 2.3 galaxy (UGC 07399, 5.4 Gyr). Higher b → younger t₅₀ (monotonic). b-Parameter Grid formalized: b < 1 → 7–9 Gyr (T3); b ~ 1 → 6.0 Gyr (Stability Block); b ~ 2 → 5.4 Gyr (Active); b > 3 → 5.0 Gyr (Active T2 Floor). HIghMass delayed-evolution population.
- Method validated: Single use (first Hallenbeck/HIghMass application)
- Known issues: Astra's data wrong. Pro's 5.5 Gyr overridden (Birthrate Inversion).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 09037 | 5.0 Gyr | T2 (Active T2 Floor / HIghMass) | MEDIUM | No | — |

---

### Method 143: MANGROVE/GLADE catalog log M★ + log SFR + SPARC gas fraction → τ-model (massive spiral) → t₅₀ [Cosmic Downsizing Anchor]
- Data source(s): MANGROVE/GLADE: log M★ = 11.29 (~1.95×10¹¹ M☉), log SFR = 0.52 (~3.31 M☉/yr), r_SFR = 4 (compilation). SPARC: L[3.6] = 374.3×10⁹ L☉, M_HI = 13.3×10⁹ M☉, i ≈ 90° (edge-on), Vflat = 269.4 km/s. M★ cross-check: SPARC log ≈ 11.27 vs MANGROVE 11.29 (0.02 dex agreement). Cross-ID: UGC 11455 = PGC 63286. Scd spiral at ~84 Mpc.
- Derived: Gas fraction ~7% (nearly exhausted), t_dbl ≈ 59 Gyr, b ≈ 0.23 (declining engine).
- Method: **Cosmic Downsizing.** At log M★ > 11, smooth secular τ-model evolution dominates (no burst traps). b = 0.23 maps to t₅₀ ≈ 10.3 Gyr. **Gas Fraction Veto:** 7% gas independently vetoes T2 — galaxy has converted 93% of baryons to stars. Even 4× dust correction (edge-on) only reaches b ≈ 0.92, still T3-adjacent. First massive spiral (log M★ > 11) in recent queue. τ-model b→t₅₀ mapping is analytically robust at high mass (no stochastic bursts).
- **Mass-Scale Protocol:** For log M★ > 10.5, prefer τ-model over bivariate dwarf interpolation.
- Method validated: Single use (first MANGROVE/GLADE application)
- Known issues: Catalog-level SFR (r_SFR = 4), not dedicated measurement. Edge-on geometry. MEDIUM confidence (Tier C).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 11455 | 10.3 Gyr | T3 (Massive Depleted Disk) | MEDIUM | No | — |

---

### Method 146: Hallenbeck+2014 HIghMass Table 1 log M★ + log SFR (column corrected) → sSFR → b → False Precision Twin (UGC 09992) → Stability Block → t₅₀ [High Spin Delayed Assembly]
- Data source(s): Hallenbeck et al. (2014) HIghMass Table 1: log M★ = 10.46, log SFR = 0.40, GF = 1.41 (gas fraction M_HI/M★), log M_HI = 10.53, D = 98.2 Mpc, R₂₅ = 40.0 kpc. i ≈ 87° (edge-on). High spin parameter λ ≈ 0.15. Cross-ID: UGC 12506.
- Derived: M★ ≈ 2.88×10¹⁰ M☉, SFR = 2.51 M☉/yr, sSFR ≈ 8.7×10⁻¹¹ yr⁻¹ (log = −10.06), b ≈ 1.20.
- Method: **Column Mixup Corrected:** Astra read GF = 1.41 as SFR → wrong sSFR (×1.8 too low) → wrong 7.8 Gyr. Pro corrected to log SFR = 0.40. **Birthrate Inversion:** Pro's 6.6 Gyr violated b > 1 physics (b > 1 cannot produce age older than 6.0 equilibrium). **False Precision Twin Match:** log sSFR = −10.06 is identical to UGC 09992 (−10.07) within 0.01 dex → same bin. **High Spin Delayed Assembly:** Despite massive disk (log M★ = 10.46), high dark matter spin parameter (λ ≈ 0.15) prevented early gas collapse → evolutionarily delayed → behaves like LSB dwarf. Edge-on (87°) creates optical dust trap — ignore visual redness.
- Method validated: Second Hallenbeck/HIghMass application
- Known issues: Astra column error. Pro's 6.6 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 12506 | 6.0 Gyr | T2 (Stability Block / High Spin) | MEDIUM-HIGH | No | — |

---



## IR / radio SFR

**Conversion:** IR or radio observations are used as the SFR tracer. We:
1) convert the tracer to SFR using the source’s stated calibration (and correct IMF if required),
2) combine with \(M_\*\) (typically SPARC) to compute sSFR, \(t_{\rm dbl}\), and \(b\),
3) map to \(t_{50}\) using the same b/sSFR→age framework documented in the audit (including dust/AGN diagnostics when IR is contaminated).

Radio special case: when a catalog reports SFR for \(M>5\,M_\odot\), we apply the audit’s IMF scaling to total SFR before computing sSFR/b.

### Worked example

**Worked example (Method 120 — radio SFR IMF scaling):**  
        This entry reports a radio SFR calibrated for \(M>5\,M_\odot\). The audit applies a fixed IMF scaling to obtain total SFR before computing sSFR and \(b\), and then maps to the recorded \(t_50\).

        ### Method 120: Wiegert PhD Table 6.9 radio SFR (M>5 M☉) + Tully+1996 Table 5 B-R + SPARC mass → Radio IMF Scaling (×3.44) → t₅₀ [Bar-Driven LSB Burst]
- Data source(s): Wiegert PhD Thesis Table 6.9: SFR(radio, M>5 M☉) = 0.14 M☉/yr (NVSS 1.4 GHz, Condon 1992 calibration). Tully et al. (1996) Table 5: (B-R) = 1.05 (inclination-corrected), LSB with small central bar. SPARC: L[3.6] = 5.298×10⁹ L☉, M_HI = 2.967×10⁹ M☉, D = 18.00 Mpc, i = 49°. M★ ≈ 2.65×10⁹ M☉ (log 9.42). Gas fraction ~112%.
- Derived: **Radio IMF Scaling:** raw 0.14 × 5.5 (Salpeter total) ÷ 1.6 (Kroupa) = 0.48 M☉/yr. Net multiplier ×3.44. b ≈ 2.5, log sSFR = −9.74, t_dbl ≈ 5.5 Gyr.
- Method: Optical (B-R = 1.05) sees quiescent LSB disk; radio sees buried nuclear burst driven by central bar funneling gas. T1 vetoed: t_dbl = 5.5 Gyr > 4.0 + LSB sub-critical disk prevents global collapse. High-b Fuel Retention Gradient: UGC 06446 (279%, 4.6) → UGC 06983 (112%, 5.1) → UGC 06923 (56%, 5.4).
- Method validated: First Wiegert PhD radio application
- Known issues: Radio IMF scaling introduces systematic uncertainty (~factor 2).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06983 | 5.1 Gyr | T2 (Bar-Driven LSB Burst) | HIGH | No | — |

---

### Method entries (verbatim from the audit)

### Method 21: SPARC mass + FIR SFR (Rossa & Dettmar 2003) → sSFR → b → anchor interpolation → t₅₀
- Data source(s): SPARC Table 1 + Rossa & Dettmar (2003) Table 1: SFR_FIR = 0.03 M☉/yr
- Method validated: Single use
- Known issues: FIR SFR = floor in dust-poor systems ("Infrared Trap")
- Correction: Deep Think claimed 0.21 (column misread — was L_FIR/D25² from Table 2). Corrected to 0.03.
- "Super-Thin Paradox": 123% gas fraction but b ≈ 0.26 (dynamically cold i = 89° disk)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 0100 | 7.5 Gyr | T2 (Inefficient Super-Thin) | MEDIUM | No | — |

---

### Method 28: SPARC mass + IRAS FIR SFR (CGS Ho+2011) → sSFR → b → anchor → t₅₀ [Hα tiebreaker]
- Data source(s): SPARC Table 1 + CGS (Ho et al. 2011): S_60μm = 0.854 Jy, S_100μm = 3.79 Jy → L_FIR → SFR ≈ 0.34 M☉/yr. Gentile et al. (2004) confirms Hα+HI rotation curve.
- Derived: log M★ ≈ 10.56, log sSFR = −11.03, b ≈ 0.13
- Method validated: Single use
- Known issues: **Boundary case** — log sSFR = −11.03 is just below T2/T3 threshold (−11.0). FIR-to-SFR calibration uncertainty ~0.2 dex. Classified T2 via Hα tiebreaker (extended Hα emission = not passive).
- Correction: Deep Think claimed flagged g-r color — not usable; IRAS recovery used instead.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 1090 | 8.0 Gyr | T2 (Inefficient Spiral) | MEDIUM | No | — |

---

### Data Gap: NGC 1705
- BCD / starburst dwarf at ~5.1 Mpc. Not in Weisz+2011 Table 2. Qualitative CMD studies exist (Annibali+2003, Tosi+2001) confirming old RGB/AGB populations and SF ≥5 Gyr, but no tabulated cumulative SFH. "≥5 Gyr of SF" does not constrain t₅₀. Hidden T3 cannot be ruled out. Astra's 6.5 Gyr unverifiable.

---

### Method 49: SPARC mass + IRAS fluxes (Chen et al. RAA) → SFR → b → t₅₀ [AGN cold dust diagnostic]
- Data source(s): SPARC Table 1: L[3.6] = 95.268×10⁹ L☉, M_HI = 2.697×10⁹ M☉. Chen et al. (RAA) IRAS: S_60 = 7.131 Jy, S_100 = 23.9 Jy. Derived: M★ ≈ 4.76×10¹⁰ M☉, SFR ≈ 1.2 M☉/yr, b ≈ 0.35, gas fraction 5.7%.
- Method: AGN contamination ruled out via cold dust diagnostic: S_100/S_60 = 3.35 >> 2.0 (AGN threshold). 60/100μm dominated by cold SF dust, not hot torus. Seyfert 1 host. "Starving" not "Choked" — low gas but actively consuming (contrast NGC 5985: more gas, engine off → T3).
- Method validated: Single use
- Known issues: Seyfert 1 AGN. Cold dust diagnostic confirms FIR traces SF, not AGN. Terminal active phase — most gas-depleted T2 in sample (5.7%). Ursa Major cluster strangulation sequence endpoint.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4051 | 7.6 Gyr | T2 (Terminal Active / Starving Spiral) | HIGH | No | — |

---

### Method 55: Tully+1996 Table 5 + Bouquin+2018 GALEX UV → Infrared Supremacy Protocol → t₅₀ [B-K overrides B-R]
- Data source(s): Tully+1996 Table 5 global: B-R = 1.17, B-K' = 3.30. Bouquin+2018: FUV-NUV = 0.61. R-K = 2.13.
- Method: Astra classified T3/8.4 based on B-R = 1.17 alone. **Infrared Supremacy Protocol** (new): When B-R (optical) and B-K' (infrared) disagree on tier, Infrared Wins. B-K' = 3.30 < 3.4 → Strong T2. R-K = 2.13 (bluer than Active benchmark ~2.2). UV confirms (FUV-NUV = 0.61 < 0.9). Red B-R attributed to bar/dust geometry.
- Method validated: Single use
- Known issues: Original T3/8.4 **REJECTED**. Barred spiral (SBbc) — bar population and dust explain red B-R without requiring old mass-weighted age. All three independent metrics (B-K', R-K, UV) unanimously T2.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 4389 | 7.2 Gyr (range 6.0–8.0) | T2 (Infrared Supremacy) | MEDIUM+ | No | — |

---

### Method 83: Gavilán+2013 Table 3 B-V + SFR + metallicity + SPARC mass → Metallicity Trap + Gas Veto → t₅₀ [Inefficient Gas-Rich Dwarf]
- Data source(s): Gavilán et al. (2013) Table 3: B-V = 0.44±0.04, SFR = 0.007218 M☉/yr, log M★ = 8.62, 12+log(O/H) = 8.02±0.03 (20% solar). SPARC: L[3.6] = 1.308×10⁹ L☉, M_HI = 0.477×10⁹ M☉, D = 10.40 Mpc. SPARC M★ = 6.54×10⁸ M☉ (adopted over Gavilán's 4.17×10⁸). Gas fraction ~73%.
- Method: Corrected b = 0.15 (forming at 15% of lifetime average). log sSFR = −10.96 (just above −11.0 threshold). t_dbl = 91 Gyr. **"Metallicity Trap":** B-V = 0.44 looks blue/young but 12+log(O/H) = 8.02 (20% solar) → metal-poor stars burn hotter/bluer at same age. Same B-V = 0.44 as NGC 2998 (b ~ 0.7, near-solar Z) but wildly different activity. **Gas Veto prevents T3:** 73% gas = abundant fuel, just processing slowly. T3 requires gas depletion or strangulation. This is inefficiency, not quiescence. Aligned with F571-8 (b = 0.12, 8.0 Gyr) with slight offset for higher gas retention.
- Method validated: Single use (third Gavilán+2013 application)
- Known issues: Mass discrepancy between Gavilán (log 8.62) and SPARC (log 8.82) — SPARC adopted for consistency. Pro initially used Gavilán mass (b = 0.24).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 02023 (DDO 25) | 7.8 Gyr | T2 (Inefficient Gas-Rich Dwarf) | HIGH | No | — |

---

### DATA GAP: UGC 02259
- **Tier: T2 CONFIRMED** via Grand Design spiral morphology + Hα velocity field.
- Grand Design structure in dwarf requires Gyr of stable rotational support → rules out T1 (young/chaotic).
- Hα velocity field requires active HII regions throughout disk → rules out T3 (passive).
- **T2 by elimination.**
- **Missing:** Tabulated SFR or integrated color. No numeric t₅₀ lockable.
- Physics estimate ~6.5 Gyr (not locked). SPARC: L[3.6] = 1.725×10⁹ L☉, M_HI = 0.494×10⁹ M☉, D = 10.50 Mpc. Gas fraction ~57%.
- Path to resolution: need tabulated Hα SFR, FUV+IR SFR, or integrated color.

---

### Method 88: IRAS RBGS Table 1 far-IR fluxes + Cold Dust Diagnostic (S100/S60) + SPARC gas fraction → Cirrus Correction → t₅₀ [Dust-Bypassed Anchor]
- Data source(s): IRAS RBGS Table 1 (F04025+6940): S60 = 6.77 Jy, S100 = 28.33 Jy, S12 = 1.08, S25 = 0.84, log L_IR = 9.97, D = 15.86 Mpc. SPARC: L[3.6] = 259.518×10⁹ L☉, M_HI = 7.678×10⁹ M☉, D = 16.50 Mpc. M★ ≈ 1.30×10¹¹ M☉ (log 11.11). Gas fraction ~5.9%. Cross-ID: UGC 02953 = IC 356 = Arp 213 = PGC 14508.
- Method: Galaxy lies in "IC 342 zone" (l,b ≈ 138°, +13°) — heavy MW foreground extinction makes optical/UV unreliable. **Far-IR bypasses extinction.** Naive IRAS SFR ≈ 1.0 M☉/yr → b ≈ 0.11 (borderline T2/T3). **Cold Dust Diagnostic:** S100/S60 = 4.18 → very cold dust (~20K) → **cirrus-dominated** (old-star heating, not SF). True b < 0.1. Gas 5.9% matches quenched population. Sab morphology + massive (log 11.1) → Downsizing. **Zone of Avoidance Protocol** established for IC 342 lane. Confidence upgraded from LOW-MEDIUM (Astra) to MEDIUM-HIGH via IRAS recovery.
- Method validated: Single use (second Cold Dust Diagnostic application after NGC 6195)
- Known issues: MEDIUM-HIGH confidence — IR calibration assumptions; cirrus correction is diagnostic, not quantitative; no direct Hα or resolved SFH.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| IC 356 (UGC 02953, Arp 213) | 8.8 Gyr | T3 (Quenched Spiral) | MEDIUM-HIGH | No | — |

---

### Method 90: IRAS RBGS Table 1 + AGN Warm Dust correction (S100/S60) + hidden H₂ via bar-driven phase transition → t₅₀ [AGN/Molecular Correction]
- Data source(s): IRAS RBGS Table 1: S60 = 6.41 Jy, S100 = 9.55 Jy, S12 = 0.44, S25 = 1.36, L_IR ≈ 1.8×10¹⁰ L☉. SPARC: L[3.6] = 101.336×10⁹ L☉, M_HI = 2.675×10⁹ M☉, D = 28.70 Mpc. M★ ≈ 5.07×10¹⁰ M☉ (log 10.70). Gas fraction (HI) ~5.3%. Cross-ID: UGC 03546 = NGC 2273. Classification: SB(r)a, Seyfert 2.
- Method: Raw b ≈ 0.54 from naive IRAS SFR. **Warm Dust Diagnostic:** S100/S60 = 1.49 → very warm → **AGN torus heating** inflates 60μm. Corrected SFR ≈ 1.0 M☉/yr → true b ≈ 0.25. **Hidden H₂ resolves gas paradox:** HI only 5.3% (would suggest T3), but SB(r)a bar acts as gravitational conveyor → funnels HI to center → HI→H₂ conversion → dense nuclear ring SF invisible to 21cm. **Low HI ≠ Quenched when bars are present.** Contrast with IC 356 (Sab, no bar, same HI%, T3) — morphology determines fate. **Warm Dust Diagnostic codified:** S100/S60 < 2.0 = AGN, 2.0–3.5 = SF-dominated, > 4.0 = cirrus.
- Method validated: Single use
- Known issues: AGN correction is factor-of-2 estimate. H₂ mass not directly measured. S100/S60 diagnostic empirically calibrated.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2273 (UGC 03546) | 7.9 Gyr | T2 (Nuclear Ring Sa) | HIGH | No | — |

---

### Method 99: van Zee (2001) Table 1 (B-V)₀ + Table 2 SFR + ⟨SFR⟩_past → dual birthrate validation → t₅₀ [Gold Standard Anchor]
- Data source(s): van Zee (2001) Table 1: (B-V)₀ = 0.37±0.03, (U-B)₀ = −0.30±0.04. Table 2: SFR(Hα) = 0.046 M☉/yr, ⟨SFR⟩_past = 0.038 M☉/yr, D = 8.23 Mpc, log M_HI = 8.60, log t_g = 9.94. SPARC: L[3.6] = 1.123×10⁹ L☉, M_HI = 0.574×10⁹ M☉, D = 9.40 Mpc. M★ ≈ 5.62×10⁸ M☉ (log 8.75). Gas fraction ~102%. Cross-ID: UGC 05414 = NGC 3104.
- Method: **"Gold Standard" — dual birthrate validation.** SPARC-based b ≈ 1.13 (sSFR × T_universe) vs van Zee internal b ≈ 1.21 (SFR/⟨SFR⟩_past). **Excellent convergence** from completely independent methods → validates entire b-parameter framework. **Only galaxy in sample with direct literature b measurement.** Galaxy forming ~21% above historical average. Gas 102% barely passes T1 gas threshold but t_dbl = 12 Gyr >> 4 Gyr → T2. (B-V)₀ = 0.37 → BC03 match at ~5.8 Gyr (mildly rising SFH bluer than 6.0 Gyr steady-state). **Color-Age anchor:** 0.37 → 5.8, 0.40 → 6.0, 0.45 → 6.5.
- Method validated: Second van Zee (2001) application (after UGC 00634/DDO 7)
- Known issues: Distance mismatch (van Zee 8.23 vs SPARC 9.40 Mpc) — minor.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3104 (UGC 05414) | 5.8 Gyr | T2 (Steady State Dwarf) | HIGH | No | — |

---

### Method 115: DATA GAP — UGC 06818 (PGC 37550)
- Data source(s): Excel sheet entry: sSFR ≈ 5.3×10⁻¹¹ yr⁻¹ → log sSFR = −10.28. **Source data (M★, SFR) NOT FOUND.** Likely source: HECATE (Kovlakas+2021) / GSWLC-2 (Salim+2018) crossmatch columns (logM_GSW, logSFR_GSW). Astra could not locate original catalog row. Pro verified methodology but not inputs.
- Method: **Ghost Data Protocol.** Derived value (log sSFR = −10.28) present but primary inputs (M★, SFR) and source locator (catalog row) missing. "Methodology ≠ Data" — correct calculation on invisible numbers does not constitute verified classification. Deep Think ruling: HOLD until source retrieved. **Pre-approved** for 7.0 Gyr (T2, MEDIUM) once raw numbers confirmed and log sSFR ∈ [−10.0, −10.6].
- Method validated: N/A (data gap)
- Known issues: Source unrecoverable. Cannot lock.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06818 | — (pre-approved 7.0 Gyr) | DATA GAP | — | N/A | Ghost Data — source traceability required |

---

### Method 117: James+2004 HaGS Table 3 SFR + Sanders (1998) B-V + SPARC mass → IMF-corrected b + Iso-Age Diagonal → t₅₀ [Late-Stage Starburst]
- Data source(s): James et al. (2004) HaGS Table 3: SFR(Salpeter) = 0.4532±0.1063 M☉/yr, EW = 3.0 nm, D = 19.8 Mpc. Sanders (1998): (B-V) = 0.42. SPARC: L[3.6] = 2.890×10⁹ L☉, M_HI = 0.809×10⁹ M☉, D = 18.00 Mpc. M★ ≈ 1.45×10⁹ M☉ (log 9.16). Gas fraction ~56%.
- Derived: SFR distance-rescaled + Kroupa = 0.234 M☉/yr, b ≈ 2.2, log sSFR = −9.79, t_dbl ≈ 6.2 Gyr.
- Method: **Pro flagged T1/3.7 Gyr (uncorrected Salpeter) — OVERRIDDEN.** IMF correction: t_dbl = 6.2 Gyr → fails T1. **Late-Stage Starburst physics:** high b (2.2) on depleted backbone (56% gas) = rejuvenation, not primordial burst. (B-V) = 0.42 confirms old+young mix (pure burst < 0.30, pure old > 0.60). Iso-Age Diagonal: UGC 04499 (142%, b=1.46, 5.4) matches UGC 06923 (56%, b=2.2, 5.4) — opposing vectors cancel.
- Method validated: Third HaGS application
- Known issues: Pro's T1/3.7 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 06923 | 5.4 Gyr | T2 (Late-Stage Starburst) | HIGH | No | — |

---

### Method 121: Wiegert+2013 Table 9 radio SFR (M>5 M☉) + Tully+1996 B-R + SPARC mass → Radio IMF Scaling (×3.44) + Bivariate Interpolation → t₅₀ [Edge-On Steady-State]
- Data source(s): Wiegert et al. (2013) Table 9: SFR(M>5 M☉) = 0.04 M☉/yr (radio continuum). Tully et al. (1996): (B-R) = 0.70 (remarkably blue for edge-on). SPARC: L[3.6] = 3.585×10⁹ L☉, M_HI = 1.214×10⁹ M☉, D = 18.00 Mpc, i = 80°. M★ ≈ 1.79×10⁹ M☉ (log 9.25). Gas fraction ~68%.
- Derived: Radio IMF: 0.04 × 3.44 = 0.138 M☉/yr (Kroupa). b ≈ 1.06, log sSFR = −10.11, t_dbl ≈ 13 Gyr.
- Method: **Pro flagged T1/3.8 Gyr (uncorrected Salpeter ×5.5) — OVERRIDDEN.** Correct scaling ×3.44 gives t_dbl = 13 Gyr → safely T2. **Bivariate Interpolation:** bounded by NGC 2552 (67%, b=1.35, 5.7) on activity axis and NGC 1003 (~50%, b=1.00, 6.0) on gas axis → 5.9 Gyr. (B-R) = 0.70 remarkably blue for i=80° → confirms intrinsic disk is active.
- Method validated: Second Wiegert radio application (after UGC 06983)
- Known issues: Pro's T1/3.8 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07089 | 5.9 Gyr | T2 (Steady-State Edge-On) | HIGH | No | — |

---

### Method 122: DATA GAP — UGC 07125 (Twin Degeneracy)
- Data source(s): SPARC: L[3.6] = 2.712×10⁹ L☉, M_HI = 4.629×10⁹ M☉, D = 19.80 Mpc, i = 90° (perfectly edge-on), Sm morphology. M★ ≈ 1.36×10⁹ M☉ (log 9.13). Gas fraction ~341%. Cross-IDs: MCG+06-27-014, PGC 38567, HIJASS J1208+36. **No SFR found** in any catalog (HÎ±, FUV, radio, IR).
- Method: **T2 confirmed via Structural Veto** — diffuse Sm morphology with sub-critical gas density (Toomre Q > 1) prevents T1 runaway collapse despite 341% gas. **Age uncalibrated due to Twin Degeneracy:** structurally identical galaxies UGC 05716 (372%, b=1.1, 5.3 Gyr) and UGC 05829 (363%, b=2.33, 4.5 Gyr) span 800 Myr range. Without measuring b, cannot resolve. Edge-on extinction wall (i=90°) prevents optical diagnostics and M★/L_B override.
- Resolution: Radio continuum or far-IR SFR measurement required. T1 threshold: SFR(Kroupa) > 0.34 M☉/yr (structurally unlikely).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07125 | — (T2 confirmed, age uncalibrated) | DATA GAP | — | N/A | Twin Degeneracy — no SFR, 4.5–5.3 Gyr range |

---

### Method 125: SFRS Table 2 UV+IR SFR + HaGS Table 3 (overcorrection demonstrated) + SPARC mass → Multiwavelength Harmony → t₅₀ [Face-On Barred Starburst]
- Data source(s): SFRS Table 2: log SFR_tot = −0.97 → 0.107 M☉/yr (UV+IR synthesis). SFRS Table 1: FUV = 15.188, NUV = 14.946, FUV-NUV = 0.24, E(B-V) = 0.034. James+2004 HaGS Table 3: SFR = 0.2835 M☉/yr (overcorrected by fixed A=1.1 mag). SPARC: L[3.6] = 1.753×10⁹ L☉, M_HI = 1.388×10⁹ M☉, D = 13.10 Mpc, i = 30° (face-on). M★ ≈ 8.77×10⁸ M☉ (log 8.94). Gas fraction ~158%. Cross-ID: NGC 4204, SFRS ID 183, SBdm.
- Derived: HaGS uncorrected: 0.2835/2.75 = 0.103 M☉/yr ≈ SFRS 0.107 → **perfect multiwavelength harmony.** SPARC-scaled SFR ≈ 0.18 M☉/yr, b ≈ 2.8, log sSFR = −9.69, t_dbl ≈ 4.87 Gyr.
- Method: Face-on (i=30°) eliminates dust — all tracers agree when HaGS fixed extinction removed. **HaGS Overcorrection Trap:** fixed A(Hα)=1.1 inflates SFR by 2.75× for dust-poor systems. Bar-driven starburst in diffuse SBdm disk. Bivariate: UGC 06446 (279%, 4.6) → UGC 07261 (158%, 4.8) → UGC 06983 (112%, 5.1).
- Method validated: First SFRS application
- Known issues: None significant.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07261 (NGC 4204) | 4.8 Gyr | T2 (Face-On Barred Starburst) | HIGH | No | — |

---

### Method 129: Hunter+2010 Table 6 SFR(FUV) + Table 2 (B-V)₀ + SPARC mass → Double Error Correction + Bivariate Interpolation → t₅₀ [Column Misread + Gas Denominator Trap]
- Data source(s): Hunter et al. (2010) Table 2: (B-V)₀ = 0.29±0.02 (BLUE), M_V = −14.85. Table 6: log SFR_FUV = −1.89 → 0.013 M☉/yr. Elmegreen & Hunter (2015): R_D = 0.87±0.03 kpc (disk scale length — NOT color). SPARC: L[3.6] = 0.109×10⁹ L☉, M_HI = 0.169×10⁹ M☉, D = 4.97 Mpc. M★ ≈ 5.5×10⁷ M☉ (log 7.74). Gas fraction ~310%. Cross-ID: DDO 126.
- Derived: SFR(Kroupa) = 0.0081 M☉/yr, b ≈ 2.03, log sSFR = −9.83, t_dbl ≈ 6.8 Gyr.
- Method: **Double error correction.** Error 1: Astra misread R_D = 0.87 kpc (disk scale length) as (B-V) = 0.87 → false T3/8.6 Gyr. Error 2: Pro used M_HI instead of M★ in sSFR denominator (Gas Denominator Trap) → b ≈ 0.9 (wrong). Corrected b = 2.03. True (B-V)₀ = 0.29 confirms vigorous blue dwarf. b/t_dbl interlock: b = 13.8/t_dbl must hold. Bivariate: UGC 05829 (363%, b=2.33, 4.5) → DDO 126 (310%, b=2.03, 4.9) → UGC 05716 (372%, b=1.1, 5.3).
- Method validated: First Hunter+2010 application
- Known issues: Both Astra and Pro errors corrected by Deep Think.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07559 (DDO 126) | 4.9 Gyr | T2 (Vigorous Sub-Critical Dwarf) | HIGH | No | — |

---

### Method 133: Dunn (2015) Table 6 UBVR at R₂₅ → B-V = 0.38±0.05 → Interpolation Principle → LSB Stability Block → t₅₀
- Data source(s): Dunn (2015) Table 6 "UBVR apparent magnitudes at R25" (extinction-corrected): U=13.07±0.05, B=13.40±0.05, V=13.02±0.02, R=12.55±0.01. R₂₅ = 52.29″. B-V = 13.40 − 13.02 = 0.38±0.05.
- Method: **Interpolation Principle.** B-V = 0.38 falls between UGC 05005 (B-V=0.35, 6.0 Gyr) and UGC 00634 (B-V=0.39, 6.0 Gyr) — both already at 6.0 Gyr. Assigning 5.9 Gyr (Astra/Pro proposal) would create "Redder is Younger" inversion relative to UGC 05005. Error bar (±0.05) spans 0.33–0.43, covering entire Blue LSB Block — 0.1 Gyr precision is statistically invalid. No Second Axis (sSFR) to justify deviation from Stability Block. Standardized to 6.0 Gyr.
- **NEW PROTOCOL: Interpolation Principle** — When candidate color falls between two existing anchors at same t₅₀, snap to anchor value. Do not interpolate. Exception: Second Axis provides independent confirmation.
- Method validated: Single use (first Dunn application)
- Known issues: Astra/Pro's 5.9 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 07690 | 6.0 Gyr | T2 (LSB Stability Block / Interpolation) | HIGH | No | — |

---

### Method 141: Gavilán+2013 MNRAS 434, Table 3 B-V + SFR/M★ → sSFR → b-parameter → Local Maximum Avoidance → Stability Block → t₅₀
- Data source(s): Gavilán et al. (2013) Table 3: B-V = 0.37±0.04, SFR = 0.0062 M☉/yr, log M★ = 7.86, D = 8.6 Mpc. Cross-ID: UGC 09992 = PGC 55809. Magellanic Irregular.
- Derived: M★ = 7.24×10⁷ M☉, sSFR = 8.6×10⁻¹¹ yr⁻¹ (log = −10.07), b ≈ 1.19.
- Method: Astra proposed 6.2 Gyr as "compromise" between B-V proxy (~4.8–5.4) and sSFR proxy (~6–7). Deep Think rejected: **Local Maximum Avoidance** — B-V = 0.37 sits between UGC 05005 (0.35, 6.0) and UGC 07690 (0.38, 6.0). Assigning 6.2 creates a "spike" — physically invalid local maximum in monotonic sequence. Also: b ≈ 1.19 (> 1) implies youthful pull, opposing older age. Standardized to 6.0 Gyr (Stability Block).
- Method validated: First Gavilán application
- Known issues: Astra's 6.2 Gyr overridden.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 09992 | 6.0 Gyr | T2 (Stability Block) | MEDIUM-HIGH | No | — |

---

### Method 142: van Zee (2001) AJ 121, Table 1 (B-V)₀ + (U-B)₀ → Second Axis Override (U-B Youth Flag) → Active T2 Cohort → t₅₀
- Data source(s): van Zee (2001) AJ 121, 2003, Table 1: (B-V)₀ = 0.39±0.05, (U-B)₀ = −0.27±0.08, m_B = 13.62±0.06, i = 31° (face-on). Table 2: D = 12.22 Mpc, SFR = 0.123 M☉/yr, ⟨SFR⟩_past = 0.064 M☉/yr. Morphology: SB(s)m. Colors corrected for Galactic extinction, NOT internal. Cross-ID: UGC 10310 = DDO 204 = Arp 2.
- Method: **Twin Color Problem:** B-V = 0.39 matches UGC 00634 (6.0 Gyr Stability Block anchor). Astra proposed 4.6 Gyr. **Second Axis Override:** U-B = −0.27 is extremely blue, proving young hot stars that UGC 00634 lacks. Twin Anchor only applies when ALL color indices match. SFR/⟨SFR⟩_past ≈ 1.9 (elevated current SF). Deep Think: False Precision prevents 4.6 vs 5.0 distinction (Δ(B-V) = 0.02 from NGC 5204 is < σ). Standardized to 5.0 Gyr (Active T2 Cohort).
- **Protocol clarified: Second Axis Override** — Multi-band evidence (U-B, NUV, verified sSFR) can override single-band Twin Anchor.
- Method validated: First van Zee (2001) application
- Known issues: Astra's 4.6 Gyr overridden. No verified sSFR (M★ not tabulated).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| UGC 10310 (DDO 204 / Arp 2) | 5.0 Gyr | T2 (Active Cohort) | MEDIUM-HIGH | No | — |

---



## Compiled / SPARC-derived

**Conversion:** These methods use SPARC as the mass backbone and compile one or more external SFR/activity indicators.  
The deterministic part is always:
\[
M_\* \text{ (from SPARC backbone)}\ \&\ \mathrm{SFR}\ \Rightarrow\ \mathrm{sSFR},\ t_{\rm dbl},\ b.
\]
The final \(t_{50}\) is assigned using the documented anchor interpolation / protocol logic (e.g., “steady-state”, “winding down”, “downsizing penalty”), with every chosen source and tie-breaker recorded in the audit entry.

### Worked example

**Worked example (Method 14 — SPARC backbone + catalog SFR):**  
        This entry uses SPARC as the mass backbone and a catalog SFR as the activity axis. The audit records any overrides/tie-breakers used (e.g., gas fraction vetoes, infrared-trap notes), and the final \(t_50\) is taken from that documented decision path.

        ### Method 14: SPARC L[3.6] + z0MGS catalog SFR → sSFR → t₅₀ [Gas fraction override]
- Data source(s): SPARC Table 1 + z0MGS (IRSA): log SFR = −1.357 ± 0.042 (WISE-only)
- SPS model: Extended SFH; M★ via M/L ≈ 0.5 at 3.6μm
- Method validated: Single use
- Known issues: WISE-only SFR = floor in dust-poor LSBs ("Infrared Trap"); boundary case resolved by gas fraction override (35% rules out T3)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| F571-8 | 8.0 Gyr | T2 (Inefficient LSB) | MEDIUM | No | — |

---

### Method entries (verbatim from the audit)

### Method 7: LITTLE THINGS b₁ birthrate parameter → τ-model → t₅₀
- Data source(s): Zhang et al. (2012) Table 3: b₁ = 0.79 ± 0.17, b₀.₁ = 1.43 ± 0.29, M★ = 5.865 × 10⁷ M☉
- SPS model: τ-model fitted to b₁; T = 13.7 Gyr
- Metallicity: Hunter et al. (2012) Table 1: 12+log(O/H) = 8.3 ± 0.07 (method unverified)
- Explicit math: b₁ = 0.79 → τ ≈ 28.1 Gyr → t₅₀ ≈ 7.68 Gyr (verified by Deep Think)
- Method validated: Single use
- Known issues: τ-model simplified; no Tier-0 CMD; metallicity method not verified
- Correction: Original T1/3.3 Gyr ("blue = young" fallacy) revised to T2/7.7 Gyr

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| DDO 168 | 7.7 (+0.7/−0.7) Gyr | T2 (borderline T3) | MEDIUM | No | — |

---

### Data Gap: DDO 170
- No global M★, no global SFR, no sSFR, no integrated B-V color found. No valid age derivation completed.

---

### Method 9: SPARC mass + compiled SFR → b-parameter → session anchor interpolation → t₅₀
- Data source(s): SPARC (Lelli+2016): L[3.6], M_HI. Compiled SFR bracket ~1.0–1.4 M☉/yr (not primary table extraction)
- SPS model: Not directly used — t₅₀ interpolated from session anchors
- Method validated: Single use
- Known issues: SFR from compiled source; t₅₀ interpolated from session anchors (correlated bias risk)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| ESO 079-G014 | 7.2 Gyr | T2 | MEDIUM | No | — |

---

### Data Gap: ESO116-G012
- No global SFR, no integrated color. Gas fraction interpolation rejected as methodological precedent risk.

---

### Data Gap: ESO563-G021
- IRAS fluxes not retrieved, no catalog SFR. "IRAS detection → T2" inference rejected — detection is not a measurement.

---

### Method 18: SPARC mass + Hα activity confirmation → age by analog
- Data source(s): SPARC Table 1 + Spekkens & Giovanelli (2006) Hα+HI rotation curve
- SPS model: Age by analogy to NGC 0891 (physical twin)
- Method validated: Single use
- Known issues: No quantitative SFR; edge-on; age adopted by analog

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| IC4202 | 7.8 Gyr | T2 (Edge-On Massive Spiral) | MEDIUM | No | — |

---

### Method 20: SFR catalog presence (K&K 2013) + SPARC mass → unverified b → session anchor interpolation → t₅₀
- Data source(s): SPARC Table 1 + Karachentsev & Kaisina (2013) Table 3 (entry confirmed, exact SFR not extracted). Dicaire et al. (2008) Hα confirms activity.
- Method validated: Single use
- Known issues: Exact SFR not extracted; b ≈ 0.7 rests on unverified claim
- Tier secured by: SFR catalog presence = definitionally Active

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 0024 | 6.5 Gyr | T2 (Steady State) | MEDIUM | No | — |

---

### Method 22: SPARC mass + dual SFR (resolved stars + LVL catalog) → sSFR → b → anchor interpolation → t₅₀
- Data source(s): SPARC Table 1 + Davidge (2021) SFR ≈ 0.1 M☉/yr (resolved MSTO) + Lehmer et al. (2010) LVL: log SFR = −1.01
- Method validated: Single use
- Known issues: Age interpolated from anchors
- Strengths: Dual SFR confirmation; TRGB distance; "Sculptor Void" explains high gas (48%) + low activity (b ≈ 0.38)

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 0247 | 7.4 Gyr | T2 (Winding Down) | HIGH | No | — |

---

### Method 27: HALOGAS SFR + SPARC/WISE M★ → sSFR → b → Main Sequence anchor → t₅₀
- Data source(s): HALOGAS: SFR = 0.34 M☉/yr. SPARC/WISE: log M★ = 9.606 ± 0.012
- Derived: log sSFR = −10.07, b ≈ 1.2
- Morphology: SAcd (negligible bulge) → age tracks disk
- Method validated: Single use
- Known issues: t₅₀ = 6.0 is standard anchor, not galaxy-specific; constant SFH assumption

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 1003 | 6.0 Gyr | T2 (Steady State) | MEDIUM | No | — |

---

### Method 34: SPARC mass + SFRS TIR (Ashby+2011) → sSFR → b → super-spiral anchor → t₅₀
- Data source(s): SPARC Table 1: L[3.6] = 319.4×10⁹ L☉. SFRS (Ashby et al. 2011): S_60μm = 1.61 Jy, S_100μm = 4.51 Jy, log L_TIR = 10.99. RC3 (NED): (B-V)_T = 0.69 ± 0.01.
- Derived: log M★ ≈ 11.20, SFR ≈ 10.5 M☉/yr (Kroupa), b ≈ 0.91
- Method validated: Single use
- Known issues: TIR required for massive dusty spirals (FIR-only misses warm dust). Deep Think column misread corrected — 2.80 and 3.63 were F100/F60 ratio and Ks-F60 color, not fluxes.
- Color paradox resolved: NGC 2955 bluer (0.69) than NGC 5371 (0.75) because more active (b = 0.91 vs 0.8).

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 2955 | 7.4 Gyr | T2 (Vigorous Super-Spiral) | HIGH | No | — |

---

### Method 43: B-V ("Nonlocal Gravity" galaxy table) + SPARC gas fraction + Hα (VS01) + Downsizing calibration → t₅₀
- Data source(s): "Nonlocal Gravity" paper galaxy table: B-V = 0.60, M_disk = 1.13×10¹⁰ M☉, M_gas = 1.7×10⁹ M☉. SPARC: L[3.6] = 21.966×10⁹ L☉, M_HI = 1.888×10⁹ M☉. VS01 (Verheijen & Sancisi 2001): Hα + HI rotation curve.
- Method: Gas fraction ~15% for Scd → "Winding Down" phase. Downsizing calibration vs NGC 5033 (same phase, higher mass): 7.4 - 0.4 = 7.0 Gyr.
- Method validated: Single use
- Known issues: Inclination 79° (moderate dust). Ursa Major cluster — "Strangulation" likely accelerated gas depletion.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3917 | 7.0 Gyr | T2 (Winding Down) | HIGH | No | — |

---

### Method 46: SPARC mass + Tucker+2024 SFR → b + VS01 Hα + morphology calibration (Cluster Twin) → t₅₀
- Data source(s): SPARC: L[3.6] = 14.353×10⁹ L☉, M_HI = 1.214×10⁹ M☉, log M★ ≈ 9.86. Tucker et al. (2024) Table 1: log SFR = -0.68 ± 0.20 → SFR = 0.21 M☉/yr. VS01: Hα + HI.
- Method: b ≈ 0.40 ("Winding Down"). Cluster Twin with NGC 3917 (same environment, same gas depletion). Morphology adjustment: Sbc (+0.3 Gyr) vs NGC 3917's Scd → 7.3 Gyr.
- Method validated: Single use
- Known issues: Tucker mass (log 9.71) differs from SPARC (log 9.86) — SPARC adopted for consistency. Ursa Major cluster strangulation.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 3972 | 7.3 Gyr | T2 (Winding Down, Cluster Twin) | HIGH | No | — |

---

### Method 66: SPARC L[3.6] → M★ + IRAS FIR → SFR → sSFR + cold dust diagnostic → t₅₀ [Fading Super-Spiral]
- Data source(s): SPARC (Lelli+2016): L[3.6] = 391.076×10⁹ L☉, M_HI = 20.907×10⁹ M☉ → M★ ≈ 1.96×10¹¹ M☉ (log M★ ≈ 11.29). NED: IRAS S_60 = 0.3069 Jy, S_100 = 1.486 Jy; B-V = 0.72; Hα+[NII] detected. Derived: SFR ≈ 2.5–4.5 M☉/yr, log sSFR ≈ −10.9 ± 0.2, b ≈ 0.18–0.32.
- Method: **Cold dust diagnostic** S_100/S_60 = 4.8 (very cold, T < 20K) confirms minimal dust heating — "Honest Red" color from old stars, not "Dusty Disguise" (contrast NGC 5371 with S_100/S_60 ≈ 2–3). Morphological quenching: massive bulge stabilizes gas disk against collapse. Hα detection keeps galaxy in T2. "Fading Super-Spiral" — burning fuel very inefficiently despite 10.7% gas fraction.
- Method validated: Single use
- Known issues: IRAS SFR range (2.5–4.5) gives b range (0.18–0.32). Near T2/T3 boundary but Hα confirms active. "Honest Red vs Dusty Disguise" principle: same B-V can arise from different causes — cold dust ratio distinguishes.

| Galaxy | t₅₀ | Tier | Confidence | Flagged? | Flag Reason |
|--------|------|------|------------|----------|-------------|
| NGC 6195 | 7.9 Gyr | T2 (Fading Super-Spiral) | HIGH | No | — |

---

