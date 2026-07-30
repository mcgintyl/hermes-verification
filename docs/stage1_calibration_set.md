# Stage-1 calibration set

In Stage 1 of the derivation, the global β parameter of the Hermes equation was
fitted against **seven** young, bulgeless galaxies, listed in
[`stage1_calibration_set.csv`](stage1_calibration_set.csv).

## The Stage-1 β fit is age-free

β is the single global parameter in a purely *geometric* boost,

    g_model(R) = g_bar(R) · [1 + β · φ(R)]

where φ(R) is built from the shear and curvature of the baryonic rotation curve.
β was fit by minimizing the unweighted mean per-galaxy χ²/N (V_obs vs V_model)
across the seven galaxies, using **only their SPARC rotation curves**. Stellar
age (t₅₀) is **not** an input to this fit (best fit β ≈ 2.807, locked at 2.81).

t₅₀ enters the derivation in two other places — never the Stage-1 β fit: it was
used to *select* a young, bulgeless sample, and its mean sets the reference age
t_ref = 3.49 Gyr used from Stage 2 onward (the age-dependent law
β_eff = 2.81·exp[(t_ref − t₅₀)/τ]).

## Historical working ages and the audit

The seven were selected on pre-audit working ages (mean = 3.49 Gyr = t_ref).
Both those values and the current post-audit ages are in
[`stage1_calibration_set.csv`](stage1_calibration_set.csv). The later full-sample
age audit revised several and dropped three from the final 133-galaxy sample
(`paper1/ages_133.csv`):

| Galaxy | Pre-audit t₅₀ | Why dropped |
|---|---|---|
| ESO 444-G084 | 3.6 | Post-audit age (≈2.2 Gyr) landed on the T1/T2 boundary; excluded. |
| UGC 06930 | 4.0 | Provisional age traced to a 21″ central-aperture SFR (Watson et al. 2012) scaled to the whole disk — overstates global SFR by ~2 dex. Rejected ("aperture fallacy"). |
| F583-4 | 4.0 | Provisional only; no reliable age could be established on audit (no clean multi-band color, UV, or SFR). |

The four that survived — D631-7, DDO 154, DDO 168, UGC 07608 — carry current
t₅₀ and g₉₈ in the published files. Reproducing the Stage-1 β fit itself needs
only the seven SPARC rotation curves (it is age-free); t₅₀ and g₉₈ are needed for
the age-dependent stages (2+) and the final validation.

## Three different "calibration" groupings — do not conflate

These are distinct sets that overlap only slightly:

| Grouping | Size | Where | Members |
|---|---|---|---|
| **Stage-1 calibration set** | 7 (4 surviving) | Derivation Summary | D631-7, DDO 154, DDO 168, ESO 444-G084, F583-4, UGC 06930, UGC 07608 |
| **Tier 1** (data-quality class) | 3 | Paper 1 §5, `Table_S1_Age_Methods.csv` | D631-7, NGC 3521, UGC 00731 |
| **Pristine-22** | 22 | `Table_S3_Pristine22_Calibration.csv` | a separate clean-photometry subset |

Only **D631-7** is common to the Stage-1 set and the Tier-1 data class. "Tier 1"
in Paper 1 is a *data-reliability* label, **not** the calibration set.

Note also that `in_calibration_sample` in
[`galaxy_age_method_index.csv`](galaxy_age_method_index.csv) marks membership in
the **final 133-galaxy validation sample** (what Config G was tested on), which
is a different thing again from the Stage-1 calibration set above.
