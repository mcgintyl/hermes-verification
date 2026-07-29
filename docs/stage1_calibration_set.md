# Stage-1 calibration set

The global β parameter of the Hermes equation was first fitted, in Stage 1 of
the derivation, against **seven** young, bulgeless galaxies. Those seven are
listed in [`stage1_calibration_set.csv`](stage1_calibration_set.csv).

**Three of the seven did not survive the later full-sample age audit** and are
therefore not in the final 133-galaxy sample (`paper1/ages_133.csv`):

| Galaxy | Why excluded |
|---|---|
| ESO 444-G084 | Age derived (t₅₀ ≈ 2.2 Gyr) but on the T1/T2 boundary; dropped from the final sample. |
| UGC 06930 | Age traced to an SFR measured in a 21″ central aperture (Watson et al. 2012) scaled to the whole disk — overstates global SFR by ~2 dex. Rejected ("aperture fallacy"). |
| F583-4 | No usable age data (no clean multi-band color, UV, or SFR). No t₅₀ was derived. |

**The four that survived — D631-7, DDO 154, DDO 168, UGC 07608 — are the
calibration galaxies present in the published files**, with `t50_gyr` and `g98`
for each. They are the ones available to reproduce or cross-validate the β fit.

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
