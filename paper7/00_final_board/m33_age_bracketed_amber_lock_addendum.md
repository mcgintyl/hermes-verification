# M33 Bracketed Amber Age Lock Addendum

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Board status:** Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete.  
**Age status:** **Amber, bracketed $t_{50}$ locked for sensitivity use.**  
**Execution status:** No Hermes run. No MOND run. No gravity comparison.

## Decision

The M33 age treatment is locked as a **bracketed Amber age anchor**, not a single Green age.

| Age lane | $t_{50}$ | Role | Use rule |
|---|---:|---|---|
| **Primary / source-native mass-weighted** | $6.0 \pm 1.2$ Gyr | Best current board-level stellar half-mass estimate using the Corbelli stellar-mass weighting kernel and Williams/Barker resolved-CMD evidence. | Must be carried into the eventual Hermes run. |
| **Conservative older-weighted bracket** | $7.1 \pm 1.3$ Gyr | Older sensitivity bracket preserving inner-disk age leverage and Barker S2 / far-outer-disk uncertainty. | Must be carried into the eventual Hermes run. |

No single global published $t_{50}$ for M33 was identified in the current audit. The two values above are transparent conversions from high-tier resolved-CMD/SFH evidence, with the difference between them treated as part of the uncertainty envelope rather than as a tuning choice.

## Rationale

The primary/source-native estimate, $t_{50}=6.0\pm1.2$ Gyr, best matches the board-level half-mass definition because it uses present-day stellar-mass weighting from the same Corbelli source-native stellar surface-density profile used in the reconstructed baryonic board.

The conservative estimate, $t_{50}=7.1\pm1.3$ Gyr, preserves older leverage from the inner disk and from the Barker S2 outer-field uncertainty. It is intentionally retained as a sensitivity bracket rather than discarded.

Digitizing or recovering numeric Williams/Barker cumulative SFH tables may refine the age estimate later, but that refinement is not treated as a prerequisite for the first execution if the PI accepts the bracketed Amber lock.

## Mandatory reporting rule

Any eventual model result for M33 must report both age lanes. The model output may not select only whichever age performs better.

Minimum future result table requirement:

| Galaxy | Age lane | $t_{50}$ used | Hermes result | MOND result | Notes |
|---|---|---:|---:|---:|---|
| M33 | Primary/source-native | $6.0\pm1.2$ Gyr | pending | pending | no run yet |
| M33 | Conservative older-weighted | $7.1\pm1.3$ Gyr | pending | pending | no run yet |

MOND does not use $t_{50}$ internally, but the same board and reporting table should carry both age lanes so that the Hermes sensitivity bracket is visible in the comparison record.

## Source classes

The age source class remains **A-method** because the primary evidence is resolved CMD/SFH. The confidence remains **Amber** because no direct global published $t_{50}$ exists and the board-level value is converted from spatial fields.

Primary evidence:

- Williams et al. 2009: resolved HST/ACS CMD/SFH fields across the M33 disk; strong inside-out growth evidence; field-based rather than global $t_{50}$ table.
- Barker et al. 2011: deep outer-disk resolved CMD/SFH fields; fills the disk-break/outer-disk gap but does not alone provide a board-wide global $t_{50}$.

Support-only sources remain PHATTER/recent-SFH, LPV/SFH, and parametric/SPS/chemical-evolution studies. These sources may support the bracket but do not replace the Williams/Barker resolved-CMD chain.

## Stop condition

This addendum updates the board-freeze package with the age decision only. It does not authorize model execution.

Still not authorized:

- no Hermes run;
- no MOND run;
- no gravity comparison;
- no interpretation of model performance.

## Updated status

**Corbelli/M33 — Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age status Amber with bracketed $t_{50}$ locked for sensitivity use; no Hermes/MOND execution authorized.**
