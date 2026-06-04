# M33-Corbelli Final Board-Freeze Package — Bracketed Age Update

**Lane:** Corbelli/M33  
**Galaxy:** M33 / NGC 598 / UGC 1117  
**Locked label:** **M33-Corbelli source-native thick-disk reconstruction**  
**Status:** **Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age status Amber with bracketed $t_{50}$ locked for sensitivity use; no Hermes/MOND execution authorized.**

## No-result warning

No Hermes run was performed. No MOND run was performed. No gravity comparison exists for this lane. This package freezes the reconstructed baryonic board and records the bracketed Amber age decision only.

## Doctrine lock

**M33 is the first real Format B candidate board, but its component velocities are reconstructed products. Therefore the disk-potential kernel is part of the data provenance and must be frozen before any physics comparison.**

## Package contents

| File | Purpose |
|---|---|
| `m33_corbelli_reconstructed_board.csv` | Reconstructed board fields: $R$, $V_{\rm obs}$, $errV$, $V_{\rm gas}$, $V_{\rm disk}$, $V_{\rm bulge}$, $V_{\rm bar}$ |
| `m33_corbelli2014_table1_import.csv` | Imported source table used to build target and surface-density inputs |
| `m33_corbelli_mass_closure_table.csv` | Stellar, H$_2$, H I, gas, and baryonic mass closure summary |
| `m33_corbelli_kernel_specification.md` | Provisionally frozen annular thick-disk kernel specification |
| `m33_corbelli_kernel_provenance_review.md` | Kernel review and comparison against Corbelli/Casertano and SPARC conventions |
| `m33_corbelli_overlay_benchmark_report.md` | Published-figure overlay benchmark report |
| `corbelli_fig12_overlay_reconstruction.png` | Visual overlay against Corbelli Figure 12 component curves |
| `m33_corbelli_reconstruction_kernel.py` | Reproducible reconstruction script |
| `m33_corbelli_remaining_limitations.md` | Limitations and warnings before physics use |
| `m33_corbelli_age_anchor_audit.md` | M33-only age-source audit, pre-lock |
| `m33_age_conversion_memo_prelock.md` | Age-conversion memo that produced the two candidate values |
| `m33_age_bracketed_amber_lock_addendum.md` | Bracketed Amber age lock decision |
| `m33_age_bracketed_amber_lock.csv` | Machine-readable age bracket table |
| `m33_corbelli_board_summary.json` | Compact board and age-status metadata |

## Source table provenance

The engineering import contains 58 radial rows with columns:

`R_kpc`, `V_r_kms`, `sigma_V_kms`, `Sigma_HI_Msun_pc2`, `Sigma_star_Msun_pc2`

Board target mapping:

| Board field | Source field |
|---|---|
| $R$ | `R_kpc` |
| $V_{\rm obs}$ | `V_r_kms` |
| $errV$ | `sigma_V_kms` |
| $\Sigma_{\rm HI}$ | `Sigma_HI_Msun_pc2` |
| $\Sigma_\star$ | `Sigma_star_Msun_pc2` |

## Reconstructed board schema

The reconstructed board has 58 radial points spanning 0.24–22.72 kpc.

Output columns:

`galaxy`, `R_kpc`, `Vobs_kms`, `errV_kms`, `Sigma_HI_Msun_pc2`, `Sigma_H2_Msun_pc2`, `Sigma_gas_total_He_Msun_pc2`, `Sigma_star_Msun_pc2`, `Vgas_kms`, `Vdisk_kms`, `Vbulge_kms`, `Vbar_kms`, `kernel`, `source_ml_convention`

## Gas profile

\[
\Sigma_{\rm H_2}(R)=10\exp(-R/2.2)\ M_\odot\ {\rm pc}^{-2}
\]

\[
\Sigma_{\rm gas}(R)=1.33\,[\Sigma_{\rm HI}(R)+\Sigma_{\rm H_2}(R)]
\]

The helium factor is **1.33**. The gas half-thickness is **0.5 kpc**.

## Stellar profile and M/L convention

The stellar profile is Corbelli's **BVIgi / pixel-SED stellar mass surface-density profile**. This is a source-native SPS / pixel-SED stellar-mass lane with radially varying M/L provenance.

It is **not** SPARC-equivalent unit-M/L photometry. The M/L convention is provenance, not a tunable Hermes parameter.

## Bulge convention

\[
V_{\rm bulge}(R)=0
\]

at every radius. This follows the source-native no-genuine-bulge/no-prominent-bar interpretation used for the board.

## Kernel specification

The board uses the provisionally frozen **annular quadrature thick-disk kernel**:

- radial integration step: 0.01 kpc;
- angular samples: 1440 per annulus;
- gas half-thickness: 0.5 kpc;
- stellar half-thickness: linear flare from 0.1 kpc to 1.0 kpc over the disk;
- same Newtonian annular force kernel applied to gas and stars;
- no M/L tuning and no fit to the observed rotation curve.

The kernel is accepted provisionally for this lane because the published-figure overlay benchmark passes against the intended Corbelli Figure 12 top-panel component target.

## Mass closure

| Quantity | Integrated value | Source reference | Ratio | Status |
|---|---:|---:|---:|---|
| Stellar mass | $4.890\times10^9\ M_\odot$ | $4.9\times10^9\ M_\odot$ | 0.998 | Pass |
| H$_2$ mass | $3.040\times10^8\ M_\odot$ | $3.0\times10^8\ M_\odot$ | 1.013 | Pass |
| H I profile mass | $1.710\times10^9\ M_\odot$ | true H I $\approx1.53\times10^9\ M_\odot$ | 1.118 | Accepted under source outer-ring filling warning |
| Total gas with He | $2.678\times10^9\ M_\odot$ | derived | — | Derived |
| Total baryonic mass | $7.569\times10^9\ M_\odot$ | derived | — | Derived |

## Overlay benchmark

The reconstructed component curves passed the source-native Figure 12 top-panel overlay benchmark:

| Component | Median offset | 90th-percentile offset | Verdict |
|---|---:|---:|---|
| Gas | 0.87 km/s | 1.58 km/s | Pass |
| Stellar disk | 0.23 km/s | 1.90 km/s | Pass |

The bottom-panel stellar offset is not a failure because the bottom panel represents heavier dynamical branches, while this board preserves the source-native SPS / pixel-SED stellar mass profile.

## Age status — Amber bracketed lock

No direct global published $t_{50}$ for M33 was identified. The primary age evidence is resolved CMD/SFH, so the source class is A-method; the confidence is Amber because the board-level age is converted from spatial fields rather than directly tabulated as a global half-mass time.

The age treatment is locked as a bracketed sensitivity pair:

| Age lane | $t_{50}$ | Role | Required future use |
|---|---:|---|---|
| **Primary / source-native mass-weighted** | $6.0\pm1.2$ Gyr | Best match to board-level stellar half-mass definition using the Corbelli stellar-mass weighting kernel | Must be reported in eventual Hermes result |
| **Conservative older-weighted bracket** | $7.1\pm1.3$ Gyr | Older sensitivity bracket preserving inner-disk and Barker S2 / outer-field age leverage | Must be reported in eventual Hermes result |

Both values must be carried into the eventual Hermes run. The model result may not select only whichever age performs better. MOND does not use $t_{50}$ internally, but both age lanes must be shown in the comparison record so the Hermes age-sensitivity envelope is explicit.

## Remaining limitations

- The board is source-native, not SPARC-equivalent.
- Component velocities are reconstructed, not directly published as rotmod-style arrays.
- The H I radial-profile mass excess is accepted but remains documented.
- The board is 1D/axisymmetric and does not implement a full 3D warp model.
- No direct global published M33 $t_{50}$ exists in the current audit.
- The age lock is Amber and bracketed, not Green and single-valued.
- No model result exists.

## Final status

**Corbelli/M33 — Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age status Amber with bracketed $t_{50}$ locked for sensitivity use; no Hermes/MOND execution authorized.**
