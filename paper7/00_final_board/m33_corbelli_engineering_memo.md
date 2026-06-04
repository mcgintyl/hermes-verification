# M33 / Corbelli Format B Reconstruction — Engineering Pass Only

**Status:** Format B candidate — source-reconstructable, not board-run-ready.  
**Action performed:** board-construction engineering only.  
**Action not performed:** no Hermes run, no MOND run, no age anchoring, no pilot board.

## Source table import

Imported 58 radial rows with columns:

`R_kpc`, `V_r_kms`, `sigma_V_kms`, `Sigma_HI_Msun_pc2`, `Sigma_star_Msun_pc2`

The target fields are mapped as:

- `R = R_kpc`
- `Vobs = V_r_kms`
- `errV = sigma_V_kms`

## Source-native profile construction

Gas profile:

`Sigma_H2(R) = 10 exp(-R / 2.2) Msun pc^-2`

`Sigma_gas,total(R) = 1.33 [Sigma_HI(R) + Sigma_H2(R)]`

Stellar profile:

`Sigma_star(R)` is the Corbelli BVIgi / pixel-SED stellar mass surface density from Table 1. This is not a SPARC-equivalent unit-M/L photometric component. It is a source-native SPS / pixel-SED stellar-mass lane.

Bulge:

`Vbulge = 0`, following the Corbelli no-genuine-bulge interpretation.

## Frozen disk-potential kernel

The reconstruction uses a uniform-thickness annular quadrature kernel in the disk plane.

- radial integration step: 0.01 kpc
- angular integration: 1440 azimuth samples per annulus
- gas half-thickness: 0.5 kpc
- stellar half-thickness: linearly flaring from 0.1 kpc at the center to 1.0 kpc at 23 kpc
- same Newtonian annular force kernel applied to gas and stars

This was selected because it is transparent, reproducible, and follows Corbelli's stated geometry choices closely enough for engineering reconstruction. It is not used as a fit lever.

## Closure checks

| Quantity | Integrated value | Source reference | Ratio |
|---|---:|---:|---:|
| Stellar mass | 4.890e9 Msun | 4.9e9 Msun | 0.998 |
| H2 mass | 3.040e8 Msun | 3.0e8 Msun | 1.013 |
| HI profile mass | 1.710e9 Msun | 1.53e9 Msun true HI | 1.118 |
| Total gas with He | 2.678e9 Msun | derived | — |
| Total baryonic mass | 7.569e9 Msun | derived | — |

The HI profile mass being above the true HI mass is expected: Corbelli states that integrating the radial profile out to 23 kpc gives a mass higher than the true HI mass because the outer HI emission does not fill the whole ring.

## Output columns

The engineering board CSV contains:

`galaxy`, `R_kpc`, `Vobs_kms`, `errV_kms`, `Sigma_HI_Msun_pc2`, `Sigma_H2_Msun_pc2`, `Sigma_gas_total_He_Msun_pc2`, `Sigma_star_Msun_pc2`, `Vgas_kms`, `Vdisk_kms`, `Vbulge_kms`, `Vbar_kms`, `kernel`, `source_ml_convention`

## Stop point

The reconstructed board passes first-order stellar and gas mass closure. It is still not authorized for Hermes or MOND until PI approval after reviewing the kernel provenance and any source-paper component-curve comparison.
