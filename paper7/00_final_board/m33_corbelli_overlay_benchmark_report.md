# M33-Corbelli Published-Figure Overlay Benchmark

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Status:** Format B reconstructed engineering board; overlay benchmark passed for the source-native/top-panel component curves; final physics execution still not authorized.

## Doctrine lock

M33 is the first real Format B candidate board, but its component velocities are reconstructed products. Therefore the disk-potential kernel is part of the data provenance and must be frozen before any physics comparison.

## Benchmark target

The overlay benchmark uses Corbelli et al. 2014 Figure 12. Figure 12 shows the M33 rotation curve and NFW dynamical models; the paper identifies the red small-dashed component as gas and the blue large-dashed component as the stellar disk. The top panel is the relevant source-native provenance target because it shows the most likely dynamical models when the rotation fit is combined with the stellar-mass-map likelihood and the cosmological C--Mh relation. The bottom panel shows dynamical best-fit branches, including heavier stellar-disk variants, and is therefore not the default no-tune board target.

## Kernel used

The current kernel is the provisional annular quadrature thick-disk kernel from the engineering reconstruction:

- gas half-thickness: 0.5 kpc;
- stellar half-thickness: flaring from 0.1 kpc at the center to 1.0 kpc at 23 kpc;
- gas profile: 1.33 × [HI + H2];
- H2 profile: 10 exp(-R/2.2) Msun pc^-2;
- stellar profile: source-native SPS / pixel-SED stellar mass surface density;
- Vbul = 0.

This is a source-native thick-disk reconstruction, not a SPARC-equivalent unit-M/L reconstruction.

## Overlay result

The reconstructed gas and stellar curves overlay the **top-panel** Corbelli Figure 12 components within plot-line tolerance. Quantitative pixel-extraction sanity checks give:

| panel   |   gas_N |   gas_median_abs_delta_kms |   gas_p90_abs_delta_kms |   stellar_N |   stellar_median_abs_delta_to_blue_median_kms |   stellar_p90_abs_delta_to_blue_median_kms |   stellar_median_abs_delta_to_blue_envelope_kms |   stellar_p90_abs_delta_to_blue_envelope_kms |   stellar_fraction_inside_blue_envelope |
|:--------|--------:|---------------------------:|------------------------:|------------:|----------------------------------------------:|-------------------------------------------:|------------------------------------------------:|---------------------------------------------:|----------------------------------------:|
| top     |      56 |                   0.866085 |                 1.58075 |          56 |                                      0.228763 |                                    1.90313 |                                         0       |                                    0.0133351 |                                0.892857 |
| bottom  |      56 |                   0.740274 |                 1.45157 |          56 |                                      5.96165  |                                    8.63348 |                                         4.82278 |                                    8.14233   |                                0        |

Interpretation:

- Gas: passes in both panels, as expected; the red gas component is effectively the same branch in the figure.
- Stellar disk: passes in the top panel. The bottom-panel offset is expected because the bottom panel includes heavier dynamical best-fit stellar branches and is not the fixed source-native board target.

## Decision

The overlay sanity check passes for the intended source-native component curves. The annular quadrature kernel can be **provisionally frozen for the M33-Corbelli lane**, subject to PI review.

No Hermes run, MOND run, age anchoring, or gravity comparison has been performed.

## Files

- `corbelli_fig12_overlay_reconstruction.png`: visual overlay benchmark.
- `corbelli_fig12_overlay_digitization_summary.csv`: quantitative pixel-extraction summary.
- `source_fig12_top_red_gas_digitized_binned.csv`, `source_fig12_top_blue_stellar_digitized_binned.csv`: digitized source component sanity curves.
- `m33_fig12_overlay_benchmark.py`: overlay/digitization script.
