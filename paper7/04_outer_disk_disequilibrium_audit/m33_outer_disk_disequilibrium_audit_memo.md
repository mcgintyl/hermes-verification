# M33 Outer-Disk Disequilibrium Audit Memo

**Lane:** M33-Corbelli source-native thick-disk reconstruction  
**Purpose:** Determine whether the locked Hermes outer-tail residual zone overlaps documented regions where M33's outer disk may violate, or at minimum complicate, the circular-equilibrium assumption.  
**Status:** Audit only. No Hermes rerun. No MOND rerun. No board edits. No refits.

## Executive finding

The spatial association is strong. The locked Hermes stress zone begins beyond roughly **10 kpc**, while the published M33 literature places the major outer-disk complications in the same radial domain: the H I warp begins near **8 kpc**, the stellar surface-brightness/age-gradient break is near **8–9 kpc**, Barker's far-outer resolved-CMD fields sit at **9.1 and 11.6 kpc**, and the extended disturbed H I features traced by Putman et al. reach **22 kpc**. The Corbelli board itself extends to **22.72 kpc**, so the Hermes penalty is almost entirely located in the part of the board most strongly flagged by the literature as warped, low-density, structurally mixed, and possibly accretion/tidal-influenced.

This does **not** erase the locked MOND result. MOND remains clean on the same frozen board. It also does **not** prove that disequilibrium causes the Hermes residual. The stronger, defensible conclusion is:

> The M33 Hermes residual is not distributed through the main disk. It is spatially concentrated in the literature-defined outer-disk transition/warp/break regime. M33 should therefore be archived as a **Hermes outer-tail stress case with a disequilibrium / multi-origin outer-disk flag**, not as a generic full-disk failure.

## Frozen inputs used

This audit uses only the already-frozen prediction files:

- Hermes primary age lane: $t_{50}=6.0$ Gyr.
- Hermes conservative age lane: $t_{50}=7.1$ Gyr.
- MOND fixed simple comparator.
- Same frozen board: source-native SPS/pixel-SED stellar mass, helium factor 1.33, H$_2$ prescription, annular thick-disk kernel, and $V_{m bul}=0$.

## Literature feature map

| feature                                                                         | radius_kpc                                                          | overlaps_locked_outer_tail_R_gt10                                                                     | source                                 |
|:--------------------------------------------------------------------------------|:--------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------|:---------------------------------------|
| Corbelli 2014 HI warp / PA twist starts                                         | R > 8                                                               | Yes; starts 2 kpc before Hermes outer-tail boundary                                                   | Corbelli et al. 2014; Koch et al. 2018 |
| M33 rotation curve becomes more perturbed and flatter                           | after ~4 kpc; out to 23 kpc                                         | Yes; outer tail lies entirely in perturbed/flatter HI regime                                          | Corbelli et al. 2014                   |
| Far outer disk break / age-gradient reversal                                    | S1 9.1 kpc, S2 11.6 kpc; break near ~9 kpc                          | Yes; most Hermes tail lies outside S1 and many points outside S2                                      | Barker et al. 2011                     |
| Extended warps, arc, diffuse gas, southern cloud/filament                       | features to 22 kpc                                                  | Yes; same 10–23 kpc radial domain                                                                     | Putman et al. 2009                     |
| Tidal radius from M31 orbit analysis                                            | <15 kpc in majority of possible orbits                              | Partial; 15–23 kpc especially exposed                                                                 | Putman et al. 2009                     |
| Kam et al. perturbed outer gas, anomalous low-brightness HI, asymmetric motions | outer gas distribution; HI disk to 23 kpc                           | Yes; supports non-circular/asymmetric outer-zone caution                                              | Kam et al. 2017                        |
| 2024 Corbelli/Burkert accretion scenario for outer misalignment                 | outer disk / HI clouds, not single radius; outer disk PA twist ~30° | Yes; interpretive support that warp may be gas-accretion/cosmic filament rather than past M31 passage | Corbelli & Burkert 2024                |

## Radial mask metrics

The table below uses the frozen per-point residuals and no refitting. Residuals are model minus observed velocity, so negative values mean underprediction.

| mask                               |   N | model                       |   chi2nu_mask |   chi2_fraction_of_model_total_pct |   mean_residual_model_minus_obs_kms |
|:-----------------------------------|----:|:----------------------------|--------------:|-----------------------------------:|------------------------------------:|
| warp_zone_R_gt_8                   |  31 | Hermes_primary_t50_6p0      |         5.292 |                             90.493 |                             -48.487 |
| warp_zone_R_gt_8                   |  31 | Hermes_conservative_t50_7p1 |         5.594 |                             90.726 |                             -49.923 |
| warp_zone_R_gt_8                   |  31 | MOND_fixed_simple           |         0.202 |                             90.693 |                              -7.843 |
| main_disk_R_le_10_locked_boundary  |  31 | Hermes_primary_t50_6p0      |         0.802 |                             13.711 |                              -7.100 |
| main_disk_R_le_10_locked_boundary  |  31 | Hermes_conservative_t50_7p1 |         0.845 |                             13.704 |                              -8.632 |
| main_disk_R_le_10_locked_boundary  |  31 | MOND_fixed_simple           |         0.024 |                             10.660 |                               0.465 |
| outer_tail_R_gt_10_locked_boundary |  27 | Hermes_primary_t50_6p0      |         5.794 |                             86.289 |                             -51.622 |
| outer_tail_R_gt_10_locked_boundary |  27 | Hermes_conservative_t50_7p1 |         6.109 |                             86.296 |                             -53.049 |
| outer_tail_R_gt_10_locked_boundary |  27 | MOND_fixed_simple           |         0.228 |                             89.340 |                              -9.394 |
| Putman_tidal_radius_R_ge_15        |  16 | Hermes_primary_t50_6p0      |         7.061 |                             62.318 |                             -60.353 |
| Putman_tidal_radius_R_ge_15        |  16 | Hermes_conservative_t50_7p1 |         7.362 |                             61.628 |                             -61.608 |
| Putman_tidal_radius_R_ge_15        |  16 | MOND_fixed_simple           |         0.322 |                             74.825 |                             -11.958 |

## Primary comparison: main disk vs outer disk

- At **$R\leq8$ kpc**, approximately the pre-warp main disk, Hermes primary has $\chi^2_
u=0.638$ and mean residual $-4.1$ km/s.
- At **$R>8$ kpc**, the warp-zone mask, Hermes primary jumps to $\chi^2_
u=5.292$ and mean residual $-48.5$ km/s; this zone contains **90.5%** of Hermes' total $\chi^2$.
- At the locked Hermes boundary **$R\leq10$ kpc**, the primary lane has $\chi^2_
u=0.802$.
- At **$R>10$ kpc**, it has $\chi^2_
u=5.794$ and mean residual $-51.6$ km/s; this is the already-locked **86.3%** outer-tail $\chi^2$ contribution.
- In the farthest Putman/tidal-radius diagnostic mask **$R\geq15$ kpc**, Hermes primary has $\chi^2_
u=7.061$ and mean residual $-60.4$ km/s.

MOND does not show the same amplitude of stress: in the locked outer zone $R>10$ kpc, fixed simple MOND has $\chi^2_
u=0.228$ and mean residual $-9.4$ km/s.

## Causal-chain assessment

### What is strongly supported

1. **The residual is spatially localized.** Hermes fits the main disk non-catastrophically and very nearly cleanly under the $R\leq10$ kpc mask, then fails primarily in the outer tail.
2. **The residual overlaps the warp/break regime.** The literature places the H I warp onset near $R\simeq8$ kpc, the stellar break near $R\simeq8$–9 kpc, and the far-outer stellar population transition at 9.1–11.6 kpc.
3. **The affected region is physically distinct.** Published work describes the M33 outskirts as warped, PA-twisted, structurally mixed, low surface density, and possibly shaped by gas accretion and/or older tidal interaction.

### What is plausible but not proven

1. **Circular-equilibrium contamination.** The H I tilted-ring target is still an observed rotation target, but the outer disk is warped and asymmetric. That makes circular-speed interpretation more fragile than in the inner disk.
2. **Multi-origin outer material.** The coexistence of young sparse stars, older metal-poor outer stars, H I clouds, diffuse gas, and possible accretion signatures makes the outer board less likely to be a single homogeneous equilibrium disk.
3. **Hermes is flagging the outer component.** The Hermes ceiling/gate structure refuses to supply the gravity needed by the outer velocities; this can be read as a stress diagnostic. But it remains a model-specific stress, not proof by itself that the velocities are non-equilibrium.

### What is not established

1. The audit does **not** prove that the M33 outer rotation curve is wrong.
2. It does **not** prove that MOND is merely overfitting tidal debris.
3. It does **not** prove that removing or downweighting the outer disk is justified for the locked board.
4. It does **not** rescue Hermes on the full frozen board.

## Paper-safe interpretation

Recommended language:

> The first external Format B board, M33-Corbelli, shows a localized Hermes stress rather than a distributed full-disk failure. Hermes fits the main disk at non-catastrophic to clean accuracy, but underpredicts the outer $R>10$ kpc disk by roughly 52–53 km/s. This residual zone overlaps the published onset of the H I warp, the stellar disk break, the Barker far-outer population fields, and the extended disturbed H I structures around M33. The association does not invalidate the MOND comparator result, which remains clean on the same frozen board, but it flags the M33 outer disk as a physically distinct regime where circular-equilibrium and single-component disk assumptions are especially load-bearing.

## Recommended archive label

**M33-Corbelli source-native thick-disk reconstruction — Hermes pass/elevated; MOND clean; Hermes outer-tail stress spatially associated with literature-defined outer-disk warp/break/accretion regime.**

## Files included

- `m33_outer_disk_mask_metrics.csv`
- `m33_outer_disk_literature_feature_map.csv`
- `m33_residuals_disequilibrium_zone_overlay.png`
- `m33_chi2_disequilibrium_zone_overlay.png`

