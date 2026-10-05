# Data sources

The two tables in `data/` are derived from two public Zenodo records, both released under the Creative Commons Attribution 4.0 license (CC BY 4.0). They are not unmodified copies: they hold a subset of the published columns, joined on PlateIFU, plus derived columns, with the changes listed below. `build_merged_table.py` rebuilds both tables from the public files and documents every step; `reproduce.py --from-public` runs that rebuild and compares it with the shipped tables (at most 2 ulps per value, identical NaN pattern).

## DynPop JAM catalog and circular-velocity curves (DynPop I and VII)

Zhu, K., Cappellari, M., Mao, S., Lu, S., Li, R., Shi, Y., Simon, D. A., Fu, Y., & Wang, X. (2026). *The dataset related to "MaNGA DynPop. VII. A Unified Bulge-Disk-Halo Model for Explaining Diversity in Circular Velocity Curves of 6000 Spiral and Early-type Galaxies"* (v2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17518315. License: CC BY 4.0.

Papers: Zhu et al. (2023), MaNGA DynPop I, *MNRAS*, 522, 6326 (the JAM catalog); Zhu et al. (2025), MaNGA DynPop VII, *ApJS*, 280, 55 (the circular-velocity curves).

Files used: `SDSSDR17_MaNGA_JAM_v2.fits` (MD5 `88f7dc25f6cadb2c5a2de4b190d00329`, 21,479,040 bytes) and `SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt` (MD5 `b7aa4a4f188196c39369f918768aa8c0`, 1,134,922 bytes).

## DynPop II stellar populations and star-formation histories

Lu, S., Zhu, K., Cappellari, M., Li, R., Mao, S., & Xu, D. (2025). *The updated catalogue of "MaNGA DynPop - II. Global stellar population, gradients, and star-formation histories from integral-field spectroscopy of 10K galaxies: link with galaxy rotation, shape, and total-density gradients"* (version 2) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15742825. License: CC BY 4.0.

Paper: Lu et al. (2023), MaNGA DynPop II, *MNRAS*, 526, 1022.

File used: `DynPop2_SP_SFH_v2.hdf5.zip` (MD5 `870368917b3a7c30f2da7598f21bdb2c`, 229,323,792 bytes), which holds `SP_SFH_v2.hdf5` (MD5 `543bc4d28dc5d33e1244a7c33bd631b6`, 2,988,377,528 bytes; the record also offers it unzipped).

## What this package changed

- **Column subset.** From the JAM catalog's global extension (HDU 1): `plateifu`, `mangaid`, `Qual`, `drp3qual`, `Lambda_Re`, `Sigma_Re`, `Eps_MGE`, `nsa_sersic_mass`, `nsa_sersic_n`, `nsa_sersic_ba`. From HDUs 2, 4 and 8 (cylindrical JAM with the mass-follows-light, NFW and gNFW mass models): the dynamical M/L, enclosed-mass, dark-matter-fraction, gNFW slope and reduced chi-squared columns, prefixed `mfl_cyl_`, `nfw_cyl_` and `gnfw_cyl_`. From the SP/SFH file: eleven scalars per galaxy (`T50`, `T90`, light- and mass-weighted age and metallicity within Re, intrinsic and observed stellar M/L, stellar mass, signal-to-noise and chi-squared), prefixed `sp_`. From the circular-velocity table: `Re`, `rmax`, `Vc_Re`, `Vc_rmax`, `Vcmax` and `Qual`, prefixed `cvc_`.
- **Join.** The three files list the same 10,296 PlateIFU entries in the same order; the build asserts that before joining.
- **Derived columns.** `DML_mfl_int` and `DML_mfl_obs` (JAM mass-follows-light log dynamical M/L minus the intrinsic or observed SPS log stellar M/L), `Dmass_mfl_Re`, `Dmass_gnfw_Re`, `gnfw_dm_to_stellar_log_Re`, `nfw_dm_to_stellar_log_Re`, `logRe` (log10 of `Re_arcsec_MGE`, an angular size in arcsec), `logSigma_Re`, `Vc_ratio_rmax_Re`, `Vc_ratio_max_Re` and `cvc_slope_rmax_Re`.
- **Cleaning.** The NSA sentinel -9999 in `nsa_sersic_n` and `nsa_sersic_ba` (53 rows) is set to NaN, and infinite values from logarithms or ratios are set to NaN. No row is dropped; the quality cut and finite-value requirements are applied by the analysis scripts.
- **Distance and subsample table.** `manga_dynpop_DA_target.csv` copies `plateifu`, `mangaid`, `DA` (adopted angular-diameter distance, Mpc) and `target` (MaNGA subsample: 0 Primary, 1 Secondary, 2 color-enhanced) unchanged from the JAM catalog's HDU 1.

## Checksums

| file | SHA-256 | bytes |
|---|---|---|
| manga_dynpop_merged_thin_firstpass.csv | 13336cd42d374514d895d51ca1334f9c32ac22d4c813360be3c6b4fa00f51f58 | 6,449,865 |
| manga_dynpop_DA_target.csv | 7076255c0b716ca7def8c65bc6de0163e84d6d3d849a78c762b2f596d49a99ec | 298,792 |

These are the hashes of the files as committed, with LF line endings.
