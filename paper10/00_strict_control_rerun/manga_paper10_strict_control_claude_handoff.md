# Paper 10 strict control rerun — Claude handoff

Primary strict control specification chosen by Pro: mass-weighted metallicity (`sp_MW_Metal_Re`), MGE ellipticity (`Eps_MGE`), and `logSigma_Re` only. This removes the v0.7 collinearity issue from including raw `Sigma_Re` and `logSigma_Re` together, and also removes duplicate metallicity and duplicate flattening terms.

Outcome: `mfl_cyl_log_ML_dyn`. Age split: `sp_T50`, larger means older. Mass bins: `nsa_sersic_mass`. Sample: `Qual >= 1` and finite required values, N = 5,952; eight equal-count mass bins of 744 galaxies; 186 young and 186 old galaxies per bin.

Controls for intrinsic-SPS set: `sp_ML_int_Re`, `sp_MW_Metal_Re`, `nsa_sersic_mass`, `logRe`, `nsa_sersic_n`, `Eps_MGE`, `Lambda_Re`, `logSigma_Re`.

Controls for observed-SPS set: `sp_ML_obs_Re`, `sp_MW_Metal_Re`, `nsa_sersic_mass`, `logRe`, `nsa_sersic_n`, `Eps_MGE`, `Lambda_Re`, `logSigma_Re`.

Controls for both-SPS set: both SPS M/L terms plus the same shared controls.

## Main result

| control_set   |    N |   young_gt_old_bins |   total_bins |   equal_bin_mean_diff_young_minus_old |   equal_bin_mean_boot95_lo |   equal_bin_mean_boot95_hi |   mass_adjusted_diff_young_minus_old |   mass_adjusted_diff_boot95_lo |   mass_adjusted_diff_boot95_hi |   shuffle_p_count_ge_real |   shuffle_p_sum_ge_real |   shuffle_p_joint_count_and_sum_ge_real |
|:--------------|-----:|--------------------:|-------------:|--------------------------------------:|---------------------------:|---------------------------:|-------------------------------------:|-------------------------------:|-------------------------------:|--------------------------:|------------------------:|----------------------------------------:|
| intSPS        | 5952 |                   7 |            8 |                               0.02378 |                    0.01318 |                    0.03366 |                              0.02247 |                        0.01139 |                        0.0342  |                   0.03679 |                  0.0002 |                                  0.0002 |
| obsSPS        | 5952 |                   7 |            8 |                               0.0393  |                    0.029   |                    0.04982 |                              0.03866 |                        0.02916 |                        0.04944 |                   0.03119 |                  0.0002 |                                  0.0002 |
| bothSPS       | 5952 |                   7 |            8 |                               0.03467 |                    0.02362 |                    0.04514 |                              0.0344  |                        0.02309 |                        0.04425 |                   0.03399 |                  0.0002 |                                  0.0002 |

## Age term / CV / rank checks

| control_set   |   cv_delta_R2_from_age |   cv_delta_R2_sd |   resid_T50_spearman_rho |   resid_T50_spearman_p |   age_T50_coef_raw_outcome_per_1sd_T50 |   age_T50_robust_se |   age_T50_robust_p |
|:--------------|-----------------------:|-----------------:|-------------------------:|-----------------------:|---------------------------------------:|--------------------:|-------------------:|
| intSPS        |               0.008582 |         0.005742 |                 -0.06491 |              5.395e-07 |                               -0.02673 |            0.003278 |          3.582e-16 |
| obsSPS        |               0.02092  |         0.009314 |                 -0.1167  |              1.669e-19 |                               -0.0404  |            0.003434 |          6.046e-32 |
| bothSPS       |               0.01921  |         0.01021  |                 -0.1033  |              1.358e-15 |                               -0.04161 |            0.00359  |          4.663e-31 |

## Per-bin note

The same mass bin (bin 5; median logM ≈ 10.663) is the lone negative bin in all three strict control sets. Its bootstrap interval crosses zero in all three cases, so the headline should be **7/8**, not 8/8, with one knife-edge/non-decisive bin.

## Sensitivity note

A quick MW/LW metallicity × Eps/MGE-vs-NSA-b/a sensitivity grid, without full bootstrap/shuffle, gives 7/8 bins for all 12 reduced-spec variants. The reduced-spec family is therefore stable at 7/8. The original kitchen-sink 8/8 result should be supplementary sensitivity only.

## Recommendation for main text

Use the strict reduced model as the primary controlled-residual result: young galaxies have higher controlled dynamical M/L residuals in **7/8 mass bins** for intrinsic, observed, and both-SPS control sets. The magnitude test is at the shuffle floor for all three; the count-only p is ~0.03–0.04, so the directional count is now significant but should still be described as modest and one-bin-sensitive.
