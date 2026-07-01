# MaNGA sledgehammer young-vs-old T50 quartile test

Conservative sample: `Qual >= 1`, finite `T50`, finite stellar mass, finite `DML_mfl_int`, finite `DML_mfl_obs`. Decision sample N = 5,952. Eight equal-count stellar-mass bins; each bin has 744 galaxies, with the lowest-T50 186 galaxies labeled young and highest-T50 186 labeled old. Ties are broken deterministically by stellar mass. Larger T50 means older assembly.

## Primary scoreboard

| outcome     | outcome_label     |   real_EP_favorable_bins |   real_mean_bin_median_diff_young_minus_old |   shuffle_p_count_ge_real |   shuffle_p_mean_diff_ge_real |   shuffle_p_joint_count_and_mean_ge_real |   mass_adjusted_median_diff_young_minus_old |   diff_bootstrap95_lo |   diff_bootstrap95_hi |
|:------------|:------------------|-------------------------:|--------------------------------------------:|--------------------------:|------------------------------:|-----------------------------------------:|--------------------------------------------:|----------------------:|----------------------:|
| DML_mfl_int | DML intrinsic SPS |                        6 |                                   0.0606433 |                0.14857    |                    0.00019996 |                               0.00019996 |                                   0.0772407 |             0.0561809 |             0.0958905 |
| DML_mfl_obs | DML observed SPS  |                        8 |                                   0.0741693 |                0.00559888 |                    0.00019996 |                               0.00019996 |                                   0.0668992 |             0.0498353 |             0.0827721 |


## Per-bin primary results

| outcome     |   mass_bin |   N_young |   N_old |   T50_young_median |   T50_old_median |   median_DML_young |   median_DML_old |   median_diff_young_minus_old |   bootstrap95_lo |   bootstrap95_hi | EP_favorable_young_gt_old   |
|:------------|-----------:|----------:|--------:|-------------------:|-----------------:|-------------------:|-----------------:|------------------------------:|-----------------:|-----------------:|:----------------------------|
| DML_mfl_int |          0 |       186 |     186 |            4.25606 |          12.8485 |        0.563867    |       0.430542   |                     0.133324  |      0.0686079   |        0.225058  | True                        |
| DML_mfl_int |          1 |       186 |     186 |            5.3     |          12.8485 |        0.309025    |       0.327705   |                    -0.0186803 |     -0.112428    |        0.109578  | False                       |
| DML_mfl_int |          2 |       186 |     186 |            5.78182 |          12.8485 |        0.218987    |       0.259625   |                    -0.0406379 |     -0.105274    |        0.0707398 | False                       |
| DML_mfl_int |          3 |       186 |     186 |            6.10303 |          12.8485 |        0.179053    |       0.143383   |                     0.0356702 |     -0.0414417   |        0.0925018 | True                        |
| DML_mfl_int |          4 |       186 |     186 |            5.62121 |          12.8485 |        0.215574    |       0.102719   |                     0.112856  |      0.0677233   |        0.161571  | True                        |
| DML_mfl_int |          5 |       186 |     186 |            6.10303 |          12.8485 |        0.135775    |       0.061404   |                     0.074371  |      0.0197121   |        0.118264  | True                        |
| DML_mfl_int |          6 |       186 |     186 |            6.10303 |          13.0091 |        0.180274    |       0.0766563  |                     0.103617  |      0.0690395   |        0.144588  | True                        |
| DML_mfl_int |          7 |       186 |     186 |            6.74545 |          12.8485 |        0.18033     |       0.0957037  |                     0.0846267 |      0.0602627   |        0.117385  | True                        |
| DML_mfl_obs |          0 |       186 |     186 |            4.25606 |          12.8485 |        0.3728      |       0.184675   |                     0.188126  |      0.119368    |        0.267763  | True                        |
| DML_mfl_obs |          1 |       186 |     186 |            5.3     |          12.8485 |        0.145721    |       0.0947855  |                     0.0509359 |     -0.0279588   |        0.137033  | True                        |
| DML_mfl_obs |          2 |       186 |     186 |            5.78182 |          12.8485 |        0.0685274   |       0.0220269  |                     0.0465006 |     -0.000609219 |        0.0848913 | True                        |
| DML_mfl_obs |          3 |       186 |     186 |            6.10303 |          12.8485 |        0.017228    |      -0.0324115  |                     0.0496395 |      0.00363844  |        0.104304  | True                        |
| DML_mfl_obs |          4 |       186 |     186 |            5.62121 |          12.8485 |        0.0537019   |      -0.0389164  |                     0.0926183 |      0.0469419   |        0.128539  | True                        |
| DML_mfl_obs |          5 |       186 |     186 |            6.10303 |          12.8485 |        0.000330569 |      -0.0523742  |                     0.0527047 |      0.0101258   |        0.0806636 | True                        |
| DML_mfl_obs |          6 |       186 |     186 |            6.10303 |          13.0091 |        0.0293466   |      -0.0285443  |                     0.0578909 |      0.0212343   |        0.0964212 | True                        |
| DML_mfl_obs |          7 |       186 |     186 |            6.74545 |          12.8485 |        0.063036    |       0.00809693 |                     0.054939  |      0.023775    |        0.0840595 | True                        |


## Overall mass-adjusted medians

| outcome     |   mass_adjusted_median_young |   mass_adjusted_median_old |   mass_adjusted_median_diff_young_minus_old |   diff_bootstrap95_lo |   diff_bootstrap95_hi |
|:------------|-----------------------------:|---------------------------:|--------------------------------------------:|----------------------:|----------------------:|
| DML_mfl_int |                    0.0757515 |                -0.00148919 |                                   0.0772407 |             0.0561809 |             0.0958905 |
| DML_mfl_obs |                    0.0544588 |                -0.0124404  |                                   0.0668992 |             0.0498353 |             0.0827721 |


## Optional paired nearest-mass check

| outcome     | outcome_label     |   n_pairs |   young_wins |   young_win_fraction |   binomial_p_greater_0p5 |   shuffle_p_win_fraction_ge_real |   shuffle_win_fraction_mean |   shuffle_win_fraction_sd |
|:------------|:------------------|----------:|-------------:|---------------------:|-------------------------:|---------------------------------:|----------------------------:|--------------------------:|
| DML_mfl_int | DML intrinsic SPS |      1488 |          881 |             0.59207  |              6.37974e-13 |                       0.00019996 |                    0.499815 |                 0.0128214 |
| DML_mfl_obs | DML observed SPS  |      1488 |          913 |             0.613575 |              8.63145e-19 |                       0.00019996 |                    0.499849 |                 0.012874  |


## Secondary fDM checks

| outcome         | outcome_label   |   real_EP_favorable_bins |   real_mean_bin_median_diff_young_minus_old |   shuffle_p_joint_count_and_mean_ge_real |
|:----------------|:----------------|-------------------------:|--------------------------------------------:|-----------------------------------------:|
| gnfw_cyl_fdm_Re | gNFW fDM(<Re)   |                        5 |                                   0.0010625 |                                 0.30074  |
| nfw_cyl_fdm_Re  | NFW fDM(<Re)    |                        4 |                                  -0.0035    |                                 0.513497 |


## Caveat

This is a blunt visualization and decision test, not proof of age-dependent gravity. DML is constructed from JAM dynamical mass-to-light minus SPS stellar mass-to-light, so age-metallicity covariance, SPS M/L behavior, IMF assumptions, dust conventions, and shared kinematic inputs remain live confounds. The result says MaNGA is worth hardening; it does not by itself settle the physics.
