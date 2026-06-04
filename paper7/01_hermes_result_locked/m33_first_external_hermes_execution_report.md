# M33 First External Hermes Execution Report

**Lane:** M33-Corbelli source-native thick-disk reconstruction  
**Execution:** Hermes only. MOND was not run.  
**Board status:** Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; final board-freeze package complete; age status Amber with bracketed t50 locked for sensitivity use.

## Run constraints preserved

- Source-native Corbelli SPS/pixel-SED stellar mass profile; radially varying M/L provenance.
- Gas profile uses helium factor 1.33 and the locked H2 prescription.
- Flaring stellar disk and gas half-thickness 0.5 kpc are preserved in the reconstructed board provenance.
- Vbul = 0.
- No changes were made to M/L, age, gas scaling, kernel parameters, distance, inclination, or uncertainties after authorization.
- The annular thick-disk kernel was not altered after the result.
- MOND was not run.

## Frozen Hermes execution settings

| Quantity | Value |
|---|---:|
| Radial points | 58 |
| g98 | 2653.950 (km/s)^2/kpc |
| a_knee | 1585 (km/s)^2/kpc |
| r_knee | 2.69 kpc |
| Gate window | 11 points |
| phi_last | 0.519 |
| median phi | 0.653 |
| sigma_int^2 | 386 |

## Bracketed-age Hermes results

| Age lane | t50 | beta | chi2 | reduced chi2 | Median abs residual | Mean signed residual | Pass/fail |
|---|---:|---:|---:|---:|---:|---:|---|
| Primary/source-native mass-weighted | 6.0 ± 1.2 Gyr | 1.834 | 181.30 | **3.126** | 29.9 km/s | -27.8 km/s | PASS / non-catastrophic |
| Conservative older-weighted bracket | 7.1 ± 1.3 Gyr | 1.699 | 191.15 | **3.296** | 31.4 km/s | -29.3 km/s | PASS / non-catastrophic |

## Residual structure

Both age lanes are non-catastrophic under the current Hermes criterion of reduced chi2 < 5. The result is not clean: the model underpredicts the outer rotation curve and the fit penalty is concentrated beyond 10 kpc.

| Age lane | Inner 0-3 kpc reduced chi2 | Mid 3-10 kpc reduced chi2 | Outer >10 kpc reduced chi2 | Fraction of total chi2 beyond 10 kpc | Outer mean residual |
|---|---:|---:|---:|---:|---:|
| Primary/source-native mass-weighted | 0.262 | 1.308 | 5.794 | 86.3% | -51.6 km/s |
| Conservative older-weighted bracket | 0.212 | 1.438 | 6.109 | 86.3% | -53.0 km/s |


## Age uncertainty sensitivity

The one-sigma age bounds were carried as diagnostic uncertainty only, not as additional tuned lanes.

| Age lane | Bound | t50 | reduced chi2 | beta |
|---|---:|---:|---:|---:|
| Primary/source-native mass-weighted | minus_1sigma | 4.8 Gyr | 2.946 | 1.992 |
| Primary/source-native mass-weighted | plus_1sigma | 7.2 Gyr | 3.311 | 1.687 |
| Conservative older-weighted bracket | minus_1sigma | 5.8 Gyr | 3.095 | 1.860 |
| Conservative older-weighted bracket | plus_1sigma | 8.4 Gyr | 3.501 | 1.549 |


## Interpretation within authorized scope

This first external Hermes execution is a **pass but elevated** result. The age bracket does not flip the classification: both age lanes remain non-catastrophic, and the younger/source-native age performs slightly better. The difference between the two age lanes is small compared with the outer-disk residual structure, so the age bracket is not the dominant source of tension.

The dominant structure is an outer-disk underprediction. Beyond 10 kpc, the primary-age lane contributes 86.3% of total chi2, with a mean residual of -51.6 km/s. The inner 0-3 kpc region is well within tolerance in both lanes.

This is not broad validation. It is one external Format B source-native reconstructed board. It supports moving M33 to the next comparator gate only if the PI authorizes it.

## Package files

- m33_hermes_age_bracket_results.csv
- m33_hermes_per_point_predictions.csv
- m33_hermes_radial_failure_structure.csv
- m33_hermes_age_uncertainty_sensitivity.csv
- m33_hermes_execution_summary.json
- m33_observed_vs_hermes_two_ages.png
- m33_hermes_residuals_vs_radius.png
- m33_baryonic_components_vbar_context.png
- m33_hermes_chi2_contribution_vs_radius.png
- m33_hermes_gate_profile.png
- m33_corbelli_reconstructed_board_input.csv
- m33_age_bracketed_amber_lock.csv

No MOND result exists for this package.
