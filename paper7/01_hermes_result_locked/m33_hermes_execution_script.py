from __future__ import annotations
import json, math, shutil, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

BASE = Path('/mnt/data/m33_first_external_hermes_execution_v2')
if BASE.exists():
    shutil.rmtree(BASE)
BASE.mkdir(parents=True)
BOARD = Path('/mnt/data/m33_corbelli_final_board_freeze_package_age_bracketed/m33_corbelli_reconstructed_board.csv')
AGECSV = Path('/mnt/data/m33_corbelli_final_board_freeze_package_age_bracketed/m33_age_bracketed_amber_lock.csv')
FREEZE_MEMO = Path('/mnt/data/m33_corbelli_final_board_freeze_package_age_bracketed/m33_corbelli_final_board_freeze_memo.md')
KERNEL_SPEC = Path('/mnt/data/m33_corbelli_final_board_freeze_package_age_bracketed/m33_corbelli_kernel_specification.md')
OVERLAY = Path('/mnt/data/m33_corbelli_final_board_freeze_package_age_bracketed/corbelli_fig12_overlay_reconstruction.png')

SIGMA_INT2 = 386.0
A_KNEE = 1585.0
C_OVER_2PI = 46654.0
BETA_MAX = math.pi - 1.0 / math.sqrt(2.0 * math.pi)
BETA_FLOOR = -1.0 / math.sqrt(2.0 * math.pi)

df = pd.read_csv(BOARD)
R = df['R_kpc'].to_numpy(float)
Vobs = df['Vobs_kms'].to_numpy(float)
errV = df['errV_kms'].to_numpy(float)
Vgas = df['Vgas_kms'].to_numpy(float)
Vdisk = df['Vdisk_kms'].to_numpy(float)
Vbul = df['Vbulge_kms'].to_numpy(float)
Vbar = np.sqrt(np.maximum(0.0, Vdisk**2 + Vbul**2 + np.sign(Vgas)*Vgas**2))
gbar = Vbar**2 / R

def gate_func():
    below = np.where(gbar < A_KNEE)[0]
    knee_idx = int(below[0]) if len(below) else int(np.argmin(np.abs(gbar - A_KNEE)))
    r_knee = float(R[knee_idx])
    x = R / r_knee
    win = 11 if len(R) >= 20 else 5
    if win >= len(R):
        win = len(R) - 1 if (len(R) - 1) % 2 else len(R) - 2
    if win % 2 == 0: win -= 1
    Vsm = savgol_filter(Vbar, win, 3, mode='interp')
    shear = np.abs(np.gradient(np.log(np.maximum(Vsm, 1e-12)), np.log(R)))
    dVdR = np.gradient(Vsm, R)
    d2VdR2 = np.gradient(dVdR, R)
    curvature = np.abs(d2VdR2) / (np.abs(dVdR) + 1e-6)
    inner = x <= 1.0
    outer = x > 1.0
    curv_med = float(np.nanmedian(curvature[inner])) if np.any(inner) else float(np.nanmedian(curvature))
    shear_med = float(np.nanmedian(shear[outer])) if np.any(outer) else float(np.nanmedian(shear))
    if not np.isfinite(curv_med) or curv_med <= 0:
        vals = curvature[curvature > 0]
        curv_med = float(np.nanmedian(vals)) if vals.size else 1.0
    if not np.isfinite(shear_med) or shear_med <= 0:
        vals = shear[shear > 0]
        shear_med = float(np.nanmedian(vals)) if vals.size else 1.0
    curv_norm = curvature / curv_med
    shear_norm = shear / shear_med
    w = 1.0 / (1.0 + np.exp(-((x - 1.40) / 0.30)))
    S_pre = (1.0 - w) * curv_norm + w * shear_norm
    S_eff = np.empty_like(S_pre)
    for i, s in enumerate(S_pre):
        S_eff[i] = s if i == 0 else 0.5 * s + 0.5 * S_eff[i-1]
    phi = 1.0 - np.exp(-S_eff)
    return dict(knee_idx=knee_idx, r_knee=r_knee, x=x, window=win, Vsm=Vsm, shear=shear,
                curvature=curvature, curv_norm=curv_norm, shear_norm=shear_norm, w=w,
                S_pre=S_pre, S_eff=S_eff, phi=phi, curv_med=curv_med, shear_med=shear_med)

gate = gate_func()
phi = gate['phi']
g98 = float(np.percentile(gbar, 98))

def beta_for_t(t50):
    psi = t50 * g98 / C_OVER_2PI
    beta = math.pi * math.exp(-psi) - 1.0 / math.sqrt(2.0 * math.pi)
    return float(psi), float(beta)

def run(t50):
    psi, beta = beta_for_t(t50)
    Vmodel = Vbar * np.sqrt(np.maximum(0.0, 1.0 + beta * phi))
    residual = Vmodel - Vobs
    chi_i = residual**2 / (errV**2 + SIGMA_INT2)
    rb_open = Vobs**2 / Vbar**2 - 1.0
    rb_actual = np.where(phi > 1e-12, rb_open / phi, np.inf)
    bands = []
    for label, lo, hi in [('inner_0_3_kpc',0,3), ('mid_3_10_kpc',3,10), ('outer_gt_10_kpc',10,np.inf)]:
        m = (R >= lo) & (R < hi)
        bands.append(dict(
            band=label, N=int(m.sum()), chi2_sum=float(chi_i[m].sum()),
            chi2_reduced=float(chi_i[m].mean()), chi2_fraction=float(chi_i[m].sum()/chi_i.sum()),
            mean_residual_kms=float(residual[m].mean()), median_abs_residual_kms=float(np.median(np.abs(residual[m])))
        ))
    summary = dict(
        t50_Gyr=float(t50), psi=psi, beta=beta, N_points=len(R), chi2=float(chi_i.sum()),
        chi2_reduced=float(chi_i.mean()), pass_fail='PASS / non-catastrophic' if chi_i.mean() < 5 else 'FAIL / catastrophic',
        median_abs_residual_kms=float(np.median(np.abs(residual))), mean_abs_residual_kms=float(np.mean(np.abs(residual))),
        rms_residual_kms=float(np.sqrt(np.mean(residual**2))), mean_signed_residual_kms=float(np.mean(residual)),
        median_signed_residual_kms=float(np.median(residual)), max_abs_residual_kms=float(np.max(np.abs(residual))),
        fraction_model_under_observed=float(np.mean(residual < 0)), outer_chi2_fraction=float(bands[2]['chi2_fraction']),
        outer_mean_residual_kms=float(bands[2]['mean_residual_kms']), max_required_beta_gate_open=float(np.nanmax(rb_open)),
        max_required_beta_actual_gate=float(np.nanmax(rb_actual[np.isfinite(rb_actual)]))
    )
    return dict(summary=summary, Vmodel=Vmodel, residual=residual, chi_i=chi_i,
                required_beta_open=rb_open, required_beta_actual_gate=rb_actual, bands=bands)

lanes = [
    ('Primary/source-native mass-weighted', 6.0, 1.2, 'primary_t60'),
    ('Conservative older-weighted bracket', 7.1, 1.3, 'conservative_t71')
]
runs = {t: run(t) for _, t, _, _ in lanes}

# write board/provenance copies
shutil.copy2(BOARD, BASE/'m33_corbelli_reconstructed_board_input.csv')
if AGECSV.exists(): shutil.copy2(AGECSV, BASE/'m33_age_bracketed_amber_lock.csv')
if FREEZE_MEMO.exists(): shutil.copy2(FREEZE_MEMO, BASE/'m33_corbelli_final_board_freeze_memo.md')
if KERNEL_SPEC.exists(): shutil.copy2(KERNEL_SPEC, BASE/'m33_corbelli_kernel_specification.md')
if OVERLAY.exists(): shutil.copy2(OVERLAY, BASE/'corbelli_fig12_overlay_reconstruction.png')

# CSVs
summary_rows=[]; band_rows=[]; uncertainty_rows=[]
per = df.copy()
per['gbar_kms2_per_kpc']=gbar
per['phi']=phi
per['x_R_over_rknee']=gate['x']
per['Vbar_smoothed_kms']=gate['Vsm']
per['gate_shear']=gate['shear']
per['gate_curvature']=gate['curvature']
per['gate_S_eff']=gate['S_eff']
for lane, t, sigma, key in lanes:
    s = runs[t]['summary'].copy()
    s.update(dict(lane_label='M33-Corbelli source-native thick-disk reconstruction', age_lane=lane, t50_sigma_Gyr=sigma,
                  g98_kms2_per_kpc=g98, r_knee_kpc=gate['r_knee'], phi_last=float(phi[-1]),
                  phi_median=float(np.median(phi)), phi_min=float(phi.min()), phi_max=float(phi.max())))
    summary_rows.append(s)
    for b in runs[t]['bands']:
        row = b.copy(); row.update(dict(age_lane=lane, t50_Gyr=t)); band_rows.append(row)
    for bound, tb in [('minus_1sigma', max(0.0,t-sigma)), ('plus_1sigma', t+sigma)]:
        bs = run(tb)['summary']; bs.update(dict(age_lane=lane, central_t50_Gyr=t, bound=bound)); uncertainty_rows.append(bs)
    per[f'Vhermes_{key}_kms'] = runs[t]['Vmodel']
    per[f'residual_model_minus_obs_{key}_kms'] = runs[t]['residual']
    per[f'chi2_contribution_{key}'] = runs[t]['chi_i']
    per[f'required_beta_open_{key}'] = runs[t]['required_beta_open']
    per[f'required_beta_actual_gate_{key}'] = runs[t]['required_beta_actual_gate']
summary_df = pd.DataFrame(summary_rows)
band_df = pd.DataFrame(band_rows)
uncert_df = pd.DataFrame(uncertainty_rows)
summary_df.to_csv(BASE/'m33_hermes_age_bracket_results.csv', index=False)
band_df.to_csv(BASE/'m33_hermes_radial_failure_structure.csv', index=False)
uncert_df.to_csv(BASE/'m33_hermes_age_uncertainty_sensitivity.csv', index=False)
per.to_csv(BASE/'m33_hermes_per_point_predictions.csv', index=False)

# plots
plt.figure(figsize=(8,5.5))
plt.errorbar(R, Vobs, yerr=errV, fmt='o', label='Observed Vobs')
plt.plot(R, Vbar, label='Vbar')
for lane, t, _, _ in lanes:
    plt.plot(R, runs[t]['Vmodel'], label=f'Hermes {t:.1f} Gyr')
plt.xlabel('R (kpc)'); plt.ylabel('Velocity (km/s)'); plt.title('M33-Corbelli: observed vs Hermes prediction')
plt.legend(); plt.tight_layout(); plt.savefig(BASE/'m33_observed_vs_hermes_two_ages.png', dpi=200); plt.close()

plt.figure(figsize=(8,5.5))
for lane, t, _, _ in lanes:
    plt.plot(R, runs[t]['residual'], marker='o', label=f'Residual {t:.1f} Gyr')
plt.axhline(0, linestyle='--')
plt.xlabel('R (kpc)'); plt.ylabel('Vmodel - Vobs (km/s)'); plt.title('M33-Corbelli: Hermes residuals vs radius')
plt.legend(); plt.tight_layout(); plt.savefig(BASE/'m33_hermes_residuals_vs_radius.png', dpi=200); plt.close()

plt.figure(figsize=(8,5.5))
plt.plot(R, Vgas, label='Vgas')
plt.plot(R, Vdisk, label='Vdisk')
plt.plot(R, Vbul, label='Vbulge = 0')
plt.plot(R, Vbar, label='Vbar')
plt.xlabel('R (kpc)'); plt.ylabel('Velocity component (km/s)'); plt.title('M33-Corbelli: baryonic components and Vbar')
plt.legend(); plt.tight_layout(); plt.savefig(BASE/'m33_baryonic_components_vbar_context.png', dpi=200); plt.close()

plt.figure(figsize=(8,5.5))
for lane, t, _, _ in lanes:
    plt.plot(R, runs[t]['chi_i'], marker='o', label=f'Chi2 contribution {t:.1f} Gyr')
plt.xlabel('R (kpc)'); plt.ylabel('Per-point chi2 contribution'); plt.title('M33-Corbelli: chi2 contribution vs radius')
plt.legend(); plt.tight_layout(); plt.savefig(BASE/'m33_hermes_chi2_contribution_vs_radius.png', dpi=200); plt.close()

plt.figure(figsize=(8,5.5))
plt.plot(R, phi, marker='o', label='phi(R)')
plt.xlabel('R (kpc)'); plt.ylabel('Gate phi'); plt.title('M33-Corbelli: Hermes gate profile')
plt.legend(); plt.tight_layout(); plt.savefig(BASE/'m33_hermes_gate_profile.png', dpi=200); plt.close()

# metadata
meta = dict(
    lane_label='M33-Corbelli source-native thick-disk reconstruction', execution='Hermes only; MOND not run',
    status='Format B reconstructed board; mass closure passed; overlay benchmark passed; annular thick-disk kernel provisionally frozen; bracketed Amber age lock',
    N_points=len(R), g98_kms2_per_kpc=g98, a_knee_kms2_per_kpc=A_KNEE, r_knee_kpc=gate['r_knee'],
    sigma_int2=SIGMA_INT2, beta_max=BETA_MAX, beta_floor=BETA_FLOOR, phi_last=float(phi[-1]),
    phi_median=float(np.median(phi)), phi_min=float(phi.min()), phi_max=float(phi.max()),
    gate_normalization='median-zone curvature/shear normalization as in Hermes equation paper',
    knee_rule='first radius where gbar drops below a_knee; nearest fallback if no crossing',
    results=summary_rows
)
(BASE/'m33_hermes_execution_summary.json').write_text(json.dumps(meta, indent=2))

# report
primary = summary_rows[0]
older = summary_rows[1]
# helper band lookup
bd = {(r['age_lane'], r['band']): r for r in band_rows}
report = f"""# M33 First External Hermes Execution Report

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
| Radial points | {len(R)} |
| g98 | {g98:.3f} (km/s)^2/kpc |
| a_knee | {A_KNEE:.0f} (km/s)^2/kpc |
| r_knee | {gate['r_knee']:.2f} kpc |
| Gate window | {gate['window']} points |
| phi_last | {phi[-1]:.3f} |
| median phi | {np.median(phi):.3f} |
| sigma_int^2 | {SIGMA_INT2:.0f} |

## Bracketed-age Hermes results

| Age lane | t50 | beta | chi2 | reduced chi2 | Median abs residual | Mean signed residual | Pass/fail |
|---|---:|---:|---:|---:|---:|---:|---|
| Primary/source-native mass-weighted | {primary['t50_Gyr']:.1f} ± {primary['t50_sigma_Gyr']:.1f} Gyr | {primary['beta']:.3f} | {primary['chi2']:.2f} | **{primary['chi2_reduced']:.3f}** | {primary['median_abs_residual_kms']:.1f} km/s | {primary['mean_signed_residual_kms']:.1f} km/s | {primary['pass_fail']} |
| Conservative older-weighted bracket | {older['t50_Gyr']:.1f} ± {older['t50_sigma_Gyr']:.1f} Gyr | {older['beta']:.3f} | {older['chi2']:.2f} | **{older['chi2_reduced']:.3f}** | {older['median_abs_residual_kms']:.1f} km/s | {older['mean_signed_residual_kms']:.1f} km/s | {older['pass_fail']} |

## Residual structure

Both age lanes are non-catastrophic under the current Hermes criterion of reduced chi2 < 5. The result is not clean: the model underpredicts the outer rotation curve and the fit penalty is concentrated beyond 10 kpc.

| Age lane | Inner 0-3 kpc reduced chi2 | Mid 3-10 kpc reduced chi2 | Outer >10 kpc reduced chi2 | Fraction of total chi2 beyond 10 kpc | Outer mean residual |
|---|---:|---:|---:|---:|---:|
"""
for lane, t, _, _ in lanes:
    inner = bd[(lane,'inner_0_3_kpc')]
    mid = bd[(lane,'mid_3_10_kpc')]
    out = bd[(lane,'outer_gt_10_kpc')]
    report += f"| {lane} | {inner['chi2_reduced']:.3f} | {mid['chi2_reduced']:.3f} | {out['chi2_reduced']:.3f} | {100*out['chi2_fraction']:.1f}% | {out['mean_residual_kms']:.1f} km/s |\n"
report += f"""

## Age uncertainty sensitivity

The one-sigma age bounds were carried as diagnostic uncertainty only, not as additional tuned lanes.

| Age lane | Bound | t50 | reduced chi2 | beta |
|---|---:|---:|---:|---:|
"""
for r in uncertainty_rows:
    report += f"| {r['age_lane']} | {r['bound']} | {r['t50_Gyr']:.1f} Gyr | {r['chi2_reduced']:.3f} | {r['beta']:.3f} |\n"
report += f"""

## Interpretation within authorized scope

This first external Hermes execution is a **pass but elevated** result. The age bracket does not flip the classification: both age lanes remain non-catastrophic, and the younger/source-native age performs slightly better. The difference between the two age lanes is small compared with the outer-disk residual structure, so the age bracket is not the dominant source of tension.

The dominant structure is an outer-disk underprediction. Beyond 10 kpc, the primary-age lane contributes {100*primary['outer_chi2_fraction']:.1f}% of total chi2, with a mean residual of {primary['outer_mean_residual_kms']:.1f} km/s. The inner 0-3 kpc region is well within tolerance in both lanes.

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
"""
(BASE/'m33_first_external_hermes_execution_report.md').write_text(report)

# include this script
shutil.copy2(Path(__file__), BASE/'m33_hermes_execution_script.py')
# zip
zip_path=Path('/mnt/data/m33_first_external_hermes_execution_v2_package.zip')
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as zf:
    for p in sorted(BASE.iterdir()):
        zf.write(p, arcname=p.name)
print(summary_df := pd.DataFrame(summary_rows).to_string(index=False))
print(zip_path)
