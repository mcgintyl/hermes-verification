"""Paper 10 size / distance / subsample sensitivity of the controlled dynamical residual.

Written October 2026 for v3 of the paper. Every fit, mass bin, quartile split, bootstrap and
shuffle comes from the package's own functions in ../00_strict_control_rerun/
run_clean_controls_fast.py (make_base, fit_controls_resid, age_coef_raw_outcome, add_mass_bins,
get_bin_arrays, eval_sledge), at full N_BOOT = 5,000 and N_SHUFFLE = 5,000.

Why: the published "size" control logRe is log10 of the MGE effective radius in ARCSEC, an
angular size. At fixed stellar mass the young quartile lies 1.2 to 1.7 times farther away than
the old quartile (median angular-diameter distance, every bin) and is 2 to 3.5 times more often a
MaNGA Secondary-sample target, so an angular size term carries distance. This script reruns the
strict (reduced) specification with

  a  published strict spec (logRe in arcsec)                       [reproduces Table 1 exactly]
  b  size in kpc: logRe_kpc = logRe + log10(DA * 1000 / 206264.806)
  c  published spec + logDA (log10 of DA in Mpc)
  d  kpc size + logDA (spans the same space as c, so the residuals are identical)
  e  MaNGA Primary-only (target = 0) and Secondary-only (target = 1) subsamples, refit within the
     subsample: published spec, kpc size, and published spec + logDA
  f  the original kitchen-sink control set of the v0.6 hardening run (adds sp_LW_Metal_Re,
     nsa_sersic_ba and raw Sigma_Re): reproduces appendix Section 8 / Table 3 (6/8, 8/8, 8/8)

and reports, for each SPS set (int, obs, both): bins with young > old, equal-bin mean and
mass-adjusted young-old residual with bootstrap 95% CIs, count / magnitude / joint shuffle p,
the fully standardized T50 coefficient (HC3), the T50 coefficient on the raw outcome per 1 SD and
per Gyr of T50 (HC3 t and p), cross-validated Delta R^2, and the Spearman rho of residual vs T50.
It also writes the per-bin median DA, Secondary fraction and SPS attenuation term
(sp_ML_obs_Re - sp_ML_int_Re = DML_int - DML_obs) of the young and old quartiles, the 12-variant
reduced-spec grid under each size treatment (bin counts plus the HC3 t and p of a linear T50
term), and one-cube-per-MaNGA-ID point estimates for every Table 1 row, with the controlled
residual in both size units (arcsec and kpc). tie_rule_redraws.csv repeats the young/old split with
the T50 ties at the quartile cuts broken at random instead of by stellar mass (200 draws per cell,
seed SEED + 9000) for the grid and the subsample rows: bin counts and mass-adjusted differences only.

Inputs: ../data/manga_dynpop_merged_thin_firstpass.csv and ../data/manga_dynpop_DA_target.csv
(plateifu, mangaid, DA, target; copied from SDSSDR17_MaNGA_JAM_v2.fits HDU 1 by
../build_merged_table.py). The JAM v2 data model defines target as "Flag for subsample of MaNGA
(Primary: 0, Secondary: 1, colour-Enhanced: 2)" and DA as the adopted angular-diameter distance
(flat Planck 2015 cosmology). Rows are joined on plateifu.

Random numbers: every variant uses the published per-set seed SEED + 1000*i (i = 0 int, 1 obs,
2 both; SEED = 20260629), i.e. common random numbers across variants, so variant a reproduces the
published scoreboard bit for bit. Subsample runs use eight near-equal mass bins (np.array_split)
and quartiles of len(bin)//4.

Usage:  python run_size_distance_sensitivity.py [--input CSV] [--da CSV] [--outdir DIR]
"""
import argparse
import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location('rc', PKG / '00_strict_control_rerun' / 'run_clean_controls_fast.py')
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

DEFAULT_IN = PKG / 'data' / 'manga_dynpop_merged_thin_firstpass.csv'
DEFAULT_DA = PKG / 'data' / 'manga_dynpop_DA_target.csv'
ARCSEC_PER_RAD = 206264.806
SPS_SETS = [('int', ['sp_ML_int_Re']), ('obs', ['sp_ML_obs_Re']), ('both', ['sp_ML_int_Re', 'sp_ML_obs_Re'])]
KITCHEN_SINK_SHARED = ['sp_MW_Metal_Re', 'sp_LW_Metal_Re', 'nsa_sersic_mass', 'logRe', 'nsa_sersic_n',
                       'nsa_sersic_ba', 'Eps_MGE', 'Lambda_Re', 'Sigma_Re', 'logSigma_Re']
SAMPLE_REQUIRED = ['sp_T50', 'nsa_sersic_mass', 'DML_mfl_int', 'DML_mfl_obs']

# key, description, sample, size column, add logDA, control family
VARIANTS = [
    ('a_published', 'published strict spec (size = logRe, arcsec)', 'all', 'logRe', False, 'strict'),
    ('b_kpc', 'size in kpc', 'all', 'logRe_kpc', False, 'strict'),
    ('c_published_plus_logDA', 'published spec + log distance', 'all', 'logRe', True, 'strict'),
    ('d_kpc_plus_logDA', 'size in kpc + log distance', 'all', 'logRe_kpc', True, 'strict'),
    ('e1_primary', 'MaNGA Primary only, published spec', 'primary', 'logRe', False, 'strict'),
    ('e2_secondary', 'MaNGA Secondary only, published spec', 'secondary', 'logRe', False, 'strict'),
    ('e3_primary_kpc', 'MaNGA Primary only, size in kpc', 'primary', 'logRe_kpc', False, 'strict'),
    ('e4_secondary_kpc', 'MaNGA Secondary only, size in kpc', 'secondary', 'logRe_kpc', False, 'strict'),
    ('e5_primary_plus_logDA', 'MaNGA Primary only, published spec + log distance', 'primary', 'logRe', True, 'strict'),
    ('e6_secondary_plus_logDA', 'MaNGA Secondary only, published spec + log distance', 'secondary', 'logRe', True, 'strict'),
    ('f_kitchen_sink', 'original kitchen-sink controls (v0.6 hardening run)', 'all', 'logRe', False, 'kitchen_sink'),
]
TARGET_OF = {'primary': 0, 'secondary': 1}


def load(input_csv, da_csv):
    df = pd.read_csv(input_csv)
    da = pd.read_csv(da_csv)
    assert da['plateifu'].is_unique and len(da) == len(df)
    m = df[['plateifu']].merge(da, on='plateifu', how='left', validate='one_to_one')
    assert m['DA'].notna().sum() > 0 and (m['mangaid'].astype(str).values == df['mangaid'].astype(str).values).all()
    df['DA'] = m['DA'].to_numpy(float)
    df['target'] = m['target'].to_numpy(int)
    with np.errstate(divide='ignore', invalid='ignore'):
        df['logDA'] = np.log10(df['DA'])
        df['logRe_kpc'] = df['logRe'] + np.log10(df['DA'] * 1000.0 / ARCSEC_PER_RAD)
    df['is_secondary'] = (df['target'] == 1).astype(float)
    return df


def controls_for(sps, size_col, add_logda, family):
    shared = rc.SHARED_CONTROLS if family == 'strict' else KITCHEN_SINK_SHARED
    c = sps + [size_col if x == 'logRe' else x for x in shared]
    return c + (['logDA'] if add_logda else [])


def run_variant(df, key, desc, sample, size_col, add_logda, family, i, sps_label, sps):
    d0 = df if sample == 'all' else df[df['target'] == TARGET_OF[sample]]
    controls = controls_for(sps, size_col, add_logda, family)
    base = rc.make_base(d0, controls)
    resid, info = rc.fit_controls_resid(base, controls)
    binned = rc.add_mass_bins(resid)
    seed = rc.SEED + 1000 * i
    score, per = rc.eval_sledge(binned, 'resid', seed=seed)
    raw = rc.age_coef_raw_outcome(base, controls)
    sd_t50 = float(base['sp_T50'].std(ddof=0))
    row = {
        'variant': key, 'description': desc, 'sample': sample, 'size_term': size_col, 'logDA_control': add_logda,
        'controls': family, 'SPS_set': sps_label, 'seed': seed, 'N': score['N'],
        'young_gt_old_bins': score['young_gt_old_bins'], 'total_bins': score['total_bins'],
        'equal_bin_mean_diff': score['equal_bin_mean_diff_young_minus_old'],
        'equal_bin_mean_ci_lo': score['equal_bin_mean_boot95_lo'], 'equal_bin_mean_ci_hi': score['equal_bin_mean_boot95_hi'],
        'mass_adjusted_diff': score['mass_adjusted_diff_young_minus_old'],
        'mass_adjusted_ci_lo': score['mass_adjusted_diff_boot95_lo'], 'mass_adjusted_ci_hi': score['mass_adjusted_diff_boot95_hi'],
        'p_count': score['shuffle_p_count_ge_real'], 'p_magnitude': score['shuffle_p_sum_ge_real'],
        'p_joint': score['shuffle_p_joint_count_and_sum_ge_real'],
        'T50_coef_standardized': info['age_T50_standardized_coef'], 'T50_coef_standardized_t': info['age_T50_robust_t'],
        'T50_coef_standardized_p_HC3': info['age_T50_robust_p'],
        'T50_coef_raw_per_1sd': float(raw.params['sp_T50']), 'T50_coef_per_Gyr': float(raw.params['sp_T50']) / sd_t50,
        'T50_coef_per_Gyr_se': float(raw.bse['sp_T50']) / sd_t50, 'T50_coef_raw_t_HC3': float(raw.tvalues['sp_T50']),
        'T50_coef_raw_p_HC3': float(raw.pvalues['sp_T50']), 'model_R2': float(raw.rsquared),
        'cv_R2_controls_only': info['cv_R2_controls_only_mean'], 'cv_R2_controls_plus_age': info['cv_R2_controls_plus_age_mean'],
        'cv_delta_R2_from_age': info['cv_delta_R2_from_age'],
        'spearman_rho_resid_T50': info['resid_T50_spearman_rho'], 'spearman_p': info['resid_T50_spearman_p'],
    }
    if add_logda:
        # unstandardized OLS of the outcome on controls + T50 (HC3): size and distance coefficients
        X = sm.add_constant(base[controls + ['sp_T50']])
        r = sm.OLS(base['mfl_cyl_log_ML_dyn'], X).fit(cov_type='HC3')
        row.update({'coef_size_term_raw': float(r.params[size_col]), 'coef_logDA_raw': float(r.params['logDA']),
                    't_logDA_HC3': float(r.tvalues['logDA'])})
    per.insert(0, 'SPS_set', sps_label)
    per.insert(0, 'variant', key)
    return row, per, resid


def distance_by_age_quartile(df):
    """Per mass bin: median DA, Secondary fraction and SPS attenuation term of the young and old quartiles (published sample)."""
    controls = controls_for(['sp_ML_obs_Re'], 'logRe', False, 'strict')
    base = rc.make_base(df, controls)
    base['att'] = base['sp_ML_obs_Re'] - base['sp_ML_int_Re']   # log10(L_int / L_obs) = DML_int - DML_obs
    b = rc.add_mass_bins(base)
    da_arr = rc.get_bin_arrays(b, 'DA')
    sec_arr = rc.get_bin_arrays(b, 'is_secondary')
    m_arr = rc.get_bin_arrays(b, 'nsa_sersic_mass')
    att_arr = rc.get_bin_arrays(b, 'att')
    rows = []
    for a, s, m, t in zip(da_arr, sec_arr, m_arr, att_arr):
        dy, do = float(np.median(a['vals'][a['young_idx']])), float(np.median(a['vals'][a['old_idx']]))
        fy, fo = float(np.mean(s['vals'][s['young_idx']])), float(np.mean(s['vals'][s['old_idx']]))
        ty, to = float(np.median(t['vals'][t['young_idx']])), float(np.median(t['vals'][t['old_idx']]))
        rows.append({'mass_bin': a['b'], 'logM_median': a['logM_median'], 'N_young': a['k'], 'N_old': a['k'],
                     'median_DA_young_Mpc': dy, 'median_DA_old_Mpc': do, 'DA_ratio_young_over_old': dy / do,
                     'frac_secondary_young': fy, 'frac_secondary_old': fo, 'secondary_ratio_young_over_old': fy / fo,
                     'median_logM_young': float(np.median(m['vals'][m['young_idx']])),
                     'median_logM_old': float(np.median(m['vals'][m['old_idx']])),
                     'median_attenuation_young_dex': ty, 'median_attenuation_old_dex': to,
                     'attenuation_old_minus_young_dex': to - ty})
    out = pd.DataFrame(rows)
    out['abs_logM_median_diff'] = (out['median_logM_young'] - out['median_logM_old']).abs()
    facts = {
        'N': len(b), 'target_counts_in_sample': {int(k): int(v) for k, v in b['target'].value_counts().sort_index().items()},
        'spearman_logDA_vs_T50': float(spearmanr(b['logDA'], b['sp_T50'])[0]),
        'spearman_logRe_arcsec_vs_logRe_kpc': float(spearmanr(b['logRe'], b['logRe_kpc'])[0]),
        'median_Re_arcsec': float(10 ** b['logRe'].median()), 'median_Re_kpc': float(10 ** b['logRe_kpc'].median()),
        'DA_ratio_range': [float(out['DA_ratio_young_over_old'].min()), float(out['DA_ratio_young_over_old'].max())],
        'secondary_ratio_range': [float(out['secondary_ratio_young_over_old'].min()), float(out['secondary_ratio_young_over_old'].max())],
        'max_abs_logM_median_diff_young_old': float(out['abs_logM_median_diff'].max()),
        'attenuation_old_minus_young_by_bin': [float(x) for x in out['attenuation_old_minus_young_dex']],
    }
    return out, facts


def grid_by_size_unit(df):
    """The 12-variant reduced-spec grid (metallicity x flattening x SPS) under each size treatment: point estimates,
    plus the HC3 t and p of a linear T50 term added to the same controls."""
    rows = []
    for size_lab, size_col, add_logda in [('arcsec (published)', 'logRe', False), ('kpc', 'logRe_kpc', False),
                                          ('arcsec + logDA', 'logRe', True), ('kpc + logDA', 'logRe_kpc', True)]:
        for met in ['sp_MW_Metal_Re', 'sp_LW_Metal_Re']:
            for fl in ['Eps_MGE', 'nsa_sersic_ba']:
                for lab, sps in SPS_SETS:
                    c = sps + [met, 'nsa_sersic_mass', size_col, 'nsa_sersic_n', fl, 'Lambda_Re', 'logSigma_Re']
                    c += ['logDA'] if add_logda else []
                    base = rc.make_base(df, c)
                    resid, _ = rc.fit_controls_resid(base, c)
                    arr = rc.get_bin_arrays(rc.add_mass_bins(resid), 'resid')
                    diffs = np.array([np.median(a['vals'][a['young_idx']]) - np.median(a['vals'][a['old_idx']]) for a in arr])
                    raw = rc.age_coef_raw_outcome(base, c)
                    rows.append({'size': size_lab, 'metallicity': met, 'flattening': fl, 'SPS_set': lab, 'N': len(base),
                                 'young_gt_old_bins': int((diffs > 0).sum()), 'mean_diff': float(diffs.mean()),
                                 'T50_t_HC3': float(raw.tvalues['sp_T50']), 'T50_p_HC3': float(raw.pvalues['sp_T50'])})
    return pd.DataFrame(rows)


def point_estimates(d_binned, col):
    arr = rc.get_bin_arrays(d_binned, col)
    diffs = np.array([np.median(a['vals'][a['young_idx']]) - np.median(a['vals'][a['old_idx']]) for a in arr])
    py = np.concatenate([a['vals'][a['young_idx']] - a['bin_median'] for a in arr])
    po = np.concatenate([a['vals'][a['old_idx']] - a['bin_median'] for a in arr])
    return int((diffs > 0).sum()), int((diffs < 0).sum()), float(diffs.mean()), float(np.median(py) - np.median(po))


def one_cube_per_mangaid(df):
    """Every Table 1 row, full sample vs one entry per MaNGA ID (highest Qual, then highest SNR); the controlled
    residual with size in arcsec (published spec) and in kpc."""
    order = df.reset_index(drop=True).reset_index().rename(columns={'index': 'row'})
    dd = order.sort_values(['Qual', 'sp_SNR_Re', 'row'], ascending=[False, False, True], kind='mergesort')
    dd = dd.drop_duplicates('mangaid').sort_values('row').drop(columns='row').reset_index(drop=True)
    rows = []
    for lab, data in [('full sample', df), ('one cube per MaNGA ID', dd)]:
        q = data[(data['Qual'] >= 1) & np.isfinite(data[SAMPLE_REQUIRED]).all(axis=1)].reset_index(drop=True)
        b = rc.add_mass_bins(q)
        for col in ['DML_mfl_obs', 'DML_mfl_int', 'mfl_cyl_log_ML_dyn', 'sp_ML_int_Re', 'sp_ML_obs_Re', 'gnfw_cyl_fdm_Re', 'nfw_cyl_fdm_Re']:
            pos, neg, mean, madj = point_estimates(b, col)
            rows.append({'sample': lab, 'N': len(q), 'row': col, 'young_gt_old_bins': pos, 'old_gt_young_bins': neg,
                         'equal_bin_mean_diff': mean, 'mass_adjusted_diff': madj})
        for size_col, size_lab in [('logRe', ''), ('logRe_kpc', ' (kpc)')]:
            for sps_label, sps in SPS_SETS:
                c = controls_for(sps, size_col, False, 'strict')
                base = rc.make_base(data, c)
                resid, _ = rc.fit_controls_resid(base, c)
                pos, neg, mean, madj = point_estimates(rc.add_mass_bins(resid), 'resid')
                rows.append({'sample': lab, 'N': len(base), 'row': f'strict residual{size_lab}, {sps_label} SPS', 'young_gt_old_bins': pos,
                             'old_gt_young_bins': neg, 'equal_bin_mean_diff': mean, 'mass_adjusted_diff': madj})
    out = pd.DataFrame(rows)
    full = out[out['sample'] == 'full sample'].set_index('row')
    dedup = out[out['sample'] != 'full sample'].set_index('row')
    summary = {'N_full': int(full['N'].iloc[0]), 'N_one_cube': int(dedup['N'].iloc[0]),
               'bin_counts_unchanged': bool((full['young_gt_old_bins'] == dedup['young_gt_old_bins']).all()),
               'max_abs_change_mass_adjusted': float((full['mass_adjusted_diff'] - dedup['mass_adjusted_diff']).abs().max())}
    return out, summary


N_TIE_DRAWS = 200
TIE_SEED = rc.SEED + 9000
SIZE_TREATMENTS = [('arcsec (published)', 'logRe', False), ('kpc', 'logRe_kpc', False), ('arcsec + logDA', 'logRe', True)]


def tie_rule_redraws(df, n_draws=N_TIE_DRAWS, seed=TIE_SEED):
    """Sensitivity of the young/old quartile split to the T50 tie rule. T50 sits on a grid about 0.16 Gyr apart, so most
    quartile cuts fall inside a group of tied values; every script orders those ties by stellar mass (lower mass sorts as
    younger). Here the ties are broken at random instead, n_draws times per cell, with one numpy default_rng(seed)
    consumed in the output's row order. Cells: the 12-variant grid in the full sample under angular size, kpc size and
    angular size + logDA (kpc + logDA gives the same residuals as angular + logDA), and the strict specification
    (mass-weighted metallicity, MGE ellipticity) in the MaNGA Primary and Secondary samples. The fits are unchanged;
    only the quartile membership of tied galaxies moves, so the linear T50 term is unaffected. Point estimates only."""
    rng = np.random.default_rng(seed)
    cells = [('all', s, met, fl) for s in SIZE_TREATMENTS for met in ('sp_MW_Metal_Re', 'sp_LW_Metal_Re')
             for fl in ('Eps_MGE', 'nsa_sersic_ba')]
    cells += [(smp, s, 'sp_MW_Metal_Re', 'Eps_MGE') for smp in ('primary', 'secondary') for s in SIZE_TREATMENTS]
    rows = []
    for sample, (size_lab, size_col, add_logda), met, fl in cells:
        d0 = df if sample == 'all' else df[df['target'] == TARGET_OF[sample]]
        for lab, sps in SPS_SETS:
            c = sps + [met, 'nsa_sersic_mass', size_col, 'nsa_sersic_n', fl, 'Lambda_Re', 'logSigma_Re']
            c += ['logDA'] if add_logda else []
            base = rc.make_base(d0, c)
            resid, _ = rc.fit_controls_resid(base, c)
            binned = rc.add_mass_bins(resid)
            pos, _, _, madj = point_estimates(binned, 'resid')
            arr = rc.get_bin_arrays(binned, 'resid')
            counts, ma = np.empty(n_draws, dtype=int), np.empty(n_draws)
            for r in range(n_draws):
                cnt, py, po = 0, [], []
                for a in arr:
                    order = np.lexsort((rng.random(len(a['vals'])), a['ages']))
                    y, o = a['vals'][order[:a['k']]], a['vals'][order[-a['k']:]]
                    cnt += int(np.median(y) - np.median(o) > 0)
                    py.append(y - a['bin_median'])
                    po.append(o - a['bin_median'])
                counts[r] = cnt
                ma[r] = np.median(np.concatenate(py)) - np.median(np.concatenate(po))
            v, f = np.unique(counts, return_counts=True)
            rows.append({'sample': sample, 'size': size_lab, 'metallicity': met, 'flattening': fl, 'SPS_set': lab, 'N': len(base),
                         'bins_mass_tie_rule': pos, 'mass_adjusted_diff_mass_tie_rule': madj, 'n_draws': n_draws,
                         'bins_min': int(counts.min()), 'bins_median': float(np.median(counts)), 'bins_max': int(counts.max()),
                         'bins_histogram': ';'.join(f'{int(a)}:{int(b)}' for a, b in zip(v, f)),
                         'mass_adjusted_min': float(ma.min()), 'mass_adjusted_median': float(np.median(ma)),
                         'mass_adjusted_max': float(ma.max())})
    return pd.DataFrame(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fmt_row(r):
    ci = f"[{r['mass_adjusted_ci_lo']:+.4f}, {r['mass_adjusted_ci_hi']:+.4f}]"
    return (f"| {r['variant']} | {r['SPS_set']} | {r['N']} | {r['young_gt_old_bins']}/8 | {r['mass_adjusted_diff']:+.4f} | {ci} | "
            f"{r['p_count']:.4f} | {r['p_magnitude']:.4f} | {r['T50_coef_standardized']:.4f} ({r['T50_coef_standardized_p_HC3']:.2g}) | "
            f"{r['T50_coef_per_Gyr']:.5f} (t {r['T50_coef_raw_t_HC3']:.2f}) | {r['spearman_rho_resid_T50']:.4f} |")


def main(argv=None):
    ap = argparse.ArgumentParser(description='Paper 10 size/distance/subsample sensitivity of the controlled residual.')
    ap.add_argument('--input', default=str(DEFAULT_IN))
    ap.add_argument('--da', default=str(DEFAULT_DA))
    ap.add_argument('--outdir', default=str(HERE / 'regenerated'))  # outputs/ holds the shipped reference run
    args = ap.parse_args(argv)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df = load(args.input, args.da)

    rows, pers, resids = [], [], {}
    for key, desc, sample, size_col, add_logda, family in VARIANTS:
        for i, (sps_label, sps) in enumerate(SPS_SETS):
            print('running', key, sps_label, flush=True)
            row, per, resid = run_variant(df, key, desc, sample, size_col, add_logda, family, i, sps_label, sps)
            rows.append(row)
            pers.append(per)
            resids[(key, sps_label)] = resid['resid'].to_numpy()
    board = pd.DataFrame(rows)
    perbin = pd.concat(pers, ignore_index=True)
    c_vs_d = max(float(np.max(np.abs(resids[('c_published_plus_logDA', s)] - resids[('d_kpc_plus_logDA', s)]))) for s, _ in SPS_SETS)

    dist, facts = distance_by_age_quartile(df)
    grid = grid_by_size_unit(df)
    gsum = (grid.groupby(['size', 'SPS_set']).agg(bins_min=('young_gt_old_bins', 'min'), bins_max=('young_gt_old_bins', 'max'),
                                                   T50_t_min=('T50_t_HC3', 'min'), T50_t_max=('T50_t_HC3', 'max')).reset_index())
    grid_counts = {s: {int(k): int(v) for k, v in grid[grid['size'] == s]['young_gt_old_bins'].value_counts().sort_index().items()}
                   for s in grid['size'].unique()}
    dedup, dedup_sum = one_cube_per_mangaid(df)
    print('running tie-rule redraws', flush=True)
    ties = tie_rule_redraws(df)

    kw = dict(index=False, lineterminator='\n')
    board.to_csv(out / 'sensitivity_scoreboard.csv', **kw)
    perbin.to_csv(out / 'sensitivity_perbin.csv', **kw)
    dist.to_csv(out / 'distance_by_age_quartile.csv', **kw)
    grid.to_csv(out / 'reduced_spec_grid_by_size_unit.csv', **kw)
    dedup.to_csv(out / 'one_cube_per_mangaid.csv', **kw)
    ties.to_csv(out / 'tie_rule_redraws.csv', **kw)
    summary = {
        'script': Path(__file__).name, 'input': Path(args.input).name, 'input_sha256': sha256(args.input),
        'da_input': Path(args.da).name, 'da_input_sha256': sha256(args.da),
        'seed_rule': 'SPS set i (0 int, 1 obs, 2 both) uses SEED + 1000*i in every variant (SEED = 20260629)',
        'n_boot': rc.N_BOOT, 'n_shuffle': rc.N_SHUFFLE, 'arcsec_per_radian': ARCSEC_PER_RAD,
        'max_abs_resid_difference_c_vs_d': c_vs_d, 'sample_facts': facts,
        'grid_bin_count_histogram_by_size': grid_counts,
        'grid_T50_t_HC3_range': [float(grid['T50_t_HC3'].min()), float(grid['T50_t_HC3'].max())],
        'grid_T50_p_HC3_max': float(grid['T50_p_HC3'].max()),
        'one_cube_per_mangaid': dedup_sum,
        'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
        'tie_rule_redraws': {'n_draws': N_TIE_DRAWS, 'seed': TIE_SEED, 'cells': int(len(ties)),
                             'rule': 'T50 ties at the quartile cuts broken at random instead of by stellar mass'},
    }
    with open(out / 'sensitivity_summary.json', 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(summary, fh, indent=2)
    md = ['# Size / distance / subsample sensitivity (generated by run_size_distance_sensitivity.py)', '',
          'Bins = mass bins with young > old controlled residual; mass-adjusted difference and its bootstrap 95% CI; '
          'count and magnitude shuffle p (5,000 shuffles, floor 0.0002); standardized T50 coefficient (HC3 p); '
          'T50 coefficient on the raw outcome in dex per Gyr (HC3 t); Spearman rho of residual vs T50.', '',
          '| Variant | SPS | N | Bins | Mass-adj. | 95% CI | Count p | Mag. p | Std. T50 coef (p) | dex/Gyr (t) | rho |',
          '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    md += [fmt_row(r) for _, r in board.iterrows()]
    md += ['', '## Young vs old quartile distance and MaNGA subsample, per mass bin', '', dist.to_markdown(index=False, floatfmt='.4g'),
           '', '## 12-variant reduced-spec grid: bins with young > old and HC3 t of a linear T50 term, by size treatment', '', gsum.to_markdown(index=False),
           '', '## One cube per MaNGA ID (point estimates)', '', dedup.to_markdown(index=False, floatfmt='.4g'), '',
           f'## T50 tie rule: young > old bins with ties at the quartile cuts broken at random ({N_TIE_DRAWS} draws per cell, '
           f'seed {TIE_SEED}; every other output breaks them by stellar mass)', '',
           ties[['sample', 'size', 'metallicity', 'flattening', 'SPS_set', 'bins_mass_tie_rule', 'bins_histogram',
                 'mass_adjusted_diff_mass_tie_rule', 'mass_adjusted_min', 'mass_adjusted_max']].to_markdown(index=False, floatfmt='.4f'), '',
           f'Sample facts: {json.dumps(facts)}', '', f'Max |residual(c) - residual(d)| = {c_vs_d:.3g}', '']
    with open(out / 'sensitivity_summary.md', 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(md))
    print('\n'.join(md[4:4 + 2 + len(board)]))
    print(json.dumps(facts))
    print(json.dumps(dedup_sum))
    print('grid', json.dumps(grid_counts))


if __name__ == '__main__':
    main()
