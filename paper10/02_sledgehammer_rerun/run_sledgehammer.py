"""Paper 10 sledgehammer rerun: raw DML split, M/L component split, secondary fDM checks,
nearest-mass paired check, sample counts and T50 grid facts.

Re-implementation written in October 2026. The code that produced the v0.6 outputs in
../01_clean_pipeline_v0_6/ was not preserved; only its reports and CSVs were. This script
rebuilds those analyses from the method text (appendix Sections 4-7) and follows the
conventions of eval_sledge() in ../00_strict_control_rerun/run_clean_controls_fast.py:

  * sample: Qual >= 1 and finite sp_T50, nsa_sersic_mass, DML_mfl_int, DML_mfl_obs
    (N = 5,952 PlateIFU entries of 5,875 distinct MaNGA IDs; repeat observations kept, as in v0.6);
  * eight equal-count bins in nsa_sersic_mass (744 each), ties in mass broken by row order;
  * within each bin the 186 youngest and 186 oldest by sp_T50, T50 ties broken by row order,
    i.e. by stellar mass (np.lexsort);
  * per-bin young-minus-old difference of medians; 5,000 bootstraps per bin (young and old
    resampled independently);
  * mass-adjusted difference: subtract each bin's median (all 744 galaxies) from its young and
    old values, pool the eight bins, take young median minus old median; 5,000 bootstraps of
    the pooled groups;
  * 5,000 within-bin shuffles, +1 correction (floor 1/5,001 = 0.0002): count p (bins in the
    tested direction >= real), magnitude p (sum of the eight bin differences >= real, or <= real
    for the old > young direction) and joint p.

Null (--null). Default value_permutation: within each bin the outcome values are permuted over
the fixed young/old positions, so the two groups are exchangeable. This is equivalent in
distribution to shuffling T50 with random tie-breaking, and it is the null whose shuffle
statistics (mean shuffled bin count) match the archived v0.6 outputs; the v0.6 code itself was
not preserved, so this is inferred. The alternative, age_permutation, is the null used by
run_clean_controls_fast.py for the controlled residual: T50 is permuted and the quartiles are
rebuilt with the same mass tie-break as the real split. Because T50 sits on a 0.16 Gyr grid,
that tie rule nudges young toward lower and old toward higher mass inside a bin in the null
too. The two nulls give the same conclusions; count p differs by up to about 0.02.

Random numbers (numpy default_rng): seed base 20260629. Outcome i (order of OUTCOMES) uses
base + 1000*i; paired outcome j uses base + 100000 + 1000*j. Within an outcome the call order
is: per-bin bootstraps, pooled bootstraps, shuffles (identical to eval_sledge), so bootstrap
intervals do not depend on the null. Deterministic quantities (counts, medians, mass-adjusted
differences, pair wins) reproduce v0.6 exactly; shuffle p-values and bootstrap intervals agree
with v0.6 to Monte Carlo precision (the v0.6 seed is unknown). sledgehammer_vs_v06.csv records
the archived and rerun values side by side.

Nearest-mass paired check: within each mass bin the 186 young and 186 old galaxies are each
sorted by stellar mass (stable, so equal masses keep T50 order) and paired rank by rank. This
rule reproduces the v0.6 counts exactly (881/1,488 and 913/1,488); greedy, with-replacement and
Hungarian matchings do not.

Usage:  python run_sledgehammer.py [--input CSV] [--outdir DIR] [--null value_permutation|age_permutation]
"""
import argparse
import hashlib
import json
import platform
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
DEFAULT_IN = PKG / 'data' / 'manga_dynpop_merged_thin_firstpass.csv'
V06 = PKG / '01_clean_pipeline_v0_6'

SEED = 20260629
N_BOOT = 5000
N_SHUFFLE = 5000
N_BINS = 8
K = 186
AGE = 'sp_T50'
MASS = 'nsa_sersic_mass'
SAMPLE_REQUIRED = [AGE, MASS, 'DML_mfl_int', 'DML_mfl_obs']
T50_MODE = 12.848485       # most common T50 value (Gyr), on the catalogue's ~0.16 Gyr grid

# (column, label, family, direction tested for the directional p)
OUTCOMES = [
    ('DML_mfl_int', 'DML intrinsic SPS', 'raw_DML', 'young_gt_old'),
    ('DML_mfl_obs', 'DML observed SPS', 'raw_DML', 'young_gt_old'),
    ('mfl_cyl_log_ML_dyn', 'JAM dynamical log(M/L)', 'component', 'young_gt_old'),
    ('sp_ML_int_Re', 'SPS intrinsic log(M*/L)', 'component', 'old_gt_young'),
    ('sp_ML_obs_Re', 'SPS observed log(M*/L)', 'component', 'old_gt_young'),
    ('gnfw_cyl_fdm_Re', 'gNFW fDM(<Re)', 'fDM', 'young_gt_old'),
    ('nfw_cyl_fdm_Re', 'NFW fDM(<Re)', 'fDM', 'young_gt_old'),
]
PAIRED_OUTCOMES = ['DML_mfl_int', 'DML_mfl_obs']


def make_sample(df):
    counts = {'joined_rows': len(df), 'joined_distinct_mangaid': int(df['mangaid'].nunique())}
    q = df[df['Qual'] >= 1]
    counts['qual_ge_1'] = len(q)
    d = q[np.isfinite(q[SAMPLE_REQUIRED]).all(axis=1)].copy().reset_index(drop=True)
    vc = d['mangaid'].value_counts()
    counts['primary_finite_sample'] = len(d)
    counts['primary_sample_distinct_mangaid'] = int(d['mangaid'].nunique())
    counts['primary_sample_mangaid_repeated'] = int((vc > 1).sum())
    counts['primary_sample_rows_of_repeated_mangaid'] = int(vc[vc > 1].sum())
    return d, counts


def add_mass_bins(d):
    d = d.copy().reset_index(drop=False).rename(columns={'index': 'orig_order'})
    d = d.sort_values([MASS, 'orig_order'], kind='mergesort').reset_index(drop=True)
    per = len(d) // N_BINS
    assert per * N_BINS == len(d), len(d)
    d['mass_bin'] = np.repeat(np.arange(N_BINS), per)
    return d


def split_young_old(ages):
    """(young_idx, old_idx) within one bin; rows are in ascending-mass order, so T50 ties go by mass."""
    order = np.lexsort((np.arange(len(ages)), ages))
    return order[:K], order[-K:]


def bin_arrays(d, outcome):
    arrays = []
    for b in range(N_BINS):
        g = d[d['mass_bin'] == b].reset_index(drop=True)
        ages = g[AGE].to_numpy(float)
        yi, oi = split_young_old(ages)
        arrays.append({
            'b': b, 'vals': g[outcome].to_numpy(float), 'ages': ages, 'mass': g[MASS].to_numpy(float),
            'young_idx': yi, 'old_idx': oi,
            'logM_min': float(g[MASS].min()), 'logM_max': float(g[MASS].max()),
            'logM_median': float(g[MASS].median()),
            'T50_young_median': float(np.median(ages[yi])), 'T50_old_median': float(np.median(ages[oi])),
            'bin_median': float(np.median(g[outcome])),
        })
    return arrays


def _p_ge(null, real, n):
    return float((np.sum(null >= real) + 1) / (n + 1))


def _p_le(null, real, n):
    return float((np.sum(null <= real) + 1) / (n + 1))


def evaluate_outcome(d, outcome, label, family, direction, seed, n_boot, n_shuffle, null):
    rng = np.random.default_rng(seed)
    arrays = bin_arrays(d, outcome)
    per_rows = []
    boot_bin_diffs = np.zeros((n_boot, N_BINS))
    pooled_y, pooled_o = [], []
    for j, a in enumerate(arrays):
        y = a['vals'][a['young_idx']]
        o = a['vals'][a['old_idx']]
        diff = float(np.median(y) - np.median(o))
        y_samp = y[rng.integers(0, len(y), size=(n_boot, len(y)))]
        o_samp = o[rng.integers(0, len(o), size=(n_boot, len(o)))]
        bdiff = np.median(y_samp, axis=1) - np.median(o_samp, axis=1)
        boot_bin_diffs[:, j] = bdiff
        lo, hi = np.percentile(bdiff, [2.5, 97.5])
        per_rows.append({
            'outcome': outcome, 'outcome_label': label, 'family': family,
            'mass_bin': a['b'], 'N_bin': len(a['vals']), 'N_young': len(y), 'N_old': len(o),
            'logM_min': a['logM_min'], 'logM_max': a['logM_max'], 'logM_median': a['logM_median'],
            'T50_young_median': a['T50_young_median'], 'T50_old_median': a['T50_old_median'],
            'median_young': float(np.median(y)), 'median_old': float(np.median(o)),
            'median_diff_young_minus_old': diff, 'bootstrap95_lo': float(lo), 'bootstrap95_hi': float(hi),
            'young_gt_old': bool(diff > 0), 'old_gt_young': bool(diff < 0),
        })
        pooled_y.append(y - a['bin_median'])
        pooled_o.append(o - a['bin_median'])
    per = pd.DataFrame(per_rows)
    diffs = per['median_diff_young_minus_old'].to_numpy()
    eq_lo, eq_hi = np.percentile(boot_bin_diffs.mean(axis=1), [2.5, 97.5])

    py = np.concatenate(pooled_y)
    po = np.concatenate(pooled_o)
    madiff = float(np.median(py) - np.median(po))
    my_boot = np.median(py[rng.integers(0, len(py), size=(n_boot, len(py)))], axis=1)
    mo_boot = np.median(po[rng.integers(0, len(po), size=(n_boot, len(po)))], axis=1)
    ma_boot = my_boot - mo_boot
    malo, mahi = np.percentile(ma_boot, [2.5, 97.5])

    count_pos = np.empty(n_shuffle, dtype=int)
    count_neg = np.empty(n_shuffle, dtype=int)
    sums = np.empty(n_shuffle)
    for s in range(n_shuffle):
        cp = cn = 0
        ssum = 0.0
        for a in arrays:
            vals = a['vals']
            if null == 'age_permutation':
                order = np.lexsort((np.arange(len(vals)), rng.permutation(a['ages'])))
                diff = float(np.median(vals[order[:K]]) - np.median(vals[order[-K:]]))
            else:
                pv = rng.permutation(vals)
                diff = float(np.median(pv[a['young_idx']]) - np.median(pv[a['old_idx']]))
            cp += diff > 0
            cn += diff < 0
            ssum += diff
        count_pos[s] = cp
        count_neg[s] = cn
        sums[s] = ssum

    real_pos = int((diffs > 0).sum())
    real_neg = int((diffs < 0).sum())
    real_sum = float(diffs.sum())
    p_pos = (_p_ge(count_pos, real_pos, n_shuffle), _p_ge(sums, real_sum, n_shuffle),
             float((np.sum((count_pos >= real_pos) & (sums >= real_sum)) + 1) / (n_shuffle + 1)))
    p_neg = (_p_ge(count_neg, real_neg, n_shuffle), _p_le(sums, real_sum, n_shuffle),
             float((np.sum((count_neg >= real_neg) & (sums <= real_sum)) + 1) / (n_shuffle + 1)))
    dir_count, dir_p = (real_pos, p_pos) if direction == 'young_gt_old' else (real_neg, p_neg)
    score = {
        'outcome': outcome, 'outcome_label': label, 'family': family, 'seed': seed,
        'N': len(d), 'n_boot': n_boot, 'n_shuffle': n_shuffle, 'null': null, 'direction_tested': direction,
        'young_gt_old_bins': real_pos, 'old_gt_young_bins': real_neg, 'directional_bins': dir_count,
        'total_bins': N_BINS,
        'equal_bin_mean_diff_young_minus_old': float(diffs.mean()),
        'equal_bin_mean_boot95_lo': float(eq_lo), 'equal_bin_mean_boot95_hi': float(eq_hi),
        'sum_of_8_bin_diffs': real_sum,
        'mass_adjusted_young_median': float(np.median(py)),
        'mass_adjusted_old_median': float(np.median(po)),
        'mass_adjusted_diff_young_minus_old': madiff,
        'mass_adjusted_diff_boot95_lo': float(malo), 'mass_adjusted_diff_boot95_hi': float(mahi),
        'p_young_gt_old_count': p_pos[0], 'p_young_gt_old_sum': p_pos[1], 'p_young_gt_old_joint': p_pos[2],
        'p_old_gt_young_count': p_neg[0], 'p_old_gt_young_sum': p_neg[1], 'p_old_gt_young_joint': p_neg[2],
        'p_directional_count': dir_p[0], 'p_directional_sum': dir_p[1], 'p_directional_joint': dir_p[2],
        # sign-test reference: P(>= real_pos of 8) for a fair coin per bin (NOT the shuffle p)
        'binomial_reference_p_young_gt_old_count': sum(comb(N_BINS, k) for k in range(real_pos, N_BINS + 1)) / 2 ** N_BINS,
        'shuffle_count_young_gt_old_mean': float(count_pos.mean()),
        'shuffle_count_young_gt_old_sd': float(count_pos.std(ddof=1)),
        'shuffle_sum_mean': float(sums.mean()), 'shuffle_sum_sd': float(sums.std(ddof=1)),
    }
    return score, per


def _pair_wins(vals, mass, yi, oi):
    ys = yi[np.argsort(mass[yi], kind='stable')]
    os_ = oi[np.argsort(mass[oi], kind='stable')]
    return int(np.sum(vals[ys] > vals[os_]))


def paired_check(d, outcome, label, seed, n_shuffle, null):
    rng = np.random.default_rng(seed)
    arrays = bin_arrays(d, outcome)
    n_pairs = N_BINS * K
    wins = sum(_pair_wins(a['vals'], a['mass'], a['young_idx'], a['old_idx']) for a in arrays)
    frac = wins / n_pairs
    nullv = np.empty(n_shuffle)
    for s in range(n_shuffle):
        w = 0
        for a in arrays:
            if null == 'age_permutation':
                yi, oi = split_young_old(rng.permutation(a['ages']))
                w += _pair_wins(a['vals'], a['mass'], yi, oi)
            else:
                w += _pair_wins(rng.permutation(a['vals']), a['mass'], a['young_idx'], a['old_idx'])
        nullv[s] = w / n_pairs
    return {
        'outcome': outcome, 'outcome_label': label, 'seed': seed, 'null': null,
        'n_pairs': n_pairs, 'young_wins': wins, 'young_win_fraction': frac,
        'binomial_p_greater_0p5': float(binomtest(wins, n_pairs, 0.5, alternative='greater').pvalue),
        'shuffle_p_win_fraction_ge_real': _p_ge(nullv, frac, n_shuffle),
        'shuffle_win_fraction_mean': float(nullv.mean()), 'shuffle_win_fraction_sd': float(nullv.std(ddof=1)),
    }


def t50_grid_check(d):
    rows = []
    cuts_inside_tie_groups = 0
    for b in range(N_BINS):
        g = d[d['mass_bin'] == b].reset_index(drop=True)
        ages = g[AGE].to_numpy(float)
        yi, oi = split_young_old(ages)
        young, old = ages[yi], ages[oi]
        s = np.sort(ages)
        young_cut_tied = bool(s[K - 1] == s[K])              # youngest-quartile boundary inside a tie group
        old_cut_tied = bool(s[-K] == s[-K - 1])              # oldest-quartile boundary inside a tie group
        cuts_inside_tie_groups += young_cut_tied + old_cut_tied
        rows.append({
            'mass_bin': b,
            'T50_young_median': float(np.median(young)), 'T50_young_min': float(young.min()),
            'T50_young_max': float(young.max()),
            'T50_young_iqr': float(np.percentile(young, 75) - np.percentile(young, 25)),
            'T50_old_median': float(np.median(old)), 'T50_old_min': float(old.min()), 'T50_old_max': float(old.max()),
            'T50_old_iqr': float(np.percentile(old, 75) - np.percentile(old, 25)),
            'n_old_at_mode': int(np.isclose(old, T50_MODE, atol=1e-4).sum()),
            'n_old_above_mode': int((old > T50_MODE + 1e-4).sum()),
            'young_cut_inside_tie_group': young_cut_tied, 'old_cut_inside_tie_group': old_cut_tied,
            'bin_T50_max': float(ages.max()),
        })
    t = d[AGE].to_numpy(float)
    grid, cnt = np.unique(np.round(t, 6), return_counts=True)
    per = pd.DataFrame(rows)
    summary = {
        'sample_T50_min': float(grid.min()), 'sample_T50_max': float(grid.max()),
        'n_distinct_T50_values': int(len(grid)), 'typical_grid_step': float(np.median(np.diff(grid))),
        'T50_mode': float(grid[np.argmax(cnt)]), 'n_sample_at_mode': int(cnt.max()),
        'n_sample_above_mode': int((t > T50_MODE + 1e-4).sum()),
        'bins_with_old_median_at_mode': int(np.isclose(per['T50_old_median'], T50_MODE, atol=1e-4).sum()),
        'old_quartile_T50_min_over_bins': float(per['T50_old_min'].min()),
        'young_quartile_T50_max_over_bins': float(per['T50_young_max'].max()),
        'old_quartile_iqr_range': [float(per['T50_old_iqr'].min()), float(per['T50_old_iqr'].max())],
        'young_quartile_iqr_range': [float(per['T50_young_iqr'].min()), float(per['T50_young_iqr'].max())],
        'old_quartile_n_above_mode_range': [int(per['n_old_above_mode'].min()), int(per['n_old_above_mode'].max())],
        'quartile_cuts_inside_tie_groups': int(cuts_inside_tie_groups), 'quartile_cuts_total': 2 * N_BINS,
    }
    return per, summary


def compare_to_v06(scoreboard, per_bin, paired):
    """Archived v0.6 values next to this rerun's values (deterministic: should agree exactly)."""
    rows = []

    def add(quantity, v06, rerun, kind):
        rows.append({'quantity': quantity, 'v06_archived': v06, 'rerun': rerun, 'abs_diff': abs(v06 - rerun), 'kind': kind})

    sc = scoreboard.set_index('outcome')
    dml = pd.read_csv(V06 / 'manga_sledgehammer_final_scoreboard.csv').set_index('outcome')
    for o in ['DML_mfl_int', 'DML_mfl_obs']:
        r, v = sc.loc[o], dml.loc[o]
        add(f'{o} young>old bins', v['real_EP_favorable_bins'], r['young_gt_old_bins'], 'deterministic')
        add(f'{o} mean bin diff', v['real_mean_bin_median_diff_young_minus_old'], r['equal_bin_mean_diff_young_minus_old'], 'deterministic')
        add(f'{o} mass-adjusted diff', v['mass_adjusted_median_diff_young_minus_old'], r['mass_adjusted_diff_young_minus_old'], 'deterministic')
        add(f'{o} mass-adjusted CI lo', v['diff_bootstrap95_lo'], r['mass_adjusted_diff_boot95_lo'], 'Monte Carlo')
        add(f'{o} mass-adjusted CI hi', v['diff_bootstrap95_hi'], r['mass_adjusted_diff_boot95_hi'], 'Monte Carlo')
        add(f'{o} count p', v['shuffle_p_count_ge_real'], r['p_young_gt_old_count'], 'Monte Carlo')
        add(f'{o} magnitude p', v['shuffle_p_mean_diff_ge_real'], r['p_young_gt_old_sum'], 'Monte Carlo')
        add(f'{o} shuffled count mean', v['shuffle_count_mean'], r['shuffle_count_young_gt_old_mean'], 'Monte Carlo')
        pr = paired.set_index('outcome').loc[o]
        add(f'{o} paired young wins', v['young_wins'], pr['young_wins'], 'deterministic')
        add(f'{o} paired binomial p', v['binomial_p_greater_0p5'], pr['binomial_p_greater_0p5'], 'deterministic')
        add(f'{o} paired shuffle p', v['shuffle_p_win_fraction_ge_real'], pr['shuffle_p_win_fraction_ge_real'], 'Monte Carlo')
    comp = pd.read_csv(V06 / 'manga_ml_split_hardening_scoreboard.csv').set_index('outcome')
    for o in ['mfl_cyl_log_ML_dyn', 'sp_ML_int_Re', 'sp_ML_obs_Re']:
        r, v = sc.loc[o], comp.loc[o]
        add(f'{o} young>old bins', v['young_gt_old_bins'], r['young_gt_old_bins'], 'deterministic')
        add(f'{o} old>young bins', v['old_gt_young_bins'], r['old_gt_young_bins'], 'deterministic')
        add(f'{o} mean bin diff', v['equal_bin_mean_diff_young_minus_old'], r['equal_bin_mean_diff_young_minus_old'], 'deterministic')
        add(f'{o} mass-adjusted diff', v['mass_adjusted_diff_young_minus_old'], r['mass_adjusted_diff_young_minus_old'], 'deterministic')
        add(f'{o} mass-adjusted CI lo', v['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_lo'], 'Monte Carlo')
        add(f'{o} mass-adjusted CI hi', v['mass_adjusted_diff_boot95_hi'], r['mass_adjusted_diff_boot95_hi'], 'Monte Carlo')
        add(f'{o} directional count p', v['p_directional_count'], r['p_directional_count'], 'Monte Carlo')
        add(f'{o} directional magnitude p', v['p_directional_sumdiff'], r['p_directional_sum'], 'Monte Carlo')
    fdm = pd.read_csv(V06 / 'manga_sledgehammer_final_secondary_fdm.csv').set_index('outcome')
    for o in ['gnfw_cyl_fdm_Re', 'nfw_cyl_fdm_Re']:
        r, v = sc.loc[o], fdm.loc[o]
        add(f'{o} young>old bins', v['real_EP_favorable_bins'], r['young_gt_old_bins'], 'deterministic')
        add(f'{o} mean bin diff', v['real_mean_bin_median_diff_young_minus_old'], r['equal_bin_mean_diff_young_minus_old'], 'deterministic')
        add(f'{o} count p', v['shuffle_p_count_ge_real'], r['p_young_gt_old_count'], 'Monte Carlo')
        add(f'{o} magnitude p', v['shuffle_p_mean_diff_ge_real'], r['p_young_gt_old_sum'], 'Monte Carlo')
        add(f'{o} joint p', v['shuffle_p_joint_count_and_mean_ge_real'], r['p_young_gt_old_joint'], 'Monte Carlo')
        add(f'{o} shuffled count mean', v['shuffle_count_mean'], r['shuffle_count_young_gt_old_mean'], 'Monte Carlo')
    # per-bin differences (deterministic) from the two archived per-bin tables
    pb = per_bin.set_index(['outcome', 'mass_bin'])['median_diff_young_minus_old']
    b1 = pd.read_csv(V06 / 'manga_sledgehammer_final_bin_results.csv')
    b2 = pd.read_csv(V06 / 'manga_ml_split_hardening_per_bin.csv')
    worst = 0.0
    for _, r in b1.iterrows():
        worst = max(worst, abs(r['median_diff_young_minus_old'] - pb.loc[(r['outcome'], r['mass_bin'])]))
    for _, r in b2[b2['outcome'].isin(['mfl_cyl_log_ML_dyn', 'sp_ML_int_Re', 'sp_ML_obs_Re'])].iterrows():
        worst = max(worst, abs(r['median_diff_young_minus_old'] - pb.loc[(r['outcome'], r['mass_bin'])]))
    add('max |per-bin young-old diff| vs v0.6 per-bin tables (DML and components)', 0.0, worst, 'deterministic')
    out = pd.DataFrame(rows)
    # Monte Carlo SE of a p-value estimated from n = 5,000 shuffles, for p rows
    is_p = out['quantity'].str.endswith(' p') & (out['kind'] == 'Monte Carlo')
    out['mc_se_of_one_estimate'] = np.nan
    pv = out.loc[is_p, 'v06_archived'].to_numpy(float)
    out.loc[is_p, 'mc_se_of_one_estimate'] = np.sqrt(pv * (1 - pv) / N_SHUFFLE)
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description='Paper 10 raw DML / component / fDM / paired rerun.')
    ap.add_argument('--input', default=str(DEFAULT_IN))
    ap.add_argument('--outdir', default=str(HERE / 'regenerated'))  # outputs/ holds the shipped reference run
    ap.add_argument('--seed', type=int, default=SEED)
    ap.add_argument('--n-boot', type=int, default=N_BOOT)
    ap.add_argument('--n-shuffle', type=int, default=N_SHUFFLE)
    ap.add_argument('--null', choices=['value_permutation', 'age_permutation'], default='value_permutation')
    args = ap.parse_args(argv)

    inp, out = Path(args.input), Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(inp)
    base, counts = make_sample(df)
    d = add_mass_bins(base)
    counts['galaxies_per_mass_bin'] = len(d) // N_BINS
    counts['galaxies_per_quartile'] = K

    scores, pers = [], []
    for i, (col, label, family, direction) in enumerate(OUTCOMES):
        print('running', col, flush=True)
        sc, per = evaluate_outcome(d, col, label, family, direction, args.seed + 1000 * i,
                                   args.n_boot, args.n_shuffle, args.null)
        scores.append(sc)
        pers.append(per)
    scoreboard = pd.DataFrame(scores)
    per_bin = pd.concat(pers, ignore_index=True)
    paired = []
    for j, col in enumerate(PAIRED_OUTCOMES):
        print('paired', col, flush=True)
        label = [o[1] for o in OUTCOMES if o[0] == col][0]
        paired.append(paired_check(d, col, label, args.seed + 100000 + 1000 * j, args.n_shuffle, args.null))
    paired = pd.DataFrame(paired)
    grid_bins, grid_summary = t50_grid_check(d)

    kw = dict(index=False, lineterminator='\n')
    scoreboard.to_csv(out / 'sledgehammer_rerun_scoreboard.csv', **kw)
    per_bin.to_csv(out / 'sledgehammer_rerun_bin_results.csv', **kw)
    paired.to_csv(out / 'sledgehammer_rerun_paired_summary.csv', **kw)
    grid_bins.to_csv(out / 'sledgehammer_rerun_t50_grid_check.csv', **kw)
    pd.DataFrame([counts]).to_csv(out / 'sledgehammer_rerun_sample_counts.csv', **kw)
    cmp = compare_to_v06(scoreboard, per_bin, paired)
    cmp.to_csv(out / 'sledgehammer_vs_v06.csv', **kw)
    manifest = {
        'script': Path(__file__).name, 'input': inp.name, 'input_sha256': sha256(inp),
        'input_bytes': inp.stat().st_size, 'seed_base': args.seed, 'null': args.null,
        'seed_rule': 'outcome i uses seed_base + 1000*i (order of OUTCOMES); paired outcome j uses seed_base + 100000 + 1000*j',
        'n_boot': args.n_boot, 'n_shuffle': args.n_shuffle, 'n_bins': N_BINS, 'k_per_quartile': K,
        'sample_counts': counts, 't50_grid_summary': grid_summary,
        'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
    }
    with open(out / 'sledgehammer_rerun_manifest.json', 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(manifest, fh, indent=2)
    show = ['outcome', 'young_gt_old_bins', 'old_gt_young_bins', 'equal_bin_mean_diff_young_minus_old',
            'mass_adjusted_diff_young_minus_old', 'mass_adjusted_diff_boot95_lo', 'mass_adjusted_diff_boot95_hi',
            'p_young_gt_old_count', 'p_young_gt_old_sum', 'p_directional_count', 'p_directional_sum']
    print(scoreboard[show].to_string(index=False))
    print(paired.to_string(index=False))
    print(json.dumps(counts))
    print(json.dumps(grid_summary))


if __name__ == '__main__':
    main()
