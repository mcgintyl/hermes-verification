"""Paper 10 one-command reproduction and number check.

    python reproduce.py                      # from the shipped joined table in data/  (about 3 minutes)
    python reproduce.py --from-public --jam SDSSDR17_MaNGA_JAM_v2.fits \
        --sfh DynPop2_SP_SFH_v2.hdf5.zip --cvc SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt
                                             # rebuild the table from the Zenodo files first

Steps: (1) check the input tables (SHA-256 of the shipped files, or rebuild them from the public
files and compare); (2) run 00_strict_control_rerun/run_clean_controls_fast.py,
02_sledgehammer_rerun/run_sledgehammer.py and 03_size_distance_sensitivity/
run_size_distance_sensitivity.py into --outdir (default reproduce_outputs/); (3) compare the
regenerated files with the archived and reference outputs; (4) check every analysis number quoted
in the v3 main paper and appendix, at the precision it is printed, and the canonical size/distance
sensitivity numbers, against the regenerated values. Literature values the paper quotes (template
ages, survey thresholds) are not analysis numbers and are not regenerated.
Prints PASS/FAIL per check and a total, writes reproduce_report.csv/.md, and exits 1 on any FAIL.

Tolerances. "det" (deterministic) and "seeded" (fixed-seed Monte Carlo from this package) values
must match the canonical value to its printed precision (half a unit of the last printed digit).
Values from the archived v0.6 run, whose seed is unknown, are checked twice: the archived file
must hold the canonical value at printed precision, and the fresh rerun must agree within
Monte Carlo tolerance: shuffle p within 4 x sqrt(2) x sqrt(p(1-p)/5000) (two independent
5,000-shuffle estimates), bootstrap 95% CI endpoints within 0.002 (about 4 x sqrt(2) times the
measured endpoint scatter of 0.0002 to 0.0004), cross-validated R^2 within 0.0075 and Delta R^2
within 0.0015 (fold-assignment scatter 0.0013 and 0.00027). "floor" p-values must equal 1/5,001.

Where the v2 text printed something other than the canonical value, the check row says so in the
v2_text column, and the report lists those rows under "Values the v2 text printed differently
(corrected in v3)": they are errata that v3 corrects, not open errors. The check itself tests the
canonical value. Locations ("main Sec 2", "app 8.4", "main T1" = main Table 1, "README" = this folder's
README.md) refer to the paper files shipped next to this script. The v3 paper exists in two wordings, one with the
angular size control in its headline and one with galaxy size in kpc; the script reads the shipped main paper to see
which, and a value that only the other wording prints is located as "not printed in this version (alternative
wording)". The check itself does not depend on the wording.

Line endings: SHA-256 values are for the files as committed (LF). If the shipped data tables were
checked out with CRLF line endings (Git core.autocrlf=true without a `paper10/** -text` attribute),
their LF-normalized bytes are hashed instead and the check note says so. The run manifests that the analysis scripts
write (sledgehammer_rerun_manifest.json, sensitivity_summary.json) record the SHA-256 of the input bytes as read, so in
such a checkout they show the CRLF files' hashes; every numeric output is the same.
"""
import argparse
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PKG = Path(__file__).resolve().parent
TABLE = 'manga_dynpop_merged_thin_firstpass.csv'
DA_TABLE = 'manga_dynpop_DA_target.csv'
SHIPPED_SHA256 = {
    TABLE: '13336cd42d374514d895d51ca1334f9c32ac22d4c813360be3c6b4fa00f51f58',
    DA_TABLE: '7076255c0b716ca7def8c65bc6de0163e84d6d3d849a78c762b2f596d49a99ec',
}
STRICT = PKG / '00_strict_control_rerun'
V06 = PKG / '01_clean_pipeline_v0_6'
SLEDGE_REF = PKG / '02_sledgehammer_rerun' / 'outputs'
SENS_REF = PKG / '03_size_distance_sensitivity' / 'outputs'
N_SHUFFLE = 5000
FLOOR = 1 / (N_SHUFFLE + 1)
SETS = {'int': 'strict_intSPS_MW_Eps_logSigma', 'obs': 'strict_obsSPS_MW_Eps_logSigma', 'both': 'strict_bothSPS_MW_Eps_logSigma'}

RESULTS = []
TEXT = {'option': None, 'main': ''}
NOT_HERE = 'not printed in this version (alternative wording)'


def detect_text():
    """Which wording of the v3 paper ships next to this script: 'B' = galaxy size in kpc in the headline and Table 1,
    'A' = angular size in the headline, None = neither (for example the v2 files of the staging copy)."""
    try:
        main = (PKG / 'manga_two_layers_paper10_final.md').read_text(encoding='utf-8')
    except OSError:
        main = ''
    TEXT['main'] = main
    TEXT['option'] = 'B' if '| 7 / 8 | +0.0188 |' in main else ('A' if 'That count depends on the size control, which is angular' in main else None)


def only(option, where, elsewhere=None):
    """Location of a value that one wording prints at `where` ('A': angular size headline, 'B': size in kpc). In the
    other wording the value is at `elsewhere`, or not printed at all."""
    if TEXT['option'] is None:
        return where + (' (angular-size wording)' if option == 'A' else ' (physical-size wording)')
    if TEXT['option'] == option:
        return where
    return NOT_HERE if elsewhere is None else elsewhere


def if_printed(marker, where, elsewhere=None):
    """Location of a value printed only when the shipped main paper contains `marker` (an optional sentence)."""
    return where if marker in TEXT['main'] else (NOT_HERE if elsewhere is None else elsewhere)


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run(script, *args):
    cmd = [sys.executable, '-B', str(script)] + [str(a) for a in args]
    print('\n$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def parse_canon(text):
    """'+0.0669' -> (0.0669, 5e-5); '59.2%' -> (0.592, 5e-4); '3.6e-16' -> (3.6e-16, 5e-18); '8' -> (8, 0)."""
    t = text.strip().replace('+', '').replace(',', '')
    pct = t.endswith('%')
    t = t.rstrip('%')
    m = re.fullmatch(r'(-?\d+)(?:\.(\d+))?(?:[eE](-?\d+))?', t)
    if not m:
        raise ValueError(text)
    dec = len(m.group(2) or '')
    exp = int(m.group(3) or 0)
    val = float(t)
    tol = 0.0 if (dec == 0 and exp == 0 and not pct) else 0.5 * 10 ** (exp - dec)
    if pct:
        val, tol = val / 100, tol / 100
    return val, tol * (1 + 1e-9) + 1e-15


def check(cid, where, what, canon, regen, kind='det', archived=None, v2=None, note='', dp=None):
    """Record one check. kind: det | seeded | mc_p | mc_ci | mc_cv_r2 | mc_cv_d | floor | count | bool.
    dp: the decimals the text prints, when it prints an integer-looking value that is a rounding (for example '2' for
    2.05): the tolerance is then half a unit of that decimal instead of exact."""
    status, detail = 'FAIL', ''
    try:
        if kind == 'bool':
            ok = bool(regen)
            target, tol = None, None
        elif kind == 'count':
            target, tol = int(canon), 0
            ok = int(regen) == target
        elif kind == 'floor':
            target, tol = FLOOR, 1e-12
            ok = abs(float(regen) - FLOOR) < 1e-12
        else:
            target, tol = parse_canon(canon)
            if dp is not None:
                tol = 0.5 * 10 ** (-dp) * (1 + 1e-9) + 1e-15
            extra = {'mc_p': 4 * math.sqrt(2) * math.sqrt(max(target * (1 - target), 0) / N_SHUFFLE),
                     'mc_ci': 0.002, 'mc_cv_r2': 0.0075, 'mc_cv_d': 0.0015}.get(kind, 0.0)
            ok = abs(float(regen) - target) <= tol + extra
            if archived is not None:
                arch_ok = abs(float(archived) - target) <= tol
                if not arch_ok:
                    detail = f'archived {archived:.6g} does not round to {canon}'
                ok = ok and arch_ok
            tol = tol + extra
        status = 'PASS' if ok else 'FAIL'
    except Exception as e:  # noqa: BLE001
        detail = f'error: {e}'
    RESULTS.append({'id': cid, 'where': where, 'quantity': what, 'canonical': 'True' if kind == 'bool' else str(canon), 'regenerated': _fmt(regen),
                    'archived_v06': '' if archived is None else _fmt(archived), 'kind': kind,
                    'tolerance': '' if kind in ('bool', 'count') else _fmt(tol if kind != 'floor' else 1e-12),
                    'status': status, 'v2_text': v2 or 'same', 'note': (note + ' ' + detail).strip()})


def ci_check(cid, where, what, canon, lo, hi, kind='det', arch=(None, None), v2=None, note=''):
    a, b = [s.strip() for s in canon.strip('[] ').split(',')]
    check(cid + '.lo', where, what + ' CI lower', a, lo, kind, arch[0], v2, note)
    check(cid + '.hi', where, what + ' CI upper', b, hi, kind, arch[1], v2, note)


def _fmt(x):
    if isinstance(x, (bool, np.bool_)):
        return str(bool(x))
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        return f'{float(x):.6g}'
    except (TypeError, ValueError):
        return str(x)


def md_table(path, heading):
    """Parse the first Markdown table after the line starting with `heading`; first column becomes the index."""
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    start = [i for i, L in enumerate(lines) if L.startswith(heading)][0]
    i = next(j for j in range(start, len(lines)) if lines[j].startswith('|'))
    cells = [c.strip() for c in lines[i].strip().strip('|').split('|')]
    rows = []
    for M in lines[i + 2:]:
        if not M.startswith('|'):
            break
        rows.append([c.strip() for c in M.strip().strip('|').split('|')])
    return pd.DataFrame(rows, columns=cells).set_index(cells[0])


def compare_csv(cid, what, ref, new, rel=1e-9, abs_=1e-12):
    try:
        a = pd.read_csv(ref, float_precision='round_trip')
        b = pd.read_csv(new, float_precision='round_trip')
        ok = list(a.columns) == list(b.columns) and a.shape == b.shape
        worst = 0.0
        if ok:
            for c in a.columns:
                if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
                    x, y = a[c].to_numpy(float), b[c].to_numpy(float)
                    if not np.array_equal(np.isnan(x), np.isnan(y)):
                        ok = False
                        continue
                    m = ~np.isnan(x)
                    d = np.abs(x[m] - y[m])
                    if d.size:
                        worst = max(worst, float(np.max(d / np.maximum(np.abs(x[m]), 1e-300))))
                        ok &= bool(np.all(d <= abs_ + rel * np.abs(x[m])))
                else:
                    ok &= bool((a[c].astype(str) == b[c].astype(str)).all())
        RESULTS.append({'id': cid, 'where': 'package files', 'quantity': what, 'canonical': Path(ref).name,
                        'regenerated': f'max rel diff {worst:.2g}', 'archived_v06': '', 'kind': 'file',
                        'tolerance': f'rel {rel:g}', 'status': 'PASS' if ok else 'FAIL', 'v2_text': '', 'note': ''})
    except Exception as e:  # noqa: BLE001
        RESULTS.append({'id': cid, 'where': 'package files', 'quantity': what, 'canonical': str(ref), 'regenerated': '',
                        'archived_v06': '', 'kind': 'file', 'tolerance': '', 'status': 'FAIL', 'v2_text': '', 'note': f'error: {e}'})


# ----------------------------------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------------------------------
def number_checks(out):
    # regenerated outputs
    st = pd.read_csv(out / 'strict' / 'manga_paper10_strict_control_residual_scoreboard.csv').set_index('control_set')
    stm = pd.read_csv(out / 'strict' / 'manga_paper10_strict_control_model_terms.csv').set_index('control_set')
    str_ = pd.read_csv(out / 'strict' / 'manga_paper10_strict_control_age_coefficients_rawY.csv').set_index('control_set')
    stp = pd.read_csv(out / 'strict' / 'manga_paper10_strict_control_residual_perbin.csv')
    stg = pd.read_csv(out / 'strict' / 'manga_paper10_strict_control_sensitivity_grid.csv')
    sh = pd.read_csv(out / 'sledgehammer' / 'sledgehammer_rerun_scoreboard.csv').set_index('outcome')
    shb = pd.read_csv(out / 'sledgehammer' / 'sledgehammer_rerun_bin_results.csv')
    shp = pd.read_csv(out / 'sledgehammer' / 'sledgehammer_rerun_paired_summary.csv').set_index('outcome')
    shm = json.loads((out / 'sledgehammer' / 'sledgehammer_rerun_manifest.json').read_text(encoding='utf-8'))
    cnt, t50 = shm['sample_counts'], shm['t50_grid_summary']
    t50b = pd.read_csv(out / 'sledgehammer' / 'sledgehammer_rerun_t50_grid_check.csv')
    se = pd.read_csv(out / 'sensitivity' / 'sensitivity_scoreboard.csv').set_index(['variant', 'SPS_set'])
    sesum = json.loads((out / 'sensitivity' / 'sensitivity_summary.json').read_text(encoding='utf-8'))
    facts = sesum['sample_facts']
    dedup = sesum['one_cube_per_mangaid']
    # archived v0.6 values
    v_dml = pd.read_csv(V06 / 'manga_sledgehammer_final_scoreboard.csv').set_index('outcome')
    v_ml = pd.read_csv(V06 / 'manga_ml_split_hardening_scoreboard.csv').set_index('outcome')
    v_ks = pd.read_csv(V06 / 'manga_ml_split_hardening_control_residual_scoreboard.csv').set_index('control_set')
    v_fdm = pd.read_csv(V06 / 'manga_sledgehammer_final_secondary_fdm.csv').set_index('outcome')
    # values archived only in the v0.6 hardening report (Markdown tables)
    rep = V06 / 'manga_ml_split_hardening_final_report.md'
    v_terms = md_table(rep, '## Controlled dynamical M/L age terms')
    v_cv = md_table(rep, '## Controlled dynamical M/L cross-validation')
    v_null = md_table(rep, '## Shuffle null')

    # ---------------- sample and design
    check('S1', 'main Sec 1 (Sample); app 4.2', 'raw three-way join, PlateIFU entries', '10296', cnt['joined_rows'], 'count')
    check('S2', 'main Sec 1 (Sample); app 4.2', 'Qual >= 1 entries', '6065', cnt['qual_ge_1'], 'count')
    check('S3', 'main abstract, Sec 1; app 4.2, 5', 'primary finite sample', '5952', cnt['primary_finite_sample'], 'count')
    check('S4', 'main Sec 1; app 5', 'galaxies per mass bin', '744', cnt['galaxies_per_mass_bin'], 'count')
    check('S5', 'main Sec 1 (Test design); app 5', 'galaxies per young/old quartile', '186', cnt['galaxies_per_quartile'], 'count')
    check('S6', 'main Sec 1 (Sample)', 'distinct MaNGA IDs in the joined table', '10160', cnt['joined_distinct_mangaid'], 'count')
    check('S7', 'main Sec 1 (Sample)', 'distinct MaNGA IDs in the 5,952 sample', '5875', cnt['primary_sample_distinct_mangaid'], 'count')
    check('S8', 'package only (app 4.2 says some galaxies were observed more than once)', 'MaNGA IDs observed more than once in the sample', '72', cnt['primary_sample_mangaid_repeated'], 'count')
    check('S9', 'package only', 'entries belonging to repeated MaNGA IDs', '149', cnt['primary_sample_rows_of_repeated_mangaid'], 'count')
    check('S10', 'main Sec 1, T1 note; app 5', 'within-bin shuffles per test (sledgehammer and sensitivity runs)', '5000',
          shm['n_shuffle'] if shm['n_shuffle'] == sesum['n_shuffle'] else -1, 'count')
    check('S11', 'app 5', 'bootstrap resamples per interval (sledgehammer and sensitivity runs)', '5000',
          shm['n_boot'] if shm['n_boot'] == sesum['n_boot'] else -1, 'count')
    # ---------------- T50 facts (main L33, app 4.3)
    v2_ceiling = 'v2 called 12.85 Gyr the "SPS grid ceiling" and the old quartile "age-saturated"; it is the most common T50 value'
    check('T1', 'main Sec 1 (Test design); app 4.3', 'mass bins whose old-quartile median T50 = 12.85 Gyr', '7', t50['bins_with_old_median_at_mode'], 'count')
    check('T2', 'main Sec 1 (Test design); app 4.3', 'old-quartile median T50 (Gyr) in those bins', '12.85', float(t50b['T50_old_median'].median()))
    check('T3', 'app 4.3', 'old-quartile median T50 in the seventh mass bin (index 6) (Gyr)', '13.01', float(t50b.loc[6, 'T50_old_median']))
    check('T4', 'main Sec 1; app 4.3', 'most common T50 value (Gyr)', '12.85', t50['T50_mode'], v2=v2_ceiling)
    check('T5', 'app 4.3', 'galaxies at the most common T50 value', '278', t50['n_sample_at_mode'], 'count')
    check('T6', 'app 4.3', 'maximum T50 in the sample (Gyr)', '14.13', t50['sample_T50_max'], v2=v2_ceiling)
    check('T7', 'app 4.3', 'sample galaxies above 12.85 Gyr', '605', t50['n_sample_above_mode'], 'count', v2=v2_ceiling)
    check('T8', 'app 4.3', 'T50 grid step (Gyr)', '0.16', t50['typical_grid_step'])
    check('T9', 'package only', 'distinct T50 values in the sample', '85', t50['n_distinct_T50_values'], 'count')
    check('T10', 'main Sec 1; app 4.3 (printed 12.0)', 'lowest old-quartile T50 over all bins (Gyr)', '12.05', t50['old_quartile_T50_min_over_bins'])
    check('T11', 'package only (not printed in v3)', 'highest young-quartile T50 over all bins (Gyr)', '8.19', t50['young_quartile_T50_max_over_bins'])
    check('T12', 'app 4.3 (printed 0.3)', 'old-quartile T50 IQR, smallest bin (Gyr)', '0.32', t50['old_quartile_iqr_range'][0])
    check('T13', 'app 4.3 (printed 0.6)', 'old-quartile T50 IQR, largest bin (Gyr)', '0.64', t50['old_quartile_iqr_range'][1])
    check('T14', 'app 4.3 (printed 1.8)', 'young-quartile T50 IQR, smallest bin (Gyr)', '1.77', t50['young_quartile_iqr_range'][0])
    check('T15', 'app 4.3 (printed 3.2)', 'young-quartile T50 IQR, largest bin (Gyr)', '3.17', t50['young_quartile_iqr_range'][1])
    check('T16', 'package only', 'old-quartile galaxies above 12.85 Gyr, fewest per bin', '62', t50['old_quartile_n_above_mode_range'][0], 'count')
    check('T17', 'package only', 'old-quartile galaxies above 12.85 Gyr, most per bin', '98', t50['old_quartile_n_above_mode_range'][1], 'count')
    check('T18', 'app 5', 'quartile cuts (of 16) falling inside a group of tied T50', '15', t50['quartile_cuts_inside_tie_groups'], 'count')

    # ---------------- main Table 1 (and appendix Sec 6, 7.1, Sec 14 Tables 1-2, which repeat it)
    W1 = 'main T1; app 6 + T1'
    for o, bins, madj, ci, cp_canon, cp_where, cp_v2_app, arch_cp, mean in [
            ('DML_mfl_obs', '8', '+0.0669', ('[+0.050, +0.083]', '[+0.0498, +0.0828]'), '0.006', 'appendix Sec 6, Table 1 and Sec 15', '0.004', 'shuffle_p_count_ge_real', '+0.0742'),
            ('DML_mfl_int', '6', '+0.0772', ('[+0.056, +0.096]', '[+0.0562, +0.0959]'), '0.15', 'appendix Sec 6 and Table 1', '0.14', 'shuffle_p_count_ge_real', '+0.0606')]:
        r, v = sh.loc[o], v_dml.loc[o]
        check(f'D.{o}.bins', W1, f'{o} bins young > old', bins, r['young_gt_old_bins'], 'count')
        check(f'D.{o}.mean', 'app 6 + T1', f'{o} mean bin young-old difference', mean, r['equal_bin_mean_diff_young_minus_old'])
        check(f'D.{o}.madj', W1, f'{o} mass-adjusted young-old difference', madj, r['mass_adjusted_diff_young_minus_old'])
        ci_check(f'D.{o}.ci3', 'main T1', f'{o} mass-adjusted', ci[0], r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'],
                 'mc_ci', (v['diff_bootstrap95_lo'], v['diff_bootstrap95_hi']))
        ci_check(f'D.{o}.ci4', 'app 6 + T1', f'{o} mass-adjusted', ci[1], r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'],
                 'mc_ci', (v['diff_bootstrap95_lo'], v['diff_bootstrap95_hi']))
        check(f'D.{o}.pcount', 'main T1; app 6, T1, Sec 15', f'{o} bin-count shuffle p', cp_canon, r['p_young_gt_old_count'], 'mc_p',
              v[arch_cp], v2=f'v2 main printed {cp_canon}; v2 {cp_where} printed {cp_v2_app}, the fair-coin '
              f'sign-test tail ({r["binomial_reference_p_young_gt_old_count"]:.4f}), not the archived shuffle p ({v[arch_cp]:.4f})')
        check(f'D.{o}.pmag', W1, f'{o} magnitude shuffle p', '0.0002', r['p_young_gt_old_sum'], 'floor')
    neg_bins = shb[(shb['outcome'] == 'DML_mfl_int') & (shb['median_diff_young_minus_old'] < 0)]['mass_bin'].tolist()
    check('D.int.negbins', 'app 6', 'DML_mfl_int: the two young<old bins are low-mass bins (index 1 and 2)', None, neg_bins == [1, 2], 'bool')
    for o, label, madj, ci3, ci4, mean in [
            ('mfl_cyl_log_ML_dyn', 'JAM dyn', '-0.1075', '[-0.120, -0.092]', '[-0.1200, -0.0920]', '-0.1045'),
            ('sp_ML_int_Re', 'SPS int', '-0.1693', '[-0.179, -0.159]', '[-0.1790, -0.1587]', '-0.1519'),
            ('sp_ML_obs_Re', 'SPS obs', '-0.1569', '[-0.168, -0.148]', '[-0.1679, -0.1478]', '-0.1579')]:
        r, v = sh.loc[o], v_ml.loc[o]
        check(f'C.{o}.yo', 'main T1; app 7.1, T2', f'{label}: bins young > old', '0', r['young_gt_old_bins'], 'count')
        check(f'C.{o}.oy', 'app 7.1, T2', f'{label}: bins old > young', '8', r['old_gt_young_bins'], 'count')
        check(f'C.{o}.mean', 'app 7.1', f'{label}: mean bin young-old difference', mean, r['equal_bin_mean_diff_young_minus_old'])
        check(f'C.{o}.madj', 'main T1; app 7.1, T2', f'{label}: mass-adjusted young-old difference', madj, r['mass_adjusted_diff_young_minus_old'])
        ci_check(f'C.{o}.ci3', 'main T1', label + ' mass-adjusted', ci3, r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'],
                 'mc_ci', (v['mass_adjusted_diff_boot95_lo'], v['mass_adjusted_diff_boot95_hi']))
        ci_check(f'C.{o}.ci4', 'app 7.1, T2', label + ' mass-adjusted', ci4, r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'],
                 'mc_ci', (v['mass_adjusted_diff_boot95_lo'], v['mass_adjusted_diff_boot95_hi']))
        check(f'C.{o}.p_yo', 'main T1 (Count p)', f'{label}: young > old count p (1 by construction at 0/8)', '1.000', r['p_young_gt_old_count'])
        arch_oy = float(v_null.loc[o, 'p_negative_count_ge_real'])
        check(f'C.{o}.p_oy', 'package only', f'{label}: old > young count p (8/8)', '0.005', r['p_old_gt_young_count'], 'mc_p', arch_oy,
              note='archived 0.0048 / 0.0052 / 0.0052; "below 0.01" is the robust wording')
        if o == 'mfl_cyl_log_ML_dyn':
            check(f'C.{o}.dirp', 'app 7.1', f'{label}: directional (young > old) magnitude p', '1.0000', r['p_directional_sum'])
        else:
            check(f'C.{o}.dirp', 'app 7.1', f'{label}: directional (old > young) magnitude p', '0.0002', r['p_directional_sum'], 'floor',
                  v2='v2 headed the column "Directional shuffle p"; the value is the magnitude p')
    W3 = only('A', 'main T1; app 8.3', 'app 8.3')
    for s, bins, madj, ci3, ci4, cp, mean in [('obs', '7', '+0.0387', '[+0.029, +0.049]', '[+0.0292, +0.0494]', '0.031', '+0.0393'),
                                              ('both', '7', '+0.0344', '[+0.023, +0.044]', '[+0.0231, +0.0442]', '0.034', '+0.0347'),
                                              ('int', '7', '+0.0225', '[+0.011, +0.034]', '[+0.0114, +0.0342]', '0.037', '+0.0238')]:
        r = st.loc[SETS[s]]
        check(f'R.{s}.bins', W3, f'strict residual ({s} SPS): bins young > old', bins, r['young_gt_old_bins'], 'count')
        check(f'R.{s}.madj', W3, f'strict residual ({s} SPS): mass-adjusted young-old', madj, r['mass_adjusted_diff_young_minus_old'], 'seeded')
        ci_check(f'R.{s}.ci3', only('A', 'main T1'), f'strict residual ({s} SPS) mass-adjusted', ci3, r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'], 'seeded')
        ci_check(f'R.{s}.ci4', 'app 8.3', f'strict residual ({s} SPS) mass-adjusted', ci4, r['mass_adjusted_diff_boot95_lo'], r['mass_adjusted_diff_boot95_hi'], 'seeded')
        check(f'R.{s}.pcount', W3, f'strict residual ({s} SPS): bin-count p', cp, r['shuffle_p_count_ge_real'], 'seeded')
        check(f'R.{s}.pmag', W3, f'strict residual ({s} SPS): magnitude p', '0.0002', r['shuffle_p_sum_ge_real'], 'floor')
        check(f'R.{s}.mean', 'app 8.3', f'strict residual ({s} SPS): mean bin young-old', mean, r['equal_bin_mean_diff_young_minus_old'], 'seeded')
    for o, name, bins, mean, pc, pm, pj in [('gnfw_cyl_fdm_Re', 'gNFW fDM', '5', '+0.0011', '0.35', '0.48', '0.3007'),
                                             ('nfw_cyl_fdm_Re', 'NFW fDM', '4', '-0.0035', '0.62', '0.65', '0.5135')]:
        r, v = sh.loc[o], v_fdm.loc[o]
        check(f'F.{o}.bins', 'main T1; app 6', f'{name}: bins young > old', bins, r['young_gt_old_bins'], 'count')
        check(f'F.{o}.mean', 'main T1; app 6', f'{name}: mean bin young-old difference', mean, r['equal_bin_mean_diff_young_minus_old'])
        check(f'F.{o}.pc', 'main T1; app 6', f'{name}: count p', pc, r['p_young_gt_old_count'], 'mc_p', v['shuffle_p_count_ge_real'])
        check(f'F.{o}.pm', 'main T1; app 6', f'{name}: magnitude p', pm, r['p_young_gt_old_sum'], 'mc_p', v['shuffle_p_mean_diff_ge_real'])
        check(f'F.{o}.pj', 'package only (archived v0.6 joint p; v3 app 6 prints count p and magnitude p)', f'{name}: joint p', pj, r['p_young_gt_old_joint'], 'mc_p', v['shuffle_p_joint_count_and_mean_ge_real'])
    check('F.gnfw.madj', 'package only', 'gNFW fDM: mass-adjusted young-old difference', '-0.00075', sh.loc['gnfw_cyl_fdm_Re', 'mass_adjusted_diff_young_minus_old'],
          note='float value -0.000749999...; the CSV prints -0.00075')
    check('F.nfw.madj', 'package only', 'NFW fDM: mass-adjusted young-old difference', '-0.0055', sh.loc['nfw_cyl_fdm_Re', 'mass_adjusted_diff_young_minus_old'])
    # paired check (app 6)
    for o, wins, frac in [('DML_mfl_int', '881', '59.2%'), ('DML_mfl_obs', '913', '61.4%')]:
        r = shp.loc[o]
        check(f'P.{o}.wins', 'app 6', f'{o} paired check: young wins of 1,488 pairs', wins, r['young_wins'], 'count')
        check(f'P.{o}.n', 'app 6', f'{o} paired check: number of pairs', '1488', r['n_pairs'], 'count')
        check(f'P.{o}.frac', 'app 6', f'{o} paired check: young win fraction', frac, r['young_win_fraction'])
        check(f'P.{o}.p', 'app 6', f'{o} paired check: shuffle p', '0.0002', r['shuffle_p_win_fraction_ge_real'], 'floor')

    # ---------------- main text
    check('M.rho', 'app 8.3', 'Spearman rho, residual vs T50, strict observed-SPS', '-0.117', stm.loc[SETS['obs'], 'resid_T50_spearman_rho'], 'det')
    check('M.rho.p', 'app 8.3', 'Spearman p, strict observed-SPS', '1.7e-19', stm.loc[SETS['obs'], 'resid_T50_spearman_p'])
    check('M.rho.both', 'app 8.3', 'Spearman rho, strict both-SPS', '-0.103', stm.loc[SETS['both'], 'resid_T50_spearman_rho'])
    check('M.rho.both.p', 'app 8.3', 'Spearman p, strict both-SPS', '1.4e-15', stm.loc[SETS['both'], 'resid_T50_spearman_p'])
    check('M.rho.int', 'app 8.3', 'Spearman rho, strict intrinsic-SPS', '-0.065', stm.loc[SETS['int'], 'resid_T50_spearman_rho'])
    check('M.rho.int.p', 'app 8.3', 'Spearman p, strict intrinsic-SPS', '5.4e-7', stm.loc[SETS['int'], 'resid_T50_spearman_p'])
    check('M.grid', 'app 8.3', 'reduced-spec grid: variants (of 12) with 7/8 bins', '12', int((stg['young_gt_old_bins'] == 7).sum()), 'count')
    check('M.range.lo', only('A', 'main Sec 2 (Layer 2)'), 'smallest strict mass-adjusted residual, 3 dp', '+0.022', st.loc[SETS['int'], 'mass_adjusted_diff_young_minus_old'], 'seeded',
          v2='v2 main Sec 2 printed +0.023 (0.022474 rounded twice); Table 1 prints +0.0225')
    check('M.range.hi', only('A', 'main Sec 2 (Layer 2)'), 'largest strict mass-adjusted residual, 3 dp', '+0.039', st.loc[SETS['obs'], 'mass_adjusted_diff_young_minus_old'], 'seeded')
    pcs = st['shuffle_p_count_ge_real']
    check('M.pc.range', only('A', 'main T1'), 'strict bin-count p all between 0.03 and 0.04', None, bool(((pcs >= 0.03) & (pcs < 0.04)).all()), 'bool')
    check('M.int.smaller', 'app 10.4', 'intrinsic-SPS residual smaller than observed-SPS', None,
          st.loc[SETS['int'], 'mass_adjusted_diff_young_minus_old'] < st.loc[SETS['obs'], 'mass_adjusted_diff_young_minus_old'], 'bool')
    b5 = stp[stp['mass_bin'] == 5]
    check('M.bin5.neg', 'app 8.3; README', 'the sixth mass bin (index 5) is the only young<old bin in all three strict sets', None,
          bool((b5['median_diff_young_minus_old'] < 0).all() and (stp[stp['mass_bin'] != 5]['median_diff_young_minus_old'] > 0).all()), 'bool')
    check('M.bin5.ci', 'app 8.3', 'its bootstrap interval crosses zero in each set', None,
          bool(((b5['bootstrap95_lo'] < 0) & (b5['bootstrap95_hi'] > 0)).all()), 'bool')
    check('M.bin5.logM', 'app 8.3 (printed 10.66); README', 'its median NSA log M (h^-2 Msun)', '10.663', float(b5['logM_median'].iloc[0]))

    # ---------------- appendix Section 8 (kitchen-sink controls, v0.6 hardening run) and the strict counterparts
    ksr = {s: se.loc[('f_kitchen_sink', s)] for s in ['int', 'obs', 'both']}
    ksn = {'int': 'controls_plus_intrinsic_SPS_ML', 'obs': 'controls_plus_observed_SPS_ML', 'both': 'controls_plus_both_SPS_ML'}
    ksv = {'int': 'intSPS_controls', 'obs': 'obsSPS_controls', 'both': 'bothSPS_controls'}
    for s, coef, p, r2, c0, c1, dr in [('int', '-0.0215', '8.1e-12', '0.4799', '0.4701', '0.4753', '+0.0052'),
                                      ('obs', '-0.0346', '4.9e-24', '0.5057', '0.4862', '0.5009', '+0.0146'),
                                      ('both', '-0.0357', '1.4e-24', '0.5058', '0.4871', '0.5008', '+0.0137')]:
        r, vt, vc = ksr[s], v_terms.loc[ksn[s]], v_cv.loc[ksn[s]]
        lab = f'kitchen-sink ({s} SPS)'
        check(f'K.{s}.N', 'app 14, T3', f'{lab}: N', '5952', r['N'], 'count')
        check(f'K.{s}.coef', 'app 14, T3 note', f'{lab}: T50 coefficient, dex per 1 SD of T50', coef, r['T50_coef_raw_per_1sd'],
              archived=float(vt['age_T50_standardized_coef']), v2='v2 app 8.1 headed the column "Standardized"; it is dex per 1 SD')
        check(f'K.{s}.p', 'archived v0.6 report (kitchen-sink; not printed in v3)', f'{lab}: HC3 p of the T50 coefficient', p, r['T50_coef_raw_p_HC3'], archived=float(vt['age_T50_robust_p']))
        check(f'K.{s}.r2', 'app 14, T3 note', f'{lab}: model R^2', r2, r['model_R2'], archived=float(vt['model_R2']))
        check(f'K.{s}.cv0', 'archived v0.6 report (kitchen-sink; not printed in v3)', f'{lab}: CV R^2 controls only', c0, r['cv_R2_controls_only'], 'mc_cv_r2', float(vc['cv_R2_controls_only_mean']),
              note='v0.6 fold seed unknown; rerun uses KFold(random_state=20260629)')
        check(f'K.{s}.cv1', 'archived v0.6 report (kitchen-sink; not printed in v3)', f'{lab}: CV R^2 controls + age', c1, r['cv_R2_controls_plus_age'], 'mc_cv_r2', float(vc['cv_R2_controls_plus_age_mean']))
        check(f'K.{s}.cvd', 'app 14, T3 note', f'{lab}: CV Delta R^2 from age', dr, r['cv_delta_R2_from_age'], 'mc_cv_d', float(vc['cv_delta_R2_from_age']))
    for s, bins, mean, madj, ci, cp, v2b, v2p, v2m in [
            ('obs', '8', '+0.0359', '+0.0358', '[+0.0270, +0.0475]', '0.004', '7-8 / 8', '0.004-0.035', None),
            ('both', '8', '+0.0342', '+0.0311', '[+0.0216, +0.0436]', '0.004', '7-8 / 8', '0.004-0.035', None),
            ('int', '6', '+0.0209', '+0.0188', '[+0.0092, +0.0312]', '0.14', None, None, '+0.0210 (0.020950 rounded twice)')]:
        r, v = ksr[s], v_ks.loc[ksv[s]]
        lab = f'kitchen-sink residual ({s} SPS)'
        mixed = 'v2 app 8.3 and Table 3 printed a two-specification range; the kitchen-sink value alone is canonical'
        check(f'K.{s}.bins', if_printed('8/8, 8/8, and 6/8 bins for the observed, both, and intrinsic sets', 'main Sec 2; app 8.3, 14 T3; README', 'app 8.3, 14 T3; README'), f'{lab}: bins young > old', bins, r['young_gt_old_bins'], 'count',
              v2=f'v2 printed {v2b} ({mixed})' if v2b else ('v2 app 8.3 and 10.4 stated 6/8 without saying it was the kitchen-sink set (strict: 7/8); '
                                                           'the 5773e94 READMEs called the kitchen-sink result 8/8' if s == 'int' else None))
        check(f'K.{s}.mean', 'app 14, T3', f'{lab}: mean bin young-old', mean, r['equal_bin_mean_diff'], v2=f'v2 printed {v2m}' if v2m else None)
        check(f'K.{s}.madj', 'app 14, T3', f'{lab}: mass-adjusted young-old', madj, r['mass_adjusted_diff'], archived=float(v['mass_adjusted_diff_young_minus_old']),
              v2=None if s == 'both' else 'v2 app 16.2 quoted +0.0188 to +0.0358 (kitchen-sink) as "the observed" residual; main Table 1 (strict) gives +0.0225 to +0.0387')
        ci_check(f'K.{s}.ci', 'app 14, T3', lab + ' mass-adjusted', ci, r['mass_adjusted_ci_lo'], r['mass_adjusted_ci_hi'], 'mc_ci',
                 (v['mass_adjusted_diff_boot95_lo'], v['mass_adjusted_diff_boot95_hi']))
        check(f'K.{s}.pc', 'app 14, T3', f'{lab}: bin-count p', cp, r['p_count'], 'mc_p', float(v['shuffle_p_count_ge_real']),
              v2=f'v2 printed {v2p} ({mixed}; 0.035 is the sign-test tail 9/256)' if v2p else None)
        check(f'K.{s}.pm', 'app 14, T3', f'{lab}: magnitude p', '0.0002', r['p_magnitude'], 'floor')
    check('K.rho', 'app 8.3', 'Spearman rho, kitchen-sink observed-SPS residual vs T50', '-0.119', ksr['obs']['spearman_rho_resid_T50'],
          v2='unlabeled in v2, while the v2 main paper quoted the strict value -0.117 under the same words')
    check('K.rho.p', 'archived v0.6 report (not printed in v3)', 'Spearman p, kitchen-sink observed-SPS', '2e-20', ksr['obs']['spearman_p'])
    check('K.rho.both', 'app 8.3', 'Spearman rho, kitchen-sink both-SPS', '-0.110', ksr['both']['spearman_rho_resid_T50'])
    check('K.rho.int', 'app 8.3', 'Spearman rho, kitchen-sink intrinsic-SPS', '-0.067', ksr['int']['spearman_rho_resid_T50'])
    # strict-spec Sections 8.1/8.2 of the v3 appendix
    for s, coef, p, r2, c0, c1, dr, std in [('int', '-0.0267', '3.6e-16', '0.4258', '0.4126', '0.4212', '+0.0086', '-0.139'),
                                            ('obs', '-0.0404', '6.0e-32', '0.4560', '0.4303', '0.4513', '+0.0209', '-0.210'),
                                            ('both', '-0.0416', '4.7e-31', '0.4561', '0.4320', '0.4512', '+0.0192', '-0.217')]:
        r, m = str_.loc[SETS[s]], stm.loc[SETS[s]]
        lab = f'strict ({s} SPS)'
        check(f'X.{s}.coef', 'app 8.1', f'{lab}: T50 coefficient, dex per 1 SD of T50', coef, r['age_T50_coef_raw_outcome_per_1sd_T50'])
        check(f'X.{s}.p', 'app 8.1', f'{lab}: HC3 p', p, r['age_T50_robust_p'])
        check(f'X.{s}.r2', 'app 8.1', f'{lab}: model R^2', r2, r['model_R2'])
        check(f'X.{s}.std', 'package only (model_terms.csv)', f'{lab}: fully standardized T50 coefficient', std, m['age_T50_standardized_coef'])
        check(f'X.{s}.cv0', 'app 8.2', f'{lab}: CV R^2 controls only', c0, m['cv_R2_controls_only_mean'], 'seeded')
        check(f'X.{s}.cv1', 'app 8.2', f'{lab}: CV R^2 controls + age', c1, m['cv_R2_controls_plus_age_mean'], 'seeded')
        check(f'X.{s}.cvd', 'app 8.2', f'{lab}: CV Delta R^2 from age', dr, m['cv_delta_R2_from_age'], 'seeded')

    # ---------------- canonical size / distance / subsample sensitivity (03_size_distance_sensitivity)
    W5 = 'sensitivity reference values (main Sec 2; app 8.4; README)'
    sens = [  # variant, set, bins, mass-adj, CI, count p, magnitude p, per-Gyr coef, t, rho
        ('a_published', 'int', '7', '+0.0225', '[+0.0114, +0.0342]', '0.037', '0.0002', '-0.00905', '-8.15', '-0.065'),
        ('a_published', 'obs', '7', '+0.0387', '[+0.0292, +0.0494]', '0.031', '0.0002', '-0.01369', '-11.76', '-0.117'),
        ('a_published', 'both', '7', '+0.0344', '[+0.0231, +0.0442]', '0.034', '0.0002', '-0.01410', '-11.59', '-0.103'),
        ('b_kpc', 'int', '3', '-0.0008', '[-0.0098, +0.0084]', '0.854', '0.525', '-0.00550', '-5.49', '-0.027'),
        ('b_kpc', 'obs', '7', '+0.0188', '[+0.0110, +0.0270]', '0.033', '0.0002', '-0.01006', '-10.17', '-0.095'),
        ('b_kpc', 'both', '6', '+0.0162', '[+0.0080, +0.0244]', '0.143', '0.0002', '-0.01059', '-10.24', '-0.083'),
        ('c_published_plus_logDA', 'int', '3', '-0.0084', '[-0.0170, +0.0003]', '0.865', '0.893', '-0.00437', '-4.19', '-0.003'),
        ('c_published_plus_logDA', 'obs', '6', '+0.0088', '[+0.0012, +0.0165]', '0.129', '0.010', '-0.00908', '-8.86', '-0.065'),
        ('c_published_plus_logDA', 'both', '3', '+0.0057', '[-0.0007, +0.0156]', '0.850', '0.055', '-0.00962', '-8.95', '-0.057'),
        ('d_kpc_plus_logDA', 'int', '3', '-0.0084', '[-0.0170, +0.0003]', '0.865', '0.893', '-0.00437', '-4.19', '-0.003'),
        ('d_kpc_plus_logDA', 'obs', '6', '+0.0088', '[+0.0012, +0.0165]', '0.129', '0.010', '-0.00908', '-8.86', '-0.065'),
        ('d_kpc_plus_logDA', 'both', '3', '+0.0057', '[-0.0007, +0.0156]', '0.850', '0.055', '-0.00962', '-8.95', '-0.057'),
        ('e1_primary', 'int', '6', '+0.0125', '[-0.0019, +0.0262]', '0.132', '0.023', '-0.00883', '-5.49', '-0.070'),
        ('e1_primary', 'obs', '6', '+0.0209', '[+0.0102, +0.0342]', '0.128', '0.0004', '-0.01175', '-5.86', '-0.103'),
        ('e1_primary', 'both', '6', '+0.0217', '[+0.0109, +0.0353]', '0.137', '0.0004', '-0.01337', '-7.05', '-0.106'),
        ('e2_secondary', 'int', '2', '-0.0192', '[-0.0327, -0.0077]', '0.963', '0.956', '-0.00141', '-0.84', '+0.032'),
        ('e2_secondary', 'obs', '4', '-0.0032', '[-0.0154, +0.0108]', '0.623', '0.426', '-0.00728', '-5.30', '-0.044'),
        ('e2_secondary', 'both', '3', '-0.0088', '[-0.0233, +0.0029]', '0.850', '0.689', '-0.00650', '-4.02', '-0.019'),
        ('e3_primary_kpc', 'int', '4', '+0.0029', '[-0.0086, +0.0166]', '0.615', '0.073', '-0.00587', '-3.74', '-0.036'),
        ('e3_primary_kpc', 'obs', '4', '+0.0110', '[+0.0005, +0.0253]', '0.618', '0.010', '-0.00887', '-4.68', '-0.077'),
        ('e3_primary_kpc', 'both', '4', '+0.0138', '[+0.0047, +0.0280]', '0.623', '0.008', '-0.01065', '-5.99', '-0.086'),
        ('e4_secondary_kpc', 'int', '2', '-0.0228', '[-0.0417, -0.0138]', '0.964', '1.000', '-0.00040', '-0.22', '+0.053'),
        ('e4_secondary_kpc', 'obs', '4', '-0.0079', '[-0.0220, +0.0056]', '0.616', '0.786', '-0.00667', '-4.54', '-0.026'),
        ('e4_secondary_kpc', 'both', '3', '-0.0133', '[-0.0269, -0.0041]', '0.842', '0.913', '-0.00583', '-3.35', '-0.001'),
        ('e5_primary_plus_logDA', 'int', '4', '+0.0033', '[-0.0087, +0.0168]', '0.611', '0.077', '-0.00624', '-3.92', '-0.039'),
        ('e5_primary_plus_logDA', 'obs', '4', '+0.0116', '[+0.0010, +0.0257]', '0.616', '0.008', '-0.00927', '-4.96', '-0.079'),
        ('e5_primary_plus_logDA', 'both', '4', '+0.0142', '[+0.0047, +0.0280]', '0.624', '0.007', '-0.01101', '-6.19', '-0.087'),
        ('f_kitchen_sink', 'int', '6', '+0.0188', '[+0.0090, +0.0312]', '0.136', '0.0002', '-0.00728', '-6.84', '-0.067'),
        ('f_kitchen_sink', 'obs', '8', '+0.0358', '[+0.0264, +0.0473]', '0.004', '0.0002', '-0.01171', '-10.11', '-0.119'),
        ('f_kitchen_sink', 'both', '8', '+0.0311', '[+0.0215, +0.0433]', '0.005', '0.0002', '-0.01209', '-10.23', '-0.110'),
    ]
    for var, s, bins, madj, ci, pc, pm, gyr, t, rho in sens:
        r = se.loc[(var, s)]
        cid = f'Z.{var}.{s}'
        lab = f'{var} ({s} SPS)'
        check(cid + '.bins', W5, f'{lab}: bins young > old', bins, r['young_gt_old_bins'], 'count')
        check(cid + '.madj', W5, f'{lab}: mass-adjusted residual', madj, r['mass_adjusted_diff'], 'seeded')
        ci_check(cid + '.ci', W5, lab + ' mass-adjusted', ci, r['mass_adjusted_ci_lo'], r['mass_adjusted_ci_hi'], 'seeded')
        if pm == '0.0002':
            check(cid + '.pm', W5, f'{lab}: magnitude p', pm, r['p_magnitude'], 'floor')
        else:
            check(cid + '.pm', W5, f'{lab}: magnitude p', pm, r['p_magnitude'], 'seeded')
        check(cid + '.pc', W5, f'{lab}: count p', pc, r['p_count'], 'seeded')
        check(cid + '.gyr', W5, f'{lab}: T50 coefficient, dex per Gyr', gyr, r['T50_coef_per_Gyr'])
        check(cid + '.t', W5, f'{lab}: HC3 t of the T50 coefficient', t, r['T50_coef_raw_t_HC3'])
        check(cid + '.rho', W5, f'{lab}: Spearman rho residual vs T50', rho, r['spearman_rho_resid_T50'])
    check('Z.c.obs.rho.p', W5, 'published + logDA (obs SPS): Spearman p', '6e-7', se.loc[('c_published_plus_logDA', 'obs'), 'spearman_p'])
    for s, p in [('obs', '2e-13'), ('both', '1e-10'), ('int', '0.04')]:
        check(f'Z.b.{s}.rho.p', 'app 8.4', f'kpc size ({s} SPS): Spearman p, residual vs T50', p,
              se.loc[('b_kpc', s), 'spearman_p'])
    check('Z.c.int.rho.p', 'app 8.4', 'published + logDA (int SPS): Spearman p, residual vs T50', '0.81',
          se.loc[('c_published_plus_logDA', 'int'), 'spearman_p'])
    t_da = [float(se.loc[('d_kpc_plus_logDA', s), 't_logDA_HC3']) for s in ['int', 'obs', 'both']]
    check('Z.d.tDA.lo', 'app 8.4', 'kpc size + logDA with T50 in the model: smallest HC3 t of logDA over the SPS sets', '3.8', min(t_da))
    check('Z.d.tDA.hi', 'app 8.4', 'kpc size + logDA with T50 in the model: largest HC3 t of logDA over the SPS sets', '4.5', max(t_da))
    full_rows = [k for k in se.index if not k[0].startswith('e')]
    check('Z.linear.all', 'main abstract, Sec 2, Sec 5; app 1, 8.4, 9, 15', 'linear T50 term significant (HC3 p < 0.001) in every full-sample '
          'variant: the 15 full-sample scoreboard fits (variants a to d and the kitchen-sink set) and the 48 size-treatment grid fits', None,
          bool((se.loc[full_rows, 'T50_coef_raw_p_HC3'] < 1e-3).all()) and sesum['grid_T50_p_HC3_max'] < 1e-3, 'bool')
    check('Z.grid.t.lo', W5, 'size-treatment grid: most negative HC3 t of the linear T50 term', '-12.7', sesum['grid_T50_t_HC3_range'][0])
    check('Z.grid.t.hi', W5, 'size-treatment grid: least negative HC3 t of the linear T50 term', '-4.2', sesum['grid_T50_t_HC3_range'][1])
    att = facts['attenuation_old_minus_young_by_bin']
    check('Z.att.low', 'app 10.4', 'SPS attenuation term (sp_ML_obs_Re - sp_ML_int_Re): old-quartile median above '
          'young-quartile median in the four lowest mass bins', None, all(x > 0 for x in att[:4]), 'bool')
    for s, size, da in [('int', '0.54', '0.61'), ('obs', '0.56', '0.62'), ('both', '0.57', '0.62')]:
        r = se.loc[('c_published_plus_logDA', s)]
        check(f'Z.coefs.{s}.size', W5, f'published + logDA ({s} SPS): raw coefficient of logRe (arcsec)', size, r['coef_size_term_raw'])
        check(f'Z.coefs.{s}.da', W5, f'published + logDA ({s} SPS): raw coefficient of logDA', da, r['coef_logDA_raw'])
    check('Z.c_eq_d', W5, 'kpc + logDA residuals identical to published + logDA (max |diff| < 1e-10)', None,
          sesum['max_abs_resid_difference_c_vs_d'] < 1e-10, 'bool')
    check('Z.da.lo', W5, 'young/old median DA ratio, smallest bin', '1.18', facts['DA_ratio_range'][0])
    check('Z.da.hi', W5, 'young/old median DA ratio, largest bin', '1.67', facts['DA_ratio_range'][1])
    check('Z.sec.lo', W5, 'young/old Secondary-sample fraction ratio, smallest bin', '2.05', facts['secondary_ratio_range'][0])
    check('Z.sec.hi', W5, 'young/old Secondary-sample fraction ratio, largest bin', '3.54', facts['secondary_ratio_range'][1])
    check('Z.logM', W5, 'young vs old median log M differ by at most (dex)', '0.03', facts['max_abs_logM_median_diff_young_old'])
    for k_, n_ in [('0', '2776'), ('1', '2337'), ('2', '839')]:
        check(f'Z.target{k_}', W5, f'sample entries with target = {k_} (0 Primary, 1 Secondary, 2 color-enhanced)', n_,
              facts['target_counts_in_sample'][k_], 'count')
    check('Z.rho.da', W5, 'Spearman rho of log DA vs T50 in the sample', '-0.180', facts['spearman_logDA_vs_T50'])
    check('Z.re.as', W5, 'median Re (arcsec)', '5.27', facts['median_Re_arcsec'])
    check('Z.re.kpc', W5, 'median Re (kpc)', '4.51', facts['median_Re_kpc'])
    gh = sesum['grid_bin_count_histogram_by_size']
    for size, hist in [('kpc', {'7': 7, '6': 3, '3': 2}), ('arcsec + logDA', {'6': 8, '5': 1, '3': 3})]:
        for b, n in hist.items():
            check(f'Z.grid.{size}.{b}', W5, f'12-variant grid with size = {size}: variants at {b}/8', str(n), gh[size].get(b, 0), 'count')
    # ---------------- the metallicity side of the size/distance grid (main Sec 2, abstract, app 8.4, README)
    gr = pd.read_csv(out / 'sensitivity' / 'reduced_spec_grid_by_size_unit.csv')
    lw, mw = gr['metallicity'] == 'sp_LW_Metal_Re', gr['metallicity'] == 'sp_MW_Metal_Re'
    kpc, dist = gr['size'] == 'kpc', gr['size'].str.contains('logDA')
    WG = 'main abstract, Sec 2; app 8.4; README'
    check('Z.grid.kpc.LW', WG, 'kpc size: all six light-weighted-metallicity variants at 7/8', None, bool((gr[kpc & lw]['young_gt_old_bins'] == 7).all()) and int((kpc & lw).sum()) == 6, 'bool')
    check('Z.grid.kpc.MW.lo', WG, 'kpc size: fewest bins among the mass-weighted variants', '3', int(gr[kpc & mw]['young_gt_old_bins'].min()), 'count')
    check('Z.grid.kpc.MW.hi', WG, 'kpc size: most bins among the mass-weighted variants', '7', int(gr[kpc & mw]['young_gt_old_bins'].max()), 'count')
    for size in ['arcsec + logDA', 'kpc + logDA']:
        g3 = gr[(gr['size'] == size) & (gr['young_gt_old_bins'] == 3)]
        check(f'Z.grid.{size}.3isMW', 'app 8.4; README', f'{size}: every 3/8 variant uses mass-weighted metallicity', None, len(g3) == 3 and bool((g3['metallicity'] == 'sp_MW_Metal_Re').all()), 'bool')
        check(f'Z.grid.{size}.ge5', 'app 8.4; README', f'{size}: variants at 5/8 or 6/8 (of 12)', '9', int(gr[(gr['size'] == size) & gr['young_gt_old_bins'].between(5, 6)].shape[0]), 'count')
    lwb = gr[(kpc | dist) & lw]['young_gt_old_bins']
    check('Z.grid.LW.lo', 'app 8.4; README', 'light-weighted metallicity, kpc size or log DA: fewest bins in any SPS set', '5', int(lwb.min()), 'count')
    check('Z.grid.LW.hi', 'app 8.4; README', 'light-weighted metallicity, kpc size or log DA: most bins in any SPS set', '7', int(lwb.max()), 'count')
    check('Z.grid.obs.ge6', 'README', 'observed SPS keeps at least 6/8 in every full-sample variant (scoreboard and grid)', None,
          bool(gr[gr['SPS_set'] == 'obs']['young_gt_old_bins'].ge(6).all())
          and bool(se.loc[[k for k in se.index if not k[0].startswith('e') and k[1] == 'obs'], 'young_gt_old_bins'].ge(6).all()), 'bool')
    # ---------------- every v3 number at the precision the text prints it; the canonical checks above use more digits
    WP = 'v3 text, printed precision'
    t1 = {('a_published', 'int'): '-8.2', ('a_published', 'obs'): '-11.8', ('a_published', 'both'): '-11.6',
          ('b_kpc', 'int'): '-5.5', ('b_kpc', 'obs'): '-10.2', ('b_kpc', 'both'): '-10.2',
          ('c_published_plus_logDA', 'int'): '-4.2', ('c_published_plus_logDA', 'obs'): '-8.9', ('c_published_plus_logDA', 'both'): '-9.0',
          ('e1_primary', 'int'): '-5.5', ('e1_primary', 'obs'): '-5.9', ('e1_primary', 'both'): '-7.0',
          ('e3_primary_kpc', 'int'): '-3.7', ('e3_primary_kpc', 'obs'): '-4.7', ('e3_primary_kpc', 'both'): '-6.0',
          ('e2_secondary', 'int'): '-0.8', ('e2_secondary', 'obs'): '-5.3', ('e2_secondary', 'both'): '-4.0',
          ('e4_secondary_kpc', 'int'): '-0.2', ('e4_secondary_kpc', 'obs'): '-4.5', ('e4_secondary_kpc', 'both'): '-3.3'}
    for (var, sps), t in t1.items():
        check(f'V.t.{var}.{sps}', 'app 8.4 table', f'{var} ({sps} SPS): HC3 t of the linear T50 term, 1 dp', t, se.loc[(var, sps), 'T50_coef_raw_t_HC3'])
    check('V.T10', 'main Sec 1; app 4.3', 'every old-quartile T50 at least 12.0 Gyr (lowest, 1 dp)', '12.0', t50['old_quartile_T50_min_over_bins'])
    check('V.T11', 'package only (not printed in v3)', 'every young-quartile T50 below 8.2 Gyr', None, t50['young_quartile_T50_max_over_bins'] < 8.2, 'bool')
    for cid, canon, val in [('V.T12', '0.3', t50['old_quartile_iqr_range'][0]), ('V.T13', '0.6', t50['old_quartile_iqr_range'][1]),
                            ('V.T14', '1.8', t50['young_quartile_iqr_range'][0]), ('V.T15', '3.2', t50['young_quartile_iqr_range'][1])]:
        check(cid, 'app 4.3', 'quartile T50 interquartile range bound, 1 dp (Gyr)', canon, val)
    check('V.da.lo', 'app 8.4; README', 'young/old median DA ratio, smallest bin, 1 dp', '1.2', facts['DA_ratio_range'][0])
    check('V.da.hi', 'app 8.4; README', 'young/old median DA ratio, largest bin, 1 dp', '1.7', facts['DA_ratio_range'][1])
    check('V.sec.lo', 'app 8.4; README', 'young/old Secondary-fraction ratio, smallest bin, printed 2', '2', facts['secondary_ratio_range'][0], dp=0)
    check('V.sec.hi', 'app 8.4; README', 'young/old Secondary-fraction ratio, largest bin, 1 dp', '3.5', facts['secondary_ratio_range'][1])
    check('V.rho.p', only('A', 'main Sec 2'), 'strict Spearman p below 1e-6 in all three SPS sets', None, bool((stm['resid_T50_spearman_p'] < 1e-6).all()), 'bool')
    check('V.coin', 'README', 'fair-coin tail P(at least 7 of 8), 9/256', '0.035', 9 / 256)
    check('V.mcnoise', 'README', 'Monte Carlo SE of a count p near 0.035 at 5,000 shuffles', '0.003', math.sqrt(9 / 256 * (1 - 9 / 256) / N_SHUFFLE))
    amp = [se.loc[('e1_primary', s_), 'mass_adjusted_diff'] / se.loc[('a_published', s_), 'mass_adjusted_diff'] for s_ in ['int', 'obs', 'both']]
    check('V.prim.amp.lo', 'README', 'Primary-only amplitude as a fraction of the full sample, smallest', '54%', min(amp))
    check('V.prim.amp.hi', 'README', 'Primary-only amplitude as a fraction of the full sample, largest', '63%', max(amp))
    kb = se.loc[('b_kpc', 'obs')], se.loc[('b_kpc', 'both')], se.loc[('b_kpc', 'int')]
    for lab_, r_, ci, pc, pm in [('obs', kb[0], '[+0.011, +0.027]', '0.033', None), ('both', kb[1], '[+0.008, +0.024]', '0.14', None),
                                 ('int', kb[2], '[-0.010, +0.008]', '0.85', '0.52')]:
        ci_check(f'V.B.T1.{lab_}.ci', only('B', 'main T1'), f'kpc controlled residual ({lab_} SPS) mass-adjusted, 3 dp', ci, r_['mass_adjusted_ci_lo'], r_['mass_adjusted_ci_hi'], 'seeded')
        check(f'V.B.T1.{lab_}.pc', only('B', 'main T1'), f'kpc controlled residual ({lab_} SPS): count p as printed', pc, r_['p_count'], 'seeded')
        if pm:
            check(f'V.B.T1.{lab_}.pm', only('B', 'main T1', 'main Sec 2; app 8.4'), f'kpc controlled residual ({lab_} SPS): magnitude p as printed', pm, r_['p_magnitude'], 'seeded')
    sbin = shb[(shb['outcome'] == 'sp_ML_int_Re') & (shb['mass_bin'] == 1)]['median_diff_young_minus_old'].iloc[0]
    check('V.sps.int.bin1', if_printed('differs by only -0.002', 'main Sec 2', 'package only'), 'SPS intrinsic M*/L: young-old median difference in mass bin index 1', '-0.002', sbin)
    check('Z.dedup.N', 'main Sec 1 (Sample)', 'one cube per MaNGA ID: N', '5875', dedup['N_one_cube'], 'count')
    check('Z.dedup.bins', 'app 4.2', 'one cube per MaNGA ID: every Table 1 bin count unchanged', None, dedup['bin_counts_unchanged'], 'bool')
    check('Z.dedup.madj', 'README', 'one cube per MaNGA ID: largest change in a Table 1 mass-adjusted difference', '0.002',
          dedup['max_abs_change_mass_adjusted'])
    # ---------------- the light-weighted side by set, the magnitude test, the subsamples and the T50 tie rule (round-2 verification)
    WR = 'main abstract, Sec 2; app 8.4; README'
    tr = pd.read_csv(out / 'sensitivity' / 'tie_rule_redraws.csv')

    def trow(size, sps, met='sp_MW_Metal_Re', fl='Eps_MGE', sample='all'):
        r = tr[(tr['sample'] == sample) & (tr['size'] == size) & (tr['metallicity'] == met) & (tr['flattening'] == fl) & (tr['SPS_set'] == sps)]
        assert len(r) == 1, (size, sps, met, fl, sample)
        return r.iloc[0]

    da_g = gr['size'] == 'arcsec + logDA'
    check('W.lw.da.obsboth', 'app 8.4; README', 'light-weighted metallicity with log DA: observed and both SPS give 6/8 with either flattening term', None,
          bool((gr[da_g & lw & gr['SPS_set'].isin(['obs', 'both'])]['young_gt_old_bins'] == 6).all()) and int((da_g & lw & gr['SPS_set'].isin(['obs', 'both'])).sum()) == 4, 'bool')
    lwi = gr[da_g & lw & (gr['SPS_set'] == 'int')]['young_gt_old_bins']
    check('W.lw.da.int.lo', 'app 8.4; README', 'light-weighted metallicity with log DA: intrinsic SPS, fewer bins of the two flattening variants', '5', int(lwi.min()), 'count')
    check('W.lw.da.int.hi', 'app 8.4; README', 'light-weighted metallicity with log DA: intrinsic SPS, more bins of the two flattening variants', '6', int(lwi.max()), 'count')
    fr = [trow('arcsec + logDA', 'int', 'sp_LW_Metal_Re', f)['mass_adjusted_diff_mass_tie_rule'] /
          trow('arcsec (published)', 'int', 'sp_LW_Metal_Re', f)['mass_adjusted_diff_mass_tie_rule'] for f in ('Eps_MGE', 'nsa_sersic_ba')]
    check('W.lw.da.int.frac', 'app 8.4; README', 'light-weighted metallicity with log DA: intrinsic-SPS mass-adjusted residual as a fraction of its angular-size value, '
          'printed "a tenth to a quarter" (smaller fraction 0.075 to 0.125, larger 0.20 to 0.30)', None,
          0.075 <= min(fr) < 0.125 and 0.20 <= max(fr) < 0.30, 'bool', note=f'MGE ellipticity {fr[0]:.3f}, NSA axis ratio {fr[1]:.3f}')
    tda = tr[(tr['sample'] == 'all') & (tr['size'] == 'arcsec + logDA')]
    check('W.da.max6', WR, 'no distance-controlled grid variant reaches 7/8, with ties broken by stellar mass or at random', None,
          int(gr[dist]['young_gt_old_bins'].max()) == 6 and int(tda['bins_max'].max()) == 6 and len(tda) == 12, 'bool')
    p6 = se[se['young_gt_old_bins'] == 6]['p_count']
    check('W.p6.lo', 'app 8.4; README', 'count p of the 6/8 rows of the sensitivity scoreboard, smallest (2 dp)', '0.13', float(p6.min()))
    check('W.p6.hi', 'app 8.4; README', 'count p of the 6/8 rows of the sensitivity scoreboard, largest (2 dp)', '0.14', float(p6.max()))
    WM = 'main Sec 2; app 8.4; README'
    for s_, pm in [('obs', '0.010'), ('both', '0.055'), ('int', '0.89')]:
        check(f'W.c.{s_}.pm', WM, f'published + logDA ({s_} SPS): magnitude p as printed', pm, se.loc[('c_published_plus_logDA', s_), 'p_magnitude'], 'seeded')
    rc_ = se.loc[('c_published_plus_logDA', 'obs')]
    ci_check('W.c.obs.ci', 'app 8.4; README', 'published + logDA (obs SPS) mass-adjusted, 3 dp', '[+0.001, +0.017]', rc_['mass_adjusted_ci_lo'], rc_['mass_adjusted_ci_hi'], 'seeded')
    rb_ = se.loc[('c_published_plus_logDA', 'both')]
    check('W.c.both.ci0', 'app 8.4', 'published + logDA (both SPS): bootstrap 95% CI crosses zero', None, rb_['mass_adjusted_ci_lo'] < 0 < rb_['mass_adjusted_ci_hi'], 'bool')
    check('W.b.int.pm', only('A', 'main Sec 2; app 8.4', 'main T1; app 8.4'), 'kpc size (int SPS): magnitude p as printed', '0.52', se.loc[('b_kpc', 'int'), 'p_magnitude'], 'seeded')
    check('W.b.floor', only('A', 'main Sec 2; app 8.4', 'main T1; app 8.4'), 'kpc size: observed- and both-SPS magnitude p at the 0.0002 floor', None,
          all(abs(se.loc[('b_kpc', s_), 'p_magnitude'] - FLOOR) < 1e-12 for s_ in ('obs', 'both')), 'bool')
    e2 = se.loc[('e2_secondary', 'int')]
    ci_check('W.e2.int.ci', 'app 8.4; README', 'Secondary only (int SPS) mass-adjusted, 3 dp', '[-0.033, -0.008]', e2['mass_adjusted_ci_lo'], e2['mass_adjusted_ci_hi'], 'seeded')
    for s_, ci in [('int', '[-0.042, -0.014]'), ('both', '[-0.027, -0.004]')]:
        e4 = se.loc[('e4_secondary_kpc', s_)]
        ci_check(f'W.e4.{s_}.ci', 'app 8.4' + ('; README' if s_ == 'int' else ''), f'Secondary only, size in kpc ({s_} SPS) mass-adjusted, 3 dp', ci,
                 e4['mass_adjusted_ci_lo'], e4['mass_adjusted_ci_hi'], 'seeded')
    check('W.sub.kpc.max4', 'main Sec 2; app 8.4; README', 'MaNGA Primary or Secondary sample alone, size in kpc: most bins young > old in any SPS set', '4',
          max(int(se.loc[(v, s_), 'young_gt_old_bins']) for v in ('e3_primary_kpc', 'e4_secondary_kpc') for s_ in ('int', 'obs', 'both')), 'count')
    check('W.prim.linear', 'main Sec 2; app 8.4; README', 'linear T50 term significant (HC3 p < 0.001) in all nine MaNGA Primary fits (angular size, kpc, plus log DA)', None,
          all(se.loc[(v, s_), 'T50_coef_raw_p_HC3'] < 1e-3 for v in ('e1_primary', 'e3_primary_kpc', 'e5_primary_plus_logDA') for s_ in ('int', 'obs', 'both')), 'bool')
    check('W.sec.neg', 'app 8.4', 'MaNGA Secondary only: young-old difference negative in all three SPS sets with either size unit', None,
          all(se.loc[(v, s_), 'mass_adjusted_diff'] < 0 for v in ('e2_secondary', 'e4_secondary_kpc') for s_ in ('int', 'obs', 'both')), 'bool')
    check('W.e4.both.rev', 'app 8.4', 'e4_secondary_kpc (both SPS): young-old difference reversed (bootstrap 95% CI below zero)', None,
          se.loc[('e4_secondary_kpc', 'both'), 'mass_adjusted_ci_hi'] < 0, 'bool')
    # the 12-variant grid as printed in appendix 8.4 (bins young > old; observed, both, intrinsic SPS)
    GRID_PRINTED = {('sp_MW_Metal_Re', 'Eps_MGE'): {'arcsec (published)': (7, 7, 7), 'kpc': (7, 6, 3), 'logDA': (6, 3, 3)},
                    ('sp_MW_Metal_Re', 'nsa_sersic_ba'): {'arcsec (published)': (7, 7, 7), 'kpc': (6, 6, 3), 'logDA': (6, 6, 3)},
                    ('sp_LW_Metal_Re', 'Eps_MGE'): {'arcsec (published)': (7, 7, 7), 'kpc': (7, 7, 7), 'logDA': (6, 6, 5)},
                    ('sp_LW_Metal_Re', 'nsa_sersic_ba'): {'arcsec (published)': (7, 7, 7), 'kpc': (7, 7, 7), 'logDA': (6, 6, 6)}}
    gi = gr.set_index(['size', 'metallicity', 'flattening', 'SPS_set'])['young_gt_old_bins']
    for (met, fl), cols in GRID_PRINTED.items():
        for size, trip in cols.items():
            for sz in (['arcsec + logDA', 'kpc + logDA'] if size == 'logDA' else [size]):
                for s_, b in zip(('obs', 'both', 'int'), trip):
                    check(f'G.{sz}.{met}.{fl}.{s_}', 'app 8.4 grid table', f'12-variant grid, size {sz}, {met}, {fl} ({s_} SPS): bins young > old', str(b),
                          int(gi[(sz, met, fl, s_)]), 'count')
    for v in ['e2_secondary', 'e4_secondary_kpc']:
        bins = [int(se.loc[(v, s_), 'young_gt_old_bins']) for s_ in ('int', 'obs', 'both')]
        check(f'W.{v}.range', 'app 8.4', f'{v}: the three SPS sets give 2/8 to 4/8', None, min(bins) == 2 and max(bins) == 4, 'bool')
        check(f'W.{v}.int.rev', WM, f'{v} (int SPS): young-old difference reversed (bootstrap 95% CI below zero)', None,
              se.loc[(v, 'int'), 'mass_adjusted_ci_hi'] < 0, 'bool')
        check(f'W.{v}.int.t', 'main Sec 2; app 8.4', f'{v} (int SPS): linear T50 term not significant (HC3 p > 0.05)', None,
              se.loc[(v, 'int'), 'T50_coef_raw_p_HC3'] > 0.05, 'bool')
    WT = 'app 8.4; README (Random seeds)'
    ang = [trow('arcsec (published)', s_) for s_ in ('int', 'obs', 'both')]
    hist = [dict((int(a), int(b)) for a, b in (x.split(':') for x in r['bins_histogram'].split(';'))) for r in ang]
    check('W.tie.ang.7', WT, 'angular size, strict specification: random tie draws at 7/8, of 600 over the three SPS sets', '598', sum(h.get(7, 0) for h in hist), 'count')
    check('W.tie.ang.8', WT, 'angular size, strict specification: the other draws are at 8/8', None, sum(h.get(8, 0) for h in hist) == 2 and all(set(h) <= {7, 8} for h in hist), 'bool')
    tk = tr[(tr['sample'] == 'all') & (tr['size'] == 'kpc')]
    check('W.tie.kpc', WT, 'kpc size: every grid count moves by at most one bin under random ties', None, len(tk) == 12 and
          bool(((tk['bins_max'] - tk['bins_min']) <= 1).all() and ((tk['bins_max'] - tk['bins_mass_tie_rule']).abs() <= 1).all()
               and ((tk['bins_min'] - tk['bins_mass_tie_rule']).abs() <= 1).all()), 'bool')
    kb7 = dict((int(a), int(b)) for a, b in (x.split(':') for x in trow('kpc', 'both')['bins_histogram'].split(';')))
    check('W.tie.kpc.both7', only('B', 'README (Random seeds)', 'package only'), 'kpc size, strict both SPS: random tie draws at 7/8 (of 200)', '15', kb7.get(7, 0), 'count')
    check('W.tie.kpc.fixed', only('B', 'README (Random seeds)', 'package only'), 'kpc size, strict observed and intrinsic SPS: the count never moves', None,
          all(trow('kpc', s_)['bins_min'] == trow('kpc', s_)['bins_max'] == trow('kpc', s_)['bins_mass_tie_rule'] for s_ in ('obs', 'int')), 'bool')
    for s_, lo, med, hi in [('obs', '4', '5', '6'), ('both', '3', '4', '6')]:
        r_ = trow('arcsec + logDA', s_)
        check(f'W.tie.da.{s_}.lo', WT, f'published + logDA ({s_} SPS): fewest bins under random ties', lo, r_['bins_min'], 'count')
        check(f'W.tie.da.{s_}.med', WT, f'published + logDA ({s_} SPS): median bins under random ties', med, r_['bins_median'], 'count')
        check(f'W.tie.da.{s_}.hi', WT, f'published + logDA ({s_} SPS): most bins under random ties', hi, r_['bins_max'], 'count')
    dev = lambda r_: max(abs(r_['mass_adjusted_min'] - r_['mass_adjusted_diff_mass_tie_rule']), abs(r_['mass_adjusted_max'] - r_['mass_adjusted_diff_mass_tie_rule']))
    dmax = max(dev(trow('arcsec + logDA', s_)) for s_ in ('int', 'obs', 'both'))
    check('W.tie.da.madj', WT, 'published + logDA: largest move of a mass-adjusted residual under random ties, printed "at most 0.004"', None, dmax <= 0.004, 'bool', note=f'{dmax:.4f}')
    dall = max(dev(r_) for _, r_ in tr.iterrows())
    check('W.tie.all.madj', 'README (Random seeds)', 'every redrawn row: largest move of a mass-adjusted difference, printed "at most 0.006"', None, dall <= 0.006, 'bool', note=f'{dall:.4f}')
    pr = [trow('arcsec (published)', s_, sample='primary') for s_ in ('int', 'obs', 'both')]
    check('W.tie.prim.lo', 'README (Random seeds)', 'Primary only, angular size: fewest bins under random ties', '4', min(r_['bins_min'] for r_ in pr), 'count')
    check('W.tie.prim.hi', 'README (Random seeds)', 'Primary only, angular size: most bins under random ties', '6', max(r_['bins_max'] for r_ in pr), 'count')
    check('W.tie.mass', 'package files', 'tie-rule file: its stellar-mass tie-rule counts equal the scoreboard and grid rows', None,
          all(int(trow(sz, s_)['bins_mass_tie_rule']) == int(se.loc[(v, s_), 'young_gt_old_bins'])
              for sz, v in (('arcsec (published)', 'a_published'), ('kpc', 'b_kpc'), ('arcsec + logDA', 'c_published_plus_logDA')) for s_ in ('int', 'obs', 'both'))
          and all(int(trow(sz, s_, sample=smp)['bins_mass_tie_rule']) == int(se.loc[(v, s_), 'young_gt_old_bins'])
                  for smp, sz, v in (('primary', 'arcsec (published)', 'e1_primary'), ('primary', 'kpc', 'e3_primary_kpc'), ('primary', 'arcsec + logDA', 'e5_primary_plus_logDA'),
                                     ('secondary', 'arcsec (published)', 'e2_secondary'), ('secondary', 'kpc', 'e4_secondary_kpc'), ('secondary', 'arcsec + logDA', 'e6_secondary_plus_logDA'))
                  for s_ in ('int', 'obs', 'both')), 'bool')
    oc = pd.read_csv(out / 'sensitivity' / 'one_cube_per_mangaid.csv')
    kp = oc[oc['row'].str.contains('kpc')]
    check('Z.dedup.kpc', 'README; app 4.2', 'one cube per MaNGA ID: the kpc controlled-residual rows are present for both samples '
          'and keep their bin counts (3/8, 7/8, 6/8)', None,
          len(kp) == 6 and kp.groupby('row')['young_gt_old_bins'].nunique().max() == 1
          and sorted(kp[kp['sample'] == 'full sample']['young_gt_old_bins']) == [3, 6, 7], 'bool')


def main(argv=None):
    ap = argparse.ArgumentParser(description='Paper 10: regenerate every output and check every analysis number the v3 paper quotes.')
    ap.add_argument('--outdir', default=str(PKG / 'reproduce_outputs'))
    ap.add_argument('--from-public', action='store_true', help='rebuild the joined table from the public files first')
    ap.add_argument('--jam', help='SDSSDR17_MaNGA_JAM_v2.fits (or set P10_JAM)')
    ap.add_argument('--sfh', help='DynPop2_SP_SFH_v2.hdf5.zip or the unzipped .hdf5 (or set P10_SFH)')
    ap.add_argument('--cvc', help='SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt (or set P10_CVC)')
    ap.add_argument('--skip-md5', action='store_true')
    ap.add_argument('--skip-run', action='store_true', help='only re-check outputs already in --outdir')
    args = ap.parse_args(argv)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    detect_text()
    print('paper text next to this script:', {'A': 'angular-size wording', 'B': 'physical-size wording'}.get(TEXT['option'], 'not a v3 text'), flush=True)

    if args.from_public:
        if not args.skip_run:
            # forward only the paths that were given, so build_merged_table.py's P10_JAM/P10_SFH/P10_CVC defaults still apply
            given = {'jam': args.jam, 'sfh': args.sfh, 'cvc': args.cvc}
            missing = [k for k, v in given.items() if not v and not os.environ.get('P10_' + k.upper())]
            if missing:
                ap.error('--from-public needs ' + ', '.join('--' + k for k in missing) + ' (or the P10_JAM, P10_SFH and P10_CVC environment variables)')
            cmd = [x for k, v in given.items() if v for x in ('--' + k, v)] + ['--outdir', out / 'data_rebuilt']
            run(PKG / 'build_merged_table.py', *(cmd + (['--skip-md5'] if args.skip_md5 else [])))
        sys.dont_write_bytecode = True
        sys.path.insert(0, str(PKG))
        import build_merged_table as bmt
        ok, rows = bmt.compare(out / 'data_rebuilt')
        res = pd.DataFrame(rows)
        n_id = int((res['status'] == 'identical').sum())
        n_ulp = int(res['status'].str.startswith('within').sum())
        # contract: every one of the 54 columns and the DA table identical or within 2 ulps, same NaN pattern.
        # The counts are information only: which columns differ by an ulp depends on the platform's log10/log.
        RESULTS.append({'id': 'B.rebuild', 'where': 'data/', 'quantity': 'public-file rebuild vs shipped joined table and DA table',
                        'canonical': 'all 54 columns and the DA table identical or within 2 ulps, same NaN pattern',
                        'regenerated': f'{n_id} identical, {n_ulp} within 2 ulps', 'archived_v06': '', 'kind': 'file', 'tolerance': '2 ulps',
                        'status': 'PASS' if ok else 'FAIL', 'v2_text': '',
                        'note': 'reference platform: 53 identical, 2 within 2 ulps (logRe, cvc_slope_rmax_Re); the build asserts 10,296 rows '
                        'in each file and row-order identity 1.000 for both joins; '
                        + '; '.join(f"{r['column']}: {r['status']}" for r in rows if r['status'] != 'identical')})
        table, da = out / 'data_rebuilt' / TABLE, out / 'data_rebuilt' / DA_TABLE
    else:
        table, da = PKG / 'data' / TABLE, PKG / 'data' / DA_TABLE
        for name, path in [(TABLE, table), (DA_TABLE, da)]:
            got, note = sha256(path), ''
            if got != SHIPPED_SHA256[name]:
                raw = Path(path).read_bytes()
                if b'\r\n' in raw and hashlib.sha256(raw.replace(b'\r\n', b'\n')).hexdigest() == SHIPPED_SHA256[name]:
                    got = SHIPPED_SHA256[name]
                    note = ('CRLF checkout: the LF-normalized bytes match the committed file. Add `paper10/** -text` to the repository '
                            '.gitattributes (or check out with core.autocrlf=false) so the files keep their committed bytes. The run manifests that the '
                            'analysis scripts write then record the CRLF hashes; the numeric outputs are unchanged')
                    print(f'NOTE: {name} has CRLF line endings; its LF-normalized SHA-256 matches', flush=True)
            RESULTS.append({'id': f'B.sha.{name}', 'where': 'data/', 'quantity': f'SHA-256 of shipped {name}', 'canonical': SHIPPED_SHA256[name][:16] + '...',
                            'regenerated': got[:16] + '...', 'archived_v06': '', 'kind': 'file', 'tolerance': 'exact (LF bytes)',
                            'status': 'PASS' if got == SHIPPED_SHA256[name] else 'FAIL', 'v2_text': '', 'note': note})

    if not args.skip_run:
        run(STRICT / 'run_clean_controls_fast.py', '--input', table, '--outdir', out / 'strict')
        run(PKG / '02_sledgehammer_rerun' / 'run_sledgehammer.py', '--input', table, '--outdir', out / 'sledgehammer')
        run(PKG / '03_size_distance_sensitivity' / 'run_size_distance_sensitivity.py', '--input', table, '--da', da,
            '--outdir', out / 'sensitivity')

    # regenerated files vs archived / reference outputs
    for f in sorted(STRICT.glob('manga_paper10_strict_control_*.csv')):
        compare_csv(f'A.strict.{f.stem[29:]}', f'strict-control output regenerates the archived {f.name}', f, out / 'strict' / f.name)
    for f in sorted(SLEDGE_REF.glob('*.csv')):
        compare_csv(f'A.sledge.{f.stem}', f'sledgehammer output matches reference {f.name}', f, out / 'sledgehammer' / f.name)
    for f in sorted(SENS_REF.glob('*.csv')):
        compare_csv(f'A.sens.{f.stem}', f'sensitivity output matches reference {f.name}', f, out / 'sensitivity' / f.name)
    cmp = pd.read_csv(out / 'sledgehammer' / 'sledgehammer_vs_v06.csv')
    det = cmp[cmp['kind'] == 'deterministic']
    RESULTS.append({'id': 'A.v06.det', 'where': '01_clean_pipeline_v0_6/', 'quantity': 'deterministic v0.6 sledgehammer values (counts, medians, '
                    'mass-adjusted differences, pair wins, per-bin differences) reproduce exactly', 'canonical': f'{len(det)} values',
                    'regenerated': f'max |diff| {det["abs_diff"].max():.2g}', 'archived_v06': '', 'kind': 'file', 'tolerance': '1e-12',
                    'status': 'PASS' if det['abs_diff'].max() <= 1e-12 else 'FAIL', 'v2_text': '', 'note': ''})

    number_checks(out)

    res = pd.DataFrame(RESULTS)
    res.to_csv(out / 'reproduce_report.csv', index=False, lineterminator='\n')
    pd.set_option('display.width', 250)
    pd.set_option('display.max_colwidth', 70)
    print()
    print(res[['id', 'status', 'canonical', 'regenerated', 'kind', 'quantity']].to_string(index=False))
    errata = res[(res['v2_text'] != 'same') & (res['v2_text'] != '')]
    n_pass, n_fail = int((res['status'] == 'PASS').sum()), int((res['status'] == 'FAIL').sum())
    lines = ['# Paper 10 reproduction report', '', f'{n_pass} PASS, {n_fail} FAIL, {len(res)} checks.', '',
             '## Values the v2 text printed differently (corrected in v3)', '',
             'These are errata in the v2 text (Zenodo 21478849) that v3 corrects; each check tests the canonical value. '
             'Locations are v3 sections.', '']
    lines += [f'- `{r.id}` ({r.where}): {r.quantity}. Canonical {r.canonical}; {r.v2_text}' for r in errata.itertuples()]
    lines += ['', '## All checks', '', res.to_markdown(index=False)]
    with open(out / 'reproduce_report.md', 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f'\nv2 errata corrected in v3: {len(errata)} (listed in {out / "reproduce_report.md"})')
    if n_fail:
        print('\nFAILED checks:')
        print(res[res['status'] == 'FAIL'][['id', 'canonical', 'regenerated', 'archived_v06', 'tolerance', 'note']].to_string(index=False))
    print(f'\nTOTAL: {n_pass} PASS, {n_fail} FAIL ({len(res)} checks)')
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
