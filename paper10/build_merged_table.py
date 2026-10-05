"""Rebuild the Paper 10 joined working table from the three public MaNGA DynPop files.

Public inputs (Zenodo; MD5 as listed by Zenodo, checked by this script):
  SDSSDR17_MaNGA_JAM_v2.fits               88f7dc25f6cadb2c5a2de4b190d00329   21,479,040 B
      https://zenodo.org/records/17518315/files/SDSSDR17_MaNGA_JAM_v2.fits
  SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt   b7aa4a4f188196c39369f918768aa8c0    1,134,922 B
      https://zenodo.org/records/17518315/files/SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt
      (record 10.5281/zenodo.17518315, DynPop VII dataset, v2, 2026-02-06; use the _v2 FITS,
       not SDSSDR17_MaNGA_JAM.fits)
  DynPop2_SP_SFH_v2.hdf5.zip               870368917b3a7c30f2da7598f21bdb2c  229,323,792 B
      https://zenodo.org/records/15742825/files/DynPop2_SP_SFH_v2.hdf5.zip
      or the unzipped DynPop2_SP_SFH_v2.hdf5 from the same record
                                           543bc4d28dc5d33e1244a7c33bd631b6  2,988,377,528 B
      (record 10.5281/zenodo.15742825, updated DynPop II catalogue, v2, 2025-06-26; the zip
       holds SP_SFH_v2.hdf5, identical to the unzipped file)

Outputs (in --outdir, default data/rebuilt/):
  manga_dynpop_merged_thin_firstpass.csv  10,296 rows x 54 columns, the joined working table
  manga_dynpop_DA_target.csv              plateifu, mangaid, DA (adopted angular-diameter distance,
                                          Mpc) and target (MaNGA subsample: 0 Primary, 1 Secondary,
                                          2 color-enhanced), copied from JAM v2 HDU 1; used by
                                          03_size_distance_sensitivity/

The joined table shipped in data/ (SHA-256 13336cd4...8f51f58) is the working table recovered
from the original analysis. A rebuild matches it in shape, column order and NaN pattern; 52 of 54
columns are bit-identical. logRe differs in 6 of 10,143 finite rows by 1 ulp, and
cvc_slope_rmax_Re (a ratio of two logarithms, unused by any analysis) in 227 of 10,127 finite rows
by 1 ulp (208 rows) or 2 ulps (19 rows); the largest difference is 2.2e-16 absolute, 3.8e-16
relative.
The cause is that log10/log round differently in different platforms' math libraries.
Which columns differ depends on the platform's math library, so the contract that --compare and
reproduce.py --from-public check is: every column identical or within 2 ulps, with an identical
NaN pattern (exit status 1 otherwise), not a particular count of differing columns.

Column notes: logRe = log10(Re_arcsec_MGE), the MGE effective radius in ARCSEC (an angular
size, not kpc); DML_mfl_int/obs = JAM MFL (cylindrical) log(M/L)_dyn minus SPS intrinsic/observed
log(M*/L) within Re; nsa_sersic_n/ba sentinels (-9999) are set to NaN.

JAM v2 HDU map (SDSSDR17_MaNGA_JAM_v2_datamodel.pdf, Table 1): 1 global, 2 JAMcyl+MFL,
4 JAMcyl+NFW, 8 JAMcyl+gNFW (others unused).

Usage:
  python build_merged_table.py --jam SDSSDR17_MaNGA_JAM_v2.fits --sfh DynPop2_SP_SFH_v2.hdf5.zip \
      --cvc SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt [--outdir data/rebuilt] [--compare]
  (paths may also come from the P10_JAM, P10_SFH and P10_CVC environment variables)
"""
import argparse
import hashlib
import os
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SHIPPED = HERE / 'data'
TABLE = 'manga_dynpop_merged_thin_firstpass.csv'
DA_TABLE = 'manga_dynpop_DA_target.csv'

MD5 = {
    'jam': '88f7dc25f6cadb2c5a2de4b190d00329',
    'cvc': 'b7aa4a4f188196c39369f918768aa8c0',
    'sfh_zip': '870368917b3a7c30f2da7598f21bdb2c',
    'sfh_hdf5': '543bc4d28dc5d33e1244a7c33bd631b6',
}
HDU = {'mfl_cyl': 2, 'nfw_cyl': 4, 'gnfw_cyl': 8}
SP_SCALARS = ['T50', 'T90', 'LW_Age_Re', 'MW_Age_Re', 'LW_Metal_Re', 'MW_Metal_Re',
              'ML_int_Re', 'ML_obs_Re', 'Mstar_Re', 'SNR_Re', 'chi2_Re']
# DynPop VII CVC table, CDS byte-by-byte layout (0-based slices)
CVC_COLSPECS = [(0, 11), (12, 21), (22, 29), (30, 36), (37, 43), (44, 49), (50, 56), (57, 63),
                (64, 72), (73, 81), (82, 90), (91, 99), (100, 106), (107, 109)]
CVC_NAMES = ['PlateIFU', 'Mangaid', 'DA', 'Re', 'MajAxis', 'FWHM', 'R_Vcmax', 'rmax',
             'Vc_Re', 'Vc_MajAxis', 'Vcmax', 'Vc_rmax', 'logMBH', 'Qual']


def md5sum(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def check_md5(path, key, skip):
    if skip:
        print(f'md5 not checked for {path}')
        return
    got = md5sum(path)
    if got != MD5[key]:
        sys.exit(f'MD5 mismatch for {path}: got {got}, Zenodo lists {MD5[key]} (use --skip-md5 to override)')
    print(f'md5 OK  {Path(path).name}  {got}')


def read_jam(path):
    from astropy.io import fits
    h = fits.open(path)
    g = h[1].data

    def s(col):
        return np.char.strip(np.asarray(g[col]).astype(str))

    jam = pd.DataFrame({
        'plateifu': s('plateifu'), 'mangaid': s('mangaid'),
        'Qual': np.asarray(g['Qual']).astype(int), 'drp3qual': np.asarray(g['drp3qual']).astype(int),
        'Lambda_Re': np.asarray(g['Lambda_Re'], float), 'Sigma_Re': np.asarray(g['Sigma_Re'], float),
        'Eps_MGE': np.asarray(g['Eps_MGE'], float), 'nsa_sersic_mass': np.asarray(g['nsa_sersic_mass'], float),
        'nsa_sersic_n': np.asarray(g['nsa_sersic_n'], float), 'nsa_sersic_ba': np.asarray(g['nsa_sersic_ba'], float),
    })
    for c in ['nsa_sersic_n', 'nsa_sersic_ba']:          # 53 rows carry the NSA sentinel -9999
        jam.loc[jam[c] <= -999, c] = np.nan
    extra = {
        'Re_arcsec_MGE': np.asarray(g['Re_arcsec_MGE'], float),
        'DA': np.asarray(g['DA'], float),
        'target': np.asarray(g['target']).astype(int),
    }
    models = {hdu: h[hdu].data for hdu in set(HDU.values())}
    return jam, extra, models


def read_sfh(path):
    import h5py
    with h5py.File(path, 'r') as f:
        keys = sorted(f.keys())
        sp = {'sp_plateIFU': [], 'sp_mangaid': []}
        for c in SP_SCALARS:
            sp['sp_' + c] = np.empty(len(keys))
        for i, k in enumerate(keys):
            grp = f[k]
            sp['sp_plateIFU'].append(grp['plateIFU'][0].decode().strip())
            sp['sp_mangaid'].append(grp['mangaid'][0].decode().strip())
            for c in SP_SCALARS:
                sp['sp_' + c][i] = grp[c][()]
    return pd.DataFrame(sp)


def read_cvc(path):
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    start = [i for i, L in enumerate(lines) if L.startswith('-----')][-1] + 1
    rows = [[L.ljust(109)[a:b].strip() for a, b in CVC_COLSPECS] for L in lines[start:] if L.strip()]
    cvc = pd.DataFrame(rows, columns=CVC_NAMES)
    for c in CVC_NAMES[2:]:
        cvc[c] = pd.to_numeric(cvc[c].replace('', np.nan), errors='coerce')
    cvc['Qual'] = cvc['Qual'].astype(int)
    return cvc


def build(jam_path, sfh_hdf5, cvc_path):
    jam, extra, models = read_jam(jam_path)
    sp = read_sfh(sfh_hdf5)
    cvc = read_cvc(cvc_path)
    n = len(jam)
    assert len(sp) == n and len(cvc) == n, (n, len(sp), len(cvc))
    id_sp = float(np.mean((sp['sp_plateIFU'].values == jam['plateifu'].values) & (sp['sp_mangaid'].values == jam['mangaid'].values)))
    id_cvc = float(np.mean((cvc['PlateIFU'].values == jam['plateifu'].values) & (cvc['Mangaid'].values == jam['mangaid'].values)))
    print(f'JAM rows {n}  SP/SFH rows {len(sp)}  CVC rows {len(cvc)}')
    print(f'row-order identity JAM<->SP/SFH {id_sp:.4f}  JAM<->CVC {id_cvc:.4f}')
    assert id_sp == 1.0 and id_cvc == 1.0 and jam['plateifu'].is_unique

    out = jam.copy()
    for c in SP_SCALARS:
        out['sp_' + c] = sp['sp_' + c].values
    for pre, hdu in HDU.items():
        d = models[hdu]
        if pre == 'mfl_cyl':
            out['mfl_cyl_log_ML_dyn'] = np.asarray(d['log_ML_dyn'], float)
            out['mfl_cyl_log_Mt_Re'] = np.asarray(d['log_Mt_Re'], float)
            out['mfl_cyl_chi2_dof'] = np.asarray(d['chi2_dof'], float)
        else:
            for c in ['log_Mt_Re', 'log_Ms_Re', 'log_Md_Re', 'fdm_Re']:
                out[f'{pre}_{c}'] = np.asarray(d[c], float)
            if pre == 'gnfw_cyl':
                out['gnfw_cyl_gamma_gNFW'] = np.asarray(d['gamma_gNFW'], float)
            out[f'{pre}_log_ML_dyn_Re'] = np.asarray(d['log_ML_dyn_Re'], float)
            out[f'{pre}_chi2_dof'] = np.asarray(d['chi2_dof'], float)
    out['cvc_Re'] = cvc['Re'].values
    out['cvc_rmax'] = cvc['rmax'].values
    out['cvc_Vc_Re'] = cvc['Vc_Re'].values
    out['cvc_Vc_rmax'] = cvc['Vc_rmax'].values
    out['cvc_Vcmax'] = cvc['Vcmax'].values
    out['cvc_CVC_Qual'] = cvc['Qual'].values

    out['DML_mfl_int'] = out['mfl_cyl_log_ML_dyn'] - out['sp_ML_int_Re']
    out['DML_mfl_obs'] = out['mfl_cyl_log_ML_dyn'] - out['sp_ML_obs_Re']
    out['Dmass_mfl_Re'] = out['mfl_cyl_log_Mt_Re'] - out['sp_Mstar_Re']
    out['Dmass_gnfw_Re'] = out['gnfw_cyl_log_Mt_Re'] - out['gnfw_cyl_log_Ms_Re']
    out['gnfw_dm_to_stellar_log_Re'] = out['gnfw_cyl_log_Md_Re'] - out['gnfw_cyl_log_Ms_Re']
    out['nfw_dm_to_stellar_log_Re'] = out['nfw_cyl_log_Md_Re'] - out['nfw_cyl_log_Ms_Re']
    with np.errstate(divide='ignore', invalid='ignore'):
        out['logRe'] = np.log10(extra['Re_arcsec_MGE'])      # ARCSEC (angular size), not kpc
        out['logSigma_Re'] = np.log10(out['Sigma_Re'])
        out['Vc_ratio_rmax_Re'] = out['cvc_Vc_rmax'] / out['cvc_Vc_Re']
        out['Vc_ratio_max_Re'] = out['cvc_Vcmax'] / out['cvc_Vc_Re']
        out['cvc_slope_rmax_Re'] = np.log(out['Vc_ratio_rmax_Re']) / np.log(out['cvc_rmax'] / out['cvc_Re'])
    out = out.replace([np.inf, -np.inf], np.nan)
    q = out['Qual'] >= 1
    prim = q & np.isfinite(out[['sp_T50', 'nsa_sersic_mass', 'DML_mfl_int', 'DML_mfl_obs']]).all(axis=1)
    print(f'Qual>=1 rows {int(q.sum())}  primary finite sample {int(prim.sum())}')
    da = pd.DataFrame({'plateifu': jam['plateifu'], 'mangaid': jam['mangaid'], 'DA': extra['DA'], 'target': extra['target']})
    return out, da


def compare(rebuilt_dir, shipped_dir=SHIPPED):
    """Compare rebuilt tables with the shipped ones. Returns (ok, list of per-column rows)."""
    a = pd.read_csv(shipped_dir / TABLE, float_precision='round_trip')
    b = pd.read_csv(rebuilt_dir / TABLE, float_precision='round_trip')
    rows, ok = [], True
    if list(a.columns) != list(b.columns) or a.shape != b.shape:
        return False, [{'column': '(shape/columns)', 'status': f'shipped {a.shape} vs rebuilt {b.shape}'}]
    for c in a.columns:
        x, y = a[c], b[c]
        if x.dtype == object or y.dtype == object:
            same = bool((x.astype(str) == y.astype(str)).all())
            rows.append({'column': c, 'identical_rows': int((x.astype(str) == y.astype(str)).sum()), 'max_ulps': 0 if same else None,
                         'status': 'identical' if same else 'DIFFERENT'})
            ok &= same
            continue
        xv, yv = x.to_numpy(float), y.to_numpy(float)
        nan_same = bool(np.array_equal(np.isnan(xv), np.isnan(yv)))
        m = ~np.isnan(xv)
        eq = int(np.sum(xv[m] == yv[m]) + np.sum(~m))
        ulps = 0
        if eq < len(xv):
            d = np.abs(xv[m] - yv[m])
            ulps = int(np.max(np.ceil(d / np.spacing(np.maximum(np.abs(xv[m]), np.abs(yv[m]))))))
        status = 'identical' if (eq == len(xv) and nan_same) else (f'within {ulps} ulp' if (nan_same and ulps <= 2) else 'DIFFERENT')
        ok &= status != 'DIFFERENT'
        rows.append({'column': c, 'identical_rows': eq, 'max_ulps': ulps, 'status': status})
    da_a = pd.read_csv(shipped_dir / DA_TABLE, float_precision='round_trip')
    da_b = pd.read_csv(rebuilt_dir / DA_TABLE, float_precision='round_trip')
    da_same = da_a.equals(da_b)
    rows.append({'column': f'{DA_TABLE} (all columns)', 'identical_rows': int(len(da_a)) if da_same else None,
                 'max_ulps': 0 if da_same else None, 'status': 'identical' if da_same else 'DIFFERENT'})
    ok &= da_same
    return ok, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description='Rebuild the Paper 10 joined table from the public DynPop files.')
    ap.add_argument('--jam', default=os.environ.get('P10_JAM'), help='SDSSDR17_MaNGA_JAM_v2.fits')
    ap.add_argument('--sfh', default=os.environ.get('P10_SFH'), help='DynPop2_SP_SFH_v2.hdf5.zip or the unzipped .hdf5')
    ap.add_argument('--cvc', default=os.environ.get('P10_CVC'), help='SDSSDR17_MaNGA_gNFW_cyl_Vcirc_ApJS.txt')
    ap.add_argument('--outdir', default=str(SHIPPED / 'rebuilt'))
    ap.add_argument('--extract-dir', default=None, help='where to unzip SP_SFH_v2.hdf5 (2.99 GB); default a temporary '
                    'folder inside --outdir that is deleted afterwards')
    ap.add_argument('--skip-md5', action='store_true')
    ap.add_argument('--compare', action='store_true', help='compare the rebuilt tables with the shipped data/ tables')
    args = ap.parse_args(argv)
    if not (args.jam and args.sfh and args.cvc):
        ap.error('give --jam, --sfh and --cvc (or set P10_JAM, P10_SFH, P10_CVC)')
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    check_md5(args.jam, 'jam', args.skip_md5)
    check_md5(args.cvc, 'cvc', args.skip_md5)
    tmp = None
    sfh = Path(args.sfh)
    if sfh.suffix == '.zip':
        check_md5(sfh, 'sfh_zip', args.skip_md5)
        if args.extract_dir:
            xdir = Path(args.extract_dir)
            xdir.mkdir(parents=True, exist_ok=True)
        else:
            tmp = tempfile.mkdtemp(prefix='sfh_unzip_', dir=outdir)
            xdir = Path(tmp)
        target = xdir / 'SP_SFH_v2.hdf5'
        if not target.exists():
            print(f'unzipping SP_SFH_v2.hdf5 to {xdir} ...', flush=True)
            with zipfile.ZipFile(sfh) as z:
                z.extract('SP_SFH_v2.hdf5', xdir)
        sfh = target
    check_md5(sfh, 'sfh_hdf5', args.skip_md5)
    try:
        table, da = build(args.jam, sfh, args.cvc)
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)
    table.to_csv(outdir / TABLE, index=False, lineterminator='\n')
    da.to_csv(outdir / DA_TABLE, index=False, lineterminator='\n')
    print('wrote', outdir / TABLE, table.shape)
    print('wrote', outdir / DA_TABLE, da.shape)
    if args.compare:
        ok, rows = compare(outdir)
        res = pd.DataFrame(rows)
        print(res.to_string(index=False))
        n_id = int((res['status'] == 'identical').sum())
        n_ulp = int(res['status'].str.startswith('within').sum())
        print(f'compare with shipped data/: {n_id} identical, {n_ulp} within 2 ulps, '
              f'{int((res["status"] == "DIFFERENT").sum())} different -> {"OK" if ok else "FAILED"}')
        if not ok:
            sys.exit(1)


if __name__ == '__main__':
    main()
