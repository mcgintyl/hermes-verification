import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
import statsmodels.api as sm

IN = Path('/mnt/data/manga_dynpop_merged_thin_firstpass.csv')
OUTDIR = Path('/mnt/data')
SEED=20260629
N_BOOT=5000
N_SHUFFLE=5000
N_BINS=8
K=186
METAL_COL='sp_MW_Metal_Re'
FLAT_COL='Eps_MGE'


def make_base(df, controls):
    cols=['Qual','sp_T50','nsa_sersic_mass','mfl_cyl_log_ML_dyn']+controls
    return df[(df['Qual']>=1) & np.isfinite(df[cols]).all(axis=1)].copy().reset_index(drop=True)

def add_mass_bins(d):
    d=d.copy().reset_index(drop=False).rename(columns={'index':'orig_order'})
    d=d.sort_values(['nsa_sersic_mass','orig_order'], kind='mergesort').reset_index(drop=True)
    per=len(d)//N_BINS
    assert per*N_BINS == len(d), len(d)
    d['mass_bin']=np.repeat(np.arange(N_BINS), per)
    return d

def fit_controls_resid(d, controls):
    X=d[controls].to_numpy(float); y=d['mfl_cyl_log_ML_dyn'].to_numpy(float)
    lr=LinearRegression().fit(X,y)
    d=d.copy(); d['resid']=y-lr.predict(X)
    # age coefficient standardized y/X with HC3
    Xz=d[controls+['sp_T50']].copy()
    Xz=(Xz-Xz.mean())/Xz.std(ddof=0)
    yz=(d['mfl_cyl_log_ML_dyn']-d['mfl_cyl_log_ML_dyn'].mean())/d['mfl_cyl_log_ML_dyn'].std(ddof=0)
    res=sm.OLS(yz, sm.add_constant(Xz)).fit(cov_type='HC3')
    # cv
    kf=KFold(n_splits=5, shuffle=True, random_state=SEED)
    Xc=X; Xa=d[controls+['sp_T50']].to_numpy(float); yy=y
    deltas=[]; r2c=[]; r2a=[]
    for tr,te in kf.split(Xc):
        m1=LinearRegression().fit(Xc[tr],yy[tr]); m2=LinearRegression().fit(Xa[tr],yy[tr])
        a=r2_score(yy[te],m1.predict(Xc[te])); b=r2_score(yy[te],m2.predict(Xa[te]))
        r2c.append(a); r2a.append(b); deltas.append(b-a)
    rho,rp=spearmanr(d['sp_T50'], d['resid'])
    info={
        'N':len(d),
        'controls':', '.join(controls),
        'age_T50_standardized_coef':float(res.params['sp_T50']),
        'age_T50_robust_se':float(res.bse['sp_T50']),
        'age_T50_robust_t':float(res.tvalues['sp_T50']),
        'age_T50_robust_p':float(res.pvalues['sp_T50']),
        'model_R2':float(res.rsquared),
        'cv_R2_controls_only_mean':float(np.mean(r2c)),
        'cv_R2_controls_plus_age_mean':float(np.mean(r2a)),
        'cv_delta_R2_from_age':float(np.mean(deltas)),
        'cv_delta_R2_sd':float(np.std(deltas,ddof=1)),
        'resid_T50_spearman_rho':float(rho),
        'resid_T50_spearman_p':float(rp),
    }
    return d, info

def get_bin_arrays(d, outcome='resid'):
    arrays=[]
    for b in range(N_BINS):
        g=d[d['mass_bin']==b].copy()
        g=g.reset_index(drop=True)
        order=np.lexsort((np.arange(len(g)), g['sp_T50'].to_numpy()))
        yi=order[:K]; oi=order[-K:]
        arrays.append({
            'b':b,
            'vals':g[outcome].to_numpy(float),
            'ages':g['sp_T50'].to_numpy(float),
            'young_idx':yi,
            'old_idx':oi,
            'logM_min':float(g['nsa_sersic_mass'].min()),
            'logM_max':float(g['nsa_sersic_mass'].max()),
            'logM_median':float(g['nsa_sersic_mass'].median()),
            'T50_young_median':float(np.median(g['sp_T50'].to_numpy()[yi])),
            'T50_old_median':float(np.median(g['sp_T50'].to_numpy()[oi])),
            'bin_median':float(np.median(g[outcome]))
        })
    return arrays

def eval_sledge(d, outcome='resid', seed=SEED):
    rng=np.random.default_rng(seed)
    arrays=get_bin_arrays(d,outcome)
    per_rows=[]
    boot_bin_diffs=np.zeros((N_BOOT,N_BINS))
    pooled_y=[]; pooled_o=[]
    for j,a in enumerate(arrays):
        vals=a['vals']; yi=a['young_idx']; oi=a['old_idx']
        y=vals[yi]; o=vals[oi]
        diff=float(np.median(y)-np.median(o))
        # vectorized bootstrap
        y_samp=y[rng.integers(0,len(y),size=(N_BOOT,len(y)))]
        o_samp=o[rng.integers(0,len(o),size=(N_BOOT,len(o)))]
        bdiff=np.median(y_samp,axis=1)-np.median(o_samp,axis=1)
        boot_bin_diffs[:,j]=bdiff
        lo,hi=np.percentile(bdiff,[2.5,97.5])
        per_rows.append({
            'mass_bin':a['b'], 'N_bin':len(vals), 'N_young':len(y), 'N_old':len(o),
            'logM_min':a['logM_min'], 'logM_max':a['logM_max'], 'logM_median':a['logM_median'],
            'T50_young_median':a['T50_young_median'], 'T50_old_median':a['T50_old_median'],
            'young_resid_median':float(np.median(y)), 'old_resid_median':float(np.median(o)),
            'median_diff_young_minus_old':diff, 'bootstrap95_lo':float(lo), 'bootstrap95_hi':float(hi),
            'young_gt_old': bool(diff>0)
        })
        # mass adjusted by subtract bin median
        pooled_y.append(y-a['bin_median']); pooled_o.append(o-a['bin_median'])
    per=pd.DataFrame(per_rows)
    eq_mean=float(per['median_diff_young_minus_old'].mean())
    eq_lo, eq_hi=np.percentile(boot_bin_diffs.mean(axis=1),[2.5,97.5])
    # mass adjusted pooled diff + bootstrap
    py=np.concatenate(pooled_y); po=np.concatenate(pooled_o)
    madiff=float(np.median(py)-np.median(po))
    py_s=py[rng.integers(0,len(py),size=(N_BOOT,len(py)))]
    po_s=po[rng.integers(0,len(po),size=(N_BOOT,len(po)))]
    ma_boot=np.median(py_s,axis=1)-np.median(po_s,axis=1)
    malo,mahi=np.percentile(ma_boot,[2.5,97.5])
    # shuffle null
    counts=np.empty(N_SHUFFLE,dtype=int); sums=np.empty(N_SHUFFLE)
    for s in range(N_SHUFFLE):
        c=0; sm=0.0
        for a in arrays:
            vals=a['vals']; ages=a['ages']
            shuffled_age=rng.permutation(ages)
            order=np.lexsort((np.arange(len(vals)), shuffled_age))
            y=vals[order[:K]]; o=vals[order[-K:]]
            diff=float(np.median(y)-np.median(o))
            c += diff>0; sm += diff
        counts[s]=c; sums[s]=sm
    real_count=int((per['median_diff_young_minus_old']>0).sum())
    real_sum=float(per['median_diff_young_minus_old'].sum())
    pcount=(np.sum(counts>=real_count)+1)/(N_SHUFFLE+1)
    psum=(np.sum(sums>=real_sum)+1)/(N_SHUFFLE+1)
    pjoint=(np.sum((counts>=real_count)&(sums>=real_sum))+1)/(N_SHUFFLE+1)
    score={
        'N':len(d), 'young_gt_old_bins':real_count, 'total_bins':N_BINS,
        'equal_bin_mean_diff_young_minus_old':eq_mean,
        'equal_bin_mean_boot95_lo':float(eq_lo), 'equal_bin_mean_boot95_hi':float(eq_hi),
        'sum_of_8_bin_diffs':real_sum,
        'mass_adjusted_diff_young_minus_old':madiff,
        'mass_adjusted_diff_boot95_lo':float(malo), 'mass_adjusted_diff_boot95_hi':float(mahi),
        'mass_adjusted_young_median':float(np.median(py)), 'mass_adjusted_old_median':float(np.median(po)),
        'shuffle_p_count_ge_real':float(pcount),
        'shuffle_p_sum_ge_real':float(psum),
        'shuffle_p_joint_count_and_sum_ge_real':float(pjoint),
        'shuffle_count_mean':float(np.mean(counts)), 'shuffle_count_sd':float(np.std(counts,ddof=1)),
        'shuffle_sum_mean':float(np.mean(sums)), 'shuffle_sum_sd':float(np.std(sums,ddof=1)),
    }
    return score,per

def main():
    df=pd.read_csv(IN)
    specs={
        'strict_intSPS_MW_Eps_logSigma':['sp_ML_int_Re',METAL_COL,'nsa_sersic_mass','logRe','nsa_sersic_n',FLAT_COL,'Lambda_Re','logSigma_Re'],
        'strict_obsSPS_MW_Eps_logSigma':['sp_ML_obs_Re',METAL_COL,'nsa_sersic_mass','logRe','nsa_sersic_n',FLAT_COL,'Lambda_Re','logSigma_Re'],
        'strict_bothSPS_MW_Eps_logSigma':['sp_ML_int_Re','sp_ML_obs_Re',METAL_COL,'nsa_sersic_mass','logRe','nsa_sersic_n',FLAT_COL,'Lambda_Re','logSigma_Re']
    }
    scoreboard=[]; perbins=[]; models=[]
    for label,controls in specs.items():
        print('running',label, flush=True)
        base=make_base(df,controls)
        resid,info=fit_controls_resid(base,controls)
        binned=add_mass_bins(resid)
        score,per=eval_sledge(binned,'resid',seed=SEED+len(scoreboard)*1000)
        score['control_set']=label
        per.insert(0,'control_set',label)
        info['control_set']=label
        scoreboard.append(score); perbins.append(per); models.append(info)
    scoredf=pd.DataFrame(scoreboard)
    scoredf=scoredf[['control_set']+[c for c in scoredf.columns if c!='control_set']]
    perdf=pd.concat(perbins,ignore_index=True)
    modeldf=pd.DataFrame(models)
    scoredf.to_csv(OUTDIR/'manga_paper10_strict_control_residual_scoreboard.csv',index=False)
    perdf.to_csv(OUTDIR/'manga_paper10_strict_control_residual_perbin.csv',index=False)
    modeldf.to_csv(OUTDIR/'manga_paper10_strict_control_model_terms.csv',index=False)
    # Quick sensitivity grid: no boot/shuffle.
    sens=[]
    for met in ['sp_MW_Metal_Re','sp_LW_Metal_Re']:
        for fl in ['Eps_MGE','nsa_sersic_ba']:
            for lab,sps in [('int',['sp_ML_int_Re']),('obs',['sp_ML_obs_Re']),('both',['sp_ML_int_Re','sp_ML_obs_Re'])]:
                controls=sps+[met,'nsa_sersic_mass','logRe','nsa_sersic_n',fl,'Lambda_Re','logSigma_Re']
                base=make_base(df,controls)
                resid,info=fit_controls_resid(base,controls)
                binned=add_mass_bins(resid)
                arr=get_bin_arrays(binned,'resid')
                diffs=[]
                for a in arr:
                    vals=a['vals']; diffs.append(float(np.median(vals[a['young_idx']])-np.median(vals[a['old_idx']])))
                sens.append({'metallicity':met,'flattening':fl,'SPS_set':lab,'N':len(binned),'young_gt_old_bins':int(np.sum(np.array(diffs)>0)),'mean_diff':float(np.mean(diffs)),'sum_diff':float(np.sum(diffs)),'age_coef':info['age_T50_standardized_coef'],'cv_delta_R2_from_age':info['cv_delta_R2_from_age'],'resid_spearman_rho':info['resid_T50_spearman_rho']})
    sensdf=pd.DataFrame(sens)
    sensdf.to_csv(OUTDIR/'manga_paper10_strict_control_sensitivity_grid.csv',index=False)
    md=[]
    md.append('# Paper 10 strict collinearity-free control rerun\n')
    md.append('Primary strict specification: mass-weighted metallicity (`sp_MW_Metal_Re`), MGE ellipticity (`Eps_MGE`), and `logSigma_Re` only. The model excludes light-weighted metallicity, duplicate flattening/axis-ratio terms, and raw `Sigma_Re`. Outcome is `mfl_cyl_log_ML_dyn`; age split uses `sp_T50`; mass bins use `nsa_sersic_mass`.\n')
    md.append('Primary controls: SPS M/L term(s), `sp_MW_Metal_Re`, `nsa_sersic_mass`, `logRe`, `nsa_sersic_n`, `Eps_MGE`, `Lambda_Re`, `logSigma_Re`.\n')
    md.append('## Residual sledgehammer scoreboard\n')
    md.append(scoredf.to_markdown(index=False, floatfmt='.6g'))
    md.append('\n## Per-bin residuals\n')
    md.append(perdf.to_markdown(index=False, floatfmt='.6g'))
    md.append('\n## Age terms / CV / Spearman\n')
    md.append(modeldf.to_markdown(index=False, floatfmt='.6g'))
    md.append('\n## Sensitivity grid, no bootstrap/shuffle\n')
    md.append(sensdf.to_markdown(index=False, floatfmt='.6g'))
    (OUTDIR/'manga_paper10_strict_control_rerun_report.md').write_text('\n'.join(md))
    print(scoredf.to_string(index=False))
    print(modeldf[['control_set','age_T50_standardized_coef','age_T50_robust_p','cv_delta_R2_from_age','resid_T50_spearman_rho','resid_T50_spearman_p']].to_string(index=False))

if __name__=='__main__':
    main()
