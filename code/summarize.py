#!/usr/bin/env python3
"""Summarize stored coefficient/covariance arrays; does not import fit code."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy.stats import norm,chi2
Z=norm.ppf(.975)

def crossed_se(a):
    # a: independent population rows x independent shared-label columns.
    M,R=a.shape
    if not np.isfinite(a).all():return np.nan
    if M==1:return np.std(a[0],ddof=1)/np.sqrt(R)
    row=a.mean(1);col=a.mean(0);res=a-row[:,None]-col[None,:]+a.mean()
    ms=np.sum(res**2)/((M-1)*(R-1))
    v=np.var(row,ddof=1)/M+np.var(col,ddof=1)/R-ms/(M*R)
    return np.sqrt(v) if v>=0 else np.nan

def summarize(root):
    meta=json.loads((root/'metadata.json').read_text());pop=np.load(root/'population.npz')
    truth=pop['truth'];meth=meta['methods'];resp=meta['responses'];coefs=meta['coefficients'];R=meta['R'];isreal=meta['dataset']=='lucas'
    rows=[];joint=[];paired=[];num=[]
    for file in sorted(root.glob('*_n*.npz')):
        tag=file.stem;sc,nn=tag.rsplit('_n',1);n=int(nn);d=np.load(file);est=d['est'];nr=est.shape[2];p=est.shape[3]
        tr=np.broadcast_to(truth,(nr,p));errors=est-tr[None,None,:,:]
        for vtype,key in [('linearization','covariance'),('jackknife','jackknife')]:
          allcov=d[key];covmask=np.zeros(est.shape,bool)
          for j,ms in enumerate(meth):
            if vtype=='jackknife' and ms['mode']=='oracle':continue
            for k,coef in enumerate(coefs):
              var=allcov[:,j,:,k,k];ee=errors[:,j,:,k]
              valid=np.isfinite(ee)&np.isfinite(var)&(var>0)
              cover=valid&(np.abs(ee)<=Z*np.sqrt(np.maximum(var,0)));covmask[:,j,:,k]=cover
              sets=[([r],resp[r]) for r in range(nr)] if isreal else [(list(range(nr)),'all_populations')]
              for rset,rname in sets:
                cv=cover[:,rset].T;va=valid[:,rset].T;er=ee[:,rset].T;vr=var[:,rset].T
                valid_est=np.isfinite(er);count=va.sum();attempts=va.size
                sd_perpop=np.nanstd(er,axis=1,ddof=1)
                bias_perpop=np.nanmean(er,axis=1)
                width=np.where(va,2*Z*np.sqrt(np.maximum(vr,0)),np.nan)
                g=d['gap'][j,rset,k];vtrue=d['oracle_covariance'][j,rset,k,k]
                lin=d['linear'][:,j,rset,k].T;rem=d['remainder'][:,j,rset,k].T
                row=dict(dataset=meta['dataset'],scenario=sc,n=n,method=ms['id'],prediction=rname,coefficient=coef,variance=vtype,
                    attempts=int(attempts),finite_intervals=int(count),covered=int(cv.sum()),coverage_all=float(cv.mean()),coverage_finite=float(cv.sum()/count) if count else np.nan,
                    coverage_mcse=float(crossed_se(cv.astype(float))),coverage_population_p10=float(np.quantile(cv.mean(1),.1)),coverage_population_p90=float(np.quantile(cv.mean(1),.9)),
                    bias=float(np.nanmean(er)),bias_rms_across_populations=float(np.sqrt(np.mean(bias_perpop**2))),bias_mcse=float(crossed_se(er)),
                    rmse=float(np.sqrt(np.nanmean(er**2))),mean_mc_sd=float(np.mean(sd_perpop)),
                    mean_se=float(np.nanmean(np.where(va,np.sqrt(np.maximum(vr,0)),np.nan))),mean_width=float(np.nanmean(width)),
                    truth_mean=float(tr[rset,k].mean()),gap_mean=float(g.mean()),gap_rms=float(np.sqrt(np.mean(g**2))),
                    oracle_linear_sd_mean=float(np.mean(np.sqrt(vtrue))),linear_rms=float(np.sqrt(np.nanmean(lin**2))),remainder_rms=float(np.sqrt(np.nanmean(rem**2))),
                    estimate_success=int(valid_est.sum()))
                rows.append(row)
            # Full coefficient and slopes-only ellipsoids: use off-diagonal covariance.
            for group,ix in [('all_coefficients',np.arange(p)),('slopes_only',np.arange(1,p))]:
                mats=allcov[:,j][:,:,ix,:][:,:,:,ix];errs=errors[:,j][:,:,ix]
                flat=mats.reshape(-1,len(ix),len(ix));eflat=errs.reshape(-1,len(ix))
                valid=np.isfinite(flat).all((1,2))&np.isfinite(eflat).all(1);stat=np.full(len(flat),np.nan)
                ids=np.flatnonzero(valid)
                if len(ids):
                    ev=np.linalg.eigvalsh(flat[ids]);v=ev[:,0]>1e-12*ev[:,-1];valid[ids]=v;ids=ids[v]
                    stat[ids]=np.einsum('ni,ni->n',eflat[ids],np.linalg.solve(flat[ids],eflat[ids,...,None])[...,0])
                cover=(valid&(stat<=chi2.ppf(.95,len(ix)))).reshape(R,nr)
                sets=[([r],resp[r]) for r in range(nr)] if isreal else [(list(range(nr)),'all_populations')]
                for rset,rname in sets:
                    cv=cover[:,rset].T;vv=valid.reshape(R,nr)[:,rset].T
                    joint.append(dict(dataset=meta['dataset'],scenario=sc,n=n,method=ms['id'],prediction=rname,group=group,variance=vtype,attempts=int(vv.size),finite_intervals=int(vv.sum()),covered=int(cv.sum()),coverage_all=float(cv.mean()),coverage_finite=float(cv.sum()/vv.sum()) if vv.sum() else np.nan,coverage_mcse=float(crossed_se(cv.astype(float)))))
          comparisons=[('target_map_spatial','map_spatial_score'),('product_score','map_spatial_score'),('soft_target_bc','soft_target'),('soft_score_bc','soft_score'),('soft_target_bc','soft_score_bc')] if isreal else [('target_B3','common_B3'),('product_B5','common_B3')]
          mids=[a['id'] for a in meth]
          for plus,minus in comparisons:
            jp,jm=mids.index(plus),mids.index(minus)
            sets=[([r],resp[r]) for r in range(nr)] if isreal else [(list(range(nr)),'all_populations')]
            for k,coef in enumerate(coefs):
             for rset,rname in sets:
                dd=(covmask[:,jp,rset,k].astype(float)-covmask[:,jm,rset,k].astype(float)).T
                sq=(errors[:,jp,rset,k]**2-errors[:,jm,rset,k]**2).T
                paired.append(dict(scenario=sc,n=n,prediction=rname,coefficient=coef,variance=vtype,plus=plus,minus=minus,coverage_difference=float(dd.mean()),coverage_difference_mcse=float(crossed_se(dd)),mse_difference=float(np.nanmean(sq)),mse_difference_mcse=float(crossed_se(sq))))
        for j,m in enumerate(meth):
            s=d['stats'][:,j];mode=m['mode']
            num.append(dict(scenario=sc,n=n,method=m['id'],attempts=R,fit_success=int(d['success'][:,j].sum()),max_condition_R=float(np.nanmax(s[:,1])),max_leverage=float(np.nanmax(s[:,0])),signed_contrast_weights=(mode=='target'),draws_with_negative_common_weights=int(np.sum(s[:,2]>0)) if mode=='common' else None,max_absolute_weight=float(np.nanmax(s[:,3])) if mode!='oracle' else None,max_total_absolute_weight=float(np.nanmax(s[:,4])) if mode!='oracle' else None))
    pd.DataFrame(rows).to_csv(root/'summary.csv',index=False,float_format='%.17g')
    pd.DataFrame(joint).to_csv(root/'joint_coverage.csv',index=False,float_format='%.17g')
    pd.DataFrame(paired).to_csv(root/'paired_comparisons.csv',index=False,float_format='%.17g')
    pd.DataFrame(num).to_csv(root/'numerical_diagnostics.csv',index=False,float_format='%.17g')
    print(root,'rows',len(rows),'joint',len(joint),'paired',len(paired))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);a=ap.parse_args();summarize(a.root)
