"""Portable saved-result reaggregation; no fits or random draws."""
from pathlib import Path
import argparse,numpy as np,pandas as pd
from scipy.stats import norm
ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);ap.add_argument('output',type=Path);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
truth=np.load(args.results/'population.npz')['truth'];rows=[]
for sc in ['uniform','map','map_east']:
 for n in [100,250,500]:
  a=np.load(args.results/f'{sc}_n{n}.npz');R=len(a['labels']);assert np.array_equal(a['sizes'],a['labels'].sum(1))
  for v,pred in enumerate(['zero','WorldCover']):
   for j,coef in [(1,'east_per100km'),(2,'north_per100km')]:
    r=dict(scenario=sc,n=n,prediction=pred,coefficient=coef,attempts=R,label_budget_empirical_mean=float(a['sizes'].mean()),max_pi=float(a['pi'].max()));inds={}
    for tag,m in [('C',3),('T',6)]:
     e=a['est'][:,m,v,j]-truth[j];vv=a['covariance'][:,m,v,j,j];finite=np.isfinite(e);valid=finite&np.isfinite(vv)&(vv>0)&a['success'][:,m]
     cover=valid&(abs(e)<=norm.ppf(.975)*np.sqrt(np.maximum(vv,0)));inds[tag]=cover.astype(int);V=float(a['oracle_covariance'][m,v,j,j]);gap=float(a['gap'][m,v,j])
     r.update({f'coverage_{tag}':100*cover.mean(),f'gap_{tag}':gap,f'RMSE_{tag}':np.sqrt(np.mean(e[finite]**2)),f'width_{tag}':np.mean(2*norm.ppf(.975)*np.sqrt(vv[valid])),f'finite_estimates_{tag}':int(finite.sum()),f'valid_intervals_{tag}':int(valid.sum()),f'fit_failures_{tag}':int((~a['success'][:,m]).sum()),f'failed_or_undefined_count_{tag}':int((~valid).sum()),f'oracle_variance_{tag}':V,f'standardized_gap_{tag}':gap/np.sqrt(V) if V>0 else np.nan,f'mean_Vhat_over_oracle_V_{tag}':vv[valid].mean()/V if V>0 else np.nan})
    d=inds['T']-inds['C'];r.update(difference_pp=100*d.mean(),paired_mcse_pp=100*d.std(ddof=1)/np.sqrt(R),abs_gap_change=abs(r['gap_T'])-abs(r['gap_C']),target_only_covered=int((d==1).sum()),common_only_covered=int((d==-1).sum()));rows.append(r)
df=pd.DataFrame(rows)
for pred,name in [('WorldCover','worldcover'),('zero','f0')]:df[df.prediction==pred].to_csv(args.output/f'lucas_crosscondition_{name}.csv',index=False,float_format='%.17g')
print('36 paired comparisons saved; fixed frame, all attempts.')
