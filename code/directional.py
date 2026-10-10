#!/usr/bin/env python3
"""Reproduce the independent directional-covariance PPI Monte Carlo study.

This is a standalone additional simulation, NOT the manuscript's original
reproduction repository. Outputs include all attempted evaluations; no
simulated failures are excluded from the denominator.
Python 3.13.5 / NumPy 2.3.5 / SciPy 1.17.0 / pandas 2.2.3 used in 2026-10-09 run.
"""
import argparse
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.special import ndtr
from scipy.stats import norm
from threadpoolctl import threadpool_limits

threadpool_limits(limits=1)

def simulate(fields: int, labels: int, seed: int, null_control: bool):
    m, k, L = 129, 2, 11
    N = m*m
    nvals = (250,500,1000)
    Ms = (1,4,16)
    delta = 0.4*k**(-1.25)
    sigma, tau = 0.4, 0.04
    zz = norm.ppf(0.975)
    a,b = np.meshgrid(np.arange(m), np.arange(m), indexing='ij')
    s1, s2 = a.ravel()/m,b.ravel()/m
    def trig(s,rmax):
        arr,keys=[np.ones(len(s))],[('0',0)]
        for q in range(1,rmax+1):
            arr.extend([np.sqrt(2)*np.cos(2*np.pi*q*s),
                        np.sqrt(2)*np.sin(2*np.pi*q*s)])
            keys.extend([('c',q),('s',q)])
        return np.column_stack(arr),keys
    tx,kx=trig(s1,4)
    ty,ky=trig(s2,4)
    B=np.einsum('ni,nj->nij',tx,ty).reshape(N,-1)
    keys=[(a,b) for a in kx for b in ky]
    degree=np.array([max(a[1],b[1]) for a,b in keys])
    protect=np.where(degree<=k)[0]
    Bh=B[:,protect]
    lam_pen=np.where(degree>k,(degree/k)**4,0.0)
    x=np.sqrt(2)*np.cos(2*np.pi*s1)
    z=np.sqrt(2)*np.cos(2*np.pi*L*s1)
    lam=1+delta*z
    ell=np.column_stack([np.ones(N),x])
    CB=B.T@ell/N
    CBh=CB[protect]
    d=np.column_stack([delta*z,delta*z*x])
    if np.abs(B.T@B/N-np.eye(B.shape[1])).max()>1e-11:
        raise RuntimeError("Fourier Gram check failed")
    if np.abs(B.T@(lam[:,None]*B)/N-np.eye(B.shape[1])).max()>1e-11:
        raise RuntimeError("Acquisition-weighted Gram check failed")
    if np.abs(B.T@d/N).max()>1e-11 or lam.min()<=0:
        raise RuntimeError("Calibration-defect / positivity check failed")
    def trigh(freqs):
        ang=2*np.pi*s1[:,None]*np.array(freqs)[None,:]
        return np.sqrt(2)*np.cos(ang),np.sqrt(2)*np.sin(ang)
    sources=(trigh([L+1+2*j for j in range(16)]),
             trigh([L+3+2*j for j in range(16)]))
    def outcomes(fid):
        rng=np.random.default_rng(np.random.SeedSequence(seed,spawn_key=(fid,0)))
        co,si=rng.normal(size=16),rng.normal(size=16)
        eps=rng.normal(size=N)
        Y=np.empty((N,6))
        for mi,M in enumerate(Ms):
            for ci,(bc,bs) in enumerate(sources):
                Y[:,2*mi+ci]=(bc[:,:M]@co[:M]+bs[:,:M]@si[:M])*np.sqrt(tau/M)+sigma*eps
        return Y
    def estimates(bs,y,C,pen,true):
        ch=cho_factor(bs.T@bs+np.diag(pen),check_finite=False)
        alpha=cho_solve(ch,bs.T@y,check_finite=False)
        U0=cho_solve(ch,C,check_finite=False)
        kappa=cho_solve(ch,pen[:,None]*alpha,check_finite=False)
        Uplus=U0+cho_solve(ch,pen[:,None]*U0,check_finite=False)
        beta=C.T@(alpha+kappa)
        err=y-bs@alpha
        fitted_r=bs@kappa
        wplus=bs@Uplus
        w0=bs@U0
        psi=wplus[:,None,:]*err[:,:,None]-w0[:,None,:]*fitted_r[:,:,None]
        v=np.sum(psi**2,axis=0).T
        if not np.isfinite(v).all() or (v<=0).any():
            raise ValueError("non-positive or non-finite variance")
        return ((np.abs(beta-true)<=zz*np.sqrt(v)).T,
                (2*zz*np.sqrt(v)).T, ((beta-true)**2).T)
    cover=np.zeros((fields,labels,3,6,2,2),bool)
    width=np.full(cover.shape,np.nan)
    sqerr=np.full(cover.shape,np.nan)
    predictions=np.full((fields,3,6,2,2),np.nan)
    gaps=np.zeros((fields,6,2))
    energy=np.zeros((fields,6))
    failures=[]
    for fid in range(fields):
        Y=outcomes(fid)
        truth=ell.T@Y/N
        gaps[fid]=(d.T@Y/N).T
        energy[fid]=np.mean(Y*Y,axis=0)
        # Oracle conditional normal approximation, with exact fixed population.
        structural=[]
        for bb,pp in ((Bh,np.zeros(len(protect))),(B,lam_pen)):
            pa=(bb.T@(lam[:,None]*Y)/N)/(1+pp[:,None])
            kap=pp[:,None]*pa/(1+pp[:,None])
            res=Y-bb@(pa+kap)
            structural.append(res[:,:,None]*ell[:,None,:])
        for ni,n in enumerate(nvals):
            pi=n/N*lam
            for ti,phi in enumerate(structural):
                v=np.einsum('i,isa->sa',pi*(1-pi),phi*phi)/n**2
                standardized=gaps[fid]/np.sqrt(v)
                predictions[fid,ni,:,ti,:]=ndtr(zz-standardized)-ndtr(-zz-standardized)
        rng=np.random.default_rng(np.random.SeedSequence(seed,spawn_key=(fid,1)))
        for li in range(labels):
            u=rng.random(N)
            for ni,n in enumerate(nvals):
                idx=np.flatnonzero(u<n/N*lam)
                for ti,(bb,cc,pp) in enumerate(((Bh,CBh,np.zeros(len(protect))),
                                                  (B,CB,n*lam_pen))):
                    try:
                        c,w,se=estimates(bb[idx],Y[idx],cc,pp,truth)
                        cover[fid,li,ni,:,ti,:]=c
                        width[fid,li,ni,:,ti,:]=w
                        sqerr[fid,li,ni,:,ti,:]=se
                    except Exception as exc:
                        failures.append((fid,li,n,ti,str(exc)))
                        # False remains in the all-attempt denominator.
    rows,field_rows=[],[]
    for ni,n in enumerate(nvals):
        for mi,M in enumerate(Ms):
            for ti,name in enumerate(("Target hard","Target soft + BC")):
                ia,io=2*mi,2*mi+1
                fa=cover[:,:,ni,ia,ti,1].mean(axis=1)
                fo=cover[:,:,ni,io,ti,1].mean(axis=1)
                diff=fo-fa
                rows.append(dict(n=n,M=M,method=name,aligned_pct=100*fa.mean(),
                    orthogonal_pct=100*fo.mean(),paired_diff_pp=100*diff.mean(),
                    paired_mcse_pp=100*diff.std(ddof=1)/np.sqrt(fields),
                    aligned_mcse_pp=100*fa.std(ddof=1)/np.sqrt(fields),
                    orthogonal_mcse_pp=100*fo.std(ddof=1)/np.sqrt(fields),
                    predicted_al_pct=100*predictions[:,ni,ia,ti,1].mean(),
                    predicted_off_pct=100*predictions[:,ni,io,ti,1].mean(),
                    aligned_width=width[:,:,ni,ia,ti,1].mean(),
                    orth_width=width[:,:,ni,io,ti,1].mean(),
                    aligned_rmse=np.sqrt(sqerr[:,:,ni,ia,ti,1].mean()),
                    orth_rmse=np.sqrt(sqerr[:,:,ni,io,ti,1].mean())))
                for fid in range(fields):
                    for ci,condition in enumerate(("aligned","orthogonal")):
                        s=2*mi+ci
                        field_rows.append(dict(field_id=fid,n=n,M=M,method=name,
                            covariance=condition,
                            slope_coverage_field=cover[fid,:,ni,s,ti,1].mean(),
                            intercept_coverage_field=cover[fid,:,ni,s,ti,0].mean(),
                            slope_oracle_coverage_pred=predictions[fid,ni,s,ti,1],
                            slope_gap=gaps[fid,s,1],intercept_gap=gaps[fid,s,0],
                            residual_energy=energy[fid,s],
                            slope_ci_width_field=width[fid,:,ni,s,ti,1].mean(),
                            slope_rmse_field=np.sqrt(sqerr[fid,:,ni,s,ti,1].mean()),
                            label_repetitions=labels))
    null_rows=[]
    if null_control:
        null_cov=np.zeros((fields,labels,2,6,2),dtype=bool)
        for fid in range(fields):
            Y=outcomes(fid)
            truth=ell.T@Y/N
            rng=np.random.default_rng(np.random.SeedSequence(seed,spawn_key=(fid,2)))
            for li in range(labels):
                u=rng.random(N)
                for ni,n in enumerate((500,1000)):
                    idx=np.flatnonzero(u<n/N)
                    try:
                        ins,_,_=estimates(B[idx],Y[idx],CB,n*lam_pen,truth)
                        null_cov[fid,li,ni]=ins
                    except Exception as exc:
                        failures.append((fid,li,n,"uniform-null",str(exc)))
        for ni,n in enumerate((500,1000)):
            for mi,M in enumerate(Ms):
                fa=null_cov[:,:,ni,2*mi,1].mean(axis=1)
                fo=null_cov[:,:,ni,2*mi+1,1].mean(axis=1)
                df=fo-fa
                null_rows.append(dict(n=n,M=M,aligned_pct=100*fa.mean(),
                    orthogonal_pct=100*fo.mean(),paired_diff_pp=100*df.mean(),
                    paired_mcse_pp=100*df.std(ddof=1)/np.sqrt(fields)))
    return pd.DataFrame(rows),pd.DataFrame(field_rows),pd.DataFrame(null_rows),failures

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields",type=int,default=128)
    ap.add_argument("--labels",type=int,default=128)
    ap.add_argument("--seed",type=int,default=2026100901)
    ap.add_argument("--with-null",action="store_true")
    ap.add_argument("--out-dir",type=Path,default=Path("."))
    args=ap.parse_args()
    if args.fields<2 or args.labels<1:
        ap.error("fields >= 2, labels >= 1 are required")
    if args.out_dir.exists():
        ap.error("output directory must not already exist")
    t0=perf_counter()
    summary,field_data,null,failed=simulate(args.fields,args.labels,args.seed,args.with_null)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    summary.to_csv(args.out_dir/"directional_summary.csv",index=False,float_format="%.10g")
    field_data.to_csv(args.out_dir/"directional_field_level.csv",index=False,float_format="%.10g")
    if args.with_null:
        null.to_csv(args.out_dir/"directional_uniform_null.csv",index=False,float_format="%.10g")
    print(summary.to_string(index=False))
    if args.with_null: print("\nUniform acquisition control:\n",null.to_string(index=False))
    print("failed evaluations:",len(failed),failed[:5])
    print("seconds elapsed",round(perf_counter()-t0,1))
if __name__=="__main__":
    main()
