#!/usr/bin/env python3
"""Independent SVD/normal-equation evaluation, numerical derivatives, exact checks."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,sys,hashlib,time
import numpy as np
import pandas as pd
import sympy as sp
from regression_core import qr_system,from_system,extract
from experiment import real_inputs
ROOT=Path(__file__).resolve().parents[1]/'inputs'
RESULTS=Path(sys.argv[1]);OUT=Path(sys.argv[2]);OUT.mkdir(parents=True,exist_ok=False)

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False))

def independent(B,z,C,pen,bc,cw=None):
    m,K=B.shape;cw=np.ones(m) if cw is None else cw
    Aug=np.vstack((np.sqrt(cw)[:,None]*B,np.diag(np.sqrt(pen))))
    v=np.vstack((np.sqrt(cw)[:,None]*z,np.zeros((K,z.shape[1]))))
    # SVD coefficients (different factorization from primary QR).
    U,s,Vt=np.linalg.svd(Aug,full_matrices=False)
    if s[-1]<1e-12*s[0]:raise ValueError('SVD rank threshold')
    be=(Vt.T/s)@(U.T@v)
    Minv=(Vt.T/s**2)@Vt
    t=Minv@C;w=B@t;e=z-B@be
    if bc:
        wr=B@(t+Minv@(pen[:,None]*t));rr=B@Minv@(pen[:,None]*be)
        val=C.T@be+t.T@(pen[:,None]*be)
        psi=wr[:,:,None]*e[:,None,:]-w[:,:,None]*rr[:,None,:]
    else:val=C.T@be;psi=w[:,:,None]*e[:,None,:]
    return val,psi

D=real_inputs();cfg=json.loads((ROOT/'protocol/protocol.json').read_text());N=D['N'];L=D['L'];Z=D['Z'];nr=2;p=3
checks=[];fd=[];deletion=[];fail=[]
for sc in cfg['lucas_scenarios']:
 for n in cfg['lucas_budgets']:
  saved=np.load(RESULTS/f'{sc}_n{n}.npz')
  for r in [0,17,199]:
   idx=np.flatnonzero(saved['labels'][r]);m=len(idx)
   for j,method in enumerate(cfg['lucas_methods']):
    if method['mode']=='oracle':continue
    B=D['bases'][method['basis']];K=B.shape[1];pen=n*method.get('rho',0)*np.r_[np.zeros(method.get('hard',K)),np.ones(K-method.get('hard',K))];bc=method.get('bc',False)
    if method['mode']=='target':C=B.T@L/N;zz=Z[idx]
    else:C=B.mean(0)[:,None];zz=(Z[idx,:,None]*L[idx,None,:]).reshape(m,nr*p)
    val,ph=independent(B[idx],zz,C,pen,bc)
    if method['mode']=='target':value=val.T;phi=ph.transpose(0,2,1)
    else:value=val[0].reshape(nr,p);phi=ph[:,0].reshape(m,nr,p)
    vc=np.einsum('mri,mrj->rij',phi,phi)
    ee=float(np.max(abs(value+D['offsets']-saved['est'][r,j])));vv=float(np.max(abs(vc-saved['covariance'][r,j])))
    checks.append(dict(scenario=sc,n=n,r=r,method=method['id'],estimate_abs=ee,covariance_abs=vv))
    if ee>1e-8 or vv>1e-8:fail.append(checks[-1])
    if sc=='map_east' and n==100 and r==17:
     fit=from_system(qr_system(B[idx],zz,pen),C,bc)
     for ii in [0,m//2,m-1]:
      h=1e-5;cw=np.ones(m);cw[ii]+=h;v1,_=independent(B[idx],zz,C,pen,bc,cw);cw[ii]-=2*h;v0,_=independent(B[idx],zz,C,pen,bc,cw)
      err=float(np.max(abs((v1-v0)/(2*h)-fit['psi'][ii])))
      fd.append(dict(method=method['id'],row=ii,max_abs=err))
      keep=np.arange(m)!=ii;vd,_=independent(B[idx][keep],zz[keep],C,pen,bc)
      err2=float(np.max(abs(fit['est']-vd-fit['deletion'][ii])))
      deletion.append(dict(method=method['id'],row=ii,max_abs=err2))
      if err>2e-7 or err2>1e-8:fail.append({'type':'derivative/deletion','method':method['id'],'fd':err,'deletion':err2})
# Independent checks of the exact finite-sample remainder for arbitrary C.
rg=np.random.default_rng(202609141599);rex=[]
for it in range(60):
 N0=29;K=5;p0=3;b=rg.normal(size=(N0,K));b[:,0]=1;d=rg.normal(size=(N0,1));ell=rg.normal(size=(N0,p0));C=b.T@ell/N0
 lam=rg.uniform(.6,1.4,N0);lam/=lam.mean();n0=15.;pi=n0*lam/N0;flag=rg.random(N0)<pi
 P=np.diag([0,0,.7,1.3,2.]);A=b.T@(lam[:,None]*b)/N0;M=A+P;E=b[flag].T@b[flag]/n0-A;MH=M+E
 aa=np.linalg.solve(M,b.T@(lam[:,None]*d)/N0);t=np.linalg.solve(M,C);tp=t+np.linalg.solve(M,P@t)
 e=d-b@aa;rvec=b@np.linalg.solve(M,P@aa)
 ge=b.T@((flag-pi)[:,None]*e)/n0;gq=E@t;gqp=E@tp;gr=b.T@((flag-pi)[:,None]*rvec)/n0
 rr=-gqp.T@np.linalg.solve(MH,ge)+gq.T@np.linalg.solve(MH,gr)-gq.T@np.linalg.solve(MH,P@np.linalg.solve(MH,ge))
 ph=(b@tp)[:,:,None]*e[:,None,:]-(b@t)[:,:,None]*rvec[:,None,:]
 lin=np.einsum('n,npr->pr',flag-pi,ph)/n0
 ab=np.linalg.solve(MH,b[flag].T@d[flag]/n0);th=np.linalg.solve(MH,C)
 sample=C.T@ab+th.T@P@ab;population=C.T@aa+t.T@P@aa
 err=float(np.max(abs(sample-population-lin-rr)));rex.append(err)
 if err>1e-10:fail.append({'type':'remainder','case':it,'error':err})
# Exact eight-point mechanism. No numerical rounding or stochastic simulation.
q=8;xx=sp.Matrix([sp.sqrt(2)*sp.sin(2*sp.pi*j/q) for j in range(q)]);gg=sp.Matrix([sp.sqrt(2)*sp.cos(2*sp.pi*j/q) for j in range(q)])
one=sp.ones(q,1);B=one.row_join(xx).row_join(gg);X=one.row_join(xx);eta=sp.Rational(3,5);la=one+eta*xx.multiply_elementwise(gg);W=sp.diag(*la)
A=sp.simplify(B.T*W*B/q);J=sp.simplify(X.T*X/q);L=X*J.inv();c=B.T*one/q;C=B.T*L/q
qm=sp.simplify(B*A.inv()*c);qt=sp.simplify(B*A.inv()*C)
common_gap=sp.simplify(L.T*(W-sp.eye(q))*gg/q);target_gap=sp.simplify((W*qt-L).T*gg/q)
assert common_gap==sp.Matrix([0,sp.Rational(3,10)])
assert target_gap==sp.zeros(2,1)
exact={'N':q,'eta':'3/5','common_mean_gap':str(common_gap[0]),'common_slope_gap':str(common_gap[1]),'target_gaps':[str(v) for v in target_gap],'lambda_min':str(min(la)),'lambda_max':str(max(la)),'A':str(A),'C':str(C)}
# Full stored table reaggregation using no summarize.py imports.
countchecks=[];df=pd.read_csv(RESULTS/'summary.csv');zcrit=1.959963984540054
for sc in cfg['lucas_scenarios']:
 for n in cfg['lucas_budgets']:
  d=np.load(RESULTS/f'{sc}_n{n}.npz');truth=np.load(RESULTS/'population.npz')['truth']
  for jj,method in enumerate(cfg['lucas_methods']):
   for vart,key in [('linearization','covariance'),('jackknife','jackknife')]:
    if vart=='jackknife' and method['mode']=='oracle':continue
    for pp,pr in enumerate(['zero','WorldCover']):
     for cc,co in enumerate(['intercept','east_per100km','north_per100km']):
      ee=d['est'][:,jj,pp,cc]-truth[cc];vv=d[key][:,jj,pp,cc,cc]
      good=np.isfinite(ee)&np.isfinite(vv)&(vv>0);cover=good&(abs(ee)<=zcrit*np.sqrt(np.maximum(0,vv)))
      row=df[(df.scenario==sc)&(df.n==n)&(df.method==method['id'])&(df.variance==vart)&(df.prediction==pr)&(df.coefficient==co)].iloc[0]
      if int(row.covered)!=int(cover.sum()) or int(row.finite_intervals)!=int(good.sum()):fail.append({'type':'table_count','method':method['id'],'n':n,'scenario':sc})
      numerr=float(max(abs(row.rmse-np.sqrt(np.mean(ee**2))),abs(row.bias-np.mean(ee))))
      countchecks.append(numerr)
      if numerr>1e-12:fail.append({'type':'table_float','error':numerr})
# Census reduction and prior intercept output: necessary index/covariance sanity check.
intchecks=[]
prior_ids=['unweighted','hard_map2','hard6','hard_map7','hard3','BCsoft_map16_h7']
new_ids=['unweighted_score','map_score','spatial_score','map_spatial_score','target_X','soft_target_bc']
# target_X intercept is NOT the previous mean hard3 when J is nonorthogonal? X centered makes Cmean, so it is.
oldmethods=json.loads((ROOT/'source_snapshot/spatial_ppi_stage11/protocol_run.json').read_text())['methods'];oi=[m['id'] for m in oldmethods];nj=[m['id'] for m in cfg['lucas_methods']]
for sc in cfg['lucas_scenarios']:
 for n in cfg['lucas_budgets']:
  old=np.load(D['old']/f'{sc}_n{n}.npz');new=np.load(RESULTS/f'{sc}_n{n}.npz')
  assert np.array_equal(old['labels'],new['labels'])
  for a,b in zip(prior_ids,new_ids):
   er=float(np.max(abs(old['est'][:,oi.index(a)]-new['est'][:,nj.index(b),:,0])))
   ev=float(np.max(abs(old['variance'][:,oi.index(a)]-new['covariance'][:,nj.index(b),:,0,0])))
   intchecks.append({'scenario':sc,'n':n,'old':a,'new':b,'estimate_abs':er,'variance_abs':ev})
   if er>1e-8 or ev>1e-8:fail.append(intchecks[-1])
# target_X is exactly sample OLS on Y or residual followed by known map slope.
olscheck=[]
for sc in cfg['lucas_scenarios']:
 for n in cfg['lucas_budgets']:
  d=np.load(RESULTS/f'{sc}_n{n}.npz')
  for rr in [0,17,199]:
   flag=d['labels'][rr];be=np.linalg.lstsq(D['X'][flag],D['Z'][flag],rcond=None)[0].T+D['offsets']
   er=float(np.max(abs(be-d['est'][rr,nj.index('target_X')])));olscheck.append(er)
   if er>1e-8:fail.append({'type':'OLS','error':er})
result={'independent_svd_cases':len(checks),'estimate_max_abs':max(x['estimate_abs'] for x in checks),'covariance_max_abs':max(x['covariance_abs'] for x in checks),'finite_difference_cases':len(fd),'finite_difference_max_abs':max(x['max_abs'] for x in fd),'delete_one_cases':len(deletion),'delete_one_max_abs':max(x['max_abs'] for x in deletion),'exact_remainder_cases':len(rex),'exact_remainder_max_abs':max(rex),'rational_example':exact,'table_rows_checked':len(countchecks),'table_float_max_abs':max(countchecks),'previous_intercept_cases':len(intchecks),'previous_intercept_max_abs':max(x['estimate_abs'] for x in intchecks),'OLS_identity_cases':len(olscheck),'OLS_identity_max_abs':max(olscheck),'failures':fail,'run_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
dump(OUT/'audit.json',result);dump(OUT/'svd_details.json',checks);dump(OUT/'finite_difference.json',fd);dump(OUT/'delete_one.json',deletion);dump(OUT/'intercept_checks.json',intchecks)
print(json.dumps(result,indent=2))

sys.exit(1 if fail else 0)
