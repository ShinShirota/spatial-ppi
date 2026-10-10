#!/usr/bin/env python3
"""Deterministic Fourier checks of the new slope-only boundary construction.
No label draws, fitted confidence intervals, or empirical coverage are computed.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import itertools,json,numpy as np
import sys
OUT=Path(sys.argv[1]);OUT.mkdir(parents=True,exist_ok=False)
spec={'dimensions':[1,2],'k_values':[1,2,3],'grid_side':65,'dictionary_radius':'2*k','omitted_frequency':'4*k+3','intensity_amplitude':'.1*k^(-s)','s':'d/2+.75','spectral_a':'d+.75','sigma':.4,'tolerance':1e-11,'no_performance_simulation':True}
(OUT/'fixed_specification.json').write_text(json.dumps(spec,indent=2))
records=[]
for d in spec['dimensions']:
 for k in spec['k_values']:
  m=spec['grid_side'];pts=np.array(list(itertools.product(np.arange(m)/m,repeat=d)))
  N=len(pts);x=np.sqrt(2)*np.cos(2*np.pi*pts[:,0]);X=np.c_[np.ones(N),x]
  L=4*k+3;z=np.sqrt(2)*np.cos(2*np.pi*L*pts[:,0]);hc=np.sqrt(2)*np.cos(2*np.pi*(L+1)*pts[:,0]);hs=np.sqrt(2)*np.sin(2*np.pi*(L+1)*pts[:,0])
  assert m>4*(L+1)
  s=d/2+.75;a=d+.75;delta=.1*k**(-s);tau=(L+1)**(-a);sig=.4
  lam=1+delta*z
  cols=[np.ones(N)];rads=[0]
  for h in itertools.product(range(-2*k,2*k+1),repeat=d):
   nz=next((v for v in h if v),0)
   if nz<=0:continue
   arg=2*np.pi*(pts@np.asarray(h))
   cols.extend([np.sqrt(2)*np.cos(arg),np.sqrt(2)*np.sin(arg)]);rads.extend([max(abs(v) for v in h)]*2)
  B=np.array(cols).T;K=B.shape[1]
  penalty=np.array([0. if r<=k else (r/k)**(2*(d+1)) for r in rads]);A=B.T@(lam[:,None]*B)/N;M=A+np.diag(penalty)
  ell=X@np.linalg.inv(X.T@X/N);C=B.T@ell/N;t=np.linalg.solve(M,C);tp=t+np.linalg.solve(M,penalty[:,None]*t);q=B@tp
  D=lam[:,None]*q-ell
  a1=D.T@hc/N;a2=D.T@hs/N
  Omega=tau*(np.outer(a1,a1)+np.outer(a2,a2))+sig**2*(D.T@D)/N**2
  expected=np.diag([delta**2*sig**2/N,delta**2*(tau/2+sig**2/N)])
  records.append({'d':d,'k':k,'N':N,'K':K,'L':L,'s':s,'a':a,'delta':delta,'tau':tau,'lambda_min':float(lam.min()),'lambda_max':float(lam.max()),'A_identity_error':float(np.max(abs(A-np.eye(K)))),'q_loading_error':float(np.max(abs(q-ell))),'cos_structural_gap_loadings':a1.tolist(),'sine_structural_gap_loadings':a2.tolist(),'covariance':Omega.tolist(),'expected_covariance':expected.tolist(),'covariance_absolute_error':float(np.max(abs(Omega-expected))),'slope_loading_error':float(abs(a1[1]-delta/np.sqrt(2))),'intercept_structural_error':float(abs(a1[0])+abs(a2[0]))})
  np.savez_compressed(OUT/f'input_d{d}_k{k}.npz',X=X,B=B,lam=lam,hc=hc,hs=hs,penalty=penalty)
res={'specification':spec,'cases':records,'all_checks_below_tolerance':all(max(r['A_identity_error'],r['q_loading_error'],r['covariance_absolute_error'],r['slope_loading_error'],r['intercept_structural_error'])<spec['tolerance'] for r in records),'max_covariance_error':max(r['covariance_absolute_error'] for r in records),'max_q_loading_error':max(r['q_loading_error'] for r in records)}
(OUT/'checks.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
assert res['all_checks_below_tolerance']
