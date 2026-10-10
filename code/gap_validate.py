import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from regression_core import population_fit
import argparse
parser=argparse.ArgumentParser();parser.add_argument('run_directory');args=parser.parse_args()
r=Path(args.run_directory);sp=json.loads((r/'spec_v1.json').read_text());N=sp['N'];s=np.arange(N)/N;k=2;L=4*k+3
lam=1+sp['c']*k**-sp['s']*np.sqrt(2)*np.cos(2*np.pi*L*s)
B=np.column_stack([np.ones(N)]+[f(2*np.pi*h*s)*np.sqrt(2) for h in range(1,2*k+1) for f in [np.cos,np.sin]])
X=B[:,:2];ell=X@np.linalg.inv(X.T@X/N);C=B.T@ell/N;y=np.load(r/'primary/field_k2_f0.npz')['y'];truth=ell.T@y/N
rng=np.random.default_rng(np.random.SeedSequence([sp['master_seed'],2,0,1]));pi=128/N*lam;ids=rng.random(N)<pi
freq=np.array([0]+[h for h in range(1,2*k+1) for _ in range(2)]);pen=np.where(freq<=k,0.,(freq/k)**4)
assert not (r/'independent_formula_checks.json').exists(), 'Validation output already exists'
checks=[]
rows=[json.loads(l) for l in (r/'primary/replicates.jsonl').open()]
for j in range(3):
 b=B[:,:2*k+1] if j==0 else B;rho=np.zeros(2*k+1) if j==0 else pen;c=C[:2*k+1] if j==0 else C
 po=population_fit(b,y[:,None],c,lam,rho,j==2);ref=next(x for x in rows if x['k']==2 and x['n']==128 and x['field']==0 and x['label']==0 and x['method']==sp['methods'][j]);ph=po['phi'][:,:,0]
 v=np.einsum('i,ij,ik->jk',pi*(1-pi)/128**2,ph,ph)
 T=b[ids].T@b[ids]+np.diag(128*rho);al=np.linalg.solve(T,b[ids].T@y[ids]);est=c.T@al
 if j==2:est+=c.T@np.linalg.solve(T,128*rho*al)
 checks.append({'method':ref['method'],'estimate_abs_error':float(np.max(abs(est-ref['estimate']))),'gap_abs_error':float(np.max(abs(po['mu'][:,0]-truth-ref['gap']))),'V_abs_error':float(np.max(abs(v-ref['V']))),'linear_abs_error':float(np.max(abs(ph[ids].sum(0)/128-ph.T@pi/128-ref['linear'])))})
assert all(v<sp['analytic_check_absolute_tolerance'] for c in checks for key,v in c.items() if key!='method')
(r/'independent_formula_checks.json').write_text(json.dumps(checks,indent=2))

print(json.dumps(checks,indent=2))
