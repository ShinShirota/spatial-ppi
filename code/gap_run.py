import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,hashlib,platform,importlib.metadata as md,contextlib,io,time
import numpy as np
from scipy.stats import norm,chi2
sys.path.insert(0,str(Path(__file__).resolve().parent))
from regression_core import qr_system,from_system
import argparse,shutil
parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
r=Path(args.output);r.mkdir(parents=True,exist_ok=False)
shutil.copy2(Path(__file__).resolve().parents[1]/'gap/spec_v1.json',r/'spec_v1.json')
spec=json.loads((r/'spec_v1.json').read_text());pilot=args.pilot
out=r/('pilot' if pilot else 'primary');out.mkdir(exist_ok=False)
buf=io.StringIO()
with contextlib.redirect_stdout(buf):np.show_config()
env={'python':platform.python_version(),'packages':{d.metadata['Name']:d.version for d in md.distributions()},'blas':buf.getvalue()}
try:
 import rasterio,pyproj
 env.update(GDAL=rasterio.__gdal_version__,PROJ=pyproj.proj_version_str)
except ImportError as e:env['geo']=str(e)
(out/'environment.json').write_text(json.dumps(env,indent=2))
(out/'hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),r/'spec_v1.json',Path(__file__).resolve().parent/'regression_core.py']},indent=2))
N=spec['N'];pos=np.arange(N)/N;X=np.column_stack([np.ones(N),np.sqrt(2)*np.cos(2*np.pi*pos)]);loading=X@np.linalg.inv(X.T@X/N);zcrit=norm.ppf(.975);jointcrit=chi2.ppf(.95,2)
records=(out/'replicates.jsonl').open('w');cells=[];attempt=0
for k in spec['k_values'][:1] if pilot else spec['k_values']:
 L=4*k+3;delta=spec['c']*k**(-spec['s']);tau=k**(-spec['a']);lam=1+delta*np.sqrt(2)*np.cos(2*np.pi*L*pos)
 hc=np.sqrt(2)*np.cos(2*np.pi*(L+1)*pos);hs=np.sqrt(2)*np.sin(2*np.pi*(L+1)*pos)
 B=np.column_stack([np.ones(N)]+[f(2*np.pi*h*pos)*np.sqrt(2) for h in range(1,2*k+1) for f in [np.cos,np.sin]])
 freq=np.array([0]+[h for h in range(1,2*k+1) for _ in range(2)]);rho=np.where(freq<=k,0.,(freq/k)**4);C=B.T@loading/N
 A=B.T@(lam[:,None]*B)/N
 err=float(np.max(abs(A-np.eye(len(freq)))));assert err<spec['analytic_check_absolute_tolerance']
 basis=[(B[:,:2*k+1],np.zeros(2*k+1),C[:2*k+1]),(B,rho,C),(B,rho,C)]
 for field in range(1 if pilot else spec['fields']):
  frng=np.random.default_rng(np.random.SeedSequence([spec['master_seed'],k,field,0]));zc,zs=frng.normal(size=2)*np.sqrt(tau);eps=frng.normal(size=N);y=zc*hc+zs*hs+spec['sigma']*eps;truth=loading.T@y/N
  np.savez_compressed(out/f'field_k{k}_f{field}.npz',y=y,zc=zc,zs=zs)
  oracle=[]
  for j,(b,pen,c) in enumerate(basis):
   al=(b.T@(lam*y)/N)/(1+pen);e=y-b@al;rr=b@(pen*al/(1+pen));ph=loading*(e-(rr if j==2 else 0))[:,None]
   g=c.T@al-truth;direct=(loading*(lam-1)[:,None]).T@y/N
   assert np.max(abs(g-direct))<spec['analytic_check_absolute_tolerance']
   oracle.append((g,ph))
  label_rng=np.random.default_rng(np.random.SeedSequence([spec['master_seed'],k,field,1]))
  for label in range(1 if pilot else spec['labels_per_field']):
   u=label_rng.random(N)
   for n in spec['n_values'][:1] if pilot else spec['n_values']:
    pi=n/N*lam;assert pi.max()<1;ind=u<pi;ids=np.flatnonzero(ind);ys=y[ind,None]
    for method in range(4):
     attempt+=1;rec={'k':k,'n':n,'field':field,'label':label,'method':spec['methods'][method],'sample_size':len(ids),'truth':truth.tolist(),'max_pi':float(pi.max()),'status':'ok'}
     try:
      if method<3:
       b,pen,c=basis[method];sysm=qr_system(b[ind],ys,n*pen);fit=from_system(sysm,c,method==2);est=fit['est'][:,0];ps=fit['psi'][:,:,0];vh=ps.T@ps;cond=fit['condition_R'];g,ph=oracle[method]
       vv=np.einsum('i,ij,ik->jk',pi*(1-pi)/n**2,ph,ph);lin=ph[ind].sum(0)/n-ph.T@pi/n
      else:
       term=loading[ind]*ys;est=(term/pi[ind,None]).sum(0)/N;vh=np.einsum('i,ij,ik->jk',(1-pi[ind])/pi[ind]**2/N**2,term,term);cond=None;g=np.zeros(2)
       ph=loading*y[:,None]/lam[:,None];vv=np.einsum('i,ij,ik->jk',pi*(1-pi)/n**2,ph,ph);lin=est-truth
      rem=est-truth-g-lin;diag=np.diag(vh);valid=np.isfinite(diag)&(diag>0);width=np.where(valid,2*zcrit*np.sqrt(np.maximum(diag,0)),np.nan);coverage=valid&(abs(est-truth)<=width/2)
      eig=np.linalg.eigvalsh(vh);jointok=bool(np.isfinite(eig).all() and eig[0]>spec['covariance_eigen_ratio_tolerance']*eig[-1] and eig[-1]>0)
      joint=bool((est-truth)@np.linalg.solve(vh,est-truth)<=jointcrit) if jointok else False
      shift=g/np.sqrt(np.diag(vv));pred=norm.cdf(zcrit-shift)-norm.cdf(-zcrit-shift)
      rec.update(estimate=est.tolist(),gap=g.tolist(),linear=lin.tolist(),remainder=rem.tolist(),Vhat=vh.tolist(),V=vv.tolist(),condition_R=cond,marginal_valid=valid.tolist(),coverage=coverage.tolist(),width=width.tolist(),joint_valid=jointok,joint_coverage=joint,noncentrality=float(g@np.linalg.solve(vv,g)),predicted_coverage=pred.tolist())
     except (np.linalg.LinAlgError,ValueError,FloatingPointError) as exc:rec.update(status='failed',reason=str(exc),coverage=[False,False],joint_coverage=False,marginal_valid=[False,False],joint_valid=False)
     records.write(json.dumps(rec,allow_nan=True)+'\n')
  records.flush();print(f'completed k={k} field={field} attempts={attempt}',flush=True)
 cells.append({'k':k,'A_identity_error':err,'gap_covariance_formula':[[delta**2*spec['sigma']**2/N,0],[0,delta**2*(tau/2+spec['sigma']**2/N)]],'K_hard':2*k+1,'K_soft':4*k+1})
records.close();(out/'population_checks.json').write_text(json.dumps(cells,indent=2));(out/'completion.json').write_text(json.dumps({'attempts':attempt,'pilot':pilot},indent=2))
if not pilot:assert attempt==spec['expected_attempts']
