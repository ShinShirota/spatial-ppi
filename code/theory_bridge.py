"""Finite deterministic identities bridging normalized theorem notation and code."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys,json,io,contextlib,platform
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import scipy,rasterio,pyproj
from regression_core import qr_system,from_system,extract
out=Path(sys.argv[1]);out.mkdir(exist_ok=False)
spec=json.loads((Path(__file__).resolve().parents[1]/'specs/theory_bridge.json').read_text());N=41;n=17;t=2*np.pi*np.arange(N)/N
B=np.column_stack([np.ones(N),np.sqrt(2)*np.cos(t),np.sqrt(2)*np.sin(t),np.sqrt(2)*np.cos(2*t),np.sqrt(2)*np.sin(2*t)])
X=np.column_stack([np.ones(N),np.cos(t)+.4*np.cos(2*t),np.sin(t)+.3*np.sin(2*t)])
f=.2*np.cos(3*t);Y=.5+np.cos(t)+.7*np.sin(2*t)+f
idx=np.flatnonzero(np.arange(N)%3!=0);rho=np.array(spec['rho']);Gamma=n*rho
R=np.array([[1,.2,-.1],[0,2,.3],[0,0,.7]])
checks=[];details=[];arrays={}
def check(name,a,b,tol=1e-10):
 err=float(np.max(np.abs(a-b)));checks.append({'name':name,'max_absolute_error':err,'atol':tol,'passed':bool(np.isfinite(err) and err<=tol)})
def fit(X,Y,f,corrected,cw=None):
 L=np.linalg.solve(X.T@X/N,X.T).T;C=B.T@L/N
 sys0=qr_system(B[idx],(Y-f)[idx,None],Gamma,cw)
 obj=from_system(sys0,C,corrected);est,cov,jk,psi,deletion=extract(obj,'target',1,3)
 return f@L/N+est[0],cov[0],psi[:,0],sys0,obj,C
for corrected in (False,True):
 est,cov,psi,s,o,C=fit(X,Y,f,corrected)
 Mh=(B[idx].T@B[idx])/n+np.diag(rho);alpha=np.linalg.solve(Mh,B[idx].T@(Y-f)[idx]/n)
 tt=np.linalg.solve(Mh,C);q=B[idx]@tt
 qp=B[idx]@(tt+np.linalg.solve(Mh,rho[:,None]*tt)) if corrected else q
 rr=B[idx]@np.linalg.solve(Mh,rho*alpha) if corrected else np.zeros(len(idx))
 phi=qp*((Y-f)[idx]-B[idx]@alpha)[:,None]-q*rr[:,None]
 name='corrected' if corrected else 'uncorrected'
 check(name+':phi=n psi',phi,n*psi)
 check(name+':covariance normalization',cov,phi.T@phi/n**2)
 h=1e-5;derivatives=[]
 for j in range(len(idx)):
  cw=np.ones(len(idx));cw[j]+=h;up=fit(X,Y,f,corrected,cw)[0];cw[j]-=2*h;down=fit(X,Y,f,corrected,cw)[0]
  derivatives.append((up-down)/(2*h))
 derivatives=np.array(derivatives);check(name+':inclusion derivative',psi,derivatives,2e-7)
 shifted=fit(X,Y,f+B[:,:3]@np.array([.3,.2,-.1]),corrected)
 check(name+':hard prediction shift estimate',est,shifted[0]);check(name+':hard prediction shift derivative',psi,shifted[2])
 recoded=fit(X@R,Y,f,corrected);Ri=np.linalg.inv(R)
 check(name+':X recoding estimate',recoded[0],Ri@est);check(name+':X recoding covariance',recoded[1],Ri@cov@Ri.T)
 Ychanged=Y.copy();outside=np.ones(N,bool);outside[idx]=False;Ychanged[outside]+=100
 changed=fit(X,Ychanged,f,corrected)
 check(name+':unobserved labels estimate',est,changed[0]);check(name+':unobserved labels covariance',cov,changed[1])
 details.append({'method':name,'omitted_second_term_max_abs_on_psi_scale':float(np.max(abs(q*rr[:,None]/n))),'finite_difference_error_if_second_term_omitted':float(np.max(abs(qp*((Y-f)[idx]-B[idx]@alpha)[:,None]/n-derivatives)))})
 arrays.update({name+'_psi':psi,name+'_phi':phi,name+'_finite_difference':derivatives,name+'_estimate':est,name+'_covariance':cov})
np.savez_compressed(out/'synthetic_arrays.npz',X=X,B=B,Y=Y,f=f,idx=idx,rho=rho,**arrays)
pirows=[]
for path in sorted((Path(__file__).resolve().parents[1]/'reference/lucas').glob('*_n*.npz')):
 with np.load(path) as z:
  nn=int(path.stem.rsplit('_n',1)[1]);pi=z['pi'];mx=float(pi.max())
  pirows.append({'cell':path.stem,'N':len(pi),'n':nn,'n_over_N':nn/len(pi),'pi_min':float(pi.min()),'pi_max':mx,'expected_labels':float(pi.sum()),'max_pi_over_one_minus_max_pi':mx/(1-mx)})
(out/'lucas_sampling_fractions.json').write_text(json.dumps(pirows,indent=2))
result={'passed':all(c['passed'] for c in checks),'checks':checks,'diagnostics':details,'N':N,'n':n,'actual_selected':len(idx),'scope':'finite synthetic checks only; not proof of asymptotic sampling or field assumptions'}
(out/'checks.json').write_text(json.dumps(result,indent=2))
buf=io.StringIO()
with contextlib.redirect_stdout(buf):np.show_config()
(out/'environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'BLAS':buf.getvalue(),'GDAL':rasterio.__gdal_version__,'rasterio_PROJ':rasterio.__proj_version__,'pyproj_PROJ':pyproj.proj_version_str,'geospatial_used':False},indent=2))
print(json.dumps(result,indent=2));sys.exit(0 if result['passed'] else 1)
