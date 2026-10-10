#!/usr/bin/env python3
from __future__ import annotations
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse,json,time,hashlib,platform,sys
from pathlib import Path
import numpy as np
import pandas as pd
from regression_core import qr_system,from_system,extract,normalize_columns,population_fit
ROOT=Path(__file__).resolve().parents[1]/'inputs'

def dump(p,obj):Path(p).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False))
def now():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())

def real_inputs():
    old=ROOT/'source_snapshot/spatial_ppi_stage11/results/primary'
    data=np.load(old/'population_inputs.npz');N=len(data['Y'])
    # Recover meters without looking at Y; all-frame projected coordinates were audited in stage 11.
    df=pd.read_csv(ROOT/'inputs/lucas_worldcover_joined.csv')
    if not np.array_equal(data['POINT_ID'],df.POINT_ID.to_numpy()):raise ValueError('Point order mismatch')
    coords=df[['easting_3035','northing_3035']].to_numpy(float)
    X=np.column_stack((np.ones(N),(coords-coords.mean(0))/100000))
    J=X.T@X/N;L=np.linalg.solve(J,X.T).T
    BM=data['Bmap'];BS=data['B'];x,y=X[:,1],X[:,2];f=data['f']
    product=np.column_stack([np.ones(N),f,x,y,x*x,x*y,y*y,x*f,y*f,x**3,x*x*y,x*y*y,y**3])
    bases={'one':BS[:,:1],'map2':BM[:,:2],'spatial6':BS[:,:6],'map7':BM[:,:7],
           'product13':normalize_columns(product),'X3':normalize_columns(X),'map16':BM}
    Z=np.column_stack((data['Y'],data['Y']-f))
    offsets=np.vstack((np.zeros(3),f@L/N))
    truth=data['Y']@L/N
    return dict(N=N,Y=data['Y'],f=f,X=X,J=J,L=L,bases=bases,Z=Z,offsets=offsets,
                truth=truth,POINT_ID=data['POINT_ID'],coords=coords,old=old)

def simulate(dataset,out,replicates=None):
    protocol=json.loads((ROOT/'protocol/protocol.json').read_text())
    if dataset=='lucas':
        D=real_inputs();methods=protocol['lucas_methods'];budgets=protocol['lucas_budgets'];scens=protocol['lucas_scenarios']
        R=2000 if replicates is None else replicates
        if not 1 <= R <= 2000: raise ValueError('replicates must be 1..2000')
        prednames=['zero','WorldCover'];coefnames=['intercept','east_per100km','north_per100km']
        labels_source=f'saved stage11, first {R} replicates per cell; unchanged order'
    elif dataset=='mechanism':
        cfg=protocol['mechanism'];G=cfg['grid_side'];N=G*G
        coord=np.arange(G)/G;ss=np.meshgrid(coord,coord,indexing='ij');x=np.sqrt(2)*np.sin(2*np.pi*ss[0].ravel());g=np.sqrt(2)*np.cos(2*np.pi*ss[0].ravel())
        X=np.column_stack((np.ones(N),x));J=X.T@X/N;L=np.linalg.solve(J,X.T).T
        bases={'one':np.ones((N,1)),'signal3':normalize_columns(np.column_stack((np.ones(N),x,g))),
               'product5':normalize_columns(np.column_stack((np.ones(N),x,g,x*x,x*g)))}
        rg=np.random.Generator(np.random.PCG64DXSM(cfg['population_seed']))
        coef=rg.normal(size=(cfg['populations'],2))*cfg['structure_sd'];eps=rg.normal(size=(N,cfg['populations']))*cfg['nugget_sd']
        Z=x[:,None]*coef[:,0]+g[:,None]*coef[:,1]+eps
        f=.7*x;Y=f[:,None]+Z
        truth=Y.T@L/N;offsets=np.tile(f@L/N,(cfg['populations'],1))
        D=dict(N=N,Y=Y,f=f,X=X,J=J,L=L,bases=bases,Z=Z,offsets=offsets,truth=truth,x=x,g=g,field_coefficients=coef,nugget=eps)
        methods=cfg['methods'];budgets=cfg['budgets'];scens=['uniform','cross_loading'];R=cfg['replicates'];prednames=[f'population_{j:02d}' for j in range(cfg['populations'])];coefnames=['intercept','sine_slope'];labels_source='new fixed PCG64DXSM streams'
    else:raise ValueError(dataset)
    out.mkdir(parents=True,exist_ok=False);N=D['N'];Z=D['Z'];L=D['L'];p=L.shape[1];nr=Z.shape[1];nm=len(methods)
    np.savez_compressed(out/'population.npz',**{k:v for k,v in D.items() if isinstance(v,np.ndarray)},**{'basis_'+k:v for k,v in D['bases'].items()})
    dump(out/'metadata.json',{'dataset':dataset,'N':N,'R':R,'responses':prednames,'coefficients':coefnames,'methods':methods,'label_source':labels_source,'started':now(),'numpy':np.__version__,'python':platform.python_version()})
    start=time.perf_counter();exceptions=[];poprows=[];detailed_checks=[]
    for si,sc in enumerate(scens):
      for ni,n in enumerate(budgets):
        if dataset=='lucas':
            with np.load(D['old']/f'{sc}_n{n}.npz') as z:
                labels=z['labels'][:R];pi=z['pi'];sizes=z['sizes'][:R]
            lam=pi*N/n
        else:
            eta=0 if sc=='uniform' else protocol['mechanism']['eta']
            lam=1+eta*D['x']*D['g'];pi=n*lam/N
            seq=np.random.SeedSequence([protocol['mechanism']['label_seed'],si,ni])
            rg=np.random.Generator(np.random.PCG64DXSM(seq));labels=rg.random((R,N))<pi;sizes=labels.sum(1)
        if pi.min()<=0 or pi.max()>=1:raise ValueError('Invalid design; no capping')
        est=np.full((R,nm,nr,p),np.nan);cov=np.full((R,nm,nr,p,p),np.nan);jk=np.full_like(cov,np.nan)
        linear=np.full_like(est,np.nan);rem=np.full_like(est,np.nan);stats=np.full((R,nm,5),np.nan);ok=np.zeros((R,nm),bool)
        config=[];popinfo=[]
        for j,sp in enumerate(methods):
            if sp['mode']=='oracle':
                config.append(None)
                gap=np.zeros((nr,p));phi=Z[:,:,None]*L[:,None,:]/lam[:,None,None]
                q=L/lam[:,None];normmatrix=None
            else:
                B=D['bases'][sp['basis']];K=B.shape[1];rho=sp.get('rho',0)*np.r_[np.zeros(sp.get('hard',K)),np.ones(K-sp.get('hard',K))]
                if sp['mode']=='target': c=B.T@L/N;zz=Z
                else:c=B.mean(0)[:,None];zz=(Z[:,:,None]*L[:,None,:]).reshape(N,nr*p)
                config.append(dict(B=B,C=c,z=zz,penalty=n*rho,bc=sp.get('bc',False)))
                po=population_fit(B,zz,c,lam,rho,sp.get('bc',False))
                if sp['mode']=='target':
                    mu=po['mu'].T;phi=po['phi'].transpose(0,2,1);q=po['q'];
                    defect=lam[:,None]*q-L
                    gdirect=Z.T@defect/N
                else:
                    mu=po['mu'][0].reshape(nr,p);phi=po['phi'][:,0,:].reshape(N,nr,p);q=po['q']
                    defect=(lam*q[:,0]-1)[:,None]*L
                    gdirect=Z.T@defect/N
                gap=D['offsets']+mu-D['truth']
                diff=float(np.max(np.abs(gdirect-gap)))
                detailed_checks.append({'scenario':sc,'n':n,'method':sp['id'],'gap_formula_difference':diff})
                if diff>1e-9:raise ArithmeticError('Population gap identity failed')
            V=np.einsum('n,nri,nrj->rij',pi*(1-pi),phi,phi)/n**2
            popinfo.append((gap,phi,V,(pi[:,None,None]*phi).sum(0)/n))
            for ri in range(nr):
             for ci in range(p):
                poprows.append(dict(scenario=sc,n=n,method=sp['id'],response=prednames[ri],coefficient=coefnames[ci],truth=float(D['truth'][ci] if D['truth'].ndim==1 else D['truth'][ri,ci]),gap=float(gap[ri,ci]),linear_design_variance=float(V[ri,ci,ci]),pi_min=float(pi.min()),pi_max=float(pi.max())))
        # Share each QR decomposition across common/target and corrected/uncorrected versions.
        for r,flag in enumerate(labels):
            idx=np.flatnonzero(flag);m=len(idx);cache={}
            for j,(sp,cf) in enumerate(zip(methods,config)):
              try:
                if sp['mode']=='oracle':
                    ss=Z[idx,:,None]*L[idx,None,:]/(N*pi[idx,None,None])
                    val=ss.sum(0);vv=np.einsum('n,nri,nrj->rij',1-pi[idx],ss,ss);vjk=np.full_like(vv,np.nan)
                    stat=[0,1.,0,0,float(np.max(np.sum(abs(ss),axis=0)))];fit=None
                else:
                    key=(sp['basis'],sp['mode'],sp.get('rho',0),sp.get('hard',D['bases'][sp['basis']].shape[1]))
                    if key not in cache:
                        cache[key]=qr_system(cf['B'][idx],cf['z'][idx],cf['penalty'])
                    fit=from_system(cache[key],cf['C'],cf['bc'])
                    val,vv,vjk,_,_=extract(fit,sp['mode'],nr,p)
                    # Target-specific slope weights are necessarily signed; not an instability count.
                    neg=float(np.sum(fit['weights']<0)) if sp['mode']=='common' else np.nan
                    stat=[fit['h'].max(),fit['condition_R'],neg,float(np.max(abs(fit['weights']))),float(np.max(np.sum(abs(fit['weights']),axis=0)))]
                est[r,j]=D['offsets']+val;cov[r,j]=vv;jk[r,j]=vjk;stats[r,j]=stat;ok[r,j]=True
                gap,phi,V,mean_phi=popinfo[j];ln=phi[idx].sum(0)/n-mean_phi
                linear[r,j]=ln;rem[r,j]=est[r,j]-D['truth']-gap-ln
              except Exception as ex:
                exceptions.append(dict(scenario=sc,n=n,r=r,method=sp['id'],m=m,error=type(ex).__name__,message=str(ex)))
            if (r+1)%500==0: print(dataset,sc,n,r+1,'seconds',round(time.perf_counter()-start,1),flush=True)
        np.savez_compressed(out/f'{sc}_n{n}.npz',est=est,covariance=cov,jackknife=jk,linear=linear,remainder=rem,success=ok,stats=stats,labels=labels,sizes=sizes,pi=pi,gap=np.stack([a[0] for a in popinfo]),oracle_covariance=np.stack([a[2] for a in popinfo]))
    pd.DataFrame(poprows).to_csv(out/'population_diagnostics.csv',index=False,float_format='%.17g')
    dump(out/'exceptions.json',exceptions);dump(out/'population_checks.json',detailed_checks)
    dump(out/'completed.json',{'dataset':dataset,'seconds':time.perf_counter()-start,'fit_attempts':len(scens)*len(budgets)*R*nm,'failures':len(exceptions),'completed':now()})
    print('COMPLETE',dataset,'seconds',time.perf_counter()-start,'failures',len(exceptions),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',choices=['lucas','mechanism'],required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replicates',type=int);args=ap.parse_args();simulate(args.dataset,args.output,args.replicates)
