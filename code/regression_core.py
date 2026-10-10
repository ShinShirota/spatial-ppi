#!/usr/bin/env python3
"""Fixed-target calibration for finite-population regression.

The scientific design matrix X and the auxiliary dictionary B are distinct.
All responses supplied here must be observed at the selected sites only.
Primary solves use augmented QR; no pseudoinverse, clipping or adaptive ridge.
"""
from __future__ import annotations
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import numpy as np
from scipy.linalg import solve_triangular


def qr_system(b: np.ndarray, z: np.ndarray, penalty: np.ndarray,
              caseweights: np.ndarray | None = None) -> dict:
    """Fit several response columns to one deterministic diagonal-penalty system."""
    m,k=b.shape
    if z.ndim!=2 or z.shape[0]!=m or penalty.shape!=(k,):
        raise ValueError('Incompatible matrix shapes')
    if m<2 or np.any(penalty<0): raise ValueError('Invalid sample or penalty')
    cw=np.ones(m) if caseweights is None else np.asarray(caseweights)
    if np.any(cw<0): raise ValueError('Negative caseweight')
    aug=np.vstack((np.sqrt(cw)[:,None]*b,np.diag(np.sqrt(penalty))))
    Q,R=np.linalg.qr(aug,mode='reduced')
    dg=np.abs(np.diag(R))
    if len(dg)!=k or dg.min()<=1e-12*dg.max():
        raise np.linalg.LinAlgError('Augmented system fails fixed rank tolerance')
    def solve(rhs):
        return solve_triangular(R,solve_triangular(R.T,rhs,lower=True,check_finite=False),check_finite=False)
    rhs=np.vstack((np.sqrt(cw)[:,None]*z,np.zeros((k,z.shape[1]))))
    be=solve_triangular(R,Q.T@rhs,check_finite=False)
    U=solve(b.T)
    e=z-b@be
    return dict(b=b,z=z,penalty=penalty,beta=be,e=e,solve=solve,U=U,
                h=np.einsum('ij,ji->i',b,U),
                condition_R=float(np.linalg.cond(R)),caseweights=cw)


def from_system(sys: dict, C: np.ndarray, corrected: bool) -> dict:
    """C is K x number_of_targets; returns target x response estimates.

    Derivative and delete-one arrays are site x target x response. The
    derivatives are with respect to each inclusion/caseweight, not Y_i.
    The delete-one formula holds at unit caseweights and fixed penalty.
    """
    b=sys['b']; P=sys['penalty']; solve=sys['solve']; be=sys['beta']; e=sys['e']
    C=np.atleast_2d(C)
    if C.shape[0]!=b.shape[1]:raise ValueError('C must be K x targets')
    t=solve(C); w0=b@t
    if corrected:
        tplus=t+solve(P[:,None]*t)
        w=b@tplus
        r=b@solve(P[:,None]*be)
        est=C.T@be+t.T@(P[:,None]*be)
        psi=w[:,:,None]*e[:,None,:]-w0[:,:,None]*r[:,None,:]
    else:
        w=w0;r=np.zeros_like(e);est=C.T@be;psi=w[:,:,None]*e[:,None,:]
    # delete-one is an exact fixed-penalty identity; no tuning in this run.
    den=1-sys['h'];deletion=np.full_like(psi,np.nan)
    if np.allclose(sys['caseweights'],1,rtol=0,atol=0) and den.min()>1e-12:
        deletion=psi/den[:,None,None]
        if corrected:
            kp=np.sum(P[:,None]*sys['U']**2,axis=0)
            deletion+=w0[:,:,None]*e[:,None,:]*(kp/den**2)[:,None,None]
    if not (np.isfinite(est).all() and np.isfinite(psi).all()):
        raise FloatingPointError('Nonfinite estimate/derivative')
    return dict(est=est,psi=psi,deletion=deletion,weights=w,w0=w0,
                condition_R=sys['condition_R'],h=sys['h'])


def extract(fit:dict,mode:str,nresponses:int,p:int) -> tuple:
    """Canonical output: response x coefficient; site x response x coefficient."""
    if mode=='target':
        est=fit['est'].T;psi=fit['psi'].transpose(0,2,1);deletion=fit['deletion'].transpose(0,2,1)
    elif mode=='common':
        est=fit['est'][0].reshape(nresponses,p)
        psi=fit['psi'][:,0,:].reshape(-1,nresponses,p)
        deletion=fit['deletion'][:,0,:].reshape(-1,nresponses,p)
    else:raise ValueError(mode)
    cov=np.einsum('mri,mrj->rij',psi,psi)
    if np.isfinite(deletion).all():
        m=len(psi);center=deletion-deletion.mean(0)
        jk=(m-1)/m*np.einsum('mri,mrj->rij',center,center)
    else:jk=np.full_like(cov,np.nan)
    return est,cov,jk,psi,deletion


def normalize_columns(raw:np.ndarray) -> np.ndarray:
    Q,R=np.linalg.qr(raw,mode='reduced');dg=np.abs(np.diag(R))
    if dg.min()<1e-12*dg.max():raise ValueError('Population dictionary is dependent')
    sign=np.where(np.diag(R)<0,-1.,1.)
    return Q*sign*np.sqrt(len(raw))


def population_fit(B,z,C,lam,rho,corrected):
    """Oracle analysis only. Never called to construct feasible sample estimates."""
    N=len(B);A=B.T@(lam[:,None]*B)/N;M=A+np.diag(rho)
    be=np.linalg.solve(M,B.T@(lam[:,None]*z)/N)
    t=np.linalg.solve(M,C);q=B@t;e=z-B@be
    if corrected:
        tp=t+np.linalg.solve(M,rho[:,None]*t)
        qp=B@tp;r=B@np.linalg.solve(M,rho[:,None]*be)
        mu=C.T@be+t.T@(rho[:,None]*be)
        phi=qp[:,:,None]*e[:,None,:]-q[:,:,None]*r[:,None,:]
    else:qp=q;r=np.zeros_like(e);mu=C.T@be;phi=q[:,:,None]*e[:,None,:]
    return dict(mu=mu,phi=phi,q=qp,q0=q,e=e,r=r,A=A,M=M,beta=be,t=t)
