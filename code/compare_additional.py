from pathlib import Path
import json,re,numpy as np,pandas as pd
import argparse
ap=argparse.ArgumentParser();ap.add_argument('results_root',type=Path);args=ap.parse_args()
P=Path(__file__).resolve().parents[1];R=args.results_root;checks=[]
def check(name,a,b,atol,rtol):
 a=np.asarray(a);b=np.asarray(b)
 ok=a.shape==b.shape
 numeric=a.dtype.kind in 'fiu'
 if ok:ok=bool(np.allclose(a,b,atol=atol,rtol=rtol,equal_nan=True)) if numeric else bool(np.array_equal(a,b))
 err=float(np.nanmax(np.abs(a-b))) if numeric and a.size and a.shape==b.shape else None
 checks.append(dict(name=name,passed=ok,max_absolute_difference=err))
for file,keys in [('summary',['n','M','method']),('field_level',['field_id','n','M','covariance','method']),('uniform_null',['n','M'])]:
 a=pd.read_csv(P/f'reference/directional/{file}.csv').rename(columns={'orth_mcse_pp':'orthogonal_mcse_pp'});b=pd.read_csv(R/f'full_directional/directional_{file}.csv')
 if 'method' in a:a['method']=a.method.replace({'hard':'Target hard','soft+BC':'Target soft + BC'})
 a=a.sort_values(keys);b=b.sort_values(keys)
 assert set(a)==set(b)
 for c in a:check('directional/'+file+'/'+c,a[c],b[c],5.1e-8,1e-10)
for f in (P/'reference/crosscondition').glob('*.csv'):
 a=pd.read_csv(f);b=pd.read_csv(R/'crosscondition'/f.name);assert list(a)==list(b)
 for c in a:check('crosscondition/'+f.name+'/'+c,a[c],b[c],1e-10,1e-8)
a=np.load(P/'inputs/source_snapshot/spatial_ppi_stage11/results/primary/population_inputs.npz');b=np.load(R/'frame_rebuilt/population_inputs.npz')
for k in a:check('frame/'+k,a[k],b[k],1e-8 if k in ['B','Bmap'] else 0,1e-10 if k in ['B','Bmap'] else 0)
a=pd.read_csv(P/'inputs/inputs/lucas_worldcover_joined.csv');b=pd.read_csv(R/'frame_rebuilt/lucas_worldcover_joined.csv')
for c in a:check('frame/'+c,a[c],b[c],0 if c=='POINT_ID' else 1e-5,0)
(R/'additional_comparison.json').write_text(json.dumps({'passed':all(c['passed'] for c in checks),'checks':checks},indent=2))
print('Additional checks:',len(checks),'failed:',[c for c in checks if not c['passed']])

raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
