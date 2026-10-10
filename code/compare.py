from pathlib import Path
import json,sys
import numpy as np
import pandas as pd
import argparse
ap=argparse.ArgumentParser();ap.add_argument('reference',type=Path);ap.add_argument('candidate',type=Path);ap.add_argument('report',type=Path);args=ap.parse_args()
old=args.reference;new=args.candidate;rows=[]
if args.report.exists():raise FileExistsError(args.report)
if not list(old.glob('*.npz')) or not list(old.glob('*.csv')):raise ValueError('Reference must contain NPZ and CSV results')
def cmp(a,b,name):
 if a.shape!=b.shape:rows.append({'name':name,'passed':False,'reason':'shape'});return
 if a.dtype.kind in 'fc':
  mask=np.array_equal(np.isnan(a),np.isnan(b)) and np.array_equal(np.isposinf(a),np.isposinf(b)) and np.array_equal(np.isneginf(a),np.isneginf(b))
  finite=np.isfinite(a)&np.isfinite(b);dif=np.abs(a[finite]-b[finite]);bad=dif>1e-8+1e-10*np.abs(a[finite])
  rows.append({'name':name,'passed':bool(mask and not bad.any()),'max_abs':float(dif.max()) if dif.size else 0,'masks_match':mask})
 else:rows.append({'name':name,'passed':bool(np.array_equal(a,b))})
for f in sorted(old.glob('*.npz')):
 with np.load(f) as a,np.load(new/f.name) as b:
  assert set(a.files)==set(b.files)
  for key in a.files:cmp(a[key],b[key],f.name+':'+key)
for f in sorted(old.glob('*.csv')):
 a=pd.read_csv(f);b=pd.read_csv(new/f.name);assert list(a.columns)==list(b.columns)
 for key in a:cmp(a[key].to_numpy(),b[key].to_numpy(),f.name+':'+key)
out={'passed':all(x['passed'] for x in rows),'checks':len(rows),'failures':[x for x in rows if not x['passed']],'details':rows}
args.report.write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='details'},indent=2));sys.exit(0 if out['passed'] else 1)
