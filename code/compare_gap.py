from pathlib import Path
import json,math
import argparse
p=argparse.ArgumentParser();p.add_argument('reference',type=Path);p.add_argument('replay',type=Path);p.add_argument('report',type=Path);args=p.parse_args()
a=args.reference;b=args.replay;count=0;maximum=0.
if args.report.exists():raise FileExistsError(args.report)
def compare(x,y,path=''):
 global count,maximum
 if isinstance(x,dict):
  assert x.keys()==y.keys(),path
  for k in x:compare(x[k],y[k],path+'/'+k)
 elif isinstance(x,list):
  assert len(x)==len(y),path
  for i,(v,w) in enumerate(zip(x,y)):compare(v,w,path+'/'+str(i))
 elif isinstance(x,float):
  count+=1
  if not math.isfinite(x):assert (math.isnan(x) and math.isnan(y)) or x==y,path
  else:
   assert math.isfinite(y) and abs(x-y)<=1e-8+1e-10*abs(x),(path,x,y)
   maximum=max(maximum,abs(x-y))
 else:assert x==y,(path,x,y)
with (a/'primary/replicates.jsonl').open() as f,(b/'primary/replicates.jsonl').open() as g:
 import itertools
 n=0
 for x,y in itertools.zip_longest(f,g):
  assert x is not None and y is not None
  compare(json.loads(x),json.loads(y),str(n));n+=1
assert n==36864
for name in ['summary.json','field_summary.json','independent_formula_checks.json','primary/population_checks.json','primary/completion.json']:
 compare(json.loads((a/name).read_text()),json.loads((b/name).read_text()),name)
result={'attempts_compared':n,'float_values_compared':count,'max_absolute_difference':maximum,'tolerance':'1e-8 + 1e-10 * abs(reference)','discrete_values_and_nonfinite_masks':'exact','passed':True,'scope':'same-seed same-environment packaging replay, not independent validation'}
args.report.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
