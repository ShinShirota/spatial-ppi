from pathlib import Path
import json,hashlib
import numpy as np
import argparse
parser=argparse.ArgumentParser();parser.add_argument('run_directory');args=parser.parse_args()
r=Path(args.run_directory);s=json.loads((r/'spec_v1.json').read_text());out=r/'primary'
assert json.loads((out/'completion.json').read_text())['attempts']==s['expected_attempts']
rows=[json.loads(l) for l in (out/'replicates.jsonl').open()];assert len(rows)==s['expected_attempts']
assert len({(x['k'],x['n'],x['field'],x['label'],x['method']) for x in rows})==len(rows)
assert not any((r/x).exists() for x in ['summary.json','field_summary.json','summary.txt','output_hashes.json']), 'Summary outputs already exist'
summ=[];fieldrows=[]
for k in s['k_values']:
 for n in s['n_values']:
  for method in s['methods']:
   rr=[x for x in rows if (x['k'],x['n'],x['method'])==(k,n,method)];ok=[x for x in rr if x['status']=='ok'];fr=[]
   for f in range(s['fields']):
    ff=[x for x in rr if x['field']==f];assert len(ff)==s['labels_per_field'];fo=[x for x in ff if x['status']=='ok']
    v={'k':k,'n':n,'method':method,'field':f,'coverage':np.mean([x['coverage'] for x in ff],axis=0).tolist(),'joint_coverage':float(np.mean([x['joint_coverage'] for x in ff])),'failures':sum(x['status']!='ok' for x in ff)}
    if fo:v.update(predicted_coverage=fo[0]['predicted_coverage'],gap=fo[0]['gap'],noncentrality=fo[0]['noncentrality'])
    fr.append(v);fieldrows.append(v)
   cov=np.array([x['coverage'] for x in fr]);jc=np.array([x['joint_coverage'] for x in fr]);pred=np.array([x['predicted_coverage'] for x in fr if 'predicted_coverage' in x])
   stat={'k':k,'n':n,'method':method,'attempts':len(rr),'fit_failures':len(rr)-len(ok),'marginal_invalid':np.sum([np.logical_not(x['marginal_valid']) for x in rr],axis=0).tolist(),'joint_invalid':sum(not x['joint_valid'] for x in rr),'coverage':cov.mean(0).tolist(),'coverage_mcse_field':(cov.std(0,ddof=1)/np.sqrt(s['fields'])).tolist(),'joint_coverage':float(jc.mean()),'joint_mcse_field':float(jc.std(ddof=1)/np.sqrt(s['fields'])),'predicted_coverage_field_mean':pred.mean(0).tolist() if len(pred) else None,'n_over_N':n/s['N'],'k2_over_n':k*k/n,'Khard2_over_n':(2*k+1)**2/n,'n_k_minus_r':n*k**-(s['a']+2*s['s'])}
   if ok:
    rem=np.array([x['remainder'] for x in ok]);ratios=np.array([np.diag(x['Vhat'])/np.diag(x['V']) for x in ok]);vgap=np.array([x['gap'] for x in fr if 'gap' in x]);
    stat.update(remainder_scaled_rms=(np.sqrt(n*np.mean(rem**2,axis=0))).tolist(),variance_ratio_mean=ratios.mean(0).tolist(),width_mean=np.nanmean([x['width'] for x in ok],axis=0).tolist(),sqrt_n_gap_rms=np.sqrt(n*np.mean(vgap**2,axis=0)).tolist())
   summ.append(stat)
(r/'summary.json').write_text(json.dumps(summ,indent=2));(r/'field_summary.json').write_text(json.dumps(fieldrows,indent=2))
lines=['Finite-budget slope-gap experiment (all-attempt coverage; field-cluster MCSE)','k n method mean_cover slope_cover slope_MCSE predicted_slope joint_cover failures joint_invalid']
for v in summ:lines.append(f"{v['k']} {v['n']} {v['method']} {v['coverage'][0]:.4f} {v['coverage'][1]:.4f} {v['coverage_mcse_field'][1]:.4f} {v['predicted_coverage_field_mean'][1]:.4f} {v['joint_coverage']:.4f} {v['fit_failures']} {v['joint_invalid']}")
(r/'summary.txt').write_text('\n'.join(lines)+'\n');print('\n'.join(lines))
files=[p for p in out.iterdir() if p.is_file()];(r/'output_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
