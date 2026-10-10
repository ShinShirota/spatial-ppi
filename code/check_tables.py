from pathlib import Path
import json,re,pandas as pd
import argparse
ap=argparse.ArgumentParser();ap.add_argument('results_root',type=Path);ap.add_argument('manuscript_directory',type=Path);args=ap.parse_args()
R=args.results_root;M=args.manuscript_directory;out=[]
def rows(s):
 return [re.sub(r'\s+','',x) for x in s.splitlines() if '&' in x and re.search(r'&\s*\$?[+-]?\d',x) and x.rstrip().endswith(r'\\')]
def verify(name,label,generated,doc):
 text=(M/doc).read_text();start=text.index('\\label{tab:'+label+'}')
 end=text.find('\\end{table}',start)
 if label=='gap-budget':end=text.find('\\end{longtable}',start)
 actual=rows(text[start:end]);expected=rows(generated)
 out.append(dict(table=name,passed=expected==actual,rows=len(expected),actual_rows=len(actual),differences=[{'generated':a,'manuscript':b} for a,b in zip(expected,actual) if a!=b]))
for name,label,file,doc in [('Main1','mechanism','mechanism','main.tex'),('Main2','east','lucas_east','main.tex'),('Main3','joint','joint','main.tex'),('S1','north','lucas_north','supplement.tex'),('S2','paired','paired','supplement.tex'),('S3','smalln','smalln','supplement.tex')]:
 verify(name,label,(R/f'tables_final/{file}_table.tex').read_text(),doc)
s=[];names={'target_hard_radius_k':'Hard','target_soft_radius_2k':'Soft','target_soft_bc_radius_2k':'Soft + BC','HT_oracle':'HT (oracle)'}
for x in json.loads((R/'full_gap/summary.json').read_text()):
 vals=[str(x['k']),str(x['n']),names[x['method']]]
 vals += [f"{100*x['coverage'][j]:.2f} ({100*x['coverage_mcse_field'][j]:.2f})" for j in range(2)]
 vals += [f"{100*x['predicted_coverage_field_mean'][1]:.2f}",f"{100*x['joint_coverage']:.2f} ({100*x['joint_mcse_field']:.2f})"]
 s.append(' & '.join(vals)+r'\\')
for name,label,lines in [('S4','gap-budget',s)]:
 (R/f'tables_final/{name}.tex').write_text('\n'.join(lines)+'\n');verify(name,label,'\n'.join(lines),'supplement.tex')
for name,label,file in [('S5','directional','summary'),('S6','directional-null','uniform_null')]:
 df=pd.read_csv(R/f'full_directional/directional_{file}.csv');keys=(['method'] if name=='S5' else [])+['n','M'];df=df.sort_values(keys);s=[]
 for _,x in df.iterrows():
  vals=([x['method']] if name=='S5' else [])+[str(int(x['n'])),str(int(x['M']))]+[f'{x[c]:.2f}' for c in ['aligned_pct','orthogonal_pct','paired_diff_pp','paired_mcse_pp']]
  s.append(' & '.join(vals)+r'\\')
 (R/f'tables_final/{name}.tex').write_text('\n'.join(s)+'\n');verify(name,label,'\n'.join(s),'supplement.tex')
df=pd.read_csv(R/'crosscondition/lucas_crosscondition_worldcover.csv');s=[]
for (scenario,n),g in df.groupby(['scenario','n'],sort=False):
 vals=[scenario.replace('_','/'),str(n)]
 for coeff in ['east_per100km','north_per100km']:
  x=g[g.coefficient==coeff].iloc[0];vals += [f"{x.coverage_C:.2f} / {x.coverage_T:.2f}",f"${x.difference_pp:+.2f}$ ({x.paired_mcse_pp:.2f})"]
 s.append(' & '.join(vals)+r'\\')
(R/'tables/S7.tex').write_text('\n'.join(s)+'\n');verify('S7','lucas-crosscondition','\n'.join(s),'supplement.tex')
(R/'table_display_comparison.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

raise SystemExit(0 if all(c['passed'] for c in out) else 1)
