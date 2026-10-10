from pathlib import Path
import pandas as pd,json,numpy as np
import argparse
ap=argparse.ArgumentParser()
ap.add_argument('--lucas',type=Path,required=True)
ap.add_argument('--mechanism',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();m=args.output;m.mkdir(parents=True,exist_ok=False)
a=pd.read_csv(args.lucas/'summary.csv');s=pd.read_csv(args.mechanism/'summary.csv');j=pd.read_csv(args.lucas/'joint_coverage.csv');p=pd.read_csv(args.lucas/'paired_comparisons.csv')
names={'unweighted_score':'Unweighted score PPI','map_score':'Common map (2)','spatial_score':'Common spatial (6)','map_spatial_score':'Common map/spatial (7)','product_score':'Common product space (13)','target_X':r'Target ${\mathbf{X}}$ (3)','target_map_spatial':'Target map/spatial (7)','soft_score':'Common soft (16)','soft_score_bc':'Common soft + BC (16)','soft_target':'Target soft (16)','soft_target_bc':'Target soft + BC (16)','HT_oracle':'Known-$\\pi$ HT (oracle)'}
def table(name,label,caption,headers,body,fmt):
 text='\\begin{table}[htbp]\n\\centering\\small\n\\caption{'+caption+'}\\label{'+label+'}\n\\begin{tabular}{'+fmt+'}\\toprule\n'+' & '.join(headers)+'\\\\\n\\midrule\n'
 text+='\n'.join(' & '.join(row)+'\\\\' for row in body)+'\n\\bottomrule\\end{tabular}\n\\end{table}\n'
 (m/name).write_text(text)
truth=np.load(args.lucas/'population.npz')['truth']
table('truth_table.tex','tab:truth','Complete-frame regression target. Coordinate slopes are fraction per 100 km; not national areal or causal effects.',['Intercept','East slope','North slope'],[[f'{x:.8f}' for x in truth]],'rrr')
for coef,tag in [('east_per100km','east'),('north_per100km','north')]:
 x=a[(a.scenario=='map_east')&(a.n==250)&(a.prediction=='WorldCover')&(a.variance=='linearization')&(a.coefficient==coef)]
 body=[]
 for _,row in x.iterrows():body.append([names[row.method],f'{100*row.coverage_all:.2f}',f'{100*row.coverage_mcse:.2f}',f'{row.rmse:.5f}',f'{row.mean_width:.5f}',f'{row.gap_mean:+.5f}'])
 table(f'lucas_{tag}_table.tex',f'tab:{tag}',f'LUCAS {tag} slope, map/east acquisition, expected $n=250$. All 12 predeclared methods; full inclusion-derivative covariance. MCSE is in percentage points.',['Method ($K$)','Coverage (\\%)','MCSE','RMSE','Width','Pop. gap'],body,'lrrrrr')
x=s[(s.scenario=='cross_loading')&(s.n==400)&(s.variance=='linearization')];body=[]
for mid,name in [('unweighted_score','Unweighted'),('common_B3','Common (3)'),('product_B5','Product (5)'),('target_B3','Target (3)'),('HT_oracle','HT oracle')]:
 u=x[(x.method==mid)&(x.coefficient=='intercept')].iloc[0];v=x[(x.method==mid)&(x.coefficient=='sine_slope')].iloc[0]
 body.append([name,f'{100*u.coverage_all:.2f}',f'{100*u.coverage_mcse:.2f}',f'{100*v.coverage_all:.2f}',f'{100*v.coverage_mcse:.2f}',f'{v.gap_rms:.4f}'])
table('mechanism_table.tex','tab:mechanism','Constructed stationary periodic mechanism, $n=400$. Crossed-population/label Monte Carlo SEs; 32 populations, 500 shared label sets.',['Method ($K$)','Mean (\\%)','MCSE','Slope (\\%)','MCSE','Slope gap RMS'],body,'lrrrrr')
x=p[(p.scenario=='map_east')&(p.n==250)&(p.prediction=='WorldCover')&(p.variance=='linearization')&(p.coefficient=='east_per100km')];body=[]
for _,row in x.iterrows():body.append([names[row['plus']]+' minus '+names[row['minus']],f'{100*row.coverage_difference:+.2f}',f'{100*row.coverage_difference_mcse:.2f}'])
table('paired_table.tex','tab:paired','Paired east-slope coverage differences at $n=250$, in percentage points. Same labels for both methods.',['Comparison','Difference','MCSE'],body,'lrr')
x=j[(j.scenario=='map_east')&(j.n==250)&(j.prediction=='WorldCover')&(j.variance=='linearization')&(j.group=='slopes_only')]
table('joint_table.tex','tab:joint','Joint 95\\% ellipsoid for both coordinate slopes, map/east acquisition, $n=250$.',['Method','Coverage (\\%)','MCSE (pp)'],[[names[row.method],f'{100*row.coverage_all:.2f}',f'{100*row.coverage_mcse:.2f}'] for _,row in x.iterrows()],'lrr')
body=[]
for mid in ['map_score','map_spatial_score','target_map_spatial','soft_target_bc']:
 x=a[(a.scenario=='map_east')&(a.n==100)&(a.prediction=='WorldCover')&(a.coefficient=='east_per100km')&(a.method==mid)]
 u=x[x.variance=='linearization'].iloc[0];v=x[x.variance=='jackknife'].iloc[0]
 body.append([names[mid],f'{100*u.coverage_all:.2f}',f'{u.mean_width:.4f}',f'{100*v.coverage_all:.2f}',f'{v.mean_width:.4f}'])
table('smalln_table.tex','tab:smalln','East slope at $n=100$: small-sample undercoverage and the width cost of the jackknife.',['Method','Linear (\\%)','Width','Jackknife (\\%)','Width'],body,'lrrrr')
