"""Rebuild only the scientific LUCAS/WorldCover frame from official source files.
Does not generate labels, fit estimators, scan all pixels, or overwrite inputs.
"""
from pathlib import Path
import argparse,csv,json,re,hashlib,platform
import numpy as np,pandas as pd,rasterio,pyproj
ap=argparse.ArgumentParser();ap.add_argument('--csv',type=Path,required=True);ap.add_argument('--tiles',type=Path,nargs=3,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=False)
rows=[];excluded=[]
with args.csv.open(encoding='utf-8-sig',newline='') as f:
 reader=csv.reader(f);header=next(reader)
 for row in reader:
  if not row or all(not x.strip() for x in row):excluded.append({'line':reader.line_num,'reason':'blank'});continue
  if len(row)==1 and re.fullmatch(r'[\d.,]+ rows selected\.',row[0].strip()):excluded.append({'line':reader.line_num,'reason':'export footer'});continue
  if len(row)!=len(header):raise ValueError('Unexpected record length')
  rows.append(row)
df=pd.DataFrame(rows,columns=header);assert len(df)==2665 and not df.POINT_ID.duplicated().any()
lon=df.TH_LONG.to_numpy(float);lat=df.TH_LAT.to_numpy(float);pred=np.zeros(len(df),int);seen=np.zeros(len(df),int)
for path in args.tiles:
 with rasterio.open(path) as ds:
  assert ds.crs.to_epsg()==4326 and ds.nodata==0 and ds.tags()['product_version']=='V2.0.0'
  ix=np.flatnonzero((lon>=ds.bounds.left)&(lon<ds.bounds.right)&(lat<=ds.bounds.top)&(lat>ds.bounds.bottom))
  for i,v in zip(ix,ds.sample(zip(lon[ix],lat[ix]))):
   if v[0]==0:raise ValueError('NoData at reference site; no deletion/imputation permitted')
   seen[i]+=1;pred[i]=int(v[0]==10)
assert (seen==1).all()
y=df.LC1.str.startswith('C').to_numpy(int);east,north=pyproj.Transformer.from_crs(4326,3035,always_xy=True).transform(lon,lat);xy=np.c_[east,north]
u=(xy-xy.min(0))/(xy.max(0)-xy.min(0));s=2*u-1;powers=[(a,d-a) for d in range(5) for a in range(d,-1,-1)]
mon=np.column_stack([s[:,0]**a*s[:,1]**b for a,b in powers])
def normal(a):
 q,r=np.linalg.qr(a,mode='reduced');return q*np.where(np.diag(r)<0,-1.,1.)*np.sqrt(len(df))
np.savez_compressed(args.output/'population_inputs.npz',Y=y,f=pred,B=normal(mon),Bmap=normal(np.column_stack((mon[:,0],pred,mon[:,1:]))),POINT_ID=df.POINT_ID.to_numpy(int))
pd.DataFrame({'POINT_ID':df.POINT_ID.to_numpy(int),'easting_3035':east,'northing_3035':north}).to_csv(args.output/'lucas_worldcover_joined.csv',index=False,float_format='%.17g')
meta={'N':len(df),'rows_dropped':0,'excluded_nonrecords':excluded,'CRS':'EPSG:3035 from EPSG:4326 TH_LONG/TH_LAT always_xy','Y':'LC1 begins C','f':'WorldCover2021 v200 native Map code10','source_hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [args.csv,*args.tiles]},'python':platform.python_version(),'numpy':np.__version__,'rasterio':rasterio.__version__,'GDAL':rasterio.__gdal_version__,'PROJ':pyproj.proj_version_str}
(args.output/'frame_provenance.json').write_text(json.dumps(meta,indent=2));print('Rebuilt',len(df),'sites; no sampling/estimation.')
