"""Exact integration of explicitly validated refined radial/shear rectangles."""
from pathlib import Path
from fractions import Fraction as F
import json
from exact import need,log_hi,log_lo,out_decimal
S=10**9

def grid(root):
 root=Path(root);meta=json.loads((root/'data/grid_meta.json').read_text());nr,nt,td=meta['radial_intervals'],meta['shear_intervals'],meta['shear_denominator']
 lines=(root/'data/grid.tsv').read_text().splitlines();den=int(lines[0]);rows=[list(map(int,l.split())) for l in lines[1:]]
 need(len(rows)==meta['cells']==nr*nt,'complete refined grid')
 for i,r in enumerate(rows):
  need(len(r)==8 and r[0]==i,'all rectangle indices')
  _,lo,hi,tl,th,a,b,d=r
  need(0<lo<hi<2*den and 0<a<10**9 and 0<d<10**9 and abs(b)<10**9,'positive arithmetic geometry')
  need(-td//2<=tl<th<=td//2,'shear coordinate domain')
  need((lo,hi)==tuple(rows[i//nt*nt][1:3]),'radial block')
  need(lo*d<=a*den<=hi*d,'centre radial containment')
  need(tl*a<=b*td<=th*a,'centre shear containment')
  if i%nt:need(tl==rows[i-1][4],'no shear gaps or overlaps')
  else:
   need(tl==-td//2 and rows[i+nt-1][4]==td//2,'entire shear interval')
   if i:need(lo==rows[i-nt][2],'no radial gaps or overlaps')
 need(3*rows[-1][2]**2>=4*den**2,'Gauss upper endpoint')
 return den,rows

def profile(p,path,root,k=1,ag=16):
 root=Path(root);meta=json.loads((root/'data/grid_meta.json').read_text());nr,nt=meta['radial_intervals'],meta['shear_intervals'];den,gr=grid(root)
 q=42*p if p>1 else 42;nl=96*(p+1) if p>1 else 96;K=q**k
 raw=[list(map(int,l.split())) for l in Path(path).read_text().splitlines()]
 need(len(raw)==len(gr),'complete profile output')
 for i,r in enumerate(raw):need(len(r)==6 and r[:3]==gr[i][:3] and 0<=r[4]<nl and -ag<=r[5]<=ag,'profile labels and cells')
 v=[F(max(z[3] for z in raw[nt*j:nt*(j+1)]),S) for j in range(nr)]
 for j in range(nr-2,-1,-1):v[j]=max(v[j],v[j+1])
 a=F(29,14);rho=F(raw[0][1],den)
 bulk=v[-1]+sum((v[j]-v[j+1])*min(F(1),a*F(raw[nt*j][2],den)) for j in range(nr-1))
 c=log_hi(q-1)+log_hi(K)/2-log_lo(rho)/2+F(5,2)*rho
 cusp=a*rho*(max(c-v[0],F(0))+F(1,2));total=bulk+cusp
 return dict(p=p,q=q,k=k,angular_nodes=2*ag+1,transfer_index=K,labels=nl,cells=len(gr),radial_intervals=nr,shear_intervals=nt,shear_denominator=meta['shear_denominator'],witnesses=len(gr)*nl*(K//q),bulk_upper=str(bulk),cusp_upper=str(cusp),total_upper=str(total),upper_decimal=out_decimal(total),target_lower=str(F(2,85)*log_lo(q)),direct_radical_gate=total<F(2,85)*log_lo(q),normalized_upper=str(max(F(0),total)/log_lo(q)),monotone_profile_ticks=[int(x*S) for x in v])
