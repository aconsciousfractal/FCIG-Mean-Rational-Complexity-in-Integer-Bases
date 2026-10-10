from pathlib import Path
from fractions import Fraction as F
import sys,math,json
from exact import log_hi,log_lo,need,out_decimal
S=10**9

def grid(root):
 lines=(Path(root)/'data/grid.tsv').read_text().splitlines();den=int(lines[0]);rows=[list(map(int,l.split())) for l in lines[1:]]
 need(len(rows)==2688 and all(len(r)==8 and r[0]==i for i,r in enumerate(rows)),'complete grid')
 for i,r in enumerate(rows):
  _,lo,hi,tl,th,a,b,d=r;t=i%21-10
  need(0<lo<hi<2*den and a>0 and d>0,'positive geometry')
  need((tl,th)==(max(-20,2*t-1),min(20,2*t+1)),'shear coverage')
  need((lo,hi)==tuple(rows[i//21*21][1:3]),'radial block')
  if i%21==0 and i:need(lo==rows[i-21][2],'radial continuity')
 need(3*rows[-1][2]**2>=4*den**2,'Gauss upper endpoint')
 return den,rows

def profile(p,path,root,k=1,ag=4):
 q=42*p if p>1 else 42;K=q**k;NL=96*(p+1) if p>1 else 96;den,grows=grid(root)
 raw=[list(map(int,l.split())) for l in Path(path).read_text().splitlines()]
 need(len(raw)==2688,'complete output')
 for i,r in enumerate(raw):need(len(r)==6 and r[:3]==grows[i][:3] and 0<=r[4]<NL and -ag<=r[5]<=ag,'profile labels and cells')
 v=[F(max(r[3] for r in raw[21*j:21*(j+1)]),S) for j in range(128)]
 for i in range(126,-1,-1):v[i]=max(v[i],v[i+1])
 a=F(29,14);rho=F(raw[0][1],den)
 bulk=v[-1]+sum((v[j]-v[j+1])*min(F(1),a*F(raw[21*j][2],den)) for j in range(127))
 c=log_hi(q-1)+log_hi(K)/2-log_lo(rho)/2+F(5,2)*rho
 cusp=a*rho*(max(c-v[0],F(0))+F(1,2));total=bulk+cusp
 return dict(p=p,q=q,k=k,angular_nodes=2*ag+1,transfer_index=K,labels=NL,cells=2688,witnesses=2688*NL*(K//q),bulk_upper=str(bulk),cusp_upper=str(cusp),total_upper=str(total),upper_decimal=out_decimal(total),target_lower=str(F(2,85)*log_lo(q)),direct_radical_gate=total<F(2,85)*log_lo(q),target_mu_4_over_25=total<F(4,25),normalized_upper=str(max(F(0),total)/log_lo(q)),monotone_profile_ticks=[int(x*S) for x in v])
if __name__=='__main__':
 root=Path(__file__).resolve().parent.parent;p=int(sys.argv[1]);r=profile(p,sys.argv[2],root);print(json.dumps(r,indent=2))
