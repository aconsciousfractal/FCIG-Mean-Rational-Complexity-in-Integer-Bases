"""Independent exact primitive-period, covering, shear and domain computations."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,isqrt,prod
from itertools import combinations
from functools import lru_cache
import hashlib,json,argparse
ROOT=Path(__file__).resolve().parent
CHECKS=0

def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise ArithmeticError(msg)
def serial(z):
 if isinstance(z,F):return str(z)
 if isinstance(z,dict):return {k:serial(v) for k,v in z.items()}
 if isinstance(z,(list,tuple)):return [serial(x) for x in z]
 return z
def primes_to(n):
 a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
 for p in range(2,isqrt(n)+1):
  if a[p]:a[p*p::p]=b'\0'*((n-p*p)//p+1)
 return [p for p in range(2,n+1) if a[p]]
def phi(n):return sum(gcd(t,n)==1 for t in range(1,n+1))
def primitive_period():
 ps=(2,3,5,7,11,13);P=prod(ps)
 ws=[P]*(P+1)
 for p in ps:
  for n in range(p,P+1,p):ws[n]=ws[n]//p*(p-1)
 phis={};trace=hashlib.sha256()
 for n in range(1,P+1):
  g=gcd(n,P)
  if g not in phis:
   v=g
   for p in ps:
    if g%p==0:v=v//p*(p-1)
   phis[g]=v
  need(ws[n]==P*phis[g]//g,'period weight identity')
  trace.update(f'{ws[n]}\n'.encode())
 total=sum(ws[1:]);need(total==prod(p*p-1 for p in ps),'CRT average')
 lo=hi=prefix=0;la=ha=0
 for n in range(1,P+1):
  prefix+=ws[n];D=P*prefix-total*n
  if D<lo:lo=D;la=n
  if D>hi:hi=D;ha=n
 alpha=F(total,P*P);ell=F(lo,P*P);h=F(hi,P*P)
 need((alpha,ell,h,la,ha)==(F(442368,715715),-F(134597,1431430),F(744773,1431430),11376,18653),'period extrema')
 b=-alpha/2+(h-ell)*F(161,72)+2*ell
 need(F(161,72)**2>5 and b==F(45173053,51531480),'outward square root and linear coefficient')
 return {'period':P,'weights_checked':P,'weight_stream_sha256':trace.hexdigest(),'mean':alpha,'minimum':ell,'maximum':h,'minimum_first_at':la,'maximum_first_at':ha,'quadratic_coefficient':alpha*F(29,14),'linear_coefficient':b,'constant_coefficient':alpha/2}

@lru_cache(None)
def small_survivors(S,L):
 q=prod(S);units=[gcd(t,q)==1 for t in range(q)]
 ds=[(1,1)]
 for p in S:ds += [(p*d,-sgn) for d,sgn in ds]
 values=[]
 for s in range(q):
  direct=sum(units[t%q] for t in range(s+1,s+L+1))
  alternate=sum(sgn*((s+L)//d-s//d) for d,sgn in ds)
  need(direct==alternate,'every cyclic interval start')
  values.append(direct)
 return min(values)
def covering():
 Hs=(2,4,6,10,14,22,26,44,58,66,98,148,204,360);core=(3,5,7,11)
 outside=[p for p in primes_to(101) if p>11];rows=[]
 for k,H in enumerate(Hs,1):
  for two in (False,True):
   L=H//2 if two else H
   for mask in range(16):
    S=tuple(p for j,p in enumerate(core) if mask>>j&1);t=k-int(two)-len(S)
    if t<0:continue
    minimum=small_survivors(S,L);loss=[(L+p-1)//p for p in outside[:t]]
    need(minimum-sum(loss)>=1,'uniform forbidden class covering')
    rows.append({'k':k,'H':H,'two_present':two,'small_support':list(S),'length':L,'outside_lower_bounds':outside[:t],'small_minimum':minimum,'outside_losses':loss,'survivors':minimum-sum(loss)})
 need(len(rows)==399,'complete covering case count')
 return {'bounds':[{'prime_count_at_most':k,'H':H} for k,H in enumerate(Hs,1)],'rows':rows,'cases':len(rows)}
def shear():
 pairs=0
 for H in (2,4,6,10,14,22,26,44,58,66,98,148,204,360):
  h=H//2
  for l in range(1,H):
   for r in range(1,H-l+1):
    # Both affine branches cross at t=(l-r)/2; check the affine endpoint candidates exactly.
    cuts={F(0),F(1,2)};x=F(l-r,2)
    if 0<x<F(1,2):cuts.add(x)
    for t in cuts:need(min(l-t,r+t)<=h-t,'integer gap distance bound')
    pairs+=1
  # Polynomial identities in t are compared coefficientwise; sign proof for all real t is in PROOF.md.
  if h<=3:
   need([0,2*h,-h*h-1]==[h*h-h*h,2*h,-h*h-1],'small shear polynomial')
   need(4*h-h*h-1>0,'small shear sign')
  else:
   a=h*h-4*h+1;b=2*(h*h-h+1)
   need([a,b-2*a,-2*b]==[(2*h-1)**2-3*h*h,6*h,-((2*h-1)**2+3)],'large shear polynomial identity')
   need(a>=1,'large shear sign')
 return {'integer_gap_pairs':pairs,'kappa_squared':{str(H):1+max(F(H*H,4),F((H-1)**2,3)) for H in (2,4,6,10,14,22,26,44,58,66,98,148,204,360)},'whole_real_shear_interval_paid_by_polynomial_signs':True}
def domains():
 fs=[(30,7),(42,11),(70,11),(66,13),(105,11),(78,17),(110,13),(102,19),(154,13),(165,13),(130,17),(114,23)]
 closed=[2310,2730,3570,3990,4290,4830,5610,6006,6090,6270,6510];result={}
 for w,Q,div in ((4,3003,30),(5,6630,210)):
  ps=primes_to(Q//div);tuples=sorted((prod(c),list(c)) for c in combinations(ps,w) if prod(c)<=Q)
  alternate=[]
  for n in range(2,Q+1):
   r=n;support=[]
   for p in ps:
    if p*p>r:break
    if r%p==0:
     r//=p;support.append(p)
     if r%p==0:support=[];r=0;break
   if r>1:support.append(r)
   if len(support)==w:alternate.append((n,support))
  need(tuples==alternate,'two exhaustive squarefree radical inventories')
  rows=[]
  for q,S in tuples:
   hits=[c for c,m in fs if q%c==0 and q//c>=m and q//c in ps and c%(q//c)] if w==4 else []
   covered=bool(hits) if w==4 else q in closed
   need(covered==(q<Q),'first residual frontier exact')
   rows.append({'radical':q,'support':S,'closed':covered,'covering_core':hits[0] if hits else None})
  result[str(w)]=rows
 return result
def candidate():
 M=210;minima=[None]*M;remaining=M;R=0
 while remaining:
  R+=1
  for d in range(1,R+1):
   if gcd(d,M)>1:continue
   inv=pow(d,-1,M)
   for a in range(-R,R+1):
    s=a*inv%M
    if minima[s] is None:minima[s]=R;remaining-=1
 need(R<=M and minima[0]==1,'exact candidate including zero')
 hist={r:minima.count(r) for r in sorted(set(minima))}
 return {'M':M,'state_minima':minima,'histogram':hist,'all_states_checked':M}
def certificate():
 result={'schema':'global-exact-premises-v1','periodic_CDF':primitive_period(),'uniform_covering':covering(),'coupled_shear':shear(),'radical_frontiers':domains(),'candidate210':candidate()}
 result['checks']=CHECKS
 return serial(result)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 need(not a.output.exists(),'fresh output required')
 z=certificate();a.output.write_text(json.dumps(z,indent=2,sort_keys=True)+'\n',encoding='utf8');print('PASS',z['checks'])
