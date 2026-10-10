"""Independent algebra, logarithm and geometric audits, and complete small banks.
The finite-state consumer is exhaustive; the second geometric audit is stratified.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,prod,isqrt
from collections import Counter
from functools import lru_cache
from exact import need,ceildiv,alternate_log,log_hi,log_lo,lower_tick,upper_tick,sieve
from assemble import grid,profile
from core import factors
FP=10**40

def reduced_fp(a,b,terms=96):
 need(b<=a<=2*b,'secondary binary normalization');n=a-b;lo=n*FP//a;hi=ceildiv(n*FP,a);s=t=0
 for j in range(1,terms+1):
  s+=lo//j;t+=ceildiv(hi,j);lo=lo*n//a;hi=ceildiv(hi*n,a)
 return s,t+ceildiv(hi*a,(terms+1)*b)
L2,U2=reduced_fp(2,1)
@lru_cache(None)
def log_fp(n,d=1):
 need(n>0 and d>0,'positive secondary logarithm');k=n.bit_length()-d.bit_length()
 a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
 if a<b:k-=1;a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
 lo,hi=reduced_fp(a,b)
 return (lo+k*L2,hi+k*U2) if k>=0 else (lo+k*U2,hi+k*L2)

def input_checks(root):
 tab=list(map(int,(root/'data/log_table.txt').read_text().split()))
 need(len(tab)==2049 and all(z==0 for z in tab[:1025]),'mantissa table layout')
 for i in range(1025,2049):
  lo,hi=log_fp(i,1024)
  need(hi*10**9<=2*tab[i]*FP and lo*10**9>=2*(tab[i]-1)*FP,'independent mantissa enclosure')
 seeds=[list(map(int,l.split())) for l in (root/'data/log_seeds.tsv').read_text().splitlines()]
 need(len(seeds)==1024,'all logarithm seeds')
 for i,(n,lo,hi) in enumerate(seeds,1):
  a,b=alternate_log(n)
  need(n==i and lo==lower_tick(a,10**17) and hi==upper_tick(b,10**17),'secondary seeds regenerated')
 for n,d in [(j+100,j+1) for j in range(120)]+[(1,j) for j in range(1,32)]:
  a,b=alternate_log(F(n,d));lo,hi=log_fp(n,d)
  need(F(lo,FP)<=a<=b<=F(hi,FP),'fixed-point secondary against exact rational series')
 return dict(status='PASS',mantissa_entries=1024,dense_log_seed_entries=1024,secondary_series_rational_checks=151)

def ps_of(p):return factors(p)

def algebra_check(p,k,prefix):
 ps=ps_of(p);q=prod(ps);K=q*k;nl=prod(a+1 for a in ps);br=[];counts=Counter();local=0
 # The output is compared to a second full HNF enumeration, using explicit
 # normalized columns and linear equations rather than the C++ inverse routine.
 for a in range(1,K+1):
  if K%a:continue
  d=K//a
  for b in range(a):
   label=0;img=[];bad=False
   for pr in ps:
    if a%pr==b%pr==d%pr==0:bad=True;break
    # Independent projective-line scan gives the single zero and all images.
    cols=[(1,z) for z in range(pr)]+[(0,1)]
    zeros=[];base=None
    for j,(x,y) in enumerate(cols):
     z=(a*x%pr,(b*x+d*y)%pr)
     if z==(0,0):zeros.append(j)
     else:
      if base is None:base=z
      need((z[0]*base[1]-z[1]*base[0])%pr==0,'same nonzero image line')
     local+=1
    need(len(zeros)==1,'one kernel line');label=label*(pr+1)+zeros[0];img.append(base)
   if bad:continue
   br.append((label,a,b,d,tuple(img)));counts[label]+=1
 br.sort(key=lambda z:(z[0],z[1:4]));need(len(br)==nl*K//q and len(counts)==nl and all(v==K//q for v in counts.values()),'all kernels and multiplicities')
 stored=[list(map(int,l.split())) for l in Path(str(prefix)+'.branches').read_text().splitlines()]
 need(len(stored)==len(br),'full algebraic branch count')
 for row,ex in zip(stored,br):
  need(row[:4]==list(ex[:4]),'HNF ordering and all representatives')
  for pr,(u,v) in zip(ps,ex[4]):need((row[4]%pr or row[5]%pr) and (row[4]*v-row[5]*u)%pr==0,'stored image mark is on the required line')
 raw=[list(map(int,l.split())) for l in Path(str(prefix)+'.tensor').read_text().splitlines()]
 need(len(raw)==nl and all(z[0]==i and z[1]==K for i,z in enumerate(raw)),'tensor constants normalized on every mark')
 def decode(j):
  out=[]
  for pr in reversed(ps):out.append(j%(pr+1));j//=pr+1
  return tuple(reversed(out))
 dec=[decode(i) for i in range(nl)];chosen=range(nl) if nl<=100 else sorted({0,1,nl//4,nl//2,3*nl//4,nl-2,nl-1});tested=0
 v1=[(i*i+17*i+13)%101-50 for i in range(nl)];v2=[i%7-3 for i in range(nl)]
 for i in chosen:
  admiss=[j for j in range(nl) if all(x!=y for x,y in zip(dec[i],dec[j]))]
  need(len(admiss)==q,'full admissible-label cardinality')
  need(raw[i][2]==(K//q)*sum(v1[j] for j in admiss) and raw[i][3]==(K//q)*sum(v2[j] for j in admiss),'two signed tensors against direct sums')
  tested+=2
 return br,dict(branches=len(br),local_projective_checks=local,tensor_constant_outputs=nl,signed_direct_outputs=tested)

def geometric_check(p,k,ag,prefix,gr,den,br):
 q=prod(ps_of(p));K=q*k
 rows=[list(map(int,l.split())) for l in Path(str(prefix)+'.audit').read_text().splitlines()]
 need(len(rows)==96,'stratified geometric audit length');logs=interiors=0
 for row in rows:
  c,h,ci,cj,u,v=row[:6];caps=row[6:];need(len(caps)==2*ag+1,'all audit angular nodes')
  _,a,b,d,image=br[h]
  need((u,v)==(ci*a+cj*b,cj*d),'coefficient transport')
  need(all((ci*x+cj*y)%pr for pr,(x,y) in zip(ps_of(p),image)),'audited unit witness')
  _,lo,hi,tl,th,*_=gr[c]
  for n,cap in zip(range(-ag,ag+1),caps):
   maxsq=F(0)
   for rho in (F(lo,den),F(hi,den)):
    for t in (F(tl,40),F(th,40)):
     xx=rho*(u+t*v);aa=ag*xx+n*v;bb=-n*xx+ag*v;sq=max(aa*aa,bb*bb)/(F(K)*rho*(ag*ag+n*n));A,B=log_fp(sq.numerator,sq.denominator)
     need(2*cap*FP>=B*10**9,'independent corner logarithm enclosure');maxsq=max(maxsq,sq);logs+=1
   for wi,wj in ((1,3),(2,2),(3,1)):
    rho=F(wi*lo+(4-wi)*hi,4*den);t=F(wj*tl+(4-wj)*th,160)
    xx=rho*(u+t*v);aa=ag*xx+n*v;bb=-n*xx+ag*v;sq=max(aa*aa,bb*bb)/(F(K)*rho*(ag*ag+n*n))
    need(sq<=maxsq,'whole-cell convexity at rational interior points');interiors+=1
  log_fp.cache_clear()
 return dict(witnesses=len(rows),independent_corner_logarithms=logs,rational_interior_checks=interiors)

@lru_cache(None)
def direct(M,q):
 vals=[];tests=0
 for s in range(M):
  best=max(1,min(s,M-s));d=1
  while d<best:
   if gcd(d,q)==1:
    a=s*d%M;best=min(best,max(d,min(a,M-a)));tests+=1
   d+=1
  vals.append(best)
 return vals,tests

def gauss(L,s):
 u=(L,0);v=(s,1)
 for _ in range(100):
  nu=sum(x*x for x in u);nv=sum(x*x for x in v)
  if nu>nv:u,v=v,u;nu=nv
  t=(2*sum(x*y for x,y in zip(u,v))+nu)//(2*nu)
  if t==0:return u,v
  v=tuple(y-t*x for x,y in zip(u,v))
 raise ArithmeticError('Gauss termination')

