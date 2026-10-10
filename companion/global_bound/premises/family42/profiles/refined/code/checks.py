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

def cases():return [(p,1,4) for p in sieve(210) if p>=37]+[(1,1,16),(1,2,16)]
def stem(p,k):return f'q{42*p}' if p>1 else f'q42_k{k}'
def ps_of(p):return (2,3,7,p) if p>1 else (2,3,7)

def algebra_check(p,k,prefix):
 ps=ps_of(p);q=prod(ps);K=q**k;nl=prod(a+1 for a in ps);br=[];counts=Counter();local=0
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
 q=prod(ps_of(p));K=q**k
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
    for t in (F(tl,80),F(th,80)):
     xx=rho*(u+t*v);aa=ag*xx+n*v;bb=-n*xx+ag*v;sq=max(aa*aa,bb*bb)/(F(K)*rho*(ag*ag+n*n));A,B=log_fp(sq.numerator,sq.denominator)
     need(2*cap*FP>=B*10**9,'independent corner logarithm enclosure');maxsq=max(maxsq,sq);logs+=1
   for wi,wj in ((1,3),(2,2),(3,1)):
    rho=F(wi*lo+(4-wi)*hi,4*den);t=F(wj*tl+(4-wj)*th,320)
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

def physical_checks(profiles,root):
 den,gr=grid(root);states=den_tests=parents=lift=0;rows=[]
 for p,k,Ls in [(1,1,(1,2,3,5,7)),(1,2,(1,2,3,5,7)),(37,1,(1,2,3,5,7,11)),(199,1,(1,2,3))]:
  q=prod(ps_of(p));K=q**k;pr=profile(p,profiles/(stem(p,k)+'.cells'),root,k,16 if p==1 else 4)
  caps=[list(map(int,l.split()))[3] for l in (profiles/(stem(p,k)+'.cells')).read_text().splitlines()]
  for L in Ls:
   M=K*L;R,cnt=direct(M,q);states+=M;den_tests+=cnt
   hh={h:log_hi(h) for h in set(R)};logM=log_lo(M)
   mean=sum(hh[h] for h in R)/M-logM/2;need(mean<=F(pr['total_upper']),'uniform mean against complete definition-based bank')
   for s in range(L):
    u,v=gauss(L,s);nu=sum(x*x for x in u);dot=sum(x*y for x,y in zip(u,v));rho=F(nu,L);t=F(dot,nu)
    indexes=[i for i,r in enumerate(gr) if F(r[1],den)<=rho<=F(r[2],den) and F(r[3],80)<=t<=F(r[4],80)]
    need(indexes,'physical geometry is inside the grid');c=indexes[0]
    child=[s+b*L for b in range(K)]
    need(len(set(child))==K and all(j%L==s for j in child),'physical fibre bijection');lift+=K
    value=sum(hh[R[j]] for j in child)/K-logM/2;need(value<=F(caps[c],10**9),'full conditional mean below continuous cell bound');parents+=1
   rows.append(dict(q=q,k=k,L=L,M=M,upper_exact_minimum_mean=str(mean),denominators_tested=cnt))
  direct.cache_clear()
 # Small-vector CDF without Gauss, checked at every jump of every L<=400.
 small=jump=0
 for L in range(1,401):
  vals=[]
  for s in range(L):
   z=L*L
   for d in range(1,isqrt(2*L)+1):a=s*d%L;a=min(a,L-a);z=min(z,a*a+d*d)
   vals.append(z);small+=1
  acc=0
  for n,h in sorted(Counter(vals).items()):acc+=h;need(14*acc<=29*n,'short-vector distribution at every jump');jump+=1
 return dict(status='PASS',complete_state_modulus_incidences=states,direct_denominator_tests=den_tests,conditional_parent_means=parents,physical_lift_states=lift,short_vector_states=small,short_vector_jumps=jump,rows=rows)

def run(profiles,root):
 den,gr=grid(root);inp=input_checks(root);totals=Counter();per=[]
 for p,k,ag in cases():
  pref=profiles/stem(p,k);br,al=algebra_check(p,k,pref);geo=geometric_check(p,k,ag,pref,gr,den,br)
  for d in (al,geo):totals.update(d)
  per.append(dict(p=p,k=k,algebra=al,geometry=geo))
  print('PASS independent p',p,'k',k,flush=True)
 return dict(status='PASS',inputs=inp,totals=dict(totals),per_profile=per,physical=physical_checks(profiles,root),not_second_exhaustive_geometric_evaluation=True)
