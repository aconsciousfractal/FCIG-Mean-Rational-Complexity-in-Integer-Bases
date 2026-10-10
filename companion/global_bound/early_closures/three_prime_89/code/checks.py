"""Independent checks of algebra, logarithms, cell geometry and physical banks."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,prod,isqrt
from collections import Counter
from functools import lru_cache
import json,sys
from certify import cases,triples,grid,stem,HERE
from exact import need,log_hi,log_lo,log_interval,ceildiv,alternate_log
FP=10**40

def raw_log_fp(a,b,terms=96):
 need(a>=b>0 and a<=2*b,'normalized secondary log')
 n=a-b;lo=n*FP//a;hi=ceildiv(n*FP,a);v=w=0
 for j in range(1,terms+1):
  v+=lo//j;w+=ceildiv(hi,j);lo=lo*n//a;hi=ceildiv(hi*n,a)
 return v,w+ceildiv(hi*a,(terms+1)*b)
L2,U2=raw_log_fp(2,1)
@lru_cache(None)
def log_fp(n,d=1):
 need(n>0 and d>0,'positive secondary log')
 k=n.bit_length()-d.bit_length();a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
 if a<b:k-=1;a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
 lo,hi=raw_log_fp(a,b)
 return (lo+k*L2,hi+k*U2) if k>=0 else (lo+k*U2,hi+k*L2)

def generate_branches(t,k):
 q=prod(t);K=q**k;items=[];counts=Counter();local=0
 for a in range(1,K+1):
  if K%a:continue
  d=K//a
  for b in range(a):
   kernel=[];image=[];bad=False
   for p in t:
    values=[(1,z) for z in range(p)]+[(0,1)]
    images=[(a*x%p,(b*x+d*y)%p) for x,y in values]
    zeros=[i for i,(x,y) in enumerate(images) if not x and not y]
    if len(zeros)>1:need(len(zeros)==p+1,'zero-rank HNF');bad=True;break
    need(len(zeros)==1,'one kernel line')
    nz=next(v for v in images if v!=(0,0))
    for x,y in images:
     if x or y:need((x*nz[1]-y*nz[0])%p==0,'all nonzero images proportional')
     local+=1
    kernel.append(zeros[0]);image.append(nz)
   if bad:continue
   label=0
   for v,p in zip(kernel,t):label=label*(p+1)+v
   counts[label]+=1;items.append((label,(a,0,b,d),tuple(image)))
 need(len(counts)==prod(p+1 for p in t) and all(n==K//q for n in counts.values()),'all projective kernels and exact multiplicities')
 items.sort(key=lambda v:(v[0],v[1]));return items,local

def tensor_direct_checks(t,branches):
 dims=[p+1 for p in t];nl=prod(dims);values=[(i*i+17*i+13)%101-50 for i in range(nl)]
 def decode(i):
  lab=[]
  for n in reversed(dims):lab.append(i%n);i//=n
  return tuple(reversed(lab))
 # Independent tensor built as three explicit sums on coordinate slices.
 x=values[:];stride=1
 for n in reversed(dims):
  out=x[:]
  for start in range(0,nl,stride*n):
   for off in range(stride):
    ids=[start+off+j*stride for j in range(n)]
    for idx in ids:out[idx]=sum(x[h] for h in ids if h!=idx)
  x=out;stride*=n
 chosen=range(nl) if nl<=300 else sorted({0,1,nl//4,nl//2,3*nl//4,nl-2,nl-1})
 decoded=[decode(i) for i in range(nl)];n=0
 for i in chosen:
  v=sum(values[h] for h in range(nl) if all(a!=b for a,b in zip(decoded[i],decoded[h])))
  need(v==x[i],'signed tensor versus direct admissible sum');n+=1
 return n

def run_checks(profile_dir,root=HERE):
 den,gr=grid(root);tab=list(map(int,(root/'data/log_table.txt').read_text().split()))
 need(len(tab)==2049 and all(x==0 for x in tab[:1025]),'short log-table layout')
 for i in range(1025,2049):
  lo,hi=log_fp(i,1024)
  need(hi*10**9<=2*tab[i]*FP and lo*10**9>=2*(tab[i]-1)*FP,'two-sided log-table enclosure')
 # Secondary fixed-point series checked against full rational sums.
 for n,d in [(j+100,j+1) for j in range(120)]+[(1,j) for j in range(1,32)]:
  lo,hi=alternate_log(F(n,d));a,b=log_fp(n,d)
  need(F(a,FP)<=lo<=hi<=F(b,FP),'fixed-point secondary versus rational series')
 branch_total=local_total=tensor_total=audits=corners=interiors=0
 for t,k in cases(root):
  q=prod(t);K=q**k;br,local=generate_branches(t,k);local_total+=local;branch_total+=len(br)
  tensor_total+=tensor_direct_checks(t,br)
  rows=[list(map(int,l.split())) for l in (Path(profile_dir)/(stem(t,k)+'.audit')).read_text().splitlines()]
  need(len(rows)==60,'audit domain complete')
  for row in rows:
   c,h,ci,cj,u,v,cap=row;_,H,image=br[h]
   need(u==ci*H[0]+cj*H[2] and v==ci*H[1]+cj*H[3],'independent parent coefficient transport')
   need(all((ci*a+cj*b)%p for p,(a,b) in zip(t,image)),'full mark admissibility')
   _,lo,hi,tl,th,*_=gr[c];mx=F(0)
   for r in (lo,hi):
    for z in (tl,th):
     sq=F(r,den*K)*F(40*u+z*v,40)**2+F(den*v*v,r*K)
     a,b=log_fp(sq.numerator,sq.denominator)
     need(2*cap*FP>=b*10**9,'C++ corner bound versus independent logarithm')
     mx=max(mx,sq);corners+=1
   # Interior rationals validate the proved two-stage convexity contract.
   for u0,v0 in ((1,3),(2,2),(3,1)):
    r=F(u0*lo+(4-u0)*hi,4);z=F(v0*tl+(4-v0)*th,4)
    sq=r/(den*K)*(F(u)+z*v/40)**2+F(den*v*v,K)/r
    need(sq<=mx,'interior convexity regression');interiors+=1
   audits+=1
 # Physical graph-lattice shortest-vector CDF, computed without Gauss reduction.
 states=0;jumps=0
 for L in range(1,401):
  vals=[]
  for s in range(L):
   val=L*L
   for d in range(1,isqrt(2*L)+1):
    a=s*d%L;a=min(a,L-a);val=min(val,a*a+d*d)
   vals.append(val);states+=1
  hist=Counter(vals);cnt=0
  for norm,count in sorted(hist.items()):cnt+=count;need(14*cnt<=29*norm,'short-vector CDF inequality at every jump');jumps+=1
 lift=0
 for t in ((2,5,7),(2,5,11),(2,5,13),(3,5,7)):
  for k in (1,2):
   K=prod(t)**k
   for L in range(1,11):
    seen=bytearray(K*L)
    for b in range(K):
     for s in range(L):
      n=s+b*L;need(n%L==s and n//L==b and not seen[n],'physical lift inverse');seen[n]=1;lift+=1
    need(all(seen),'all states of physical lift')
 # No equality is assumed between infinity and Euclidean heights.
 def height(M,q,s,squared=False):
  best=M*M+1 if squared else M
  for d in range(1,M+1):
   if gcd(d,q)>1:continue
   a=s*d%M;a=min(a,M-a);h=a*a+d*d if squared else max(a,d)
   best=min(best,h)
  return best
 need(height(70,70,2)==2 and height(70,70,2,True)==5,'Euclidean relaxation is strict')
 need(height(210,70,70)==3 and height(210,210,70)==70,'extra-prime relaxed model is not original')
 return dict(status='PASS',log_table_entries=1024,secondary_rational_series_cases=151,
  algebraic_branches=branch_total,local_projective_checks=local_total,signed_tensor_direct_outputs=tensor_total,
  audit_witnesses=audits,independent_corner_logs=corners,interior_rational_tests=interiors,
  shortest_vector_bank_states=states,shortest_vector_CDF_jumps=jumps,physical_lift_states=lift,
  euclidean_equals_original_refuted=True,extra_prime_identification_refuted=True,
  independent_checks_not_a_second_full_replay_of_every_geometric_cost=True)

if __name__=='__main__':
 print(json.dumps(run_checks(sys.argv[1]),indent=2))
