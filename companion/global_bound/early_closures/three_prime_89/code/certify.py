"""Certificate assembly. Uses exact rational enclosures; never accepts a float."""
from pathlib import Path
from fractions import Fraction as F
from math import prod,gcd,isqrt
from functools import reduce
import csv,json,sys,hashlib,heapq
from exact import need,log_hi,log_lo,upper_tick,out_decimal,log_interval,alternate_log
S=10**9
HERE=Path(__file__).resolve().parent.parent
SPECIAL={(2,5,7):12000,(2,5,11):700,(2,5,13):300,(3,5,7):300}

def triples(root=HERE):
 rows=[tuple(int(r[k]) for k in ('r','s','p')) for r in csv.DictReader(open(root/'data/remaining_89.csv'))]
 need(len(rows)==89 and len(set(rows))==89,'exactly 89 distinct residual triples')
 for t in rows:need(2<=t[0]<t[1]<t[2]<=47 and all(all(p%d for d in range(2,isqrt(p)+1)) for p in t),'prime triples')
 return rows

def cases(root=HERE):return [(t,1) for t in triples(root)]+[(t,2) for t in SPECIAL]
def stem(t,k):return f'q{prod(t)}'+('' if k==1 else f'_k{k}')
def grid(root=HERE):
 rows=[list(map(int,l.split())) for l in (root/'data/grid.tsv').read_text().splitlines()];den=rows[0][0];rows=rows[1:]
 need(len(rows)==2688 and all(len(r)==8 and r[0]==i for i,r in enumerate(rows)),'grid complete')
 for i,r in enumerate(rows):
  _,lo,hi,tl,th,a,b,d=r
  need(0<lo<hi<2*den and a>0 and d>0,'positive intervals and seed bases')
  t=i%21-10;need(tl==max(-20,2*t-1) and th==min(20,2*t+1),'whole shear coverage')
  need((lo,hi)==tuple(rows[21*(i//21)][1:3]),'radial block identity')
  if i%21==0 and i:need(lo==rows[i-21][2],'radial intervals contiguous')
 need(3*rows[-1][2]**2>=4*den**2,'covers Gauss upper endpoint')
 return den,rows

def profile(t,k,profile_path,root=HERE):
 q=prod(t);K=q**k;nl=prod(p+1 for p in t);den,grows=grid(root)
 raw=[list(map(int,l.split())) for l in Path(profile_path).read_text().splitlines()]
 need(len(raw)==2688,'all cell bounds')
 for i,r in enumerate(raw):need(len(r)==5 and r[:3]==grows[i][:3] and 0<=r[4]<nl,'cell domain and mark')
 v=[F(max(r[3] for r in raw[21*j:21*(j+1)]),S) for j in range(128)]
 for i in range(126,-1,-1):v[i]=max(v[i],v[i+1])
 a=F(29,14);rho=F(raw[0][1],den)
 bulk=v[-1]+sum((v[j]-v[j+1])*min(F(1),a*F(raw[21*j][2],den)) for j in range(127))
 c=log_hi(q-1)+log_hi(K)/2-log_lo(rho)/2+F(5,2)*rho
 cusp=a*rho*(max(c-v[0],F(0))+F(1,2));total=bulk+cusp
 return dict(primes=list(t),q=q,k=k,transfer_index=K,labels=nl,primitive_branches=(K//q)*nl,
  cells=2688,witnesses=2688*(K//q)*nl,bulk_upper=str(bulk),cusp_upper=str(cusp),total_upper=str(total),
  upper_decimal=out_decimal(total),target_lower=str(F(2,85)*log_lo(q)),
  direct_radical_gate=total<F(2,85)*log_lo(q),monotone_profile_ticks=[int(x*S) for x in v])

def finite_domain(root=HERE):
 allrows=[]
 for t,cut in SPECIAL.items():
  q=prod(t);vals=[]
  def walk(i,m,es):
   if i==3:vals.append((m,tuple(es)));return
   m*=t[i];e=1
   while m<cut:walk(i+1,m,es+[e]);m*=t[i];e+=1
  walk(0,1,[]);vals.sort()
  seen={q};heap=[q];ind=[]
  while heap:
   m=heapq.heappop(heap);ind.append(m)
   for p in t:
    n=m*p
    if n<cut and n not in seen:seen.add(n);heapq.heappush(heap,n)
  need([m for m,e in vals]==ind,'two independent finite inventories')
  for m,e in vals:allrows.append(dict(M=m,q=q,r=t[0],s=t[1],p=t[2],exponents=':'.join(map(str,e)),base=m,length=1,base_cutoff=cut))
 return sorted(allrows,key=lambda r:r['M'])

def candidate():
 M=210;R=[];w=[]
 for s in range(M):
  best=max(1,min(s,M-s));aa=s if s<=M//2 else s-M;dd=1
  for d in range(1,M//2+1):
   if gcd(d,M)>1:continue
   a=s*d%M;a=a if a<=M//2 else a-M
   h=max(abs(a),d)
   if h<best:best=h;aa=a;dd=d
  R.append(best);w.append([s,aa,dd,best])
 need(all(d>0 and gcd(d,M)==1 and (a-s*d)%M==0 and h==max(abs(a),d) for s,a,d,h in w),'candidate witnesses')
 P=prod(R);need(P**170>210**18690,'candidate >2/85 exact integer gate')
 # Directly check every smaller possible denominator against the attained height.
 cnt=0
 for s,a,d,h in w:
  for e in range(1,h):
   if gcd(e,M)==1:cnt+=1;x=s*e%M;need(max(e,min(x,M-x))>=h,'candidate minimality')
 z=P;fac={}
 for p in range(2,max(R)+1):
  e=0
  while z%p==0:z//=p;e+=1
  if e:fac[str(p)]=e
 need(z==1,'complete product factorization')
 return dict(status='PASS',M=210,records=w,product_factors=fac,lower_denominator_checks=cnt,integer_comparison='P^170 > 210^18690',comparison_passed=True)

def finite_certificate(work,root=HERE):
 rows=finite_domain(root);out=[];logs_used=set()
 for r in rows:
  M=r['M'];hist=[tuple(map(int,l.split())) for l in (Path(work)/f'M{M}.hist').read_text().splitlines()]
  need(sum(n for h,n in hist)==M and len({h for h,n in hist})==len(hist),'full finite histogram mass')
  need(all(h>=1 and n>0 for h,n in hist),'positive heights and multiplicities')
  u=sum(n*log_hi(h) for h,n in hist)/M-log_lo(M)/2
  logs_used.update(h for h,n in hist);logs_used.add(M)
  need(u<F(2,85)*log_lo(M),'finite base-length gate');need(u<0,'negative complete witness mean')
  out.append({**r,'mean_upper':str(u),'upper_decimal':out_decimal(u),'normalized_upper':str(max(u,F(0))/log_lo(M)),'histogram':hist})
 checks=0
 for n in sorted(logs_used):
  lo,hi=alternate_log(n)
  need(lo>=log_lo(n) and hi<=log_hi(n),'secondary rational log interval')
  checks+=1
 return dict(status='PASS',all_means_negative=True,moduli=len(rows),representatives=sum(r['M']//2+1 for r in rows),states=sum(r['M'] for r in rows),independently_checked_integer_logs=checks,rows=out)

def assemble(profile_dir,finite,root=HERE):
 profiles=[];routes=[]
 for t,k in cases(root):
  pr=profile(t,k,Path(profile_dir)/(stem(t,k)+'.cells'),root);profiles.append(pr)
  u=F(pr['total_upper']);q=pr['q']
  if k==1 and t in SPECIAL:
   base_min=SPECIAL[t];ratio=max(F(0),u)/log_lo(base_min)
   routes.append(dict(route='large_base',q=q,k=k,base_min=base_min,normalized_upper=str(ratio)))
  else:
   need(pr['direct_radical_gate'],'uniform radical exclusion')
   ratio=max(F(0),u)/log_lo(q)
   routes.append(dict(route='all_moduli' if k==1 else 'N_ge_2',q=q,k=k,base_min=q,normalized_upper=str(ratio)))
 for r in finite['rows']:routes.append(dict(route='finite_N1',q=r['q'],M=r['M'],base_min=r['M'],normalized_upper=r['normalized_upper']))
 cap=F(233,10000);gamma=max(F(r['normalized_upper']) for r in routes)
 need(all(F(r['normalized_upper'])<cap for r in routes),'common 233/10000 gate')
 need(cap<F(2,85),'strict candidate gap')
 need(sum(p['direct_radical_gate'] for p in profiles if p['k']==1)==85,'85 immediate uniform closures')
 return dict(status='PASS',new_supports_closed=89,one_level_closed=85,special_supports=[prod(t) for t in SPECIAL],
  finite_profiles=93,cell_profiles=93*2688,witness_records=sum(p['witnesses'] for p in profiles),
  corner_checks=4*sum(p['witnesses'] for p in profiles),common_bound=str(gamma),common_bound_decimal=out_decimal(gamma),
  simple_bound=str(cap),candidate_lower='2/85',worst_routes=[r for r in routes if F(r['normalized_upper'])==gamma],
  all_new_supports_closed_all_valuations=True,all_three_prime_supports_closed_with_declared_imports=True,
  previous_exclusions_replayed=False,global_Cstar_identified=False,profiles=profiles,routes=routes)

def write_json(path,data):Path(path).write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
