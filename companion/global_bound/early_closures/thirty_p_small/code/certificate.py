"""Rational integration and complete covering of the selected four-prime sectors."""
from pathlib import Path
from fractions import Fraction as F
from math import prod,gcd,isqrt
from collections import Counter
import struct,heapq,json
from exact import need,log_lo,log_hi,out_decimal,alternate_log
from new_math import large_prime_gates
PS=(11,13,17,19,23,29,31)
ETA=F(233,10000);S=10**9;DEN=2**52

def configurations():return [(5,False,16,3)]+[(p,True,8,2) for p in PS]
def stem(p,four):return f'q{(30 if four else 6)*p}'

def profile(p,four,k,prefix):
 q=(30 if four else 6)*p;nl=(72 if four else 12)*(p+1)
 c=[list(map(int,l.split())) for l in Path(str(prefix)+f'.cells{k}').read_text().splitlines()]
 r=[list(map(int,l.split())) for l in Path(str(prefix)+f'.radial{k}').read_text().splitlines()]
 need(len(c)==2688 and len(r)==128,'all profile cells')
 for i,z in enumerate(c):need(len(z)==4 and z[0]==i and 0<=z[2]<nl and 0<=z[3]<33,'cell metadata')
 for j,z in enumerate(r):
  need(len(z)==5 and z[0]==j-125 and 0<z[1]<z[2] and z[3]==max(x[1] for x in c[j*21:(j+1)*21]),'radial maximum')
  need(z[4] in [x[0] for x in c[j*21:(j+1)*21] if x[1]==z[3]],'actual maximizing cell')
  if j:need(z[1]==r[j-1][2],'no radial gaps')
 need(3*r[-1][2]**2>=4*DEN**2,'Gauss domain covered')
 v=[F(z[3],S) for z in r]
 for j in range(126,-1,-1):v[j]=max(v[j],v[j+1])
 a=F(29,14);rho=F(r[0][1],DEN)
 bulk=v[-1]+sum((v[j]-v[j+1])*min(1,a*F(r[j][2],DEN)) for j in range(127))
 t=log_hi(q-1)+k*log_hi(q)/2-log_lo(rho)/2+4*rho
 cusp=a*rho*(max(t-v[0],F(0))+F(1,2));total=bulk+cusp
 return dict(p=p,four=four,q=q,k=k,labels=nl,directions=32 if four else 64,
             bulk_upper=str(bulk),cusp_upper=str(cusp),total_upper=str(total),upper_decimal=out_decimal(total),
             monotone_profile_ticks=[int(x*S) for x in v],radical_gate=total<ETA*log_lo(q))

def collect(profiles_dir):
 rows=[]
 for p,four,_,mk in configurations():
  for k in range(1,mk+1):rows.append(profile(p,four,k,Path(profiles_dir)/stem(p,four)))
 a={1:F(57,200),2:F(1,8),3:F(31,500)}
 for z in rows:
  if not z['four']:need(F(z['total_upper'])<a[z['k']],'rerun marked-30 cap')
 return rows

def inventory(profiles):
 by={(z['q'],z['k']):F(z['total_upper']) for z in profiles};rows=[];routes=[];cuts={}
 for p in PS:
  q=30*p;u=by[q,1]
  if u<ETA*log_lo(q):routes.append(dict(p=p,q=q,route='all_N',base_min=q,normalized_upper=str(max(0,u)/log_lo(q))));continue
  v=by[q,2];need(v<ETA*log_lo(q),'N>=2 uniform gate')
  routes.append(dict(p=p,q=q,route='N>=2',base_min=q,normalized_upper=str(max(0,v)/log_lo(q))))
  cut=q
  while not u<ETA*log_lo(cut):cut*=2
  cuts[q]=cut
  routes.append(dict(p=p,q=q,route='N=1_large_B',base_min=cut,normalized_upper=str(max(0,u)/log_lo(cut))))
  def rec(j,n,ee):
   if j==4:rows.append(dict(M=n,q=q,p=p,base=n,N=1,exponents=ee,base_cutoff=cut));return
   pr=(2,3,5,p)[j];e=1;n*=pr
   while n<cut:rec(j+1,n,ee+[e]);n*=pr;e+=1
  rec(0,1,[])
  expected=sorted(z['M'] for z in rows if z['q']==q)
  scanned=[]
  for n in range(q,cut,q):
   z=n
   for pr in (2,3,5,p):
    while z%pr==0:z//=pr
   if z==1:scanned.append(n)
  need(expected==scanned,'independent finite-inventory scan')
 need(len({z['M'] for z in rows})==len(rows),'disjoint radicals have disjoint moduli')
 return sorted(rows,key=lambda z:z['M']),routes,cuts

def finite_certificate(folder,rows):
 out=[];logs=set()
 for z in rows:
  M=z['M'];hs=[tuple(map(int,l.split())) for l in (Path(folder)/f'M{M}.hist').read_text().splitlines()]
  need(sum(n for h,n in hs)==M and len({h for h,n in hs})==len(hs),'complete finite mass')
  need(all(1<=h<=M//2 and n>0 for h,n in hs),'histogram domain')
  u=sum(n*log_hi(h) for h,n in hs)/M-log_lo(M)/2;ratio=max(F(0),u)/log_lo(M)
  need(ratio<ETA,'finite base comparison')
  out.append({**z,'mean_upper':str(u),'mean_upper_decimal':out_decimal(u),'normalized_upper':str(ratio),'histogram':hs})
  logs.add(M);logs.update(h for h,n in hs)
 for n in sorted(logs):
  lo,hi=alternate_log(n);need(log_lo(n)<=lo<=hi<=log_hi(n),'all finite integer logarithms independently enclosed')
 return dict(status='PASS',moduli=len(rows),states=sum(z['M'] for z in rows),representatives=sum(z['M']//2+1 for z in rows),
             independent_log_arguments=len(logs),rows=out)

def physical_bank(q,prefix,refdir,maxk):
 import arithmetic as ar
 vals={};roots=[];points=conditional=0
 for k in range(1,min(maxk,2)+1):
  M=q**k;raw=(Path(refdir)/f'M{M}_q{q}.min').read_bytes();need(len(raw)==4*M,'definition-based complete minimum array')
  R=struct.unpack('<'+'I'*M,raw);need(all(1<=h<=max(1,M//2) for h in R),'reference height domain');vals[k]=R
 if q==30:
  raw=(Path(refdir)/'M27000_q30.min').read_bytes();need(len(raw)==4*27000,'cubic mark30 reference');vals[3]=struct.unpack('<'+'I'*27000,raw)
 ticks={h:ar.logs(h)[1] for R in vals.values() for h in set(R)}
 for line in Path(str(prefix)+'.points').read_text().splitlines():
  k,bound=map(int,line.split())
  if k==0:need(bound>=0,'identity observable');continue
  if k not in vals:continue
  M=q**k;u=F(sum(ticks[h] for h in vals[k]),M*ar.LS)-F(ar.logs(M)[0],2*ar.LS)
  need(u<=F(bound,ar.S),'root envelope against every exact modular minimum');roots.append(dict(k=k,M=M,mean_upper=str(u)))
 for line in Path(str(prefix)+'.bank').read_text().splitlines():
  k,s,bound=map(int,line.split());need(0<=s<q,'bank state')
  if k==0:
   u=F(ticks[vals[1][s]],ar.LS)-F(ar.logs(q)[0],2*ar.LS);points+=1
  elif k==1:
   R=vals[2];u=F(sum(ticks[R[s+q*b]] for b in range(q)),q*ar.LS)-F(ar.logs(q*q)[0],2*ar.LS);conditional+=1
  else:continue
  need(u<=F(bound,ar.S),'conditional envelope against complete descendants')
 return dict(q=q,root_checks=roots,pointwise_checks=points,conditional_checks=conditional,
             independently_enumerated_modular_states=sum(len(R) for R in vals.values()))

def closure(profiles,finite):
 rows,routes,cuts=inventory(profiles);need([r['M'] for r in rows]==[r['M'] for r in finite['rows']],'exact finite residual')
 routes.extend(dict(q=r['q'],route='finite_N1',M=r['M'],base_min=r['M'],normalized_upper=r['normalized_upper']) for r in finite['rows'])
 gamma=max(F(r['normalized_upper']) for r in routes);need(gamma<ETA<F(2,85),'strict common gap for seven sectors')
 large=large_prime_gates();globalbound=max([gamma]+[F(r['delta_upper']) for r in large['new_unbounded_regions']]+[F(large['large_base_gate']['delta_upper'])])
 need(globalbound<ETA,'all new regions share the declared gap')
 return dict(status='PASS',closed_primes=list(PS),closed_radicals=[30*p for p in PS],
             all_valuations_and_lengths_for_selected_primes=True,base_cutoffs={str(k):v for k,v in cuts.items()},
             finite_moduli=len(rows),finite_states=finite['states'],profiles=profiles,routes=routes,
             selected_sectors_common_upper=str(gamma),selected_sectors_common_upper_decimal=out_decimal(gamma),
             all_new_regions_common_upper=str(globalbound),all_new_regions_common_upper_decimal=out_decimal(globalbound),
             simple_common_upper='233/10000',candidate_lower='2/85',large_prime_regions=large,
             entire_30p_family_closed=False,global_Cstar_identified=False,
             historical_three_prime_closure_used=False,historical_210_closure_used=False)
