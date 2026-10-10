"""Exact finite residual after the uniform one-step block and prime-erasure bound."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,gcd
from functools import lru_cache
import csv,json,sys
from exact import need,sieve,prime_trial,phi_hi,log_lo,out_decimal
BMAX=4_000_000; PMAX=20_000; PMIN=211
A1=F(7,25);A2=F(61,500);DELTA=F(2,85); ERR=F(7304,399)

def upper(M,p):return A1+phi_hi(F(8,3*p))+ERR/isqrt(M)

def make():
 ps=[p for p in sieve(PMAX-1) if p>=PMIN]
 need(ps==[p for p in range(PMIN,PMAX) if prime_trial(p)],'prime list second enumeration')
 co=[]
 def sm(i,m,ee):
  if i==3:
   if m*PMIN<BMAX:co.append((m,ee))
   return
  pr=(2,3,5)[i];n=m*pr;e=1
  while n*PMIN<BMAX:sm(i+1,n,ee+[e]);n*=pr;e+=1
 sm(0,1,[]);co.sort();rows=[]
 for p in ps:
  pp=p;e=1
  while 30*pp<BMAX:
   for s,ee in co:
    M=s*pp
    if M>=BMAX:break
    cap=upper(M,p);gate=DELTA*log_lo(M);typ='analytic' if cap<gate else 'witness'
    rows.append(dict(M=M,p=p,q=30*p,exponents=':'.join(map(str,ee+[e])),route=typ,mean_upper=str(cap),normalized_upper=str(cap/log_lo(M))))
   pp*=p;e+=1
 rows.sort(key=lambda r:r['M']);need(len({r['M'] for r in rows})==len(rows),'no duplicate modulus')
 # Independent scan of every multiple of thirty below BMAX.
 check=[]
 for M in range(30*PMIN,BMAX,30):
  v=M
  for s in (2,3,5):
   while v%s==0:v//=s
  if v==1:continue
  p=v
  for s in range(2,isqrt(v)+1):
   if v%s==0:p=s;break
  if not PMIN<=p<PMAX or not prime_trial(p):continue
  while v%p==0:v//=p
  if v==1:check.append((M,p))
 need(check==[(r['M'],r['p']) for r in rows],'power products and cofactor scan agree')
 gates=[]
 for name,p,M,base,a in [('large_p',20000,30*20000,30*20000,A1),('N_ge_2',211,(30*211)**2,30*211,A2),('large_base',211,BMAX,BMAX,A1)]:
  u=a+phi_hi(F(8,3*p))+ERR/isqrt(M);ratio=u/log_lo(base)
  need(ratio<DELTA,'infinite region gate '+name)
  gates.append(dict(name=name,p_min=p,M_min=M,base_min=base,mark30_cap=str(a),mean_upper=str(u),normalized_upper=str(ratio),normalized_upper_decimal=out_decimal(ratio),margin_to_delta=str(DELTA-ratio)))
 witnesses=[r for r in rows if r['route']=='witness']
 return rows,dict(status='PASS',base_upper_exclusive=BMAX,prime_range=[PMIN,PMAX],all_moduli=len(rows),analytic_moduli=len(rows)-len(witnesses),witness_moduli=len(witnesses),witness_states=sum(r['M'] for r in witnesses),witness_representatives=sum(r['M']//2+1 for r in witnesses),maximum_witness_modulus=max(r['M'] for r in witnesses),infinite_regions=gates)

if __name__=='__main__':
 root=Path(__file__).resolve().parent.parent;rows,meta=make()
 for name,chosen in [('finite_inventory.csv',rows),('witness_inventory.csv',[r for r in rows if r['route']=='witness'])]:
  with (root/'results'/name).open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(chosen)
 (root/'results/inventory_certificate.json').write_text(json.dumps(meta,sort_keys=True,indent=2)+'\n')
 print(json.dumps(meta,indent=2),flush=True)
