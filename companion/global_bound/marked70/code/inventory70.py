"""Exhaustive N=1 residual for rad(B)=70p; all decisions use rational intervals."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import csv, json
from exact import need, sieve, prime_trial, phi_hi, log_lo, out_decimal

BMAX=3_000_000
PMAX=1000
PMIN=73
DELTA=F(2,85)
A1=F(7,50)
A2=F(0)
ERR=F(51480,2107)

@lru_cache(None)
def penalty(p):
    return phi_hi(F(24,7*p))

def upper(M,p):
    # This extra use of k=2 includes N=1 with unequal positive valuations.
    a=A2 if M%(70**2)==0 else A1
    return a+penalty(p)+ERR/isqrt(M)

def make():
    ps=[p for p in sieve(PMAX-1) if p>=PMIN]
    need(ps==[p for p in range(PMIN,PMAX) if prime_trial(p)],'prime list by trial division')
    cofactors=[]
    def rec(i,m,ee):
        if i==3:
            cofactors.append((m,ee));return
        pr=(2,5,7)[i];n=m*pr;e=1
        while n*PMIN<BMAX:
            rec(i+1,n,ee+[e]);n*=pr;e+=1
    rec(0,1,[]);cofactors.sort();rows=[]
    for p in ps:
        pp=p;e=1
        while 70*pp<BMAX:
            for s,ee in cofactors:
                M=s*pp
                if M>=BMAX:break
                u=upper(M,p);gate=DELTA*log_lo(M)
                rows.append(dict(M=M,p=p,q=70*p,exponents=':'.join(map(str,ee+[e])),
                    profile_k=2 if M%4900==0 else 1,
                    route='analytic' if u<gate else 'witness',
                    mean_upper=str(u),normalized_upper=str(u/log_lo(M))))
            pp*=p;e+=1
    rows.sort(key=lambda r:r['M'])
    need(len({r['M'] for r in rows})==len(rows),'unique modulus')
    # Independently scan every multiple of 70, instead of multiplying powers.
    check=[]
    for M in range(70*PMIN,BMAX,70):
        v=M
        for s in (2,5,7):
            while v%s==0:v//=s
        if v==1:continue
        p=v
        for s in range(2,isqrt(v)+1):
            if v%s==0:p=s;break
        if not PMIN<=p<PMAX or not prime_trial(p):continue
        while v%p==0:v//=p
        if v==1:check.append((M,p))
    need(check==[(r['M'],r['p']) for r in rows],'independent cofactor scan')
    gates=[]
    for name,p,M,b,a in [
        ('large_p',PMAX,70*PMAX,70*PMAX,A1),
        ('N_ge_2',PMIN,(70*PMIN)**2,70*PMIN,A2),
        ('large_base',PMIN,BMAX,BMAX,A1),
    ]:
        u=a+penalty(p)+ERR/isqrt(M);ratio=u/log_lo(b)
        need(ratio<DELTA,'infinite region '+name)
        gates.append(dict(name=name,p_min=p,M_min=M,base_min=b,
            relaxed_cap=str(a),mean_upper=str(u),normalized_upper=str(ratio),
            normalized_upper_decimal=out_decimal(ratio),margin_to_delta=str(DELTA-ratio)))
    witnesses=[r for r in rows if r['route']=='witness']
    meta=dict(status='PASS',base_upper_exclusive=BMAX,prime_range=[PMIN,PMAX],
        prime_count=len(ps),cofactor_count=len(cofactors),all_moduli=len(rows),
        analytic_moduli=len(rows)-len(witnesses),witness_moduli=len(witnesses),
        witness_states=sum(r['M'] for r in witnesses),
        witness_representatives=sum(r['M']//2+1 for r in witnesses),
        maximum_witness_modulus=max(r['M'] for r in witnesses),infinite_regions=gates,
        second_profile_used_at_N1=sum(r['profile_k']==2 for r in rows),
        finite_moduli_are_bases_at_N1=True,other_four_prime_supports_closed=False,
        all_five_prime_supports_closed=False)
    return rows,meta

def write(root):
    rows,meta=make();out=root/'results';out.mkdir(parents=True,exist_ok=True)
    for name,chosen in [('finite_inventory.csv',rows),('witness_inventory.csv',[r for r in rows if r['route']=='witness'])]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(chosen)
    (out/'inventory_certificate.json').write_text(json.dumps(meta,sort_keys=True,indent=2)+'\n')
    (out/'finite_domain.txt').write_text(''.join(f"{r['M']} {r['p']}\n" for r in rows if r['route']=='witness'))
    return meta

if __name__=='__main__':
    print(json.dumps(write(Path(__file__).resolve().parent.parent),indent=2),flush=True)
