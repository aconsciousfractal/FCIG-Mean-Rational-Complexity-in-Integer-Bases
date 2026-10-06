"""Fixed new portable definition/oracle smoke, not historical job replay."""
from collections import Counter
import json
from math import gcd,isqrt
from central_oracle import central_counts

MODULI=(2,3,9,6,10,44,352,900)


def minimum(M,s):
    best=max(1,min(s,M-s));d=2
    while d<best:
        if gcd(d,M)==1:
            A=s*d%M;best=min(best,max(d,min(A,M-A)))
        d+=1
    return best


def factors(M):
    out=[];p=2
    while p*p<=M:
        if M%p==0:
            e=0
            while M%p==0:M//=p;e+=1
            out.append((p,e))
        p+=1
    if M>1:out.append((M,1))
    return out


def main():
    rows=[]
    for M in MODULI:
        heights=[minimum(M,s) for s in range(M)];hist=sorted(Counter(heights).items());fs=factors(M)
        checks=0
        if len(fs)==2:
            (p,a),(q,b)=fs
            for T in range(1,isqrt(M)+1):
                actual=central_counts(M,p,q,a,b,T)
                expected=[[g,sum(h<=T and gcd(s,M)==g for s,h in enumerate(heights))] for g,c in actual['strata']]
                if actual['cdf']!=sum(h<=T for h in heights) or actual['strata']!=expected or actual['zero_count']!=1:
                    raise ValueError('fixed central definition mismatch')
                checks+=1
        if M==900:
            product=1
            for h,c in hist:product*=h**c
            if product*product<=M**M:raise ValueError('exact B30,N2 counterexample missing')
        rows.append(dict(M=M,states=M,central_threshold_checks=checks,histogram=hist))
    print(json.dumps(dict(status='PORTABLE_FIXED_DEFINITION_SMOKE_PASS',fixed_moduli=list(MODULI),definition_states=sum(MODULI),rows=rows,large_domain_or_historical_jobs_replayed=False),sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
