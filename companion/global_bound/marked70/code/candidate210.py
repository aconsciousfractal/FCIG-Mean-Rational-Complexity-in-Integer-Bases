"""Exact definition-based certificate for Delta(210,1)>2/85."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from math import gcd
import json
from exact import need, log_lo, log_hi, out_decimal

ROOT=Path(__file__).resolve().parent.parent
FACTORS={2:76,3:57,5:25,7:21,11:24,13:26,17:20,19:26,23:20,
         29:4,31:4,37:8,41:2,53:2}

def certificate():
    M=210
    # d=1 always gives height <=105, so d>105 cannot improve the answer.
    ds=[d for d in range(1,106) if gcd(d,M)==1]
    values=[]
    for s in range(M):
        values.append(min(max(d,min((s*d)%M,M-(s*d)%M)) for d in ds))
    hist=dict(sorted(Counter(values).items()))
    P=1
    for h,count in hist.items():P*=h**count
    factored=1
    for p,e in FACTORS.items():factored*=p**e
    need(P==factored,'factorized product of all exact minima')
    need(P**170 > M**18690,'exact candidate exceeds 2/85')
    lo=sum(F(n)*log_lo(h) for h,n in hist.items())/M-log_hi(M)/2
    hi=sum(F(n)*log_hi(h) for h,n in hist.items())/M-log_lo(M)/2
    need(lo>0,'positive candidate mean')
    clo=lo/log_hi(M);chi=hi/log_lo(M)
    need(F(2,85)<clo<=chi<F(3,125),'candidate interval')
    ticks=(clo*10**15).numerator//(clo*10**15).denominator
    lower_decimal=f'{ticks//10**15}.{ticks%10**15:015d}'
    upper_decimal=out_decimal(chi)
    need(F(lower_decimal)<=clo<=chi<=F(upper_decimal),
         'outward displayed decimal interval')
    return dict(status='PASS',M=M,states=M,denominators=ds,
        heights=values,histogram={str(k):v for k,v in hist.items()},
        product_factors={str(k):v for k,v in FACTORS.items()},
        integer_comparison='P^170 > 210^18690',integer_comparison_pass=True,
        mu_lower=str(lo),mu_upper=str(hi),c0_lower=str(clo),c0_upper=str(chi),
        c0_lower_decimal=lower_decimal,c0_upper_decimal=upper_decimal,
        method='All 210 states; every unit denominator <=105; exact integer minima.')

if __name__=='__main__':
    r=certificate()
    (ROOT/'results/candidate210.json').write_text(json.dumps(r,indent=2)+'\n')
    print('PASS candidate 210: all exact minima; P^170 > 210^18690; c0 interval',
          r['c0_lower_decimal'],r['c0_upper_decimal'])
