"""Finite physical banks against denominator enumeration, not Gauss."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import arithmetic as ar
import excision

def run(p,prefix):
    q=6*p;R,_=excision.direct(q,q);sq=None
    if p in (5,7,11,13):sq=excision.direct(q*q,q)[0]
    roots=0;point=0;conditional=0
    for line in Path(str(prefix)+'.points').read_text().splitlines():
        k,bound=map(int,line.split())
        vals=[1] if k==0 else R if k==1 else sq if k==2 else None
        if vals is None:continue
        M=q**k;upper=sum(ar.logs(x)[1] for x in vals)/F(M*ar.LS)-F(ar.logs(M)[0],2*ar.LS)
        ar.need(upper<=F(bound,ar.S),'complete root bank versus direct denominators');roots+=len(vals)
    for line in Path(str(prefix)+'.bank').read_text().splitlines():
        k,s,bound=map(int,line.split());ar.need(0<=s<q,'bank source state')
        if k==0:
            upper=F(ar.logs(R[s])[1],ar.LS)-F(ar.logs(q)[0],2*ar.LS);point+=1
        elif k==1 and sq is not None:
            upper=sum(ar.logs(sq[s+q*b])[1] for b in range(q))/F(q*ar.LS)-F(ar.logs(q*q)[0],2*ar.LS);conditional+=1
        else:continue
        ar.need(upper<=F(bound,ar.S),'conditional marked bank versus all direct descendants')
    return {'p':p,'q':q,'pointwise_F_banks':point,'conditional_U_banks':conditional,'total_root_states':roots,'maximum_direct_modulus':q*q if sq is not None else q,'minima_method':'every admissible denominator up to half the modulus','cubic_moduli_not_enumerated':True}
