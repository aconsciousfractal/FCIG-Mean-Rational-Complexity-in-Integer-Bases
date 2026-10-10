#!/usr/bin/env python3
"""Fresh standard-library rational arithmetic for the105 component.

No executable helper is imported from the42/70/66 deliveries. Logarithms
use the positive -log(1-z) series after binary range reduction; table
construction follows the structural odd-prime proof rather than a search.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd

CHECKS=0

def require(condition,label):
    global CHECKS
    CHECKS+=1
    if not condition: raise AssertionError(label)

def record(x):
    x=F(x)
    return dict(numerator=x.numerator,denominator=x.denominator)

def fraction(x):
    return F(x['numerator'],x['denominator'])

@lru_cache(None)
def log_interval(x):
    x=F(x);require(x>0,'logarithm domain')
    exponent=0
    while x>=2:x/=2;exponent+=1
    while x<1:x*=2;exponent-=1
    def series(y):
        z=1-1/y;term=z;value=F(0)
        for n in range(1,129):
            value+=term/n;term*=z
        return value,value+term/(129*(1-z))
    lower,upper=series(x);low2,high2=series(F(2))
    if exponent>=0:return lower+exponent*low2,upper+exponent*high2
    return lower+exponent*high2,upper+exponent*low2

@lru_cache(None)
def phi_interval(x):
    x=F(x);require(x>0,'Phi domain')
    a,b=log_interval(1+x);c,d=log_interval(1+1/x)
    return a+x*c,b+x*d

def decimal_bound(x,up=True,digits=15):
    x=F(x);s=10**digits
    n=-((-x.numerator*s)//x.denominator) if up else (x.numerator*s)//x.denominator
    sign='-' if n<0 else '';n=abs(n)
    return f'{sign}{n//s}.{n%s:0{digits}d}'

def structural_table(q,radius):
    """Positive-half table from the odd three/four-prime argument."""
    require(radius in (3,5),'supported structural radius')
    table={}
    for z in range(q):
        if gcd(z,q)==1:rows=[(F(0),F(1,2),1,0)]
        else:
            usable=[i for i in range(1,4 if radius==3 else 5) if gcd(i*z-1,q)==1]
            if usable:rows=[(F(0),F(1,2),usable[0],-1)]
            else:
                require(radius==5,'three-prime candidate exists')
                rows=[(F(0),F(2,5),5,-1),(F(1,3),F(1,2),3,-2)]
        table[z]=rows
    return table

def periodic(q):
    require(q>1,'unit-sum symmetry excludes q1')
    residues=[j for j in range(1,q) if gcd(j,q)==1]
    n=len(residues);rho=F(n,q);A=F(5,2)*rho
    require(sum(residues)==F(q*n,2),'unit-sum identity')
    rows=[];count=total=0;C=F(0)
    for j in range(q):
        if j in residues:count+=1;total+=j
        error=count-rho*j;residual=2*j*count+total-A*j*j
        slope=3*error;endpoint=residual/j if j else F(0)
        C=max(C,slope,endpoint)
        rows.append(dict(j=j,U=count,S=total,slope=record(slope),endpoint=record(endpoint)))
    for row in rows:
        require(fraction(row['slope'])<=C,'all-block slope')
        require(fraction(row['endpoint'])<=C,'all-block endpoint')
    return rho,A,C,rows

def direct_minima(M,q):
    """Definition-level enumeration with the admissible d=1 upper bound."""
    answer=[]
    for s in range(M):
        a=s if 2*s<=M else s-M
        best=(max(1,abs(a)),1,a)
        d=2
        while d<=best[0]:
            if gcd(d,q)==1:
                a=s*d%M
                if 2*a>M:a-=M
                candidate=(max(d,abs(a)),d,a)
                if candidate<best:best=candidate
            d+=1
        H,d,a=best;answer.append((H,a,d))
    return answer
