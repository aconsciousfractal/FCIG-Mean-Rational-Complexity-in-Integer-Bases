"""Rational logarithms and outward integer ticks. No floating arithmetic."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
LS = 10**15

def need(value, message):
    if not value:
        raise ArithmeticError(message)

def ceildiv(a, b):
    need(b > 0, 'positive divisor')
    return -((-a)//b)

def atanh_log(y, terms=28):
    y = F(y)
    need(1 <= y <= 2, 'range-reduced logarithm')
    z=(y-1)/(y+1);z2=z*z
    power=z; total=F(0)
    for j in range(terms):
        total += 2*power/(2*j+1)
        power *= z2
    err=2*power/((2*terms+1)*(1-z2))
    return total,total+err

@lru_cache(None)
def log_interval(x, terms=28):
    x=F(x)
    need(x>0,'log domain')
    if x<1:
        lo,hi=log_interval(1/x,terms)
        return -hi,-lo
    e=x.numerator.bit_length()-x.denominator.bit_length()
    if x<F(2)**e:e-=1
    y=x/F(2)**e
    need(1<=y<2,'range reduction exact')
    lo,hi=atanh_log(y,terms)
    a,b=atanh_log(F(2),terms)
    return lo+e*a,hi+e*b

@lru_cache(None)
def ticks(x):
    lo,hi=log_interval(F(x))
    return lo.numerator*LS//lo.denominator,ceildiv(hi.numerator*LS,hi.denominator)

def log_lo(x):return F(ticks(F(x))[0],LS)
def log_hi(x):return F(ticks(F(x))[1],LS)

def alternate_log(x, terms=96):
    """Independent series: -log(1-z)=sum z^j/j, scale reduction by 2."""
    x=F(x)
    if x<1:
        a,b=alternate_log(1/x,terms);return -b,-a
    y=x;e=0
    while y>=2:y/=2;e+=1
    def series(t):
        z=1-1/t;power=z;v=F(0)
        for j in range(1,terms+1):v+=power/j;power*=z
        return v,v+power/((terms+1)*(1-z))
    a,b=series(y);c,d=series(F(2))
    return a+e*c,b+e*d

def phi_hi(p):
    need(p>=5,'prime parameter lower range')
    x=F(5,3*p)
    return log_hi(1+x)+x*log_hi(1+1/x)

def analytic_upper(M,p):
    """Theorem A: replace sqrt(M) by its lower integer, in decreasing errors."""
    need(M>=6*p and p>=5,'analytic theorem domain')
    need(all(p%d for d in range(2,isqrt(p)+1)), 'prime parameter required')
    need(M%2==0 and M%3==0 and M%p==0, 'all three prime factors required')
    r=M
    for t in (2,3,p):
        while r%t==0:r//=t
    need(r==1,'radical must be exactly 6p')
    a=isqrt(M); l=log_hi(a)
    return phi_hi(p)-F(1,200)+F(17,30*p)+(20+F(31,20)*l)/a+(F(13,4)+F(4,9)*l)/(a*a)

def out_decimal(q, digits=15):
    q=F(q);v=ceildiv(q.numerator*10**digits,q.denominator)
    sign='-' if v<0 else ''; v=abs(v)
    return sign+str(v//10**digits)+'.'+str(v%10**digits).zfill(digits)
