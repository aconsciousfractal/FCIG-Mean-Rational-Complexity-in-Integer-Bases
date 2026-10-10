"""Exact rational interval arithmetic. No floating-point decisions."""
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
LS=10**15
CS=10**12

def need(value,message):
    if not value:raise ArithmeticError(message)
def ceildiv(a,b):
    need(b>0,'positive divisor');return -((-a)//b)
def upper_tick(v,scale=CS):
    v=F(v);return ceildiv(v.numerator*scale,v.denominator)
def lower_tick(v,scale=CS):
    v=F(v);return (v.numerator*scale)//v.denominator

def atanh_log(y,terms=28):
    y=F(y);need(1<=y<=2,'reduced logarithm domain')
    z=(y-1)/(y+1);z2=z*z;p=z;v=F(0)
    for j in range(terms):v+=2*p/(2*j+1);p*=z2
    return v,v+2*p/((2*terms+1)*(1-z2))
@lru_cache(None)
def log_interval(x):
    x=F(x);need(x>0,'positive logarithm')
    if x<1:
        a,b=log_interval(1/x);return -b,-a
    e=x.numerator.bit_length()-x.denominator.bit_length()
    if x<F(2)**e:e-=1
    y=x/F(2)**e;a,b=atanh_log(y);c,d=atanh_log(F(2))
    return a+e*c,b+e*d
@lru_cache(None)
def ticks(x):
    a,b=log_interval(F(x));return lower_tick(a,LS),upper_tick(b,LS)
def log_lo(x):return F(ticks(F(x))[0],LS)
def log_hi(x):return F(ticks(F(x))[1],LS)
def phi_hi(x):
    x=F(x);need(x>0,'positive entropy argument')
    return log_hi(1+x)+x*log_hi(1+1/x)
@lru_cache(None)
def alternate_log(x,terms=96):
    x=F(x);need(x>0,'positive logarithm')
    if x<1:
        a,b=alternate_log(1/x,terms);return -b,-a
    y=x;e=0
    while y>=2:y/=2;e+=1
    def series(t):
        z=1-1/t;p=z;v=F(0)
        for j in range(1,terms+1):v+=p/j;p*=z
        return v,v+p/((terms+1)*(1-z))
    a,b=series(y);c,d=series(F(2));return a+e*c,b+e*d

def alpha_interval(n=1000):
    v=sum((F(1,i*i) for i in range(1,n+1)),F(0))
    zlo=v+F(1,n+1)+F(1,2*(n+1)**2)
    zhi=v+F(2,2*n+1)
    need(zlo<zhi,'zeta trapezoid and midpoint enclosure')
    return 1/zhi,1/zlo

def out_decimal(v,digits=15):
    v=F(v);z=upper_tick(v,10**digits);sign='-' if z<0 else '';z=abs(z)
    return sign+str(z//10**digits)+'.'+str(z%10**digits).zfill(digits)

def sieve(n):
    a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if a[p]:a[p*p::p]=b'\0'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if a[i]]
def prime_trial(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
