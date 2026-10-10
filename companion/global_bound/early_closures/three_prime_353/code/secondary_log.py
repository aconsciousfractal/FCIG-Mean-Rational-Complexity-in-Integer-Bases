"""Independent log enclosures from -log(1-z); outward integer arithmetic.\nNo calls to the primary atanh-log implementation.\n"""
from functools import lru_cache
from exact import need
PREC=10**45
TERMS=96

def ceildiv(a,b):return (a+b-1)//b

def series(num,den):
    # 0<=z=num/den<=1/2. Invariant: pl/PREC <= z^j <= pu/PREC.
    need(den>0 and 0<=2*num<=den,'secondary series domain')
    pl=PREC*num//den;pu=ceildiv(PREC*num,den);sl=su=0
    for j in range(1,TERMS+1):
        sl+=pl//j;su+=ceildiv(pu,j)
        pl=pl*num//den;pu=ceildiv(pu*num,den)
    # tail <= z^(TERMS+1)/[(TERMS+1)*(1-z)].
    su+=ceildiv(pu*den,(TERMS+1)*(den-num))
    return sl,su

@lru_cache(None)
def log2():return series(1,2)

def enclosure(x):
    need(isinstance(x,int) and x>=1,'positive integer secondary log input')
    e=x.bit_length()-1;base=1<<e
    a,b=series(x-base,x);c,d=log2()
    return a+e*c,b+e*d

def contained(x,lo,hi,primary_scale):
    a,b=enclosure(x)
    return lo*PREC<=a*primary_scale and b*primary_scale<=hi*PREC
