"""Independent integer/rational checks; no `assert` is used.
Squarefree-mark HNF branches and rational logarithm series, independent of C++ production.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import gcd,prod
from pathlib import Path
from struct import Struct,iter_unpack
from itertools import product
import mmap,hashlib,json

Q=30; PRIMES=(2,3,5); S=10**9

def need(condition:bool, message:str)->None:
    if not condition: raise ArithmeticError(message)

def det(B):return B[0]*B[3]-B[1]*B[2]
def mm(A,B):
    return (A[0]*B[0]+A[1]*B[2],A[0]*B[1]+A[1]*B[3],
            A[2]*B[0]+A[3]*B[2],A[2]*B[1]+A[3]*B[3])
def mv(A,r):return (A[0]*r[0]+A[1]*r[1],A[2]*r[0]+A[3]*r[1])
def mark(r):
    return tuple((b*pow(a,-1,p))%p if a else p
                 for p in PRIMES for a,b in [(r[0]%p,r[1]%p)])
def idx(l):
    z=0
    for t,p in zip(l,PRIMES):z=z*(p+1)+t
    return z
def crt(vals):return sum(v*(Q//p)*pow(Q//p,-1,p) for p,v in zip(PRIMES,vals))%Q
def rep(l):return (crt([int(t<p) for t,p in zip(l,PRIMES)]),
                   crt([t if t<p else 1 for t,p in zip(l,PRIMES)]))
LABELS=list(product(*(range(p+1) for p in PRIMES)))

def branches():
    out=[None]*len(LABELS)
    for a in range(1,Q+1):
        if Q%a:continue
        d=Q//a
        for b in range(a):
            H=(a,0,b,d); ker=[]; image=[]
            for p in PRIMES:
                lines=[(1,t) for t in range(p)]+[(0,1)]
                killed=[t for t,r in enumerate(lines) if all(x%p==0 for x in mv(H,r))]
                need(len(killed)==1,'rank-one kernel')
                ker.append(killed[0])
                nonzero=next(mv(H,r) for r in lines if any(x%p for x in mv(H,r)))
                image.append(tuple(x%p for x in nonzero))
            rh=(crt([v[0] for v in image]),crt([v[1] for v in image]))
            j=idx(ker);need(out[j] is None,'branch duplicate')
            out[j]=(H,rh)
    need(all(x is not None for x in out),'branch coverage')
    return out

BR=branches()

def canonical(B):
    """Return Y,V,P with Y=V*B*P. V integral unimodular; P signed permutation."""
    V=(1,0,0,1);P=(1,0,0,1);Y=tuple(B)
    while True:
        u,v=Y[:2],Y[2:];nu=sum(x*x for x in u);nv=sum(x*x for x in v)
        if nu>nv:
            A=(0,1,1,0);Y=mm(A,Y);V=mm(A,V);u,v=Y[:2],Y[2:];nu=nv
        dot=sum(x*y for x,y in zip(u,v));t=(2*dot+nu)//(2*nu)
        if not t:break
        A=(1,0,-t,1);Y=mm(A,Y);V=mm(A,V)
    for A in [( -1,0,0,1) if Y[0]<0 else (1,0,0,1)]:Y=mm(Y,A);P=mm(P,A)
    for A in [(1,0,0,-1) if Y[1]<0 else (1,0,0,1)]:Y=mm(Y,A);P=mm(P,A)
    if Y[0]<Y[1]:A=(0,1,1,0);Y=mm(Y,A);P=mm(P,A)
    if det(Y)<0:A=(1,0,0,-1);Y=mm(A,Y);V=mm(A,V)
    need(Y==mm(mm(V,B),P),'basis transformation')
    need(abs(det(V))==1 and abs(det(P))==1,'unimodular transformations')
    need(Y[0]>=Y[1]>=0 and det(Y)>0,'canonical orientation')
    return Y,V,P

def partition_key(Y,part):
    den,rb,ab=part;D=det(Y);nu=Y[0]**2+Y[1]**2
    kl=next(i-160 for i,t in enumerate(rb) if nu*den<=t*D)
    ka=next(i for i,t in enumerate(ab) if Y[1]*den<=Y[0]*t)
    dot=Y[0]*Y[2]+Y[1]*Y[3];kt=(40*dot+nu)//(2*nu)
    return (kl,kt,ka)

def load_partition(path):
    data=list(map(int,Path(path).read_text().split()));need(len(data)==186,'partition data count')
    return data[0],data[1:165],data[165:]

def setup(p, four=True):
    global Q,PRIMES,LABELS,BR
    need(p>=5 and all(p%d for d in range(2,__import__('math').isqrt(p)+1)), 'prime parameter')
    need((four and 11<=p<=31) or (not four and p==5), 'selected parameter domain')
    Q=(30 if four else 6)*p
    PRIMES=(2,3,5,p) if four else (2,3,p)
    LABELS=list(product(*(range(x+1) for x in PRIMES)))
    BR=branches()

def algebra_checks():
    local=0;resets=0;outputs=0;allowed_pairs=0
    for lab in LABELS:need(mark(rep(lab))==lab,'CRT label inverse')
    for j,(H,rh) in enumerate(BR):
        for ell in PRIMES:
            im=mark(rh)[PRIMES.index(ell)]
            for t,r in enumerate([(1,z) for z in range(ell)]+[(0,1)]):
                v=mv(H,r);a,b=v[0]%ell,v[1]%ell
                killed=(a==0 and b==0)
                need(killed==(t==LABELS[j][PRIMES.index(ell)]),'local kernel filter')
                if not killed:need((b*pow(a,-1,ell)%ell if a else ell)==im,'local image reset');resets+=1
                local+=1
    allowed=[]
    for lab in LABELS:
        row=[idx(t) for t in product(*(tuple(x for x in range(p+1) if x!=l) for p,l in zip(PRIMES,lab)))]
        need(len(row)==Q and len(set(row))==Q,'all admissible children and their weights')
        allowed.append(row);allowed_pairs+=len(row)
    for test in range(7):
        values=[((k+test+1)**3%101)-50 for k in range(len(LABELS))]
        x=values[:];stride=1
        for p in reversed(PRIMES):
            n=p+1
            for start in range(0,len(LABELS),stride*n):
                for off in range(stride):
                    ids=[start+off+t*stride for t in range(n)];total=sum(x[i] for i in ids)
                    for i in ids:x[i]=total-x[i]
            stride*=n
        for j,ids in enumerate(allowed):need(x[j]==sum(values[i] for i in ids),'signed tensor/direct sum');outputs+=1
    return dict(local_prime_line_checks=local,nonzero_image_checks=resets,labels=len(LABELS),branches=len(BR),admissible_pairs=allowed_pairs,signed_tensor_outputs=outputs)


def rational_log(n:int,d:int=1,terms:int=24):
    """Full rational atanh sum, independent of the original fixed-point generator."""
    need(n>0 and d>0,'positive logarithm argument')
    k=n.bit_length()-d.bit_length()
    a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
    if a<b:k-=1;a,b=(n,d<<k) if k>=0 else (n<<(-k),d)
    def series(a,b):
        z=F(a-b,a+b);v=z;ss=F(0)
        for j in range(terms):ss+=2*v/(2*j+1);v*=z*z
        return ss,ss+2*v/((2*terms+1)*(1-z*z))
    lo,hi=series(a,b);l2,u2=series(2,1)
    return (lo+k*l2,hi+k*u2) if k>=0 else (lo+k*u2,hi+k*l2)

def log_checks(path):
    tab=list(map(int,Path(path).read_text().split()))
    need(len(tab)==12545,'table length')
    for i,u in enumerate(tab):
        if i<=1024:need(u==0,'unused compact prefix');continue
        lo,hi=rational_log(i,1024)
        need(hi*S/2<=u,'upper log table')
        if i<=2048:need(lo*S/2>=u-1,'lower mantissa table')
    lo,hi=rational_log(210);need(lo>F(16,3),'log 210 lower bound')
    return tab,{'upper_log_entries':11520,'two_sided_mantissa_entries':1024,'log210_gt_16_over_3':True}
