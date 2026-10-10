"""Independent integer/rational checks; no `assert` is used.
The large membership pass streams the C++ witness dump and verifies every pair.
The small end-to-end pass checks all 44,100 actual states against the definition.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import gcd,prod
from pathlib import Path
from struct import Struct,iter_unpack
from itertools import product
import mmap,hashlib,json

Q=210; PRIMES=(2,3,5,7); S=10**9

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
    out=[None]*576
    for a in range(1,211):
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

def algebra_checks():
    transfer_cases=0
    for lab in LABELS:
        r=rep(lab);n=0
        for k,(H,rh) in enumerate(BR):
            t=mv(H,r);ok=gcd(gcd(t[0],t[1]),Q)==1
            equivalent=all(a!=b for a,b in zip(lab,LABELS[k]))
            need(ok==equivalent,'kernel-filter equivalence')
            if ok:
                need(mark(t)==mark(rh),'image reset');n+=1;transfer_cases+=1
        need(n==210,'probability normalisation')
    for test in range(7):
        values=[((k+test+1)**3%101)-50 for k in range(576)]
        x=values[:];stride=1
        for p in reversed(PRIMES):
            n=p+1
            for start in range(0,576,stride*n):
                for off in range(stride):
                    ids=[start+off+t*stride for t in range(n)]
                    total=sum(x[i] for i in ids)
                    for i in ids:x[i]=total-x[i]
            stride*=n
        for j,lab in enumerate(LABELS):
            slow=sum(v for k,v in enumerate(values) if all(a!=b for a,b in zip(lab,LABELS[k])))
            need(slow==x[j],'tensor/direct sum mismatch')
    return {'rank_one_branches':576,'projective_input_labels':576,
            'admitted_label_branch_pairs':transfer_cases,'tensor_signals':7,'tensor_outputs_checked':7*576}


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

def egcd(a,b):
    x0,x1,y0,y1=1,0,0,1
    while b:q,a,b=a//b,b,a%b;x0,x1=x1,x0-q*x1;y0,y1=y1,y0-q*y1
    if a<0:a,x0,y0=-a,-x0,-y0
    return a,x0,y0

def verify_dump(path:Path, part, direct_R, table):
    """Every stored witness is independently checked, using arbitrary-size Python integers."""
    head=Struct('<iiiqqqq');pair=Struct('<qq');record_size=head.size+576*pair.size
    keys={};checked=0;maxw=0
    with path.open('rb') as f:
        n=int.from_bytes(f.read(8),'little',signed=True)
        for row in range(n):
            raw=f.read(record_size);need(len(raw)==record_size,'truncated witness data')
            a=head.unpack_from(raw);key=tuple(a[:3]);B=tuple(a[3:]);need(det(B)>0,'centre determinant')
            need(key not in keys,'duplicate witness centre');keys[key]=8+row*record_size+head.size
            for k,(wx,wy) in enumerate(iter_unpack('<qq',raw[head.size:])):
                H,rh=BR[k];a,_,b,d=H
                nx=wx*d-wy*b;ny=wy*a
                need(nx%Q==0 and ny%Q==0,'witness membership')
                need(gcd((nx//Q)*rh[0]+(ny//Q)*rh[1],Q)==1,'witness primitive mark')
                maxw=max(maxw,abs(wx),abs(wy));checked+=1
        need(not f.read(1),'trailing witness data')
    # Audit every q^2 physical state, with actual denominators from an independent program.
    M=Q*Q;seen=bytearray(M);maxgap=0;sum_log_upper=0
    with path.open('rb') as f:
        mmf=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)
        for s in range(Q):
            HB=(Q,0,s,1);Y,V,P=canonical(HB);r=mv(V,(0,1));lab=mark(r);key=partition_key(Y,part)
            need(key in keys,'missing bank centre')
            for k,(H,rh) in enumerate(BR):
                if any(a==b for a,b in zip(lab,LABELS[k])):continue
                wx,wy=pair.unpack_from(mmf,keys[key]+k*16)
                Z=mm(mm(H,V),HB);g,x,y=egcd(Z[1],Z[3]);need(g==1,'physical denominator projection')
                state=(x*Z[0]+y*Z[2])%M
                need(not seen[state],'duplicate state leaf');seen[state]=1
                W=mm(V,HB);aa=wx*W[0]+wy*W[2];dd=wx*W[1]+wy*W[3]
                need(gcd(dd,M)==1 and (aa-state*dd)%M==0,'end-to-end witness congruence')
                R=max(abs(aa),abs(dd));need(R>=direct_R[state],'upper witness below true minimum')
                maxgap=max(maxgap,R-direct_R[state])
                # The selected height need NOT equal the true minimum.
        mmf.close()
    need(sum(seen)==M,'physical bank mass')
    return {'witness_centres':n,'integer_witnesses_checked':checked,'maximum_coefficient':maxw,
            'complete_physical_bank_modulus':M,'complete_physical_bank_states':sum(seen),
            'maximum_selected_height_minus_true_height':maxgap}


def direct_result(path:Path,M:int):
    raw=path.read_bytes();need(len(raw)==4*M,'direct-definition size')
    R=[x[0] for x in iter_unpack('<i',raw)]
    need(R[0]==1 and min(R)>=1,'zero and positivity')
    need(all(R[s]==R[-s%M] for s in range(M)),'sign symmetry')
    return R

def candidate_check(R):
    need(len(R)==210,'candidate mass')
    P=prod(R);need(P**170>210**18690,'candidate inequality')
    lo,hi=rational_log(P);lb,ub=rational_log(210)
    c_lo=lo/(210*ub)-F(1,2);c_hi=hi/(210*lb)-F(1,2)
    decimal_scale=10**24
    c_lo=F(c_lo.numerator*decimal_scale//c_lo.denominator,decimal_scale)
    c_hi=F(-((-c_hi.numerator*decimal_scale)//c_hi.denominator),decimal_scale)
    return {'mass':210,'sum_heights':sum(R),'product_sha256':hashlib.sha256(str(P).encode()).hexdigest(),
            'P_power_170_gt_210_power_18690':True,'c0_lower':str(c_lo),'c0_upper':str(c_hi)}
