"""Exact prime-repair certificate for mark 30. General proofs are supplied in the main manuscript."""
from fractions import Fraction as F
from math import gcd, isqrt, prod
from functools import lru_cache
from exact import need, log_lo, log_hi, phi_hi, out_decimal, sieve, prime_trial
BREAKS=(F(-1,2),F(-1,3),F(-1,4),F(0),F(1,4),F(1,3),F(1,2))

def repair_table():
    rows=[]
    for eta in range(30):
        for lo,hi in zip(BREAKS,BREAKS[1:]):
            choices=[(i,j) for i in range(1,5) for j in range(-4,5)
                     if gcd(j+i*eta,30)==1 and max(abs(j+i*lo),abs(j+i*hi))<=1]
            need(choices,'uncovered real repair interval')
            i,j=min(choices)
            rows.append(dict(eta=eta,lo=str(lo),hi=str(hi),i=i,j=j))
    bad=[]
    for i in range(1,4):
        for j in range(-4,5):
            if abs(F(j)+F(i,4))<=1:bad.append([i,j,gcd(j+3*i,30)])
    need(all(z[2]>1 for z in bad),'four-coefficient necessity within specified class')
    need(gcd(-1+4*3,30)==1 and -1+4*F(1,4)==0,'four-coefficient repair')
    return dict(status='PASS',closed_real_intervals=180,rows=rows,
                bounded_negative=dict(eta=3,tau='1/4',coefficient_range=[1,3],admissible_longitudinal_pairs=bad,
                                      conclusion='no pair with i<=3 and |j+i*tau|<=1; not a lower bound on the true minimum'))

TABLE=repair_table()['rows']
def choose(eta,tau):
    need(0<=eta<30 and -F(1,2)<=tau<=F(1,2),'repair parameter domain')
    z=next(r for r in TABLE if r['eta']==eta and F(r['lo'])<=tau<=F(r['hi']))
    return z['i'],z['j']

def periodic():
    u=F(4,15);U=S=0;rows=[];b=F(0)
    for j in range(30):
        if j and gcd(j,30)==1:U+=1;S+=j
        eps=U-u*j
        a=3*eps
        c=F(2*j*U+S,j)-F(5,2)*u*j if j else F(0)
        b=max(b,a,c)
        rows.append(dict(j=j,U=U,S=S,eps=str(eps),linear_coefficient=str(a),constant_coefficient=str(c)))
    need(b==F(166,57),'exact periodic linear coefficient')
    # Coefficientwise comparison proves all n>=0 for A=30n+j, not an extrapolation.
    for r in rows:
        j=r['j'];eps=F(r['eps']);D=F(2*j*r['U']+r['S'])-F(2,3)*j*j
        need(90*eps<=30*b and D<=b*j,'all-period polynomial coefficients')
    U=S=checks=0
    for A in range(1,3001):
        if gcd(A,30)==1:U+=1;S+=A
        n,j=divmod(A,30);r=rows[j]
        pu=8*n+r['U'];ps=120*n*n+30*n*r['U']+r['S']
        need((U,S)==(pu,ps),'independent periodic count and sum')
        need(2*A*U+S<=F(2,3)*A*A+b*A,'periodic bound bank')
        checks+=1
    return dict(status='PASS',density='4/15',quadratic_coefficient='2/3',linear_coefficient=str(b),
                coefficient_tests=60,integer_bank_max=3000,integer_bank_tests=checks,rows=rows)

@lru_cache(None)
def direct(M,q):
    """Definition-based exact minima; d=1 supplies a finite initial bound."""
    need(M>=1 and q>=1,'positive model')
    R=[];W=[];tests=0
    for s in range(M):
        a=s if 2*s<=M else s-M;best=max(1,abs(a));ba,bd=a,1
        d=2
        while d<best:
            if gcd(d,q)==1:
                x=s*d%M;x=x if 2*x<=M else x-M;h=max(abs(x),d);tests+=1
                if h<best:best=h;ba,bd=x,d
            d+=1
        R.append(best);W.append((ba,bd))
    return R,W,tests

def egcd(a,b):
    oa,ob=a,b;x0,x1,y0,y1=1,0,0,1
    while b:
        k=a//b;a,b=b,a-k*b;x0,x1=x1,x0-k*x1;y0,y1=y1,y0-k*y1
    if a<0:a,x0,y0=-a,-x0,-y0
    need(x0*oa+y0*ob==a,'Bezout identity')
    return a,x0,y0

def repair(M,s,p,a,d):
    need(prime_trial(p) and p>=7 and M%p==0,'additional prime domain')
    need(d>0 and gcd(d,30)==1 and (a-s*d)%M==0 and d%p==0,'damaged witness')
    n=(a-s*d)//M;g,x,y=egcd(n,d);need(g==1,'primitive in lattice coordinates')
    b,e=-y*M+x*s,x
    need(a*e-b*d==M and gcd(d,e)==1,'oriented basis completion')
    norm=a*a+d*d;t=(2*(a*b+d*e)+norm)//(2*norm);b-=t*a;e-=t*d
    tau=F(a*b+d*e,norm);eta=e*pow(d,-1,30)%30;i,j=choose(eta,tau)
    aa=j*a+i*b;dd=j*d+i*e
    if dd<0:aa,dd=-aa,-dd
    need(dd>0 and gcd(dd,30*p)==1 and (aa-s*dd)%M==0,'repaired full mark')
    H=max(abs(a),d);K=max(abs(aa),dd)
    need(K*H<=H*H+4*M,'pointwise repair budget')
    return aa,dd,K

def physical_checks():
    ps=[7,11,13,17,19,23,29,31,37,41,43]
    domain=sorted((p*m,p) for p in ps for m in (1,2,3,5,6,10,30))
    count=badn=density_checks=0;rows=[];b=F(166,57)
    for M,p in domain:
        R,W,nt=direct(M,30);T,_,nf=direct(M,30*p)
        bad={}
        for s,((a,d),h,k) in enumerate(zip(W,R,T)):
            need(h<=k,'mark relaxation ordering')
            if d%p:
                need(h==k,'undamaged exact witness')
            else:
                _,_,v=repair(M,s,p,a,d);need(k<=v,'repaired witness upper, not presumed minimum');badn+=1
                bad[h]=bad.get(h,0)+1
            count+=1
        acc=0
        for H in sorted(bad):
            acc+=bad[H]
            need(acc<=min(F(M),F(2*H*H,3*p)+b*H),'all damaged-height jumps')
            density_checks+=1
        diff=sum(log_hi(y)-log_lo(x) for x,y in zip(R,T))/M
        cap=phi_hi(F(8,3*p))+F(7304,399*isqrt(M))
        need(diff<=cap,'mean repair finite regression')
        rows.append(dict(M=M,p=p,states=M,damaged=sum(bad.values()),upper_excision_mean=str(diff),upper_formula=str(cap),denominators_tested=nt+nf))
    R,W,_=direct(330,30);T,_,_=direct(330,330)
    need(R[30]==11 and T[30]==30 and W[30]==(0,11),'relaxed/full separation fixture')
    need((0-30*1)%330!=0,'illegal division by ambient gcd detected')
    return dict(status='PASS',pairs=len(domain),states=count,damaged_witnesses_repaired=badn,
                density_jump_checks=density_checks,domain=rows,
                separation_fixture=dict(M=330,state=30,relaxed_height=11,full_height=30,relaxed_witness=[0,11]))

def candidate():
    R,W,nt=direct(210,210);P=prod(R);need(P**170>210**18690,'candidate exact integer gate')
    lower=0
    for s,h in enumerate(R):
        for d in range(1,h):
            if gcd(d,210)==1:
                a=s*d%210;need(max(d,min(a,210-a))>=h,'candidate minimality');lower+=1
    factors={};z=P
    for p in sieve(max(R)):
        e=0
        while z%p==0:e+=1;z//=p
        if e:factors[str(p)]=e
    need(z==1,'complete candidate product')
    return dict(status='PASS',M=210,minima=R,witnesses=W,product=str(P),factorization=factors,
                lower_denominator_checks=lower,integer_comparison='P^170 > 210^18690',candidate_lower='2/85')

def simple_table():
    """An independent two-half-interval implementation of the same repair lemma."""
    rows=[];groups={}
    for eta in range(30):
        if gcd(eta,30)==1:i,j=1,0
        else:
            i=next(i for i in range(1,5) if gcd(i*eta-1,30)==1);j=-1
        need(gcd(j+i*eta,30)==1 and max(abs(F(j)),abs(F(j)+F(i,2)))<=1,'positive half-interval')
        groups.setdefault(f'{i},{j}',[]).append(eta)
        rows.append(dict(eta=eta,tau_lo='0',tau_hi='1/2',i=i,j=j))
        ee=(-eta)%30
        need(gcd(-j+i*ee,30)==1 and max(abs(F(-j)),abs(F(-j)-F(i,2)))<=1,'negative half-interval')
        rows.append(dict(eta=ee,tau_lo='-1/2',tau_hi='0',i=i,j=-j))
    return dict(status='PASS',closed_half_intervals=60,positive_half_groups=groups,rows=rows)
