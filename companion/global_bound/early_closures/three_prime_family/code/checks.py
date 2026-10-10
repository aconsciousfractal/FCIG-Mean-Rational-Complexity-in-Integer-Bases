"""Finite independent controls for the general proof. Tests are not a general theorem."""
from fractions import Fraction as F
from math import gcd,isqrt
from collections import Counter
from functools import lru_cache
from exact import need,log_lo,log_hi,analytic_upper,phi_hi,alternate_log
import json
from pathlib import Path

@lru_cache(None)
def phi(n):return sum(gcd(n,j)==1 for j in range(1,n+1))
def U(n):return n-n//2-n//3+n//6

def direct(M,q):
    """All denominators smaller than the current upper height; no lattice algorithm."""
    out=[]
    for s in range(M):
        a=s if 2*s<=M else s-M;d0=1;H=max(1,abs(a));d=2
        while d<H:
            if gcd(d,q)==1:
                b=s*d%M;b=b if 2*b<=M else b-M
                h=max(abs(b),d)
                if h<H:H=h;a=b;d0=d
            d+=1
        need(d0>0 and gcd(d0,q)==1 and (a-s*d0)%M==0,'direct witness')
        out.append((H,a,d0))
    return out

def factors_6p(M):
    if M%6:return None
    r=M
    for p in (2,3):
        while r%p==0:r//=p
    if r<5:return None
    p=next((d for d in range(5,isqrt(r)+1,2) if r%d==0),r)
    while r%p==0:r//=p
    return p if r==1 else None

def far_holes(M,T):
    need(T*T>2*M,'strict far-tail start requires T^2>2M')
    return sum((F(phi(d))*(1+F(2*(M-T*d),d*(T-2*d))) for d in range(1,(M-1)//T+1) if gcd(d,6)>1),F(0))

def kernel(eta,t):
    if eta==0:return (1,-1) if t>=0 else (1,1)
    if eta in (1,5):return (1,0)
    if eta==2:return (1,-1) if t>=0 else (2,1)
    if eta==4:return (1,1) if t<=0 else (2,-1)
    return (2,-1) if t>=0 else (2,1)

def expect_failure(fn,message):
    try:fn()
    except (ArithmeticError,ValueError) as e:return str(e)
    raise ArithmeticError('negative control accepted: '+message)

def run():
    count_period=0
    # Each infinite residue-class polynomial is verified by its coefficient formula.
    polynomial=[]
    for r in range(6):
        # Difference bound minus S_6(6n+r) is affine in n.
        val=[]
        for n in (0,1,2):
            A=6*n+r;u=2*n+int(r>=1)+int(r>=5)
            S=6*n*n+(6*n+1)*int(r>=1)+(6*n+5)*int(r>=5)
            need(U(A)==u,'unit count period')
            val.append(F(A*A,6)+F(2*A,3)+F(1,6)-S)
        need(val[0]>=0 and val[1]>=val[0] and val[2]-val[1]==val[1]-val[0],'unit sum polynomial')
        polynomial.append({'r':r,'constant':str(val[0]),'slope':str(val[1]-val[0])})
    S=0
    for A in range(1,5001):
        if gcd(A,6)==1:S+=A
        need(U(A)<=F(A,3)+F(2,3),'unit count upper')
        need(S<=F(A*A,6)+F(2*A,3)+F(1,6),'unit sum upper')
        need(2*A*U(A)+S<=F(5*A*A,6)+F(13*A,6),'new incidence polynomial')
        count_period+=1
    table=0
    for eta in range(6):
        for n in range(-256,257):
            t=F(n,512);i,j=kernel(eta,t)
            need(1<=i<=2 and gcd(j+i*eta,6)==1 and abs(j+i*t)<=1,'repair table')
            table+=1
    ms=[M for M in range(30,1201,6) if factors_6p(M)]
    ms+= [6*127,12*127,18*127,36*127,6*131,12*131,6*149,6*167,6*257,36*5*5,12*7*7]
    ms=sorted(set(ms));states=0;badstates=0;bad_cdf=0;far_tests=0;central_tests=0
    total_means=0;small_margin=None;missing_atom=None
    for M in ms:
        p=factors_6p(M);need(p is not None,'bank support')
        A=direct(M,6);B=direct(M,6*p);states+=M
        bad=Counter()
        for s,((H,a,d),(R,aa,dd)) in enumerate(zip(A,B)):
            need(H<=R,'relaxation direction')
            if d%p:
                need(R==H,'already p-unit')
            else:
                bad[H]+=1;badstates+=1
                need(a%p==0,'bad numerator divisibility')
                need(R*H<=H*H+2*M,'pointwise repair bound')
        cum=0
        for T in range(1,max(x[0] for x in A)+1):
            cum+=bad[T]
            need(6*p*cum<=5*T*T+13*p*T,'refined bad-state CDF')
            bad_cdf+=1
        root=isqrt(M)
        # Exact primitive images compared with the full-state height CDF.
        for T in sorted(set([1,max(1,root//2),isqrt((M-1)//2),root])):
            P6=0;Pq=0;pos={};neg={}
            for d in range(1,T+1):
                if gcd(d,6)>1:continue
                for a in range(1,T+1):
                    if gcd(a,d)>1:continue
                    P6+=1
                    if gcd(d,M)==1:
                        Pq+=1;v=a*pow(d,-1,M)%M
                        need(v not in pos,'same sign collision unexpected');pos[v]=(a,d)
                        v=(-v)%M;need(v not in neg,'negative sign collision unexpected');neg[v]=(a,d)
            need(Pq>=P6-T*U(T//p),'deleted primitive denominators')
            J=len(set(pos)&set(neg));C=1+len(set(pos)|set(neg))
            trueC=sum(h<=T for h,a,d in B)
            need(C==trueC and C==1+2*Pq-J,'exact central CDF')
            if 2*T*T<=M:need(J==0,'central collision-free boundary')
            else:
                L=2*F(T)-F(M,T)
                need(J<=L*L/18+8*L/9+F(16,9),'collision triangle')
            central_tests+=1
        for T in sorted(set([isqrt(2*M)+1,2*root+1,3*root+1,max(isqrt(2*M)+1,M//4)])):
            cap=far_holes(M,T);Strue=sum(h>T for h,a,d in A)
            need(Strue<=cap,'Farey survival bound for marked6 with extra prime')
            noatoms=cap-sum(phi(d) for d in range(1,(M-1)//T+1) if gcd(d,6)>1)
            if missing_atom is None and Strue>noatoms:missing_atom={'M':M,'T':T,'survival':Strue,'length_only_upper':str(noatoms)}
            far_tests+=1
        v=sum(log_hi(h) for h,a,d in B)/M-log_lo(M)/2
        hi=analytic_upper(M,p)
        need(v<hi,'analytic bound versus complete small-model witness upper')
        total_means+=1
    # Exact deleted-density identity with nonintegral endpoints; a tighter alpha interval only in this bank.
    z=sum(F(1,n*n) for n in range(1,101))
    alo=1/(z+F(2,201));ahi=1/(z+F(1,101))
    def mob(n):
        count=0;p=2
        while p*p<=n:
            if n%p==0:
                n//=p;count+=1
                if n%p==0:return 0
            p+=1
        if n>1:count+=1
        return -1 if count%2 else 1
    density=0
    for z in range(4,2001):
        Q=F(z,4);n=Q.numerator//Q.denominator
        a=sum((F(phi(d),d) for d in range(1,n+1) if gcd(d,6)>1),F(0))
        b=sum((F(mob(d),d)*(F(1,2)*(Q//(2*d))+F(2,3)*(Q//(3*d))-F(5,6)*(Q//(6*d))) for d in range(1,n+1) if gcd(d,6)==1),F(0))
        need(a==b,'deleted density floor identity')
        E=1+F(2,9)*log_lo(Q)
        need(a-alo*Q/2<=E and ahi*Q/2-a<=E,'density discrepancy enclosure')
        density+=1
    a=direct(30,6)[6][0];b=direct(30,30)[6][0]
    need((a,b)==(5,6),'relaxation inequivalence fixture')
    neg={
        'extra_prime_radical':expect_failure(lambda:analytic_upper(210,5),'wrong radical'),
        'composite_parameter':expect_failure(lambda:analytic_upper(6*169,169),'composite p'),
        'strict_far_boundary':expect_failure(lambda:far_holes(450,30),'far equality endpoint'),
        'false_model_equality':expect_failure(lambda:need(a==b,'R_30,6(6)=5 differs from R_30,30(6)=6'),'relaxation equality'),
    }
    return {'status':'PASS','unit_prefix_checks':count_period,'unit_sum_residue_polynomials':polynomial,'repair_table_checks':table,'complete_moduli':ms,'state_modulus_incidents':states,'bad_states_repaired':badstates,'bad_CDF_checks':bad_cdf,'central_CDF_checks':central_tests,'marked6_Farey_checks':far_tests,'complete_mean_checks':total_means,'density_real_cutoff_checks':density,'missing_integer_atom_fixture':missing_atom,'negative_controls':neg,'finite_banks_are_not_general_proof':True}

if __name__=='__main__':
    import sys
    obj=run();Path(sys.argv[1]).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    print(json.dumps(obj,indent=2),flush=True)
