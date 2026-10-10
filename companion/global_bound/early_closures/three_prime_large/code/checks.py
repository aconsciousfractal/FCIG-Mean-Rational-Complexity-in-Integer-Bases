"""Independent finite checks. None of these ranges is a substitute for a theorem."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd,isqrt
from functools import lru_cache
from collections import Counter
import csv,json
from exact import *
from family import parameters,bound,DELTA

@lru_cache(None)
def direct(M,q):
    out=[]
    for s in range(M):
        a=s if 2*s<=M else s-M;d0=1;H=max(1,abs(a));d=2
        while d<H:
            if gcd(d,q)==1:
                a2=s*d%M;a2=a2 if 2*a2<=M else a2-M
                h=max(abs(a2),d)
                if h<H:H=h;a=a2;d0=d
            d+=1
        need(d0>0 and gcd(d0,q)==1 and (a-s*d0)%M==0,'direct witness')
        out.append((H,a,d0))
    return out

def egcd(a,b):
    if b==0:return (abs(a),1 if a>=0 else -1,0)
    g,x,y=egcd(b,a%b);return g,y,x-(a//b)*y

def repair(M,s,a,d,q,p):
    need(p>=3 and prime_trial(p) and q%p!=0,'odd additional prime')
    f=[];z=q
    for r in range(2,isqrt(q)+1):
        if z%r==0:
            f.append(r)
            while z%r==0:z//=r
    if z>1:f.append(z)
    need(len(f)==2,'repair requires exactly two mark primes')
    need(d>0 and gcd(d,q)==1 and d%p==0 and (a-s*d)%M==0,'problematic witness input')
    k=(a-s*d)//M;g,x,y=egcd(k,d);need(g==1,'lattice primitivity')
    b=-M*y+s*x;e=x;nu=a*a+d*d;dot=a*b+d*e
    shift=(2*dot+nu)//(2*nu);b-=shift*a;e-=shift*d;t=F(a*b+d*e,nu)
    need(-F(1,2)<=t<=F(1,2),'reduced shear')
    eta=e*pow(d,-1,q)%q;sg=1 if t>=0 else -1
    if gcd(eta,q)==1:i,j=1,0
    elif gcd(eta-sg,q)==1:i,j=1,-sg
    else:i,j=2,-sg
    need(gcd(j+i*eta,q)==1 and abs(j+i*t)<=1,'general repair coefficients')
    aa=j*a+i*b;dd=j*d+i*e
    if dd<0:aa=-aa;dd=-dd
    H=max(abs(a),d);R=max(abs(aa),dd)
    need(dd>0 and gcd(dd,q*p)==1 and (aa-s*dd)%M==0,'repaired witness full mark')
    need(R*H<=H*H+2*M,'repaired height bound')
    return i,j,t,aa,dd

def euler(n):return sum(gcd(n,i)==1 for i in range(1,n+1))
def mob(n):
    sign=1;p=2
    while p*p<=n:
        if n%p==0:
            n//=p;sign=-sign
            if n%p==0:return 0
        p+=1
    if n>1:sign=-sign
    return sign

def holes(M,T,q):
    need(T*T>2*M,'strict Farey tail T^2>2M')
    return sum((F(euler(d))*(1+F(2*(M-T*d),d*(T-2*d))) for d in range(1,(M-1)//T+1) if gcd(d,q)>1),F(0))

def expect(fn):
    try:fn()
    except (ArithmeticError,ValueError) as e:return str(e)
    raise ArithmeticError('a deliberately invalid input was accepted')

def parameter_checks(rows):
    residues=0;polynomial=0;harmonic=0
    for row in rows:
        r=int(row['r']);s=int(row['s']);q=r*s;u=F(row['u']);dp=F(row['Dplus']);dm=F(row['Dminus']);K=F(row['primitive_floor_K']);span=F(row['interval_span']);bb=F(row['bad_incidence_b']);dd=F(row['deleted_discrepancy']);c=F(row['deleted_c']);h=F(row['harmonic_h'])
        A=F(r-1,r);B=F(s-1,s);C=1-F(1,q)
        for j in range(q):
            U=j-j//r-j//s+j//q
            S=j*(j+1)//2-r*(j//r)*(j//r+1)//2-s*(j//s)*(j//s+1)//2+q*(j//q)*(j//q+1)//2
            delta=U-u*j
            need(delta<=dp and u-delta<=dm and delta<=K and 2*u-delta<=K,'whole-residue floor bound')
            need(3*delta<=bb,'bad-count infinite slope')
            need(2*j*U+S-F(5,2)*u*j*j<=bb*j,'bad-count residue intercept')
            z=A*(j//r)+B*(j//s)-C*(j//q)-c*j
            need(-dd<=z-c<=z<=dd,'whole-residue deleted discrepancy')
            # Algebraic quadratic block law, two nonconstant coefficients.
            ph=(r-1)*(s-1)
            need(F(5,2)*q*ph==F(5,2)*u*q*q,'quadratic leading coefficient')
            need(3*q*U+2*j*ph-5*ph*j==3*q*delta,'quadratic linear coefficient')
            residues+=1;polynomial+=2
        H=F(0);U=0
        for j in range(1,17):
            if gcd(j,q)==1:H+=F(1,j);U+=1
            lo,hi=alternate_log(j)
            need(H-u*lo<=h,'independent harmonic prefix bound');harmonic+=1
        lo,hi=alternate_log(16)
        need(H-u*lo+(16*u-U+dp)/16<=h,'harmonic continuation beyond 16');harmonic+=1
        need(span==dp+dm,'closed interval discrepancy')
    return {'pair_rows':len(rows),'complete_residue_classes':residues,'polynomial_coefficient_checks':polynomial,'independent_harmonic_checks':harmonic}

def run(results):
    results=Path(results)
    rows=list(csv.DictReader((results/'pair_profiles.csv').open()))
    pa=parameter_checks(rows)
    small=[(2,5,7),(2,5,11),(2,7,11),(3,5,7),(3,5,11),(3,7,11),(5,7,11)]
    mods=[]
    for r,s,p in small:
        for fac in (1,r,s):mods.append((r*s*p*fac,r,s,p))
    mods+= [(2*5*p,2,5,p) for p in (53,59,113,127,131,137,149,257)]
    mods+= [(3*5*p,3,5,p) for p in (53,59,113,127,131,149)]
    mods+= [(5*7*127,5,7,127),(7*11*127,7,11,127),(2*5*127*2,2,5,127)]
    mods=sorted(set(mods));states=bad=central=tail=cdf=0;float_free=0
    missing_atom=None;first_bad=None
    for M,r,s,p in mods:
        q=r*s;par=parameters(r,s);u=F(par['u']);b=F(par['bad_incidence_b']);Kc=F(par['interval_span'])
        A=direct(M,q);B=direct(M,q*p);states+=M;bad_hist=Counter()
        for st,((H,a,d),(R,aa,dd)) in enumerate(zip(A,B)):
            need(H<=R,'monotone admissibility')
            if d%p:need(H==R,'unaffected minimum')
            else:
                bad_hist[H]+=1;bad+=1;rep=repair(M,st,a,d,q,p)
                need(R*H<=H*H+2*M,'true height versus repair')
                if first_bad is None:first_bad={'M':M,'q':q,'p':p,'state':st,'relaxed':H,'full':R,'chosen':[a,d],'repaired':[rep[3],rep[4]]}
        ct=0
        for T in range(1,max(x[0] for x in A)+1):
            ct+=bad_hist[T]
            need(ct<=F(5,2)*u*T*T/p+b*T,'problematic CDF exact bound')
            cdf+=1
        for T in sorted(set([1,isqrt((M-1)//2),isqrt(M)])):
            Pq=Pfull=0;pos=set();neg=set()
            for d in range(1,T+1):
                if gcd(d,q)>1:continue
                for a in range(1,T+1):
                    if gcd(a,d)>1:continue
                    Pq+=1
                    if gcd(d,M)==1:
                        Pfull+=1;v=a*pow(d,-1,M)%M
                        need(v not in pos,'same-sign collision');pos.add(v)
                        v=(-v)%M;need(v not in neg,'same-sign negative collision');neg.add(v)
            J=len(pos&neg);need(1+2*Pfull-J==sum(H<=T for H,a,d in B),'full primitive central identity')
            need(Pfull>=Pq-u*T*T/p-F(par['Dplus'])*T,'additional-prime central loss')
            if 2*T*T<=M:need(J==0,'collision-free endpoint')
            else:
                L=2*F(T)-F(M,T)
                need(J<=u*u*L*L/2+2*u*Kc*L+Kc*Kc,'general closed collision triangle')
            central+=1
        for T in sorted(set([isqrt(2*M)+1,2*isqrt(M)+1,3*isqrt(M)+1])):
            h=holes(M,T,q);truth=sum(H>T for H,a,d in A)
            need(truth<=h,'general two-prime Farey coverage')
            lens=h-sum(euler(d) for d in range(1,(M-1)//T+1) if gcd(d,q)>1)
            if truth>lens and missing_atom is None:missing_atom={'M':M,'T':T,'q':q,'survival':truth,'length_only':str(lens)}
            tail+=1
        if p>=53:
            mu=sum((log_hi(H) for H,a,d in B),F(0))/M-log_lo(M)/2
            need(mu<bound(M,r,s,p),'analytic formula versus full exact small model');float_free+=1
    # The Farey statement applies to arbitrary M, not only to multiples of the mark.
    other=0
    for q in (10,14,15,21,35):
        for M in range(2,101):
            vals=direct(M,q)
            for T in (isqrt(2*M)+1,2*isqrt(M)+1):
                need(sum(H>T for H,a,d in vals)<=holes(M,T,q),'arbitrary-initial-modulus Farey check');other+=1
    # Two-prime repair table checked across every eta and rational shears, independent of M.
    table=0
    for q in (10,14,15,21,22,26,33,35,55,77,143):
        for eta in range(q):
            for z in (-128,-1,0,1,128):
                t=F(z,256);sg=1 if t>=0 else -1
                if gcd(eta,q)==1:i,j=1,0
                elif gcd(eta-sg,q)==1:i,j=1,-sg
                else:i,j=2,-sg
                need(gcd(j+i*eta,q)==1 and abs(j+i*t)<=1,'two-prime repair table');table+=1
    # Exact deleted-density identity at real endpoints.
    den=0
    for r,s in ((2,5),(2,7),(3,5),(3,7),(5,7),(7,11)):
        q=r*s;par=parameters(r,s);c=F(par['deleted_c']);err=F(par['deleted_discrepancy']);h=F(par['harmonic_h']);u=F(par['u']);rho=F(par['rho'])
        from family import AL,AH
        for j in range(4,401):
            Q=F(j,4);n=Q//1
            val=sum((F(euler(d),d) for d in range(1,n+1) if gcd(d,q)>1),F(0))
            identity=sum((F(mob(g),g)*(F(r-1,r)*(Q//(r*g))+F(s-1,s)*(Q//(s*g))-(1-F(1,q))*(Q//(q*g))) for g in range(1,n+1) if gcd(g,q)==1),F(0))
            need(val==identity,'two-prime deleted-density floor identity')
            E=err*(h+u*log_lo(Q))+2*c
            need(val-AL*(1-rho)*Q<=E and AH*(1-rho)*Q-val<=E,'deleted-density error')
            den+=1
    negative={
        'four_prime_support':expect(lambda:bound(2*3*5*127,2,5,127)),
        'composite_third_parameter':expect(lambda:bound(2*5*129,2,5,129)),
        'duplicate_mark_primes':expect(lambda:parameters(5,5)),
        'previous_family_not_new':expect(lambda:parameters(2,3)),
        'strict_Farey_boundary':expect(lambda:holes(450,30,10)),
        'old_mod6_count':expect(lambda:need(3-3//2-3//5+3//10<=F(3,3)+F(2,3),'U_10(3)=2 exceeds old U_6 bound 5/3')),
        'three_mark_prime_repair':expect(lambda:need(any(gcd(j+3*i,30)==1 and abs(F(j)+F(i,4))<=1 for i in (1,2) for j in range(-3,4)),'q=30, eta=3, tau=1/4 has no admissible coefficient i<=2')),
    }
    cert={'status':'PASS','parameter_checks':pa,'complete_moduli':mods,'state_modulus_incidents':states,'bad_minima_repaired':bad,'bad_CDF_thresholds':cdf,'exact_central_counts':central,'three_prime_Farey_checks':tail,'arbitrary_modulus_Farey_checks':other,'repair_residue_shear_checks':table,'deleted_density_real_cutoffs':den,'analytic_small_mean_checks':float_free,'first_problematic_fixture':first_bad,'integer_atom_fixture':missing_atom,'negative_controls':negative,'finite_tests_not_general_proof':True}
    return cert
if __name__=='__main__':
    import sys
    obj=run(sys.argv[1]);Path(sys.argv[2]).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n');print(json.dumps(obj,indent=2),flush=True)
