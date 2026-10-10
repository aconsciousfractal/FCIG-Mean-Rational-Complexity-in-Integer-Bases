"""Independent finite controls. Their finite ranges are explicit, not general proofs."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd,isqrt,prod
from pathlib import Path
import csv,json
from central import *
from inventory import write_json,power_tuples,heap_moduli,primitive_base

def mobius(n):
    z=n;k=0
    for p in range(2,isqrt(n)+1):
        if z%p==0:
            z//=p;k+=1
            if z%p==0:return 0
    if z>1:k+=1
    return -1 if k%2 else 1

def divisors_signed(prs):
    return [(prod(t),(-1)**len(t)) for j in range(len(prs)+1) for t in combinations(prs,j)]
def units_by_floors(x,prs):return sum(sg*(x//d) for d,sg in divisors_signed(prs))
def roots_by_search(M):
    out=[]
    for e in range(1,M.bit_length()):
        lo=1;hi=1<<((M.bit_length()+e-1)//e)
        while lo<hi:
            mid=(lo+hi)//2
            if mid**e<M:lo=mid+1
            else:hi=mid
        if lo>=2 and lo**e==M:out.append((lo,e))
    return sorted(out)

def direct_heights(M,q):
    hs=[];ds=[]
    for s in range(M):
        H=max(1,min(s,M-s));best=1
        for d in range(1,M//2+1):
            if gcd(d,q)!=1:continue
            a=s*d%M;h=max(d,min(a,M-a))
            if h<H:H=h;best=d
        hs.append(H);ds.append(best)
    return hs,ds

def run(out):
    out=Path(out);rows=list(csv.DictReader((out/'profiles_442.csv').open()))
    residues=polys=harmonic=primitive=0
    seenpairs=set();pair_residues=0
    for row in rows:
        prs=tuple(int(row[k]) for k in ('r','s','p'));R=prod(prs);lc=central_constants(*prs);u=lc['u']
        ds=divisors_signed(prs);low=high=0;ph=int(u*R);Knum=int(lc["K"]*R)
        need(F(Knum,R)==lc["K"],"integer-scaled floor constant")
        for n in range(R):
            cn=sum(sg*(n//d) for d,sg in ds);eps=R*cn-ph*n
            low=min(low,eps);high=max(high,eps);residues+=1
            # The two endpoint expressions for floor(x) U(floor(x))-u*x*x.
            need(abs(eps)<=Knum,'left floor-product endpoint')
            right=R*n*cn-ph*(n+1)*(n+1)
            need(abs(right)<=Knum*(n+1),'right floor-product endpoint');polys+=2
        need(F(high,R)==lc['dp'] and u-F(low,R)==lc['dminus'],'independent full period extrema')
        hh=F(0);un=0
        for n in range(1,17):
            if all(n%p for p in prs):hh+=F(1,n);un+=1
            llo,lhi=alternate_log(n)
            need(hh-u*llo<=lc['h'],'independent harmonic prefix');harmonic+=1
        llo,lhi=alternate_log(16)
        need(hh-u*llo+(16*u-un+lc['dp'])/16<=lc['h'],'harmonic tail constant');harmonic+=1
        for T in (1,2,3,7,16,23,40):
            count=sum(1 for d in range(1,T+1) if gcd(d,R)==1 for a in range(1,T+1) if gcd(a,d)==1)
            inv=sum(mobius(e)*(T//e)*units_by_floors(T//e,prs) for e in range(1,T+1) if gcd(e,R)==1)
            need(count==inv,'primitive pair Mobius identity')
            err=(lc['K']*lc['h']+u)*T+lc['K']*u*T*log_hi(T)
            need(AL*lc['rho']*T*T-err<=count<=AH*lc['rho']*T*T+err,'primitive pair bounds');primitive+=1
        # Two-prime quantities inherited but reconstructed by independent floor formulas.
        pair=prs[:2]
        if pair not in seenpairs:
            seenpairs.add(pair);r,s=pair;q=r*s;U=S=0;dp=F(0);mn=F(0);b=F(0);dw=F(0);d=parameters(*pair);uq=F(d['u']);cr=F(r-1,r*r)+F(s-1,s*s)-F(q-1,q*q)
            for n in range(q):
                U=n-n//r-n//s+n//q
                S=n*(n+1)//2-r*(n//r)*(n//r+1)//2-s*(n//s)*(n//s+1)//2+q*(n//q)*(n//q+1)//2
                eps=U-uq*n;dp=max(dp,eps);mn=min(mn,eps)
                if n:b=max(b,3*eps,(2*n*U+S-F(5,2)*uq*n*n)/n)
                w=F(r-1,r)*(n//r)+F(s-1,s)*(n//s)-F(q-1,q)*(n//q)-cr*n
                dw=max(dw,abs(w),abs(w-cr));pair_residues+=1
            need(dw==F(d['deleted_discrepancy']) and cr==F(d['deleted_c']),'independent deleted-density constants')
            need(dp==F(d['Dplus']) and uq-mn==F(d['Dminus']) and b==F(d['bad_incidence_b']),'independent two-prime discrepancy')
    # Every modulus below a selected tail has all its exact roots checked by bisection.
    allrows=list(csv.DictReader((out/'finite_domain.csv').open()))+list(csv.DictReader((out/'analytic_moduli.csv').open()))
    root_tests=representations=0
    for row in allrows:
        M=int(row['M']);prs=tuple(int(row[k]) for k in ('r','s','p'));B,g,es=primitive_base(M,prs)
        rr=roots_by_search(M)
        expected=sorted((B**(g//e),e) for e in range(1,g+1) if g%e==0)
        need(rr==expected,'all base-length representations by integer root search')
        root_tests+=1;representations+=len(rr)
    # Small complete physical laboratory, not a sample used as an infinite theorem.
    moduli=[70,110,130,154,182,238,286,374,430,518,555,595,665,715,935,1001]
    states=central_tests=repairs=bad_cdf_tests=0
    for M in moduli:
        fac=[];z=M
        for p in range(2,M+1):
            if z%p==0:
                fac.append(p)
                while z%p==0:z//=p
        need(len(fac)==3 and fac[:2]!=[2,3],'small complete three-prime laboratory')
        r,s,p=fac;R=r*s*p;q=r*s
        full,_=direct_heights(M,R);relaxed,den=direct_heights(M,q)
        need(all(a<=b for a,b in zip(relaxed,full)),'marked/full model order')
        for j,H in enumerate(relaxed):
            if den[j]%p==0:
                need(full[j]*H<=H*H+2*M,'pointwise repair');repairs+=1
        lc=central_constants(r,s,p);u=lc['u'];ki=lc['ki']
        for T in range(1,isqrt(M)+1):
            P=sum(1 for d in range(1,T+1) if gcd(d,R)==1 for a in range(1,T+1) if gcd(a,d)==1)
            actual=sum(H<=T for H in full);J=1+2*P-actual
            L=max(F(0),2*T-F(M,T))
            Jupper=F(0) if 2*T*T<=M else u*u*L*L/2+2*u*ki*L+ki*ki
            need(0<=J<=Jupper,'full-unit collision bound')
            er=(lc['K']*lc['h']+u)*T+lc['K']*u*T*log_hi(T)
            lower=1+2*(AL*lc['rho']*T*T-er)-Jupper
            need(lower<=actual<=sum(H<=T for H in relaxed),'central lower CDF and relaxation');central_tests+=1
        # Integrated complete mean bounded independently by the new formula.
        mu=sum(log_hi(h) for h in full)/M-log_lo(M)/2
        need(mu<=true_upper(M,r,s,p),'small-model complete mean majorant')
        states+=M
    from secondary_log import enclosure,PREC
    secondary_inputs=list(range(1,513))+[999,1000,1001,999999,170953875,74490839040,2**120-1,2**120,2**120+1]
    for x in secondary_inputs:
        a,b=alternate_log(x);lo,hi=enclosure(x)
        need(F(lo,PREC)<=a<=b<=F(hi,PREC),'fixed-point secondary vs exact Fraction series')
    result=dict(secondary_log_rational_controls=len(secondary_inputs),status='PASS',full_period_classes=residues,floor_product_endpoint_checks=polys,full_harmonic_checks=harmonic,primitive_pair_checks=primitive,two_prime_pairs=len(seenpairs),two_prime_period_classes=pair_residues,all_complement_root_tests=root_tests,exact_base_length_representations=representations,small_physical_moduli=moduli,small_physical_state_modulus_incidences=states,full_central_CDF_checks=central_tests,pointwise_repairs=repairs,finite_tests_are_not_the_general_proof=True)
    write_json(out/'independent_checks.json',result);print(json.dumps(result,indent=2),flush=True)
    return result
if __name__=='__main__':
    import sys
    run(sys.argv[1])
