"""Prime excision: finite checks, independent definitions and exact family gates.
The all-modulus argument is supplied in the main manuscript, not inferred from tests.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import gcd,prod,isqrt
from collections import Counter
import arithmetic as ar
import reference as ref
need=ar.need


def choose(eta:int,tau:F)->tuple[int,int]:
    need(0<=eta<6 and -F(1,2)<=tau<=F(1,2),'repair table domain')
    if eta==0:return (1,-1) if tau>=0 else (1,1)
    if eta in (1,5):return (1,0)
    if eta==2:return (1,-1) if tau>=0 else (2,1)
    if eta==4:return (1,1) if tau<=0 else (2,-1)
    return (2,-1) if tau>=0 else (2,1)


def direct(M:int,q:int):
    """All minima and witnesses, from d=1,...,M//2 (d=1 for M=1).
    d=1 already gives height <= max(1,M//2), so larger d cannot improve it.
    Ties retain the first denominator, for a reproducible chosen minimizer.
    """
    need(M>=1 and q>=1,'positive model parameters')
    cap=max(1,M//2)
    R=[max(1,min(s,M-s)) for s in range(M)]
    W=[(s if 2*s<=M else s-M,1) for s in range(M)]
    for d in range(2,cap+1):
        if gcd(d,q)!=1:continue
        for s in range(M):
            if d>=R[s]:continue
            a=s*d%M
            if 2*a>M:a-=M
            h=max(abs(a),d)
            if h<R[s]:R[s]=h;W[s]=(a,d)
    return R,W


def logratio_bounds(n:int,d:int):
    a,b=ar.logs(n);c,e=ar.logs(d)
    return F(a-e,ar.LS),F(b-c,ar.LS)


def psi_bounds(p:int):
    need(p>=5,'p range')
    l1,u1=logratio_bounds(p+6,p);l2,u2=logratio_bounds(p+6,6)
    return l1+F(6,p)*l2,u1+F(6,p)*u2


def egcd(a,b):
    oldr,r,oldx,x,oldy,y=a,b,1,0,0,1
    while r:
        q=oldr//r;oldr,r=r,oldr-q*r;oldx,x=x,oldx-q*x;oldy,y=y,oldy-q*y
    if oldr<0:oldr,oldx,oldy=-oldr,-oldx,-oldy
    return oldr,oldx,oldy


def repair_witness(M,s,p,a,d):
    need(d>0 and gcd(d,6)==1 and (a-s*d)%M==0 and d%p==0,'bad chosen witness')
    n=(a-s*d)//M
    g,x,y=egcd(n,d);need(g==1,'primitive in physical lattice, not ambient gcd')
    # n*i-d*j=1, hence v=j*(M,0)+i*(s,1).
    i,j=x,-y;b,e=j*M+i*s,i
    need(a*e-b*d==M and gcd(d,e)==1,'oriented completion')
    nu=a*a+d*d;t=(2*(a*b+d*e)+nu)//(2*nu)
    b-=t*a;e-=t*d
    tau=F(a*b+d*e,nu);need(-F(1,2)<=tau<=F(1,2),'centred shear')
    eta=e*pow(d,-1,6)%6;i,j=choose(eta,tau)
    aa=j*a+i*b;dd=j*d+i*e
    need(dd!=0,'unit denominator is nonzero')
    if dd<0:aa,dd=-aa,-dd
    H=max(abs(a),d);K=max(abs(aa),dd)
    need(gcd(dd,6*p)==1 and (aa-s*dd)%M==0,'repaired admissibility')
    need(K*H<=H*H+2*M,'pointwise repair upper bound')
    return K


def candidate():
    R,W=direct(210,210);P=prod(R)
    need(P**170>210**18690,'candidate integer comparison')
    expected={2:76,3:57,5:25,7:21,11:24,13:26,17:20,19:26,23:20,29:4,31:4,37:8,41:2,53:2}
    need(prod(p**e for p,e in expected.items())==P,'candidate factorization')
    lower,upper=ref.rational_log(P,terms=48);lb,ub=ref.rational_log(210,terms=48)
    # Positive log(P) permits these endpoint choices.
    clo=lower/(210*ub)-F(1,2);chi=upper/(210*lb)-F(1,2)
    z=10**24
    return {'status':'PASS','level':'CL3','mass':210,'R':R,'height_sum':sum(R),'factorization':{str(p):e for p,e in expected.items()},
            'P':str(P),'P170_gt_21018690':True,'c0_lower':str(F(clo.numerator*z//clo.denominator,z)),
            'c0_upper':str(F(ar.ceildiv(chi.numerator*z,chi.denominator),z)),'all_minima_from_definition':True}


def finite_checks():
    # Exact table endpoints and many rational interior parameters; proof covers the continuum.
    tc=0
    for denominator in range(1,65):
        for numerator in range(-denominator//2-1,denominator//2+2):
            t=F(numerator,denominator)
            if not -F(1,2)<=t<=F(1,2):continue
            for eta in range(6):
                i,j=choose(eta,t)
                need(1<=i<=2 and gcd(j+i*eta,6)==1 and abs(j+i*t)<=1,'repair table exact check');tc+=1
    cases=[]
    for M in range(2,601):
        for p in (5,7,11,13):
            if M%p==0:cases.append((M,p))
    cases+= [(6*p, p) for p in (17,19,23,29,31,37)]+[(210,5),(210,7),(462,7),(990,5),(1050,7),(1296*5,5)]
    cases=sorted(set(cases))
    states=badstates=thresholds=incidencechecks=0;rows=[]
    cache={}
    for M,p in cases:
        if M not in cache:cache[M]=direct(M,6)
        R,W=cache[M];S,V=direct(M,6*p)
        bad=Counter();numbad=0
        for s,(h,k,(a,d)) in enumerate(zip(R,S,W)):
            need(h<=k,'relaxation order')
            if d%p:
                need(h==k,'good chosen witness retains exact minimum')
            else:
                numbad+=1;bad[h]+=1
                repaired=repair_witness(M,s,p,a,d)
                need(k<=repaired,'repaired candidate compared with exact minimum')
            states+=1
        count=0
        for H,n in sorted(bad.items()):
            count+=n;need(p*count<=3*H*H,'bad CDF bound');thresholds+=1
        upper=sum(logratio_bounds(k,h)[1] for h,k in zip(R,S))/M
        psi_l,psi_u=psi_bounds(p)
        need(upper<psi_l,'finite mean difference below analytic penalty')
        badstates+=numbad
        rows.append({'M':M,'p':p,'bad_chosen_states':numbad,'mean_log_ratio_upper':str(upper),'psi_lower':str(psi_l)})
        for A in range(1,min(9,M//p)+1):
            actual=0;m=M//p
            for e in range(1,A+1):
                for b in range(-A,A+1):
                    actual+=sum((s*e-b)%m==0 for s in range(M))
            formula=p*sum(gcd(e,m)*(2*(A//gcd(e,m))+1) for e in range(1,A+1))
            need(actual==formula<=3*p*A*A,'all-state incidence formula and upper bound');incidencechecks+=1
    r,w=direct(30,6);s,_=direct(30,30)
    need(r[6]==5 and s[6]==6 and w[6]==(0,5),'ambient nonprimitive regression')
    need((0-6*1)%30!=0,'dividing ambient gcd changes the state')
    # A genuine missing-denominator failure, caught rather than reported as a proof.
    try:
        need(s[6]<=r[6],'FALSE equality after deleting the prime condition')
    except ArithmeticError as e:negative=str(e)
    else:raise ArithmeticError('negative regression did not fail')
    wheel_checks=[]
    for p in range(5,127):
        if any(p%d==0 for d in range(2,isqrt(p)+1)):continue
        units=[n for n in range(6*p) if gcd(n,6*p)==1]
        J=max(b-a for a,b in zip(units,units[1:]+[units[0]+6*p]))
        need(J==6,'complete finite wheel gap 6p')
        wheel_checks.append({'p':p,'J_6p':J})
    return {'status':'PASS','level':'CL3 finite domains only','wheel_checks':wheel_checks,'table_rational_tests':tc,'modulus_prime_pairs':len(cases),
            'all_states_in_declared_pairs':states,'bad_chosen_witnesses_repaired':badstates,'bad_CDF_thresholds':thresholds,
            'incidence_formulas_checked_by_all_states':incidencechecks,
            'negative_regression':{'M':30,'s':6,'p':5,'relaxed_height':5,'full_height':6,'chosen_relaxed_vector':[0,5],
             'lattice_coordinates':[-1,5],'ambient_gcd':5,'output':'ArithmeticError: '+negative},'rows':rows}


def family_gates(uniform):
    thresholds=[3000000,1471,277,157,127]
    rows=[]
    for row,p in zip(uniform['rows'],thresholds):
        k=row['k'];cap=F(row['rational_cap']);pl,pu=psi_bounds(p)
        lb,ub=logratio_bounds(6*p,1)
        margin=F(2,85)*lb-cap-pu
        need(margin>0,'family gate k='+str(k))
        # Independent full rational series for all six log arguments of each gate.
        L1,U1=ref.rational_log(p+6,p,48);L2,U2=ref.rational_log(p+6,6,48);L3,U3=ref.rational_log(6*p,1,48)
        need(pl<=L1+F(6,p)*L2<=U1+F(6,p)*U2<=pu,'penalty independent series containment')
        need(F(2,85)*L3-cap-U1-F(6,p)*U2>0,'independent exact family gate')
        rows.append({'k':k,'condition':'min(v2(B^N),v3(B^N)) >= k','p_threshold':p,'nu6_cap':str(cap),
                     'psi_upper':str(pu),'log_6p_lower':str(lb),'positive_margin_lower':str(margin),'all_p_above_threshold':True})
    primes=[n for n in range(5,127) if all(n%d for d in range(2,isqrt(n)+1))]
    sieve=[True]*127;sieve[0]=sieve[1]=False
    for d in range(2,12):
        if sieve[d]:
            for n in range(d*d,127,d):sieve[n]=False
    need(primes==[n for n in range(5,127) if sieve[n]],'two prime enumerations')
    need(len(primes)==28 and primes[-1]==113,'remaining small primes')
    pl,pu=psi_bounds(113);lb,ub=logratio_bounds(678,1)
    need(F(-19,500)+pl-F(2,85)*ub>0,'p113 chosen majorant does not pass')
    # Complete finite prime frontiers: sieve and a separate segmented construction.
    import hashlib,bisect
    limit=thresholds[0]
    sieve=bytearray(b'\x01')*limit;sieve[0:2]=b'\x00\x00'
    for d in range(2,isqrt(limit-1)+1):
        if sieve[d]:sieve[d*d:limit:d]=b'\x00'*(((limit-1-d*d)//d)+1)
    allprimes=[n for n in range(5,limit) if sieve[n]]
    basic=[n for n in range(2,isqrt(limit-1)+1) if all(n%d for d in range(2,isqrt(n)+1))]
    other=[]
    for lo in range(5,limit,32768):
        hi=min(limit,lo+32768);block=bytearray(b'\x01')*(hi-lo)
        for p0 in basic:
            first=max(p0*p0,ar.ceildiv(lo,p0)*p0)
            if first<hi:block[first-lo:hi-lo:p0]=b'\x00'*(((hi-1-first)//p0)+1)
        other.extend(lo+i for i,b in enumerate(block) if b)
    need(allprimes==other,'complete and segmented prime frontiers agree')
    for r in rows:
        n=bisect.bisect_left(allprimes,r['p_threshold'])
        r['remaining_prime_count_for_N_ge_k']=n
        r['largest_unexcluded_prime_for_N_ge_k']=allprimes[n-1]
    digest=hashlib.sha256(('\n'.join(map(str,allprimes))+'\n').encode()).hexdigest()
    return {'status':'PASS','level':'CL3 exact comparisons with general monotonic extension proved in text',
            'prime_frontier_full_range':[5,limit-1],'prime_frontier_methods':['Eratosthenes','segmented sieve with trial-division base primes'],
            'full_prime_list_sha256':digest,
            'rows':rows,'N_ge_5_unexcluded_primes_by_this_filter':primes,
            'negative_gate':'FAIL: p=113,k=5: this majorant does not prove Delta < 2/85; not a counterexample',
            'entire_6p_family_closed':False,'global_Cstar_identified':False,'formalized':False}


def bank_checks(prefix):
    from pathlib import Path
    points=[tuple(map(int,x.split())) for x in Path(str(prefix)+'.points').read_text().splitlines()]
    roots=0
    for k,bound in points:
        if k>4:continue
        M=6**k;R,_=direct(M,6)
        upper=sum(ar.logs(h)[1] for h in R)/F(M*ar.LS)-F(ar.logs(M)[0],2*ar.LS)
        need(upper<=F(bound,ar.S),'complete root bank against direct definition');roots+=M
    R6,_=direct(6,6);R36,_=direct(36,6);conditional=0
    for line in Path(str(prefix)+'.bank').read_text().splitlines():
        k,s,bound=map(int,line.split())
        if k==0:
            upper=F(ar.logs(R6[s])[1],ar.LS)-F(ar.logs(6)[0],2*ar.LS)
        else:
            need(k==1,'bank horizon domain')
            upper=sum(ar.logs(R36[s+6*b])[1] for b in range(6))/F(6*ar.LS)-F(ar.logs(36)[0],2*ar.LS)
        need(upper<=F(bound,ar.S),'conditional all-successor bank');conditional+=1
    return {'root_bank_lengths':[0,1,2,3,4],'root_bank_states':roots,
            'conditional_bank_comparisons':conditional,'definition_not_gauss_used_for_minima':True}
