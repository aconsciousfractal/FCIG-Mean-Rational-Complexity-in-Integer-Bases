#!/usr/bin/env python3
"""Portable global scalar bound; fixed S in both integrals and the logarithmic divisor."""
import argparse
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
from math import isqrt, prod
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
SOURCE = Path(__file__).resolve().parent
TARGET = SOURCE / 'certificate.json'
PS = (2,3,5,7,11,13,17,19)
CLOSED = (2310,2730,3570,3990,4290,4830,5610,6006,6090,6270,6510)
GRID = 10**40
START = time.monotonic()
CHECKS = 0

def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ArithmeticError(message)
    if time.monotonic()-START > 600:
        raise TimeoutError('scalar runtime budget 600 seconds exceeded')

def canonical(obj):
    return (json.dumps(obj,sort_keys=True,indent=2)+'\n').encode()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def module(path,name):
    spec = importlib.util.spec_from_file_location(name,path)
    obj = importlib.util.module_from_spec(spec)
    sys.modules[name] = obj
    spec.loader.exec_module(obj)
    return obj

def floor_grid(x):
    return F(x.numerator*GRID//x.denominator,GRID)

def ceil_grid(x):
    return F(-(-x.numerator*GRID//x.denominator),GRID)

def decimal_up(x,digits=15):
    x=F(x);n=-(-x.numerator*10**digits//x.denominator)
    need(n>=0,'nonnegative display')
    return str(n//10**digits)+'.'+str(n%10**digits).zfill(digits)

@lru_cache(None)
def log_interval(x):
    x=F(x);need(x>0,'positive logarithm argument')
    if x<1:
        lo,hi=log_interval(1/x)
        return -hi,-lo
    e=x.numerator.bit_length()-x.denominator.bit_length()
    if x<F(2)**e:e-=1
    y=x/F(2)**e
    need(e>=0 and 1<=y<2,'exact logarithm range reduction')
    def series(t):
        z=(t-1)/(t+1);z2=z*z;p=z;total=F(0)
        for j in range(80):
            total+=2*p/(2*j+1);p*=z2
        return total,total+2*p/(161*(1-z2))
    lo,hi=series(y);l2,h2=series(F(2))
    return floor_grid(lo+e*l2),ceil_grid(hi+e*h2)

@lru_cache(None)
def sqrt_interval(x):
    x=F(x);need(x>=0,'nonnegative square root')
    n=isqrt(x.numerator*GRID**2//x.denominator)
    lo,hi=F(n,GRID),F(n+1,GRID)
    need(lo*lo<=x<hi*hi,'exact square root enclosure')
    return lo,hi

def interval(values):
    return [str(v) for v in values]

def positive(Q,k2,periodic):
    a,b,d=(F(periodic[k]) for k in ('quadratic_coefficient','linear_coefficient','constant_coefficient'))
    coefficient=1-d/Q;need(coefficient>0,'positive cut equation coefficient')
    sl,sh=sqrt_interval(k2/Q)
    dl,dh=sqrt_interval(b*b*k2/Q+4*coefficient*a*k2)
    cutlo,cuthi=(b*sl+dl)/(2*coefficient),(b*sh+dh)/(2*coefficient)
    z=F(-(-(cuthi*10**6).numerator//(cuthi*10**6).denominator),10**6)
    need(1<cutlo<=cuthi<=z and z-cuthi<F(1,10**6),'outward stationary cut')
    cQ=k2*Q/z**2;ll,lh=log_interval(cQ)
    E=a/2+b;tl,th=sqrt_interval(k2*Q)
    margin=b*b*k2*Q-4*z*z*E*E
    derivative=E-b*tl/(2*z)
    need(cQ>1 and ll>1 and margin>0 and derivative<0,'all-real-M derivative guards')
    zl,zh=log_interval(z)
    corrlo=a*k2/(2*z*z)+b*sl/z+d*ll/(2*Q)-E/Q
    corrhi=a*k2/(2*z*z)+b*sh/z+d*lh/(2*Q)-E/Q
    lo,hi=zl+corrlo,zh+corrhi
    need(0<corrlo<=corrhi and 0<lo<=hi,'positive finite endpoint integral enclosure')
    return dict(cut_z=str(z),cut_selection=dict(stationary_cut_interval=interval((cutlo,cuthi)),rational_grid_denominator=10**6,exact_optimal_cut_claimed=False),
        positive_interval=interval((lo,hi)),positive_correction_interval=interval((corrlo,corrhi)),
        derivative_guards=dict(cutoff_ratio_squared_at_Q=str(cQ),log_cutoff_ratio_squared_interval=interval((ll,lh)),
        E=str(E),sqrt_k2Q_interval=interval((tl,th)),squared_gate_margin=str(margin),derivative_M_upper_at_Q=str(derivative),
        scaled_derivative_formula='E-b*sqrt(k2*M)/(2*z)+d/2*(1-ln(k2*M/z^2))',all_real_M_ge_Q_S=True,
        finite_endpoint_subtractions_included=True))

def validate_profile(p):
    scale=p['tick_scale'];bs=p['breakpoints_radius_tick']
    need(isinstance(scale,int) and scale>0 and bs[0]==[0,0],'envelope initial value and exact scale')
    need(all(isinstance(r,int) and isinstance(t,int) and r>=0 and t>=0 for r,t in bs),'integer nonnegative envelope ticks')
    need(all(a[0]<b[0] and a[1]<=b[1] for a,b in zip(bs,bs[1:])),'ordered nondecreasing complete envelope')

def h(p,r):
    bs=p['breakpoints_radius_tick'];i=bisect_right([b[0] for b in bs],r)-1
    need(i>=0,'envelope radius covered')
    return F(bs[i][1],p['tick_scale'])

def strict_radius(Q,alpha):
    x=Q*alpha;r=isqrt(x.numerator//x.denominator)
    if r*r==x:r-=1
    need(r*r<x<=(r+1)**2,'strict injective boundary')
    return r

def negative(square,cross,Q):
    modes=[(square,F(1,2))]
    if cross is not None:modes.append((cross,F(4,5)))
    cuts={F(0),F(1)}
    for p,alpha in modes:
        validate_profile(p);cuts.add(alpha)
        cuts.update(F(r*r,Q) for r,t in p['breakpoints_radius_tick'] if 0<F(r*r,Q)<alpha)
    def affine(p,alpha,u):
        if u<alpha:
            x=Q*u;r=isqrt(x.numerator//x.denominator)
            return h(p,r),F(0)
        return F(0),alpha*h(p,strict_radius(Q,alpha))
    bands=[];lo=hi=F(0);cuts=sorted(cuts)
    for l,r in zip(cuts,cuts[1:]):
        coefficients=[affine(p,alpha,(l+r)/2) for p,alpha in modes]
        split=[l,r]
        if len(coefficients)==2:
            (m1,b1),(m2,b2)=coefficients
            if m1!=m2:
                switch=(b2-b1)/(m1-m2)
                if l<switch<r:split=[l,switch,r]
        for x,y in zip(split,split[1:]):
            mid=(x+y)/2;m,b=max(coefficients,key=lambda mb:mb[0]*mid+mb[1])
            need(0<=x<y<=1 and m>=0 and b>=0,'valid nonnegative fixed-support integral band')
            base=m*(y-x)/2;lo+=base;hi+=base
            if b:
                need(x>0,'positive plateau logarithm endpoint')
                ll,lh=log_interval(y/x);lo+=b*ll/2;hi+=b*lh/2
            bands.append(interval((x,y,m,b)))
    need(bands[0][0]=='0' and bands[-1][1]=='1' and all(a[1]==b[0] for a,b in zip(bands,bands[1:])),'complete disjoint partition of u in [0,1]')
    need(0<=lo<=hi,'negative integral interval')
    return dict(negative_interval=interval((lo,hi)),negative_bands=bands,
        negative_method='square' if cross is None else 'pointwise maximum of square and cross',
        strict_square_radius=strict_radius(Q,F(1,2)),strict_cross_radius=None if cross is None else strict_radius(Q,F(4,5)),
        negative_bound_for_all_real_M_ge_Q_S=True)

def sieve(n):
    a=bytearray(b'\1')*(n+1);a[:2]=b'\0\0'
    for p in range(2,isqrt(n)+1):
        if a[p]:a[p*p::p]=b'\0'*((n-p*p)//p+1)
    return [i for i in range(2,n+1) if a[i]]

def build():
    old=read(SOURCE/'inputs/global_sectors.json');periodic=old['periodic_CDF']
    square_doc=read(SOURCE/'inputs/square.json')
    crossdoc=read(SOURCE/'inputs/cross.json')
    need([F(periodic[k]) for k in ('quadratic_coefficient','linear_coefficient','constant_coefficient')]==[F(6414336,5010005),F(45173053,51531480),F(221184,715715)],'PCDF coefficients')
    Hs={r['prime_count_at_most']:r['H'] for r in old['new_uniform_covering']}
    sectors={};comparisons=[]
    for w,Q in ((4,3003),(5,6630),(6,30030)):
        H=Hs[w-1];k2=1+max(F(H*H,4),F((H-1)**2,3))
        squares=square_doc['support_envelopes'][str(w)]
        expected=[mask for mask in range(256) if mask.bit_count()<=w and not (w==5 and mask.bit_count()==w and prod(p for i,p in enumerate(PS) if mask&(1<<i)) in CLOSED)]
        need([p['mask'] for p in squares]==expected,'complete exact mask inventory')
        crosses={} if w==4 else {p['mask']:p for p in crossdoc['support_envelopes'][str(w)]}
        if w!=4:need(list(crosses)==expected,'same square/cross support inventory')
        rows=[]
        for square in squares:
            mask=square['mask'];marked=[p for i,p in enumerate(PS) if mask&(1<<i)]
            radical=prod(marked);t=w-len(marked);outside=[]
            for p in sieve(101):
                if p<=19 or (w==5 and t==1 and radical*p in CLOSED):continue
                if len(outside)<t:outside.append(p)
            need(len(outside)==t and square['kernel_radical']==radical and square['external_prime_lower_bounds']==outside,'supportwise external ordered lower bounds')
            physical=radical*prod(outside);QS=max(Q,physical)
            cross=crosses.get(mask)
            if cross is not None:
                need(cross['mask']==mask and cross['omega']==w and cross['outside_prime_lower_bounds']==outside and cross['tick_scale']==square['tick_scale'],'identical physical support across both shapes')
            row=dict(omega=w,mask=mask,marked_kernel_primes=marked,kernel_radical=radical,outside_prime_lower_bounds=outside,physical_base_min=physical,Q_S=QS,
                square_profile=square,cross_profile=cross)
            row.update(positive(QS,k2,periodic));row.update(negative(square,cross,QS))
            ql,qh=log_interval(QS);U=F(row['positive_interval'][1])-F(row['negative_interval'][0])
            need(ql>0,'positive logarithmic divisor')
            upper=U/ql if U>0 else F(0)
            row.update(ln_Q_interval=interval((ql,qh)),U_S=str(U),numerator_positive=U>0,final_upper=str(upper),final_upper_decimal=decimal_up(upper),
                nonpositive_numerator_uses_zero=U<=0,support_fixed_through_positive_negative_division=True)
            rows.append(row)
        winner=max(rows,key=lambda r:F(r['final_upper']));upper=F(winner['final_upper'])
        sectors[str(w)]=dict(omega=w,Q=Q,H=H,kappa_squared=str(k2),support_count=len(rows),rows=rows,
            max_after_division=True,selected_support_mask=winner['mask'],selected_Q_S=winner['Q_S'],upper=str(upper),upper_decimal=decimal_up(upper))
        comparisons.append(dict(sector='omega='+str(w),upper=str(upper),upper_decimal=decimal_up(upper),provenance='complete supportwise scalar calculation'))
    inherited=[r for r in old['finite_sectors'] if r['omega']>=7]
    need([r['omega'] for r in inherited]==list(range(7,16)),'all unchanged finite sectors seven through fifteen')
    for r in inherited:
        comparisons.append(dict(sector='omega='+str(r['omega']),upper=r['final_upper'],upper_decimal=r['final_upper_decimal'],provenance='finite scalar sector'))
    tail=old['infinite_tail'];candidate=old['candidate']
    comparisons.extend([dict(sector='omega>=16',upper=tail['upper'],upper_decimal=tail['upper_decimal'],provenance='classical tail'),
        dict(sector='closed/candidate sectors',upper=candidate['c0_upper'],upper_decimal=decimal_up(F(candidate['c0_upper'])),provenance='closed domains and exact candidate')])
    winner=max(comparisons,key=lambda r:F(r['upper']));cap=F(winner['upper']);oldcap=F(old['global_upper'])
    names=['inputs/square.json','inputs/cross.json','inputs/global_sectors.json','producer.py','check.py','check_premises.py','check_closure_tables.py','check_radical_closures.py','verify_closures.py','reconstruct_closures.py','closure_input_catalog.json','finite_certificate_index.json','proofs/PROOF.md','proofs/NEGATIVE_PROFILES.md']
    bindings=[dict(path=name,bytes=(SOURCE/name).stat().st_size,sha256=sha(SOURCE/name)) for name in names]
    return dict(schema='global-supportwise-v1',status='SCALAR_CERTIFICATE',normalization='Delta(B,N)=mu(B^N)/ln(B)',
        all_N_ge_1=True,all_positive_valuations=True,all_states_including_zero=True,kernel_primes=list(PS),conditional_excluded_five_radicals=list(CLOSED),
        source_bindings=bindings,
        periodic_CDF=periodic,sectors=sectors,inherited_finite_sectors=inherited,infinite_tail=tail,candidate=candidate,
        declared_four_prime_families=old['declared_four_prime_families'],every_conditional_exclusion_separately_paid=old['every_conditional_exclusion_separately_paid'],
        all_sector_comparison=comparisons,dominant_sector=winner['sector'],global_upper=str(cap),global_upper_decimal=decimal_up(cap),
        previous_global_upper=str(oldcap),previous_global_upper_decimal=old['global_upper_decimal'],strict_improvement=cap<oldcap,exact_improvement=str(oldcap-cap),
        log_terms=80,log_method='atanh with positive geometric remainder and outward rational grid10^-40',sqrt_decimal_places=40,
        own_explicit_guard_count=CHECKS,stock_count_producers_replayed=False,stock_global_verifier_replayed=False,historical_large_raw_replayed=False,
        optimization_mode_scope='own scalar guards and compressed-profile API authentication only',
        
        exact_Cstar_equals_c0_proved=False,offset_sign_classification_claimed=False,public_or_manuscript_changes=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true')
    parser.add_argument('--output',type=Path)
    global TARGET
    args=parser.parse_args()
    if args.output: TARGET=args.output.resolve()
    result=build();data=canonical(result)
    need(len(data)<=30*1024*1024,'compressed evidence output below 30 MiB')
    if args.build:
        need(not TARGET.exists(),'build refuses to overwrite a frozen certificate')
        TARGET.write_bytes(data)
    else:need(TARGET.read_bytes()==data,'entire own scalar certificate regenerated identically')
    print(json.dumps(dict(status='PASS',certificate_sha256=sha(TARGET),certificate_bytes=len(data),global_upper_decimal=result['global_upper_decimal'],
        dominant_sector=result['dominant_sector'],strict_improvement=result['strict_improvement'],support_counts={w:s['support_count'] for w,s in result['sectors'].items()},
        runtime_seconds=round(time.monotonic()-START,3)),sort_keys=True))

if __name__=='__main__':main()
