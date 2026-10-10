#!/usr/bin/env python3
"""Every prime-count sector: paid6270/130p, coupled W, and uniform negative log."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import isqrt, prod
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
from common import need, sha, read, module, canonical, check_imports, manifest
from exact import alternate_log, out_decimal

OLD = module(ROOT/'global/previous/verify_global.py', 'p72_previous_global')
AR = OLD.AR


def log(x):
    return alternate_log(F(x), terms=160)


def sqrt_interval(x):
    x = F(x); scale = 10**40
    scaled = x*scale*scale
    n = isqrt(scaled.numerator//scaled.denominator)
    lo, hi = F(n,scale), F(n+1,scale)
    need(lo*lo <= x < hi*hi, 'independent square-root interval')
    return lo, hi


def choose_cut(Q,k2,periodic):
    """An exact rational ceiling of the stationary cut; optimality is unnecessary."""
    a,b,d = [F(periodic[k]) for k in
             ('quadratic_coefficient','linear_coefficient','constant_coefficient')]
    coefficient = 1-d/Q
    need(coefficient>0, 'positive cut quadratic coefficient')
    sl,sh = sqrt_interval(F(k2)/Q)
    dl,dh = sqrt_interval(b*b*k2/Q+4*coefficient*a*k2)
    lo,hi = (b*sl+dl)/(2*coefficient),(b*sh+dh)/(2*coefficient)
    scaled = hi*10**6
    z = F(-(-scaled.numerator//scaled.denominator),10**6)
    need(1<lo<=hi<=z and z-hi<F(1,10**6), 'outward rational cut selection')
    return z, dict(stationary_cut_lower=str(lo),stationary_cut_upper=str(hi),
                   rational_grid_denominator=10**6,exact_optimal_cut_claimed=False)


def independent_positive(row,periodic):
    Q,k2,z = F(row['base_min']),F(row['kappa_squared']),F(row['cut_z'])
    a,b,d = [F(periodic[k]) for k in
             ('quadratic_coefficient','linear_coefficient','constant_coefficient')]
    sl,sh = sqrt_interval(k2/Q); ll,lh = log(k2*Q/z**2); zl,zh = log(z)
    E = a/2+b
    lo = zl+a*k2/(2*z*z)+b*sl/z+d*ll/(2*Q)-E/Q
    hi = zh+a*k2/(2*z*z)+b*sh/z+d*lh/(2*Q)-E/Q
    need(F(row['mean_lower'])<=lo<=hi<=F(row['mean_upper']),
         'independent160-term positive integral inside primary enclosure')
    need(ll>1 and b*b*k2*Q>4*z*z*E*E,
         'positive integral decreases for every real M>=Q')
    return lo,hi


def independent_negative(profile,Q):
    """Integrate directly in t^2 using stored, fully replayed envelope segments."""
    scale=int(profile['tick_scale'])
    breaks=[(int(i),F(int(v),scale)) for i,v in profile['breakpoints_radius_tick']]
    need(breaks[0]==(0,0) and all(breaks[i][0]<breaks[i+1][0]
         and breaks[i][1]<=breaks[i+1][1] for i in range(len(breaks)-1)),
         'complete ordered nondecreasing envelope')
    radius=isqrt(Q//2)
    if 2*radius*radius==Q:radius-=1
    need(2*radius*radius<Q<=2*(radius+1)**2, 'strict negative-log boundary')
    rational=F(0);endpoint=F(0)
    for i,(start,value) in enumerate(breaks):
        stop=breaks[i+1][0] if i+1<len(breaks) else None
        left=F(start*start,Q)
        right=min(F(stop*stop,Q),F(1,2)) if stop is not None else F(1,2)
        if left<right:rational+=value*(right-left)/2
        if start<=radius and (stop is None or radius<stop):endpoint=value
    coefficient=endpoint/4
    need(rational>=0 and coefficient>0, 'nonnegative exact integration coefficients')
    l2,h2=log(2)
    return rational+coefficient*l2,rational+coefficient*h2,rational,coefficient


def classical_tail(periodic,gates):
    """The new coupled coefficient is intentionally NOT scaled in this tail."""
    omega=16;ps=OLD.primes_to(101);P=prod(ps[:omega]);nextp=ps[omega]
    H=2**(omega-1);k2=F(H*H+3,3)
    z,cut=choose_cut(P,k2,periodic)
    row=AR.integrated_cap(P,k2,z,periodic,gates,'tail16')
    independent_positive(row,periodic)
    l2,h2=AR.log_bounds(2);lp,hp=AR.log_bounds(P);lr,hr=AR.log_bounds(nextp)
    Cl,Ch=F(row['mean_lower']),F(row['mean_upper'])
    derivative=h2*hp-lr*Cl
    gates.require(derivative<0, 'infinite tail derivative for every real v>=0')
    il2,ih2=log(2);ipl,iph=log(P);irl,irh=log(nextp)
    need(ih2*iph-irl*Cl<0, 'independent all-real-v derivative sign')
    need(nextp==59 and ps[:omega][-1]==53 and H==32768,
         'sixteen-prime tail anchor and complete prime catalogue')
    row.update(from_omega=omega,v='omega(B)-16',primorial_primes=ps[:omega],
        next_prime_lower=nextp,H_at_anchor=H,base_growth_lower='59^v',
        kappa_growth_upper='sqrt(kappa_squared)*2^v',cut_growth='cut_z*2^v',
        ratio_kappa_over_cut_is_constant=True,
        function='(C+v*log(2))/(log(base_min)+v*log(next_prime_lower))',
        C_lower=str(Cl),C_upper=str(Ch),derivative_v_upper=str(derivative),
        monotone_for_every_real_v_ge_0=True,
        Kanold_universal_bound_is_a_declared_classical_import=True,
        coefficient_source='classical W, kappa_squared=(H^2+3)/3',
        improved_coupled_coefficient_scaled_from_anchor=False,
        finite_cardinality_scan_used_as_proof=False,cut_selection=cut)
    return row


def build_certificate(root=ROOT):
    root=Path(root);bindings=check_imports(root);gates=AR.Gates()
    OLD.check_manifest(root/'global/previous')
    previous_component=OLD.build_certificate(root/'global/previous')
    need(previous_component==read(root/'global/previous/results/global_certificate.json'),
         'complete preceding arithmetic component regenerated without changes')
    old=read(root/'imports/global_6090_165p.json')
    old_fixed=read(root/'imports/closure6090.json');old_family=read(root/'imports/family165.json')
    need(old['new_6090_premise_paid'] and old['new_165p_family_premise_paid'],
         'historical master paid both immediately preceding premises')
    need(old['new_6090_certificate_sha256']==sha(root/'imports/closure6090.json')
         and old['new_165p_certificate_sha256']==sha(root/'imports/family165.json'),
         'authenticated historical closure identities')
    need(old_fixed['status']==old_family['status']=='PASS'
         and old_fixed['radical']==6090 and old_family['core']==165
         and old_family['prime_minimum']==13,
         'historical closed domains')
    for obj,gamma in [(old_fixed,'Gamma6090'),(old_family,'Gamma165')]:
        need(obj['all_positive_valuations'] and obj['all_N_ge_1']
             and obj['completed_full_replay_verified'] and F(obj[gamma])<F(2,85),
             'historical all-length sector with paid raw replay')
    need(old['first_open_four_prime_radical']==2210
         and old['first_open_five_prime_radical']==6270, 'immediate preceding frontiers')
    fixed=read(root/'results/closure6270.json');family=read(root/'results/family130_certificate.json')
    need(fixed['status']==family['status']=='PASS' and fixed['radical']==6270
         and family['core']==130 and family['prime_minimum']==17, 'both new sector statements')
    for name,obj,gamma in [('6270',fixed,'Gamma6270'),('130p',family,'Gamma130')]:
        need(obj['all_positive_valuations'] and obj['all_N_ge_1']
             and obj['normalization']=='Delta(B,N)=mu(B^N)/ln(B)'
             and F(obj[gamma])<F(2,85) and obj['completed_full_replay_verified'],
             'complete new geometric premise '+name)
    need(family['completed_independent_physical_consumption_verified'],
         'all new family physical records independently consumed')
    families=[(r['core'],r['new_prime_minimum']) for r in old['declared_four_prime_families']]
    need(families==[(30,7),(42,11),(70,11),(66,13),(105,11),(78,17),(110,13),
                   (102,19),(154,13),(165,13)], 'all ten preceding four-prime families')
    families+=[(130,17)]
    closed5=old['five_prime_radicals_closed_below_2_over_85']+[6270]
    need(closed5==[2310,2730,3570,3990,4290,4830,5610,6006,6090,6270],
         'all ten closed five-prime supports')
    Q4,Q5=2622,6510;inventories={}
    for k,Q,divisor in [(4,Q4,30),(5,Q5,210)]:
        ps=OLD.primes_to(Q//divisor)
        tuples=sorted((prod(c),list(c)) for c in combinations(ps,k) if prod(c)<=Q)
        integers=[]
        for n in range(2,Q+1):
            fs=OLD.factor(n)
            if len(fs)==k and all(e==1 for p,e in fs):integers.append((n,[p for p,e in fs]))
        need(tuples==integers, 'two exhaustive support inventories at k='+str(k))
        rows=[]
        for radical,primes in tuples:
            if k==4:
                hits=[c for c,pmin in families if radical%c==0
                      and OLD.prime(radical//c) and radical//c>=pmin and c%(radical//c)!=0]
                covered=bool(hits)
            else:hits=[];covered=radical in closed5
            need(covered==(radical<Q), 'exact first uncovered support')
            rows.append(dict(radical=radical,prime_support=primes,covered=covered,
                             covering_family_core=hits[0] if hits else None))
        inventories[str(k)]=dict(frontier=Q,count=len(rows),rows=rows,
                                 largest_prime_enumerated=Q//divisor)
    need(inventories['4']['rows'][-1]['prime_support']==[2,3,19,23]
         and inventories['5']['rows'][-1]['prime_support']==[2,3,5,7,31],
         'new frontiers include every combination')
    cover_module=module(root/'covering/verify_covering.py','p72_current_covering')
    cover_module.check_manifest();cover=cover_module.build_certificate()
    need(canonical(cover)==(root/'covering/certificate.json').read_bytes(),
         'entire covering certificate reconstructed')
    Hs={r['prime_count_at_most']:r['H'] for r in cover['bounds']}
    need([Hs[k] for k in range(1,15)]==[2,4,6,10,14,22,26,44,58,66,98,148,204,360]
         and cover['case_count']==399 and cover['all_prime_supports']
         and cover['all_arbitrary_forbidden_classes'] and cover['all_integer_intervals']
         and not cover['primorial_extremality_assumed'], 'uniform covering scope through14 residual primes')
    manifest(root/'geometry');manifest(root/'negative')
    geometry=module(root/'geometry/verify_geometry.py','p72_current_geometry')
    gc=read(root/'geometry/certificate.json')
    need(geometry.kappa_squared(6)==10 and geometry.kappa_squared(10)==28
         and 4*geometry.kappa_squared(10)<geometry.kappa_squared(20),
         'coupled endpoints retained and invalid anchor scaling rejected')
    neg=module(root/'negative/uniform_negative.py','p72_current_negative')
    nc=read(root/'negative/certificate.json')
    need(nc['status']=='PASS_UNIFORM_NEGATIVE_PART_FOUR_FIVE_SIX_PRIMES'
         and nc['checks']==47759 and not nc['finite_examples_extrapolated'],
         'uniform negative proof scope, not finite examples')
    periodic=AR.period_certificate(gates)
    need(periodic==previous_component['periodic_CDF']==old['periodic_CDF'],
         'periodic CDF unchanged and completely recomputed')
    primes=OLD.primes_to(101);finite=[]
    for omega in range(4,16):
        Q=Q4 if omega==4 else Q5 if omega==5 else prod(primes[:omega])
        H=Hs[omega-1];k2=geometry.kappa_squared(H);z,cut=choose_cut(Q,k2,periodic)
        row=AR.integrated_cap(Q,k2,z,periodic,gates,'current'+str(omega))
        independent_positive(row,periodic)
        lower=F(0);negative_interval=None;coefficients=None
        if omega in (4,5,6):
            lower=F(neg.lower_bound(omega,Q));primary=tuple(map(F,neg.lower_interval(omega,Q)))
            ownlo,ownhi,A,B=independent_negative(nc['envelopes'][str(omega)],Q)
            declared=neg.decompose(omega,Q)
            need(A==declared['rational'] and B==declared['log2_coefficient']
                 and primary[0]==lower and 0<lower<=ownlo<=ownhi<=primary[1],
                 'independent exact bands and negative-log interval')
            negative_interval=list(map(str,primary))
            coefficients=dict(rational=str(A),log2_coefficient=str(B),
                              strict_radius=declared['strict_radius'])
        ql,qh=AR.log_bounds(Q);numerator=F(row['mean_upper'])-lower
        need(numerator>0, 'positive final numerator divided by lnB then lnQ')
        upper=numerator/ql
        row.update(omega=omega,H=H,coefficient_source='coupled-shear W',
            classical_kappa_squared=str(F(H*H+3,3)),cut_selection=cut,
            uniform_negative_lower=str(lower),negative_formula_interval=negative_interval,
            negative_exact_coefficients=coefficients,final_mean_upper=str(numerator),
            final_upper=str(upper),final_upper_decimal=out_decimal(upper),
            positive_part_only_upper=row['upper'],positive_part_only_upper_decimal=row['upper_decimal'])
        finite.append(row)
    tail=classical_tail(periodic,gates)
    candidate=read(root/'results/candidate210.json')
    need(candidate['status']=='PASS' and candidate['M']==210
         and candidate['integer_comparison_pass'] and F(2,85)<F(candidate['c0_lower']),
         'candidate and exact intermediate exclusion threshold')
    bounds=[(F(r['final_upper']),'omega='+str(r['omega'])) for r in finite]
    bounds += [(F(tail['upper']),'omega>=16'),(F(candidate['c0_upper']),'closed/candidate sectors')]
    cap,dominant=max(bounds);oldcap=F(old['global_upper'])
    need(cap<oldcap and dominant=='omega=5', 'strict global improvement with all sectors included')
    sources=['global/PROOF.md','global/verify_global.py','code/common.py','code/exact.py',
        'results/closure6270.json','results/family130_certificate.json','results/candidate210.json',
        'negative/PROOF.md','negative/uniform_negative.py','negative/certificate.json',
        'geometry/PROOF.md','geometry/verify_geometry.py','geometry/certificate.json',
        'covering/PROOF.md','covering/verify_covering.py','covering/certificate.json',
        'imports/global_6090_165p.json','global/previous/results/global_certificate.json']
    return dict(status='PASS_WITH_DECLARED_HISTORICAL_DEPENDENCIES',
        normalization='Delta(B,N)=mu(B^N)/ln(B)',all_positive_valuations=True,all_N_ge_1=True,
        all_states_including_zero=True,first_open_four_prime_radical=Q4,
        first_open_five_prime_radical=Q5,radical_inventories=inventories,
        declared_four_prime_families=[dict(core=c,new_prime_minimum=p) for c,p in families],
        five_prime_radicals_closed_below_2_over_85=closed5,
        new_6270_premise_paid=True,new_130p_family_premise_paid=True,
        new_6270_certificate_sha256=sha(root/'results/closure6270.json'),
        new_130p_certificate_sha256=sha(root/'results/family130_certificate.json'),
        historical_30p_includes_candidate_and_is_not_below_2_over_85=True,
        uniform_negative_used_at_prime_counts=[4,5,6],periodic_CDF=periodic,
        finite_sectors=finite,infinite_tail=tail,dominant_sector=dominant,
        global_upper=str(cap),global_upper_decimal=out_decimal(cap),
        previous_global_upper=str(oldcap),previous_global_upper_decimal=out_decimal(oldcap),
        exact_improvement_lower=str(oldcap-cap),global_candidate_lower=candidate['c0_lower'],
        candidate=candidate,new_uniform_covering=cover['bounds'],covering_cases=cover['case_count'],
        coupled_shear_bound_used_at_prime_counts=list(range(4,16)),
        negative_examples_used_as_uniform_constants=False,
        continuous_M_and_prime_count_monotonicity_checked=True,
        independent_log_series_terms=160,independent_sqrt_decimal_places=40,
        historical_large_raw_replayed_this_task=False,
        new_external_asymptotic_theorem_used_as_numerical_premise=False,
        inherited_42pr_regions=old['inherited_42pr_regions'],full_42pr_family_closed=False,
        inherited_far_tail=old['inherited_far_tail'],
        remaining_domain_for_a_value_greater_than_c0=old['remaining_domain_for_a_value_greater_than_c0'],
        all_four_prime_supports_closed=False,all_five_prime_supports_closed=False,
        exact_Cstar_equals_c0_proved=False,offset_sign_classification_claimed=False,
        internal_AI_review_is_not_specialist_human_review=True,
        primary_arithmetic_assertions=gates.count,authenticated_imports=len(bindings),
        inputs=[dict(path=n,sha256=sha(root/n)) for n in sources])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true');parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    if args.report:need(not args.report.resolve().is_relative_to(ROOT), 'external report required')
    started=time.monotonic();result=build_certificate();path=ROOT/'global/certificate.json'
    if args.build:path.write_bytes(canonical(result))
    else:need(path.read_bytes()==canonical(result), 'entire global certificate regenerated identically')
    report=dict(status='PASS',certificate_sha256=sha(path),global_upper=result['global_upper'],
        global_upper_decimal=result['global_upper_decimal'],
        first_open_four_prime_radical=result['first_open_four_prime_radical'],
        first_open_five_prime_radical=result['first_open_five_prime_radical'],
        runtime_seconds=round(time.monotonic()-started,3))
    if args.report:args.report.write_bytes(canonical(report))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
