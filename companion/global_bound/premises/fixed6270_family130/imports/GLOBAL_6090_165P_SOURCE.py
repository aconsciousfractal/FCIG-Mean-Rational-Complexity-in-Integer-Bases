#!/usr/bin/env python3
"""All prime-count sectors with paid 6090/165p and uniform negative log."""
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
from common import need, sha, read, module, canonical, check_imports
from exact import alternate_log, out_decimal

OLD = module(ROOT/'global/previous/verify_global.py', 'p72_previous_global')
AR = OLD.AR


def log(x):
    return alternate_log(F(x), terms=160)


def sqrt_interval(x):
    x = F(x)
    scale = 10**40
    n = isqrt((x*scale*scale).numerator//(x*scale*scale).denominator)
    lo, hi = F(n,scale), F(n+1,scale)
    need(lo*lo <= x < hi*hi, 'independent square-root interval')
    return lo, hi


def independent_positive(row, periodic):
    Q,k2,z = F(row['base_min']),F(row['kappa_squared']),F(row['cut_z'])
    a,b,d = [F(periodic[k]) for k in
             ('quadratic_coefficient','linear_coefficient','constant_coefficient')]
    sl,sh = sqrt_interval(k2/Q)
    ll,lh = log(k2*Q/z**2)
    zl,zh = log(z)
    E = a/2+b
    lo = zl+a*k2/(2*z*z)+b*sl/z+d*ll/(2*Q)-E/Q
    hi = zh+a*k2/(2*z*z)+b*sh/z+d*lh/(2*Q)-E/Q
    need(F(row['mean_lower']) <= lo <= hi <= F(row['mean_upper']),
         'independent160-term positive integral inside primary enclosure')
    need(ll>1 and b*b*k2*Q>4*z*z*E*E,
         'positive integral decreases for every real M>=Q')
    return lo,hi


def independent_negative(omega, Q):
    """Independent evaluation of the proved continuous minorant, not its proof."""
    c = {5:F(4,25),6:F(3,25)}[omega]
    Q = F(Q)
    need(Q>=882, 'declared uniform negative-log domain')
    sl,sh = sqrt_interval(F(2)/Q)
    l2,h2 = log(2); lq,hq = log(Q); lr,hr = log(20)
    # Regroup the log2 coefficient, which is strictly positive in this domain.
    need(F(1,4)-sh/2 > 0, 'positive log2 coefficient')
    lo = 2*c*(F(1,4)-sh+(F(1,4)-sh/2)*l2
              +(F(-160)+lq/2-hr)/Q)
    hi = 2*c*(F(1,4)-sl+(F(1,4)-sl/2)*h2
              +(F(-160)+hq/2-lr)/Q)
    need(0<lo<=hi, 'strictly positive uniform subtraction')
    return lo,hi


def build_certificate(root=ROOT):
    root = Path(root)
    bindings = check_imports(root)
    gates = AR.Gates()
    OLD.check_manifest(root/'global/previous')
    previous_component = OLD.build_certificate(root/'global/previous')
    need(previous_component == read(root/'global/previous/results/global_certificate.json'),
         'complete preceding arithmetic component regenerated without changes')
    old = read(root/'imports/coverage_6006_154p.json')
    old_fixed = read(root/'imports/closure6006.json')
    old_family = read(root/'imports/family154.json')
    need(old['new_6006_premise_paid'] and old['new_154p_family_premise_paid'],
         'historical master paid both preceding premises')
    need(old['new_fixed_certificate_sha256']==sha(root/'imports/closure6006.json')
         and old['new_family_certificate_sha256']==sha(root/'imports/family154.json'),
         'historical closed-sector certificates have exact identities')
    need(old_fixed['status']==old_family['status']=='PASS'
         and old_fixed['radical']==6006 and old_family['core']==154
         and old_family['prime_minimum']==13
         and old_fixed['all_positive_valuations'] and old_fixed['all_N_ge_1']
         and old_family['all_positive_valuations'] and old_family['all_N_ge_1']
         and F(old_fixed['Gamma6006'])<F(2,85) and F(old_family['Gamma154'])<F(2,85),
         'historical sector domains and strict margins')
    need(old['first_open_four_prime_radical']==2145
         and old['first_five_prime_radical_not_closed_vs_c0']==6090,
         'immediate preceding frontiers')
    fixed = read(root/'results/closure6090.json')
    family = read(root/'results/family165_certificate.json')
    need(fixed['status']==family['status']=='PASS' and fixed['radical']==6090
         and family['core']==165 and family['prime_minimum']==13,
         'both exact new sector statements')
    for name,obj,gamma in [('6090',fixed,'Gamma6090'),('165p',family,'Gamma165')]:
        need(obj['all_positive_valuations'] and obj['all_N_ge_1']
             and obj['normalization']=='Delta(B,N)=mu(B^N)/ln(B)'
             and F(obj[gamma])<F(2,85) and obj['completed_full_replay_verified'],
             'complete new geometric premise '+name)
    need(family['completed_independent_physical_consumption_verified'],
         'new odd-family physical records independently consumed')
    families = [(r['core'],r['new_prime_minimum']) for r in old['declared_four_prime_families']]
    need(families==[(30,7),(42,11),(70,11),(66,13),(105,11),(78,17),(110,13),(102,19),(154,13)],
         'all nine preceding four-prime families')
    families += [(165,13)]
    closed5 = old['five_prime_radicals_closed_below_2_over_85']+[6090]
    need(closed5==[2310,2730,3570,3990,4290,4830,5610,6006,6090],
         'all nine closed five-prime supports')
    Q4,Q5=2210,6270
    inventories={}
    for k,Q,divisor in [(4,Q4,30),(5,Q5,210)]:
        ps=OLD.primes_to(Q//divisor)
        tuples=sorted((prod(c),list(c)) for c in combinations(ps,k) if prod(c)<=Q)
        integers=[]
        for n in range(2,Q+1):
            fs=OLD.factor(n)
            if len(fs)==k and all(e==1 for p,e in fs):
                integers.append((n,[p for p,e in fs]))
        need(tuples==integers, 'two exhaustive support inventories at k='+str(k))
        rows=[]
        for radical,primes in tuples:
            if k==4:
                hits=[c for c,pmin in families if radical%c==0
                      and OLD.prime(radical//c) and radical//c>=pmin
                      and c%(radical//c)!=0]
                covered=bool(hits)
            else:
                hits=[]; covered=radical in closed5
            need(covered==(radical<Q), 'exact first uncovered support')
            rows.append(dict(radical=radical,prime_support=primes,covered=covered,
                             covering_family_core=hits[0] if hits else None))
        inventories[str(k)]=dict(frontier=Q,count=len(rows),rows=rows,
                                  largest_prime_enumerated=Q//divisor)
    need(inventories['4']['rows'][-1]['prime_support']==[2,5,13,17]
         and inventories['5']['rows'][-1]['prime_support']==[2,3,5,11,19],
         'new frontiers include all combinations')

    cover_module=module(root/'covering/verify_covering.py','p72_current_covering')
    cover_module.check_manifest()
    cover=cover_module.build_certificate()
    need(canonical(cover)==(root/'covering/certificate.json').read_bytes(),
         'complete new covering certificate reconstructed')
    Hs={r['prime_count_at_most']:r['H'] for r in cover['bounds']}
    need(cover['all_prime_supports'] and cover['all_arbitrary_forbidden_classes']
         and cover['all_integer_intervals'] and not cover['primorial_extremality_assumed'],
         'covering has exactly the uniform scope needed by W')
    neg_module=module(root/'negative/uniform_negative.py','p72_current_negative')
    periodic=AR.period_certificate(gates)
    need(periodic==previous_component['periodic_CDF'], 'periodic CDF unchanged and recomputed')
    primes=OLD.primes_to(47)
    zs={4:F(4),5:F(20,3),6:F(37,4),7:F(29,2),8:F(17),
        9:F(29),10:F(38),11:F(43),12:F(64),13:F(97)}
    finite=[]
    for omega,z in zs.items():
        Q=Q4 if omega==4 else Q5 if omega==5 else prod(primes[:omega])
        H=Hs[omega-1]
        row=AR.integrated_cap(Q,F(H*H+3,3),z,periodic,gates,'current'+str(omega))
        independent_positive(row,periodic)
        lower=F(0); negative_interval=None
        if omega in (5,6):
            lower=F(neg_module.lower_bound(omega,Q))
            primary_interval=tuple(map(F,neg_module.lower_interval(omega,Q)))
            ownlo,ownhi=independent_negative(omega,Q)
            need(primary_interval[0] == lower and 0<lower<=ownlo<=ownhi<=primary_interval[1],
                 'independent uniform negative-log integration')
            negative_interval=list(map(str,primary_interval))
        ql,qh=AR.log_bounds(Q)
        numerator=F(row['mean_upper'])-lower
        need(numerator>0, 'positive final numerator divided by lnB then lnQ')
        upper=numerator/ql
        row.update(omega=omega,H=H,uniform_negative_lower=str(lower),
                   negative_formula_interval=negative_interval,
                   final_mean_upper=str(numerator),final_upper=str(upper),
                   final_upper_decimal=out_decimal(upper),
                   positive_part_only_upper=row['upper'],
                   positive_part_only_upper_decimal=row['upper_decimal'])
        finite.append(row)
    tail=AR.monotone_tail(14,5358,periodic,gates)
    independent_positive(tail,periodic)
    need(tail['primorial_primes']==primes[:14] and tail['next_prime_lower']==47
         and tail['H_at_anchor']==8192 and tail['monotone_for_every_real_v_ge_0'],
         'complete infinite prime-count tail starts at14')
    l2,h2=log(2); lp,hp=log(tail['base_min']); lr,hr=log(47)
    need(h2*hp-lr*F(tail['C_lower'])<0, 'independent all-real-v derivative sign')
    candidate=read(root/'results/candidate210.json')
    need(candidate['status']=='PASS' and candidate['M']==210
         and candidate['integer_comparison_pass'] and F(2,85)<F(candidate['c0_lower']),
         'candidate and intermediate exclusion threshold')
    bounds=[(F(r['final_upper']),'omega='+str(r['omega'])) for r in finite]
    bounds += [(F(tail['upper']),'omega>=14'),(F(candidate['c0_upper']),'closed/candidate sectors')]
    cap,dominant=max(bounds)
    need(cap<F(old['unrounded_integration_upper']), 'strict improvement over preceding global cap')
    need(max(F(r['final_upper']) for r in finite if r['omega']>=6)<=cap
         and F(tail['upper'])<=cap, 'every finite and infinite prime-count sector included')
    oldcap=F(old['unrounded_integration_upper'])
    onlyfront=max([F(r['upper']) for r in finite if r['omega']<=8]
                  +[F(previous_component['infinite_tail']['upper'])])
    sources=['global/PROOF.md','global/verify_global.py','code/common.py','code/exact.py',
             'results/closure6090.json','results/family165_certificate.json','results/candidate210.json',
             'negative/PROOF.md','negative/uniform_negative.py','negative/certificate.json',
             'covering/PROOF.md','covering/verify_covering.py','covering/certificate.json',
             'imports/coverage_6006_154p.json','global/previous/results/global_certificate.json']
    inputs=[dict(path=n,sha256=sha(root/n)) for n in sources]
    return dict(status='PASS_WITH_DECLARED_HISTORICAL_DEPENDENCIES',
                normalization='Delta(B,N)=mu(B^N)/ln(B)',all_positive_valuations=True,all_N_ge_1=True,
                all_states_including_zero=True,first_open_four_prime_radical=Q4,
                first_open_five_prime_radical=Q5,radical_inventories=inventories,
                declared_four_prime_families=[dict(core=c,new_prime_minimum=p) for c,p in families],
                five_prime_radicals_closed_below_2_over_85=closed5,
                new_6090_premise_paid=True,new_165p_family_premise_paid=True,
                new_6090_certificate_sha256=sha(root/'results/closure6090.json'),
                new_165p_certificate_sha256=sha(root/'results/family165_certificate.json'),
                historical_30p_includes_candidate_and_is_not_below_2_over_85=True,
                uniform_negative_used_at_prime_counts=[5,6],periodic_CDF=periodic,
                finite_sectors=finite,infinite_tail=tail,dominant_sector=dominant,
                global_upper=str(cap),global_upper_decimal=out_decimal(cap),
                previous_global_upper=str(oldcap),previous_global_upper_decimal=out_decimal(oldcap),
                exact_improvement_lower=str(oldcap-cap),
                frontier_only_with_old_tail_upper=str(onlyfront),
                frontier_only_with_old_tail_upper_decimal=out_decimal(onlyfront),
                global_candidate_lower=candidate['c0_lower'],candidate=candidate,
                new_uniform_covering=cover['bounds'],covering_cases=cover['case_count'],
                negative_examples_used_as_uniform_constants=False,
                continuous_M_and_prime_count_monotonicity_checked=True,
                independent_log_series_terms=160,independent_sqrt_decimal_places=40,
                historical_large_raw_replayed_this_task=False,
                new_external_asymptotic_theorem_used_as_numerical_premise=False,
                inherited_42pr_regions=old['new_42pr_infinite_regions'],full_42pr_family_closed=False,
                inherited_far_tail=old['inherited_far_tail'],
                remaining_domain_for_a_value_greater_than_c0=old['remaining_domain_for_a_value_greater_than_c0'],
                all_four_prime_supports_closed=False,all_five_prime_supports_closed=False,
                exact_Cstar_equals_c0_proved=False,offset_sign_classification_claimed=False,
                internal_AI_review_is_not_specialist_human_review=True,
                primary_arithmetic_assertions=gates.count,authenticated_imports=len(bindings),inputs=inputs)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    if args.report:
        need(not args.report.resolve().is_relative_to(ROOT), 'external report required')
    started=time.monotonic()
    result=build_certificate()
    path=ROOT/'global/certificate.json'
    if args.build:
        path.write_bytes(canonical(result))
    else:
        need(path.read_bytes()==canonical(result), 'entire global certificate regenerated identically')
    report=dict(status='PASS',certificate_sha256=sha(path),
                global_upper=result['global_upper'],global_upper_decimal=result['global_upper_decimal'],
                first_open_four_prime_radical=result['first_open_four_prime_radical'],
                first_open_five_prime_radical=result['first_open_five_prime_radical'],
                runtime_seconds=round(time.monotonic()-started,3))
    if args.report:args.report.write_bytes(canonical(report))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
