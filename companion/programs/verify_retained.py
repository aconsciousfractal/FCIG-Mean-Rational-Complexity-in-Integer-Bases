"""Portable retained-evidence verification, not replay of historical science.

The central intermediate arrays remain a producer/proof dependency.
E032 has no retained histogram: its recorded Q is source-bound, not rebuilt.
"""
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import gcd
from archive import PACKAGE, data, entry, manifest, read_object

S = 2**64
W = 2**192


def need(ok, message):
    if not ok:
        raise ValueError(message)


def domain_rows(rows, expected, keys, fixed_primes=None):
    """Bind certificate obligations to a distinct, ordered exact domain."""
    actual = []
    for row in rows:
        need(type(row) is dict and all(key in row for key in keys), 'complete domain profile fields')
        values = [row[key] for key in keys]
        if fixed_primes is not None:
            values = [values[0], *fixed_primes, *values[1:]]
        need(all(type(value) is int for value in values), 'exact integer domain profile')
        actual.append(values)
    need(all(type(profile) is list and all(type(value) is int for value in profile)
             for profile in expected), 'exact independent domain profiles')
    need(len({tuple(profile) for profile in actual}) == len(actual)
         and len({tuple(profile) for profile in expected}) == len(expected)
         and actual == expected, 'complete ordered domain identity')


def hist_check(hist, M):
    need(type(hist) is list and all(type(v) is list and len(v)==2
        and all(type(x) is int for x in v) and v[0]>=1 and v[1]>0 for v in hist), 'histogram shape')
    need([h for h,c in hist]==sorted({h for h,c in hist}), 'ordered unique heights')
    need(sum(c for h,c in hist)==M, 'original histogram mass')


def series32(a,b):
    if a==0:return 0,0
    lower=0;numerator,denominator=a,b
    for j in range(32):
        lower+=2*numerator*W//((2*j+1)*denominator)
        numerator*=a*a;denominator*=b*b
    rn=2*numerator*b*b*W;rd=65*denominator*(b*b-a*a)
    return lower,lower+32+-((-rn)//rd)


def log32(n):
    need(type(n) is int and n>0,'positive log integer')
    k=n.bit_length()-1;p=2**k
    L,U=series32(n-p,n+p);l2,u2=series32(1,3)
    return (L+k*l2)*S//W,-((-((U+k*u2)*S))//W)


@lru_cache(maxsize=None)
def recipe16(n):
    k=n.bit_length()-1;p=2**k
    def series(t):
        total=Fraction(0);power=t
        for j in range(16):
            total+=2*power/(2*j+1);power*=t*t
        return total,total+2*power/(33*(1-t*t))
    L,U=series(Fraction(n-p,n+p));l2,u2=series(Fraction(1,3))
    L=(L+k*l2)*S;U=(U+k*u2)*S
    return L.numerator//L.denominator,-((-U.numerator)//U.denominator)


def early(round_number,name,count):
    payload=data(round_number,name);rows=payload['rows'];need(len(rows)==count,'early complete domain count')
    if round_number == 26:
        domain_rows(rows, [[n, 6**n] for n in range(1, 9)], ('N', 'M'))
    elif round_number == 27:
        domain_rows(rows, [[b, 6, b**6, rad] for b, rad in ((10,10),(12,6),(14,14))],
                    ('base', 'N', 'M', 'radical'))
    else:
        raise ValueError('unsupported early domain')
    for row in rows:
        M=row['M'];hist_check(row['histogram'],M)
        den=row['upper_centered_log_denominator'];factor=den//(M*S)
        need(factor in (1,2) and den==factor*M*S and (factor*M)%2==0,'early exact centered denominator')
        q=factor*sum(c*recipe16(h)[1] for h,c in row['histogram'])-(factor*M//2)*recipe16(M)[0]
        need(q==row['upper_centered_log_numerator'] and q<0 and row['strict_negative'] is True,'early outward Q reconstruction')
    return dict(rows=len(rows),mass=sum(r['M'] for r in rows),log_recipe='exact K16 reconstruction from histogram')


def aggregate(round_number,name,count):
    payload=data(round_number,name);rows=payload['rows'];need(len(rows)==count,'complete aggregate count')
    domain_name={29:'RADICAL6_DOMAIN.json',30:'RADICAL10_DOMAIN.json',31:'PARAMETRIC_DOMAIN.json',32:'LENGTH_FOUR_DOMAIN.json',33:'LENGTH_THREE_DOMAIN.json',35:'LENGTH_TWO_DOMAIN.json',39:'BRANCH_EXECUTION_DOMAIN.json',40:'BRANCH_DOMAIN.json',41:'BRANCH_DOMAIN.json'}[round_number]
    need(payload['domain_sha256']==entry(round_number,domain_name)['sha256'],'aggregate original domain byte binding')
    domain = data(round_number, domain_name)
    expected = domain['profiles' if round_number in (39, 40, 41) else 'rows']
    need(len(expected) == count, 'independent aggregate domain count')
    if round_number in (29, 30):
        domain_rows(rows, expected, ('M', 'a', 'b'), (2, 3 if round_number == 29 else 5))
    else:
        domain_rows(rows, expected, ('M', 'p', 'r', 'a', 'b'))
    ledger={n:(L,U) for n,L,U in payload['log_tick_ledger']}
    need(len(ledger)==len(payload['log_tick_ledger']) and list(ledger)==sorted(ledger),'unique sorted log ledger')
    for n,(L,U) in ledger.items():
        l,u=log32(n);need(L<=l<=u<=U,'distinct K32 log containment')
    used=set();paid=[]
    for r in rows:
        M=r['M'];hist=r['upper_histogram'];hist_check(hist,M)
        need(all(h in ledger for h,c in hist) and M in ledger,'log closure')
        q=2*sum(c*ledger[h][1] for h,c in hist)-M*ledger[M][0]
        need(q==r['upper_centered_numerator'] and r['upper_centered_denominator']==2*M*S and q<0 and r['status']=='PAID','exact outward Q')
        used.add(M);used.update(h for h,c in hist);paid.append(M)
    # Earlier producers share a ledger with their retained small controls.
    # These controls are not additional main obligations or mandatory PAID rows.
    validation=payload.get('validation',[])
    for r in validation:
        M=r['M'];hist=r['upper_histogram'];hist_check(hist,M)
        need(M in ledger and all(h in ledger for h,c in hist),'retained validation log closure')
        q=2*sum(c*ledger[h][1] for h,c in hist)-M*ledger[M][0]
        need(q==r['upper_centered_numerator'] and r['upper_centered_denominator']==2*M*S and r['status']==('PAID' if q<0 else 'UNPAID'),'retained validation outward Q/status')
        used.add(M);used.update(h for h,c in hist)
    need(set(ledger)==used and payload['paid_moduli']==paid and payload['unpaid_moduli']==[],'complete strict sign partition')
    return dict(rows=len(rows),mass=sum(r['M'] for r in rows),consumed_log_enclosures=len(ledger),retained_validation_Q_checks=len(validation),all_intermediate_central_values_independently_recounted=False)


def prime_power():
    ternary = data(23,'FINITE_BRIDGE.json')['rows']
    general = data(23,'GENERAL_PRIME_BRIDGE.json')['rows']
    domain_rows(ternary, [[n, 3**n] for n in range(1, 9)], ('N', 'M'))
    cases = ((5,3),(5,4),(5,5),(5,6),(7,3),(7,4),(7,5),(11,3),(11,4),
             (13,3),(17,3),(19,3),(23,3))
    domain_rows(general, [[p, n, p**n] for p, n in cases], ('p', 'N', 'M'))
    rows=ternary+general
    need(len(rows)==21,'complete prime-power bridge')
    for r in rows:
        hist_check(r['exact_histogram'],r['M'])
        h=sum(t*c for t,c in r['exact_histogram'])
        need(h==r['height_sum'] and r['M']**3-h*h==r['integer_Jensen_margin']>0,'integer Jensen certificate')
    return dict(rows=21,mass=sum(r['M'] for r in rows),cdf_formulas_replayed=False)


def analytic_projection(round_number,name,total,paid,left):
    x=data(round_number,name);need(len(x['rows'])==total and x['analytic_paid']==paid and len(x['remaining_profiles'])==left,'analytical complete partition cardinality')
    need(len(x['analytic_paid_profiles'])==paid,'analytical paid-list count')
    yes=[];no=[]
    for r in x['rows']:
        profile=r['profile'] if round_number==34 else [r[k] for k in ('M','p','q','a','b')]
        value=r['mu_upper'] if round_number==34 else r['upper']
        need(type(value) is list and len(value)==2 and all(type(v) is int for v in value) and value[1]>0,'retained rational analytical upper')
        status='ANALYTIC_PAID' if value[0]<0 else 'NEEDS_FINITE'
        need(r['status']==status,'strict recorded analytical classification')
        (yes if value[0]<0 else no).append(profile)
    need(yes==x['analytic_paid_profiles'] and no==x['remaining_profiles'],'complete recorded analytical profile projection')
    return dict(rows=total,analytic_paid=paid,remaining=left,analytical_bounds_recomputed=False)


def main():
    m=manifest();seen=set()
    for e in m['entries']:
        need(e['locator'] not in seen,'unique provenance locator');seen.add(e['locator'])
        raw=read_object(e)
        need(len(raw)==e['bytes'] and hashlib.sha256(raw).hexdigest()==e['sha256'],'complete original byte object')
    results={'prime_power':prime_power(),'radical6_diagonal':early(26,'FINITE_BRIDGE.json',8),'sixth_length':early(27,'SIXTH_LENGTH_BRIDGE.json',3)}
    for n,name,count in ((29,'RADICAL6_AGGREGATE.json',143),(30,'RADICAL10_AGGREGATE.json',97),(31,'PARAMETRIC_AGGREGATE.json',132),(32,'LENGTH_FOUR_AGGREGATE.json',3),(33,'LENGTH_THREE_AGGREGATE.json',35),(35,'LENGTH_TWO_AGGREGATE.json',99),(39,'BRANCH_AGGREGATE.json',1087),(40,'BRANCH_AGGREGATE.json',1422),(41,'BRANCH_AGGREGATE.json',1901)):
        results['aggregate_'+str(n)]=aggregate(n,name,count)
    square=data(28,'SEMIPRIME_BRIDGE.json');need(len(square['rows'])==4596,'squarefree complete bridge')
    need(square['domain_sha256']==entry(28,'SEMIPRIME_DOMAIN.json')['sha256'],'squarefree domain byte binding')
    need([[r[k] for k in ('M','p','r')] for r in square['rows']]==[[r[k] for k in ('M','p','r')] for r in data(28,'SEMIPRIME_DOMAIN.json')['rows']],'squarefree complete original model order')
    need(all(r['upper_histogram_mass']==r['M'] and r['upper_centered_log_numerator']<0 and r['strict_negative'] is True and 'histogram' not in r for r in square['rows']),'squarefree recorded source-bound signs')
    results['squarefree']=dict(rows=4596,histogram_and_log_ledger_retained=False,Q_rebuilt_from_retained_data=False,reconstruction_command='python programs/reconstruct_semiprime.py --all')
    results['analytic_N2']=analytic_projection(34,'ANALYTIC_N2.json',454,355,99)
    results['analytic_N1']=analytic_projection(37,'ANALYTIC_N1.json',11382,6972,4410)
    generator=data(36,'inventory_generate_normal1.json');factor=data(36,'inventory_factor_normal1.json')
    for key in ('all_rows','retained_rows','active_rows'):
        need(generator[key]==factor[key],'two complete original metadata inventories')
    need(len(generator['retained_rows'])==19239 and len(generator['active_rows'])==11382,'complete conservative and ratio-filtered inventory')
    analytical=data(37,'ANALYTIC_N1.json')
    need([[r[k] for k in ('M','p','q','a','b')] for r in generator['active_rows']]==[[r[k] for k in ('M','p','q','a','b')] for r in analytical['rows']],'complete analytical input order')
    del generator,factor,analytical
    results['inventories']=dict(conservative=19239,active=11382,distinct_retained_algorithms_compared=True,inventory_algorithms_rerun=False)
    remaining=data(37,'ANALYTIC_N1.json')['remaining_profiles'];combined=[]
    for n in (39,40,41):
        name='BRANCH_EXECUTION_DOMAIN.json' if n==39 else 'BRANCH_DOMAIN.json'
        domain=data(n,name);payload=data(n,'BRANCH_AGGREGATE.json')
        need(payload['domain_sha256']==entry(n,name)['sha256'],'actual final domain digest')
        profiles=domain['profiles'];need(profiles==[[r[k] for k in ('M','p','r','a','b')] for r in payload['rows']],'complete final profile order')
        combined.extend(profiles)
    need(sorted(combined)==remaining and len(combined)==4410,'complete final 4410 union')
    need(data(39,'REMAINING_N1.json')['remaining_profiles']==[r for r in remaining if r[1]<5],'p2/p3 projection')
    need(data(40,'REMAINING_N1.json')['remaining_profiles']==[r for r in remaining if r[1]==2],'p2 projection')
    need(data(41,'REMAINING_N1.json')['remaining_profiles']==[],'empty two-prime obligations')
    physical={}
    for n,count,types in ((39,2361,6068480),(40,2505,13888896),(41,2600,5116390)):
        x=data(n,'verify_local_controls_normal1.json');records=x['complete_curve_records']
        need(len(records)==x['curves_checked']==count and all(r['status']=='PASS' for r in records),'complete retained physical policy')
        need(sum(r['signed_types_checked'] for r in records)==x['signed_types_checked']==types,'retained signed-type totals')
        need(sum(r['step_tests'] for r in records)==x['test_steps'] and sum(r['marginal_comparisons'] for r in records)==x['marginal_comparisons']==2*x['test_steps'],'retained physical counters')
        if n==41:
            singleton=[r for r in records if r['profile']['kind']=='m2']
            need(len(singleton)==x['m2_singletons']==650 and all(r['signed_types_checked']==1 and r['step_tests']==3 and r['missing_prime_repairs']==0 for r in singleton),'undoubled physical singleton records')
        physical[str(n)]=dict(curves=count,signed_types=types,physical_types_reexecuted=False)
    results['physical_records']=physical
    counter=next(r for r in data(24,'mixed_normal.json','reports')['rows'] if r['base']==30 and r['N']==2)
    hist_check(counter['histogram'],900);product=1
    for h,c in counter['histogram']:product*=h**c
    need(product*product>900**900,'exact B30 N2 product counterexample')
    print(json.dumps(dict(status='PORTABLE_RETAINED_ARITHMETIC_PASS',objects_checked=len({e['sha256'] for e in m['entries']}),locators_checked=len(m['entries']),stages=results,final_bridge_rows=4410,two_prime_remaining=0,universal_zero_offset=False,central_producer_and_proofs_required=True,large_minima_or_CDFs_replayed=False,historical_jobs_replayed=False),sort_keys=True,separators=(',',':')))


if __name__=='__main__':main()
