"""Exact assembly for all selected full marks and the two marked-core caps."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json
from exact import need, sieve, log_lo, out_decimal, prime_trial
from assemble import profile
from core import CORE, CORE_PRIMES, MIN_NEW_PRIME, mark
ROOT = Path(__file__).resolve().parents[1]

def assemble(root=ROOT):
    contract = json.loads((root / 'profile_contract.json').read_text())
    plan = json.loads((root / 'profiles_plan.json').read_text())
    need(contract['core'] == CORE and tuple(contract['core_primes']) == CORE_PRIMES, 'declared core factors')
    need(contract['first_new_prime'] == MIN_NEW_PRIME, 'first new prime')
    last = contract['last_small_prime']
    tail = contract['tail_first_prime']
    need(prime_trial(last) and prime_trial(tail) and (last >= MIN_NEW_PRIME) and (tail > last), 'prime cutoffs')
    need([p for p in sieve(tail) if last < p] == [tail], 'no omitted prime at the tail boundary')
    primes = [p for p in sieve(last) if p >= MIN_NEW_PRIME]
    need(sorted((z['p'] for z in plan if z['p'] > 1)) == primes, 'every declared small prime exactly once')
    need(sorted(((z['p'], z['k']) for z in plan if z['p'] == 1)) == [(1, 1), (1, 2)], 'both marked-core divisibilities')
    rows = []
    routes = []
    marks = []
    for z in plan:
        p, k, ag = (z['p'], z['k'], z['ag'])
        q, _ = mark(p)
        need(z['stem'] == (f'q{q}' if p > 1 else f'q{CORE}_k{k}'), 'unambiguous result stems')
        got = profile(p, root / 'results' / (z['stem'] + '.cells'), root, k, ag)
        stored = json.loads((root / 'results' / (z['stem'] + '.json')).read_text())
        need(set(stored) == set(got) | {'witness_bytes', 'witness_sha256'}, 'profile schema')
        need(all((stored[x] == v for x, v in got.items())), 'exact reintegration of every stored field')
        need(stored['witness_bytes'] == 8 * got['witnesses'] and len(stored['witness_sha256']) == 64, 'witness metadata only')
        cap = F(got['total_upper'])
        rows.append(got)
        if p == 1:
            simple = F(z['simple_cap'])
            need(simple == F(contract['marked_simple_caps'][str(k)]), 'declared core input')
            need(cap < simple, 'strict marked-core simple cap')
            marks.append({'q': CORE, 'k': k, 'divisor': CORE ** k, 'exact_upper': str(cap), 'upper_decimal': out_decimal(cap), 'simple_cap': str(simple), 'valid_for': f'nu_{CORE}(M) for every positive M divisible by divisor'})
        else:
            need(k == 1 and ag in (4, 8, 16), 'full mark one-step configuration')
            norm = max(F(0), cap) / log_lo(q)
            need(norm < F(2, 85), 'strict all-exponents and all-lengths route')
            routes.append({'p': p, 'q': q, 'domain': 'B=2^a 3^b 5^c 7^d p^e; a,b,c,d,e,N>=1', 'required_divisor_of_M': q, 'base_min': q, 'mean_upper': str(cap), 'normalized_upper': str(norm), 'normalized_upper_decimal': out_decimal(norm)})
    gamma = max((F(z['normalized_upper']) for z in routes))
    return {'status': 'PASS', 'closed_primes': primes, 'closed_radicals': [CORE * p for p in primes], 'all_positive_prime_exponents': True, 'all_lengths_N_ge_1': True, f'marked{CORE}_caps': marks, 'routes': routes, 'common_upper': str(gamma), 'common_upper_decimal': out_decimal(gamma), 'comparison_threshold': '2/85', 'strict_margin': str(F(2, 85) - gamma), 'direct_profiles': len(rows), 'geometric_witnesses_in_direct_profiles': sum((z['witnesses'] for z in rows)), 'corner_evaluations_in_direct_profiles': sum((4 * z['angular_nodes'] * z['witnesses'] for z in rows)), 'historical_profile_caps_assumed': False, f'entire_{CORE}p_family_closed_by_this_component_alone': False, 'remaining_tail_starts_at_prime': tail, 'global_Cstar_identified': False, 'second_geometric_audit_exhaustive': False}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    z = assemble()
    s = json.dumps(z, indent=2) + '\n'
    out = ROOT / 'results/small_prime_closure.json'
    if a.check:
        need(out.read_text() == s, 'exact closure report')
    else:
        out.write_text(s, newline='\n')
    print('PASS small-prime closure', z['common_upper_decimal'], flush=True)
if __name__ == '__main__':
    main()
