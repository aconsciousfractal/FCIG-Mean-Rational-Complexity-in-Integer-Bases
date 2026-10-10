"""Reintegrate every declared profile and assemble its exact marked scope."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json
from exact import need, log_lo, out_decimal, prime_trial, sieve
from assemble import profile
from core import CORE, CORE_PRIMES, MIN_NEW_PRIME, mark
ROOT = Path(__file__).resolve().parents[1]

def assemble(root=ROOT):
    plan = json.loads((root / 'profiles_plan.json').read_text())
    contract = json.loads((root / 'profile_contract.json').read_text())
    need(CORE == 114 and CORE_PRIMES == (2, 3, 19) and (MIN_NEW_PRIME == 23), 'fixed scientific core')
    need(contract['core'] == CORE and tuple(contract['core_primes']) == CORE_PRIMES, 'declared core match')
    need(contract['cdf_alpha'] == '3/2', 'new universal CDF used')
    need(contract['first_new_prime'] == MIN_NEW_PRIME, 'first admissible new prime')
    tail = contract['tail_first_prime']
    need(prime_trial(tail) and tail >= MIN_NEW_PRIME, 'next prime starts the external tail')
    need(len({r['stem'] for r in plan}) == len(plan), 'unique profile outputs')
    configured = [p for p in contract['small_primes']]
    need(configured == sorted(set(configured)) and all((prime_trial(p) and p >= 23 for p in configured)), 'declared small prime set')
    need(sorted((r['p'] for r in plan if r['p'] > 1)) == configured, 'every selected small prime present once')
    need(configured == [p for p in sieve(tail) if MIN_NEW_PRIME <= p < tail], 'complete small prime interval through declared tail')
    need(sorted((r['multiplier'] for r in plan if r['p'] == 1)) == [1, 2, 3, 19], 'all four declared marked indices')
    need(set(contract['simple_marked_caps']) == {'114', '228', '342', '2166'}, 'exact simple-cap index domain')
    rows = []
    marked = []
    small = []
    for z in plan:
        p, m, ag = (z['p'], z['multiplier'], z['ag'])
        q, labels = mark(p)
        need(ag in (4, 8, 16), 'audited angular rule')
        need(m in (1, 2, 3, 19) if p == 1 else m == 1, 'supported index')
        need(z['stem'] == (f'q114_K{114 * m}' if p == 1 else f'q{q}'), 'unambiguous stem')
        got = profile(p, root / 'results' / (z['stem'] + '.cells'), root, m, ag)
        stored = json.loads((root / 'results' / (z['stem'] + '.json')).read_text())
        need(set(stored) == set(got) | {'witness_bytes', 'witness_sha256'}, 'exact stored schema')
        need(all((stored[k] == v for k, v in got.items())), 'all exact integral fields')
        need(stored['witness_bytes'] == 8 * got['witnesses'] and len(stored['witness_sha256']) == 64, 'record count and digest schema')
        cap = F(got['total_upper'])
        rows.append(got)
        if p == 1:
            simple = F(contract['simple_marked_caps'][str(114 * m)])
            need(cap < simple, 'strict simple marked cap')
            marked.append({'q': 114, 'required_divisor': 114 * m, 'exact_upper': str(cap), 'upper_decimal': out_decimal(cap), 'simple_cap': str(simple), 'scope': 'all positive M divisible by required_divisor, with extra primes allowed; marked mean nu114 only'})
        else:
            norm = max(F(0), cap) / log_lo(q)
            need(norm < F(2, 85), 'strict small-radical closure')
            small.append({'p': p, 'q': q, 'mean_upper': str(cap), 'normalized_upper': str(norm), 'normalized_upper_decimal': out_decimal(norm), 'scope': 'B=2^a*3^b*19^c*p^d, a,b,c,d,N>=1; normalization ln B'})
    result = {'status': 'PASS', 'cdf_alpha': '3/2', 'marked_caps': marked, 'small_prime_routes': small, 'closed_primes': configured, 'profiles': len(rows), 'geometric_witnesses': sum((r['witnesses'] for r in rows)), 'geometric_corner_evaluations': sum((4 * r['angular_nodes'] * r['witnesses'] for r in rows)), 'arbitrary_prime_valuations': True, 'all_N_positive': True, 'geometric_witnesses_replayed_here': False, 'entire_114p_family_closed_here': False, 'global_Cstar_identified': False}
    if small:
        gamma = max((F(v['normalized_upper']) for v in small))
        result.update(small_prime_common_upper=str(gamma), small_prime_common_upper_decimal=out_decimal(gamma), small_prime_strict_margin=str(F(2, 85) - gamma))
    return result
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    z = assemble()
    s = json.dumps(z, indent=2) + '\n'
    p = ROOT / 'results/profile_certificate.json'
    if a.check:
        need(p.read_text() == s, 'exact profile certificate')
    else:
        p.write_text(s, newline='\n')
    print('PASS profiles assembly', z['profiles'], 'profiles')
