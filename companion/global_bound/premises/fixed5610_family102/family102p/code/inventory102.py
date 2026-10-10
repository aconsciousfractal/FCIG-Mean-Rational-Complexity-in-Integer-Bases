"""Exhaustive all-valuation residual for rad(B)=102p, with rational gates.

The numerical anchor caps below are proof obligations discharged separately
by the complete continuous profiles. This module only enumerates and assembles.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from math import isqrt
import csv
import io
import json
from exact import need, sieve, prime_trial, phi_hi, log_lo, out_decimal
CORE = 102
PMIN = 37
PMAX = 160
BMAX = 200000
SMALL = [19, 23, 29, 31]
DELTA = F(2, 85)
ANCHORS = {102: F(6, 125), 204: F(9, 200), 10404: F(-1, 100)}
K_GEOMETRY = F(46339, 250000)
J4_UPPER = F(51, 10)
A_DENSITY = F(40, 51)
C_DENSITY = F(14744, 4029)
PER_LINEAR = F(92678, 159375)
PER_ERROR = F(7372, 395)
FREE_LINEAR = F(139017, 62500)

@lru_cache(None)
def periodic_part(p):
    need(p >= 19, 'scalar prime lower bound')
    return phi_hi(F(80, 51 * p)) + PER_LINEAR / p

@lru_cache(None)
def m_free(p):
    need(p >= 19, 'M-free integration domain')
    return phi_hi(F(6, p)) + FREE_LINEAR / p

def penalty(p, M):
    need(M >= 1, 'positive modulus')
    periodic = periodic_part(p) + PER_ERROR / isqrt(M)
    free = m_free(p)
    return (periodic, 'periodic') if periodic <= free else (free, 'M_free')

def upper(M, p):
    available = [(a, K) for K, a in ANCHORS.items() if M % K == 0]
    need(bool(available), 'at least one certified anchor divisor')
    anchor, K = min(available)
    cost, method = penalty(p, M)
    return (anchor + cost, K, method)

def normalized(value, B):
    return max(F(0), value) / log_lo(B)

def serialize(rows):
    need(bool(rows), 'nonempty inventory')
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return out.getvalue().encode()

def make():
    need(SMALL == [p for p in sieve(PMIN - 1) if p >= 19], 'complete small-prime prefix with no gap')
    need(PER_LINEAR == 4 * A_DENSITY * K_GEOMETRY, 'periodic linear coefficient from sharp geometry')
    need(PER_ERROR == C_DENSITY * J4_UPPER, 'periodic square-root coefficient')
    need(FREE_LINEAR == 12 * K_GEOMETRY, 'M-free sharp coefficient')
    primes = [p for p in sieve(PMAX - 1) if p >= PMIN]
    need(primes == [p for p in range(PMIN, PMAX) if prime_trial(p)], 'independent residual prime list')
    cofactors = []

    def rec(i, m, exponents):
        if i == 3:
            cofactors.append((m, exponents))
            return
        prime = (2, 3, 17)[i]
        n, exponent = (m * prime, 1)
        while n * PMIN < BMAX:
            rec(i + 1, n, exponents + [exponent])
            n *= prime
            exponent += 1
    rec(0, 1, [])
    cofactors.sort()
    rows = []
    for p in primes:
        power, exponent = (p, 1)
        while CORE * power < BMAX:
            for cofactor, exponents in cofactors:
                M = cofactor * power
                if M >= BMAX:
                    break
                value, K, method = upper(M, p)
                ratio = normalized(value, M)
                rows.append(dict(M=M, p=p, q=CORE * p, exponents=':'.join(map(str, exponents + [exponent])), anchor_divisor=K, penalty_method=method, route='analytic' if ratio < DELTA else 'witness', mean_upper=str(value), normalized_upper=str(ratio)))
            power *= p
            exponent += 1
    rows.sort(key=lambda r: r['M'])
    need(len({r['M'] for r in rows}) == len(rows), 'unique finite modulus')
    check = []
    for M in range(CORE * PMIN, BMAX, CORE):
        remainder = M
        for p in (2, 3, 17):
            while remainder % p == 0:
                remainder //= p
        if remainder == 1:
            continue
        p = remainder
        for d in range(2, isqrt(remainder) + 1):
            if remainder % d == 0:
                p = d
                break
        if not PMIN <= p < PMAX or not prime_trial(p):
            continue
        while remainder % p == 0:
            remainder //= p
        if remainder == 1:
            check.append((M, p))
    need(check == [(r['M'], r['p']) for r in rows], 'independent complete cofactor factorization')
    gates = []
    for name, p, M, B, K in [('large_p', PMAX, CORE * PMAX, CORE * PMAX, CORE), ('N_ge_2', PMIN, (CORE * PMIN) ** 2, CORE * PMIN, CORE ** 2), ('large_base', PMIN, BMAX, BMAX, CORE)]:
        cost, method = penalty(p, M)
        value = ANCHORS[K] + cost
        ratio = normalized(value, B)
        need(0 < ratio < DELTA, 'infinite region ' + name)
        gates.append(dict(name=name, p_scalar_min=p, M_min=M, base_min=B, anchor_divisor=K, relaxed_cap=str(ANCHORS[K]), penalty_method=method, mean_upper=str(value), normalized_upper=str(ratio), normalized_upper_decimal=out_decimal(ratio), margin_to_delta=str(DELTA - ratio)))
    witnesses = [r for r in rows if r['route'] == 'witness']
    meta = dict(status='PASS', computational_scope='exact gates and complete residual inventory; anchors are separate obligations', core=CORE, base_upper_exclusive=BMAX, prime_range=[PMIN, PMAX], small_primes=SMALL, prime_count=len(primes), cofactor_count=len(cofactors), all_moduli=len(rows), analytic_moduli=len(rows) - len(witnesses), witness_moduli=len(witnesses), witness_states=sum((r['M'] for r in witnesses)), witness_representatives=sum((r['M'] // 2 + 1 for r in witnesses)), maximum_witness_modulus=max((r['M'] for r in witnesses), default=0), infinite_regions=gates, anchor_usage_at_N1={str(K): sum((r['anchor_divisor'] == K for r in rows)) for K in ANCHORS}, M_free_chosen=sum((r['penalty_method'] == 'M_free' for r in rows)), finite_moduli_are_bases_at_N1=True, all_positive_prime_valuations=True, other_four_prime_supports_closed=False, all_five_prime_supports_closed=False, minimum_penalty_monotonic=True, normalization='mu(B^N)/ln(B)')
    return (rows, meta)

def write(root):
    rows, meta = make()
    output = root / 'results'
    output.mkdir(parents=True, exist_ok=True)
    for name, selected in [('finite_inventory.csv', rows), ('witness_inventory.csv', [r for r in rows if r['route'] == 'witness'])]:
        (output / name).write_bytes(serialize(selected))
    (output / 'inventory_certificate.json').write_text(json.dumps(meta, sort_keys=True, indent=2) + '\n', newline='\n')
    (output / 'finite_domain.txt').write_text(''.join((f"{r['M']} {r['p']}\n" for r in rows if r['route'] == 'witness')), newline='\n')
    return meta
if __name__ == '__main__':
    print(json.dumps(write(Path(__file__).resolve().parents[1]), indent=2), flush=True)
