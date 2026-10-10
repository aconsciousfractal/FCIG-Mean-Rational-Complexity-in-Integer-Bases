"""Exhaustive residual for rad(B)=66p; exact rational route selection."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import csv, io, json
from exact import need, sieve, prime_trial, phi_hi, log_lo, out_decimal
BMAX = 1000000
PMAX = 350
PMIN = 47
CORE = 66
DELTA = F(2, 85)
A1 = F(127, 1000)
A2 = F(2, 125)
ERR = F(84371299, 4278582)

@lru_cache(None)
def periodic_part(p):
    return phi_hi(F(50, 33 * p)) + F(25, 44 * p)

@lru_cache(None)
def m_free(p):
    need(p >= 11, 'sharp M-free domain')
    return phi_hi(F(6, p)) + F(9, 4 * p)

def penalty(p, M):
    periodic = periodic_part(p) + ERR / isqrt(M)
    free = m_free(p)
    return (periodic, 'periodic') if periodic <= free else (free, 'M_free')

def upper(M, p):
    k = 2 if M % CORE ** 2 == 0 else 1
    cost, method = penalty(p, M)
    return ((A2 if k == 2 else A1) + cost, k, method)

def serialize(rows):
    out = io.StringIO(newline='')
    w = csv.DictWriter(out, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    return out.getvalue().encode()

def make():
    ps = [p for p in sieve(PMAX - 1) if p >= PMIN]
    need(ps == [p for p in range(PMIN, PMAX) if prime_trial(p)], 'independent prime list')
    cofactors = []

    def rec(i, m, ee):
        if i == 3:
            cofactors.append((m, ee))
            return
        pr = (2, 3, 11)[i]
        n, e = (m * pr, 1)
        while n * PMIN < BMAX:
            rec(i + 1, n, ee + [e])
            n *= pr
            e += 1
    rec(0, 1, [])
    cofactors.sort()
    rows = []
    for p in ps:
        pp, e = (p, 1)
        while CORE * pp < BMAX:
            for s, ee in cofactors:
                M = s * pp
                if M >= BMAX:
                    break
                value, k, method = upper(M, p)
                rows.append(dict(M=M, p=p, q=CORE * p, exponents=':'.join(map(str, ee + [e])), profile_k=k, penalty_method=method, route='analytic' if value < DELTA * log_lo(M) else 'witness', mean_upper=str(value), normalized_upper=str(value / log_lo(M))))
            pp *= p
            e += 1
    rows.sort(key=lambda r: r['M'])
    need(len({r['M'] for r in rows}) == len(rows), 'unique finite modulus')
    check = []
    for M in range(CORE * PMIN, BMAX, CORE):
        v = M
        for s in (2, 3, 11):
            while v % s == 0:
                v //= s
        if v == 1:
            continue
        p = v
        for s in range(2, isqrt(v) + 1):
            if v % s == 0:
                p = s
                break
        if not PMIN <= p < PMAX or not prime_trial(p):
            continue
        while v % p == 0:
            v //= p
        if v == 1:
            check.append((M, p))
    need(check == [(r['M'], r['p']) for r in rows], 'independent cofactor scan')
    gates = []
    for name, p, M, B, cap in [('large_p', PMAX, CORE * PMAX, CORE * PMAX, A1), ('N_ge_2', PMIN, (CORE * PMIN) ** 2, CORE * PMIN, A2), ('large_base', PMIN, BMAX, BMAX, A1)]:
        cost, method = penalty(p, M)
        value = cap + cost
        ratio = value / log_lo(B)
        need(0 < ratio < DELTA, 'infinite region ' + name)
        gates.append(dict(name=name, p_min=p, M_min=M, base_min=B, relaxed_cap=str(cap), penalty_method=method, mean_upper=str(value), normalized_upper=str(ratio), normalized_upper_decimal=out_decimal(ratio), margin_to_delta=str(DELTA - ratio)))
    witnesses = [r for r in rows if r['route'] == 'witness']
    meta = dict(status='PASS', core=CORE, base_upper_exclusive=BMAX, prime_range=[PMIN, PMAX], prime_count=len(ps), cofactor_count=len(cofactors), all_moduli=len(rows), analytic_moduli=len(rows) - len(witnesses), witness_moduli=len(witnesses), witness_states=sum((r['M'] for r in witnesses)), witness_representatives=sum((r['M'] // 2 + 1 for r in witnesses)), maximum_witness_modulus=max((r['M'] for r in witnesses), default=0), infinite_regions=gates, second_profile_used_at_N1=sum((r['profile_k'] == 2 for r in rows)), M_free_chosen=sum((r['penalty_method'] == 'M_free' for r in rows)), finite_moduli_are_bases_at_N1=True, all_positive_prime_valuations=True, other_four_prime_supports_closed=False, all_five_prime_supports_closed=False)
    return (rows, meta)

def write(root):
    rows, meta = make()
    out = root / 'results'
    out.mkdir(parents=True, exist_ok=True)
    for name, chosen in [('finite_inventory.csv', rows), ('witness_inventory.csv', [r for r in rows if r['route'] == 'witness'])]:
        (out / name).write_bytes(serialize(chosen))
    (out / 'inventory_certificate.json').write_text(json.dumps(meta, sort_keys=True, indent=2) + '\n', newline='\n')
    (out / 'finite_domain.txt').write_text(''.join((f"{r['M']} {r['p']}\n" for r in rows if r['route'] == 'witness')), newline='\n')
    return meta
if __name__ == '__main__':
    print(json.dumps(write(Path(__file__).resolve().parent.parent), indent=2), flush=True)
