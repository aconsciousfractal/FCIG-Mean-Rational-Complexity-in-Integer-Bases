"""Complete finite residual for rad(B)=105p; exact interval route selection.

This file only accepts its fixed cutoffs after checking both fresh geometric
anchors. The two enumerations cover all positive, independently varying
valuations. All normalized upper bounds are nonnegative to preserve the
direction of division when a marked mean upper bound is negative.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import csv, hashlib, io, json
from exact import need, sieve, prime_trial, phi_hi, log_lo, out_decimal
ROOT = Path(__file__).resolve().parent.parent
CORE = 105
FACTORS = (3, 5, 7)
BMAX = 1000000
PMIN = 31
PMAX = 125
DELTA = F(2, 85)
A1 = F(6, 125)
A2 = F(-2, 25)
ERR = F(90985, 3948)

@lru_cache(None)
def periodic_part(p):
    return phi_hi(F(12, 7 * p)) + F(9, 14 * p)

@lru_cache(None)
def m_free(p):
    need(p >= 7, 'radius-three sharp M-free integration domain')
    return phi_hi(F(9, 2 * p)) + F(27, 16 * p)

def penalty(p, M):
    periodic = periodic_part(p) + ERR / isqrt(M)
    free = m_free(p)
    return (periodic, 'periodic') if periodic <= free else (free, 'M_free')

def upper(M, p):
    k = 2 if M % CORE ** 2 == 0 else 1
    cost, method = penalty(p, M)
    return ((A2 if k == 2 else A1) + cost, k, method)

def normalized(value, B):
    return max(F(0), value) / log_lo(B)

def anchors(root):
    rows = []
    for k, cap in ((1, A1), (2, A2)):
        path = root / f'profiles/results/q105_k{k}.json'
        z = json.loads(path.read_text())
        value = F(z['total_upper'])
        need(z['q'] == CORE and z['k'] == k and (z['transfer_index'] == CORE ** k), 'fresh marked105 divisor')
        need(value < cap, 'fresh marked105 cap required before domain acceptance')
        need(value > 0 if k == 1 else value < 0, 'actual sign of marked105 anchor')
        rows.append(dict(k=k, divisor=CORE ** k, exact_upper=str(value), relaxed_cap=str(cap), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return rows

def serialize(rows):
    need(bool(rows), 'nonempty finite inventory')
    out = io.StringIO(newline='')
    w = csv.DictWriter(out, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
    return out.getvalue().encode()

def make(root=ROOT):
    accepted_anchors = anchors(root)
    ps = [p for p in sieve(PMAX - 1) if p >= PMIN]
    need(ps == [p for p in range(PMIN, PMAX) if prime_trial(p)], 'independent prime list')
    cofactors = []

    def rec(i, m, ee):
        if i == len(FACTORS):
            cofactors.append((m, ee))
            return
        pr = FACTORS[i]
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
                ratio = normalized(value, M)
                rows.append(dict(M=M, p=p, q=CORE * p, exponents=':'.join(map(str, ee + [e])), profile_k=k, penalty_method=method, route='analytic' if ratio < DELTA else 'witness', mean_upper=str(value), normalized_upper=str(ratio)))
            pp *= p
            e += 1
    rows.sort(key=lambda r: r['M'])
    need(len({r['M'] for r in rows}) == len(rows), 'unique finite modulus')
    need(all((r['M'] % 2 for r in rows)), 'all bases in this family are odd')
    check = []
    for M in range(CORE * PMIN, BMAX, CORE):
        v = M
        for s in FACTORS:
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
        ratio = normalized(value, B)
        need(0 < ratio < DELTA, 'infinite region ' + name)
        gates.append(dict(name=name, p_min=p, M_min=M, base_min=B, relaxed_cap=str(cap), penalty_method=method, mean_upper=str(value), normalized_upper=str(ratio), normalized_upper_decimal=out_decimal(ratio), margin_to_delta=str(DELTA - ratio)))
    witnesses = [r for r in rows if r['route'] == 'witness']
    meta = dict(status='PASS', core=CORE, core_factors=list(FACTORS), base_upper_exclusive=BMAX, prime_range=[PMIN, PMAX], prime_count=len(ps), cofactor_count=len(cofactors), marked_anchors=accepted_anchors, all_moduli=len(rows), analytic_moduli=len(rows) - len(witnesses), witness_moduli=len(witnesses), witness_states=sum((r['M'] for r in witnesses)), witness_representatives=sum((r['M'] // 2 + 1 for r in witnesses)), maximum_witness_modulus=max((r['M'] for r in witnesses), default=0), infinite_regions=gates, second_profile_used_at_N1=sum((r['profile_k'] == 2 for r in rows)), M_free_chosen=sum((r['penalty_method'] == 'M_free' for r in rows)), finite_moduli_are_bases_at_N1=True, all_positive_prime_valuations=True, normalization='log B, including N>=2', nonnegative_normalized_uppers=True, other_four_prime_supports_closed=False, all_five_prime_supports_closed=False)
    return (rows, meta)

def write(root=ROOT):
    rows, meta = make(root)
    out = root / 'results'
    out.mkdir(parents=True, exist_ok=True)
    for name, chosen in [('finite_inventory.csv', rows), ('witness_inventory.csv', [r for r in rows if r['route'] == 'witness'])]:
        (out / name).write_bytes(serialize(chosen))
    (out / 'inventory_certificate.json').write_text(json.dumps(meta, sort_keys=True, indent=2) + '\n', newline='\n')
    (out / 'finite_domain.txt').write_text(''.join((f"{r['M']} {r['p']}\n" for r in rows if r['route'] == 'witness')), newline='\n')
    return meta
if __name__ == '__main__':
    print(json.dumps(write(), indent=2), flush=True)
