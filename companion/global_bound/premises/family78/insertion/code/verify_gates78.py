#!/usr/bin/env python3
"""Exact80-term crosscheck of the fixed78p gates and residual arithmetic.

Input anchor inequalities are declared assumptions here, discharged by
the separate continuous profiles in the complete project.  This script
does not promote arithmetic agreement into proof of those geometries.
"""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from math import gcd, isqrt
from pathlib import Path
import insertion78 as E

CHECKS = 0
DELTA = F(2, 85)
CAPS = {78: F(3, 40), 156: F(7, 100), 6084: F(1, 250)}


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def upper_fraction(x, scale=10**18):
    x = F(x)
    return F(-((-x.numerator * scale) // x.denominator), scale)


def verify(inputs):
    paths = {name: inputs / name for name in (
        'results__inventory_certificate.json', 'results__finite_inventory.csv',
        'code__inventory78.py', 'code__exact.py', 'PROVENANCE.json')}
    blobs = {name: path.read_bytes() for name, path in paths.items()}
    provenance = json.loads(blobs['PROVENANCE.json'])
    for row in provenance['sources']:
        name = Path(row['copy']).name
        need(name in blobs and sha256(blobs[name]).hexdigest() == row['sha256']
             and len(blobs[name]) == row['bytes'], 'snapshot source byte identity')
    root = json.loads(blobs['results__inventory_certificate.json'])
    need(root['status'] == 'PASS' and root['core'] == 78, 'root inventory identity')
    need(root['small_primes'] == [17, 19, 23, 29, 31, 37], 'small-prime prefix')
    gates = []
    specifications = [
        ('large_p', 200, 15600, 15600, 78, 'M_free'),
        ('N_ge_2', 41, 3198, 10227204, 6084, 'periodic'),
        ('large_base', 41, 250000, 250000, 78, 'periodic'),
    ]
    source_gates = {g['name']: g for g in root['infinite_regions']}
    need(set(source_gates) == {s[0] for s in specifications}, 'complete three-gate inventory')
    for name, p, B, M, K, expected in specifications:
        cap = CAPS[K]
        periodic, free = E.periodic_upper(p, M), E.free_upper(p)
        chosen = 'periodic' if periodic <= free else 'M_free'
        cost = min(periodic, free)
        mean = cap + cost
        ratio = E.normalized_upper(p, B, M, cap)
        need(chosen == expected and 0 < ratio < DELTA, 'fixed infinite gate')
        old = source_gates[name]
        need((old['p_scalar_min'], old['base_min'], old['M_min'], old['anchor_divisor'])
             == (p, B, M, K), 'source gate domain')
        need(F(old['relaxed_cap']) == cap and old['penalty_method'] == chosen,
             'source cap and penalty choice')
        root_mean, root_ratio = F(old['mean_upper']), F(old['normalized_upper'])
        need(0 <= root_mean - mean < F(1, 10**12), 'root rounded mean above80-term upper')
        need(0 <= root_ratio - ratio < F(1, 10**12), 'root rounded gate above80-term upper')
        gates.append(dict(name=name, p_scalar=p, B_min=B, M_min=M, anchor_divisor=K,
                          anchor_assumed_upper=E.fraction_json(cap), selected=chosen,
                          mean_upper=E.fraction_json(mean), normalized_upper=E.fraction_json(ratio),
                          normalized_upper_decimal=E.decimal_ceiling(ratio),
                          margin_to_delta=E.fraction_json(DELTA-ratio),
                          rounded_root_difference_upper=E.fraction_json(upper_fraction(root_ratio-ratio))))
    need([p for p in range(17, 41) if E.is_prime(p)] == root['small_primes'], 'lower prime gap and prefix')
    primes = [p for p in range(41, 200) if E.is_prime(p)]
    need(len(primes) == 34 and [p for p in range(199, 212) if E.is_prime(p)] == [199, 211],
         'upper prime gap and residual primes')
    rows = list(csv.DictReader(io.StringIO(blobs['results__finite_inventory.csv'].decode())))
    need(len(rows) == 486 == root['all_moduli'], 'all486 residual rows')
    seen = set()
    analytic = witnesses = 0
    maximum_mean_error = maximum_ratio_error = F(0)
    summaries = []
    for row in rows:
        M, p, q, K = (int(row[k]) for k in ('M', 'p', 'q', 'anchor_divisor'))
        need(M not in seen and M < 250000 and p in primes and q == 78*p,
             'residual row domain and uniqueness')
        seen.add(M)
        need(E.validate_regime(p, M, True), 'residual four-prime physical support')
        powers = [int(v) for v in row['exponents'].split(':')]
        need(len(powers) == 4 and all(v >= 1 for v in powers), 'positive valuations')
        reconstructed = 1
        for prime, exponent in zip((2, 3, 13, p), powers):
            reconstructed *= prime**exponent
        need(M == reconstructed, 'valuation product equality')
        available = [(cap, divisor) for divisor, cap in CAPS.items() if M % divisor == 0]
        cap, divisor = min(available)
        need(divisor == K, 'best admissible anchor divisor')
        periodic, free = E.periodic_upper(p, M), E.free_upper(p)
        method = 'periodic' if periodic <= free else 'M_free'
        mean = cap + min(periodic, free)
        ratio = E.normalized_upper(p, M, M, cap)
        root_mean, root_ratio = F(row['mean_upper']), F(row['normalized_upper'])
        need(0 <= root_mean-mean < F(1, 10**12), 'residual mean rounding direction')
        need(0 <= root_ratio-ratio < F(1, 10**12), 'residual ratio rounding direction')
        need(row['penalty_method'] == method, 'residual penalty choice')
        route = 'analytic' if ratio < DELTA else 'witness'
        need(row['route'] == route, 'residual gate classification stable under log refinement')
        analytic += route == 'analytic'
        witnesses += route == 'witness'
        maximum_mean_error = max(maximum_mean_error, root_mean-mean)
        maximum_ratio_error = max(maximum_ratio_error, root_ratio-ratio)
        summaries.append(dict(M=M, p=p, anchor_divisor=K, route=route,
                              normalized_upper_decimal=E.decimal_ceiling(ratio)))
    need((analytic, witnesses) == (304, 182) == (root['analytic_moduli'], root['witness_moduli']),
         'unchanged residual route totals')
    # Independent factorization scan: every multiple78, no recursive exponent tree.
    independently_seen = set()
    for M in range(78*41, 250000, 78):
        remainder = M
        for prime in (2, 3, 13):
            while remainder % prime == 0:
                remainder //= prime
        for p in primes:
            if remainder % p == 0:
                remaining = remainder
                while remaining % p == 0:
                    remaining //= p
                if remaining == 1:
                    independently_seen.add(M)
                break
    need(independently_seen == seen, 'independent residual completeness scan')
    return dict(status='PASS', checks=CHECKS, anchors={str(k): E.fraction_json(v) for k, v in CAPS.items()},
                anchors_are_geometric_proof_obligations=True, profile_geometry_verified_here=False,
                infinite_gates=gates, residual_rows=len(rows), analytic_rows=analytic, witness_rows=witnesses,
                complete_all_positive_valuations_at_N1=True,
                maximum_mean_rounding_difference_upper=E.fraction_json(upper_fraction(maximum_mean_error)),
                maximum_ratio_rounding_difference_upper=E.fraction_json(upper_fraction(maximum_ratio_error)),
                residual_summaries=summaries,
                input_files=[dict(file=name, bytes=len(raw), sha256=sha256(raw).hexdigest())
                             for name, raw in sorted(blobs.items())],
                claim_scope='exact insertion constants, three conditional anchor gates, and residual arithmetic; finite physical witnesses are a separate component',
                family78p_closed_here=False, Cstar_exact_proved_here=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(result, indent=2) + '\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps(dict(status='PASS', checks=CHECKS,
                         gates={g['name']:g['normalized_upper_decimal'] for g in result['infinite_gates']},
                         residual_rows=result['residual_rows'], analytic_rows=result['analytic_rows'],
                         witness_rows=result['witness_rows'], output_bytes=len(raw),
                         output_sha256=sha256(raw).hexdigest()), indent=2))


if __name__ == '__main__':
    main()
