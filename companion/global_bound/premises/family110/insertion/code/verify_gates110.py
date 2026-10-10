#!/usr/bin/env python3
"""Exact scalar gates and an independently generated multiplicative closure.

The fixed anchor values are geometric proof obligations.  This checker
reconstructs their arithmetic consequences; it does not certify profiles
or consume the physical witness streams.  Root code is hashed, never run.
"""
import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import insertion110 as E

CHECKS = 0
DELTA = F(2, 85)
PMIN, PMAX, BMAX = 37, 150, 300000
SMALL = [13, 17, 19, 23, 29, 31]
CAPS = {110: F(17, 500), 220: F(23, 1000), 12100: F(-6, 125)}


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def upper_fraction(x, scale=10**18):
    x = F(x)
    return F(-((-x.numerator * scale) // x.denominator), scale)


def multiplicative_closure(primes):
    """Start at110p and multiply by2,5,11,p while remaining below BMAX.

    Every target base occurs by induction on the excess exponent sum;
    every generated value has exactly the specified four-prime support.
    This uses neither recursive exponent products nor a factorization scan.
    """
    result = set()
    per_prime = []
    transitions = 0
    for p in primes:
        first = 110*p
        seen = {first} if first < BMAX else set()
        queue = list(seen)
        while queue:
            B = queue.pop()
            for prime in (2, 5, 11, p):
                nxt = B*prime
                if nxt < BMAX:
                    transitions += 1
                    if nxt not in seen:
                        seen.add(nxt)
                        queue.append(nxt)
        need(not (result & seen), 'closures for distinct new primes are disjoint')
        result.update(seen)
        per_prime.append(dict(p=p, seed=first, vertices=len(seen)))
    return result, dict(algorithm='multiplicative reachability from each110p seed',
                        vertices=len(result), admissible_directed_edges=transitions,
                        prime_components=per_prime,
                        completeness_proof='remove one excess prime factor to reduce any positive-valuation base toward its seed')


def verify(inputs):
    filenames = ['results__inventory_certificate.json', 'results__finite_inventory.csv',
                 'code__inventory110.py', 'code__exact.py', 'PROVENANCE.json']
    blobs = {name: (inputs/name).read_bytes() for name in filenames}
    provenance = json.loads(blobs['PROVENANCE.json'])
    copied = {Path(row['copy']).name for row in provenance['sources']}
    need(copied == set(filenames)-{'PROVENANCE.json'}, 'exact four-source snapshot inventory')
    for row in provenance['sources']:
        name = Path(row['copy']).name
        need(sha256(blobs[name]).hexdigest() == row['sha256'] and len(blobs[name]) == row['bytes'],
             'snapshot source byte identity')
    root = json.loads(blobs['results__inventory_certificate.json'])
    need(root['status'] == 'PASS' and root['core'] == 110, 'root inventory identity')
    need(root['small_primes'] == SMALL and root['prime_range'] == [PMIN, PMAX]
         and root['base_upper_exclusive'] == BMAX, 'fixed exhaustive partition parameters')
    need(root['normalization'] == 'mu(B^N)/ln(B)', 'full normalization declaration')
    need(root['finite_moduli_are_bases_at_N1'] and root['all_positive_prime_valuations'],
         'finite domain is N1 with all positive valuations')
    specifications = [
        ('large_p', PMAX, 110*PMAX, 110*PMAX, 110, 'M_free'),
        ('N_ge_2', PMIN, 110*PMIN, (110*PMIN)**2, 12100, 'periodic'),
        ('large_base', PMIN, BMAX, BMAX, 110, 'periodic'),
    ]
    source_gates = {g['name']: g for g in root['infinite_regions']}
    need(set(source_gates) == {s[0] for s in specifications} and len(root['infinite_regions']) == 3,
         'three distinct infinite gates')
    gates = []
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
             'source signed anchor and penalty choice')
        root_mean, root_ratio = F(old['mean_upper']), F(old['normalized_upper'])
        need(0 <= root_mean - mean < F(1, 10**12), 'root mean above80-term upper')
        need(0 <= root_ratio - ratio < F(1, 10**12), 'root gate above80-term upper')
        gates.append(dict(name=name, p_scalar=p, B_min=B, M_min=M, anchor_divisor=K,
                          anchor_assumed_upper=E.fraction_json(cap), selected=chosen,
                          mean_upper=E.fraction_json(mean), normalized_upper=E.fraction_json(ratio),
                          normalized_upper_decimal=E.decimal_ceiling(ratio),
                          margin_to_delta=E.fraction_json(DELTA-ratio),
                          rounded_root_difference_upper=E.fraction_json(upper_fraction(root_ratio-ratio))))
    need([p for p in range(13, PMIN) if E.is_prime(p)] == SMALL, 'complete lower prime prefix')
    primes = [p for p in range(PMIN, PMAX) if E.is_prime(p)]
    need(len(primes) == 24 == root['prime_count'], 'complete residual prime set')
    need([p for p in range(149, 152) if E.is_prime(p)] == [149, 151], 'upper prime gap')
    rows = list(csv.DictReader(io.StringIO(blobs['results__finite_inventory.csv'].decode())))
    need(len(rows) == 292 == root['all_moduli'], 'all292 residual rows')
    seen = set()
    analytic = witnesses = witness_states = witness_representatives = 0
    maximum_witness = M_free = 0
    anchor_usage = {K: 0 for K in CAPS}
    max_mean_error = max_ratio_error = F(0)
    summaries = []
    for row in rows:
        M, p, q, K = (int(row[k]) for k in ('M', 'p', 'q', 'anchor_divisor'))
        need(M not in seen and M < BMAX and p in primes and q == 110*p,
             'residual row domain and uniqueness')
        seen.add(M)
        need(E.validate_regime(p, M, True), 'full residual four-prime support')
        powers = [int(v) for v in row['exponents'].split(':')]
        need(len(powers) == 4 and all(v >= 1 for v in powers), 'positive valuations')
        reconstructed = 1
        for prime, exponent in zip((2, 5, 11, p), powers):
            reconstructed *= prime**exponent
        need(M == reconstructed, 'valuation product equality')
        cap, divisor = min((cap, divisor) for divisor, cap in CAPS.items() if M % divisor == 0)
        need(divisor == K, 'best admissible signed anchor')
        periodic, free = E.periodic_upper(p, M), E.free_upper(p)
        method = 'periodic' if periodic <= free else 'M_free'
        mean = cap + min(periodic, free)
        ratio = E.normalized_upper(p, M, M, cap)
        root_mean, root_ratio = F(row['mean_upper']), F(row['normalized_upper'])
        need(0 <= root_mean-mean < F(1, 10**12), 'residual mean rounding direction')
        need(0 <= root_ratio-ratio < F(1, 10**12), 'residual ratio rounding direction')
        need(row['penalty_method'] == method, 'residual penalty selection')
        route = 'analytic' if ratio < DELTA else 'witness'
        need(row['route'] == route, 'stable route under refined logarithm bounds')
        analytic += route == 'analytic'
        witnesses += route == 'witness'
        if route == 'witness':
            witness_states += M
            witness_representatives += M//2+1
            maximum_witness = max(maximum_witness, M)
        M_free += method == 'M_free'
        anchor_usage[K] += 1
        max_mean_error = max(max_mean_error, root_mean-mean)
        max_ratio_error = max(max_ratio_error, root_ratio-ratio)
        summaries.append(dict(M=M, p=p, anchor_divisor=K, route=route,
                              normalized_upper_decimal=E.decimal_ceiling(ratio)))
    need((analytic, witnesses) == (185, 107) == (root['analytic_moduli'], root['witness_moduli']),
         'unchanged residual route totals')
    need((witness_states, witness_representatives) == (4855510, 2427862)
         == (root['witness_states'], root['witness_representatives']), 'physical domain totals')
    need(maximum_witness == root['maximum_witness_modulus'], 'maximum physical modulus')
    need(anchor_usage == {110: 98, 220: 194, 12100: 0}
         and {str(K): n for K, n in anchor_usage.items()} == root['anchor_usage_at_N1'], 'anchor usage')
    need(M_free == 32 == root['M_free_chosen'], 'M-free route count')
    independent, closure = multiplicative_closure(primes)
    need(independent == seen, 'independent multiplicative closure equals the full inventory')
    return dict(status='PASS', checks=CHECKS, anchors={str(k): E.fraction_json(v) for k, v in CAPS.items()},
                anchors_are_geometric_proof_obligations=True, profile_geometry_verified_here=False,
                infinite_gates=gates, residual_rows=len(rows), analytic_rows=analytic, witness_rows=witnesses,
                witness_states=witness_states, witness_representatives=witness_representatives,
                complete_all_positive_valuations_at_N1=True, independent_multiplicative_closure=closure,
                maximum_mean_rounding_difference_upper=E.fraction_json(upper_fraction(max_mean_error)),
                maximum_ratio_rounding_difference_upper=E.fraction_json(upper_fraction(max_ratio_error)),
                residual_summaries=summaries,
                input_files=[dict(file=name, bytes=len(raw), sha256=sha256(raw).hexdigest())
                             for name, raw in sorted(blobs.items())],
                claim_scope='exact insertion constants, three conditional anchor gates, and residual arithmetic; profiles and physical streams are separate components',
                family110p_closed_here=False, Cstar_exact_proved_here=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(result, indent=2)+'\n').encode()
    args.output.write_bytes(raw)
    print(json.dumps(dict(status='PASS', checks=CHECKS,
                         gates={g['name']:g['normalized_upper_decimal'] for g in result['infinite_gates']},
                         residual_rows=result['residual_rows'], analytic_rows=result['analytic_rows'],
                         witness_rows=result['witness_rows'], output_bytes=len(raw),
                         output_sha256=sha256(raw).hexdigest()), indent=2))


if __name__ == '__main__':
    main()
