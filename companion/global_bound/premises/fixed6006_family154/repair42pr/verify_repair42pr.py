"""Exact certificate for joint repair of the marked denominator 42*p*r.

Python standard library only. The finite residue table proves an infinite-prime
statement by its nonzero determinants; finite-prime and physical checks below
are regressions, not a replacement for that proof. Historical marked42 profile
bounds are explicit imported premises. No claim of a complete 42pr closure.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import time
ROOT = Path(__file__).resolve().parent
A, C = (F(5, 7), F(753, 217))
K_UP, J4_UP, J6_UP = (F(46339, 250000), F(51, 10), F(25, 4))
DELTA = F(2, 85)
INTERVALS = ((F(0), F(1, 3)), (F(1, 3), F(2, 5)), (F(2, 5), F(1, 2)))
GROUPS = (((0, 6, 12, 24, 30), (((1, -1), (2, -1)),) * 3), ((1, 13, 19, 31, 37), (((1, 0), (2, -1)),) * 3), ((2,), (((1, -1), (3, -1)),) * 3), ((3,), (((2, -1), (4, -1)), ((2, -1), (5, -2)), ((2, -1), (5, -2)))), ((4,), (((5, -1), (6, -1)), ((3, -1), (5, -1)), ((3, -1), (5, -3)))), ((5, 11, 17, 29, 41), (((1, 0), (4, -1)),) * 3), ((7, 9), (((2, -1), (6, -1)), ((2, -1), (3, -2)), ((2, -1), (3, -2)))), ((8, 14, 18, 20, 32, 38), (((3, -1), (4, -1)),) * 3), ((10, 16, 22, 28, 34, 36), (((2, -1), (3, -1)),) * 3), ((15, 21, 27, 33), (((2, -1), (4, -1)), ((2, -1), (3, -2)), ((2, -1), (3, -2)))), ((23, 25), (((1, 0), (6, -1)), ((1, 0), (3, -2)), ((1, 0), (3, -2)))), ((26,), (((4, -1), (6, -1)), ((1, -1), (4, -1)), ((1, -1), (5, -3)))), ((35, 39), (((4, -1), (6, -1)), ((4, -1), (5, -2)), ((4, -1), (5, -2)))), ((40,), (((5, -1), (6, -1)), ((2, -1), (5, -1)), ((2, -1), (5, -3)))))
TABLE = {z: row for zs, row in GROUPS for z in zs}
CHECKS = Counter()

def require(test, label):
    CHECKS[label] += 1
    if not test:
        raise AssertionError(label)

def enc(x):
    x = F(x)
    return f'{x.numerator}/{x.denominator}'

def upward(x, digits=35):
    x, s = (F(x), 10 ** digits)
    return F(-(-x.numerator * s // x.denominator), s)

def downward(x, digits=35):
    x, s = (F(x), 10 ** digits)
    return F(x.numerator * s // x.denominator, s)

def decimal_upper(x, digits=15):
    s = 10 ** digits
    k = -(-F(x).numerator * s // F(x).denominator)
    return f'{k // s}.{k % s:0{digits}d}'

@lru_cache(None)
def log_interval(x):
    """Directed rational artanh series, 80 terms, explicit positive tail."""
    x = F(x)
    if x <= 0:
        raise ValueError('positive logarithm argument required')
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def part(z):
        s, t, z2 = (F(0), z, z * z)
        for n in range(80):
            s += 2 * t / (2 * n + 1)
            t *= z2
        return (s, s + 2 * t / (161 * (1 - z2)))
    lo, hi = part((x - 1) / (x + 1))
    l2, h2 = part(F(1, 3))
    return (lo + k * l2, hi + k * h2) if k >= 0 else (lo + k * h2, hi + k * l2)

def phi_upper(x):
    x = F(x)
    return log_interval(1 + x)[1] + x * log_interval(1 + 1 / x)[1]

def atan_interval(x, even_last=80):
    x = F(x)
    if not 0 <= x <= F(3, 4):
        raise ValueError('arctangent series domain')
    a = sum(((-1) ** n * x ** (2 * n + 1) / F(2 * n + 1) for n in range(even_last + 1)), F(0))
    return (a - x ** (2 * even_last + 3) / F(2 * even_last + 3), a)

def G_upper(y):
    """Upper for the exact capped-linear integral for every y>0."""
    y = F(y)
    if y <= 0:
        raise ValueError('positive argument required')
    if y <= F(4, 3):
        return phi_upper(y / 2) + y * K_UP
    return log_interval(1 + y * y)[1] / 2 + y * atan_interval(1 / y)[1]

def periodic_upper(p, r, M):
    """Scalar upper, including non-prime monotone boundaries p,r>=11."""
    if min(p, r) < 11 or M < 1:
        raise ValueError('invalid scalar boundary')
    lam = F(1, p) + F(1, r)
    require(6 * A * lam <= F(4, 3), 'periodic integral domain')
    return phi_upper(3 * A * lam) + 6 * A * K_UP * lam + 2 * C * J6_UP / isqrt(M)

def free_upper(p, r):
    return G_upper(18 * (F(1, p) + F(1, r)))

def residue_certificate():
    require(sorted(TABLE) == list(range(42)), 'all42 residue classes')
    require(sum((len(zs) for zs, _ in GROUPS)) == 42, 'residue groups disjoint')
    rows, dets = ([], Counter())
    for z in range(42):
        for h, (left, right) in enumerate(INTERVALS):
            pair = TABLE[z][h]
            for i, j in pair:
                require(1 <= i <= 6, 'radius6 coefficient')
                require(gcd(i * z + j, 42) == 1, 'core42 unit')
                for t in (left, right):
                    require(abs(j + i * t) <= 1, 'closed shear endpoint')
            (i, j), (a, b) = pair
            det = abs(i * b - a * j)
            require(1 <= det <= 4, 'nonzero determinant at most4')
            dets[det] += 1
            rows.append({'z': z, 'interval': [enc(left), enc(right)], 'pairs': [list(x) for x in pair], 'abs_determinant': det})
    radius4 = []
    for z in range(42):
        if gcd(z, 42) == 1:
            i, j = (1, 0)
        else:
            opts = (1, 2, 3) if z % 2 == 0 else (2, 4)
            i, j = next(((i, -1) for i in opts if gcd(i * z - 1, 42) == 1))
        require(gcd(i * z + j, 42) == 1 and 1 <= i <= 4, 'radius4 core unit')
        require(max(abs(F(j)), abs(j + F(i, 2))) <= 1, 'radius4 full positive half')
        radius4.append({'z': z, 'pair': [i, j]})
    return {'positive_closed_intervals': 3, 'positive_cells': 126, 'negative_cells_by_reflection': 126, 'groups': 14, 'max_abs_determinant': max(dets), 'determinant_histogram': dict(sorted(dets.items())), 'rows': rows, 'both_bad_radius4_rows': radius4}

def independent_interval_cover():
    """Independent brute candidate union, without using TABLE/GROUPS.

    Each candidate's exact shear interval is intersected with [-1/2,1/2].
    Union coverage, including all endpoints, is checked for all residues modulo
    42*r for eight primes. The infinite-prime proof uses determinants above.
    """
    total, lengths = (0, {})
    for prime in (11, 13, 17, 19, 23, 29, 31, 37):
        q = 42 * prime
        for z in range(q):
            intervals = []
            for i in range(1, 7):
                for j in range(-4, 5):
                    if gcd(i * z + j, q) == 1:
                        a = max(F(-1, 2), F(-1 - j, i))
                        b = min(F(1, 2), F(1 - j, i))
                        if a <= b:
                            intervals.append((a, b))
            end = F(-1, 2)
            for a, b in sorted(intervals):
                if b < end:
                    continue
                require(a <= end, 'independent interval no gap')
                end = max(end, b)
            require(end == F(1, 2), 'independent full shear coverage')
            total += 1
        lengths[str(prime)] = q
    return {'method': 'union of exact intervals from all admissible coefficient pairs', 'primes': [11, 13, 17, 19, 23, 29, 31, 37], 'residues_by_prime': lengths, 'total_residue_classes': total, 'universal_prime_proof': False}

def obstructions():
    rows = []
    for radius, q, z in ((4, 462, 275), (4, 546, 65), (5, 462, 259)):
        tau = F(1, 4)
        pairs = [(i, j, gcd(i * z + j, q)) for i in range(1, radius + 1) for j in range(-radius, 2) if abs(j + i * tau) <= 1]
        require(all((g > 1 for i, j, g in pairs)), 'declared radius obstruction')
        rows.append({'radius': radius, 'preserved_mark': q, 'z': z, 'tau': enc(tau), 'eligible_pairs_and_gcd': [list(x) for x in pairs]})
    ca = [(i, j) for i in range(1, 6) for j in range(-5, 2) if abs(F(j) + F(i, 4)) <= 1 and gcd(i * 7 + j, 42) == 1]
    require(ca == [(2, -1)], 'symbolic radius5 obstruction sole candidate')
    return {'scope': 'coefficient ansatz only; not a physical shortest-vector counterexample', 'examples': rows, 'all_primes_argument': 'z=7 mod42 and 2z=1 mod r; CRT applies for every prime r>=11'}

def periodic_certificate():
    rows = []
    for j in range(42):
        units = [e for e in range(1, j + 1) if gcd(e, 42) == 1]
        U, S = (len(units), sum(units))
        mob = {1: 1, 2: -1, 3: -1, 7: -1, 6: 1, 14: 1, 21: 1, 42: -1}
        U2 = sum((mu * (j // d) for d, mu in mob.items()))
        S2 = sum((mu * d * (j // d) * (j // d + 1) // 2 for d, mu in mob.items()))
        require((U, S) == (U2, S2), 'periodic inclusion exclusion')
        eps = U - F(2, 7) * j
        D = 2 * j * U + S - A * j * j
        require(3 * eps <= C, 'periodic linear slope')
        require(D <= C * j, 'periodic initial intercept')
        rows.append({'j': j, 'U': U, 'S': S, 'epsilon': enc(eps), 'D': enc(D)})
    require(F(rows[31]['D']) / 31 == C, 'periodic sharp endpoint31')
    return {'mark': 42, 'rho': '2/7', 'A': enc(A), 'C': enc(C), 'sharp_endpoint_for_fixed_A': 31, 'rows': rows}

def integral_certificate():
    lo, hi = atan_interval(F(3, 4), 40)
    L, H = log_interval(F(5, 2))
    k_lo, k_hi = (lo - H / 2, hi - L / 2)
    require(0 < k_lo <= k_hi < K_UP, 'certified K upper')
    a0, a = (F(866025, 10 ** 6), F(866026, 10 ** 6))
    b0, b = (F(1154700, 10 ** 6), F(1154701, 10 ** 6))
    require(a0 * a0 < F(3, 4) < a * a and a < 1, 'J endpoint a')
    require(b0 * b0 < F(4, 3) < b * b and b * b < 2, 'J endpoint b')

    def first(x, n):
        return sum(((-1) ** j * x ** (4 * j + 1) / F(4 * j + 1) for j in range(n + 1)), F(0))

    def second(x, n):
        return sum(((-1) ** j * x ** (2 * j + 1) / F(2 ** (j + 1) * (2 * j + 1)) for j in range(n + 1)), F(0))
    jl = 4 * (first(a0, 21) + second(b0, 41))
    jh = 4 * (first(a, 20) + second(b, 40))
    require(0 < jl < jh < J4_UP, 'certified J4 upper')
    require(F(3, 2) * J4_UP * J4_UP < J6_UP * J6_UP, 'J6 scaled upper25over4')
    require(6 * A * K_UP == F(139017, 175000), 'joint periodic linear coefficient')
    require(2 * C * J6_UP == F(18825, 434), 'joint root coefficient')
    require(1 + F(4, 3) ** 2 == F(25, 9), 'G boundary first logarithm')
    require(1 + F(3, 2) == F(5, 2), 'G boundary Phi cancellation')
    return {'K_interval': [enc(downward(k_lo)), enc(upward(k_hi))], 'K_upper': enc(K_UP), 'J4_interval': [enc(downward(jl)), enc(upward(jh))], 'J4_upper': enc(J4_UP), 'J6_upper': enc(J6_UP), 'periodic_Phi_coefficient': '15/7', 'periodic_lambda_coefficient': '139017/175000', 'periodic_Msqrt_coefficient': '18825/434', 'general_G_branches': 'Phi(y/2)+yK for y<=4/3; log(1+y*y)/2+y*atan(1/y) for y>=4/3'}

def minimum(M, s, q):
    a = s if 2 * s <= M else s - M
    best = (max(abs(a), 1), a, 1)
    d = 2
    while d <= best[0]:
        if gcd(d, q) == 1:
            a = s * d % M
            if 2 * a > M:
                a -= M
            h = max(abs(a), d)
            if h < best[0]:
                best = (h, a, d)
        d += 1
    return best

def construct(M, s, a, d, p, r):
    k = (a - s * d) // M
    require(a - s * d == M * k and gcd(k, d) == 1, 'primitive lattice coordinates')
    x = pow(k, -1, d)
    y = (1 - k * x) // d
    b, e = (s * x - M * y, x)
    den = a * a + d * d
    tau = F(a * b + d * e, den)
    shift = (tau + F(1, 2)).numerator // (tau + F(1, 2)).denominator
    b -= shift * a
    e -= shift * d
    tau -= shift
    require(a * e - d * b == M and abs(tau) <= F(1, 2), 'reduced companion determinant')
    require(gcd(d, e) == 1, 'companion denominator coprime')
    sign = 1 if tau >= 0 else -1
    z = sign * e * pow(d, -1, 42) % 42
    if d % p == 0 and d % r == 0:
        if gcd(z, 42) == 1:
            opts = [(1, 0)]
        else:
            inds = (1, 2, 3) if z % 2 == 0 else (2, 4)
            opts = [(i, -1) for i in inds if gcd(i * z - 1, 42) == 1]
        radius = 4
    else:
        h = next((h for h, (a0, b0) in enumerate(INTERVALS) if a0 <= abs(tau) <= b0))
        opts = TABLE[z][h]
        radius = 6
    for i, j0 in opts:
        j = sign * j0
        wa, wd = (j * a + i * b, j * d + i * e)
        if gcd(wd, 42 * p * r) == 1:
            if wd < 0:
                wa, wd = (-wa, -wd)
            require(wd > 0 and (wa - s * wd) % M == 0, 'repaired vector in lattice')
            require(abs(j + i * tau) <= 1 and 1 <= i <= radius, 'repaired vector coefficients')
            return (max(abs(wa), wd), radius)
    raise AssertionError('no repaired vector')

def physical_regressions():
    cases = []
    for p, r, M in ((11, 13, 143), (11, 13, 286), (11, 13, 429), (11, 13, 1573), (11, 13, 1859), (11, 13, 3146), (11, 13, 6006), (11, 17, 187), (13, 17, 221)):
        require(M % (p * r) == 0 and p != r, 'regression theorem domain')
        bad = Counter()
        both = 0
        zeros = 0
        worst_radius = 0
        for s in range(M):
            H, a, d = minimum(M, s, 42)
            new, an, dn = minimum(M, s, 42 * p * r)
            require(H <= new, 'marked minimum monotonicity')
            if d % p and d % r:
                require(new == H, 'no bad prime no cost')
                continue
            bad[H] += 1
            zeros += a == 0
            both += d % (p * r) == 0
            witness, R = construct(M, s, a, d, p, r)
            worst_radius = max(worst_radius, R)
            require(new <= witness, 'actual minimum below repair witness')
            if 4 * H * H <= 3 * R * M:
                require(witness * witness * H * H <= H ** 4 + R * R * M * M, 'pointwise low-X geometric bound')
            else:
                require(2 * witness * H <= 2 * H * H + R * M, 'pointwise high-X geometric bound')
        cumulative = 0
        lam = F(1, p) + F(1, r)
        for T in range(max(bad, default=0) + 1):
            cumulative += bad[T]
            require(cumulative <= A * lam * T * T + 2 * C * T, 'observed joint periodic CDF')
            require(cumulative <= 3 * lam * T * T, 'observed joint free CDF')
        cases.append({'p': p, 'r': r, 'M': M, 'states': M, 'bad_selected_minima': sum(bad.values()), 'both_bad': both, 'zero_numerator_bad_minima': zeros, 'max_repair_radius_used': worst_radius, 'bad_height_histogram': dict(sorted(bad.items()))})
    return {'method': 'all states, exact denominator enumeration plus constructive lattice repair', 'cases': cases, 'states': sum((x['states'] for x in cases)), 'does_not_replace_universal_proof': True}

def source_and_gate_certificate():
    receipt = json.loads((ROOT / 'imports/SOURCE_RECEIPT.json').read_text())
    require(receipt['status'] == 'PASS_SELECTIVE_SOURCE_BYTES', 'source receipt status')
    expected_names = {'PROFILE42_PROOF_SOURCE.md', 'q42_k1.json', 'q42_k2.json', 'small_prime_closure.json', 'PROFILE42_MANIFEST_SOURCE.json', 'ROOT42_MANIFEST_SOURCE.json', 'INTEGRAL_CONSTANTS_SOURCE.md'}
    copied = [r['copied_as'] for r in receipt['sources']]
    require(len(copied) == len(set(copied)) == len(expected_names) and set(copied) == {'imports/' + n for n in expected_names}, 'source receipt exact unique inventory')
    require(len({r['source'] for r in receipt['sources']}) == len(expected_names), 'source receipt original paths unique')
    require({p.name for p in (ROOT / 'imports').iterdir()} == expected_names | {'SOURCE_RECEIPT.json'}, 'source directory exact inventory')
    for row in receipt['sources']:
        p = ROOT / row['copied_as']
        raw = p.read_bytes()
        require(hashlib.sha256(raw).hexdigest() == row['sha256'] and len(raw) == row['bytes'], 'imported source byte identity')
    mf = json.loads((ROOT / 'imports/PROFILE42_MANIFEST_SOURCE.json').read_text())
    root_manifest_bytes = (ROOT / 'imports/ROOT42_MANIFEST_SOURCE.json').read_bytes()
    require(hashlib.sha256(root_manifest_bytes).hexdigest() == '1340e72a05355e39a767212c2bae9e464be16c8f9a0bcea466e393f0baf21409', 'historical root42 manifest expected identity')
    root_manifest = json.loads(root_manifest_bytes)
    root_files = {row['path']: row for row in root_manifest['files']}
    require(root_manifest['algorithm'] == 'SHA-256' and len(root_files) == len(root_manifest['files']), 'historical root manifest schema and unique inventory')
    links = {'PROFILE42_PROOF_SOURCE.md': 'profiles/PROOF.md', 'q42_k1.json': 'profiles/results/q42_k1.json', 'q42_k2.json': 'profiles/results/q42_k2.json', 'small_prime_closure.json': 'profiles/results/small_prime_closure.json', 'PROFILE42_MANIFEST_SOURCE.json': 'profiles/MANIFEST_SHA256.json'}
    for name, old_path in links.items():
        raw = (ROOT / 'imports' / name).read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        require(root_files[old_path]['sha256'] == sha and root_files[old_path]['bytes'] == len(raw), 'selected anchor source linked to frozen root42 manifest')
        if name != 'PROFILE42_MANIFEST_SOURCE.json':
            require(sha == mf[old_path.removeprefix('profiles/')], 'historical profile manifest link')
    a1 = F(json.loads((ROOT / 'imports/q42_k1.json').read_text())['total_upper'])
    a2 = F(json.loads((ROOT / 'imports/q42_k2.json').read_text())['total_upper'])
    require(a1 < F(107, 500), 'historical anchor42 below107over500')
    require(0 < a2 < F(67, 1000), 'historical anchor42square positive below67over1000')
    gates = []
    for P, N, anchor, K in ((229, 1, F(107, 500), 42), (83, 2, F(67, 1000), 1764)):
        B = 42 * P * P
        M = B ** N
        cost = periodic_upper(P, P, M)
        ub = (anchor + cost) / log_interval(B)[0]
        q = upward(ub)
        require(q < DELTA, 'new infinite joint prime gate')
        require(M % K == 0, 'gate anchor divisor')
        gates.append({'p_and_r_minimum': P, 'all_N_at_least': N, 'all_positive_valuations': True, 'B_minimum': B, 'M_minimum': M, 'anchor_divisor': K, 'anchor_upper': enc(anchor), 'insertion_cost_upper': enc(upward(cost)), 'normalized_upper': enc(q), 'normalized_upper_decimal': decimal_upper(q), 'gap_to_2over85_lower': enc(downward(DELTA - q)), 'uses_distinctness_to_tighten_scalar_boundary': False})
    return {'anchor_status': 'inherited proven marked42 profiles; no large-stream replay in this module', 'anchor_exact_upper_42': enc(a1), 'anchor_exact_upper_1764': enc(a2), 'anchor_1764_is_negative': False, 'gates': gates, 'full_family_closed': False}

def build():
    CHECKS.clear()
    out = {'schema': 'P72_JOINT42PR_RADIUS6_V1', 'status': 'PASS_UNIFORM_REPAIR_AND_TWO_INFINITE_GATES_WITH_DECLARED_HISTORICAL_ANCHORS', 'theorem_domain': 'distinct primes p,r>=11 dividing M; no requirement42|M for insertion; physical specialization rad(B)=42*p*r,M=B**N', 'radius6': residue_certificate(), 'independent_finite_prime_cover': independent_interval_cover(), 'radius_obstructions': obstructions(), 'periodic_count': periodic_certificate(), 'integrals': integral_certificate(), 'physical_regressions': physical_regressions(), 'infinite_gates': source_and_gate_certificate()}
    out['checks_by_type'] = dict(sorted(CHECKS.items()))
    out['exact_assertions'] = sum(CHECKS.values())
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--report', required=True, type=Path)
    ap.add_argument('--write-certificate', action='store_true')
    args = ap.parse_args()
    if args.report.resolve().is_relative_to(ROOT):
        raise ValueError('report path must be outside the immutable component')
    start = time.monotonic()
    out = build()
    payload = json.dumps(out, indent=2, sort_keys=True) + '\n'
    certificate = ROOT / 'certificate.json'
    if args.write_certificate:
        certificate.write_text(payload, newline='\n')
    else:
        require(certificate.read_text() == payload, 'stored certificate exact byte replay')
    report = {'status': 'PASS', 'certificate_sha256': hashlib.sha256(certificate.read_bytes()).hexdigest(), 'exact_assertions': out['exact_assertions'], 'physical_states': out['physical_regressions']['states'], 'positive_residue_shear_cells': 126, 'independent_finite_prime_residue_classes': out['independent_finite_prime_cover']['total_residue_classes'], 'infinite_gates': [{'p_and_r_minimum': g['p_and_r_minimum'], 'N_minimum': g['all_N_at_least'], 'upper': g['normalized_upper_decimal']} for g in out['infinite_gates']['gates']], 'elapsed_seconds': round(time.monotonic() - start, 6), 'scope': 'uniform radius6 theorem; two corner exclusions; full42pr family remains open; imported42 anchors not geometrically rerun'}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + '\n', newline='\n')
    print(json.dumps(report))
if __name__ == '__main__':
    main()
