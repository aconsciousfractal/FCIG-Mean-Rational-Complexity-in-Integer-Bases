"""Exact finite certificates supporting R42, P42, and PE42.

The infinite conclusions depend on the proofs in REPAIR42_PROOF.md; these
finite checks are not presented as a proof by sampling. No third-party package.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
import json
import sys
ROOT = Path(__file__).resolve().parent
Q = 42
GROUPS = [([1, 5, 11, 13, 17, 19, 23, 25, 29, 31, 37, 41], 1, 0), ([0, 2, 6, 12, 14, 18, 20, 24, 26, 30, 32, 38], 1, -1), ([3, 7, 9, 10, 15, 16, 21, 22, 27, 28, 33, 34, 36, 40], 2, -1), ([4, 8], 3, -1), ([35, 39], 4, -1)]

def need(test, label):
    if not test:
        raise RuntimeError(label)

def fs(x):
    x = F(x)
    return {'numerator': x.numerator, 'denominator': x.denominator}

def log_bounds(x, terms=80):
    """Rational enclosure for natural log by binary reduction + atanh."""
    x = F(x)
    need(x > 0, 'positive logarithm')
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def atanh_bounds(z):
        z = F(z)
        zz = z * z
        power = z
        total = F(0)
        for j in range(terms):
            total += 2 * power / (2 * j + 1)
            power *= zz
        tail = 2 * power / ((2 * terms + 1) * (1 - zz))
        return (total, total + tail)
    lo, hi = atanh_bounds((x - 1) / (x + 1))
    l2, h2 = atanh_bounds(F(1, 3))
    return (lo + k * l2, hi + k * h2) if k >= 0 else (lo + k * h2, hi + k * l2)

def phi_bounds(x):
    a, b = log_bounds(1 + x)
    c, d = log_bounds(1 + 1 / x)
    return (a + x * c, b + x * d)

def oracle(M, q):
    """Definition directly: all d <= M//2 (d=1 already bounds height)."""
    cap = max(1, M // 2)
    out = []
    ds = [d for d in range(1, cap + 1) if gcd(d, q) == 1]
    for s in range(M):
        best = None
        for d in ds:
            if best is not None and d > best[0]:
                break
            a = s * d % M
            if 2 * a > M:
                a -= M
            candidate = (max(abs(a), d), d, a)
            if best is None or candidate < best:
                best = candidate
        need(best is not None, 'denominator 1 exists')
        H, d, a = best
        out.append((H, a, d))
    return out

def bezout(a, b):
    old_r, r, old_s, s, old_t, t = (a, b, 1, 0, 0, 1)
    while r:
        c = old_r // r
        old_r, r = (r, old_r - c * r)
        old_s, s = (s, old_s - c * s)
        old_t, t = (t, old_t - c * t)
    if old_r < 0:
        return (-old_r, -old_s, -old_t)
    return (old_r, old_s, old_t)

def main():
    table = {}
    for zs, i, j in GROUPS:
        for z in zs:
            need(z not in table, 'no repeated repair class')
            need(gcd(j + i * z, Q) == 1, 'q-unit witness')
            need(max(abs(F(j)), abs(F(j) + F(i, 2))) <= 1, 'closed half-interval')
            table[z] = (i, j)
    need(set(table) == set(range(Q)), 'all 42 classes')
    obstructions = []
    for z in (35, 39):
        tau = F(1, 4)
        inspected = []
        for i in range(1, 4):
            for j in range(-4, 5):
                if abs(j + i * tau) <= 1:
                    g = gcd(j + i * z, Q)
                    inspected.append({'i': i, 'j': j, 'gcd': g})
                    need(g > 1, 'restricted i<=3 obstruction')
        obstructions.append({'z': z, 'tau': fs(tau), 'pairs': inspected})
    u = F(2, 7)
    c = F(753, 217)
    U = S = 0
    period = []
    max_seen = F(0)
    for j in range(Q):
        if j and gcd(j, Q) == 1:
            U += 1
            S += j
        eps = U - u * j
        D = 2 * j * U + S - F(5, 7) * j * j
        need(3 * eps <= c, 'P42 coefficient of n')
        need(D <= c * j, 'P42 coefficient at n=0')
        if j:
            max_seen = max(max_seen, 3 * eps, D / j)
        period.append({'j': j, 'U': U, 'S': S, 'epsilon': fs(eps), 'D': fs(D), 'slope': fs(3 * eps)})
    need(max_seen == c, 'sharp periodic c at j=31')
    U = S = 0
    for A in range(4201):
        if A and gcd(A, Q) == 1:
            U += 1
            S += A
        need(2 * A * U + S <= F(5, 7) * A * A + c * A, 'direct P42 bank')
    pairs = [(11, 11), (22, 11), (55, 11), (77, 11), (110, 11), (121, 11), (143, 11), (154, 11), (242, 11), (462, 11), (13, 13), (39, 13), (143, 13), (169, 13), (286, 13), (17, 17), (187, 17), (289, 17), (19, 19), (361, 19)]
    states = bad = 0
    model_records = []
    for M, p in pairs:
        marked = oracle(M, Q)
        full = oracle(M, Q * p)
        bad_M = 0
        for s, (H, a, d) in enumerate(marked):
            Hf = full[s][0]
            need(H <= Hf, 'relaxation direction')
            if d % p:
                need(H == Hf, 'no repair for p-unit chosen minimizer')
                continue
            bad += 1
            bad_M += 1
            need(a % p == 0, 'bad numerator divisible by p')
            k = (a - s * d) // M
            g, x, y = bezout(k, d)
            need(g == 1, 'primitive in lattice coordinates')
            e = x
            b = s * x - M * y
            need(a * e - d * b == M, 'positive basis determinant')
            norm2 = a * a + d * d
            tau = F(a * b + d * e, norm2)
            shift = (tau + F(1, 2)).numerator // (tau + F(1, 2)).denominator
            b -= shift * a
            e -= shift * d
            tau -= shift
            need(abs(tau) <= F(1, 2), 'reduced shear')
            z = e * pow(d, -1, Q) % Q
            if tau >= 0:
                i, j = table[z]
            else:
                i, j = table[-z % Q]
                j = -j
            aa, dd = (j * a + i * b, j * d + i * e)
            if dd < 0:
                aa, dd = (-aa, -dd)
            need(dd > 0 and gcd(dd, Q * p) == 1, 'fully admissible repaired denominator')
            need((aa - s * dd) % M == 0, 'repaired congruence')
            Hr = max(abs(aa), dd)
            need(F(Hr) <= H + F(4 * M, H), 'R42 height bound')
            need(Hf <= Hr, 'full minimum below witness')
        states += M
        model_records.append({'M': M, 'p': p, 'states': M, 'bad': bad_M})
    example = (oracle(462, 42)[42], oracle(462, 462)[42])
    need(example[0][0] == 11 and example[1][0] == 42, 'relaxation counterexample')
    need(-42 % 462 != 0, 'division by ambient gcd invalid')
    gates = []
    routes = [('p_ge_3400_all_N', 3400, 42 * 3400, 42 * 3400, F(107, 500)), ('p_ge_107_N_ge_2', 107, 42 * 107, (42 * 107) ** 2, F(67, 1000)), ('p_ge_107_B_ge_4000000', 107, 4000000, 4000000, F(107, 500)), ('p_ge_83_N1_core_exponents_ge2', 83, 42 * 42 * 83, 42 * 42 * 83, F(67, 1000))]
    for label, p, B, M, anchor in routes:
        root = isqrt(M)
        need(root * root <= M < (root + 1) * (root + 1), 'integer square-root enclosure')
        _, phi_hi = phi_bounds(F(20, 7 * p))
        lnB_lo, _ = log_bounds(B)
        upper = (anchor + phi_hi + F(33132, 1519 * root)) / lnB_lo
        need(upper < F(2, 85), 'strict certified gate: ' + label)
        gates.append({'route': label, 'p_scalar': p, 'B_min': B, 'M_min': M, 'anchor': fs(anchor), 'sqrt_lower': root, 'upper': fs(upper), 'upper_descriptive': float(upper)})
    p = 103
    phi_lo, _ = phi_bounds(F(20, 7 * p))
    _, lnB_hi = log_bounds(42 * p)
    lower_of_bound = (F(67, 1000) + phi_lo + F(33132, 1519 * 42 * p)) / lnB_hi
    need(lower_of_bound > F(2, 85), 'record limited failure of chosen sufficient bound')
    output = {'status': 'PASS', 'q': Q, 'repair_classes': 42, 'closed_half_intervals': 84, 'repair_coefficient': 4, 'periodic_coefficient': fs(c), 'quadratic_density_coefficient': fs(F(5, 7)), 'PE_argument': '20/(7p)', 'PE_sqrtM_coefficient': fs(F(33132, 1519)), 'periodic_rows': period, 'repair_table': [{'z': z, 'i': table[z][0], 'j': table[z][1]} for z in range(Q)], 'restricted_repair_obstructions': obstructions, 'direct_P42_A_range': [0, 4200], 'direct_model_states': states, 'direct_bad_minima_repaired': bad, 'model_records': model_records, 'relaxation_example': {'M': 462, 's': 42, 'marked': example[0], 'full': example[1]}, 'gates': gates, 'limited_negative': {'p': 103, 'N': 2, 'sufficient_bound_exceeds_delta': True, 'lower_of_bound': fs(lower_of_bound), 'interpretation': 'not a counterexample for the true mean'}, 'claim_scope': 'finite arithmetic certificate plus finite model checks; general proofs separate'}
    path = ROOT / 'repair42_certificate.json'
    path.write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': 'PASS', 'certificate': path.name, 'states': states, 'bad_minima_repaired': bad, 'gate_uppers': {x['route']: x['upper_descriptive'] for x in gates}}, indent=2))
if __name__ == '__main__':
    main()
