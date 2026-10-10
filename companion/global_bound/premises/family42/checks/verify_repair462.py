"""Finite certificate for PE462 and a five-prime tail, no external packages.

Imports the basic rational log/definition oracle from verify_repair42.py.
The q=462 marked-profile anchor is a frozen input from the parallel profile
lane; checking its rational value does not replay the geometric producer.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path
from hashlib import sha256
import json
from verify_repair42 import need, fs, phi_bounds, log_bounds, oracle, bezout
ROOT = Path(__file__).resolve().parent
Q = 462
R = 6

def make_table(q, r):
    table = {}
    for z in range(q):
        candidates = []
        for i in range(1, r + 1):
            for j in range(-r, 2):
                if gcd(j + i * z, q) != 1:
                    continue
                lo = max(F(0), F(-1 - j, i))
                hi = min(F(1, 2), F(1 - j, i))
                if lo <= hi:
                    candidates.append((lo, hi, i, j))
        reach = F(0)
        intervals = []
        while reach < F(1, 2):
            allowed = [x for x in candidates if x[0] <= reach and x[1] > reach]
            need(bool(allowed), f'continuous coverage for q={q}, z={z}')
            row = max(allowed, key=lambda x: (x[1], -x[2], -abs(x[3])))
            intervals.append(row)
            reach = row[1]
        table[z] = intervals
    return table

def main():
    table = make_table(Q, R)
    total_intervals = 0
    serialized_table = []
    for z, intervals in table.items():
        previous = F(0)
        serialized_rows = []
        for lo, hi, i, j in intervals:
            need(lo <= previous < hi, 'no uncovered interval')
            need(0 <= lo <= hi <= F(1, 2), 'closed half-interval domain')
            need(gcd(j + i * z, Q) == 1 and 1 <= i <= R, 'repair mark and radius')
            need(max(abs(j + i * lo), abs(j + i * hi)) <= 1, 'continuous affine bound')
            previous = hi
            serialized_rows.append({'lo': fs(lo), 'hi': fs(hi), 'i': i, 'j': j})
        need(previous == F(1, 2), 'right endpoint included')
        total_intervals += len(intervals)
        serialized_table.append({'z': z, 'intervals': serialized_rows})
    units = [j for j in range(1, Q) if gcd(j, Q) == 1]
    need(len(units) == 120 and sum(units) == 27720, 'q462 unit period')
    u = F(20, 77)
    Aq = F(50, 77)
    Cq = F(492, 77)
    U = S = 0
    periodic_rows = []
    maxC = F(0)
    for j in range(Q):
        if j and gcd(j, Q) == 1:
            U += 1
            S += j
        eps = U - u * j
        D = 2 * j * U + S - Aq * j * j
        need(3 * eps <= Cq, 'P462 coefficient of n')
        need(D <= Cq * j, 'P462 n=0 endpoint')
        maxC = max(maxC, 3 * eps, D / j if j else F(0))
        periodic_rows.append({'j': j, 'U': U, 'S': S, 'epsilon': fs(eps), 'D': fs(D)})
    need(maxC == Cq, 'periodic C462 exact')
    need(3 * F(periodic_rows[115]['epsilon']['numerator'], periodic_rows[115]['epsilon']['denominator']) == Cq, 'C462 attained in slope j115')
    U = S = 0
    for A in range(4621):
        if A and gcd(A, Q) == 1:
            U += 1
            S += A
        need(2 * A * U + S <= Aq * A * A + Cq * A, 'direct periodic sanity check')
    restricted_negative = []
    z, tau = (245, F(1, 20))
    for i in range(1, 6):
        for j in range(-6, 3):
            if abs(j + i * tau) <= 1:
                g = gcd(j + i * z, Q)
                need(g > 1, 'restricted r<=5 fails at chosen point')
                restricted_negative.append({'i': i, 'j': j, 'gcd': g})
    model_rows = []
    total_states = total_bad = 0
    for M, p in [(6006, 13), (7854, 17), (169, 13)]:
        marked = oracle(M, Q)
        full = oracle(M, Q * p)
        bad = 0
        for s, (H, a, d) in enumerate(marked):
            need(H <= full[s][0], 'marked relaxation')
            if d % p:
                need(H == full[s][0], 'already full admissible')
                continue
            bad += 1
            need(a % p == 0, 'bad numerator')
            k = (a - s * d) // M
            g, x, y = bezout(k, d)
            need(g == 1, 'lattice primitivity')
            b, e = (s * x - M * y, x)
            need(a * e - d * b == M, 'basis determinant')
            norm2 = a * a + d * d
            tau = F(a * b + d * e, norm2)
            t = tau + F(1, 2)
            shift = t.numerator // t.denominator
            b -= shift * a
            e -= shift * d
            tau -= shift
            need(abs(tau) <= F(1, 2), 'shear range')
            z = e * pow(d, -1, Q) % Q
            if tau < 0:
                lookup_z, lookup_tau = (-z % Q, -tau)
            else:
                lookup_z, lookup_tau = (z, tau)
            applicable = [row for row in table[lookup_z] if row[0] <= lookup_tau <= row[1]]
            need(bool(applicable), 'table branch present at boundary')
            _, _, i, j = applicable[0]
            if tau < 0:
                j = -j
            aa, dd = (j * a + i * b, j * d + i * e)
            if dd < 0:
                aa, dd = (-aa, -dd)
            need(dd > 0 and gcd(dd, Q * p) == 1, 'fully unit positive denominator')
            need((aa - s * dd) % M == 0, 'congruence')
            Hr = max(abs(aa), dd)
            need(F(Hr) <= H + F(R * M, H), 'R462 height bound')
            need(full[s][0] <= Hr, 'definition minimum below repaired witness')
        model_rows.append({'M': M, 'p': p, 'states': M, 'bad_repaired': bad})
        total_states += M
        total_bad += bad
    anchor_path = ROOT / 'imports' / 'q462_anchor_33nodes.json'
    anchor_bytes = anchor_path.read_bytes()
    anchor_source = json.loads(anchor_bytes)
    need(anchor_source['q'] == 462 and anchor_source['k'] == 1, 'imported anchor mark and horizon')
    anchor_exact = F(anchor_source['total_upper'])
    need(anchor_exact < F(17, 100), 'source anchor below rational cap')
    need(F(49, 20) ** 2 > 6, 'sqrt6 strictly below 49/20')
    need(F(22, 7) * F(49, 20) * Cq == F(246, 5), 'rational PE462 coefficient')
    gates = []
    for p in [650, 700]:
        B = M = Q * p
        root = isqrt(M)
        need(root * root <= M < (root + 1) * (root + 1), 'exact integer-root lower bound')
        _, phi_hi = phi_bounds(F(300, 77 * p))
        loglo, _ = log_bounds(B)
        upper = (F(17, 100) + phi_hi + F(246, 5 * root)) / loglo
        need(upper < F(2, 85), 'five-prime tail gate')
        gates.append({'p_scalar': p, 'B_min': B, 'M_min': M, 'sqrt_lower': root, 'upper': fs(upper), 'upper_descriptive': float(upper)})
    report = {'status': 'PASS', 'q': Q, 'repair_radius': R, 'positive_half_intervals': total_intervals, 'reflected_closed_intervals': 2 * total_intervals, 'periodic_Aq': fs(Aq), 'periodic_Cq': fs(Cq), 'PE_argument': '300/(77p)', 'PE_sqrtM_coefficient': fs(F(246, 5)), 'repair_table': serialized_table, 'periodic_rows': periodic_rows, 'direct_P462_A_range': [0, 4620], 'model_states': total_states, 'model_bad_repaired': total_bad, 'model_rows': model_rows, 'restricted_negative': {'z': 245, 'tau': fs(F(1, 20)), 'radius_max': 5, 'pairs': restricted_negative, 'scope': 'this longitudinal-bound repair class only'}, 'anchor_import': {'file': str(anchor_path.relative_to(ROOT)), 'sha256': sha256(anchor_bytes).hexdigest(), 'source_total_upper': fs(anchor_exact), 'cap': fs(F(17, 100)), 'geometric_producer_replayed_here': False}, 'gates': gates, 'claim_scope': 'general insertion proof + finite table + separate profile-anchor dependency'}
    outpath = ROOT / 'repair462_certificate.json'
    outpath.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': 'PASS', 'certificate': outpath.name, 'closed_intervals': 2 * total_intervals, 'model_states': total_states, 'bad_repaired': total_bad, 'gate_uppers': {str(g['p_scalar']): g['upper_descriptive'] for g in gates}}, indent=2))
if __name__ == '__main__':
    main()
