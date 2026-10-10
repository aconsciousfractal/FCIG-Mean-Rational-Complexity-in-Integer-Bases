"""Exact radius-six insertion at mark770 and the five-prime scalar tail300.

The imported q770 marked-profile bound applies to every multiple of770.
Reading its rational summary is not a replay of the geometric certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json
from verify_repair70 import need, fs, oracle, bezout, periodic_certificate, log_bounds, phi_bounds
from verify_gates70 import decimal_ceiling, decimal_floor, prime
ROOT = Path(__file__).resolve().parent
Q = 770
R = 6

def make_table(q, radius):
    table = {}
    for z in range(q):
        candidates = []
        for i in range(1, radius + 1):
            for j in range(-radius, 2):
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
            need(bool(allowed), f'continuous table coverage q{q},z{z}')
            row = max(allowed, key=lambda x: (x[1], -x[2], -abs(x[3])))
            intervals.append(row)
            reach = row[1]
        table[z] = intervals
    return table

def main():
    table = make_table(Q, R)
    serialized_table = []
    total_intervals = 0
    special = []
    for z, intervals in table.items():
        previous = F(0)
        rows = []
        for lo, hi, i, j in intervals:
            need(lo <= previous < hi, 'no open coverage gap')
            need(0 <= lo <= hi <= F(1, 2), 'closed table domain')
            need(gcd(j + i * z, Q) == 1 and 1 <= i <= R, 'old-mark unit and radius')
            need(max(abs(j + i * lo), abs(j + i * hi)) <= 1, 'affine interval bound')
            need(gcd(-j + i * -z, Q) == 1, 'reflected old-mark unit')
            need(max(abs(-j - i * lo), abs(-j - i * hi)) <= 1, 'reflected affine interval bound')
            previous = hi
            rows.append({'lo': fs(lo), 'hi': fs(hi), 'i': i, 'j': j})
        need(previous == F(1, 2), 'right endpoint covered')
        total_intervals += len(intervals)
        serialized_table.append({'z': z, 'intervals': rows})
        if len(intervals) > 1:
            special.append(z)
    need(total_intervals == 776, 'number of half-intervals')
    need(special == [25, 259, 319, 553, 583, 765], 'special residue classes')
    negative = []
    z, tau = (25, F(1, 10))
    for i in range(1, 6):
        for j in range(-6, 3):
            if abs(j + i * tau) <= 1:
                g = gcd(j + i * z, Q)
                need(g > 1, 'radius-five restricted obstruction')
                negative.append({'i': i, 'j': j, 'gcd': g})
    need(len(negative) == 10, 'all possible obstruction pairs')
    rho, A, C, periods = periodic_certificate(Q)
    need(rho == F(24, 77) and A == F(60, 77) and (C == F(453, 77)), 'reconstructed periodic coefficients770')
    slope = periods[183]['slope']
    need(F(slope['numerator'], slope['denominator']) == C, 'periodic slope attains C770 at j183')
    need(sum((j for j in range(Q) if gcd(j, Q) == 1)) == 92400, 'unit residue sum770')
    U = S = 0
    for T in range(7701):
        if T and gcd(T, Q) == 1:
            U += 1
            S += T
        need(2 * T * U + S <= A * T * T + C * T, 'direct periodic770 regression')
    model_rows = []
    states = bad_total = density_points = 0
    for M, p in [(10010, 13), (13090, 17), (169, 13)]:
        need(p > R and Q % p and (M % p == 0), 'insertion hypotheses770')
        marked, full = (oracle(M, Q), oracle(M, Q * p))
        bad = 0
        heights = {}
        for s, (H, a, d) in enumerate(marked):
            need(H <= full[s][0], 'marked relaxation770')
            if d % p:
                need(H == full[s][0], 'already full admissible770')
                continue
            bad += 1
            heights[H] = heights.get(H, 0) + 1
            need(a % p == 0, 'bad numerator770')
            k = (a - s * d) // M
            g, x, y = bezout(k, d)
            need(g == 1, 'primitive lattice coordinates770')
            b, e = (s * x - M * y, x)
            need(a * e - d * b == M, 'basis determinant770')
            need(gcd(d, e) == 1 and e % p, 'second denominator excludes p770')
            norm2 = a * a + d * d
            tau = F(a * b + d * e, norm2)
            shifted = tau + F(1, 2)
            shift = shifted.numerator // shifted.denominator
            b -= shift * a
            e -= shift * d
            tau -= shift
            need(abs(tau) <= F(1, 2), 'reduced shear770')
            z = e * pow(d, -1, Q) % Q
            zz, tt = (-z % Q, -tau) if tau < 0 else (z, tau)
            matches = [row for row in table[zz] if row[0] <= tt <= row[1]]
            need(bool(matches), 'chosen shear covered770')
            _, _, i, j = matches[0]
            if tau < 0:
                j = -j
            aa, dd = (j * a + i * b, j * d + i * e)
            if dd < 0:
                aa, dd = (-aa, -dd)
            need(dd > 0 and gcd(dd, Q * p) == 1, 'positive fully unit denominator770')
            need((aa - s * dd) % M == 0, 'lattice membership770')
            Hr = max(abs(aa), dd)
            need(F(Hr) <= H + F(R * M, H), 'radius-six height bound770')
            need(full[s][0] <= Hr, 'witness above full minimum770')
        cumulative = 0
        for H in sorted(heights):
            cumulative += heights[H]
            need(cumulative <= min(F(M), A * H * H / p + C * H), 'D770 finite density')
            density_points += 1
        states += M
        bad_total += bad
        model_rows.append({'M': M, 'p': p, 'states': M, 'bad_repaired': bad, 'density_breakpoints': len(heights)})
    anchor_path = ROOT / 'imports' / 'q770.json'
    anchor_bytes = anchor_path.read_bytes()
    anchor_source = json.loads(anchor_bytes)
    need(anchor_source['q'] == 770 and anchor_source['k'] == 1 and (anchor_source['transfer_index'] == 770), 'marked770 import contract')
    anchor_exact = F(anchor_source['total_upper'])
    cap = F(9, 80)
    need(anchor_exact < cap, 'exact cap770 below9/80')
    need(F(49, 20) ** 2 > 6, 'sqrt6 upper by integer squares')
    need(F(22, 7) * F(49, 20) == F(77, 10), 'pi sqrt6 rational coefficient')
    need(C * F(77, 10) == F(453, 10), 'PE770 rational coefficient')
    p = 300
    B = M = Q * p
    root = isqrt(M)
    need(root * root <= M < (root + 1) * (root + 1), 'gate integer square root770')
    arg = R * A / p
    need(arg == F(360, 77 * p), 'PE770 argument')
    phi_lo, phi_hi = phi_bounds(arg)
    loglo, loghi = log_bounds(B)
    upper = (cap + phi_hi + F(453, 10 * root)) / loglo
    margin = F(2, 85) - upper
    need(margin > 0, 'all-exponent all-length five-prime tail300')
    neighbors = [v for v in range(293, 308) if prime(v)]
    need(neighbors == [293, 307], 'prime gap around scalar300')
    output = {'status': 'PASS', 'q': Q, 'repair_radius': R, 'positive_half_intervals': total_intervals, 'reflected_closed_intervals': 2 * total_intervals, 'special_classes': special, 'repair_table': serialized_table, 'rho': fs(rho), 'periodic_Aq': fs(A), 'periodic_Cq': fs(C), 'periodic_rows': periods, 'direct_periodic_range': [0, 7700], 'PE_argument': '360/(77*ell)', 'PE_sqrtM_coefficient': fs(F(453, 10)), 'restricted_radius_five_obstruction': {'z': 25, 'tau': fs(F(1, 10)), 'pairs': negative, 'scope': 'only the ansatz abs(j+i*tau)<=1; not global geometric optimality'}, 'model_rows': model_rows, 'model_instances': states, 'bad_minima_repaired': bad_total, 'density_breakpoints_checked': density_points, 'anchor_import': {'file': anchor_path.relative_to(ROOT).as_posix(), 'bytes': len(anchor_bytes), 'sha256': sha256(anchor_bytes).hexdigest(), 'total_upper': fs(anchor_exact), 'cap': fs(cap), 'scope': 'nu770(M) for every multiple M of770, including additional prime divisors', 'geometric_producer_replayed_here': False}, 'gate': {'ell_scalar': p, 'first_actual_prime': 307, 'B_min': B, 'M_min': M, 'sqrt_lower': root, 'Phi_argument': fs(arg), 'Phi_lower': fs(phi_lo), 'Phi_upper': fs(phi_hi), 'lnB_lower': fs(loglo), 'lnB_upper': fs(loghi), 'upper': fs(upper), 'upper_decimal_ceiling': decimal_ceiling(upper), 'margin_from_delta': fs(margin), 'margin_decimal_floor': decimal_floor(margin), 'prime_neighbors': neighbors, 'domain': 'rad(B)=770*ell, ell>=300 prime; all positive valuations and N>=1', 'upper_is_value_of_Delta': False}, 'claim_scope': 'analytic insertion plus finite table/period certificates and separate marked-profile input', 'all_five_prime_cases_closed': False}
    path = ROOT / 'repair770_certificate.json'
    path.write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': 'PASS', 'certificate': path.name, 'closed_intervals': 2 * total_intervals, 'model_instances': states, 'bad_minima_repaired': bad_total, 'density_breakpoints_checked': density_points, 'tail_upper': decimal_ceiling(upper), 'margin_floor': decimal_floor(margin)}, indent=2))
if __name__ == '__main__':
    main()
