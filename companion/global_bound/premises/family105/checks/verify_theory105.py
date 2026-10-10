"""Exact tables/counting/penalties; finite model regressions stay separate."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
import exact105 as E
from exact105 import require, record, fraction, periodic, structural_table, direct_minima, log_interval
HERE = Path(__file__).resolve().parent

def verify_table(q, radius):
    table = structural_table(q, radius)
    out = []
    intervals = 0
    special = []
    for z, rows in table.items():
        reached = F(0)
        encoded = []
        for lo, hi, i, j in rows:
            require(0 <= lo <= reached < hi <= F(1, 2), 'continuous positive-half coverage')
            require(1 <= i <= radius and gcd(j + i * z, q) == 1, 'unit and radius')
            require(abs(j + i * lo) <= 1 and abs(j + i * hi) <= 1, 'affine endpoints')
            require(gcd(-j - i * z, q) == 1, 'reflected unit')
            require(abs(-j - i * lo) <= 1 and abs(-j - i * hi) <= 1, 'reflected endpoints')
            reached = hi
            intervals += 1
            encoded.append(dict(lo=record(lo), hi=record(hi), i=i, j=j))
        require(reached == F(1, 2), 'closed final endpoint')
        if len(rows) > 1:
            special.append(z)
        out.append(dict(z=z, intervals=encoded))
    return (table, out, intervals, special)

def regression(M, q, p, radius, table, A, C):
    require(M % p == 0 and q % p and (p > radius), 'model hypotheses')
    marked = direct_minima(M, q)
    full = direct_minima(M, q * p)
    heights = {}
    bad = 0
    zero = 0
    for s, (H, a, d) in enumerate(marked):
        require(H <= full[s][0], 'restriction direction')
        if d % p:
            require(H == full[s][0], 'chosen full-admissible minimum')
            continue
        bad += 1
        zero += a == 0
        heights[H] = heights.get(H, 0) + 1
        require(a % p == 0, 'bad numerator including zero')
        k = (a - s * d) // M
        require(gcd(k, d) == 1, 'lattice primitivity')
        x = pow(k, -1, d)
        y = (1 - k * x) // d
        b = s * x - M * y
        e = x
        require(a * e - d * b == M and gcd(d, e) == 1, 'oriented lattice basis')
        norm = a * a + d * d
        dot = a * b + d * e
        shift = (2 * dot + norm) // (2 * norm)
        b -= shift * a
        e -= shift * d
        tau = F(a * b + d * e, norm)
        require(abs(tau) <= F(1, 2), 'reduced shear')
        z = e * pow(d, -1, q) % q
        key = -z % q if tau < 0 else z
        t = abs(tau)
        choices = [r for r in table[key] if r[0] <= t <= r[1]]
        require(bool(choices), 'branch at closed boundary')
        _, _, i, j = choices[0]
        if tau < 0:
            j = -j
        aa = j * a + i * b
        dd = j * d + i * e
        if dd < 0:
            aa = -aa
            dd = -dd
        height = max(abs(aa), dd)
        require(dd > 0 and gcd(dd, q * p) == 1, 'repaired mark and sign')
        require((aa - s * dd) % M == 0 and full[s][0] <= height, 'repaired lattice and height')
        require(F(height) <= H + F(radius * M, H), 'triangular repair bound')
        eta = F(min(abs(a), d), H)
        v = F(radius * M, H * H)
        coordinate = max(1 + v * eta / (1 + eta * eta), eta + v / (1 + eta * eta))
        envelope = (1 + v / 2) ** 2 if v <= F(4, 3) else 1 + v * v
        require(F(height, H) <= coordinate and coordinate ** 2 <= envelope, 'sharp coordinate model')
    cumulative = 0
    for H in sorted(heights):
        cumulative += heights[H]
        require(cumulative <= min(F(M), A * H * H / p + C * H), 'periodic CDF breakpoint')
        require(cumulative <= min(F(M), F(3 * H * H, p)), 'M-free CDF breakpoint')
    return dict(M=M, q=q, p=p, states=M, bad=bad, zero_numerator_bad=zero, breakpoints=len(heights), p_squared_divides_M=M % (p * p) == 0, model_scope='finite regression, not universal quantifier')

def main():
    marks = {}
    tables = {}
    specs = [(105, 3, F(16, 35), F(8, 7), F(1550, 329), 25, F(1, 4), 105), (1155, 5, F(32, 77), F(80, 77), F(582, 77), 235, F(1, 8), 1161)]
    for q, radius, erho, eA, eC, z, tau, count in specs:
        table, encoded, intervals, special = verify_table(q, radius)
        tables[q] = table
        rho, A, C, rows = periodic(q)
        require((rho, A, C) == (erho, eA, eC), 'fresh exact periodic coefficients')
        require(intervals == count, 'positive interval count')
        smaller = []
        for i in range(1, radius):
            candidates = [j for j in range(-radius - 2, radius + 3) if abs(j + i * tau) <= 1]
            require(candidates == [-1, 0], 'complete smaller-radius affine set')
            for j in candidates:
                g = gcd(j + i * z, q)
                require(g > 1, 'strict smaller-radius ansatz obstruction')
                smaller.append(dict(i=i, j=j, gcd=g))
        targets = [row['j'] for row in rows if max(fraction(row['slope']), fraction(row['endpoint'])) == C]
        marks[str(q)] = dict(q=q, radius=radius, rho=record(rho), A=record(A), C=record(C), positive_intervals=intervals, reflected_intervals=2 * intervals, special_classes=special, repair_table=encoded, periodic_rows=rows, C_maximizers=targets, obstruction=dict(z=z, tau=record(tau), radius=radius - 1, pairs=smaller, scope='unit/affine ansatz only'))
    require(marks['1155']['special_classes'] == [235, 418, 763, 952, 985, 1012], 'six exceptional1155 residues')
    require(marks['105']['C_maximizers'] == [47], 'periodic105 endpoint47')
    require(marks['1155']['C_maximizers'] == [634], 'periodic1155 slope634')
    atan_upper = sum(((-1) ** k * F(3, 4) ** (2 * k + 1) / F(2 * k + 1) for k in range(41)))
    atan_lower = atan_upper - F(3, 4) ** 83 / 83
    lo, hi = log_interval(F(5, 2))
    Klo = atan_lower - hi / 2
    Khi = atan_upper - lo / 2
    require(0 < Klo <= Khi < F(3, 16), 'integrated sharp remainder')
    log2lo, log2hi = log_interval(F(2))
    require(log2hi < F(25, 36), 'log2 rational upper')
    require(F(26, 15) ** 2 > 3 and F(49, 20) ** 2 > 6, 'radius3 square-root bounds')
    require(F(19, 6) ** 2 > 10 and F(31, 8) ** 2 > 15 and (F(9, 4) ** 2 > 5), 'radius5 square-root bounds')
    J3 = F(22, 7) * F(49, 40) + F(3, 2) * F(25, 36)
    J5 = F(22, 7) * F(19, 12) + F(31, 16) * F(25, 36)
    require(J3 == F(587, 120) and J5 == F(25489, 4032), 'exact sharp geometric factors')
    penalties = {}
    for q, R, J, basefactor in [(105, 3, J3, F(22, 7) * F(26, 15)), (1155, 5, J5, F(22, 7) * F(9, 4))]:
        A = fraction(marks[str(q)]['A'])
        C = fraction(marks[str(q)]['C'])
        penalties[str(q)] = dict(base_Phi_numerator=record(R * A), base_sqrtM_coefficient=record(basefactor * C), sharp_Phi_numerator=record(R * A / 2), sharp_linear_numerator=record(F(3, 16) * R * A), sharp_sqrtM_coefficient=record(J * C), sharp_geometry=record(J), M_free_Phi_numerator=record(F(3 * R, 2)), M_free_linear_numerator=record(F(9 * R, 16)), M_free_p_scalar_min=record(F(9 * R, 4)))
    require(fraction(penalties['105']['base_sqrtM_coefficient']) == F(177320, 6909), 'PE105 base coefficient')
    require(fraction(penalties['105']['sharp_sqrtM_coefficient']) == F(90985, 3948), 'PE105 sharp coefficient')
    banks = []
    for factors, R in [([3, 5, 11], 3), ([5, 7, 11], 3), ([3, 5, 7, 13], 5), ([3, 7, 11, 13], 5), ([5, 7, 11, 13], 5)]:
        q = 1
        for p in factors:
            q *= p
        table, _, n, special = verify_table(q, R)
        if 3 not in factors and len(factors) == 4:
            require(not special and all((row[2] <= 4 for rows in table.values() for row in rows)), 'r4 for four odd primes excluding3')
        banks.append(dict(factors=factors, q=q, radius=R, positive_intervals=n, special_classes=len(special), scope='finite regression of analytic universal lemma'))
    models = {105: [(121, 11), (1155, 11), (12705, 11), (169, 13), (1365, 13), (221, 13), (1785, 17)], 1155: [(169, 13), (15015, 13), (221, 17)]}
    for q, configs in models.items():
        row = marks[str(q)]
        A = fraction(row['A'])
        C = fraction(row['C'])
        result = [regression(M, q, p, row['radius'], tables[q], A, C) for M, p in configs]
        row['regressions'] = result
        row['regression_states'] = sum((x['states'] for x in result))
        row['repaired_bad'] = sum((x['bad'] for x in result))
    out = dict(status='PASS', checks=E.CHECKS, marks=marks, penalties=penalties, sharp_K_interval=dict(lower=record(Klo), upper=record(Khi), rational_cap=record(F(3, 16))), universal_odd_core_bank=banks, normalization='Delta(B,N)=mu(B^N)/ln(B)', M_free_density='min(1,3*X/p)', M_free_includes_zero=True, p_squared_in_M_allowed=True, profile_anchors_proved_here=False, finite_regressions_are_universal_proofs=False, code_reused_from_prior_deliveries=False)
    (HERE / 'theory105_certificate.json').write_text(json.dumps(out, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status='PASS', checks=E.CHECKS, marks={q: {k: v[k] for k in ['radius', 'positive_intervals', 'regression_states', 'repaired_bad']} for q, v in marks.items()}, penalties=penalties), indent=2))
if __name__ == '__main__':
    main()
