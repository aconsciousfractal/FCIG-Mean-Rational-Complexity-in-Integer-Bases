"""Exact insertion certificates for marks6,66,858 and uniform6r counting.

Infinite proofs are separate from finite model regressions. Arithmetic
helpers are transparently adapted from frozen70; this is not a second
implementation-independent audit. Only the Python standard library.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json
import exact_support as E
from exact_support import need, fs, log_bounds, periodic_certificate, make_table, oracle, bezout
ROOT = Path(__file__).resolve().parent
CONFIG = {6: dict(radius=2, rho=F(1, 3), A=F(5, 6), C=F(13, 6), intervals=6, obstruction=(3, F(1, 4)), example_p=5, models=[(5, 5), (10, 5), (25, 5), (30, 5), (50, 5), (7, 7), (42, 7), (49, 7), (66, 11), (121, 11)]), 66: dict(radius=4, rho=F(10, 33), A=F(25, 33), C=F(5645, 1617), intervals=66, obstruction=(11, F(1, 4)), example_p=13, models=[(13, 13), (26, 13), (65, 13), (91, 13), (130, 13), (169, 13), (221, 13), (286, 13), (338, 13), (858, 13), (17, 17), (51, 17), (221, 17), (289, 17), (374, 17), (1122, 17), (19, 19), (247, 19), (361, 19), (1254, 19), (23, 23), (529, 23)]), 858: dict(radius=6, rho=F(40, 143), A=F(100, 143), C=F(633619, 120263), intervals=866, obstruction=(267, F(1, 10)), example_p=17, models=[(6006, 7), (14586, 17), (289, 17)])}

def certify_mark(q, cfg):
    radius = cfg['radius']
    table = make_table(q, radius)
    serialized = []
    count = 0
    special = []
    for z, rows in table.items():
        reach = F(0)
        out = []
        for lo, hi, i, j in rows:
            need(0 <= lo <= reach < hi <= F(1, 2), 'continuous positive coverage')
            need(1 <= i <= radius and gcd(j + i * z, q) == 1, 'radius and old mark')
            need(max(abs(j + i * lo), abs(j + i * hi)) <= 1, 'positive affine endpoints')
            need(gcd(-j - i * z, q) == 1 and max(abs(-j - i * lo), abs(-j - i * hi)) <= 1, 'reflected affine and unit conditions')
            reach = hi
            out.append(dict(lo=fs(lo), hi=fs(hi), i=i, j=j))
        need(reach == F(1, 2), 'closed half endpoint')
        count += len(rows)
        if len(rows) > 1:
            special.append(z)
        serialized.append(dict(z=z, intervals=out))
    need(count == cfg['intervals'], 'interval count')
    z, tau = cfg['obstruction']
    obstruct = []
    for i in range(1, radius):
        for j in range(-radius, 3):
            if abs(j + i * tau) <= 1:
                g = gcd(j + i * z, q)
                need(g > 1, 'restricted smaller-radius obstruction')
                obstruct.append(dict(i=i, j=j, gcd=g))
    need(len(obstruct) == 2 * (radius - 1), 'all possible obstruction pairs')
    rho, A, C, period = periodic_certificate(q)
    need((rho, A, C) == (cfg['rho'], cfg['A'], cfg['C']), 'reconstructed exact periodic coefficients')
    U = S = 0
    for T in range(10 * q + 1):
        if T and gcd(T, q) == 1:
            U += 1
            S += T
        need(2 * T * U + S <= A * T * T + C * T, 'finite periodic regression')
    model_rows = []
    states = bad_total = breakpoints = 0
    gap = None
    for M, p in cfg['models']:
        need(p > radius and q % p and (M % p == 0), 'model insertion hypotheses')
        marked, full = (oracle(M, q), oracle(M, q * p))
        if M == q * cfg['example_p'] and p == cfg['example_p']:
            need(marked[q] == (p, 0, p) and full[q][0] == q, 'strict marked/full model gap')
            need(-q % M != 0, 'ambient gcd division fails lattice membership')
            gap = dict(M=M, s=q, marked=marked[q], full=full[q])
        bad = 0
        heights = {}
        for s, (H, a, d) in enumerate(marked):
            need(H <= full[s][0], 'restriction direction')
            if d % p:
                need(H == full[s][0], 'no cost for chosen p-unit minimum')
                continue
            bad += 1
            heights[H] = heights.get(H, 0) + 1
            need(a % p == 0, 'bad numerator includes zero')
            k = (a - s * d) // M
            g, x, y = bezout(k, d)
            need(g == 1, 'lattice-coordinate primitivity')
            b, e = (s * x - M * y, x)
            need(a * e - d * b == M and gcd(d, e) == 1 and e % p, 'basis and second denominator')
            norm2 = a * a + d * d
            tau = F(a * b + d * e, norm2)
            shifted = tau + F(1, 2)
            shift = shifted.numerator // shifted.denominator
            b -= shift * a
            e -= shift * d
            tau -= shift
            need(abs(tau) <= F(1, 2), 'reduced shear')
            z = e * pow(d, -1, q) % q
            zz, tt = (-z % q, -tau) if tau < 0 else (z, tau)
            choices = [v for v in table[zz] if v[0] <= tt <= v[1]]
            need(bool(choices), 'closed interval branch present')
            _, _, i, j = choices[0]
            if tau < 0:
                j = -j
            aa, dd = (j * a + i * b, j * d + i * e)
            if dd < 0:
                aa, dd = (-aa, -dd)
            Hr = max(abs(aa), dd)
            need(dd > 0 and gcd(dd, q * p) == 1 and ((aa - s * dd) % M == 0), 'repaired full admissibility')
            need(F(Hr) <= H + F(radius * M, H) and full[s][0] <= Hr, 'repaired height and full minimum')
            eta = F(min(abs(a), d), H)
            v = F(radius * M, H * H)
            coordinate = max(1 + v * eta / (1 + eta * eta), eta + v / (1 + eta * eta))
            need(F(Hr, H) <= coordinate, 'coordinate-sensitive estimate on actual witness')
            envelope = (1 + v / 2) ** 2 if v <= F(4, 3) else 1 + v * v
            need(coordinate ** 2 <= envelope, 'sharp geometric envelope on actual model')
        cumulative = 0
        for H in sorted(heights):
            cumulative += heights[H]
            need(cumulative <= min(F(M), A * H * H / p + C * H), 'periodic density at every breakpoint')
            need(cumulative <= min(F(M), F(3 * H * H, p)), 'mark-independent M-free density')
            breakpoints += 1
        states += M
        bad_total += bad
        model_rows.append(dict(M=M, p=p, states=M, bad=bad, density_breakpoints=len(heights)))
    need(gap is not None, 'explicit model-separation example included')
    return dict(q=q, radius=radius, rho=fs(rho), A=fs(A), C=fs(C), positive_half_intervals=count, reflected_closed_intervals=2 * count, special_classes=special, repair_table=serialized, periodic_rows=period, restricted_obstruction=dict(z=cfg['obstruction'][0], tau=fs(cfg['obstruction'][1]), smaller_radius=radius - 1, pairs=obstruct, scope='only abs(j+i*tau)<=1 ansatz'), periodic_regression_range=[0, 10 * q], model_instances=states, bad_minima_repaired=bad_total, density_breakpoints=breakpoints, model_records=model_rows, model_gap=gap)

def main():
    marks = {str(q): certify_mark(q, cfg) for q, cfg in CONFIG.items()}
    uniform = []
    for r in [5, 7, 11, 13, 17, 19]:
        q = 6 * r
        count = 0
        for n in range(24 * r + 1):
            if n and n % 2 == 0 and (gcd(n // 2, q) == 1):
                count += 1
            t = F(n, 2)

            def U6(x):
                return (x + 5) // 6 + (x + 1) // 6
            first = U6(t)
            second = U6(t / r)
            need(abs(first - t / 3) <= F(2, 3), 'U6 real-count error')
            need(count == first - second, 'exact U6r inclusion/exclusion')
            rho = F(r - 1, 3 * r)
            need(abs(count - rho * t) <= F(4, 3), 'uniform U6r real-count error')
        rho, A, C, _ = periodic_certificate(q)
        need(A <= F(5, 6) and C <= F(16, 3), 'uniform periodic bank')
        uniform.append(dict(r=r, q=q, A=fs(A), C=fs(C), half_grid_max=12 * r))
    for s in [3, 5, 7, 11, 13]:
        for z, rows in make_table(2 * s, 2).items():
            need(len(rows) == 1, 'two-prime analytic radius-two bank')
    for A in range(1001):
        need(A == 0 or F(5 * A * A + A, 2) <= 3 * A * A, 'M-free integer floor inequality')
    atan_hi = sum(((-1) ** j * F(3, 4) ** (2 * j + 1) / F(2 * j + 1) for j in range(21)))
    atan_lo = atan_hi - F(3, 4) ** 43 / 43
    l, h = log_bounds(F(5, 2))
    Klo = atan_lo - h / 2
    Khi = atan_hi - l / 2
    need(0 < Klo <= Khi < F(3, 16), 'sharp integrated linear remainder')
    lo2, hi2 = log_bounds(2)
    need(F(2, 3) < lo2 <= hi2 < F(25, 36), 'both rational log2 bounds')
    need(F(99, 70) ** 2 > 2 and F(26, 15) ** 2 > 3 and (F(49, 20) ** 2 > 6), 'rational square roots')
    geometry = {2: F(22, 7) + F(49, 40) * F(25, 36), 4: F(22, 7) * F(99, 70) + F(26, 15) * F(25, 36), 6: F(22, 7) * F(26, 15) + F(3, 2) * F(99, 70) * F(25, 36)}
    need(geometry == {2: F(8051, 2016), 4: F(74731, 13230), 6: F(11627, 1680)}, 'sharp geometric coefficients')
    basegeometry = {2: F(22, 7) * F(99, 70), 4: F(44, 7), 6: F(77, 10)}
    penalties = {}
    for q, cfg in CONFIG.items():
        r, A, C = (cfg['radius'], cfg['A'], cfg['C'])
        penalties[str(q)] = dict(base_Phi_numerator=fs(r * A), base_sqrtM_coefficient=fs(basegeometry[r] * C), sharp_Phi_numerator=fs(r * A / 2), sharp_linear_numerator=fs(F(3, 16) * r * A), sharp_sqrtM_coefficient=fs(geometry[r] * C), M_free_Phi_numerator=fs(F(3 * r, 2)), M_free_linear_numerator=fs(F(9 * r, 16)), M_free_condition_p_ge=fs(F(9 * r, 4)))
    need(F(44, 7) * CONFIG[66]['C'] == F(248380, 11319) == F(22580, 1029), 'base66 coefficient identity')
    need(geometry[4] * CONFIG[66]['C'] == F(84371299, 4278582), 'sharp66 coefficient identity')
    need(F(3, 16) * 4 * CONFIG[66]['A'] == F(25, 44), 'sharp66 linear identity')
    need(geometry[6] * CONFIG[858]['C'] == F(95676469, 2623920), 'sharp858 coefficient identity')
    need(basegeometry[6] * CONFIG[858]['C'] == F(4435333, 109330), 'base858 coefficient identity')
    need(F(3, 16) * 6 * CONFIG[858]['A'] == F(225, 286), 'sharp858 linear identity')
    need(geometry[2] * CONFIG[6]['C'] == F(104663, 12096), 'sharp6 coefficient identity')
    need(geometry[4] * F(16, 3) == F(597848, 19845), 'uniform6r coefficient identity')
    output = dict(status='PASS', checks=E.CHECKS, marks=marks, penalties=penalties, sharp_linear_remainder=dict(lower=fs(Klo), upper=fs(Khi), rational_cap=fs(F(3, 16))), sharp_geometry={str(k): fs(v) for k, v in geometry.items()}, M_free=dict(mark_independent=True, density='min(1,3*X/p)', base='Phi(3*radius/p)', sharp='Phi(3*radius/(2*p))+9*radius/(16*p)', condition='p>=9*radius/4', integer_regression_range=[0, 1000], p_squared_in_M_allowed=True, numerator_zero_included=True), uniform6r=dict(r_prime_min=5, rho='(r-1)/(3*r)', A_upper=fs(F(5, 6)), C_upper=fs(F(16, 3)), sharp_Phi_numerator=fs(F(5, 3)), sharp_linear_numerator=fs(F(5, 8)), sharp_sqrtM_coefficient=fs(F(597848, 19845)), finite_count_bank=uniform), claim_scope='exact arithmetic and finite regressions; universal arguments in written proofs', profile_caps_proved_here=False)
    (ROOT / 'theory66_certificate.json').write_text(json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status='PASS', checks=E.CHECKS, marks={q: {k: v[k] for k in ['radius', 'reflected_closed_intervals', 'model_instances', 'bad_minima_repaired']} for q, v in marks.items()}), indent=2))
if __name__ == '__main__':
    main()
