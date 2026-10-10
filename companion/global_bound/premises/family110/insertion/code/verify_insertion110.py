#!/usr/bin/env python3
"""Reconstruct the exact insertion110 certificate (standard library only).

Finite model checks are regressions, not a substitute for the universal
proof.  All decisive inequalities here use integers or Fraction.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import insertion110 as E

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def fs(x):
    return E.fraction_json(x)


def table110():
    table = {}
    grouped = {}
    for z in range(110):
        if gcd(z, 110) == 1:
            i, j = 1, 0
        else:
            candidates = [i for i in range(1, 5) if gcd(i * z - 1, 110) == 1]
            need(bool(candidates), "whole-half repair candidate exists")
            i, j = candidates[0], -1
        need(1 <= i <= 4 and gcd(j + i * z, 110) == 1, "positive table admissibility")
        need(max(abs(F(j)), abs(j + F(i, 2))) <= 1, "closed positive-half affine bound")
        need(gcd(-j - i * z, 110) == 1, "reflected residue admissibility")
        need(max(abs(F(-j)), abs(-j - F(i, 2))) <= 1, "closed reflected-half affine bound")
        table[z] = (i, j)
        grouped.setdefault((i, j), []).append(z)
    obstructions = []
    for z in (33, 105):
        obstruction = []
        for i in range(1, 4):
            for j in range(-4, 3):
                if abs(j + F(i, 4)) <= 1:
                    common = gcd(j + z * i, 110)
                    need(common > 1, "radius-three obstruction at tau1/4")
                    obstruction.append({"i": i, "j": j, "gcd": common})
        need(len(obstruction) == 6, "all possible smaller-radius pairs")
        obstructions.append({"z": z, "tau": fs(F(1, 4)), "pairs": obstruction})
    need(grouped[(4, -1)] == [33, 105], "reconstructed radius-four classes")
    return table, {
        "mark": 110, "radius": 4, "positive_half_intervals": 110,
        "closed_intervals_both_halves": 220,
        "grouped_table": [dict(i=i, j=j, z_values=zs) for (i, j), zs in sorted(grouped.items())],
        "restricted_radius_three_obstructions": obstructions,
        "radius_minimality_scope": "only the ansatz 1<=i<=R and abs(j+i*tau)<=1",
    }


def periodic110():
    units = [r for r in range(110) if gcd(r, 110) == 1]
    phi = len(units)
    need(phi == 40 and sum(units) == 110 * phi // 2, "unit count and symmetry sum")
    rho, A = F(phi, 110), F(5 * phi, 2 * 110)
    U = S = 0
    rows = []
    C = F(0)
    for j in range(110):
        if j and gcd(j, 110) == 1:
            U += 1
            S += j
        epsilon = U - rho * j
        D = 2 * j * U + S - A * j * j
        endpoint = D / j if j else F(0)
        C = max(C, 3 * epsilon, endpoint)
        rows.append((j, U, S, epsilon, D, endpoint))
    need((rho, A, C) == (E.RHO, E.PERIODIC_A, E.PERIODIC_C), "independently reconstructed coefficients")
    serialized = []
    maximizers = []
    for j, U, S, epsilon, D, endpoint in rows:
        need(3 * epsilon <= C and D <= C * j, "all-block slope and intercept inequalities")
        # Confirm both polynomial coefficients, then several exact block substitutions.
        need(F(5, 2) * 110 * phi - A * 110 * 110 == 0, "quadratic cancellation")
        need(2 * j * phi + 3 * 110 * U - 2 * A * 110 * j == 3 * 110 * epsilon,
             "linear coefficient identity")
        for n in (0, 1, 2, 100):
            T = 110 * n + j
            UT = phi * n + U
            ST = 110 * phi * n * n // 2 + 110 * n * U + S
            need(2 * T * UT + ST - A * T * T == 3 * 110 * epsilon * n + D,
                 "exact periodic polynomial regression")
        if endpoint == C:
            maximizers.append(j)
        serialized.append(dict(j=j, U=U, S=S, epsilon=fs(epsilon), D=fs(D),
                               slope=fs(3 * epsilon), endpoint=fs(endpoint)))
    need(maximizers == [73], "periodic endpoint maximum")
    j, Uj, Sj, _, _, _ = rows[73]
    need(2*j*Uj + Sj == A*j*j + C*j, "C is sharp for the fixed-A polynomial count")
    # A second counting identity avoids the residue-prefix implementation:
    # 1_(gcd(t,110)=1) = sum_(d|gcd(t,110)) mu(d).
    signed_divisors = [(1, 1), (2, -1), (5, -1), (10, 1),
                       (11, -1), (22, 1), (55, 1), (110, -1)]
    U = S = 0
    for T in range(1101):
        if T and gcd(T, 110) == 1:
            U += 1
            S += T
        ie_U = sum(sign * (T // d) for d, sign in signed_divisors)
        ie_S = sum(sign * d * (T // d) * (T // d + 1) // 2
                   for d, sign in signed_divisors)
        need((U, S) == (ie_U, ie_S), "independent inclusion-exclusion unit count and sum")
        need(2 * T * U + S <= A * T * T + C * T, "positive-integer periodic regression")
    for L in range(4097):
        need(L == 0 or F(5 * L * L + L, 2) <= 3 * L * L,
             "mark-free floor inequality, including L0 separately")
    return dict(rho=fs(rho), A=fs(A), C=fs(C), maximizer_j=maximizers,
                period_rows=serialized, regression_T_range=[0, 1100],
                independent_count="Mobius inclusion-exclusion over the eight squarefree divisors of 110",
                inclusion_exclusion_checks=1101, free_floor_L_range=[0, 4096],
                C_optimal_for_fixed_A_polynomial_count=True, C_optimality_witness_T=73,
                actual_bad_minimum_CDF_claimed_optimal=False)


def integral_constants():
    # Both positive finite integrals have a convergent alternating expansion.
    alo, ahi = F(866025, 10**6), F(866026, 10**6)
    blo, bhi = F(1154700, 10**6), F(1154701, 10**6)
    need(alo * alo < F(3, 4) < ahi * ahi < 1, "directed sqrt3/2 bounds")
    need(blo * blo < F(4, 3) < bhi * bhi < 2, "directed 2/sqrt3 bounds")

    def first(x, last):
        return sum((-1)**j * x**(4*j + 1) / (4*j + 1) for j in range(last + 1))

    def second(x, last):
        return sum((-1)**j * x**(2*j + 1) / (2**(j + 1) * (2*j + 1)) for j in range(last + 1))

    Jlo = 4 * (first(alo, 21) + second(blo, 41))
    Jhi = 4 * (first(ahi, 20) + second(bhi, 40))
    need(5 < Jlo < Jhi < E.J4_UPPER, "exact J4 below51/10")
    # Arctangent alternating series; an even last index gives the upper.
    atan_hi = sum((-1)**j * F(3, 4)**(2*j + 1) / (2*j + 1) for j in range(41))
    atan_lo = atan_hi - F(3, 4)**83 / 83
    loglo, loghi = E.log_bounds(F(5, 2))
    Klo, Khi = atan_lo - loghi / 2, atan_hi - loglo / 2
    need(0 < Klo < Khi < E.K_UPPER < F(3, 16), "exact linear integral constant")
    need(E.J4_UPPER < F(74731, 13230), "improved J coefficient below inherited bound")
    need(4 * E.PERIODIC_A / 2 == E.PERIODIC_PHI_NUMERATOR, "periodic Phi coefficient")
    need(4 * E.PERIODIC_A * E.K_UPPER == E.PERIODIC_LINEAR_NUMERATOR, "periodic linear coefficient")
    need(E.J4_UPPER * E.PERIODIC_C == E.PERIODIC_ROOT_COEFFICIENT, "periodic radical coefficient")
    need(F(12, 2) == E.FREE_PHI_NUMERATOR and 12 * E.K_UPPER == E.FREE_LINEAR_NUMERATOR,
         "free coefficients")
    need(F(74731, 13230) * E.PERIODIC_C == F(44315483, 1770615), "legacy root crosscheck")
    need(F(3, 4) * E.PERIODIC_A == F(15, 22), "legacy linear crosscheck")
    return {
        "J4": dict(lower=fs(Jlo), upper=fs(Jhi), simple_upper=fs(E.J4_UPPER),
                   upper_decimal=E.decimal_ceiling(Jhi),
                   a_interval=[fs(alo), fs(ahi)], b_interval=[fs(blo), fs(bhi)],
                   first_series_last_indices={"lower": 21, "upper": 20},
                   second_series_last_indices={"lower": 41, "upper": 40},
                   identity="4*int_0^(sqrt(3)/2)dt/(1+t^4)+4*int_0^(2/sqrt(3))dt/(t^2+2)"),
        "K": dict(lower=fs(Klo), upper=fs(Khi), simple_upper=fs(E.K_UPPER),
                  upper_decimal=E.decimal_ceiling(Khi), arctan_upper_last_index=40,
                  log_terms=80, identity="atan(3/4)-ln(5/2)/2"),
    }


def bezout(a, b):
    r0, r1, x0, x1, y0, y1 = a, b, 1, 0, 0, 1
    while r1:
        quotient = r0 // r1
        r0, r1 = r1, r0 - quotient * r1
        x0, x1 = x1, x0 - quotient * x1
        y0, y1 = y1, y0 - quotient * y1
    return (r0, x0, y0) if r0 > 0 else (-r0, -x0, -y0)


def oracle(M, q):
    """Definition-level enumeration with positive denominators; no geometry."""
    cap = max(1, M // 2)  # d=1 always supplies this cap.
    denominators = [d for d in range(1, cap + 1) if gcd(d, q) == 1]
    result = []
    for s in range(M):
        best = None
        for d in denominators:
            if best is not None and d > best[0]:
                break
            a = (s * d) % M
            if 2 * a > M:
                a -= M
            candidate = (max(abs(a), d), d, a)
            if best is None or candidate < best:
                best = candidate
        need(best is not None, "d1 supplies every marked minimum")
        H, d, a = best
        result.append((H, a, d))
    return result


def model_checks(table):
    models = [(13, 13), (26, 13), (39, 13), (65, 13), (143, 13), (169, 13),
              (286, 13), (338, 13), (1430, 13), (2860, 13), (4290, 13),
              (18590, 13), (17, 17), (289, 17), (1870, 17),
              (19, 19), (361, 19), (2090, 19)]
    records = []
    total_states = total_bad = total_zero = total_breakpoints = 0
    gap = None
    saw_p_square_with_full_mark = False
    for M, p in models:
        need(E.validate_regime(p, M), "actual model insertion regime")
        old, new = oracle(M, 110), oracle(M, 110 * p)
        bad = zero = 0
        heights = {}
        for s, (H, a, d) in enumerate(old):
            need(H <= new[s][0], "restriction cannot reduce height")
            if d % p:
                need(H == new[s][0], "chosen p-unit minimum incurs zero cost")
                continue
            bad += 1
            zero += (a == 0)
            heights[H] = heights.get(H, 0) + 1
            need(a % p == 0, "bad numerator divisible by p, zero allowed")
            k = (a - s * d) // M
            need(a == s * d + M * k and gcd(k, d) == 1, "lattice-coordinate primitivity")
            common, x, y = bezout(k, d)
            need(common == 1 and k * x + d * y == 1, "exact Bezout identity")
            b, e = s * x - M * y, x
            need(a * e - d * b == M and gcd(d, e) == 1 and e % p != 0,
                 "primitive lattice basis and nonzero new-prime residue")
            tau = F(a * b + d * e, a * a + d * d)
            shifted = tau + F(1, 2)
            n = shifted.numerator // shifted.denominator
            b -= n * a
            e -= n * d
            tau -= n
            need(abs(tau) <= F(1, 2), "centered shear")
            z = e * pow(d, -1, 110) % 110
            i, j = table[-z % 110] if tau < 0 else table[z]
            if tau < 0:
                j = -j
            aa, dd = j * a + i * b, j * d + i * e
            if dd < 0:
                aa, dd = -aa, -dd
            need(dd > 0 and gcd(dd, 110 * p) == 1 and (aa - s * dd) % M == 0,
                 "repaired marked denominator and lattice membership")
            Hr = max(abs(aa), dd)
            need(new[s][0] <= Hr and F(Hr) <= H + F(4 * M, H), "height repair and oracle comparison")
            eta = F(min(abs(a), d), H)
            v = F(4 * M, H * H)
            coordinate = max(1 + v * eta / (1 + eta * eta), eta + v / (1 + eta * eta))
            envelope_square = (1 + v / 2)**2 if v <= F(4, 3) else 1 + v * v
            need(F(Hr, H) <= coordinate and coordinate**2 <= envelope_square,
                 "structured coordinate geometry on actual repaired vector")
            if M % (110 * p) == 0:
                residual = M
                for prime in (2, 5, 11, p):
                    while residual % prime == 0:
                        residual //= prime
                if residual == 1:
                    need(gcd(dd, M) == 1, "full physical denominator at target radical")
        cumulative = 0
        for H in sorted(heights):
            cumulative += heights[H]
            need(cumulative <= min(F(M), E.PERIODIC_A * H * H / p + E.PERIODIC_C * H),
                 "periodic density at every bad-minimum breakpoint")
            need(cumulative <= min(F(M), F(3 * H * H, p)), "free density at every breakpoint")
        if M == 1430:
            need(old[110] == (13, 0, 13) and new[110][0] == 110, "marked/physical separation at zero numerator")
            need((-110) % M != 0, "ambient division leaves the lattice")
            gap = dict(M=M, s=110, old_minimum=list(old[110]), new_minimum=list(new[110]),
                       invalid_ambient_division=[0, 1])
        if M == 18590:
            saw_p_square_with_full_mark = M % (p * p) == 0 and zero > 0 and bad > 0
        total_states += M
        total_bad += bad
        total_zero += zero
        total_breakpoints += len(heights)
        records.append(dict(M=M, p=p, states=M, bad_minima=bad, zero_numerators=zero,
                            density_breakpoints=len(heights)))
    need(gap is not None and saw_p_square_with_full_mark, "zero and p-squared regimes explicitly exercised")
    return dict(total_states=total_states, repaired_bad_minima=total_bad,
                bad_zero_numerators=total_zero, density_breakpoints=total_breakpoints,
                models=records, marked_physical_gap=gap,
                claim_scope="finite regressions; universal validity comes from the written proof")


def geometry_grid():
    values = [F(k, 16) for k in range(1, 129)] + [F(4, 3), F(1, 1000), F(1000)]
    count = 0
    for eta in (F(k, 64) for k in range(65)):
        for v in values:
            first = 1 + v * eta / (1 + eta * eta)
            second = eta + v / (1 + eta * eta)
            envelope_square = (1 + v / 2)**2 if v <= F(4, 3) else 1 + v * v
            need(max(first, second)**2 <= envelope_square, "rational geometry grid")
            count += 1
    return dict(exact_pairs=count, eta_grid_denominator=64,
                scope="regression only; continuous proof is in REPAIR110_PROOF.md")


def consumer_checks():
    positives = [(lambda: E.validate_regime(13, 13)),
                 (lambda: E.validate_regime(13, 18590, True)),
                 (lambda: E.validate_regime(13, 4290, False)),
                 (lambda: E.periodic_upper(350, 1000000) > 0),
                 (lambda: E.normalized_upper(13, 1430, 1430, -100) == 0)]
    for fn in positives:
        need(bool(fn()), "positive consumer fixture")
    negatives = [
        ("nonprime theorem input", lambda: E.validate_regime(15, 1650)),
        ("new prime below domain", lambda: E.validate_regime(11, 1210)),
        ("new prime does not divide M", lambda: E.validate_regime(13, 1691)),
        ("extra prime in physical identification", lambda: E.validate_regime(13, 4290, True)),
        ("missing target primes", lambda: E.validate_regime(13, 169, True)),
        ("zero modulus", lambda: E.periodic_upper(13, 0)),
        ("boolean scalar prime", lambda: E.free_upper(True)),
        ("integration domain not paid", lambda: E.free_upper(5)),
        ("negative logarithm argument", lambda: E.log_bounds(-1)),
        ("unit normalization base", lambda: E.normalized_upper(13, 1, 1430, 0)),
        ("M below B", lambda: E.normalized_upper(13, 1430, 13, 0)),
    ]
    rejected = []
    for name, fn in negatives:
        try:
            fn()
        except ValueError:
            rejected.append(name)
        else:
            raise RuntimeError("negative fixture accepted: " + name)
        need(True, "negative consumer fixture rejected")
    # N=2 normalization must retain ln(B).  Clipping acts on the total numerator.
    p, B, M = 13, 1430, 1430**2
    cost = E.insertion_upper(p, M)
    anchor = -cost / 2
    got = E.normalized_upper(p, B, M, anchor)
    expected = cost / (2 * E.log_bounds(B)[0])
    need(got == expected and got > cost / (2 * E.log_bounds(M)[0]), "signed anchor and lnB normalization")
    need(got < E.insertion_upper(p, M) / E.log_bounds(B)[0], "anchor itself is not clipped")
    samples = []
    for p in (13, 17, 19, 47, 101, 350, 1000):
        previous = None
        for M in (110 * p, (110 * p)**2, 10**12):
            per, free = E.periodic_upper(p, M), E.free_upper(p)
            need(per < E.legacy_periodic_upper(p, M) and free < E.legacy_free_upper(p), "improved bounds beat legacy bounds")
            if previous is not None:
                need(per <= previous, "periodic scalar monotonicity in M")
            previous = per
            samples.append(dict(p=p, M=M, periodic_upper_decimal=E.decimal_ceiling(per),
                                free_upper_decimal=E.decimal_ceiling(free)))
    monotone_checks = 0
    for anchor in (F(-1), F(-1, 100), F(0), F(1, 10)):
        p0, B0, M0 = 13, 1430, 1430
        u0 = E.normalized_upper(p0, B0, M0, anchor)
        for p1, B1, M1 in ((17, 1430, 1430), (13, 2000, 2000), (13, 1430, 1430**2)):
            need(E.normalized_upper(p1, B1, M1, anchor) <= u0, "signed clipped scalar monotonicity")
            monotone_checks += 1
    return dict(positive_fixtures=len(positives), negative_fixtures=len(rejected), rejected=rejected,
                signed_anchor_and_lnB_checks=2, scalar_monotonicity_checks=monotone_checks,
                scalar_samples=samples, scalar_samples_are_not_Delta_values=True)


def build_certificate():
    table, repair = table110()
    periodic = periodic110()
    integrals = integral_constants()
    geometry = geometry_grid()
    models = model_checks(table)
    consumers = consumer_checks()
    return dict(status="PASS", checks=CHECKS, mark=110, new_prime_domain="p >= 13 prime, p | M",
                q_divides_M_required_for_insertion=False, p_squared_in_M_allowed=True,
                repair=repair, periodic=periodic, integral_constants=integrals,
                penalties=dict(
                    periodic_Phi_numerator=fs(E.PERIODIC_PHI_NUMERATOR),
                    periodic_linear_numerator=fs(E.PERIODIC_LINEAR_NUMERATOR),
                    periodic_sqrtM_coefficient=fs(E.PERIODIC_ROOT_COEFFICIENT),
                    free_Phi_numerator=fs(E.FREE_PHI_NUMERATOR),
                    free_linear_numerator=fs(E.FREE_LINEAR_NUMERATOR)),
                model_regressions=models, continuous_geometry_regressions=geometry,
                consumer_checks=consumers,
                normalization="max(0, anchor + insertion_upper) / log_lower(B)",
                profile_caps_proved_here=False, family110p_closed_here=False,
                universal_proof="REPAIR110_PROOF.md and INTEGRAL_CONSTANTS.md")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(result, indent=2) + "\n").encode()
    args.output.write_bytes(raw)
    print(json.dumps(dict(status="PASS", checks=CHECKS,
                         model_states=result["model_regressions"]["total_states"],
                         repaired_bad_minima=result["model_regressions"]["repaired_bad_minima"],
                         J4_upper=result["integral_constants"]["J4"]["upper_decimal"],
                         K_upper=result["integral_constants"]["K"]["upper_decimal"],
                         output_bytes=len(raw), output_sha256=sha256(raw).hexdigest()), indent=2))


if __name__ == "__main__":
    main()
