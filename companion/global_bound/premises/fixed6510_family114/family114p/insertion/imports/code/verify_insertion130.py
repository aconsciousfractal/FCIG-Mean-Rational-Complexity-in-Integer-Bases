#!/usr/bin/env python3
"""Reconstruct the exact insertion130 certificate (standard library only).

Finite model checks are regressions, not a substitute for the universal
proof.  All decisive inequalities here use integers or Fraction.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd
from pathlib import Path
import insertion130 as E

CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def fs(x):
    return E.fraction_json(x)


def table130():
    table = {}
    grouped = {}
    for z in range(130):
        if gcd(z, 130) == 1:
            i, j = 1, 0
        else:
            candidates = [i for i in range(1, 5) if gcd(i * z - 1, 130) == 1]
            need(bool(candidates), "whole-half repair candidate exists")
            i, j = candidates[0], -1
        need(1 <= i <= 4 and gcd(j + i * z, 130) == 1, "positive table admissibility")
        need(max(abs(F(j)), abs(j + F(i, 2))) <= 1, "closed positive-half affine bound")
        need(gcd(-j - i * z, 130) == 1, "reflected residue admissibility")
        need(max(abs(F(-j)), abs(-j - F(i, 2))) <= 1, "closed reflected-half affine bound")
        table[z] = (i, j)
        grouped.setdefault((i, j), []).append(z)
    obstructions = []
    for z in (13, 85):
        obstruction = []
        for i in range(1, 4):
            for j in range(-4, 3):
                if abs(j + F(i, 4)) <= 1:
                    common = gcd(j + z * i, 130)
                    need(common > 1, "radius-three obstruction at tau1/4")
                    obstruction.append({"i": i, "j": j, "gcd": common})
        need(len(obstruction) == 6, "all possible smaller-radius pairs")
        obstructions.append({"z": z, "tau": fs(F(1, 4)), "pairs": obstruction})
    need(grouped[(4, -1)] == [13, 85], "reconstructed radius-four classes")
    return table, {
        "mark": 130, "radius": 4, "positive_half_intervals": 130,
        "closed_intervals_both_halves": 260,
        "grouped_table": [dict(i=i, j=j, z_values=zs) for (i, j), zs in sorted(grouped.items())],
        "restricted_radius_three_obstructions": obstructions,
        "radius_minimality_scope": "only the ansatz 1<=i<=R and abs(j+i*tau)<=1",
    }


def periodic130():
    units = [r for r in range(130) if gcd(r, 130) == 1]
    phi = len(units)
    need(phi == 48 and sum(units) == 130 * phi // 2, "unit count and symmetry sum")
    rho, A = F(phi, 130), F(5 * phi, 2 * 130)
    U = S = 0
    rows = []
    C = F(0)
    for j in range(130):
        if j and gcd(j, 130) == 1:
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
        need(F(5, 2) * 130 * phi - A * 130 * 130 == 0, "quadratic cancellation")
        need(2 * j * phi + 3 * 130 * U - 2 * A * 130 * j == 3 * 130 * epsilon,
             "linear coefficient identity")
        for n in (0, 1, 2, 100):
            T = 130 * n + j
            UT = phi * n + U
            ST = 130 * phi * n * n // 2 + 130 * n * U + S
            need(2 * T * UT + ST - A * T * T == 3 * 130 * epsilon * n + D,
                 "exact periodic polynomial regression")
        if endpoint == C:
            maximizers.append(j)
        serialized.append(dict(j=j, U=U, S=S, epsilon=fs(epsilon), D=fs(D),
                               slope=fs(3 * epsilon), endpoint=fs(endpoint)))
    need(maximizers == [83], "periodic endpoint maximum")
    j, Uj, Sj, _, _, _ = rows[83]
    need(2*j*Uj + Sj == A*j*j + C*j, "C is sharp for the fixed-A polynomial count")
    # A second counting identity avoids the residue-prefix implementation:
    # 1_(gcd(t,130)=1) = sum_(d|gcd(t,130)) mu(d).
    signed_divisors = [(1, 1), (2, -1), (5, -1), (10, 1),
                       (13, -1), (26, 1), (65, 1), (130, -1)]
    U = S = 0
    for T in range(1301):
        if T and gcd(T, 130) == 1:
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
                period_rows=serialized, regression_T_range=[0, 1300],
                independent_count="Mobius inclusion-exclusion over the eight squarefree divisors of 130",
                inclusion_exclusion_checks=1301, free_floor_L_range=[0, 4096],
                C_optimal_for_fixed_A_polynomial_count=True, C_optimality_witness_T=83,
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
    models = [(17,17),(34,17),(85,17),(221,17),(51,17),(289,17),
              (578,17),(2210,17),(4420,17),(6630,17),(37570,17),
              (19,19),(361,19),(2470,19)]
    records = []
    total_states = total_bad = total_zero = total_breakpoints = 0
    gap = None
    saw_p_square_with_full_mark = False
    for M, p in models:
        need(E.validate_regime(p, M), "actual model insertion regime")
        old, new = oracle(M, 130), oracle(M, 130 * p)
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
            z = e * pow(d, -1, 130) % 130
            i, j = table[-z % 130] if tau < 0 else table[z]
            if tau < 0:
                j = -j
            aa, dd = j * a + i * b, j * d + i * e
            if dd < 0:
                aa, dd = -aa, -dd
            need(dd > 0 and gcd(dd, 130 * p) == 1 and (aa - s * dd) % M == 0,
                 "repaired marked denominator and lattice membership")
            Hr = max(abs(aa), dd)
            need(new[s][0] <= Hr and F(Hr) <= H + F(4 * M, H), "height repair and oracle comparison")
            eta = F(min(abs(a), d), H)
            v = F(4 * M, H * H)
            coordinate = max(1 + v * eta / (1 + eta * eta), eta + v / (1 + eta * eta))
            envelope_square = (1 + v / 2)**2 if v <= F(4, 3) else 1 + v * v
            need(F(Hr, H) <= coordinate and coordinate**2 <= envelope_square,
                 "structured coordinate geometry on actual repaired vector")
            if M % (130 * p) == 0:
                residual = M
                for prime in (2, 5, 13, p):
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
        if M == 2210:
            need(old[130] == (17, 0, 17) and new[130][0] == 130, "marked/physical separation at zero numerator")
            need((-130) % M != 0, "ambient division leaves the lattice")
            gap = dict(M=M, s=130, old_minimum=list(old[130]), new_minimum=list(new[130]),
                       invalid_ambient_division=[0, 1])
        if M == 37570:
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
                scope="regression only; continuous proof is in REPAIR130_PROOF.md")


def consumer_checks():
    positives = [lambda:E.validate_regime(17,17), lambda:E.validate_regime(17,37570,True),
                 lambda:E.validate_regime(17,6630,False), lambda:E.normalized_upper(17,2210,2210,-100)==0]
    for fn in positives: need(bool(fn()), "positive consumer fixture")
    negatives = [
      ("nonprime theorem input", lambda:E.validate_regime(25,3250)),
      ("new prime below domain", lambda:E.validate_regime(13,1690)),
      ("new prime does not divide M", lambda:E.validate_regime(17,2211)),
      ("extra prime in physical identification", lambda:E.validate_regime(17,6630,True)),
      ("missing target primes", lambda:E.validate_regime(17,289,True)),
      ("zero modulus", lambda:E.periodic_upper(17,0)),
      ("boolean scalar prime", lambda:E.free_upper(True)),
      ("negative logarithm", lambda:E.log_bounds(-1)),
      ("unit base", lambda:E.normalized_upper(17,1,2210,0)),
      ("M below B", lambda:E.normalized_upper(17,2210,17,0))]
    rejected=[]
    for name,fn in negatives:
      try:fn()
      except ValueError:rejected.append(name)
      else:raise RuntimeError("negative fixture accepted: "+name)
      need(True,"negative fixture rejected")
    p,B,M=17,2210,2210**2;cost=E.insertion_upper(p,M);anchor=-cost/2
    got=E.normalized_upper(p,B,M,anchor)
    need(got==cost/(2*E.log_bounds(B)[0]) and got>cost/(2*E.log_bounds(M)[0]),"signed anchor and lnB")
    need(got<E.insertion_upper(p,M)/E.log_bounds(B)[0],"anchor itself not clipped")
    checks=0
    for anchor in (F(-1),F(-1,100),F(0),F(1,10)):
      initial=E.normalized_upper(17,2210,2210,anchor)
      for p,B,M in [(19,2210,2210),(17,2500,2500),(17,2210,2210**2)]:
       need(E.normalized_upper(p,B,M,anchor)<=initial,"signed clipped scalar monotonicity");checks+=1
    return dict(positive_fixtures=len(positives),negative_fixtures=len(rejected),rejected=rejected,
                signed_anchor_and_lnB_checks=2,scalar_monotonicity_checks=checks)


def build_certificate():
    table, repair = table130()
    periodic = periodic130()
    integrals = integral_constants()
    geometry = geometry_grid()
    models = model_checks(table)
    consumers = consumer_checks()
    return dict(status="PASS", checks=CHECKS, mark=130, new_prime_domain="p >= 17 prime, p | M",
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
                profile_caps_proved_here=False, family130p_closed_here=False,
                universal_proof="REPAIR130_PROOF.md and INTEGRAL_CONSTANTS.md")


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
