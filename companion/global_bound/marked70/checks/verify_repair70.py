#!/usr/bin/env python3
"""Exact certificates for R70/P70/D70/PE70 and finite model regressions.

The infinite conclusions are proved in REPAIR70_PROOF.md.  Finite model
regressions are explicitly separate from the proof of a universal assertion.
The arithmetic helpers adapt the archived verify_repair42.py; this is not
an implementation-independent audit.  Only the Python standard library.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
Q = 70
RHO = F(12, 35)
A70 = F(6, 7)
C70 = F(1170, 301)
PE70 = F(51480, 2107)
CHECKS = 0


def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def fs(x):
    x = F(x)
    return {"numerator": x.numerator, "denominator": x.denominator}


def log_bounds(x, terms=80):
    """Rational enclosure of ln(x), binary reduction and an explicit tail."""
    x = F(x)
    need(x > 0, "positive logarithm")
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def atanh_bounds(z):
        z = F(z)
        zz = z*z
        power = z
        total = F(0)
        for j in range(terms):
            total += 2*power/(2*j+1)
            power *= zz
        tail = 2*power/((2*terms+1)*(1-zz))
        return total, total+tail

    lo, hi = atanh_bounds((x-1)/(x+1))
    l2, h2 = atanh_bounds(F(1, 3))
    return (lo+k*l2, hi+k*h2) if k >= 0 else (lo+k*h2, hi+k*l2)


def phi_bounds(x):
    x = F(x)
    need(x > 0, "positive Phi argument")
    a, b = log_bounds(1+x)
    c, d = log_bounds(1+1/x)
    return a+x*c, b+x*d


def repair_pair(z, q):
    """Analytic radius-four rule for q=2*r*s, distinct odd primes r,s."""
    z %= q
    if z % 2 == 0:
        for i in (1, 2, 3):
            if gcd(i*z-1, q) == 1:
                return i, -1
    elif gcd(z, q) == 1:
        return 1, 0
    else:
        for i in (2, 4):
            if gcd(i*z-1, q) == 1:
                return i, -1
    raise RuntimeError("radius-four rule used outside its proved core class")


def table_for(q):
    table = {}
    for z in range(q):
        i, j = repair_pair(z, q)
        need(1 <= i <= 4, "repair radius")
        need(gcd(i*z+j, q) == 1, "old-mark unit")
        need(max(abs(F(j)), abs(F(j)+F(i, 2))) <= 1,
             "entire closed positive half-interval")
        # Reflection is algebraic; check its modular and endpoint conditions.
        ni, nj = repair_pair(-z, q)
        nj = -nj
        need(gcd(ni*z+nj, q) == 1, "reflected old-mark unit")
        need(max(abs(F(nj)-F(ni, 2)), abs(F(nj))) <= 1,
             "entire closed negative half-interval")
        table[z] = (i, j)
    return table


def periodic_certificate(q):
    phi = sum(gcd(j, q) == 1 for j in range(q))
    rho = F(phi, q)
    quadratic = F(5, 2)*rho
    need(sum(j for j in range(q) if gcd(j, q) == 1) == q*phi//2,
         "unit residue sum")
    U = S = 0
    rows = []
    C = F(0)
    for j in range(q):
        if j and gcd(j, q) == 1:
            U += 1
            S += j
        eps = U-rho*j
        D = 2*j*U+S-quadratic*j*j
        C = max(C, 3*eps, D/j if j else F(0))
        rows.append((j, U, S, eps, D))
    output = []
    for j, U, S, eps, D in rows:
        need(3*eps <= C, "periodic slope coefficient")
        need(D <= C*j, "periodic endpoint coefficient")
        output.append({"j": j, "U": U, "S": S, "epsilon": fs(eps),
                       "D": fs(D), "slope": fs(3*eps),
                       "endpoint": fs(D/j) if j else None})
    return rho, quadratic, C, output


def oracle(M, q):
    """Definition itself: d=1 bounds the search by max(1,floor(M/2))."""
    cap = max(1, M//2)
    ds = [d for d in range(1, cap+1) if gcd(d, q) == 1]
    out = []
    for s in range(M):
        best = None
        for d in ds:
            if best is not None and d > best[0]:
                break
            a = (s*d) % M
            if 2*a > M:
                a -= M
            candidate = max(abs(a), d), d, a
            if best is None or candidate < best:
                best = candidate
        need(best is not None, "denominator one exists")
        H, d, a = best
        out.append((H, a, d))
    return out


def bezout(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        c = old_r//r
        old_r, r = r, old_r-c*r
        old_s, s = s, old_s-c*s
        old_t, t = t, old_t-c*t
    if old_r < 0:
        return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def main():
    table = table_for(Q)
    groups = []
    for i, j in ((1, 0), (1, -1), (2, -1), (3, -1), (4, -1)):
        groups.append({"i": i, "j": j,
                       "classes": [z for z in range(Q) if table[z] == (i, j)]})
    obstructions = []
    for z in (25, 63):
        tau = F(1, 4)
        pairs = []
        for i in range(1, 4):
            for j in range(-4, 5):
                if abs(j+i*tau) <= 1:
                    g = gcd(j+i*z, Q)
                    need(g > 1, "restricted radius-three obstruction")
                    pairs.append({"i": i, "j": j, "gcd": g})
        obstructions.append({"z": z, "tau": fs(tau), "pairs": pairs})

    rho, quadratic, C, period = periodic_certificate(Q)
    need(rho == RHO and quadratic == A70 and C == C70, "P70 exact constants")
    need(C70 == F(period[43]["D"]["numerator"],
                  period[43]["D"]["denominator"])/43,
         "C70 attained at j=43")
    need(PE70 == F(44, 7)*C70, "rational pi coefficient")
    U = S = 0
    for T in range(7001):
        if T and gcd(T, Q) == 1:
            U += 1
            S += T
        need(2*T*U+S <= A70*T*T+C70*T, "direct P70 regression")

    # This bank corroborates the analytic theorem, not its quantifier over r,s.
    core_bank = []
    for r, s in ((3, 5), (3, 7), (3, 11), (5, 7),
                 (3, 13), (5, 11), (7, 11), (11, 13)):
        q = 2*r*s
        table_for(q)
        rr, aa, cc, _ = periodic_certificate(q)
        need(rr == F((r-1)*(s-1), 2*r*s), "three-prime density")
        need(cc <= 8, "uniform periodic coefficient regression")
        inv = pow(2*r, -1, s)
        odd_multiplier = inv if inv % 2 else inv+s
        z = r*odd_multiplier % q
        need(z % 2 == 1 and z % r == 0 and (2*z-1) % s == 0,
             "CRT obstruction construction")
        for i in (1, 2, 3):
            for j in (-1, 0):
                need(gcd(i*z+j, q) > 1, "generic restricted radius-three failure")
        count = 0
        for n in range(4*q+1):
            if n and n % 2 == 0 and gcd(n//2, q) == 1:
                count += 1
            def odd_count(d):
                return (n+2*d)//(4*d)
            inclusion = odd_count(1)-odd_count(r)-odd_count(s)+odd_count(r*s)
            need(count == inclusion, "odd-multiple inclusion-exclusion")
            need(abs(F(count)-rr*F(n, 2)) <= 2, "uniform real-count error bank")
        core_bank.append({"r": r, "s": s, "q": q, "rho": fs(rr),
                          "A": fs(aa), "C_periodic": fs(cc),
                          "radius_three_obstruction_z": z,
                          "half_grid_range": [0, 2*q]})

    model_pairs = [
        (11, 11), (22, 11), (55, 11), (77, 11), (110, 11),
        (121, 11), (143, 11), (154, 11), (242, 11), (770, 11),
        (13, 13), (39, 13), (143, 13), (169, 13), (286, 13),
        (910, 13), (17, 17), (187, 17), (289, 17), (1190, 17),
        (19, 19), (361, 19),
    ]
    states = bad = density_checks = 0
    model_records = []
    for M, p in model_pairs:
        need(M % p == 0 and p > 4 and Q % p != 0, "PE model hypotheses")
        marked, full = oracle(M, Q), oracle(M, Q*p)
        bad_M = 0
        bad_heights = {}
        for s, (H, a, d) in enumerate(marked):
            Hf = full[s][0]
            need(H <= Hf, "restriction direction")
            if d % p:
                need(H == Hf, "chosen p-unit minimizer needs no repair")
                continue
            bad += 1
            bad_M += 1
            bad_heights[H] = bad_heights.get(H, 0)+1
            need(a % p == 0, "bad numerator also divisible by p")
            k = (a-s*d)//M
            g, x, y = bezout(k, d)
            need(g == 1, "lattice-coordinate primitive minimum")
            e, b = x, s*x-M*y
            need(a*e-d*b == M, "lattice basis determinant")
            need(gcd(d, e) == 1 and e % p != 0, "second denominator avoids p")
            norm2 = a*a+d*d
            tau = F(a*b+d*e, norm2)
            shifted = tau+F(1, 2)
            shift = shifted.numerator//shifted.denominator
            b -= shift*a
            e -= shift*d
            tau -= shift
            need(abs(tau) <= F(1, 2), "reduced lattice shear")
            z = e*pow(d, -1, Q) % Q
            if tau >= 0:
                i, j = table[z]
            else:
                i, j = table[-z % Q]
                j = -j
            aa, dd = j*a+i*b, j*d+i*e
            if dd < 0:
                aa, dd = -aa, -dd
            need(dd > 0 and gcd(dd, Q*p) == 1, "repaired admissible denominator")
            need((aa-s*dd) % M == 0, "repaired lattice membership")
            Hr = max(abs(aa), dd)
            need(F(Hr) <= H+F(4*M, H), "R70 height estimate")
            need(Hf <= Hr, "repaired witness bounds the full minimum")
        cumulative = 0
        for H in sorted(bad_heights):
            cumulative += bad_heights[H]
            need(cumulative <= min(F(M), A70*H*H/p+C70*H),
                 "D70 density at every bad-height breakpoint")
            density_checks += 1
        states += M
        model_records.append({"M": M, "p": p, "states": M,
                              "bad": bad_M, "distinct_bad_heights": len(bad_heights)})

    marked = oracle(770, 70)[70]
    full = oracle(770, 770)[70]
    need(marked == (11, 0, 11) and full[0] == 70, "strict marked/full gap")
    need((-70) % 770 != 0, "ambient gcd division does not preserve lattice")
    output = {
        "status": "PASS", "q": Q, "checks": CHECKS,
        "repair_radius": 4, "classes": 70, "closed_half_intervals": 140,
        "repair_groups": groups,
        "repair_table": [{"z": z, "i": table[z][0], "j": table[z][1]}
                         for z in range(Q)],
        "restricted_radius_three_obstructions": obstructions,
        "rho": fs(RHO), "A70": fs(A70), "C70": fs(C70),
        "periodic_rows": period,
        "PE_argument": "24/(7p)", "PE_sqrtM_coefficient": fs(PE70),
        "uniform_2rs": {
            "repair_radius": 4, "unit_count_error": 2,
            "periodic_linear_coefficient": 8,
            "PE_argument": "10*rho/p < 5/p",
            "PE_sqrtM_coefficient": fs(F(352, 7)),
            "core_regression_bank": core_bank,
            "scope": "infinite theorem established by proof, bank is finite regression only"
        },
        "direct_periodic_integer_range": [0, 7000],
        "model_pairs": model_records, "model_instances": states,
        "bad_minima_repaired": bad, "density_breakpoints_checked": density_checks,
        "model_gap": {"M": 770, "s": 70, "marked": marked, "full": full},
        "claim_scope": "exact arithmetic certificates and finite model regressions; analytic proofs separate",
        "profile_caps": "separate dependency, not proved or imported by this script",
    }
    path = ROOT/"repair70_certificate.json"
    path.write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "PASS", "certificate": path.name,
                      "checks": CHECKS, "model_instances": states,
                      "bad_minima_repaired": bad,
                      "density_breakpoints_checked": density_checks}, indent=2))


if __name__ == "__main__":
    main()
