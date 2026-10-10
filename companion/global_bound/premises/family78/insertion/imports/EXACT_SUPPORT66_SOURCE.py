#!/usr/bin/env python3
"""Exact shared helpers, adapted transparently from the frozen70 lane.

This is not an independent arithmetic implementation. New scientific
constants and tables are reconstructed by verify_theory66.py.
"""
from fractions import Fraction as F
from math import gcd, isqrt
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


def make_table(q, radius):
    table = {}
    for z in range(q):
        candidates = []
        for i in range(1, radius+1):
            for j in range(-radius, 2):
                if gcd(j+i*z, q) != 1:
                    continue
                lo = max(F(0), F(-1-j, i))
                hi = min(F(1, 2), F(1-j, i))
                if lo <= hi:
                    candidates.append((lo, hi, i, j))
        reach = F(0)
        intervals = []
        while reach < F(1, 2):
            allowed = [x for x in candidates if x[0] <= reach and x[1] > reach]
            need(bool(allowed), f"continuous table coverage q{q},z{z}")
            row = max(allowed, key=lambda x: (x[1], -x[2], -abs(x[3])))
            intervals.append(row)
            reach = row[1]
        table[z] = intervals
    return table
