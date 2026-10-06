"""E030: three analytically forced N6 exceptions, generic two-prime upper heights."""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import sys
K, SCALE = (16, 1 << 64)

def require(ok, message):
    if not ok:
        raise ValueError(message)

def ceildiv(a, b):
    return -(-a // b)

def atanh_log(t, terms):
    require(F(0) <= t <= F(1, 3), 'range-reduced atanh argument')
    lower, power = (F(0), t)
    for j in range(terms):
        lower += 2 * power / (2 * j + 1)
        power *= t * t
    return (lower, lower + 2 * power / ((2 * terms + 1) * (1 - t * t)))

def log_interval(n, terms=K):
    if type(n) is not int or n < 1:
        raise ValueError('log input must be a positive exact integer')
    k = n.bit_length() - 1
    t = F(n - (1 << k), n + (1 << k))
    lo, hi = atanh_log(t, terms)
    l2, u2 = atanh_log(F(1, 3), terms)
    return (k * l2 + lo, k * u2 + hi)

def log_ticks(n):
    lo, hi = log_interval(n)
    return (lo.numerator * SCALE // lo.denominator, ceildiv(hi.numerator * SCALE, hi.denominator))

def certify(hist, m):
    require(sum(hist.values()) == m and all((type(h) is int and 1 <= h <= m // 2 and (type(c) is int) and (c > 0) for h, c in hist.items())), 'histogram mass and range')
    upper_sum = sum((c * log_ticks(h)[1] for h, c in hist.items()))
    lower_m = log_ticks(m)[0]
    surplus = 2 * upper_sum - m * lower_m
    require(surplus < 0, 'finite upper-mean sign')
    return dict(log_scale=SCALE, upper_centered_log_numerator=surplus, upper_centered_log_denominator=2 * SCALE * m, strict_negative=True, histogram=[[h, c] for h, c in sorted(hist.items())])

def mobius_table(limit):
    mu, primes, composite = ([0] * (limit + 1), [], [False] * (limit + 1))
    mu[1] = 1
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for p in primes:
            if n * p > limit:
                break
            composite[n * p] = True
            if n % p == 0:
                mu[n * p] = 0
                break
            mu[n * p] = -mu[n]
    return mu

def central(m, rad):
    l = isqrt(m)
    heights = [0] * m
    heights[0] = 1
    units = [b for b in range(1, l + 1) if gcd(b, rad) == 1]
    pairs = 0
    for b in units:
        inv = pow(b, -1, m)
        for u in range(1, l + 1):
            if gcd(u, b) != 1:
                continue
            pairs += 1
            h = max(u, b)
            r = u * inv % m
            require((u - r * b) % m == 0, 'positive central witness')
            require((-u - -r % m * b) % m == 0, 'negative central witness')
            for residue in (r, -r % m):
                if not heights[residue] or h < heights[residue]:
                    heights[residue] = h
    collision_bins = Counter()
    for b in units:
        for d in units:
            if l * (b + d) < m or gcd(b, d) != 1:
                continue
            u0 = m * pow(d, -1, b) % b if b > 1 else 0
            lower = max(1, ceildiv(m - l * b, d))
            upper = min(l, (m - b) // d)
            u = u0 + b * ceildiv(lower - u0, b)
            if u > upper:
                continue
            e = (m - u * d) // b
            require(1 <= u <= l and 1 <= e <= l and (u * d + e * b == m), 'central collision diophantine equation')
            require(u + b > upper, 'one central collision per ordered pair')
            collision_bins[max(u, e, b, d)] += 1
    mu = mobius_table(l)
    unit_cumulative = [0] * (l + 1)
    for t in range(1, l + 1):
        unit_cumulative[t] = unit_cumulative[t - 1] + (gcd(t, rad) == 1)
    hist = Counter((h for h in heights if h))
    actual_cdf = collisions = 0
    for t in range(1, l + 1):
        p = sum((mu[g] * (t // g) * unit_cumulative[t // g] for g in range(1, t + 1) if gcd(g, rad) == 1))
        collisions += collision_bins[t]
        actual_cdf += hist[t]
        require(actual_cdf == 1 + 2 * p - collisions, 'every exact central histogram bin')
    require(p == pairs, 'primitive enumeration versus Mobius count')
    return (heights, hist, dict(central_cutoff=l, primitive_positive_pairs=pairs, central_residues=sum(hist.values()), central_thresholds=l, ordered_collisions=sum(collision_bins.values())))

def envelope(m, rad):
    exact, exact_hist, evidence = central(m, rad)
    t = isqrt(2 * m)
    t += 1
    require(t * t > 2 * m, 'strict far-cover threshold')
    tail = [0] * m
    hole_count = incidences = candidate_checks = 0
    for q in range(2, (m - 1) // t + 1):
        if gcd(q, rad) == 1:
            continue
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            inv = pow(a, -1, q)
            bs = []
            for v in (inv, q - inv):
                b = v + q * ((t - v) // q)
                if gcd(b, rad) != 1:
                    b -= q
                require(0 < b <= t and gcd(b, rad) == 1 and (b > t - 2 * q), 'low-center enclosing unit neighbor')
                if q % rad == 0:
                    require(b > t - q, 'radical-divisible refinement')
                bs.append(b)
            left = F(m * a, q) - F(m - q * t, q * bs[0])
            right = F(m * a, q) + F(m - q * t, q * bs[1])
            lo, hi = (left.numerator // left.denominator + 1, ceildiv(right.numerator, right.denominator) - 1)
            require(0 <= lo and hi < m, 'unclipped open integer hole')
            if lo > hi:
                continue
            hole_count += 1
            for r in range(lo, hi + 1):
                incidences += 1
                delta = m * a - q * r
                d = abs(delta)
                v = inv if delta >= 0 else q - inv
                b0 = v + q * ((m - v * (q + d)) // (q * (q + d)))
                candidates = []
                below = b0
                if below > 0 and gcd(below, rad) != 1:
                    below -= q
                if below > 0:
                    candidates.append(below)
                above = b0 + q
                if gcd(above, rad) != 1:
                    above += q
                require(above > 0, 'positive upper bracketing denominator')
                candidates.append(above)
                best = max(1, min(r, m - r))
                for b in candidates:
                    require(b > 0 and gcd(b, rad) == 1, 'positive unit tail denominator')
                    numerator = m - d * b
                    require(numerator % q == 0, 'integer determinant-one numerator')
                    num = numerator // q
                    if delta < 0:
                        num = -num
                    require((num - r * b) % m == 0, 'tail witness congruence')
                    best = min(best, max(abs(num), b))
                    candidate_checks += 1
                if not tail[r] or best < tail[r]:
                    tail[r] = best
    upper = [0] * m
    for r in range(m):
        if exact[r]:
            upper[r] = exact[r]
        elif tail[r]:
            require(tail[r] > evidence['central_cutoff'], 'central completeness versus tail witness')
            upper[r] = tail[r]
        else:
            upper[r] = min(t, max(1, min(r, m - r)))
    hist = Counter(upper)
    evidence.update(far_threshold=t, nonempty_holes=hole_count, integer_hole_incidences=incidences, tail_witness_checks=candidate_checks, distinct_in_holes=sum((h > 0 for h in tail)), exact_minima_scope='central window only')
    require({h: c for h, c in hist.items() if h <= evidence['central_cutoff']} == dict(exact_hist), 'full envelope preserves exact central histogram')
    return (upper, hist, exact, evidence)

def defining_minimum(m, r, rad):
    best = max(1, min(r, m - r))
    for b in range(1, best + 1):
        if b >= best:
            break
        if gcd(b, rad) != 1:
            continue
        residue = r * b % m
        best = min(best, max(b, min(residue, m - residue)))
    return best

def farey(t):
    a, b, c, d = (0, 1, 1, t)
    yield (a, b)
    while True:
        yield (c, d)
        if c == d == 1:
            break
        k = (t + b) // d
        a, b, c, d = (c, d, k * c - a, k * d - b)

def farey_cdf(m, rad, t):
    intervals = sorted(((max(0, ceildiv(m * a - t, b)), min(m - 1, (m * a + t) // b)) for a, b in farey(t) if gcd(b, rad) == 1))
    covered, left, right = (0, 0, -1)
    for lo, hi in intervals:
        if lo > hi:
            continue
        if lo > right + 1:
            covered += max(0, right - left + 1)
            left, right = (lo, hi)
        else:
            right = max(right, hi)
    return covered + max(0, right - left + 1)

def premise_tests():
    rejected = 0
    for value in (0, -1, True, 1.0, F(1)):
        try:
            log_interval(value)
        except ValueError:
            rejected += 1
        else:
            raise ValueError('invalid exact-log input accepted')
    for n in (1, 2, 3, 6, 36, 1296):
        lo, hi = log_interval(n)
        finer_lo, finer_hi = log_interval(n, 32)
        lower, upper = log_ticks(n)
        require(F(lower, SCALE) <= lo <= finer_lo <= finer_hi <= hi <= F(upper, SCALE), 'nested outward logarithm bounds')
    require(F(1691, 14400) > F(1, 10) and F(263, 4500) > F(1, 20), 'parameter margins')
    error3000 = (22 + 15 * F(81, 10)) / 3000 + (18 + 4 * F(81, 10)) / 3000 ** 2
    require(F(1, 20) - error3000 == F(2026, 937500) > 0, 'even-other finite error')
    error2048 = (22 + 15 * F(77, 10)) / 2048 + (18 + 4 * F(77, 10)) / 2048 ** 2
    require(error2048 < F(7, 100), 'odd finite error')
    require(6 ** 9 > 9000000 and 6 ** 9 > 2 ** 22 and (15 ** 6 > 6 ** 9) and (12 ** 7 > 6 ** 9) and (10 ** 7 > 9000000), 'uniform domain reduction')
    small = []
    for b in range(2, 15):
        factors, value = ([], b)
        for p in (2, 3, 5, 7, 11, 13):
            if value % p == 0:
                factors.append(p)
                while value % p == 0:
                    value //= p
        require(value == 1, 'small integer prime cover')
        if len(factors) == 2:
            small.append(b)
    require(small == [6, 10, 12, 14], 'exact small-base exceptions')
    return dict(invalid_log_inputs_rejected=rejected, exact_integer_log_cases=6, atanh_terms=K, scale=SCALE, even_other_error_at_3000=str(error3000), odd_error_at_2048=str(error2048), small_bases_with_two_primes=small, analytic_lemmas_proved_by_checker=False)
