"""Uniform negative contribution with unit-aware sieve and integer-radius envelope.

See PROOF.md. All decisions are integer/rational. Full replay regenerates every
finite support/radius value and the complete certificate; the public API reads
only the certified compressed envelope, authenticated by the component manifest.
"""
from array import array
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import gcd, isqrt, prod
from pathlib import Path
import argparse
import hashlib
import json
import time
ROOT = Path(__file__).resolve().parent
KERNEL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)
NEXT_PRIMES = (23, 29, 31, 37, 41, 43)
COMMON_CUTOFF = 199
PRIME_SQUARE_UPPER = F(4523, 10000)
COEFFICIENTS = {4: F(11, 50), 5: F(1, 5), 6: F(37, 200)}
MIN_Q = 2
DECIMAL_SCALE = 10 ** 36

def need(condition, message):
    if not condition:
        raise ArithmeticError(message)

def ceildiv(a, b):
    need(b > 0, 'positive denominator')
    return -(-a // b)

def tick_lower(value, scale=DECIMAL_SCALE):
    value = F(value)
    return F(value.numerator * scale // value.denominator, scale)

def tick_upper(value, scale=DECIMAL_SCALE):
    value = F(value)
    return F(ceildiv(value.numerator * scale, value.denominator), scale)

def decimal_lower(value, digits=15):
    value = F(value)
    count = value.numerator * 10 ** digits // value.denominator
    sign = '-' if count < 0 else ''
    count = abs(count)
    return sign + str(count // 10 ** digits) + '.' + str(count % 10 ** digits).zfill(digits)

@lru_cache(None)
def log_interval(value):
    """atanh series, 48 terms; geometric positive remainder, outward rounded."""
    value = F(value)
    need(value > 0, 'positive logarithm argument')
    if value < 1:
        lo, hi = log_interval(1 / value)
        return (-hi, -lo)
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if value < F(2) ** exponent:
        exponent -= 1
    reduced = value / F(2) ** exponent
    need(1 <= reduced < 2 and exponent >= 0, 'log range reduction')

    def series(x):
        z = (x - 1) / (x + 1)
        z2 = z * z
        power = z
        total = F(0)
        for j in range(48):
            total += 2 * power / (2 * j + 1)
            power *= z2
        return (total, total + 2 * power / (97 * (1 - z2)))
    lo, hi = series(reduced)
    l2, h2 = series(F(2))
    return (tick_lower(lo + exponent * l2), tick_upper(hi + exponent * h2))

def sqrt_interval(value):
    value = F(value)
    need(value >= 0, 'nonnegative square root')
    t = isqrt(value.numerator * DECIMAL_SCALE ** 2 // value.denominator)
    lo, hi = (F(t, DECIMAL_SCALE), F(t + 1, DECIMAL_SCALE))
    need(lo * lo <= value < hi * hi, 'integer square root enclosure')
    return (lo, hi)

def sieve(n):
    prime = bytearray(b'\x01') * (n + 1)
    prime[:2] = b'\x00\x00'
    for p in range(2, isqrt(n) + 1):
        if prime[p]:
            prime[p * p::p] = b'\x00' * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if prime[p]]

def prime_trial(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

def prime_factors(n):
    answer = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            answer.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        answer.append(n)
    return answer

def mask_product(mask):
    return prod((p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i))

def signed_divisors(primes):
    return [(prod(subset), (-1) ** size) for size in range(len(primes) + 1) for subset in combinations(primes, size)]

def square_density(mask):
    return prod((F(p - 1, p) if mask & 1 << i else F(p * p - 1, p * p) for i, p in enumerate(KERNEL_PRIMES)))

def square_error(mask):
    """One-sided floor-product error, keeping only positive Mobius terms."""
    count = mask.bit_count()
    plus_rem = prod((F(p + 1, p) for i, p in enumerate(KERNEL_PRIMES) if not mask & 1 << i))
    minus_rem = prod((F(p - 1, p) for i, p in enumerate(KERNEL_PRIMES) if not mask & 1 << i))
    plus_all = prod((F(p + 1, p) for p in KERNEL_PRIMES))
    minus_all = prod((F(p - 1, p) for p in KERNEL_PRIMES))
    return (2 ** count * plus_rem + (minus_rem if count == 0 else 0)) / 2 + (plus_all + minus_all) / 2

def square_error_direct(mask):
    marked = [p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i]
    unmarked = [p for p in KERNEL_PRIMES if p not in marked]
    return sum((F(1, k) + F(1, k * d) for k, mu_k in signed_divisors(unmarked) for d, mu_d in signed_divisors(marked) if mu_k * mu_d > 0), F(0))

def count_square_mobius(mask, R):
    marked = [p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i]
    unmarked = [p for p in KERNEL_PRIMES if p not in marked]
    return sum((mu_k * (R // k) * sum((mu_d * (R // (k * d)) for d, mu_d in signed_divisors(marked))) for k, mu_k in signed_divisors(unmarked)))

class Checks:

    def __init__(self):
        self.n = 0
        self.categories = Counter()

    def require(self, condition, message, category):
        self.n += 1
        self.categories[category] += 1
        need(condition, message)

def prime_square_certificate(checks):
    cutoff, scale = (20000, 10 ** 12)
    primes = sieve(cutoff)
    by_trial = [n for n in range(2, cutoff + 1) if prime_trial(n)]
    checks.require(primes == by_trial, 'independent primality enumeration', 'prime_square')
    upper_tick = sum((ceildiv(scale, p * p) for p in primes)) + ceildiv(scale, cutoff)
    bound = F(upper_tick, scale)
    checks.require(bound < PRIME_SQUARE_UPPER, 'prime-square series explicit upper', 'prime_square')
    checks.require(tuple(primes[:8]) == KERNEL_PRIMES and tuple(primes[8:14]) == NEXT_PRIMES, 'kernel and first permissible residual primes', 'prime_square')
    return dict(cutoff=cutoff, rounding_scale=scale, prime_count=len(primes), finite_sum_plus_integer_tail_upper=str(bound), claimed_all_prime_square_sum_upper=str(PRIME_SQUARE_UPPER), tail_proof='sum_{n>T} 1/n^2 <= integral_T^infty t^-2 dt = 1/T', primality_crosscheck='complete sieve versus trial division up to 20000')

def primitive_count(M, R):
    return 1 + 2 * sum((gcd(a, d) == 1 for d in range(1, R + 1) if gcd(d, M) == 1 for a in range(1, R + 1)))

def physical_negative_part(M):
    """Direct full-state/denominator definition; no inverses or moment identities."""
    upper_radius = isqrt(M - 1)
    denominators = [d for d in range(1, upper_radius + 1) if gcd(d, M) == 1]
    histogram = Counter()
    for s in range(M):
        best = upper_radius + 1
        for d in denominators:
            if d >= best:
                break
            residue = s * d % M
            radius = max(d, min(residue, M - residue))
            if radius < best:
                best = radius
        if best <= upper_radius:
            histogram[best] += 1
    lo = hi = F(0)
    for radius, count in histogram.items():
        lower, upper = log_interval(F(M, radius * radius))
        lo += F(count, 2 * M) * lower
        hi += F(count, 2 * M) * upper
    return (histogram, tick_lower(lo), tick_upper(hi))

def alternative_log_interval(value):
    """Independent 160-term Mercator series after repeated halving."""
    value = F(value)
    need(value > 0, 'positive alternate log')
    if value < 1:
        lo, hi = alternative_log_interval(1 / value)
        return (-hi, -lo)
    exponent = 0
    while value >= 2:
        value /= 2
        exponent += 1

    def series(x):
        z = 1 - 1 / x
        power = z
        total = F(0)
        for j in range(1, 161):
            total += power / j
            power *= z
        return (total, total + power / (161 * (1 - z)))
    lo, hi = series(value)
    l2, h2 = series(F(2))
    return (lo + exponent * l2, hi + exponent * h2)

def unit_density(mask):
    return prod((F(p - 1, p) for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i))

def unit_upper_error(mask):
    size = mask.bit_count()
    return 2 ** (size - 1) if size else 0

def unit_count_mobius(mask, R):
    return sum((sign * (R // d) for d, sign in signed_divisors([p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i])))

def integration_coefficients(profile, Q):
    """Exact A,B in L(Q)=A+B*ln(2); compressed finite bands plus infinite plateau."""
    need(isinstance(Q, int) and Q >= MIN_Q, 'integer scalar Q >= 2')
    breaks = [(int(row[0]), int(row[1])) for row in profile['breakpoints_radius_tick']]
    need(breaks[0] == (0, 0), 'no vector at radius zero')
    radius = isqrt((Q - 1) // 2)
    total_tick = 0
    endpoint = 0
    for i, (left, tick) in enumerate(breaks):
        if left > radius:
            break
        right = breaks[i + 1][0] if i + 1 < len(breaks) else radius + 1
        end = min(right, radius)
        if end > left:
            total_tick += tick * (end * end - left * left)
        endpoint = tick
    rational = F(total_tick, 2 * Q * DECIMAL_SCALE) + F(endpoint, DECIMAL_SCALE) * (F(1, 4) - F(radius * radius, 2 * Q))
    return dict(rational=rational, log2_coefficient=F(endpoint, 4 * DECIMAL_SCALE), strict_radius=radius, endpoint=F(endpoint, DECIMAL_SCALE), weighted_prefix_tick=total_tick, tick_scale=DECIMAL_SCALE)

def lower_from_profile(profile, Q, alternate=False):
    """Integrate the nondecreasing integer-radius envelope, including strict edge."""
    need(isinstance(Q, int) and Q >= MIN_Q, 'integer scalar Q >= 2')
    if not alternate:
        terms = integration_coefficients(profile, Q)
        l2, h2 = log_interval(2)
        low = terms['rational'] + terms['log2_coefficient'] * l2
        high = terms['rational'] + terms['log2_coefficient'] * h2
        return (tick_lower(low), tick_upper(high))
    radius = isqrt(Q // 2)
    if 2 * radius * radius == Q:
        radius -= 1
    need(2 * radius * radius < Q <= 2 * (radius + 1) ** 2, 'independent strict square boundary')
    integral = F(0)
    end_value = F(0)
    breaks = [(int(row[0]), int(row[1])) for row in profile['breakpoints_radius_tick']]
    for i, (start, tick) in enumerate(breaks):
        value = F(tick, DECIMAL_SCALE)
        stop = breaks[i + 1][0] if i + 1 < len(breaks) else None
        bottom = F(start * start, Q)
        top = F(1, 2) if stop is None else min(F(stop * stop, Q), F(1, 2))
        if bottom < top:
            integral += value * (top - bottom) / 2
        if start <= radius and (stop is None or radius < stop):
            end_value = value
    l2, h2 = alternative_log_interval(2)
    return (integral + end_value * l2 / 4, integral + end_value * h2 / 4)

@lru_cache(None)
def certified_profiles():
    data = json.loads((ROOT / 'certificate.json').read_text())
    return data['envelopes']

def lower_interval(omega, Q):
    need(omega in COEFFICIENTS, 'proved prime count is 4, 5 or 6')
    return lower_from_profile(certified_profiles()[str(omega)], Q)

def lower_bound(omega, Q):
    return lower_interval(omega, Q)[0]

def decompose(omega, Q):
    """Expose the exact rational and log(2) terms for independent global checks."""
    need(omega in COEFFICIENTS, 'proved prime count is 4, 5 or 6')
    return integration_coefficients(certified_profiles()[str(omega)], Q)

def sieve_certificate(checks):
    prime_square = prime_square_certificate(checks)
    medium = [p for p in sieve(COMMON_CUTOFF) if p > KERNEL_PRIMES[-1]]
    medium_square = sum((F(1, p * p) for p in medium), F(0))
    medium_harmonic = sum((F(1, p) for p in medium), F(0))
    tail_square = PRIME_SQUARE_UPPER - sum((F(1, p * p) for p in sieve(COMMON_CUTOFF)), F(0))
    checks.require(tail_square > 0, 'positive proved prime-square tail after 199', 'asymptotic')
    rows = []
    for omega, coefficient in COEFFICIENTS.items():
        for mask in range(256):
            if mask.bit_count() > omega:
                continue
            residual = omega - mask.bit_count()
            rho = unit_density(mask)
            density = square_density(mask)
            outside = sum((F(1, p) for p in NEXT_PRIMES[:residual]), F(0))
            gamma = density - rho * (outside + medium_square) - tail_square
            error = square_error(mask) + unit_upper_error(mask) * (residual + medium_harmonic)
            checks.require(square_error(mask) == square_error_direct(mask), 'one-sided square error equals independent signed sum', 'asymptotic')
            checks.require(gamma > coefficient, 'positive coefficient margin', 'asymptotic')
            ratio = error / (gamma - coefficient)
            cutoff = max(4, ceildiv(ratio.numerator, ratio.denominator))
            checks.require(gamma - error / cutoff >= coefficient, 'every real radius at or above analytic cutoff is covered', 'asymptotic')
            rows.append(dict(omega=omega, mask=mask, marked_kernel_primes=[p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i], residual_prime_count=residual, unit_density=str(rho), square_density=str(density), gamma=str(gamma), floor_error=str(error), coefficient=str(coefficient), analytic_cutoff=cutoff))
    cutoffs = {w: max((row['analytic_cutoff'] for row in rows if row['omega'] == w)) for w in COEFFICIENTS}
    max_radius = max(cutoffs.values())
    support_lists = {w: [k for k in range(256) if k.bit_count() <= w] for w in COEFFICIENTS}
    active_masks = support_lists[max(COEFFICIENTS)]
    units = [array('I', [0]) for _ in range(256)]
    squares = [0] * 256
    medium_loss = [0] * 256
    medium_divisors = [[] for _ in range(max_radius + 1)]
    tail_loss = [0] * (max_radius + 1)
    primes = sieve(max_radius)
    for p in primes:
        if KERNEL_PRIMES[-1] < p <= COMMON_CUTOFF:
            for R in range(p, max_radius + 1, p):
                medium_divisors[R].append(p)
        elif p > COMMON_CUTOFF:
            for R in range(p, max_radius + 1, p):
                tail_loss[R] += 2 * (R // p) - 1
    for R in range(1, max_radius + 1):
        tail_loss[R] += tail_loss[R - 1]
    finite_minima = {w: [] for w in COEFFICIENTS}
    minimum_digests = {w: hashlib.sha256() for w in COEFFICIENTS}
    comparison_radii = set(range(1, 25)) | {28, 32, 40, 57, 80, 100, 136, 199, 220, 257, 586, 997, max_radius - 1, max_radius} | set(cutoffs.values())
    support_radius_evaluations = {w: 0 for w in COEFFICIENTS}
    small_floor_gates = 0
    for R in range(1, max_radius + 1):
        divisor_mask = sum((1 << i for i, p in enumerate(KERNEL_PRIMES) if R % p == 0))
        old_units = units[divisor_mask][-1]
        for mask in range(256):
            units[mask].append(units[mask][-1] + int(mask & divisor_mask == 0))
        for mask in active_masks:
            squares[mask] += units[mask | divisor_mask][-1] + (old_units if mask & divisor_mask == 0 else 0)
            for p in medium_divisors[R]:
                quotient = R // p
                medium_loss[mask] += quotient * units[mask][quotient] - (quotient - 1) * units[mask][quotient - 1]
        if R in comparison_radii:
            direct_tail = sum(((R // p) ** 2 for p in primes if COMMON_CUTOFF < p <= R))
            checks.require(tail_loss[R] == direct_tail, 'prime-square floor tail recurrence', 'count_crosscheck')
            for mask in active_masks:
                checks.require(units[mask][R] == unit_count_mobius(mask, R), 'unit-count recurrence versus Mobius formula', 'count_crosscheck')
                checks.require(squares[mask] == count_square_mobius(mask, R), 'primitive-kernel square recurrence versus Mobius formula', 'count_crosscheck')
                direct_medium = sum((R // p * unit_count_mobius(mask, R // p) for p in medium))
                checks.require(medium_loss[mask] == direct_medium, 'unit-aware common-prime loss versus direct Mobius calculation', 'count_crosscheck')
        if R <= 24:
            for mask in active_masks:
                q = mask_product(mask)
                direct = sum((gcd(d, q) == 1 and gcd(gcd(a, d), prod(KERNEL_PRIMES)) == 1 for d in range(1, R + 1) for a in range(1, R + 1)))
                checks.require(squares[mask] == direct, 'independent two-dimensional kernel predicate', 'direct_definition')
                rho, err = (unit_density(mask), unit_upper_error(mask))
                checks.require(units[mask][R] <= rho * R + err, 'unit-count upper error at small integer radii', 'direct_definition')
                small_floor_gates += 1
        for omega in COEFFICIENTS:
            if R >= cutoffs[omega]:
                continue
            minimum, worst = (R * R, 0)
            for mask in support_lists[omega]:
                count = omega - mask.bit_count()
                outside_loss = R * sum((units[mask][R // p] for p in NEXT_PRIMES[:count]))
                bound = squares[mask] - outside_loss - medium_loss[mask] - tail_loss[R]
                if bound < minimum:
                    minimum, worst = (bound, mask)
                support_radius_evaluations[omega] += 1
            final = max(R, minimum)
            finite_minima[omega].append(final)
            minimum_digests[omega].update(f'{R},{minimum},{worst},{final}\n'.encode())
    envelopes = {}
    for omega, coefficient in COEFFICIENTS.items():
        cutoff = cutoffs[omega]
        checks.require(len(finite_minima[omega]) == cutoff - 1, 'every radius from one to cutoff minus one is present', 'inventory')
        checks.require(support_radius_evaluations[omega] == (cutoff - 1) * len(support_lists[omega]), 'complete support by radius Cartesian product', 'inventory')
        checks.require(2 * coefficient * cutoff > 1, 'analytic tail ratio is increasing throughout the infinite tail', 'asymptotic')
        tail_ratio = F(1 + 2 * coefficient * cutoff * cutoff, (cutoff + 1) ** 2)
        values = [0] * (cutoff + 1)
        suffix = tail_ratio
        values[cutoff] = tick_lower(suffix).numerator * (DECIMAL_SCALE // tick_lower(suffix).denominator)
        for R in range(cutoff - 1, 0, -1):
            candidate = F(1 + 2 * finite_minima[omega][R - 1], (R + 1) ** 2)
            suffix = min(suffix, candidate)
            values[R] = suffix.numerator * DECIMAL_SCALE // suffix.denominator
        checks.require(all((values[R] <= values[R + 1] for R in range(cutoff))), 'nondecreasing integer-radius lower envelope', 'envelope')
        changes = [[0, 0]]
        for R in range(1, cutoff + 1):
            if values[R] != values[R - 1]:
                changes.append([R, values[R]])
        envelopes[str(omega)] = dict(coefficient=str(coefficient), tail_begin=cutoff, tail_ratio_exact=str(tail_ratio), tick_scale=DECIMAL_SCALE, finite_radius_count=cutoff - 1, support_case_count=len(support_lists[omega]), support_radius_evaluations=support_radius_evaluations[omega], complete_finite_minimum_trace_sha256=minimum_digests[omega].hexdigest(), envelope_breakpoint_count=len(changes), breakpoints_radius_tick=changes, initial_exact_minima=finite_minima[omega][:128], tail_constant_starts=changes[-1][0], tail_is_infinite=True, support_extremizer_assumed=False)
    return (dict(kernel_primes=list(KERNEL_PRIMES), kernel_product=prod(KERNEL_PRIMES), first_residual_primes=list(NEXT_PRIMES), unit_aware_common_primes=medium, common_prime_cutoff=COMMON_CUTOFF, prime_square_series=prime_square, prime_square_tail_after_cutoff_upper=str(tail_square), largest_computed_radius=max_radius, support_case_counts={str(w): len(support_lists[w]) for w in COEFFICIENTS}, analytic_cutoffs={str(w): cutoffs[w] for w in COEFFICIENTS}, support_radius_evaluations={str(w): support_radius_evaluations[w] for w in COEFFICIENTS}, count_crosscheck_radii=sorted(comparison_radii), rows=rows), envelopes)

def definition_regressions(checks, envelopes):
    support_count = 0
    radii = (12, 20, 28, 40, 57)
    for omega in COEFFICIENTS:
        bounds = envelopes[str(omega)]['initial_exact_minima']
        for mask in range(256):
            if mask.bit_count() > omega:
                continue
            support = [p for i, p in enumerate(KERNEL_PRIMES) if mask & 1 << i] + list(NEXT_PRIMES[:omega - mask.bit_count()])
            modulus = prod(support)
            checks.require(len(prime_factors(modulus)) == omega, 'representative has exactly the declared support size', 'primitive_definition')
            for R in radii:
                count = primitive_count(modulus, R)
                checks.require(count >= 1 + 2 * bounds[R - 1], 'finite all-support lower count below direct primitive definition', 'primitive_definition')
            support_count += 1
    power_rows = []
    for B in (2210, 2622, 2310, 6270, 6510, 30030):
        for R in (12, 28, 57):
            value = primitive_count(B, R)
            checks.require(value == primitive_count(B * B, R), 'full-unit primitive condition unaffected by prime powers', 'prime_powers')
            power_rows.append(dict(B=B, N=2, M=B * B, R=R, A=value))
    small_orbit_count = 0
    boundary_failures = []
    for M in range(2, 81):
        for R in sorted({1, max(1, isqrt((M - 1) // 2)), max(1, isqrt(M // 2))}):
            states = Counter((a * pow(d, -1, M) % M for d in range(1, R + 1) if gcd(d, M) == 1 for a in range(-R, R + 1) if gcd(a, d) == 1))
            physical = sum((any((min(s * d % M, M - s * d % M) <= R for d in range(1, R + 1) if gcd(d, M) == 1)) for s in range(M)))
            checks.require(len(states) == physical, 'primitive states agree with full physical CDF including zero', 'finite_orbit')
            if 2 * R * R < M:
                checks.require(sum(states.values()) == physical, 'strict determinant collision-free gate', 'finite_orbit')
            if 2 * R * R == M and max(states.values()) > 1:
                boundary_failures.append(dict(M=M, R=R, maximum_fibre=max(states.values())))
            small_orbit_count += 1
    checks.require(bool(boundary_failures), 'closed square boundary has explicit counterexamples', 'finite_orbit')
    physical_rows = []
    for M in (2210, 2622, 2310, 6270, 6510, 30030):
        omega = len(prime_factors(M))
        histogram, lo, hi = physical_negative_part(M)
        proved, unused_upper = lower_from_profile(envelopes[str(omega)], M)
        checks.require(lo >= proved, 'new uniform lower bound below complete physical negative part', 'negative_definition')
        checks.require(histogram[1] >= 1, 'zero numerator is included', 'negative_definition')
        physical_rows.append(dict(M=M, omega=omega, complete_state_count=M, direct_negative_lower=str(lo), direct_negative_upper=str(hi), uniform_lower=str(proved), histogram={str(r): histogram[r] for r in sorted(histogram)}))
    return dict(representative_supports=support_count, representative_radii=list(radii), representative_examples_are_not_uniform_premises=True, prime_power_rows=power_rows, small_complete_orbits=small_orbit_count, forbidden_boundary_examples=boundary_failures, direct_physical_negative_examples=physical_rows)

def build_certificate():
    checks = Checks()
    sieve_data, envelopes = sieve_certificate(checks)
    regressions = definition_regressions(checks, envelopes)
    uniform_rows = []
    declared = ((4, 2210), (4, 2622), (5, 6270), (5, 6510), (6, 30030))
    for omega, Q in declared:
        profile = envelopes[str(omega)]
        terms = integration_coefficients(profile, Q)
        lo, hi = lower_from_profile(profile, Q)
        ilo, ihi = lower_from_profile(profile, Q, alternate=True)
        checks.require(lo <= ilo <= ihi <= hi, 'independent band integration and 160-term Mercator enclosure', 'integration')
        uniform_rows.append(dict(omega=omega, Q=Q, lower=str(lo), upper=str(hi), lower_decimal_truncated=decimal_lower(lo), strict_radius=terms['strict_radius'], rational_term=str(terms['rational']), log2_coefficient=str(terms['log2_coefficient']), endpoint_h=str(terms['endpoint']), weighted_prefix_tick=terms['weighted_prefix_tick'], independent_integral_crosscheck=True))
    boundary_Q = {2, 3, 4, 5, 2210, 2622, 6270, 6510, 30030, 10 ** 80}
    for r in (1, 2, 10, 31, 55, 57, 80, 122, 1000):
        boundary_Q.update((2 * r * r - 1, 2 * r * r, 2 * r * r + 1))
    for omega, profile in envelopes.items():
        for Q in sorted((q for q in boundary_Q if q >= 2)):
            lo, hi = lower_from_profile(profile, Q)
            ilo, ihi = lower_from_profile(profile, Q, alternate=True)
            checks.require(lo <= ilo <= ihi <= hi, 'strict boundary and large-scalar independent integration', 'integration')
            next_lo, next_hi = lower_from_profile(profile, Q + 1)
            checks.require(next_lo >= lo, 'boundary regression of proved monotonicity', 'integration')
    return dict(status='PASS_UNIFORM_NEGATIVE_PART_FOUR_FIVE_SIX_PRIMES', theorem_scope='every integer M>=Q>=2, omega(M) in {4,5,6}; arbitrary positive valuations; scalar Q integer', auxiliary_envelope_scope='pointwise monotonicity in the real comparison parameter x>=Q', all_states_including_zero=True, zero_numerator_allowed=True, physical_norm='infinity', logarithms='natural', full_mark=True, normalization='Delta(B,N)=mu(B^N)/ln(B); never divide by ln(M) when N>1', finite_examples_extrapolated=False, primorial_extremality_assumed=False, external_numerical_or_asymptotic_premises=False, mathematical_novelty_claimed=False, uniform_rows=uniform_rows, envelopes=envelopes, sieve=sieve_data, regressions=regressions, checks=checks.n, check_categories=dict(sorted(checks.categories.items())))

def manifest_check():
    path = ROOT / 'MANIFEST_SHA256.json'
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    inventory = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p.name != 'MANIFEST_SHA256.json' and ('__pycache__' not in p.parts)}
    need(inventory == set(data['files']), 'negative-component exact inventory')
    for name, info in data['files'].items():
        raw = (ROOT / name).read_bytes()
        need(len(raw) == info['bytes'] and hashlib.sha256(raw).hexdigest() == info['sha256'], 'negative-component source hash: ' + name)
    return data

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-certificate', action='store_true')
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    need(ROOT not in args.report.resolve().parents, 'report must stay outside frozen component')
    start = time.monotonic()
    before = manifest_check()
    certificate = build_certificate()
    raw = (json.dumps(certificate, indent=2, sort_keys=True) + '\n').encode()
    if args.write_certificate:
        need(before is None, 'cannot overwrite frozen certificate')
        (ROOT / 'certificate.json').write_bytes(raw)
    else:
        need((ROOT / 'certificate.json').read_bytes() == raw, 'complete deterministic certificate replay')
    after = manifest_check()
    need(before == after, 'frozen negative-component bytes unchanged')
    result = dict(status='PASS', certificate_status=certificate['status'], certificate_sha256=hashlib.sha256(raw).hexdigest(), checks=certificate['checks'], check_categories=certificate['check_categories'], support_radius_evaluations=certificate['sieve']['support_radius_evaluations'], uniform_rows=certificate['uniform_rows'], frozen_inventory_unchanged=before is not None, payload_file_count=len(before['files']) if before else 0, runtime_seconds=round(time.monotonic() - start, 6))
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', newline='\n')
    print(json.dumps(result, sort_keys=True))
if __name__ == '__main__':
    main()
