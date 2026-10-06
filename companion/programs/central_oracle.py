"""Exact original-model central counts by half-Farey/divisibility inversion."""
from functools import lru_cache
import hashlib
from math import gcd
R35_ANTECEDENT_PATH = 'results/batch_35/verify_large.py'
R35_ANTECEDENT_SHA256 = '13b44f61425cc6b3d978909a3b2b61be9a71cc63c2f1c9a08ecf5bda8d7bccbb'

def _need(ok, message):
    if not ok:
        raise ValueError(message)

@lru_cache(maxsize=None)
def _prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

def _profile(p, q, a, b):
    _need(all((type(v) is int for v in (p, q, a, b))), 'exact integer profile; bool rejected')
    _need(1 < p < q and a > 0 and (b > 0), 'ordered primes and positive valuations')
    _need(_prime(p) and _prime(q), 'prime premises')
    return p ** a * q ** b

def _cells(p, q, a, b):
    return sorted(((p ** i * q ** j, i, j) for i in range(a + 1) for j in range(b + 1)))

def invert_divisibility(p, q, a, b, prefix):
    """Invert a complete capped divisibility grid, INCLUDING terminal g=M.

    prefix is an exact dict mapping every divisor of p**a*q**b to a
    nonnegative integer cumulative population. The returned sorted cells
    describe that population's exact gcd classes. Arbitrary interval inputs,
    including zero, are supported; central geometry is not assumed here.
    """
    M = _profile(p, q, a, b)
    cells = _cells(p, q, a, b)
    _need(type(prefix) is dict, 'complete divisor-prefix dict required')
    _need(all((type(g) is int for g in prefix)), 'exact divisor keys; bool rejected')
    _need(set(prefix) == {g for g, i, j in cells}, 'complete saturated divisor grid')
    _need(all((type(v) is int and v >= 0 for v in prefix.values())), 'nonnegative exact prefix populations')
    result = []
    for g, i, j in cells:
        value = prefix[g]
        if i < a:
            value -= prefix[g * p]
        if j < b:
            value -= prefix[g * q]
        if i < a and j < b:
            value += prefix[g * p * q]
        _need(value >= 0, 'invalid cumulative grid: negative exact cell')
        result.append([g, value])
    _need(result[-1][0] == M and sum((value for g, value in result)) == prefix[1], 'saturated inversion partition')
    return result

def _half_farey(T, counters):
    """Yield complete reduced Farey order through 1/2, before unit filtering."""
    counters['half_Farey_emissions'] += 1
    yield (0, 1)
    if T == 1:
        return
    a, b, c, d = (0, 1, 1, T)
    while 2 * c <= d:
        _need(b * c - a * d == 1 and 1 <= d <= T and (b + d > T), 'Farey neighbor invariant')
        counters['half_Farey_emissions'] += 1
        yield (c, d)
        z = (T + b) // d
        counters['Farey_recurrence_floor_divisions'] += 1
        a, b, c, d = (c, d, z * c - a, z * d - b)
    _need(2 * a == b and 2 * c > d, 'complete half-order stopping boundary')
    counters['stopping_probes'] += 1

def central_counts(M, p, q, a, b, T):
    """Return exact CDF and ALL proper gcd strata at integer 1<=T<=sqrt(M).

    Original states are 0..M-1; physical denominators are positive and
    coprime to the full M; height is max(abs(A),d). Zero is covered once.
    No original-state or defining-height enumeration is used.
    """
    _need(type(M) is int and type(T) is int, 'exact M and T; bool rejected')
    expected_M = _profile(p, q, a, b)
    _need(M == expected_M, 'original modulus/valuation identity')
    _need(1 <= T and T * T <= M, 'central integer threshold scope')
    cells = _cells(p, q, a, b)
    divisors = [g for g, i, j in cells]
    small = [g for g in divisors if 1 < g <= T]
    prefix = {g: 0 for g in divisors}
    counters = dict(half_Farey_emissions=0, stopping_probes=0, half_boundary_floor_divisions=1, Farey_recurrence_floor_divisions=0, unit_gcd_tests=0, full_unit_centres=0, endpoint_floor_divisions=0, nonempty_integer_intervals=0, merged_intervals=0, divisor_queries=0, class_floor_divisions=0, divisibility_differences=0, divisibility_accumulator_updates=0, interval_mass_additions=0, proper_reconstruction_signed_references=0, reconstruction_additions_subtractions=0, terminal_grid_references=1, proper_stratum_doublings=0, proper_strata_output=(a + 1) * (b + 1) - 1, supported_divisor_cells=len(small), forced_zero_strata=sum((g > T for g in divisors[:-1])))
    half_end = (M - 1) // 2
    digest = hashlib.sha256()
    active = previous = None

    def finish(L, U):
        counters['merged_intervals'] += 1
        counters['interval_mass_additions'] += 1
        prefix[1] += U - L + 1
        digest.update((str(L) + ',' + str(U) + '\n').encode())
        for g in small:
            prefix[g] += U // g - (L - 1) // g
            counters['divisor_queries'] += 1
            counters['class_floor_divisions'] += 2
            counters['divisibility_differences'] += 1
            counters['divisibility_accumulator_updates'] += 1
    for c, d in _half_farey(T, counters):
        counters['unit_gcd_tests'] += 1
        if gcd(d, M) != 1:
            continue
        counters['full_unit_centres'] += 1
        L = max(1, -((T - M * c) // d))
        U = min(half_end, (M * c + T) // d)
        counters['endpoint_floor_divisions'] += 2
        if L > U:
            continue
        counters['nonempty_integer_intervals'] += 1
        if previous is not None:
            _need(L >= previous[0] and U >= previous[1], 'ordered clipped integer endpoints')
        previous = (L, U)
        if active is None:
            active = (L, U)
        elif L <= active[1] + 1:
            active = (active[0], max(active[1], U))
        else:
            finish(*active)
            active = (L, U)
    if active is not None:
        finish(*active)
    exact = invert_divisibility(p, q, a, b, prefix)
    _need(exact[-1] == [M, 0], 'positive terminal population must vanish')
    strata = []
    for (g, i, j), (exact_g, value) in zip(cells[:-1], exact[:-1]):
        _need(g == exact_g and (g <= T or value == 0), 'complete supported stratum alignment')
        terms = 1 + int(i < a) + int(j < b) + int(i < a and j < b)
        counters['proper_reconstruction_signed_references'] += terms
        counters['reconstruction_additions_subtractions'] += terms - 1
        if g <= T:
            strata.append([g, 2 * value])
            counters['proper_stratum_doublings'] += 1
        else:
            strata.append([g, 0])
    cdf = 1 + 2 * prefix[1]
    _need(1 + sum((count for g, count in strata)) == cdf, 'zero plus all proper strata mass')
    _need(counters['class_floor_divisions'] == 2 * len(small) * counters['merged_intervals'], 'exact merged-component query accounting')
    _need(counters['merged_intervals'] <= counters['nonempty_integer_intervals'] <= counters['full_unit_centres'] <= counters['half_Farey_emissions'], 'exact streaming counter hierarchy')
    counters['candidate_inspections'] = counters['half_Farey_emissions'] + counters['stopping_probes']
    counters['principal_floor_divisions'] = counters['Farey_recurrence_floor_divisions'] + counters['endpoint_floor_divisions'] + counters['class_floor_divisions']
    return dict(cdf=cdf, strata=strata, zero_count=1, counters=counters, interval_sha256=digest.hexdigest())
