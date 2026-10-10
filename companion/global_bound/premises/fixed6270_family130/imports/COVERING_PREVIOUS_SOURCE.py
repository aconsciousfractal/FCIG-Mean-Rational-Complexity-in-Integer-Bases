#!/usr/bin/env python3
"""Uniform forbidden-class covering via an exact small-prime sieve.

Every integer interval and every prime support of the declared size is covered.
No assumption that the worst support is a primorial is used.
"""
from pathlib import Path
from math import gcd, prod, isqrt
from itertools import combinations
from functools import lru_cache
import argparse
import hashlib
import json
import time

ROOT = Path(__file__).resolve().parent
CORE = (3, 5, 7, 11)
BOUNDS = {1: 2, 2: 4, 3: 6, 4: 10, 5: 14, 6: 22,
          7: 26, 8: 44, 9: 58, 10: 66, 11: 98, 12: 148}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()


def primes_sieve(limit):
    sieve = bytearray(b'\1') * (limit + 1)
    sieve[:2] = b'\0\0'
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = b'\0' * ((limit-p*p)//p + 1)
    return [p for p in range(2, limit+1) if sieve[p]]


def divisors_signed(support):
    return [(prod(sub), (-1)**len(sub))
            for size in range(len(support)+1)
            for sub in combinations(support, size)]


@lru_cache(None)
def interval_counts(support, length):
    """Intervals (s,s+length], all s modulo product(support)."""
    modulus = prod(support)
    values = [int(gcd(t, modulus) == 1) for t in range(modulus)]
    complete, remainder = divmod(length, modulus)
    base = complete * sum(values)
    current = base + sum(values[t % modulus] for t in range(1, remainder+1))
    primary = []
    for s in range(modulus):
        primary.append(current)
        current += values[(s+remainder+1) % modulus] - values[(s+1) % modulus]
    # A separate integer inclusion-exclusion computation retains every start,
    # not only the minimum. It works for incomplete and multiple periods.
    signed = divisors_signed(support)
    alternate = [sum(sign * ((s+length)//d - s//d) for d, sign in signed)
                 for s in range(modulus)]
    require(primary == alternate, 'all interval starts agree between two algorithms')
    require(current == primary[0], 'cyclic sliding count returns to its first value')
    require(sum(values) == prod(p-1 for p in support), 'exact totient of small support')
    return tuple(primary)


def build_certificate():
    require(__debug__, 'assertions must be enabled')
    ps = primes_sieve(101)
    trial = [p for p in range(2,102)
             if all(p % d for d in range(2,isqrt(p)+1))]
    require(ps == trial, 'sieve and trial primality agree')
    large = [p for p in ps if p > max(CORE)]
    require(len(large) >= max(BOUNDS), 'enough exact lower-bound primes')
    supports = [tuple(p for i,p in enumerate(CORE) if mask >> i & 1)
                for mask in range(1 << len(CORE))]
    require(len(set(supports)) == 16, 'all sixteen small-prime subsets')
    cases, checked_starts = [], 0
    # k=0 is the separate empty-support convention J(1)=1.
    for k, H in BOUNDS.items():
        require(H % 2 == 0, 'even H leaves exactly H/2 allowed parity positions')
        for has_two in (False, True):
            L = H//2 if has_two else H
            residual = k-int(has_two)
            for support in supports:
                if len(support) > residual:
                    continue
                missing = [p for p in CORE if p not in support]
                outside_count = residual-len(support)
                lower_primes = large[:outside_count]
                require(len(lower_primes) == outside_count,
                        'the support outside the core is counted exactly')
                vals = interval_counts(support, L)
                checked_starts += len(vals)
                minimum = min(vals)
                losses = [(L+p-1)//p for p in lower_primes]
                survivors = minimum-sum(losses)
                require(survivors >= 1, 'every support case has an allowed integer')
                row = dict(prime_count_at_most=k, H=H, contains_two=has_two,
                           indexed_interval_length=L, small_primes_present=list(support),
                           small_primes_absent=missing, small_modulus=prod(support),
                           all_interval_starts=len(vals), minimum_small_sieve_survivors=minimum,
                           first_minimizing_start=vals.index(minimum),
                           all_interval_counts_sha256=sha_bytes(canonical(list(vals))),
                           outside_prime_count=outside_count,
                           ordered_outside_prime_lower_bounds=lower_primes,
                           outside_union_upper_terms=losses,
                           surviving_positions_lower=survivors)
                cases.append(row)
    # Check the schema covers exactly the parity/subset cases, with no omissions.
    for k,H in BOUNDS.items():
        got = [r for r in cases if r['prime_count_at_most'] == k]
        expected = {(has_two, support) for has_two in (False, True)
                    for support in supports if len(support)+int(has_two) <= k}
        actual = {(r['contains_two'],tuple(r['small_primes_present'])) for r in got}
        require(actual == expected and len(got) == len(expected), 'complete unique case split')
        require(min(r['surviving_positions_lower'] for r in got) >= 1,
                'universal positive residual at every declared H')
    # Regression against direct forbidden-class masks for short prime supports.
    # CRT translations enumerate every choice of forbidden classes on the core.
    regression_checks = 0
    for support in [(), (3,), (5,), (3,5), (3,7), (3,5,7), (3,5,7,11)]:
        q = prod(support)
        for length in [1,2,7,11,26,74,148, q+3]:
            vals = interval_counts(support, length)
            for shift in range(q):
                direct = sum(all((j+shift) % p for p in support)
                             for j in range(1,length+1))
                require(vals[shift] == direct, 'direct forbidden-class translation check')
                regression_checks += 1
    return dict(status='PASS', theorem='J(q)<=H_k for every omega(q)<=k',
                convention='Every interval of H consecutive integers contains an allowed integer; J(1)=1',
                all_prime_supports=True, all_arbitrary_forbidden_classes=True,
                all_integer_intervals=True, arbitrary_prime_powers_irrelevant=True,
                small_core=list(CORE), maximum_small_modulus=prod(CORE),
                empty_support_H=1,
                bounds=[dict(prime_count_at_most=k,H=H) for k,H in BOUNDS.items()],
                cases=cases, case_count=len(cases),
                all_interval_starts_compared=checked_starts,
                direct_forbidden_class_regression_checks=regression_checks,
                primorial_extremality_assumed=False, optimality_claimed=False,
                external_new_theorem_used=False,
                proof='PROOF.md', exact_Cstar_equals_c0_proved=False)


def check_manifest():
    mf = ROOT/'MANIFEST_SHA256.json'
    entries = json.loads(mf.read_text())['files']
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and p != mf}
    require(actual == set(entries), 'complete component inventory')
    for rel, meta in entries.items():
        path = ROOT/rel
        require(path.resolve().is_relative_to(ROOT), 'contained component path')
        data = path.read_bytes()
        require(len(data) == meta['bytes'] and sha_bytes(data) == meta['sha256'],
                'immutable component bytes '+rel)
    return len(entries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', action='store_true', help='development only: write certificate before freezing')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    if args.report:
        require(not args.report.resolve().is_relative_to(ROOT), 'external report required')
    started = time.monotonic()
    count = None if args.build else check_manifest()
    cert = build_certificate()
    path = ROOT/'certificate.json'
    if args.build:
        path.write_bytes(canonical(cert))
    else:
        require(canonical(cert) == path.read_bytes(), 'complete regenerated certificate byte identity')
        require(check_manifest() == count, 'component remains unchanged')
    result = dict(status='PASS', payload_files=count,
                  certificate_sha256=sha_bytes(path.read_bytes()),
                  case_count=cert['case_count'],
                  all_interval_starts_compared=cert['all_interval_starts_compared'],
                  direct_forbidden_class_regression_checks=cert['direct_forbidden_class_regression_checks'],
                  bounds=cert['bounds'], component_bytes_unchanged=not args.build,
                  runtime_seconds=round(time.monotonic()-started,3))
    if args.report:
        args.report.write_bytes(canonical(result))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
