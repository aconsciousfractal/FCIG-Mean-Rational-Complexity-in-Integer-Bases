"""Full independent physical raw consumer; no lattice search or minimality.

The raw directory is an explicit argument and can be a fresh master replay.
All source state weights and logarithm sums are independently reconstructed.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter
from math import gcd, isqrt
import argparse, csv, hashlib, json, struct
CHECKS = 0

def need(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)

def floor_fraction(x):
    x = F(x)
    return x.numerator // x.denominator

def ceil_fraction(x):
    x = F(x)
    return -(-x.numerator // x.denominator)

def prime(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

def prime_factors(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        out.append(n)
    return out

@lru_cache(None)
def fast_log_interval(x, scale=10 ** 40, terms=128):
    x = F(x)
    need(x > 0, 'positive independent logarithm argument')
    exponent = 0
    while x >= 2:
        x /= 2
        exponent += 1
    while x < 1:
        x *= 2
        exponent -= 1

    def series(t):
        numerator, denominator = (t.numerator, t.denominator)
        need(0 <= 2 * numerator <= denominator, 'independent Mercator convergence interval')
        lo = hi = scale
        lower = upper = 0
        for n in range(1, terms + 1):
            lo = lo * numerator // denominator
            hi = -(-hi * numerator // denominator)
            lower += lo // n
            upper += -(-hi // n)
        next_hi = -(-hi * numerator // denominator)
        tail = -(-next_hi * denominator // ((terms + 1) * (denominator - numerator)))
        return (lower, upper + tail)
    lo, hi = series(1 - 1 / x)
    lo2, hi2 = series(F(1, 2))
    if exponent >= 0:
        lo += exponent * lo2
        hi += exponent * hi2
    else:
        lo += exponent * hi2
        hi += exponent * lo2
    return (F(lo, scale), F(hi, scale))
S = 10 ** 9
T = 10 ** 18

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def certificate(root, rawdir):
    meta = json.loads((root / 'results/inventory_certificate.json').read_text())
    need(meta['status'] == 'PASS' and meta['core'] == 102 and meta['finite_moduli_are_bases_at_N1'] and meta['all_positive_prime_valuations'] and (meta['normalization'] == 'mu(B^N)/ln(B)'), 'explicit residual N1 domain; completeness is audited separately')
    case = dict(pmin=meta['prime_range'][0], pmax=meta['prime_range'][1], bmax=meta['base_upper_exclusive'], physical_states=meta['witness_states'], physical_records=meta['witness_representatives'])
    expected_moduli = meta['witness_moduli']
    rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader((root / 'results/finite_bounds.csv').open())]
    specified = [tuple(map(int, l.split())) for l in (root / 'results/finite_domain.txt').read_text().splitlines()]
    streams = json.loads((root / 'results/finite_stream_manifest.json').read_text())
    receipt = json.loads((root / 'results/finite_receipt.json').read_text())
    need([(r['M'], r['p']) for r in rows] == specified == [(r['M'], r['p']) for r in streams], 'every specified physical stream exactly once')
    need(len(rows) == expected_moduli and len(set(specified)) == expected_moduli, 'all unique physical modules in the frozen partition')
    table = list(map(int, (root / 'data/log_table.txt').read_text().split()))
    need(len(table) == 2049 and all((t == 0 for t in table[:1025])), 'primary logarithm table shape')

    @lru_cache(None)
    def logupper(n):
        k = n.bit_length() - 1
        shift = max(k - 10, 0)
        mantissa = n >> shift if k >= 10 else n << 10 - k
        rest = n - mantissa * (1 << shift) if k >= 10 else 0
        base = mantissa * (1 << shift) if k >= 10 else 1
        declared = 2 * table[mantissa] + 2 * k * table[2048] + ceil_fraction(F(rest * S, base))
        low, high = fast_log_interval(F(n))
        independent = ceil_fraction(high * T)
        need(F(independent, T) <= F(declared, S), 'separate rational series validates every actually used primary upper')
        return (declared, independent)
    results = []
    records = states = 0
    negative = 0
    for row, stream in zip(rows, streams):
        M, p = (row['M'], row['p'])
        need(prime(p) and case['pmin'] <= p < case['pmax'] and (prime_factors(M) == [2, 3, 17, p]), 'true complete support of every finite base')
        need(1938 <= M < case['bmax'] and M % 2 == 0, 'physical modulus and even reflection domain')
        path = rawdir / f'M{M}.bin'
        raw = path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        expected = M // 2 + 1
        need(len(raw) == 8 * expected == stream['bytes'] and expected == stream['records'] and (digest == stream['sha256']), 'raw length and exact stream hash')
        heights = Counter()
        mass = 0
        maxheight = 0
        for s, (a, d) in enumerate(struct.iter_unpack('<ii', raw)):
            height = max(abs(a), d)
            need(0 < d <= M // 2 and gcd(d, M) == 1 and ((a - s * d) % M == 0) and (1 <= height <= M // 2), 'every independent physical witness has positive full-unit denominator, correct state and height')
            weight = 1 if s == 0 or 2 * s == M else 2
            heights[height] += weight
            mass += weight
            maxheight = max(maxheight, height)
            records += 1
        need(mass == M and row['states'] == M and (row['representatives'] == expected), 'all states including zero and midpoint have their exact mass')
        total = alternate = 0
        for height, weight in heights.items():
            value, alt = logupper(height)
            total += weight * value
            alternate += weight * alt
        need(total == row['sum_log_upper_ticks'] and maxheight == row['max_height'], 'full declared weighted logarithm sum reconstructed')
        numerator = 2 * total - M * row['log_M_lower_ticks']
        need(numerator == row['mean_numerator'] and row['stencil'] in (0, 1), 'signed finite mean and declared refinement mode')
        need(row['stencil'] == 0 and row['first_numerator'] == numerator or (row['stencil'] == 1 and 85 * max(0, row['first_numerator']) >= 4 * M * row['log_M_lower_ticks']), 'first-pass equality or explicitly failed first gate')
        low, high = fast_log_interval(F(M))
        independent_loglow = floor_fraction(low * T)
        need(F(row['log_M_lower_ticks'], S) <= F(independent_loglow, T), 'separate lower log of the actual base M=B at N1')
        alt_numerator = 2 * alternate - M * independent_loglow
        primary_mean = F(numerator, 2 * M * S)
        alternate_mean = F(alt_numerator, 2 * M * T)
        need(alternate_mean <= primary_mean, 'independent full consumer improves the declared upper mean')
        normalized = F(max(0, numerator), 2 * M * row['log_M_lower_ticks'])
        need(normalized < F(2, 85), 'finite normalized gate with ln B=ln M only because N=1')
        negative += primary_mean < 0
        states += M
        results.append(dict(M=M, p=p, states=M, records=expected, sha256=digest, primary_mean_upper=str(primary_mean), independent_mean_upper=str(alternate_mean), normalized_upper=str(normalized), max_height=maxheight, distinct_heights=len(heights), height_histogram={str(h): count for h, count in sorted(heights.items())}))
    need(states == receipt['states'] == case['physical_states'] and records == receipt['representatives'] == case['physical_records'], 'all physical states and raw records consumed independently')
    need(receipt['status'] == 'PASS' and receipt['all_upper_means_negative'] == (negative == len(rows)) and receipt['inputs_unchanged'] and (receipt['stream_bytes'] == 8 * records) and (receipt['refined_moduli'] == sum((r['stencil'] for r in rows))), 'complete generation receipt arithmetic and scope')
    need(receipt['output_sha256'] == sha(root / 'results/finite_bounds.csv'), 'finite result table checksum')
    return dict(status='PASS', checks=CHECKS, moduli=len(rows), states=states, fully_consumed_records=records, raw_bytes=8 * records, negative_upper_means=negative, independent_full_raw_consumer=True, uses_lattice_reduction=False, assumes_minimality=False, log_series='independent directed integer Mercator128,scale10^40,explicit tail', declared_domain_completeness_audited_separately=True, consumer_has_no_producer_imports=True, independently_checked_height_logs=logupper.cache_info().currsize, finite_scope_N=1, all_N_claimed=False, results=results, inputs={name: sha(root / name) for name in ['code/finite.cpp', 'data/log_table.txt', 'data/log_seeds.tsv', 'results/finite_bounds.csv', 'results/finite_domain.txt', 'results/finite_stream_manifest.json', 'results/finite_receipt.json', 'results/inventory_certificate.json']})

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--raw-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    z = certificate(a.root.resolve(), a.raw_dir.resolve())
    a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status=z['status'], checks=z['checks'], records=z['fully_consumed_records'], states=z['states'], moduli=z['moduli']), indent=2))
if __name__ == '__main__':
    main()
