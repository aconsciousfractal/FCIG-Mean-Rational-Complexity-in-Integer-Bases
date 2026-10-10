"""Independent algebra, logarithm and geometric audits, and complete small banks.
The finite-state consumer is exhaustive; the second geometric audit is stratified.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd, prod, isqrt
from collections import Counter
from functools import lru_cache
from exact import need, ceildiv, alternate_log, log_hi, log_lo, lower_tick, upper_tick, sieve
FP = 10 ** 40

def reduced_fp(a, b, terms=96):
    need(b <= a <= 2 * b, 'secondary binary normalization')
    n = a - b
    lo = n * FP // a
    hi = ceildiv(n * FP, a)
    s = t = 0
    for j in range(1, terms + 1):
        s += lo // j
        t += ceildiv(hi, j)
        lo = lo * n // a
        hi = ceildiv(hi * n, a)
    return (s, t + ceildiv(hi * a, (terms + 1) * b))
L2, U2 = reduced_fp(2, 1)

@lru_cache(None)
def log_fp(n, d=1):
    need(n > 0 and d > 0, 'positive secondary logarithm')
    k = n.bit_length() - d.bit_length()
    a, b = (n, d << k) if k >= 0 else (n << -k, d)
    if a < b:
        k -= 1
        a, b = (n, d << k) if k >= 0 else (n << -k, d)
    lo, hi = reduced_fp(a, b)
    return (lo + k * L2, hi + k * U2) if k >= 0 else (lo + k * U2, hi + k * L2)

def input_checks(root):
    tab = list(map(int, (root / 'data/log_table.txt').read_text().split()))
    need(len(tab) == 2049 and all((z == 0 for z in tab[:1025])), 'mantissa table layout')
    for i in range(1025, 2049):
        lo, hi = log_fp(i, 1024)
        need(hi * 10 ** 9 <= 2 * tab[i] * FP and lo * 10 ** 9 >= 2 * (tab[i] - 1) * FP, 'independent mantissa enclosure')
    seeds = [list(map(int, l.split())) for l in (root / 'data/log_seeds.tsv').read_text().splitlines()]
    need(len(seeds) == 1024, 'all logarithm seeds')
    for i, (n, lo, hi) in enumerate(seeds, 1):
        a, b = alternate_log(n)
        need(n == i and lo == lower_tick(a, 10 ** 17) and (hi == upper_tick(b, 10 ** 17)), 'secondary seeds regenerated')
    for n, d in [(j + 100, j + 1) for j in range(120)] + [(1, j) for j in range(1, 32)]:
        a, b = alternate_log(F(n, d))
        lo, hi = log_fp(n, d)
        need(F(lo, FP) <= a <= b <= F(hi, FP), 'fixed-point secondary against exact rational series')
    return dict(status='PASS', mantissa_entries=1024, dense_log_seed_entries=1024, secondary_series_rational_checks=151)
if __name__ == '__main__':
    import json
    root = Path(__file__).resolve().parent.parent
    result = input_checks(root)
    (root / 'results/log_input_checks.json').write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(result))
