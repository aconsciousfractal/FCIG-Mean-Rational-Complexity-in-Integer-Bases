"""Independent two-chart P^1(Z/KZ) enumeration for the new unequal indices.
No imports from the producer/HNF implementation. Standard library only.
"""
from pathlib import Path
from itertools import product
from math import gcd, prod
from collections import Counter
import argparse, json, hashlib
PS = (2, 5, 13)
Q = 130
MULTIPLIERS = (1, 2, 5, 13)
CHECKS = 0

def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ArithmeticError(message)

def crt(rs, ms):
    k = prod(ms)
    return sum((r * (k // m) * pow(k // m, -1, m) for r, m in zip(rs, ms))) % k

def label(x, y):
    v = 0
    for p in PS:
        need(x % p or y % p, 'primitive local column')
        j = y * pow(x % p, -1, p) % p if x % p else p
        v = v * (p + 1) + j
    return v

def enumerate_index(m):
    K = Q * m
    mods = [p * p if m % p == 0 else p for p in PS]
    charts = []
    for p, mod in zip(PS, mods):
        local = [(1, t) for t in range(mod)] + [(p * t, 1) for t in range(mod // p)]
        need(len(local) == (p + 1) * (mod // p), 'two-chart count')
        counts = Counter()
        for x, y in local:
            need(gcd(gcd(x, y), mod) == 1, 'surjective local form')
            counts[y * pow(x, -1, p) % p if x % p else p] += 1
        need(counts == Counter({j: mod // p for j in range(p + 1)}), 'uniform local lifts')
        charts.append(local)
    lifted = {}
    multiplicities = Counter()
    for row in product(*charts):
        x = crt([r[0] for r in row], mods)
        y = crt([r[1] for r in row], mods)
        need(gcd(gcd(x, y), K) == 1, 'global primitive form')
        g = gcd(x, K)
        a = K // g
        d = g
        b = -y * pow(x // g, -1, a) % a if a > 1 else 0
        need(a * x % K == 0 and (b * x + d * y) % K == 0, 'annihilated rows')
        need(a * d == K and 0 <= b < a, 'canonical HNF index')
        h = (a, b, d)
        lab = label(x, y)
        need(h not in lifted, 'distinct local classes have distinct HNF')
        lifted[h] = lab
        multiplicities[lab] += 1
    reference = {}
    for a in range(1, K + 1):
        if K % a:
            continue
        d = K // a
        for b in range(a):
            if gcd(gcd(a, b), d) != 1:
                continue
            labs = []
            for p in PS:
                zero = [(x, y) for x, y in [(1, t) for t in range(p)] + [(0, 1)] if a * x % p == 0 and (b * x + d * y) % p == 0]
                need(len(zero) == 1, 'single rank-one kernel')
                x, y = zero[0]
                labs.append(y if x else p)
            lab = 0
            for p, j in zip(PS, labs):
                lab = lab * (p + 1) + j
            reference[a, b, d] = lab
    need(lifted == reference, 'two independent complete HNF constructions')
    nl = prod((p + 1 for p in PS))
    need(len(lifted) == m * nl, 'total primitive sublattices')
    need(multiplicities == Counter({j: m for j in range(nl)}), 'every label has multiplier lifts')

    def decode(j):
        z = []
        for p in reversed(PS):
            z.append(j % (p + 1))
            j //= p + 1
        return tuple(reversed(z))
    dec = [decode(j) for j in range(nl)]
    vals1 = [(j * j + 17 * j + 13) % 101 - 50 for j in range(nl)]
    vals2 = [j % 7 - 3 for j in range(nl)]
    tensor = []
    for i in range(nl):
        admitted = [j for j in range(nl) if all((x != y for x, y in zip(dec[i], dec[j])))]
        need(len(admitted) * m == K, 'all parents have K admitted branches')
        tensor.append([i, K, m * sum((vals1[j] for j in admitted)), m * sum((vals2[j] for j in admitted))])
    return ({'K': K, 'multiplier': m, 'local_prime_powers': mods, 'branches': len(lifted), 'labels': nl, 'lifts_per_label': m, 'all_parent_tensor_rows': nl, 'hnf_sha256': hashlib.sha256(json.dumps(sorted(((lab, *h) for h, lab in lifted.items()))).encode()).hexdigest()}, lifted, tensor)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    rows = []
    for m in MULTIPLIERS:
        row, hs, tensor = enumerate_index(m)
        if a.results:
            p = a.results / f'q130_K{Q * m}'
            stored = [list(map(int, l.split())) for l in Path(str(p) + '.branches').read_text().splitlines()]
            need([r[:4] for r in stored] == [list((lab, *h)) for h, lab in sorted(hs.items(), key=lambda kv: (kv[1], kv[0]))], 'all stored branches compared to local-chart enumeration')
            got = [list(map(int, l.split())) for l in Path(str(p) + '.tensor').read_text().splitlines()]
            need(got == tensor, 'every signed tensor row checked directly')
        rows.append(row)
    out = {'status': 'PASS', 'core': Q, 'core_primes': list(PS), 'checks': CHECKS, 'stored_profile_algebra_compared': bool(a.results), 'enumerations': rows, 'scope': 'proof plus independent complete local quotient and HNF construction, no numerical profile assumed'}
    a.output.write_text(json.dumps(out, indent=2) + '\n', newline='\n')
    print('PASS mixed indices', CHECKS, 'checks')
if __name__ == '__main__':
    main()
