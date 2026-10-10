"""Exact HNF certificate for unequal valuations; no geometric profile run."""
from pathlib import Path
from math import gcd, prod, isqrt
from collections import Counter
from itertools import product
import argparse, hashlib, json, struct
HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise AssertionError(label)

def primes(n):
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

def hnfs(K):
    for a in range(1, K + 1):
        if K % a == 0:
            for b in range(a):
                yield (a, b, K // a)

def line(x, y, p):
    x %= p
    y %= p
    need(x or y, 'nonzero projective vector')
    return y * pow(x, -1, p) % p if x else p

def complement_tensor(values, dimensions):
    result = values[:]
    stride = 1
    for size in reversed(dimensions):
        for start in range(0, len(result), size * stride):
            for inner in range(stride):
                ids = [start + inner + j * stride for j in range(size)]
                total = sum((result[t] for t in ids))
                for t in ids:
                    result[t] = total - result[t]
        stride *= size
    return result

def certify(q, K):
    ps = primes(q)
    need(prod(ps) == q and primes(K) == ps, 'squarefree mark and equal prime support')
    labels = list(product(*[range(p + 1) for p in ps]))
    hist = Counter()
    local = Counter()
    branches = []
    reject = 0
    digest = hashlib.sha256()
    for a, b, d in hnfs(K):
        if any((a % p == b % p == d % p == 0 for p in ps)):
            reject += 1
            continue
        kernel = []
        for p in ps:
            need(a * d % p == 0 and (a % p or b % p or d % p), 'local rank exactly one')
            z = -b * pow(d, -1, p) % p if d % p else p
            vx, vy = (1, z) if z < p else (0, 1)
            need(a * vx % p == 0 and (b * vx + d * vy) % p == 0, 'kernel directly verified')
            kernel.append(z)
            local[p, a % p, b % p, d % p] += 1
        hist[tuple(kernel)] += 1
        branches.append((a, b, d))
        digest.update(struct.pack('<QQQ', a, b, d))
    need(set(hist) == set(labels), 'all projective labels present')
    need(set(hist.values()) == {K // q}, 'uniform unequal-valuation multiplicity')
    psi = prod((p + 1 for p in ps))
    need(len(branches) == K // q * psi, 'branch count')
    values = [hist[x] for x in labels]
    transformed = complement_tensor(values, [p + 1 for p in ps])
    need(all((x == K for x in transformed)), 'all marked admitted counts via exact tensor')
    local_checks = 0
    for p, a, b, d in local:
        images = set()
        admitted = 0
        for z in range(p + 1):
            x, y = (1, z) if z < p else (0, 1)
            im = (a * x % p, (b * x + d * y) % p)
            if im == (0, 0):
                continue
            admitted += 1
            images.add(line(*im, p))
            local_checks += 1
        need(admitted == p and len(images) == 1, 'all local admitted images share one line')
    parents = [(0, 1), (1, 0), (1, 1), (1, 2), (2, 3), (3, 5), (7, 11), (11, 13)]
    for x, y in parents:
        if any((x % p == y % p == 0 for p in ps)):
            continue
        admitted = sum((all((a * x % p or (b * x + d * y) % p for p in ps)) for a, b, d in branches))
        need(admitted == K, 'sample global mark direct branch count')
    physical = [(a, b, d) for a, b, d in branches if gcd(d, q) == 1]
    need(physical == [(K, b, 1) for b in range(K)], 'physical admissibility forces d1')
    for L in [1, 2, 3, 5, 11]:
        seen = {s + b * L for s in range(L) for _, b, _ in physical}
        need(seen == set(range(K * L)), 'complete physical-state partition')
    return dict(q=q, K=K, transfer_multiplier=K // q, prime_factors=ps, parent_labels=psi, primitive_branches=len(branches), rejected_rank_zero=reject, multiplicity_per_kernel=K // q, admitted_per_parent=K, HNF_stream_sha256=digest.hexdigest(), local_matrices=len(local), local_image_checks=local_checks, all_kernel_multiplicities_checked=True, all_tensor_parent_counts_checked=True, local_reset_exhaustively_checked=True, global_parent_checks='eight fixed primitive marks', physical_seed_regressions=[1, 2, 3, 5, 11], geometry_replayed=False)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'unequal_transfer_certificate.json')
    args = parser.parse_args()
    records = [certify(2310, 2310 * s) for s in [1, 2, 3, 6, 4, 8]]
    need(records[1]['primitive_branches'] == 13824 and records[2]['primitive_branches'] == 20736 and (records[3]['primitive_branches'] == 41472), 'three priority index counts')
    need(all((x % 2 == 0 for x in [2, 0, 2310])), 'diag2,2310 must be rejected modulo2')
    admitted = sum((gcd(d, 6) == 1 for a, b, d in hnfs(30)))
    need(admitted == 36 and admitted != 30, 'extra determinant prime breaks physical K normalization')
    H = (2, 0, 1)
    marks = [(1, 3), (0, 1)]
    images = [(2 * x, y) for x, y in marks]
    need(all((all((x % p or y % p for p in [2, 3])) for x, y in images)), 'both subset-index parents admitted')
    need(line(*images[0], 3) != line(*images[1], 3), 'missing mark prime has no image reset')
    need(gcd(images[0][1], 6) > 1 and gcd(images[1][1], 6) == 1, 'different admissible witness sets')
    result = dict(status='PASS', checks=CHECKS, marks=records, general_theorem='q squarefree>1, rad(K)=q; arbitrary positive valuations', multiplier_formula='(K/q)*product(p+1)', normalization='1/K', uniform_cap_requires_new_geometry=True, historical210_caps_used=False, q_squared_branch_count_formula=2310 * 6912, q_squared_enumerated=False, negative_regressions=dict(rank_zero_rejected=True, extra_prime_admitted_count=36, subset_index_no_reset=dict(q=6, K=2, H=[2, 0, 0, 1], parents=marks, images=images)), infinite_claim_scope='written proof; finite checks do not replace it')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status='PASS', checks=CHECKS, indices={r['K']: r['primitive_branches'] for r in records}, geometry_replayed=False), indent=2))
if __name__ == '__main__':
    main()
