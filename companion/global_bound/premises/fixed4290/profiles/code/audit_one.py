"""Second Python arithmetic/algebra audit; geometric coverage is stratified.
The C++ consumer is exhaustive. Complete definition-based physical banks are separate
cross-checks and do not replace the proof of continuous coverage.
"""
from pathlib import Path
from fractions import Fraction as F
import sys, json, argparse
from checks import input_checks, algebra_check, geometric_check, direct, gauss
from assemble import grid, profile
from exact import need, log_hi, log_lo
from core import mark
ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('ag', type=int)
    ap.add_argument('--k', type=int, default=1)
    ap.add_argument('--prefix', type=Path)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    p, ag, k = (args.p, args.ag, args.k)
    q, _ = mark(p)
    K = q * k
    prefix = args.prefix or ROOT / f'results/q{q}'
    den, gr = grid(ROOT)
    pr = profile(p, Path(str(prefix) + '.cells'), ROOT, k, ag)
    alldata, algebra = algebra_check(p, k, prefix)
    geo = geometric_check(p, k, ag, prefix, gr, den, alldata)
    td = 40
    if (ROOT / 'data/grid_meta.json').exists():
        td = json.loads((ROOT / 'data/grid_meta.json').read_text())['shear_denominator']
    cases = (1, 2, 3, 5, 7)
    physical = []
    states = tests = parents = lift = 0
    caps = [list(map(int, l.split()))[3] for l in Path(str(prefix) + '.cells').read_text().splitlines()]
    for L in cases:
        M = K * L
        R, nt = direct(M, q)
        states += M
        tests += nt
        hh = {h: log_hi(h) for h in set(R)}
        logM = log_lo(M)
        mean = sum((hh[h] for h in R)) / M - logM / 2
        need(mean <= F(pr['total_upper']), 'complete definition bank versus integrated upper')
        for s in range(L):
            u, v = gauss(L, s)
            nu = sum((x * x for x in u))
            dot = sum((x * y for x, y in zip(u, v)))
            rho = F(nu, L)
            t = F(dot, nu)
            indexes = [i for i, r in enumerate(gr) if F(r[1], den) <= rho <= F(r[2], den) and F(r[3], td) <= t <= F(r[4], td)]
            need(indexes, 'physical geometry inside grid')
            c = indexes[0]
            child = [s + b * L for b in range(K)]
            need(len(set(child)) == K and all((j % L == s for j in child)), 'physical fibre cardinality and congruence')
            lift += K
            value = sum((hh[R[j]] for j in child)) / K - logM / 2
            need(value <= F(caps[c], 10 ** 9), 'complete conditional physical mean versus cell upper')
            parents += 1
        physical.append(dict(L=L, M=M, mean_upper=str(mean), denominator_tests=nt))
    direct.cache_clear()
    result = {'status': 'PASS', 'p': p, 'q': q, 'algebra': algebra, 'stratified_geometry': geo, 'physical': {'state_modulus_incidences': states, 'denominator_tests': tests, 'conditional_parent_means': parents, 'fibre_states': lift, 'rows': physical}, 'geometric_second_audit_exhaustive': False}
    out = args.output or Path(str(prefix) + '.independent.json')
    out.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print('PASS independent', p, 'algebra branches', algebra['branches'], 'physical states', states, flush=True)
if __name__ == '__main__':
    main()
