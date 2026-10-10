"""Independent Python checks. Geometry is sampled; C++ consumption is exhaustive."""
from pathlib import Path
from fractions import Fraction as F
import json, sys, argparse
from assemble_refined import R, refined_grid, assemble
from checks import algebra_check, geometric_check, direct, gauss
from exact import need, log_hi, log_lo

def audit(prefix):
    den, old, gr, mapping = refined_grid()
    pr = assemble(prefix)
    br, alg = algebra_check(11, 1, prefix, 6)
    geo = geometric_check(11, 1, 16, prefix, gr, den, br, 6)
    raw = [list(map(int, x.split())) for x in Path(str(prefix) + '.cells').read_text().splitlines()]
    coarse = [list(map(int, x.split())) for x in (R / 'imports/q2310_K13860_AG4.cells').read_text().splitlines()]
    combined = old[:1890] + [[1890 + i, *r[1:]] for i, r in enumerate(gr)]
    caps = [r[3] for r in coarse[:1890]] + [r[7] for r in raw]
    rows = []
    parents = states = tests = fibres = 0
    for L in [1, 2, 3, 5, 7]:
        M = 13860 * L
        vals, nt = direct(M, 2310)
        hh = {h: log_hi(h) for h in set(vals)}
        lm = log_lo(M)
        mean = sum((hh[h] for h in vals)) / M - lm / 2
        need(mean <= F(pr['total_upper']), 'complete physical mean versus hybrid upper')
        for s in range(L):
            u, v = gauss(L, s)
            nu = sum((x * x for x in u))
            dot = sum((x * y for x, y in zip(u, v)))
            rho = F(nu, L)
            t = F(dot, nu)
            ix = [i for i, r in enumerate(combined) if F(r[1], den) <= rho <= F(r[2], den) and F(r[3], 40) <= t <= F(r[4], 40)]
            need(ix, 'physical parent in hybrid cover')
            children = [s + b * L for b in range(13860)]
            need(len(set(children)) == 13860 and all((j % L == s for j in children)), 'physical fibre')
            value = sum((hh[vals[j]] for j in children)) / 13860 - lm / 2
            need(value <= F(caps[ix[0]], 10 ** 9), 'conditional physical mean versus hybrid cell')
            parents += 1
            fibres += 13860
        rows.append(dict(L=L, M=M, mean_upper=str(mean), denominator_tests=nt))
        states += M
        tests += nt
    return dict(status='PASS', algebra=alg, stratified_geometry=geo, physical=dict(rows=rows, state_modulus_incidences=states, denominator_tests=tests, conditional_parent_means=parents, fibre_states=fibres), all_refined_records_consumed_by_this_python_audit=False, consumer_CPP_exhaustive_separately=True, geometric_second_audit_exhaustive=False, coarse_certificate_imported=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', type=Path, default=R / 'results/profile')
    ap.add_argument('--output', type=Path, default=R / 'results/profile.independent.json')
    a = ap.parse_args()
    z = audit(a.prefix)
    a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps({k: v for k, v in z.items() if k != 'physical'}, indent=2))
