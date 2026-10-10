"""Independent HNF enumeration, rational geometry, and complete small state banks."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json
from checks import algebra_check, log_fp, FP, ps_of, direct, gauss
from integrate import rectangles, integrate
from exact import need, log_hi, log_lo
R = Path(__file__).resolve().parents[1]

def audit(p, m, style, ag, prefix):
    den, old, gr = rectangles(R, style, ag)
    pr = integrate(p, m, style, str(prefix) + '.cells', R, ag)
    q = pr['q']
    K = pr['transfer_index']
    br, algebra = algebra_check(p, 1, prefix, m)
    rows = [list(map(int, l.split())) for l in Path(str(prefix) + '.audit').read_text().splitlines()]
    need(len(rows) == 96, 'stratified audit length')
    logs = interiors = 0
    for row in rows:
        c, h, ci, cj, u, v, a = row[:7]
        caps = row[7:]
        need(len(caps) == 2 * a + 1 and a == gr[c][7], 'audit angular nodes')
        _, ha, hb, hd, image = br[h]
        need((u, v) == (ci * ha + cj * hb, cj * hd), 'independent coefficient transport')
        need(all(((ci * x + cj * y) % pp for pp, (x, y) in zip(ps_of(p), image))), 'all marked primes admitted')
        _, _, _, lo, hi, tl, th, _ = gr[c]
        for n, cap in zip(range(-a, a + 1), caps):
            maxsq = F(0)
            for rho in (F(lo, den), F(hi, den)):
                for t in (F(tl, 40), F(th, 40)):
                    xx = rho * (u + t * v)
                    aa = a * xx + n * v
                    bb = -n * xx + a * v
                    sq = max(aa * aa, bb * bb) / (F(K) * rho * (a * a + n * n))
                    _, upper = log_fp(sq.numerator, sq.denominator)
                    need(2 * cap * FP >= upper * 10 ** 9, 'independent outward corner logarithm')
                    maxsq = max(maxsq, sq)
                    logs += 1
            for wi, wj in ((1, 3), (2, 2), (3, 1)):
                rho = F(wi * lo + (4 - wi) * hi, 4 * den)
                t = F(wj * tl + (4 - wj) * th, 160)
                xx = rho * (u + t * v)
                aa = a * xx + n * v
                bb = -n * xx + a * v
                sq = max(aa * aa, bb * bb) / (F(K) * rho * (a * a + n * n))
                need(sq <= maxsq, 'rational interior versus rectangle vertices')
                interiors += 1
        log_fp.cache_clear()
    raw = [list(map(int, l.split())) for l in Path(str(prefix) + '.cells').read_text().splitlines()]
    states = tests = parents = lifts = 0
    physical = []
    for L in (1, 2, 3, 5, 7):
        M = K * L
        values, nt = direct(M, q)
        states += M
        tests += nt
        hh = {h: log_hi(h) for h in set(values)}
        logM = log_lo(M)
        mean = sum((hh[h] for h in values)) / M - logM / 2
        need(mean <= F(pr['total_upper']), 'complete definition bank versus integrated profile')
        for s in range(L):
            u, v = gauss(L, s)
            nu = sum((x * x for x in u))
            dot = sum((x * y for x, y in zip(u, v)))
            rho = F(nu, L)
            t = F(dot, nu)
            ix = [i for i, r in enumerate(gr) if F(r[3], den) <= rho <= F(r[4], den) and F(r[5], 40) <= t <= F(r[6], 40)]
            need(ix, 'physical parent inside certified domain')
            c = ix[0]
            child = [s + b * L for b in range(K)]
            need(len(set(child)) == K and all((j % L == s for j in child)), 'exact physical fibre')
            value = sum((hh[values[j]] for j in child)) / K - logM / 2
            need(value <= F(raw[c][8], 10 ** 9), 'conditional physical mean versus whole cell cap')
            parents += 1
            lifts += K
        physical.append(dict(L=L, M=M, mean_upper=str(mean), denominator_tests=nt))
    direct.cache_clear()
    return dict(status='PASS', q=q, transfer_index=K, algebra=algebra, stratified_geometry=dict(witnesses=96, independent_corner_logarithms=logs, rational_interior_checks=interiors), physical=dict(state_modulus_incidences=states, denominator_tests=tests, conditional_parent_means=parents, fibre_states=lifts, rows=physical), alternative_geometric_audit_exhaustive=False)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('multiplier', type=int)
    ap.add_argument('style')
    ap.add_argument('--angular', type=int, default=8)
    ap.add_argument('--prefix', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    r = audit(a.p, a.multiplier, a.style, a.angular, a.prefix)
    a.output.write_text(json.dumps(r, indent=2) + '\n', newline='\n')
    print(json.dumps(r, indent=2))
