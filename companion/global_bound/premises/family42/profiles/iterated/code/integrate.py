from pathlib import Path
from fractions import Fraction as F
from exact import need, log_hi, log_lo, out_decimal
S = 10 ** 9
DEN = 2 ** 52

def profile(root, k=2, prefix=None):
    root = Path(root)
    p = 11
    q = 462
    nl = 1152
    prefix = Path(prefix) if prefix is not None else root / 'results/q462'
    c = [list(map(int, l.split())) for l in Path(str(prefix) + f'.cells{k}').read_text().splitlines()]
    r = [list(map(int, l.split())) for l in Path(str(prefix) + f'.radial{k}').read_text().splitlines()]
    need(len(c) == 2688 and len(r) == 128, 'all profile cells')
    for i, z in enumerate(c):
        need(len(z) == 4 and z[0] == i and (0 <= z[2] < nl) and (0 <= z[3] < 33), 'cell metadata')
    for j, z in enumerate(r):
        need(len(z) == 5 and z[0] == j - 125 and (0 < z[1] < z[2]) and (z[3] == max((x[1] for x in c[j * 21:(j + 1) * 21]))), 'radial maximum')
        need(z[4] in [x[0] for x in c[j * 21:(j + 1) * 21] if x[1] == z[3]], 'actual maximizing cell')
        if j:
            need(z[1] == r[j - 1][2], 'no radial gaps')
    need(3 * r[-1][2] ** 2 >= 4 * DEN ** 2, 'Gauss domain covered')
    v = [F(z[3], S) for z in r]
    for j in range(126, -1, -1):
        v[j] = max(v[j], v[j + 1])
    a = F(29, 14)
    rho = F(r[0][1], DEN)
    bulk = v[-1] + sum(((v[j] - v[j + 1]) * min(1, a * F(r[j][2], DEN)) for j in range(127)))
    t = log_hi(q - 1) + k * log_hi(q) / 2 - log_lo(rho) / 2 + 4 * rho
    cusp = a * rho * (max(t - v[0], F(0)) + F(1, 2))
    total = bulk + cusp
    return dict(p=p, q=q, k=k, labels=nl, directions=32, cells=2688, bulk_upper=str(bulk), cusp_upper=str(cusp), total_upper=str(total), upper_decimal=out_decimal(total), target_lower=str(F(2, 85) * log_lo(q)), radical_gate=total < F(2, 85) * log_lo(q), normalized_upper=str(max(F(0), total) / log_lo(q)), monotone_profile_ticks=[int(x * S) for x in v])
if __name__ == '__main__':
    import json, argparse
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', type=Path)
    ap.add_argument('--output', type=Path, default=root / 'results/integrated.json')
    args = ap.parse_args()
    rows = [profile(root, k, args.prefix) for k in (1, 2)]
    args.output.write_text(json.dumps(rows, indent=2) + '\n', newline='\n')
    for r in rows:
        print('ITERATED', r['k'], r['upper_decimal'], 'PASS', r['radical_gate'], flush=True)
