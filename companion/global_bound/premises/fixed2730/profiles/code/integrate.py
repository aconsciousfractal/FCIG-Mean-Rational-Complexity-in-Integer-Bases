"""Exact integration of an autonomous coarse or hybrid marked profile."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import json, sys
from assemble import grid
from core import mark
from exact import need, log_hi, log_lo, out_decimal
S = 10 ** 9
CROP = 5200308914369309

def rectangles(root, style, ag=8):
    den, old = grid(root)
    need(style in ('coarse', 'hybrid'), 'grid style')
    need(3 * (CROP - 1) ** 2 < 4 * den ** 2 <= 3 * CROP ** 2, 'minimal rational Gauss endpoint')
    need(old[126 * 21][1] < CROP < old[126 * 21][2] and old[127 * 21][1] > CROP, 'Gauss restriction')
    rows = []
    for band in range(127):
        np = 4 if style == 'hybrid' and band >= 90 else 1
        for part in range(np):
            for t in range(21):
                c = band * 21 + t
                _, lo, hi, tl, th, *_ = old[c]
                width = min(hi, CROP) - lo
                nlo, nhi = (lo + width * part // np, lo + width * (part + 1) // np)
                need(old[c][1] <= nlo < nhi <= old[c][2], 'subrectangle inclusion')
                rows.append([len(rows), c, part, nlo, nhi, tl, th, ag if np == 4 else 4])
    for i, r in enumerate(rows):
        if i % 21 == 0 and i:
            need(r[3] == rows[i - 21][4], 'complete radial continuity')
    need(rows[0][3] == old[0][1] and rows[-1][4] == CROP, 'full continuous coverage')
    return (den, old, rows)

def integrate(p, multiplier, style, path, root, ag=8):
    need((multiplier, style) in ((1, 'coarse'), (6, 'hybrid')), 'declared index and grid')
    q, nl = mark(p)
    K = q * multiplier
    den, old, gr = rectangles(root, style, ag)
    raw = [list(map(int, l.split())) for l in Path(path).read_text().splitlines()]
    need(len(raw) == len(gr), 'all continuous output rectangles')
    for row, g in zip(raw, gr):
        need(len(row) == 15 and row[:8] == g, 'exact profile geometry')
        _, c, _, _, _, _, _, a, *_ = row
        embed = style == 'hybrid' and c >= 90 * 21
        need(row[8] == min(row[9], row[10]), 'minimum of complete continuous caps')
        need(all((0 <= row[i] < nl for i in (11, 13))) and all((-a <= row[i] <= a for i in (12, 14))), 'angular and mark labels')
        need(embed and row[14] % 2 == 0 or (not embed and row[8] == row[9] == row[10] and (row[11:13] == row[13:15])), 'embedded cap domain')
    pieces = [(F(gr[j][4], den), F(max((r[8] for r in raw[j:j + 21])), S)) for j in range(0, len(gr), 21)]
    v = [x[1] for x in pieces]
    for i in range(len(v) - 2, -1, -1):
        v[i] = max(v[i], v[i + 1])
    alpha = F(29, 14)
    rho = F(old[0][1], den)
    bulk = v[-1] + sum(((v[i] - v[i + 1]) * min(1, alpha * pieces[i][0]) for i in range(len(v) - 1)))
    cusp_at_rho = log_hi(q - 1) + log_hi(K) / 2 - log_lo(rho) / 2 + F(5, 2) * rho
    cusp = alpha * rho * (max(cusp_at_rho - v[0], 0) + F(1, 2))
    total = bulk + cusp
    nb = nl * multiplier
    fine_sources = 777 if style == 'hybrid' else 0
    coarse_sources = 2667 - fine_sources
    return dict(status='PASS', p=p, q=q, transfer_multiplier=multiplier, transfer_index=K, required_divisor_of_M=K, labels=nl, branches=nb, style=style, fine_angular_parameter=ag if fine_sources else None, source_cells=2667, output_cells=len(gr), radial_bands=len(pieces), gauss_upper=str(F(CROP, den)), geometric_witness_records=2667 * nb, witness_bytes=2667 * nb * 8, geometric_corner_evaluations=nb * (coarse_sources * 9 * 4 + fine_sources * 4 * (2 * ag + 1) * 4), embedded_grid_win_cells=sum((row[10] < row[9] for row in raw)), bulk_upper=str(bulk), cusp_upper=str(cusp), total_upper=str(total), upper_decimal=out_decimal(total), normalized_upper=str(max(total, F(0)) / log_lo(q)), normalized_upper_decimal=out_decimal(max(total, F(0)) / log_lo(q)), strict_2_over_85_gate=total < F(2, 85) * log_lo(q), all_positive_multiples_including_extra_primes=True, all_residues_including_zero=True, normalization='log B', monotone_profile_ticks=[int(x * S) for x in v])
if __name__ == '__main__':
    root = Path(__file__).resolve().parents[1]
    p, m, ag = map(int, sys.argv[1:4])
    style = sys.argv[4]
    pre = Path(sys.argv[5])
    out = Path(sys.argv[6])
    r = integrate(p, m, style, str(pre) + '.cells', root, ag)
    out.write_text(json.dumps(r, indent=2) + '\n', newline='\n')
    print(json.dumps({k: v for k, v in r.items() if k != 'monotone_profile_ticks'}, indent=2))
