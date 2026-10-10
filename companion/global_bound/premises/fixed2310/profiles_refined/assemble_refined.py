"""Exact hybrid integration: certified coarse bands0..89 plus all refined90..126."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json, sys
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
from exact import need, log_hi, log_lo, out_decimal
from assemble import grid, profile
S = 10 ** 9

def refined_grid():
    den, old = grid(R)
    con = json.loads((R / 'CONTRACT.json').read_text())
    hi = con['gauss_crop_hi']
    need(con['mark'] == 2310 and con['transfer_index'] == 13860 and (con['source_bands'] == [90, 126]) and (con['radial_parts'] == 4) and (con['angular_grid'] == 16), 'fixed refinement contract')
    need(3 * (hi - 1) ** 2 < 4 * den ** 2 <= 3 * hi ** 2, 'minimal outward Gauss endpoint')
    need(old[126 * 21][1] < hi < old[126 * 21][2] and old[127 * 21][1] > hi, 'tail restriction')
    new = []
    mapping = []
    for b in range(90, 127):
        for part in range(4):
            for t in range(21):
                c = b * 21 + t
                _, lo, top, tl, th, a, v, d = old[c]
                top = min(hi, top)
                width = top - lo
                nlo = lo + width * part // 4
                nhi = lo + width * (part + 1) // 4
                need(0 < nlo < nhi <= old[c][2] and old[c][1] <= nlo, 'subrectangle inclusion')
                new.append([len(new), nlo, nhi, tl, th, a, v, d])
                mapping.append((c, part))
    for i, row in enumerate(new):
        if i % 21 == 0 and i:
            need(row[1] == new[i - 21][2], 'continuous refined radial coverage')
    need(new[0][1] == old[89 * 21][2] and new[-1][2] == hi, 'hybrid joins and coverage')
    return (den, old, new, mapping)

def assemble(prefix=None):
    prefix = Path(prefix) if prefix else R / 'results/profile'
    den, old, gr, mapping = refined_grid()
    oldpath = R / 'imports/q2310_K13860_AG4.cells'
    oldcaps = [list(map(int, l.split())) for l in oldpath.read_text().splitlines()]
    recomputed = profile(11, oldpath, R, 1, 4, 6)
    record = json.loads((R / 'imports/q2310_K13860_AG4.json').read_text())
    need(all((record[k] == v for k, v in recomputed.items())), 'imported coarse integration')
    raw = [list(map(int, l.split())) for l in Path(str(prefix) + '.cells').read_text().splitlines()]
    need(len(raw) == 3108, 'all new cells')
    for i, row in enumerate(raw):
        c, p = mapping[i]
        _, lo, hi, tl, th, *_ = gr[i]
        need(len(row) == 14 and row[:7] == [i, c, p, lo, hi, tl, th], 'subcell output exact geometry')
        need(row[7] == min(row[8], row[9]) and 0 <= row[10] < 6912 and (-16 <= row[11] <= 16) and (0 <= row[12] < 6912) and (-16 <= row[13] <= 16) and (row[13] % 2 == 0), 'complete-cell choice and angular labels')
    pieces = []
    for b in range(90):
        pieces.append((F(old[b * 21][2], den), F(max((r[3] for r in oldcaps[21 * b:21 * (b + 1)])), S)))
    for b in range(148):
        pieces.append((F(gr[b * 21][2], den), F(max((r[7] for r in raw[21 * b:21 * (b + 1)])), S)))
    v = [p[1] for p in pieces]
    for i in range(len(v) - 2, -1, -1):
        v[i] = max(v[i], v[i + 1])
    alpha = F(29, 14)
    rho = F(old[0][1], den)
    bulk = v[-1] + sum(((v[i] - v[i + 1]) * min(1, alpha * pieces[i][0]) for i in range(len(v) - 1)))
    cusp_at_rho = log_hi(2309) + log_hi(13860) / 2 - log_lo(rho) / 2 + F(5, 2) * rho
    cusp = alpha * rho * (max(cusp_at_rho - v[0], 0) + F(1, 2))
    total = bulk + cusp
    return dict(status='PASS', q=2310, transfer_index=13860, transfer_multiplier=6, valid_for_every_positive_multiple_including_extra_primes=True, labels=6912, branches=41472, coarse_cells_retained=1890, new_source_cells=777, new_subcells=3108, hybrid_radial_bands=238, hybrid_cells=4998, angular_grid=16, embedded_grid=8, embedded_grid_win_cells=sum((r[9] < r[8] for r in raw)), coarse_component_dependency='baseline marked2310 K13860 AG4', coarse_cells_sha256=sha256(oldpath.read_bytes()).hexdigest(), bulk_upper=str(bulk), cusp_upper=str(cusp), total_upper=str(total), upper_decimal=out_decimal(total), target_lower=str(F(2, 85) * log_lo(2310)), direct_radical_gate=total < F(2, 85) * log_lo(2310), normalized_upper=str(max(total, F(0)) / log_lo(2310)), required_divisor_of_M=13860, monotone_profile_ticks=[int(x * S) for x in v], rho_upper=str(F(gr[-1][2], den)), geometric_witness_records=777 * 41472, geometric_corner_evaluations=777 * 41472 * 4 * 33 * 4)
if __name__ == '__main__':
    prefix = Path(sys.argv[1]) if len(sys.argv) > 1 else R / 'results/profile'
    z = assemble(prefix)
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else R / 'results/profile.json'
    out.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps({k: v for k, v in z.items() if k != 'monotone_profile_ticks'}, indent=2))
