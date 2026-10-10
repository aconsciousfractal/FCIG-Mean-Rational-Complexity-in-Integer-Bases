"""Same consumed cells with an explicit alternative universal CDF coefficient.

This arithmetic certificate is conditional until the named CDF lemma is proved
and bound to the main deliverable. It never changes the geometric contract.
"""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import argparse, json
from integrate import integrate, rectangles
from exact import need, log_hi, log_lo, out_decimal
R = Path(__file__).resolve().parents[1]

def assemble_alpha(p, m, style, ag, prefix, alpha):
    need(0 < alpha <= F(29, 14), 'positive sharpened CDF coefficient')
    base = integrate(p, m, style, str(prefix) + '.cells', R, ag)
    den, old, gr = rectangles(R, style, ag)
    v = [F(x, 10 ** 9) for x in base['monotone_profile_ticks']]
    ends = [F(gr[j][4], den) for j in range(0, len(gr), 21)]
    rho = F(old[0][1], den)
    bulk = v[-1] + sum(((v[i] - v[i + 1]) * min(1, alpha * ends[i]) for i in range(len(v) - 1)))
    q, K = (base['q'], base['transfer_index'])
    c = log_hi(q - 1) + log_hi(K) / 2 - log_lo(rho) / 2 + F(5, 2) * rho
    cusp = alpha * rho * (max(c - v[0], 0) + F(1, 2))
    total = bulk + cusp
    need(total <= F(base['total_upper']), 'sharper CDF cannot worsen a nonincreasing profile')
    return dict(status='ARITHMETIC_PASS_CONDITIONAL_ON_CDF_LEMMA', q=q, transfer_index=K, alpha=str(alpha), required_CDF_lemma=f'For every L, Pr(rho<=x)<=min(1,{alpha}*x) for all positive x', source_cells_sha256=sha256(Path(str(prefix) + '.cells').read_bytes()).hexdigest(), base_total_upper=base['total_upper'], bulk_upper=str(bulk), cusp_upper=str(cusp), total_upper=str(total), upper_decimal=out_decimal(total), normalized_upper=str(max(total, F(0)) / log_lo(q)), normalized_upper_decimal=out_decimal(max(total, F(0)) / log_lo(q)), strict_2_over_85_gate=total < F(2, 85) * log_lo(q), unchanged_geometric_contract=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('p', type=int)
    ap.add_argument('multiplier', type=int)
    ap.add_argument('style')
    ap.add_argument('--angular', type=int, default=8)
    ap.add_argument('--alpha', default='3/2')
    ap.add_argument('--prefix', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    z = assemble_alpha(a.p, a.multiplier, a.style, a.angular, a.prefix, F(a.alpha))
    a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print(json.dumps(z, indent=2))
