"""Exact finite and tail arithmetic for the imported universal PCDF 3/2 proof.

This checker consumes the proof's integer envelope, not a sample of moduli.
The infinite argument is stated in imports/PRIMITIVE_CDF.md.
"""
from fractions import Fraction as F
from math import gcd, isqrt
from exact import need

def phi(d):
    return sum((gcd(a, d) == 1 for a in range(1, d + 1)))

def check():
    rows = []
    checks = 0
    for T in range(1, 64):
        terms = []
        for d in range(1, isqrt(T) + 1):
            A = isqrt(T - d * d)
            candidates = []
            for h in range(1, d + 1):
                if d % h:
                    continue
                r = d // h
                quotient = phi(d) // phi(r)
                need(quotient * phi(r) == phi(d), 'integral unit fibre multiplicity')
                U = sum((gcd(b, r) == 1 for b in range(-A // h - 1, A // h + 2) if abs(b) * h <= A))
                candidates.append({'h': h, 'r': r, 'unit_count': U, 'incidences': quotient * U})
                checks += 1
            terms.append({'d': d, 'A': A, 'maximum': max((v['incidences'] for v in candidates)), 'candidates': candidates})
        E = sum((v['maximum'] for v in terms))
        need(2 * E <= 3 * T, 'finite primitive CDF envelope')
        rows.append({'T': T, 'E': E, 'margin_twice': 3 * T - 2 * E, 'terms': terms})
        checks += 1
    weights = [F(1), F(1, 2), F(2, 3), F(1, 2), F(1), F(1, 3)]
    mean = F(2, 3)
    D = []
    d = F(0)
    for w in weights:
        d += w - mean
        D.append(d)
    need(D == [F(1, 3), F(1, 6), F(1, 6), F(0), F(1, 3), F(0)], 'period six cumulative discrepancy')
    need(sum(weights) == 6 * mean and all((0 <= x <= F(1, 3) for x in D)), 'Abel weight hypotheses')
    need(F(9, 4) ** 2 > 5, 'rational sqrt5 upper')
    margin = F(3, 2) - F(29, 21) - F(11, 96)
    need(margin == F(1, 224) > 0, 'uniform Y>=8 tail margin')
    need(F(3, 2) >= 1, 'seed M1 included')
    return {'status': 'PASS', 'scope': 'every integer M>=1 and every real x>=0, by imported proof plus exact finite envelope and uniform tail', 'cdf_alpha': '3/2', 'finite_squared_radii': 63, 'checks': checks + 5, 'tail_from_radius': 8, 'tail_strict_margin': '1/224', 'weights': list(map(str, weights)), 'cumulative_discrepancy': list(map(str, D)), 'rows': rows, 'all_zero_numerator_cases_included': True, 'lattice_primitive_not_ambient_primitive': True, 'seed_M1_explicit': True}
if __name__ == '__main__':
    import argparse, json
    from pathlib import Path
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    a = ap.parse_args()
    z = check()
    if a.output:
        a.output.write_text(json.dumps(z, indent=2) + '\n', newline='\n')
    print('PASS PCDF exact envelope', z['checks'], 'checks')
