"""Rational gate checker with fresh marked-profile inputs.

Negative caps use max(0,cap+cost) before division, preserving monotonicity.
This verifies the profile contract/bytes/numerical use, not its geometry.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
from hashlib import sha256
import json
from exact105 import require, record, fraction, log_interval, phi_interval, decimal_bound
HERE = Path(__file__).resolve().parent
DELTA = F(2, 85)
C105 = F(90985, 3948)
C1155 = F(2472433, 51744)

def prime(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

def lower_sqrt(M):
    z = isqrt(M)
    require(z * z <= M < (z + 1) ** 2, 'directed integer square root')
    return z

def periodic105(p, M):
    require(F(24, 7 * p) <= F(4, 3), 'periodic105 integration domain')
    return phi_interval(F(12, 7 * p))[1] + F(9, 14 * p) + C105 / lower_sqrt(M)

def free105(p):
    require(F(9, p) <= F(4, 3), 'M-free105 integration domain')
    return phi_interval(F(9, 2 * p))[1] + F(27, 16 * p)

def normalized105(p, B, M, cap):
    require(B > 1 and M >= B and (p > 0), 'normalization domain')
    costs = dict(periodic=periodic105(p, M), M_free=free105(p))
    loglower = log_interval(F(B))[0]
    bounds = {name: max(F(0), cap + cost) / loglower for name, cost in costs.items()}
    choice = min(bounds, key=bounds.get)
    return (costs, bounds, choice, bounds[choice])

def anchor(name, q, k, cap):
    path = HERE / 'imports' / name
    raw = path.read_bytes()
    j = json.loads(raw)
    provenance = json.loads((HERE / 'PROVENANCE.json').read_text())
    bridges = [v for v in provenance['sources'] if v['local_copy'] == 'imports/' + name]
    require(len(bridges) == 1, 'unique source bridge for marked anchor')
    require(bridges[0]['sha256'] == sha256(raw).hexdigest() and bridges[0]['bytes'] == len(raw), 'anchor byte identity with recorded fresh source')
    require((j['q'], j['k'], j['transfer_index']) == (q, k, q ** k), 'anchor mark/level/divisor')
    total = F(j['total_upper'])
    require(total < cap, 'fresh anchor strictly below simple cap')
    return dict(file='imports/' + name, q=q, k=k, divisor=q ** k, total_upper=record(total), simple_cap=record(cap), sha256=sha256(raw).hexdigest(), bytes=len(raw), valid_for_every_multiple_including_extra_primes=True, geometric_replay_here=False)

def main():
    theory = json.loads((HERE / 'theory105_certificate.json').read_text())
    require(theory['status'] == 'PASS', 'theory status')
    for q, R, A, C in [(105, 3, F(8, 7), F(1550, 329)), (1155, 5, F(80, 77), F(582, 77))]:
        row = theory['marks'][str(q)]
        require(row['radius'] == R and fraction(row['A']) == A and (fraction(row['C']) == C), 'exact theory contract')
    require(fraction(theory['penalties']['105']['sharp_sqrtM_coefficient']) == C105, 'sharp105 bridge')
    require(fraction(theory['penalties']['1155']['sharp_sqrtM_coefficient']) == C1155, 'sharp1155 bridge')
    anchors = [anchor('q105_k1.json', 105, 1, F(6, 125)), anchor('q105_k2.json', 105, 2, -F(2, 25)), anchor('q1155.json', 1155, 1, F(6, 125))]
    gates = []
    for name, p, B, M, cap, expected in [('large_p', 125, 13125, 13125, F(6, 125), 'M_free'), ('N_ge_2', 31, 3255, 3255 ** 2, -F(2, 25), 'periodic'), ('large_base', 31, 1000000, 1000000, F(6, 125), 'periodic')]:
        costs, choices, selected, upper = normalized105(p, B, M, cap)
        require(selected == expected and upper < DELTA, 'strict fixed105 gate')
        require(cap + costs[selected] > 0, 'positive corner numerator despite signed cap')
        gates.append(dict(name=name, p_scalar=p, B_min=B, M_min=M, sqrt_lower=lower_sqrt(M), cap=record(cap), costs={k: record(v) for k, v in costs.items()}, choices={k: record(v) for k, v in choices.items()}, selected=selected, upper=record(upper), upper_decimal_ceiling=decimal_bound(upper), margin=record(DELTA - upper), margin_decimal_floor=decimal_bound(DELTA - upper, False), normalization_uses_positive_part=True))
    small = [p for p in range(11, 31) if prime(p)]
    finite = [p for p in range(31, 125) if prime(p)]
    require(small == [11, 13, 17, 19, 23, 29] and len(finite) == 20, 'prime partition105')
    require([p for p in range(29, 32) if prime(p)] == [29, 31], 'lower prime boundary')
    require([p for p in range(113, 128) if prime(p)] == [113, 127], 'upper prime boundary')
    ell = 150
    B = M = 1155 * ell
    rt = lower_sqrt(M)
    require(F(400, 77 * ell) <= F(4, 3) and F(15, ell) <= F(4, 3), 'both1155 integration domains')
    costs = dict(periodic=phi_interval(F(200, 77 * ell))[1] + F(75, 77 * ell) + C1155 / rt, M_free=phi_interval(F(15, 2 * ell))[1] + F(45, 16 * ell))
    choices = {name: (F(6, 125) + cost) / log_interval(F(B))[0] for name, cost in costs.items()}
    selected = min(choices, key=choices.get)
    upper = choices[selected]
    require(selected == 'periodic' and upper < DELTA, 'tail1155 gate150')
    require([p for p in range(149, 152) if prime(p)] == [149, 151], 'first actual prime151')
    tail = dict(ell_scalar=ell, first_actual_prime=151, B_min=B, M_min=M, sqrt_lower=rt, cap=record(F(6, 125)), costs={k: record(v) for k, v in costs.items()}, choices={k: record(v) for k, v in choices.items()}, selected=selected, upper=record(upper), upper_decimal_ceiling=decimal_bound(upper), margin=record(DELTA - upper), margin_decimal_floor=decimal_bound(DELTA - upper, False), scope='rad(B)=1155*ell, ell>=150 prime', all_positive_valuations=True, all_N_ge_1=True)
    output = dict(status='PASS', normalization='Delta(B,N)=mu(B^N)/ln(B)', anchors=anchors, gates105=gates, tail1155=tail, small_primes=small, finite_primes=finite, finite_complement=dict(N=1, p_min=31, p_upper_exclusive=125, B_upper_exclusive=1000000), signed_cap_rule='max(0,cap+cost)/ln(B); positive boundary numerator checked', upper_values_are_not_values_of_Delta=True, exact_Cstar_proved=False, all_four_prime_supports_closed=False, all_five_prime_supports_closed=False, geometric_replay_here=False)
    (HERE / 'gates105_certificate.json').write_text(json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status='PASS', gates={v['name']: v['upper_decimal_ceiling'] for v in gates}, margins={v['name']: v['margin_decimal_floor'] for v in gates}, tail1155_upper=tail['upper_decimal_ceiling'], tail1155_margin=tail['margin_decimal_floor']), indent=2))
if __name__ == '__main__':
    main()
