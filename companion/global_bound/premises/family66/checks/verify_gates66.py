"""Exact gates66, tail858 and uniform6rp, using fresh frozen marked caps.

This checks numerical use and byte identities, not the profile geometries.
Every decisive logarithm/root/margin computation is rational or integral.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
import json
from exact_support import need, fs, log_bounds, phi_bounds
ROOT = Path(__file__).resolve().parent
DELTA = F(2, 85)
SHARP66_C = F(84371299, 4278582)
SHARP858_C = F(95676469, 2623920)

def fraction(v):
    return F(v['numerator'], v['denominator'])

def directed_decimal(x, up=True, places=15):
    x = F(x)
    scale = 10 ** places
    n = -(-x.numerator * scale // x.denominator) if up else x.numerator * scale // x.denominator
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n // scale}.{n % scale:0{places}d}'

def prime(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

@lru_cache(None)
def phi_upper(x):
    return phi_bounds(F(x))[1]

def root_lower(M):
    r = isqrt(M)
    need(r * r <= M < (r + 1) ** 2, 'directed integer square root')
    return r

def periodic_sharp66(p, M):
    need(p > 0 and F(100, 33 * p) <= F(4, 3), 'periodic sharp integration range66')
    return phi_upper(F(50, 33 * p)) + F(25, 44 * p) + SHARP66_C / root_lower(M)

def mfree_sharp66(p):
    need(p > 0 and F(12, p) <= F(4, 3), 'M-free integration range66')
    return phi_upper(F(6, p)) + F(9, 4 * p)

def bound66(p, B, M, cap):
    cap = F(cap)
    need(B > 1 and M >= B and (cap >= 0), 'positive normalization domain66')
    L = log_bounds(B)[0]
    choices = {'periodic': (cap + periodic_sharp66(p, M)) / L, 'M_free': (cap + mfree_sharp66(p)) / L}
    selected = min(choices, key=choices.get)
    return (choices, selected, choices[selected])

def marked_anchor(name, q, k, cap):
    path = ROOT / 'imports' / name
    raw = path.read_bytes()
    j = json.loads(raw)
    need((j['q'], j['k'], j['transfer_index']) == (q, k, q ** k), 'marked anchor domain')
    value = F(j['total_upper'])
    need(value < cap, 'fresh exact anchor below simple rational cap')
    return dict(file=path.relative_to(ROOT).as_posix(), q=q, k=k, divisor=q ** k, total_upper=fs(value), simple_cap=fs(cap), bytes=len(raw), sha256=sha256(raw).hexdigest(), every_multiple_including_extra_primes=True, geometric_replay_here=False)

def main():
    theory = json.loads((ROOT / 'theory66_certificate.json').read_text())
    need(theory['status'] == 'PASS', 'theory status')
    for q, r, A, C in [(6, 2, F(5, 6), F(13, 6)), (66, 4, F(25, 33), F(5645, 1617)), (858, 6, F(100, 143), F(633619, 120263))]:
        row = theory['marks'][str(q)]
        need(row['radius'] == r and fraction(row['A']) == A and (fraction(row['C']) == C), 'exact theory bridge')
    need(fraction(theory['penalties']['66']['sharp_sqrtM_coefficient']) == SHARP66_C, 'sharp66 bridge')
    need(fraction(theory['penalties']['858']['sharp_sqrtM_coefficient']) == SHARP858_C, 'sharp858 bridge')
    anchors = [marked_anchor('q66_k1.json', 66, 1, F(127, 1000)), marked_anchor('q66_k2.json', 66, 2, F(2, 125)), marked_anchor('q858.json', 858, 1, F(1, 8)), marked_anchor('q6.json', 6, 1, F(79, 200))]
    regions = []
    specs = [('large_p', 350, 23100, 23100, F(127, 1000), 'M_free'), ('N_ge_2', 47, 3102, 3102 ** 2, F(2, 125), 'periodic'), ('large_base', 47, 1000000, 1000000, F(127, 1000), 'periodic')]
    for label, p, B, M, cap, expected in specs:
        choices, selected, upper = bound66(p, B, M, cap)
        need(selected == expected and upper < DELTA, 'strict fixed gate66')
        margin = DELTA - upper
        regions.append(dict(name=label, p_scalar=p, B_min=B, M_min=M, cap=fs(cap), sqrt_lower=root_lower(M), selected=selected, choices={k: fs(v) for k, v in choices.items()}, upper=fs(upper), upper_decimal_ceiling=directed_decimal(upper), margin=fs(margin), margin_decimal_floor=directed_decimal(margin, False)))
    small = [p for p in range(13, 47) if prime(p)]
    finite = [p for p in range(47, 350) if prime(p)]
    need(small == [13, 17, 19, 23, 29, 31, 37, 41, 43] and len(finite) == 56, 'prime partition66')
    need([p for p in range(43, 48) if prime(p)] == [43, 47], 'lower prime gap66')
    need([p for p in range(349, 354) if prime(p)] == [349, 353], 'upper prime gap66')
    ell = 200
    B = M = 858 * ell
    rt = root_lower(M)
    periodic = phi_upper(F(300, 143 * ell)) + F(225, 286 * ell) + SHARP858_C / rt
    need(F(18, ell) <= F(4, 3), 'M-free range858')
    free = phi_upper(F(9, ell)) + F(27, 8 * ell)
    L = log_bounds(B)[0]
    candidates = {'periodic': (F(1, 8) + periodic) / L, 'M_free': (F(1, 8) + free) / L}
    selected = min(candidates, key=candidates.get)
    upper858 = candidates[selected]
    need(selected == 'periodic' and upper858 < DELTA, 'tail858 scalar200')
    need([p for p in range(199, 212) if prime(p)] == [199, 211], 'tail858 first prime')
    tail858 = dict(ell_scalar=ell, first_actual_prime=211, B_min=B, M_min=M, sqrt_lower=rt, cap=fs(F(1, 8)), selected=selected, choices={k: fs(v) for k, v in candidates.items()}, upper=fs(upper858), upper_decimal_ceiling=directed_decimal(upper858), margin=fs(DELTA - upper858), margin_decimal_floor=directed_decimal(DELTA - upper858, False), all_positive_valuations=True, all_N_ge_1=True, scope='rad(B)=858*ell, ell>=200 prime')
    r = p = 3000
    B = M = 6 * r * p
    rt = root_lower(M)
    first = phi_upper(F(5, 6 * r)) + F(5, 16 * r) + F(104663, 12096 * rt)
    second = phi_upper(F(5, 3 * p)) + F(5, 8 * p) + F(597848, 19845 * rt)
    need(F(5, 3 * r) <= F(4, 3) and F(10, 3 * p) <= F(4, 3), 'double insertion sharp ranges')
    upper6rp = (F(79, 200) + first + second) / log_bounds(B)[0]
    need(upper6rp < DELTA, 'uniform 6rp scalar3000')
    pair = [n for n in range(3000, 3020) if prime(n)][:2]
    need(pair == [3001, 3011], 'first two actual primes above uniform scalar')
    uniform = dict(r_scalar=3000, p_scalar=3000, actual_primes_required_distinct=True, ordered_convention='r<p', B_min=B, M_min=M, sqrt_lower=rt, cap6=fs(F(79, 200)), first_penalty_upper=fs(first), second_penalty_upper=fs(second), upper=fs(upper6rp), upper_decimal_ceiling=directed_decimal(upper6rp), margin=fs(DELTA - upper6rp), margin_decimal_floor=directed_decimal(DELTA - upper6rp, False), first_actual_pair=pair, all_positive_valuations=True, all_N_ge_1=True, scope='rad(B)=6*r*p, distinct primes r,p>=3000', intermediate_mark='6r', final_mark='6rp')
    output = dict(status='PASS', normalization='Delta(B,N)=mu(B^N)/ln(B)', anchors=anchors, gates66=regions, small_primes=small, finite_primes=finite, finite_complement=dict(N=1, p_min=47, p_upper_exclusive=350, B_upper_exclusive=1000000), tail858=tail858, uniform6rp=uniform, exact_Cstar_proved=False, all_four_prime_supports_closed=False, all_five_prime_supports_closed=False, upper_values_are_not_values_of_Delta=True, claim_scope='exact arithmetic with separately verified fresh marked-profile inputs')
    (ROOT / 'gates66_certificate.json').write_text(json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status='PASS', gate66_uppers={r['name']: r['upper_decimal_ceiling'] for r in regions}, gate66_margin_floors={r['name']: r['margin_decimal_floor'] for r in regions}, tail858_upper=tail858['upper_decimal_ceiling'], tail858_margin_floor=tail858['margin_decimal_floor'], uniform6rp_upper=uniform['upper_decimal_ceiling'], uniform6rp_margin_floor=uniform['margin_decimal_floor']), indent=2))
if __name__ == '__main__':
    main()
