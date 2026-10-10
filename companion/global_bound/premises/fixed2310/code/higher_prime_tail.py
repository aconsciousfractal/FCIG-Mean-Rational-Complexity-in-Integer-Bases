"""Uniform improvement for every support with at least six primes.

Uses the existing SV/W/MEAN proof, an elementary H76 union bound for
seven remaining primes, and the declared classical Kanold tail import.
No geometric profile or numerical sampling is assumed.
"""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import argparse, json
from exact import need, log_lo, log_hi, out_decimal
ROOT = Path(__file__).resolve().parent.parent

def certificate():
    H = 76
    odd = (3, 5, 7, 11, 13, 17, 19)
    odd_cover = [(H + p - 1) // p for p in odd]
    n = H // 2
    even_cover = [(n + p - 1) // p for p in odd[:-1]]
    need(odd_cover == [26, 16, 11, 7, 6, 5, 4] and sum(odd_cover) == 75 < 76, 'seven odd primes leave a point in every76 consecutive integers')
    need(even_cover == [13, 8, 6, 4, 3, 3] and sum(even_cover) == 37 < 38, 'six odd primes leave an odd point in every76 consecutive integers')
    rows = []
    for k, B, h in [(6, 30030, 18), (7, 510510, 38), (8, 9699690, 76), (9, 223092870, 256)]:
        K = F(29 * (h * h + 3), 42)
        lo = (1 + log_lo(K)) / (2 * log_hi(B))
        hi = (1 + log_hi(K)) / (2 * log_lo(B))
        rows.append(dict(omega=k, base_min=B, H=h, K=str(K), lower=str(lo), upper=str(hi), upper_decimal=out_decimal(hi)))
    need(prod((2, 3, 5, 7, 11, 13, 17, 19, 23)) == 223092870, 'P9 exact primorial')
    K9 = F(1900631, 42)
    P9 = 223092870
    der = log_hi(4) * log_hi(P9) - log_lo(29) * (1 + log_lo(K9))
    need(der < 0, 'F9 decreasing on every real v>=0')
    l6 = F(rows[0]['lower'])
    u6 = F(rows[0]['upper'])
    need(all((F(r['upper']) < l6 for r in rows[1:])), 'six-prime bound dominates all k>=6')
    left = (49 * 3161) ** 8
    right = 30030 ** 5 * (18 * 14) ** 8
    need(left < right and u6 < F(5, 16), 'independent integer and logarithmic5/16 gates')
    return dict(status='PASS_WITH_DECLARED_IMPORTS', theorem='For every B with omega(B)>=6 and every N>=1, Delta(B,N)<=U6<5/16', exact_upper_expression='(1+log(3161/14))/(2*log(30030))', upper=str(u6), upper_decimal=out_decimal(u6), simple_strict_upper='5/16', all_positive_valuations=True, all_N_ge_1=True, all_prime_counts_ge_6=True, below_c0_exclusion_claimed=False, new_geometric_profile_used=False, H76=dict(max_remaining_primes=7, block=76, odd_cover=odd_cover, odd_sum=75, odd_subsequence_length=38, even_cover=even_cover, even_sum=37, optimal_Jacobsthal_constant_claimed=False), sectors=rows, infinite_tail=dict(omega_from=9, P9=P9, K9=str(K9), v='omega(B)-9', base_growth_lower='29^v', K_growth_upper='4^v', derivative_upper=str(der), decreasing_for_every_real_v_ge_0=True, finite_k_scan_is_not_the_proof=True), integer_gate=dict(expression='(49*3161)^8 < 30030^5*(18*14)^8', left=str(left), right=str(right), positive_difference=str(right - left)), imports=['SV/W/MEAN, sources/GLOBAL_CAP_SOURCE.md sections3-4', 'elementary H18 andH38 bounds, same source section2', 'Kanold J(q)<=2^omega(q), only for omega(B)>=9'], formalized=False)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--check', action='store_true')
    a = p.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    f = ROOT / 'results/higher_prime_tail.json'
    if a.check:
        need(f.read_text() == s, 'frozen high-prime certificate')
    else:
        f.write_text(s, newline='\n')
    print(json.dumps({k: z[k] for k in ('status', 'upper_decimal', 'simple_strict_upper')}, indent=2))
