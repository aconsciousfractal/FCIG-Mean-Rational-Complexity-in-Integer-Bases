"""Exact finite gate and monotonicity for the whole omega(B)>=2000 tail.

Stevens' universal Jacobsthal bound is a declared classical import. This
certificate checks the new deduction, not the proof of the imported theorem.
No floating point value decides a theorem or inequality.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import argparse, json
from exact import need, log_lo, log_hi, out_decimal
ROOT = Path(__file__).resolve().parent.parent

def mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c

def certificate():
    lower = sum((F(1, factorial(i)) for i in range(7)), F(0))
    upper = sum((F(1, factorial(i)) for i in range(5)), F(0)) + F(1, 120) * F(6, 5)
    need(lower == F(1957, 720) and upper == F(1631, 600) < F(49, 18), 'partial exponential series and strict geometric remainder')
    comparisons = [('ln2 > 69/100', 49 ** 69, 2 ** 100 * 18 ** 69), ('ln2 < 7/10', 2 ** 10 * 24 ** 7, 65 ** 7), ('ln2000 > 38/5', 1631 ** 38, 2000 ** 5 * 600 ** 38), ('ln2000 < 761/100', 2000 ** 100 * 720 ** 761, 1957 ** 761)]
    gates = []
    for name, left, right in comparisons:
        need(left < right, name)
        gates.append(dict(claim=name, left=str(left), right=str(right), positive_difference=str(right - left)))
    P = [F(6, 5), F(2), F(49, 9)]
    D = mul([P[1], 2 * P[2]], [-F(1, 3), F(1)])
    Q = mul(P, [F(2, 3), F(1)])
    D += [F(0)] * (len(Q) - len(D))
    derivative = [a - b for a, b in zip(D, Q)]
    need(derivative == [-F(22, 15), -F(562, 135), F(142, 27), -F(49, 9)], 'symbolic polynomial numerator of the derivative')
    need(derivative[3] + derivative[2] == -F(5, 27) and derivative[0] < 0 and (derivative[1] < 0), 'derivative is negative for every real t>=1')
    tlo = F(38, 5)
    thi = F(761, 100)
    gamma = (P[0] + P[1] * thi + P[2] * thi ** 2) / (2000 * (tlo - F(1, 3)))
    need(gamma == F(29854729, 1308000000) < F(2, 85), 'uniform exact 2000-prime exclusion')
    need(tlo < log_lo(2000) < log_hi(2000) < thi, 'independent atanh logarithm check')
    tu = log_hi(2000)
    tl = log_lo(2000)
    tighter = (P[0] + P[1] * tu + P[2] * tu ** 2) / (2000 * (tl - F(1, 3)))
    need(tighter < gamma, 'separate endpoint enclosure contained in integer gate')
    gamma2048 = (P[0] + P[1] * F(77, 10) + P[2] * F(77, 10) ** 2) / (2048 * (F(759, 100) - F(1, 3)))
    need(gamma2048 == F(305461, 13375488) < F(2, 85), 'additional 2048 integer endpoint')
    mean_all_small_k = P[0] + P[1] * thi + P[2] * thi ** 2
    size_gamma = mean_all_small_k / (20480 * F(69, 100))
    need(mean_all_small_k == F(29854729, 90000) and size_gamma == F(29854729, 1271808000) < F(2, 85), 'every base at least 2^20480 is excluded, using the two prime-count regions')
    return dict(status='PASS_WITH_DECLARED_IMPORTS', theorem='For every base B with omega(B)>=2000 and every N>=1, Delta(B,N)<29854729/1308000000<2/85<c0', omega_from=2000, upper=str(gamma), upper_decimal=out_decimal(gamma), strict_upper=True, delta='2/85', strict_margin_to_delta=str(F(2, 85) - gamma), all_prime_supports=True, all_positive_valuations=True, all_N_ge_1=True, normalization='Delta(B,N)=mu(B^N)/ln(B)', auxiliary_function='G(k)=(6/5+2*ln(k)+(49/9)*ln(k)^2)/(k*(ln(k)-1/3))', monotonicity=dict(variable='t=ln(k)', domain='every real t>=1', derivative_numerator_coefficients_ascending=[str(v) for v in derivative], cubic_and_quadratic_upper='-(5/27)*t^2', strict_decrease_on_entire_half_line=True, finite_scan_not_the_proof=True), exponential_bounds=dict(lower=str(lower), upper=str(upper), upper_tail='sum(n=5..infinity)1/n! < (1/120)/(1-1/6)'), integer_gates=gates, independent_log_endpoint_upper=str(tighter), additional_2048_endpoint_upper=str(gamma2048), all_large_bases=dict(base_lower='2^20480', all_N_ge_1=True, omega_ge_2000_route='main infinite-tail theorem', omega_below_2000_route_mean_upper=str(mean_all_small_k), uniform_upper=str(max(gamma, size_gamma)), uniform_upper_decimal=out_decimal(max(gamma, size_gamma)), strict_margin_to_delta=str(F(2, 85) - max(gamma, size_gamma)), zero_remaining_primes_handled_by_J1_equals_1=True, finite_base_region_for_any_competitor='2<=B<2^20480 and omega(B)<=1999', N_still_unbounded=True, practical_exhaustive_scan_claimed=False), prime_product_lower='ln(B)>=ln(2)+k*(ln(2*k)-1)>k*(ln(k)-1/3)', primorial_maximality_assumed=False, external_import='Stevens: J(n)<=2*r^(2+2*e*ln(r)), r=omega(n)>=1, every n', external_theorem_proof_reconstructed=False, internal_import='SV/W/MEAN: mu(B^N)<=(1+ln(29*(H^2+3)/42))/2', source_check='sources/STEVENS_SOURCE_CHECK.json', new_geometric_profile_used=False, all_five_prime_supports_excluded=False, exact_global_Cstar_proved=False, formalized=False)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    z = certificate()
    s = json.dumps(z, indent=2) + '\n'
    p = ROOT / 'results/stevens_tail.json'
    if a.check:
        need(p.read_text() == s, 'frozen Stevens tail certificate')
    else:
        p.write_text(s, newline='\n')
    print(json.dumps({k: z[k] for k in ['status', 'omega_from', 'upper_decimal', 'strict_margin_to_delta']}, indent=2))
