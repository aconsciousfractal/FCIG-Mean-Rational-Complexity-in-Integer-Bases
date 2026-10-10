"""Integer-only verification of the new deduction from the imported Stevens bound."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb
import argparse, hashlib, json
HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise AssertionError(message)

def subtract(a, b):
    return [(a[i] if i < len(a) else F()) - (b[i] if i < len(b) else F()) for i in range(max(len(a), len(b)))]

def multiply(a, b):
    out = [F()] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out

def upper_decimal(x, digits=15):
    scale = 10 ** digits
    numerator = (x.numerator * scale + x.denominator - 1) // x.denominator
    return f'{numerator // scale}.{numerator % scale:0{digits}d}'

def lower_decimal(x, digits=15):
    scale = 10 ** digits
    numerator = x.numerator * scale // x.denominator
    return f'{numerator // scale}.{numerator % scale:0{digits}d}'

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'stevens2000_certificate.json')
    args = parser.parse_args()
    provenance = json.loads((HERE / 'PROVENANCE.json').read_text())
    for source in provenance['sources']:
        raw = (HERE / source['local_copy']).read_bytes()
        need(len(raw) == source['bytes'] and hashlib.sha256(raw).hexdigest() == source['sha256'], 'source byte bridge')
    classical = json.loads((HERE / 'imports/STEVENS_SOURCE_CHECK.json').read_text())
    need(classical['natural_logarithms'] and classical['primary_author_paper_used']['universal_gn_scope_explicit'], 'import explicitly covers all n with natural logarithms')
    need(classical['primary_author_paper_used']['identifier'] == 'arXiv:1208.5342v2' and classical['original_bibliography']['full_original_proof_retrieved'] is False, 'source version and original-proof limit')
    s4 = sum((F(1, factorial(i)) for i in range(5)), F())
    s6 = sum((F(1, factorial(i)) for i in range(7)), F())
    u6 = s6 + F(1, 6 * factorial(6))
    need(s4 == F(65, 24), 'exponential lower through fourth term')
    need(s6 == F(1957, 720), 'exponential lower through sixth term')
    need(u6 == F(11743, 4320) < F(1631, 600) < F(49, 18), 'exponential upper chain')
    gates = []
    for name, left, right, direction, conclusion in [('ln2_lower', 49 ** 69, 2 ** 100 * 18 ** 69, '<', 'ln2>69/100'), ('ln2_upper', 65 ** 7, 2 ** 10 * 24 ** 7, '>', 'ln2<7/10'), ('ln2000_lower', 1631 ** 38, 2000 ** 5 * 600 ** 38, '<', 'ln2000>38/5'), ('ln2000_upper', 1957 ** 761, 2000 ** 100 * 720 ** 761, '>', 'ln2000<761/100')]:
        need(left < right if direction == '<' else left > right, name)
        gates.append(dict(name=name, left=left, right=right, comparison=direction, conclusion=conclusion))
    need(F(69, 100) > F(2, 3), 'positive lower denominator gain')
    need(13 * 3 ** 2 > 87, 'H>=3 implies K_H<=H^2')
    need(29 * (3 ** 2 + 3) <= 42 * 3 ** 2, 'direct H3 polynomial check')
    A = [F(6, 5), F(2), F(49, 9)]
    Ap = [A[1], 2 * A[2]]
    Q = subtract(multiply(subtract(Ap, A), [-F(1, 3), F(1)]), A)
    expected = [-F(22, 15), -F(562, 135), F(142, 27), -F(49, 9)]
    need(Q == expected, 'differentiation polynomial identity')
    translated = [sum((Q[j] * comb(j, i) for j in range(i, len(Q))), F()) for i in range(len(Q))]
    need(translated == [-F(157, 27), -F(449, 45), -F(299, 27), -F(49, 9)], 'Q(1+u) exact coefficients')
    need(all((x < 0 for x in translated)), 'strict negativity for every u>=0')
    need(-F(49, 9) + F(142, 27) == -F(5, 27) < 0, 'independent grouping sign at t>=1')
    need(all((x > 0 for x in A)), 'numerator positive and increasing at positive t')
    tlo = F(38, 5)
    thi = F(761, 100)
    need(1 < tlo < thi and tlo - F(1, 3) > 0, 'G domain at lower bound2000')
    scalar = (A[0] + A[1] * thi + A[2] * thi ** 2) / (2000 * (tlo - F(1, 3)))
    need(scalar == F(29854729, 1308000000) < F(2, 85), 'main gate2000')
    margin = F(2, 85) - scalar
    need(margin == F(15669607, 22236000000) > 0, 'exact positive margin')
    need(F(upper_decimal(scalar)) >= scalar, 'outward upper decimal')
    need(F(lower_decimal(margin)) <= margin, 'outward lower margin')
    need(2048 == 2 ** 11, 'power-of-two crosscheck')
    alo = 11 * F(69, 100)
    ahi = 11 * F(7, 10)
    need((alo, ahi) == (F(759, 100), F(77, 10)), 'ln2048 interval from ln2')
    alt = (A[0] + A[1] * ahi + A[2] * ahi ** 2) / (2048 * (alo - F(1, 3)))
    need(alt == F(305461, 13375488) < F(2, 85), 'redundant gate2048')
    pupper = A[0] + A[1] * thi + A[2] * thi ** 2
    need(pupper == F(29854729, 90000), 'uniform mean cap below2000 prime factors')
    base_upper = pupper / (20480 * F(69, 100))
    need(base_upper == F(29854729, 1271808000) < F(2, 85), 'large-base corollary')
    base_margin = F(2, 85) - base_upper
    need(base_margin == F(1192807, 21620736000) > 0, 'large-base exact margin')
    result = dict(status='PASS_WITH_DECLARED_IMPORTS', checks=CHECKS, theorem=dict(minimum_distinct_prime_count=2000, all_larger_prime_counts=True, all_positive_valuations=True, all_N_ge1=True, normalization='mu(B^N)/ln(B)', exact_uniform_upper=str(scalar), strict_upper_decimal=upper_decimal(scalar), threshold='2/85', exact_margin=str(margin), margin_lower_decimal=lower_decimal(margin), strict_below_threshold=True, below_c0_uses_previous_candidate_comparison=True), sources=provenance['sources'], imported_Stevens=dict(statement='J(n)<=2*r^(2+2*e*ln(r)), r=omega(n)>=1', all_prime_supports=True, original_proof_reconstructed=False, not_assuming_primorial_extremality=True), imported_internal=['SV/W/MEAN with K_H=29*(H^2+3)/42', '2/85<c0 for candidate210'], e_bounds=dict(strict_lower=str(s6), series_upper=str(u6), upper_2000='1631/600', upper_numerator='49/18', lower_ln2='65/24'), integer_log_gates=gates, numerator_coefficients_ascending=list(map(str, A)), denominator_lower='ln(B)>k*(ln(k)-1/3)', derivative_polynomial_coefficients_ascending=list(map(str, Q)), derivative_at_one_plus_u_coefficients_ascending=list(map(str, translated)), monotonicity_scope='all real k>=e; proved by negative polynomial, not sampled', crosscheck2048=dict(exact_upper=str(alt), used_for_main_claim=False), large_base_corollary=dict(base_lower_bound='2^20480', all_prime_counts=True, all_positive_valuations=True, all_N_ge1=True, case_k_ge2000='main theorem', case_k_below2000_mean_upper=str(pupper), exact_uniform_upper=str(base_upper), threshold='2/85', exact_margin=str(base_margin), strict_below_threshold=True, necessary_condition_for_Delta_above_c0=dict(B_strictly_below='2^20480', omega_B_at_most=1999, N_still_unbounded=True), finite_enumeration_of_pairs_claimed=False, practical_enumeration_claimed=False), no_claim_for_remaining_5_to_1999_sectors=True, exact_Cstar_not_identified=True, geometric_profiles_required=False, geometry_replayed=False, finite_checks_replace_imported_or_infinite_proofs=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status=result['status'], checks=CHECKS, threshold_k=2000, upper=str(scalar), upper_decimal=upper_decimal(scalar), margin=str(margin)), indent=2))
if __name__ == '__main__':
    main()
