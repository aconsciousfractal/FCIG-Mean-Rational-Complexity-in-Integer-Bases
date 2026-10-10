"""Integer global gate C*<=U2730<1/3, with explicit imported premises.

Only the fresh K=q marked2310 cap is needed. Exclusion of that radical
below2/85, higher-index profiles, and previous210 numeric caps are unused.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import prod, factorial, isqrt
import argparse, hashlib, json
HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(v, message):
    global CHECKS
    CHECKS += 1
    if not v:
        raise AssertionError(message)

def prime(n):
    return n >= 2 and all((n % d for d in range(2, isqrt(n) + 1)))

def upper_gate(K, B):
    left = (49 * K.numerator) ** 50
    right = B ** 33 * (18 * K.denominator) ** 50
    need(left < right, 'strict integer sector gate below33/100')
    return dict(K=str(K), B_min=B, strict_upper='33/100', integer_left=left, integer_right=right, positive_integer_margin=right - left, inequality='(49*Knum)^50 < B^33*(18*Kden)^50')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE / 'global_third_certificate.json')
    args = parser.parse_args()
    provenance = json.loads((HERE / 'PROVENANCE.json').read_text())
    inputs = []
    for name in ['q2310.json', 'PREVIOUS_GLOBAL_CAP_PROOF.md']:
        path = HERE / 'imports' / name
        raw = path.read_bytes()
        row = next((x for x in provenance['sources'] if x['local_copy'] == 'imports/' + name))
        need(len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'], 'import/source hash')
        inputs.append(dict(file='imports/' + name, bytes=len(raw), sha256=row['sha256']))
    anchor = json.loads((HERE / 'imports/q2310.json').read_text())
    A = F(anchor['total_upper'])
    need((anchor['q'], anchor['k'], anchor['transfer_index']) == (2310, 1, 2310), 'fresh anchor domain')
    need(0 < A < F(4, 15), 'fresh marked2310 cap below4/15')
    e_upper = sum((F(1, factorial(n)) for n in range(9)), F(0)) + F(1, 8 * factorial(8))
    need(F(8, 3) < e_upper < F(49, 18), 'rational exponential sandwich')
    need(sum((F(1, factorial(n)) for n in range(4)), F(0)) == F(8, 3), 'strict partial sum lower')
    need(F(2310) > F(49, 18) ** 7, 'ln2310>7')
    need(F(4, 105) < F(33, 100), 'whole2310 family below separator')
    first = sorted(((prod(c), c) for c in combinations([p for p in range(2, 14) if prime(p)], 5)))
    need(first[:2] == [(2310, (2, 3, 5, 7, 11)), (2730, (2, 3, 5, 7, 13))], 'first two five-prime radicals')
    need(2730 // (2 * 3 * 5 * 7) == 13, 'all possible fifth primes at the boundary included')
    small_four = [prod(c) for c in combinations([p for p in range(2, 26) if prime(p)], 4) if prod(c) < 770]
    small_four.sort()
    need(small_four == [210, 330, 390, 462, 510, 546, 570, 690, 714], 'four-prime frontier770')
    need(all((n % 30 == 0 and prime(n // 30) and (n // 30 >= 7) or (n % 42 == 0 and prime(n // 42) and (n // 42 >= 11)) for n in small_four)), 'imported30p/42p cover frontier')
    need(769 // 30 == 25, 'four-prime enumeration complete')
    sectors = {}
    for k, B, H in [(4, 770, 6), (6, 30030, 18), (7, 510510, 38), (8, 9699690, 128)]:
        K = F(29 * (H * H + 3), 42)
        sectors[str(k)] = upper_gate(K, B)
    K8 = F(475223, 42)
    P8 = 9699690
    need(F(4) ** 2 < F(8, 3) ** 3, 'ln4<3/2')
    need(F(P8) < F(8, 3) ** 17, 'lnP8<17')
    need(F(23) > F(49, 18) ** 3, 'ln23>3')
    need(K8 > F(49, 18) ** 9, 'lnK8>9')
    derivative_upper = F(3, 2) * 17 - 3 * (1 + 9)
    need(derivative_upper == -F(9, 2) < 0, 'F8 derivative strictly negative')
    need(prod([2, 3, 5, 7, 11, 13, 17, 19]) == P8 and prime(23), 'tail minimum product and added-prime floor')
    left = (49 * 2987) ** 3
    right = 2730 ** 2 * (18 * 42) ** 3
    need(left < right, 'U2730<one third by integers')
    lower_left = (8 * 2987) ** 50
    lower_right = 2730 ** 33 * (3 * 42) ** 50
    need(lower_left > lower_right, 'U2730>33/100 by integers')
    need(F(2, 85) < F(1, 40) < F(33, 100), 'imported small-sector separators')
    output = dict(status='PASS_WITH_DECLARED_IMPORTS', checks=CHECKS, inputs=inputs, theorem='C* <= (1+ln(2987/42))/(2ln2730) < 1/3', exact_global_Cstar_proved=False, strong_delta_exclusion_of_2310_required=False, strong_delta_exclusion_of_2310_claimed=False, higher_index_profiles_used=False, anchor2310=dict(exact_cap=str(A), simple_cap='4/15', divisor=2310, all_positive_multiples=True, full_radical_family_normalized_upper='4/105', all_positive_valuations=True, all_N_ge_1=True), e_upper=str(e_upper), e_relaxation='49/18', e_lower='8/3', separator='33/100', five_prime_boundary=dict(first=2310, next=2730, next_primes=[2, 3, 5, 7, 13], no_five_prime_radical_between=True), residual_four_prime_frontier=770, residual_four_prime_radicals_below=small_four, sector_gates=sectors, tail8=dict(base=P8, K=str(K8), scaling_B='23^v', scaling_K='4^v', derivative_strict_upper=str(derivative_upper), every_real_v_ge_0=True, classical_Kanold_import=True), U2730_lower_gate=dict(strict_lower='33/100', integer_left=lower_left, integer_right=lower_right, comparison='left>right'), U2730_upper_gate=dict(strict_upper='1/3', integer_left=left, integer_right=right, positive_integer_margin=right - left, comparison='left<right'), imports_declared=['omega(B)<=2: Delta<=0', 'omega(B)=3: Delta<=2/85', 'complete30p theorem, candidate210 exception with c0<1/40', 'complete42p theorem for p>=11', 'SV/W/MEAN and elementary Jacobsthal bounds6,10,18,38', 'Kanold J(q)<=2^omega(q) only for k>=8', 'fresh marked2310 profile valid on every multiple2310'], previous70_66_105_closures_required=False, geometry_replayed_here=False, imported_theorem_replays_performed_here=False, claim_classes='new internal raccordo and integer checks; historical and classical imports explicit')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps(dict(status=output['status'], checks=CHECKS, theorem=output['theorem'], strong2310_exclusion_required=False, integer_margin=right - left), indent=2))
if __name__ == '__main__':
    main()
